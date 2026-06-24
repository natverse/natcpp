#include <Rcpp.h>
#include <algorithm>
#include <atomic>
#include <memory>
#include <vector>
#include <RcppThread.h>

using namespace Rcpp;

// Progress bar that suppresses output for the first `delaySecs` seconds,
// then shows a bar with time estimate (using RcppThread's Rcout).
class DelayedProgressBar : public RcppThread::ProgressPrinter {
public:
  DelayedProgressBar(size_t numIt, size_t printEvery, size_t delaySecs,
                     bool display = true)
    : ProgressPrinter(numIt, printEvery), delaySecs_(delaySecs),
      display_(display) {}

private:
  size_t delaySecs_;
  bool display_;
  bool started_{false};

  void printProgress() override {
    if (!display_ || isDone_) return;
    if (it_ == numIt_) isDone_ = true;

    using namespace std::chrono;
    auto elapsed = duration<float>(steady_clock::now() - startTime_).count();
    if (!started_ && elapsed < delaySecs_) return;
    started_ = true;

    double pct = std::round(it_ * 100.0 / numIt_);
    std::ostringstream msg;
    msg << "\rComputing: " << makeBar(pct) << progressString();
    RcppThread::Rcout << msg.str();
  }

  std::string makeBar(size_t pct, size_t numBars = 40) {
    std::ostringstream msg;
    msg << "[";
    size_t i = 0;
    for (; i < pct / 100.0 * numBars; i++) msg << "=";
    for (; i < numBars; i++) msg << " ";
    msg << "] ";
    return msg.str();
  }
};

//' @noRd
// [[Rcpp::export]]
NumericVector weighted_jaccard_sparse_fill(
    const S4& x, const S4& pattern, bool transpose = false,
    bool display_progress = true, int threads = 4) {

  // Original matrix slots (dgCMatrix: nr x nc)
  IntegerVector dims = x.slot("Dim");
  IntegerVector xi = x.slot("i");
  IntegerVector xp = x.slot("p");
  NumericVector xx = x.slot("x");

  const int nr = dims[0];
  const int nc = dims[1];
  const int nfeat = transpose ? nc : nr;

  // Pattern matrix slots (dgCMatrix: ncomp x ncomp)
  IntegerVector Ai = pattern.slot("i");
  IntegerVector Ap = pattern.slot("p");
  IntegerVector Adims = pattern.slot("Dim");
  const int ncomp = Adims[0];
  const int annz = Ai.size();

  // Output: min_sums values aligned with pattern's (i, p) structure
  NumericVector out(annz, 0.0);

  // Build CSR for feature dimension: for each feature f, the items (columns
  // or rows being compared) that are nonzero in that feature, with their values.
  std::vector<int> feat_count(nfeat, 0);
  if (!transpose) {
    for (int col = 0; col < nc; ++col)
      for (int idx = xp[col]; idx < xp[col + 1]; ++idx)
        feat_count[xi[idx]]++;
  } else {
    for (int col = 0; col < nc; ++col)
      feat_count[col] = xp[col + 1] - xp[col];
  }

  std::vector<int> feat_offsets(nfeat + 1, 0);
  for (int f = 0; f < nfeat; ++f)
    feat_offsets[f + 1] = feat_offsets[f] + feat_count[f];

  const int total_nnz = feat_offsets[nfeat];
  std::vector<int> feat_items(total_nnz);   // which item (col/row being compared)
  std::vector<double> feat_vals(total_nnz);  // value of x at that position
  std::vector<int> feat_fill(feat_offsets.begin(), feat_offsets.end());

  if (!transpose) {
    for (int col = 0; col < nc; ++col) {
      for (int idx = xp[col]; idx < xp[col + 1]; ++idx) {
        const int row = xi[idx];
        const int dest = feat_fill[row]++;
        feat_items[dest] = col;
        feat_vals[dest] = xx[idx];
      }
    }
  } else {
    for (int col = 0; col < nc; ++col) {
      for (int idx = xp[col]; idx < xp[col + 1]; ++idx) {
        const int dest = feat_fill[col]++;
        feat_items[dest] = xi[idx];
        feat_vals[dest] = xx[idx];
      }
    }
  }

  // Also need item->features CSR: for each compared item, which features
  // it participates in and with what value.
  std::vector<int> item_count(ncomp, 0);
  for (int f = 0; f < nfeat; ++f) {
    for (int idx = feat_offsets[f]; idx < feat_offsets[f + 1]; ++idx)
      item_count[feat_items[idx]]++;
  }

  std::vector<int> item_offsets(ncomp + 1, 0);
  for (int c = 0; c < ncomp; ++c)
    item_offsets[c + 1] = item_offsets[c] + item_count[c];

  std::vector<int> item_feats(total_nnz);   // which feature
  std::vector<double> item_vals(total_nnz);  // value
  std::vector<int> item_fill(item_offsets.begin(), item_offsets.end());

  for (int f = 0; f < nfeat; ++f) {
    for (int idx = feat_offsets[f]; idx < feat_offsets[f + 1]; ++idx) {
      const int item = feat_items[idx];
      const int dest = item_fill[item]++;
      item_feats[dest] = f;
      item_vals[dest] = feat_vals[idx];
    }
  }

  // Raw pointer for thread-safe access (Rcpp NumericVector not thread-safe)
  double* out_ptr = out.begin();
  const int* Ai_ptr = &Ai[0];
  const int* Ap_ptr = &Ap[0];
  const int* item_feats_ptr = item_feats.data();
  const double* item_vals_ptr = item_vals.data();
  const int* feat_items_ptr = feat_items.data();
  const double* feat_vals_ptr = feat_vals.data();
  const int* feat_offsets_ptr = feat_offsets.data();
  const int* item_offsets_ptr = item_offsets.data();

  DelayedProgressBar bar(ncomp, 1, 2, display_progress);

  const size_t nThreads = (threads > 0) ? static_cast<size_t>(threads)
                                       : std::thread::hardware_concurrency();

  // Pre-allocate per-thread row_to_pos vectors (freed when function exits)
  std::vector<std::vector<int>> per_thread_r2p(nThreads, std::vector<int>(ncomp, -1));

  // Thread-local index into per_thread_r2p, assigned once per thread via atomic
  static thread_local int tl_tid = -1;
  static thread_local int tl_gen = -1;
  static std::atomic<int> generation{0};
  const int cur_gen = ++generation;
  std::atomic<int> next_tid{0};

  // Parallel over output columns — each column writes to its own slice of out[]
  RcppThread::parallelFor(0, ncomp, [&](int cb) {
    if (tl_gen != cur_gen) {
      tl_tid = next_tid++;
      tl_gen = cur_gen;
    }
    auto& row_to_pos = per_thread_r2p[tl_tid];

    // Populate row_to_pos for pattern column cb
    for (int idx = Ap_ptr[cb]; idx < Ap_ptr[cb + 1]; ++idx)
      row_to_pos[Ai_ptr[idx]] = idx;

    // For each feature f where cb is nonzero
    for (int fi = item_offsets_ptr[cb]; fi < item_offsets_ptr[cb + 1]; ++fi) {
      const int f = item_feats_ptr[fi];
      const double vb = item_vals_ptr[fi];

      // For each other item ca also in this feature, accumulate min
      for (int idx = feat_offsets_ptr[f]; idx < feat_offsets_ptr[f + 1]; ++idx) {
        const int ca = feat_items_ptr[idx];
        const int pos = row_to_pos[ca];
        if (pos >= 0)
          out_ptr[pos] += std::min(feat_vals_ptr[idx], vb);
      }
    }

    // Clear row_to_pos
    for (int idx = Ap_ptr[cb]; idx < Ap_ptr[cb + 1]; ++idx)
      row_to_pos[Ai_ptr[idx]] = -1;

    ++bar;
  }, nThreads);

  return out;
}

// Helper: build feature CSR and totals from dgCMatrix slots.
// Shared by both dense matrix and dense triangle functions.
struct FeatureCSR {
  std::vector<double> totals;
  std::vector<int> feat_offsets;
  std::vector<int> feat_items;
  std::vector<double> feat_vals;
  int ncomp;
  int nfeat;
  int total_nnz;
};

static FeatureCSR build_feature_csr(const IntegerVector& xi, const IntegerVector& xp,
                                    const NumericVector& xx, int nr, int nc,
                                    bool transpose) {
  FeatureCSR csr;
  csr.ncomp = transpose ? nr : nc;
  csr.nfeat = transpose ? nc : nr;
  csr.totals.resize(csr.ncomp, 0.0);

  std::vector<int> feat_count(csr.nfeat, 0);
  if (!transpose) {
    for (int col = 0; col < nc; ++col) {
      for (int idx = xp[col]; idx < xp[col + 1]; ++idx) {
        csr.totals[col] += xx[idx];
        feat_count[xi[idx]]++;
      }
    }
  } else {
    for (int col = 0; col < nc; ++col) {
      feat_count[col] = xp[col + 1] - xp[col];
      for (int idx = xp[col]; idx < xp[col + 1]; ++idx) {
        csr.totals[xi[idx]] += xx[idx];
      }
    }
  }

  csr.feat_offsets.resize(csr.nfeat + 1, 0);
  for (int f = 0; f < csr.nfeat; ++f)
    csr.feat_offsets[f + 1] = csr.feat_offsets[f] + feat_count[f];

  csr.total_nnz = csr.feat_offsets[csr.nfeat];
  csr.feat_items.resize(csr.total_nnz);
  csr.feat_vals.resize(csr.total_nnz);
  std::vector<int> feat_fill(csr.feat_offsets.begin(), csr.feat_offsets.end());

  if (!transpose) {
    for (int col = 0; col < nc; ++col) {
      for (int idx = xp[col]; idx < xp[col + 1]; ++idx) {
        const int row = xi[idx];
        const int dest = feat_fill[row]++;
        csr.feat_items[dest] = col;
        csr.feat_vals[dest] = xx[idx];
      }
    }
  } else {
    for (int col = 0; col < nc; ++col) {
      for (int idx = xp[col]; idx < xp[col + 1]; ++idx) {
        const int dest = feat_fill[col]++;
        csr.feat_items[dest] = xi[idx];
        csr.feat_vals[dest] = xx[idx];
      }
    }
  }

  return csr;
}

//' Dense weighted Jaccard similarity via C++
//'
//' Compute the full weighted Jaccard similarity matrix for a dgCMatrix using
//' an adaptive dense accumulation strategy. For small output matrices, uses a
//' feature-oriented loop; for large outputs, switches to a column-oriented
//' loop for better cache performance.
//'
//' @param x A dgCMatrix (sparse column-compressed matrix)
//' @param transpose If \code{FALSE}, compare columns; if \code{TRUE}, compare
//'   rows
//' @param threads Number of threads (default 4). Set to 0 for all cores.
//' @param triangle If \code{TRUE}, return only the lower triangle as a flat
//'   numeric vector in \code{\link{dist}} layout. If \code{FALSE} (default),
//'   return a full square matrix.
//' @param distance If \code{TRUE}, return distance (\code{1 - similarity})
//'   instead of similarity. Default \code{FALSE}.
//' @return A dense numeric similarity matrix, or a numeric vector in dist
//'   layout when \code{triangle = TRUE}.
//' @export
// [[Rcpp::export]]
SEXP c_weighted_jaccard_dense(const S4& x, bool transpose = false,
                              int threads = 4, bool triangle = false,
                              bool distance = false) {
  IntegerVector dims = x.slot("Dim");
  IntegerVector xi = x.slot("i");
  IntegerVector xp = x.slot("p");
  NumericVector xx = x.slot("x");

  const int nr = dims[0];
  const int nc = dims[1];

  FeatureCSR csr = build_feature_csr(xi, xp, xx, nr, nc, transpose);
  const int ncomp = csr.ncomp;
  const int nfeat = csr.nfeat;
  const size_t nT = (threads > 0) ? static_cast<size_t>(threads)
                                   : std::thread::hardware_concurrency();

  if (triangle) {
    // --- Triangle path: output is a flat vector in dist layout ---
    // dist position for pair (r,c) with r>c (1-based):
    //   pos = (c-1)*n - c*(c-1)/2 + (r-c)   [1-based]
    // We use 0-based: for (ca, cb) with ca > cb:
    //   pos = cb*n - cb*(cb+1)/2 + (ca-cb) - 1  [0-based]
    // which simplifies to: cb*(ncomp-1) - cb*(cb-1)/2 + (ca-1)  ... let's use a
    // direct formula.

    const std::size_t dist_len = static_cast<std::size_t>(ncomp) * (ncomp - 1) / 2;
    NumericVector out(dist_len, 0.0);
    double* out_ptr = out.begin();
    const double* totals_ptr = csr.totals.data();

    // Build item->features CSR (always needed for colwise triangle)
    std::vector<int> item_count(ncomp, 0);
    for (int f = 0; f < nfeat; ++f)
      for (int idx = csr.feat_offsets[f]; idx < csr.feat_offsets[f + 1]; ++idx)
        item_count[csr.feat_items[idx]]++;

    std::vector<int> item_offsets(ncomp + 1, 0);
    for (int c = 0; c < ncomp; ++c)
      item_offsets[c + 1] = item_offsets[c] + item_count[c];

    std::vector<int> item_feats(csr.total_nnz);
    std::vector<double> item_vals(csr.total_nnz);
    std::vector<int> item_fill(item_offsets.begin(), item_offsets.end());

    for (int f = 0; f < nfeat; ++f) {
      for (int idx = csr.feat_offsets[f]; idx < csr.feat_offsets[f + 1]; ++idx) {
        const int item = csr.feat_items[idx];
        const int dest = item_fill[item]++;
        item_feats[dest] = f;
        item_vals[dest] = csr.feat_vals[idx];
      }
    }

    const int* feat_items_ptr = csr.feat_items.data();
    const double* feat_vals_ptr = csr.feat_vals.data();
    const int* feat_offsets_ptr = csr.feat_offsets.data();
    const int* item_feats_ptr = item_feats.data();
    const double* item_vals_ptr = item_vals.data();
    const int* item_offsets_ptr = item_offsets.data();

    // For column cb, accumulate min-sums for rows ca > cb only,
    // writing directly into dist layout.
    // dist index for (ca, cb) with ca > cb (0-based items, 0-based output):
    //   pos = cb * ncomp - cb*(cb+1)/2 + (ca - cb - 1)
    RcppThread::parallelFor(0, ncomp - 1, [=](int cb) {
      const std::size_t base = static_cast<std::size_t>(cb) * ncomp
                             - static_cast<std::size_t>(cb) * (cb + 1) / 2
                             - static_cast<std::size_t>(cb) - 1;
      // Thread-local accumulator for this column's min-sums (only rows > cb)
      std::vector<double> col_ms(ncomp, 0.0);

      for (int fi = item_offsets_ptr[cb]; fi < item_offsets_ptr[cb + 1]; ++fi) {
        const int f = item_feats_ptr[fi];
        const double vb = item_vals_ptr[fi];
        for (int idx = feat_offsets_ptr[f]; idx < feat_offsets_ptr[f + 1]; ++idx) {
          const int ca = feat_items_ptr[idx];
          if (ca > cb)
            col_ms[ca] += std::min(feat_vals_ptr[idx], vb);
        }
      }

      // Convert to similarity or distance and store
      const double tb = totals_ptr[cb];
      for (int ca = cb + 1; ca < ncomp; ++ca) {
        const double ms = col_ms[ca];
        const double denom = totals_ptr[ca] + tb - ms;
        double sim = (denom > 0.0) ? ms / denom : 0.0;
        out_ptr[base + ca] = distance ? (1.0 - sim) : sim;
      }
    }, nT);

    // Return as a dist object
    out.attr("class") = "dist";
    out.attr("Size") = ncomp;
    out.attr("Diag") = false;
    out.attr("Upper") = false;
    return out;

  } else {
    // --- Full matrix path (existing) ---
    NumericMatrix out(ncomp, ncomp);

    // Threshold: ncomp^2 * 8 bytes > ~12MB (typical L3 share per core)
    const bool use_colwise = (static_cast<std::size_t>(ncomp) * ncomp * 8 > 12 * 1024 * 1024);

    if (use_colwise) {
      std::vector<int> item_count(ncomp, 0);
      for (int f = 0; f < nfeat; ++f)
        for (int idx = csr.feat_offsets[f]; idx < csr.feat_offsets[f + 1]; ++idx)
          item_count[csr.feat_items[idx]]++;

      std::vector<int> item_offsets(ncomp + 1, 0);
      for (int c = 0; c < ncomp; ++c)
        item_offsets[c + 1] = item_offsets[c] + item_count[c];

      std::vector<int> item_feats(csr.total_nnz);
      std::vector<double> item_vals(csr.total_nnz);
      std::vector<int> item_fill(item_offsets.begin(), item_offsets.end());

      for (int f = 0; f < nfeat; ++f) {
        for (int idx = csr.feat_offsets[f]; idx < csr.feat_offsets[f + 1]; ++idx) {
          const int item = csr.feat_items[idx];
          const int dest = item_fill[item]++;
          item_feats[dest] = f;
          item_vals[dest] = csr.feat_vals[idx];
        }
      }

      double* out_ptr = &out[0];
      const int* feat_items_ptr = csr.feat_items.data();
      const double* feat_vals_ptr = csr.feat_vals.data();
      const int* feat_offsets_ptr = csr.feat_offsets.data();
      const int* item_feats_ptr = item_feats.data();
      const double* item_vals_ptr = item_vals.data();
      const int* item_offsets_ptr = item_offsets.data();

      RcppThread::parallelFor(0, ncomp, [&](int cb) {
        double* col = out_ptr + static_cast<std::size_t>(cb) * ncomp;
        for (int fi = item_offsets_ptr[cb]; fi < item_offsets_ptr[cb + 1]; ++fi) {
          const int f = item_feats_ptr[fi];
          const double vb = item_vals_ptr[fi];
          for (int idx = feat_offsets_ptr[f]; idx < feat_offsets_ptr[f + 1]; ++idx) {
            col[feat_items_ptr[idx]] += std::min(feat_vals_ptr[idx], vb);
          }
        }
      }, nT);
    } else {
      // Feature-oriented: iterate features, scatter into output
      for (int f = 0; f < nfeat; ++f) {
        const int start = csr.feat_offsets[f];
        const int end = csr.feat_offsets[f + 1];
        const int k = end - start;

        for (int a = 0; a < k; ++a) {
          const int ca = csr.feat_items[start + a];
          const double va = csr.feat_vals[start + a];
          out[static_cast<std::size_t>(ca) * ncomp + ca] += va;
          for (int b = a + 1; b < k; ++b) {
            const int cb = csr.feat_items[start + b];
            const double mv = std::min(va, csr.feat_vals[start + b]);
            out[static_cast<std::size_t>(ca) * ncomp + cb] += mv;
            out[static_cast<std::size_t>(cb) * ncomp + ca] += mv;
          }
        }
      }
    }

    // Convert min_sums to similarity (or distance)
    {
      double* out_ptr = &out[0];
      const double* totals_ptr = csr.totals.data();
      RcppThread::parallelFor(0, ncomp, [=](int cb) {
        for (int ca = 0; ca < ncomp; ++ca) {
          const std::size_t pos = static_cast<std::size_t>(cb) * ncomp + ca;
          if (ca == cb) {
            out_ptr[pos] = distance ? 0.0 : 1.0;
          } else {
            const double ms = out_ptr[pos];
            const double denom = totals_ptr[ca] + totals_ptr[cb] - ms;
            double sim = (denom > 0.0) ? ms / denom : 0.0;
            out_ptr[pos] = distance ? (1.0 - sim) : sim;
          }
        }
      }, nT);
    }

    return out;
  }
}
