#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <progress.hpp>

using namespace Rcpp;

//' @noRd
// [[Rcpp::export]]
NumericVector weighted_jaccard_sparse_fill(
    const S4& x, const S4& pattern, bool transpose = false,
    bool display_progress = true) {

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

  // row_to_pos[row] = position in out[] for entry (row, current_col).
  // Populated once per output column, then cleared.
  std::vector<int> row_to_pos(ncomp, -1);

  Progress prog(ncomp, display_progress);

  // Iterate over output columns
  for (int cb = 0; cb < ncomp; ++cb) {
    if (Progress::check_abort()) return NumericVector(annz);
    prog.increment();
    // Populate row_to_pos for pattern column cb
    for (int idx = Ap[cb]; idx < Ap[cb + 1]; ++idx)
      row_to_pos[Ai[idx]] = idx;

    // For each feature f where cb is nonzero
    for (int fi = item_offsets[cb]; fi < item_offsets[cb + 1]; ++fi) {
      const int f = item_feats[fi];
      const double vb = item_vals[fi];

      // For each other item ca also in this feature, accumulate min into (ca, cb)
      for (int idx = feat_offsets[f]; idx < feat_offsets[f + 1]; ++idx) {
        const int ca = feat_items[idx];
        const int pos = row_to_pos[ca];
        if (pos >= 0)
          out[pos] += std::min(feat_vals[idx], vb);
      }
    }

    // Clear row_to_pos
    for (int idx = Ap[cb]; idx < Ap[cb + 1]; ++idx)
      row_to_pos[Ai[idx]] = -1;
  }

  return out;
}

//' Dense weighted Jaccard similarity via C++
//'
//' Compute the full weighted Jaccard similarity matrix for a dgCMatrix using
//' an adaptive dense accumulation strategy. For small output matrices, uses a
//' feature-oriented loop; for large outputs, switches to a column-oriented
//' loop for better cache performance. Intended for internal use by
//' \code{coconat::jaccard_sim}. See also \code{\link{c_weighted_jaccard_sparse}}.
//'
//' @param x A dgCMatrix (sparse column-compressed matrix)
//' @param transpose If \code{FALSE}, compare columns; if \code{TRUE}, compare
//'   rows
//' @return A dense numeric similarity matrix
//' @export
// [[Rcpp::export]]
NumericMatrix c_weighted_jaccard_dense(const S4& x, bool transpose = false) {
  IntegerVector dims = x.slot("Dim");
  IntegerVector xi = x.slot("i");
  IntegerVector xp = x.slot("p");
  NumericVector xx = x.slot("x");

  const int nr = dims[0];
  const int nc = dims[1];
  const int ncomp = transpose ? nr : nc;
  const int nfeat = transpose ? nc : nr;

  // Totals per compared item
  std::vector<double> totals(ncomp, 0.0);

  // Build feature CSR: for each feature f, which items are nonzero
  std::vector<int> feat_count(nfeat, 0);
  if (!transpose) {
    for (int col = 0; col < nc; ++col) {
      for (int idx = xp[col]; idx < xp[col + 1]; ++idx) {
        totals[col] += xx[idx];
        feat_count[xi[idx]]++;
      }
    }
  } else {
    for (int col = 0; col < nc; ++col) {
      feat_count[col] = xp[col + 1] - xp[col];
      for (int idx = xp[col]; idx < xp[col + 1]; ++idx) {
        totals[xi[idx]] += xx[idx];
      }
    }
  }

  std::vector<int> feat_offsets(nfeat + 1, 0);
  for (int f = 0; f < nfeat; ++f)
    feat_offsets[f + 1] = feat_offsets[f] + feat_count[f];

  const int total_nnz = feat_offsets[nfeat];
  std::vector<int> feat_items(total_nnz);
  std::vector<double> feat_vals(total_nnz);
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

  NumericMatrix out(ncomp, ncomp);

  // Threshold: ncomp^2 * 8 bytes > ~12MB (typical L3 share per core)
  const bool use_colwise = (static_cast<std::size_t>(ncomp) * ncomp * 8 > 12 * 1024 * 1024);

  if (use_colwise) {
    // Column-oriented: build item->features CSR, iterate output columns
    std::vector<int> item_count(ncomp, 0);
    for (int f = 0; f < nfeat; ++f)
      for (int idx = feat_offsets[f]; idx < feat_offsets[f + 1]; ++idx)
        item_count[feat_items[idx]]++;

    std::vector<int> item_offsets(ncomp + 1, 0);
    for (int c = 0; c < ncomp; ++c)
      item_offsets[c + 1] = item_offsets[c] + item_count[c];

    std::vector<int> item_feats(total_nnz);
    std::vector<double> item_vals(total_nnz);
    std::vector<int> item_fill(item_offsets.begin(), item_offsets.end());

    for (int f = 0; f < nfeat; ++f) {
      for (int idx = feat_offsets[f]; idx < feat_offsets[f + 1]; ++idx) {
        const int item = feat_items[idx];
        const int dest = item_fill[item]++;
        item_feats[dest] = f;
        item_vals[dest] = feat_vals[idx];
      }
    }

    for (int cb = 0; cb < ncomp; ++cb) {
      double* col = &out[static_cast<std::size_t>(cb) * ncomp];
      for (int fi = item_offsets[cb]; fi < item_offsets[cb + 1]; ++fi) {
        const int f = item_feats[fi];
        const double vb = item_vals[fi];
        for (int idx = feat_offsets[f]; idx < feat_offsets[f + 1]; ++idx) {
          col[feat_items[idx]] += std::min(feat_vals[idx], vb);
        }
      }
    }
  } else {
    // Feature-oriented: iterate features, scatter into output
    for (int f = 0; f < nfeat; ++f) {
      const int start = feat_offsets[f];
      const int end = feat_offsets[f + 1];
      const int k = end - start;

      for (int a = 0; a < k; ++a) {
        const int ca = feat_items[start + a];
        const double va = feat_vals[start + a];
        out[static_cast<std::size_t>(ca) * ncomp + ca] += va;
        for (int b = a + 1; b < k; ++b) {
          const int cb = feat_items[start + b];
          const double mv = std::min(va, feat_vals[start + b]);
          out[static_cast<std::size_t>(ca) * ncomp + cb] += mv;
          out[static_cast<std::size_t>(cb) * ncomp + ca] += mv;
        }
      }
    }
  }

  // Convert min_sums to similarity
  for (int cb = 0; cb < ncomp; ++cb) {
    for (int ca = 0; ca < ncomp; ++ca) {
      if (ca == cb) {
        out[static_cast<std::size_t>(cb) * ncomp + ca] = 1.0;
      } else {
        const std::size_t pos = static_cast<std::size_t>(cb) * ncomp + ca;
        const double ms = out[pos];
        const double denom = totals[ca] + totals[cb] - ms;
        out[pos] = (denom > 0.0) ? ms / denom : 0.0;
      }
    }
  }

  return out;
}
