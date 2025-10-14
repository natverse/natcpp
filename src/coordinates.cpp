// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

// Define a lightweight accessor template that provides (i,j) access

// Generic accessor for NumericMatrix (no bounds checking needed)
template <typename T>
struct PtAccessor {
  int ncol, nrow;
  const T& xyz;

  PtAccessor(const T& xyz_) : xyz(xyz_) {
    nrow = xyz_.nrow();
    ncol = xyz_.ncol();
  }
  inline double operator()(int i, int j) const {
    return xyz(i, j); // NumericMatrix handles row & column bounds
  }
};

// Specialization for DataFrame
template <>
struct PtAccessor<DataFrame> {
  int ncol, nrow;
  std::vector<NumericVector> cols;

  PtAccessor(const DataFrame& df) {
    ncol = df.size();
    nrow = df.nrows();
    cols.reserve(ncol);
    for (int j = 0; j < ncol; ++j)
      cols.push_back(as<NumericVector>(df[j])); // shallow copy of each column
  }

  inline double operator()(int i, int j) const {
    if (j < 0 || j >= ncol)
      Rcpp::stop("Column index out of range");
    return cols[j][i]; // row bounds handled by NumericVector
  }
};


// Core computation, templated on accessor type
template <typename XYZ>
IntegerMatrix ijkpos_core(const PtAccessor<XYZ>& xyz,
                          const NumericVector& origin,
                          const NumericVector& voxdims,
                          const IntegerVector& dims,
                          bool clamp) {
  int n = xyz.nrow;
  int p = xyz.ncol;

  if (p != 3 || origin.size() != 3 || voxdims.size() != 3 || dims.size() != 3)
    stop("All coordinate-related vectors/matrices must have length or ncol = 3");

  IntegerMatrix ijk(n, 3);
  bool warn = false;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < 3; ++j) {
      double val = (xyz(i, j) - origin[j]) / voxdims[j];
      int idx = static_cast<int>(std::round(val)) + 1;

      if (idx < 1 || idx > dims[j]) {
        warn = true;
        if (clamp) {
          if (idx < 1) idx = 1;
          else if (idx > dims[j]) idx = dims[j];
        }
      }

      ijk(i, j) = idx;
    }
  }

  if (warn)
    Rcpp::warning("pixel coordinates outside image data");

  return ijk;
}


//' Convert physical coordinates to pixel coordinates
//'
//' @param xyz Nx3 matrix of physical coordinates
//' @param dims Integer dimensions of the 3d image array
//' @param origin Numeric: 3d coordinates of the origin
//' @param voxdims Numeric: 3 numbers describing the voxel dimensions
//' @param clamp Logical: whether or not to clamp values within the pixel
//'   boundaries of the image.
//' @return Nx3 integer matrix of pixel coordinates
//' @export
// [[Rcpp::export]]
IntegerMatrix c_ijkpos(SEXP xyz,
                NumericVector origin,
                NumericVector voxdims,
                IntegerVector dims,
                bool clamp = false) {
  // Dispatch on type without copying
  if (is<NumericMatrix>(xyz)) {
    NumericMatrix mat = as<NumericMatrix>(xyz);
    return ijkpos_core(PtAccessor<NumericMatrix>(mat), origin, voxdims, dims, clamp);
  } else if (is<DataFrame>(xyz)) {
    DataFrame df = as<DataFrame>(xyz);
    return ijkpos_core(PtAccessor<DataFrame>(df), origin, voxdims, dims, clamp);
  } else {
    stop("xyz must be a numeric matrix or a data frame with 3 numeric columns");
  }
}


template <typename T>
NumericVector sub2ind_core(IntegerVector dims, const T& indices) {
  int n = indices.nrow();
  int d = indices.ncol();

  if (dims.size() != d)
    stop("indices must have the same number of columns as dimensions in dims");

  // Compute stride multipliers (cumprod(c(1, dims[-length(dims)])))
  std::vector<double> stride(d);
  stride[0] = 1.0;
  for (int i = 1; i < d; ++i)
    stride[i] = stride[i - 1] * dims[i - 1];

  NumericVector ndx(n);
  for (int r = 0; r < n; ++r) {
    double idx = 1.0; // R is 1-based
    for (int j = 0; j < d; ++j) {
      double v = indices(r, j);
      if (v < 1 || v > dims[j])
        Rcpp::warning("index out of range");
      idx += (v - 1.0) * stride[j];
    }
    ndx[r] = idx;
  }

  return ndx;
}

//' Find 1D index given n-dimensional indices
//' @param dims Integer dimensions of the array (usually 3d)
//' @param indices Nx3 integer matrix of pixel coordinates
//' @return numeric vector of linear indices into the array
//' @export
// [[Rcpp::export]]
NumericVector c_sub2ind(IntegerVector dims, IntegerMatrix indices) {
  return sub2ind_core(dims, indices);
}

NumericVector c_sub2ind(IntegerVector dims, NumericMatrix indices) {
  return sub2ind_core(dims, indices);
}

template <typename coords>
NumericVector coords21dindex_core(const PtAccessor<coords>& xyz,
                            const NumericVector& origin,
                            const NumericVector& voxdims,
                            const IntegerVector& dims,
                            bool clamp = false) {
  int n = xyz.nrow;
  int d = xyz.ncol;
  if (origin.size() != d || voxdims.size() != d || dims.size() != d)
    Rcpp::stop("origin, voxdims, and dims must have same length as number of columns in xyz");

  // Stride for linear indexing. Note R has column-major orientation
  std::vector<double> stride(d);
  stride[0] = 1.0;
  for (int j = 1; j < d; ++j)
    stride[j] = stride[j-1] * dims[j-1];

  NumericVector linidx(n);
  bool warn_flag = false;

  for (int i = 0; i < n; ++i) {
    double idx = 1.0; // 1-based R index
    for (int j = 0; j < d; ++j) {
      double pix = (xyz(i,j) - origin[j]) / voxdims[j] + 1.0;
      int pixidx = static_cast<int>(std::round(pix));

      if (pixidx < 1 || pixidx > dims[j]) {
        warn_flag = true;
        if (clamp) {
          if (pixidx < 1) pixidx = 1;
          else if (pixidx > dims[j]) pixidx = dims[j];
        }
      }

      idx += (pixidx - 1.0) * stride[j];
    }
    linidx[i] = idx;
  }

  if (warn_flag)
    Rcpp::warning("pixel coordinates outside image data");

  return linidx;
}

//' Convert physical coordinates to 1d indices into image array
//'
//' @param xyz Nx3 matrix or data.frame of physical coordinates
//' @param dims Integer dimensions of the 3d image array
//' @param origin Numeric: 3d coordinates of the origin
//' @param voxdims Numeric: 3 numbers describing the voxel dimensions
//' @param clamp Logical: whether or not to clamp values within the pixel
//'   boundaries of the image.
//' @return Nx3 integer matrix of pixel coordinates
//' @export
// [[Rcpp::export]]
NumericVector c_coords21dindex(SEXP xyz,
                  NumericVector origin,
                  NumericVector voxdims,
                  IntegerVector dims,
                  bool clamp = false) {
  if (Rcpp::is<NumericMatrix>(xyz)) {
    NumericMatrix mat = as<NumericMatrix>(xyz);
    return coords21dindex_core(PtAccessor<NumericMatrix>(mat), origin, voxdims, dims, clamp);
  } else if (Rcpp::is<DataFrame>(xyz)) {
    DataFrame df = as<DataFrame>(xyz);
    return coords21dindex_core(PtAccessor<DataFrame>(df), origin, voxdims, dims, clamp);
  } else {
    Rcpp::stop("xyz must be a numeric matrix or a data frame");
  }
}
