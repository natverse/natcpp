# natcpp 0.3.0

* add weighted Jaccard similarity functions for sparse data
  - `c_weighted_jaccard_sparse()` returns a sparse similarity (or distance) matrix
  - `c_weighted_jaccard_dense()` returns a dense matrix or a `dist`-layout vector
  - both run in parallel via RcppThread, with an optional progress bar (RcppProgress)
  - `triangle=` switches between symmetric upper-triangle output and a full matrix
  - `distance=` returns `1 - similarity`
* portable `src/Makevars` using `$(SHLIB_OPENMP_CXXFLAGS)` in place of a GNU make
  `$(shell ...)` callout to `RcppThread::LdFlags()`

**Full Changelog**: https://github.com/natverse/natcpp/compare/v0.2...v0.3.0

# natcpp 0.2

* add fast coordinate / index conversion functions in https://github.com/natverse/natcpp/pull/2
  - `c_coords21dindex()` 
  - `c_ijkpos()`
  - `c_sub2ind()`
  to make `nat::coord2ind()` and friends much more memory efficient

**Full Changelog**: https://github.com/natverse/natcpp/compare/v0.1.1...v0.2

# natcpp 0.1.1

* adds `c_ListofMatrixRows()`

# natcpp 0.1.0

* First version of the package with basic functions for manipulating segment
  lists and calculating cable lengths.
