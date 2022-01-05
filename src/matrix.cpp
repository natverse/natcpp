#include "natcpp.h"

template <typename T>
List c_ListofMatrixRowsT(const T &m) {
  List L(m.nrow());

  for(int i=0; i<(m.nrow()); i++) {
    // NumericVector v = m( 0 , _ );
    L[i]=m.row(i);
  }
  return L;
}


//' Convert a matrix into list of row vectors
//'
//' @details Typically this will be for 3D coordinates but there are no limits
//'   on row length.
//' @param object An integer, numeric, character or logical matrix of N rows and
//'   M columns
//' @return a list containing N vectors of length M corresponding to the rows of
//'   \code{object}.
//' @export
//' @examples
//' \dontrun{
//' library(nat)
//' xyz=xyzmatrix(Cell07PNs)
//' mat2list = function(m) {
//' um=unname(m)
//' lapply(1:nrow(um), function(i) um[i,])
//' }
//' bench::mark(rcpp=c_ListofMatrixRows(xyz), r=mat2list(xyz))
//' }
//' @export
// [[Rcpp::export]]
List c_ListofMatrixRows(const SEXP &object) {
  if (!Rf_isMatrix(object))
    stop("Please provide a matrix!");

  switch (TYPEOF(object))
  {
  case INTSXP: return c_ListofMatrixRowsT<IntegerMatrix>(object);
  case REALSXP: return c_ListofMatrixRowsT<NumericMatrix>(object);
  case CHARSXP: return c_ListofMatrixRowsT<CharacterMatrix>(object);
  case LGLSXP: return c_ListofMatrixRowsT<LogicalMatrix>(object);
  }
  stop("Unimplemented matrix type!");
}
