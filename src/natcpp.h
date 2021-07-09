#include <Rcpp.h>
using namespace Rcpp;

// Check that we really have an integer vector inside our seglist
// and cast if necessary
inline IntegerVector check_segvec(const SEXP &x) {
  switch(TYPEOF(x)) {
  case INTSXP:
    return IntegerVector(x);
    break;
  case REALSXP:
    return as<IntegerVector>(x);
    break;
  }
  stop("seglist must contain integer (or numeric) vectors!");
}
