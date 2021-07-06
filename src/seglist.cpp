#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector c_listlengths(const List &L) {
  IntegerVector lens(L.size());

  for (int i=0; i<L.size(); i++) {
    SEXP x = L[i];
    lens[i] = Rf_length(x);
  }
  return lens;
}

// [[Rcpp::export]]
IntegerMatrix c_topntail(const List &L) {
  IntegerMatrix m( 2 , L.length() );

  for (int i=0; i<L.size(); i++) {
    Rcpp::IntegerVector V = L[i];
    if(V.length()<1)
      continue;
    m(0,i)=V[0];
    m(1,i)=V[V.length()-1];
  }
  return m;
}

// [[Rcpp::export]]
IntegerMatrix c_EdgeListFromSegList(const List &L) {

  IntegerVector lsl=c_listlengths(L);
  // length of segments
  lsl=lsl[lsl>1];
  int totseglen=sum(lsl-1);

  IntegerMatrix el(totseglen, 2);
  // count elements so far
  int nelems=0;

  for (int i=0; i<L.size(); i++) {
    Rcpp::IntegerVector V = L[i];
    int nv=V.length();
    if(nv>1) {
      for(int j=0;j<(nv-1);j++) {
        el(nelems+j,0)=V[j];
        el(nelems+j,1)=V[j+1];
      }
      nelems=nelems+nv-1;
    }
  }
  return el;
}
