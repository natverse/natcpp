#include <Rcpp.h>
using namespace Rcpp;

//' A simple function to compute the lengths of the elements of an R list
//'
//' @details This is equivalent to the \code{base::lengths} however it it much
//' faster for long lists (and somewhat slower for short ones).
//' @param L a list
//' @return An integer vector containing the length of each element of \code{L}
//' @export
// [[Rcpp::export]]
IntegerVector c_listlengths(const List &L) {
  IntegerVector lens(L.size());

  for (int i=0; i<L.size(); i++) {
    SEXP x = L[i];
    lens[i] = Rf_length(x);
  }
  return lens;
}

//' Find the first and last element of a list
//'
//' @param L a list containing integer vectors
//' @return An integer vector containing the length of each element of \code{L}
//' @export
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

//' Turn a segment list into an edgelist suitable for constructing an ngraph
//' @details It is up to the caller to generate the seglist.
//' Note that isolated points will be dropped since they have no edges.
//' @param L a list containing integer vectors from \code{as.seglist}
//' @return An integer matrix of
//' @export
//' @examples
//'
//' \dontrun{
//' library(nat)
//' # make a neuron with multiple subtrees
//' n=prune_vertices(Cell07PNs[[1]], 48L)
//' # Must use flatten=T if including all subtrees
//' sl=as.seglist(n, all = TRUE, flatten = TRUE)
//' c_EdgeListFromSegList(sl)
//' }
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
