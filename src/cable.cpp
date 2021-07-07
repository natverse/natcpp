#include <Rcpp.h>
using namespace Rcpp;

//' Compute summed segment lengths or total cable
//' @description \code{c_seglengths} comutes the summed segment length equivalent
//'  to \code{nat::seglengths(sumsegment = T)}
//' @param sl A \code{seglist} with 1-indices into vectors x,y,z
//' @param x,y,z Numeric vectors with 3D coordinate data
//' (which could be columns from a data frame)
//' @export
// [[Rcpp::export]]
NumericVector c_seglengths(const List &sl, const NumericVector &x,
                           const NumericVector &y, const NumericVector &z) {
  NumericVector lens(sl.size());
  for (int i=0; i<sl.size(); i++) {
    // nb must decrement by 1 for C style indices
    const Rcpp::IntegerVector idxs = sl[i];
    const int nv=idxs.length();

    double sd=0.0;
    for(int j=0; j<(nv-1); j++) {
      double dx=x[idxs[j+1]-1]-x[idxs[j]-1];
      double dy=y[idxs[j+1]-1]-y[idxs[j]-1];
      double dz=z[idxs[j+1]-1]-z[idxs[j]-1];
      double d=std::sqrt(dx*dx+dy*dy+dz*dz);
      if(!ISNAN(d)) {
        // this deals with NA at the level of each point
        sd += d;
      }
    }
    lens[i]=sd;
  }
  return lens;
}

//' @description \code{c_total_cable} computes the summed total cable for a
//' whole neuron. It's intended use is the \code{nat::summary.neuron} function.
//' @rdname c_seglengths
// [[Rcpp::export]]
double c_total_cable(const List &sl, const NumericVector &x,
                        const NumericVector &y, const NumericVector &z) {
  double total_cable=0.0;
  for (int i=0; i<sl.size(); i++) {
    // nb must decrement by 1 for C style indices
    const Rcpp::IntegerVector idxs = sl[i];
    const int nv=idxs.length();

    double sd=0.0;
    for(int j=0; j<(nv-1); j++) {
      double dx=x[idxs[j+1]-1]-x[idxs[j]-1];
      double dy=y[idxs[j+1]-1]-y[idxs[j]-1];
      double dz=z[idxs[j+1]-1]-z[idxs[j]-1];
      double d=std::sqrt(dx*dx+dy*dy+dz*dz);
      if(!ISNAN(d)) {
        // this deals with NA at the level of each point
        sd += d;
      }
    }
    total_cable += sd;
  }
  return total_cable;
}

// This is slower than c_total_cable and more importantly uses much more
// memory on the R side to compute the parent indices
// expects parent indices and columns of xyz data
// (which could be columns from a data frame)
// parent_idx is -1
double c_total_cablepi(const IntegerVector parent_idx,
                     const NumericVector &x,
                     const NumericVector &y,
                     const NumericVector &z) {
  double total_cable=0.0;

  for (int i=0; i<parent_idx.size(); i++) {
    // nb must decrement by 1 for C style indices
    if(parent_idx[i]<1) {
      continue;
    }
    double dx=x[parent_idx[i]-1]-x[i];
    double dy=y[parent_idx[i]-1]-y[i];
    double dz=z[parent_idx[i]-1]-z[i];
    double d = std::sqrt(dx*dx+dy*dy+dz*dz);
    if(!ISNAN(d)) {
      total_cable += d;
    }
  }
  return total_cable;
}
