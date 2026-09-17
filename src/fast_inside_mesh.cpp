// Fast generalised winding number via libigl's "Fast Winding Numbers for Soups
// and Clouds" (Barill et al. 2018). A bounding-volume hierarchy is built once
// over the mesh, then each query point is evaluated in O(log F). This is the
// accelerated back end intended for large point sets on large meshes; the
// brute-force c_mesh_winding_number() is the O(P*F) reference.
//
// We call libigl's per-point primitive (UT_SolidAngle::computeSolidAngle, the
// same call libigl's batch overload makes internally) directly inside an
// RcppThread::parallelFor. This keeps the core count controllable via the
// `threads` argument -- so CRAN's check-farm limit can be respected -- while
// distributing points across threads ourselves rather than via libigl's own
// std::thread pool (whose size is fixed by a process-global singleton and so
// cannot honour a per-call thread count). Routing each query through libigl's
// templated fast_winding_number(bvh, acc, p) wrapper instead is dramatically
// slower, so we avoid it.
//
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <thread>
#include <RcppThread.h>
#include "vendor/igl/fast_winding_number.h"

using namespace Rcpp;
typedef igl::FastWindingNumber::HDK_Sample::UT_Vector3T<float> UTVec3f;

// [[Rcpp::export]]
NumericVector c_fast_mesh_winding_number(NumericMatrix points,
                                         NumericMatrix vertices,
                                         IntegerMatrix faces,
                                         int threads = 4,
                                         double accuracy = 2.0) {
  const int np = points.nrow();
  const int nv = vertices.nrow();
  const int nf = faces.nrow();
  if (points.ncol() != 3)   stop("points must be an Nx3 matrix");
  if (vertices.ncol() != 3) stop("vertices must be an Nx3 matrix");
  if (faces.ncol() != 3)    stop("faces must be an Nx3 matrix");

  Eigen::MatrixXd V(nv, 3);
  for (int i = 0; i < nv; ++i) {
    V(i, 0) = vertices(i, 0); V(i, 1) = vertices(i, 1); V(i, 2) = vertices(i, 2);
  }
  Eigen::MatrixXi F(nf, 3);
  for (int i = 0; i < nf; ++i) {
    const int a = faces(i, 0) - 1, b = faces(i, 1) - 1, c = faces(i, 2) - 1;
    if (a < 0 || b < 0 || c < 0 || a >= nv || b >= nv || c >= nv)
      stop("faces contains a vertex index outside [1, nrow(vertices)]");
    F(i, 0) = a; F(i, 1) = b; F(i, 2) = c;
  }

  // Precompute the BVH once (Taylor expansion order 2, as in libigl examples).
  igl::FastWindingNumberBVH bvh;
  igl::fast_winding_number(V, F, 2, bvh);

  // Copy query points into a plain buffer so the worker threads touch no R data
  // structures and no Eigen expression templates.
  std::vector<float> Q(static_cast<size_t>(np) * 3);
  for (int i = 0; i < np; ++i) {
    Q[3 * i + 0] = static_cast<float>(points(i, 0));
    Q[3 * i + 1] = static_cast<float>(points(i, 1));
    Q[3 * i + 2] = static_cast<float>(points(i, 2));
  }

  NumericVector out(np);
  double* out_ptr = out.begin();           // NumericVector is not thread-safe
  const float acc = static_cast<float>(accuracy);
  const size_t nThreads = (threads > 0) ? static_cast<size_t>(threads)
                                        : std::thread::hardware_concurrency();

  // computeSolidAngle is a const, read-only query on the shared BVH, so
  // concurrent single-point evaluation is thread-safe.
  RcppThread::parallelFor(0, np, [&](int i) {
    UTVec3f p;
    p[0] = Q[3 * i + 0]; p[1] = Q[3 * i + 1]; p[2] = Q[3 * i + 2];
    out_ptr[i] = bvh.ut_solid_angle.computeSolidAngle(p, acc) / (4.0 * igl::PI);
  }, nThreads);

  return out;
}
