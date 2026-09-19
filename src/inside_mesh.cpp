#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <thread>
#include <RcppThread.h>

using namespace Rcpp;

// Generalised (solid-angle) winding number of a triangle mesh at a set of query
// points. For a closed, consistently oriented mesh the winding number is ~ +/-1
// for interior points and ~0 for exterior points, regardless of orientation, so
// thresholding |w| > 0.5 gives a robust inside/outside test that -- unlike a
// closest-point signed distance -- has no dependence on local surface normals
// and no ray-casting tie-breaking. Signed solid angle per triangle uses the
// Van Oosterom & Strackee (1983) formula.

// [[Rcpp::export]]
NumericVector c_mesh_winding_number(NumericMatrix points,
                                    NumericMatrix vertices,
                                    IntegerMatrix faces,
                                    int threads = 4) {
  const int np = points.nrow();
  const int nv = vertices.nrow();
  const int nf = faces.nrow();
  if (points.ncol() != 3)   stop("points must be an Nx3 matrix");
  if (vertices.ncol() != 3) stop("vertices must be an Nx3 matrix");
  if (faces.ncol() != 3)    stop("faces must be an Nx3 matrix");

  // Copy into contiguous std::vectors for thread-safe, cache-friendly reads.
  std::vector<double> vx(nv), vy(nv), vz(nv);
  for (int i = 0; i < nv; ++i) {
    vx[i] = vertices(i, 0); vy[i] = vertices(i, 1); vz[i] = vertices(i, 2);
  }
  std::vector<int> fa(nf), fb(nf), fc(nf);
  for (int i = 0; i < nf; ++i) {
    const int a = faces(i, 0) - 1, b = faces(i, 1) - 1, c = faces(i, 2) - 1; // 1-based -> 0-based
    if (a < 0 || b < 0 || c < 0 || a >= nv || b >= nv || c >= nv)
      stop("faces contains a vertex index outside [1, nrow(vertices)]");
    fa[i] = a; fb[i] = b; fc[i] = c;
  }
  std::vector<double> px(np), py(np), pz(np);
  for (int i = 0; i < np; ++i) {
    px[i] = points(i, 0); py[i] = points(i, 1); pz[i] = points(i, 2);
  }

  NumericVector out(np);
  double* out_ptr = out.begin();          // NumericVector is not thread-safe; write via raw pointer
  const double inv4pi = 1.0 / (4.0 * M_PI);

  const size_t nThreads = (threads > 0) ? static_cast<size_t>(threads)
                                        : std::thread::hardware_concurrency();

  RcppThread::parallelFor(0, np, [&](int i) {
    const double qx = px[i], qy = py[i], qz = pz[i];
    double omega = 0.0;
    for (int f = 0; f < nf; ++f) {
      // triangle vertices relative to the query point
      const double ax = vx[fa[f]] - qx, ay = vy[fa[f]] - qy, az = vz[fa[f]] - qz;
      const double bx = vx[fb[f]] - qx, by = vy[fb[f]] - qy, bz = vz[fb[f]] - qz;
      const double cx = vx[fc[f]] - qx, cy = vy[fc[f]] - qy, cz = vz[fc[f]] - qz;
      const double la = std::sqrt(ax*ax + ay*ay + az*az);
      const double lb = std::sqrt(bx*bx + by*by + bz*bz);
      const double lc = std::sqrt(cx*cx + cy*cy + cz*cz);
      // numerator: scalar triple product a . (b x c)
      const double num = ax*(by*cz - bz*cy) - ay*(bx*cz - bz*cx) + az*(bx*cy - by*cx);
      // denominator
      const double den = la*lb*lc
        + (ax*bx + ay*by + az*bz) * lc
        + (bx*cx + by*cy + bz*cz) * la
        + (cx*ax + cy*ay + cz*az) * lb;
      omega += 2.0 * std::atan2(num, den);   // signed solid angle of this triangle
    }
    out_ptr[i] = omega * inv4pi;
  }, nThreads);

  return out;
}
