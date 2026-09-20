# Test which points lie inside a triangle mesh (generalised winding number)

Robust point-in-mesh test based on the generalised (solid-angle) winding
number. For a closed mesh the winding number is approximately \\\pm 1\\
for interior points and \\0\\ for exterior points, so `abs(w) > 0.5`
classifies points as inside. Unlike a closest-point signed-distance test
it does not depend on surface normals and has no ray-casting
tie-breaking, so it does not produce the spurious "outside point
classified as inside" results that normal-based tests can give near thin
protrusions or sharp features.

## Usage

``` r
c_pointsinside(
  points,
  vertices,
  faces,
  method = c("auto", "bvh", "bruteforce"),
  threads = NULL,
  accuracy = 2
)
```

## Arguments

- points:

  An Nx3 matrix of query point coordinates (or anything coercible with
  `as.matrix`).

- vertices:

  An Nx3 matrix of mesh vertex coordinates.

- faces:

  An Nx3 integer matrix of 1-based vertex indices (one triangle per
  row), e.g. `t(mesh$it)` for an rgl `mesh3d`.

- method:

  Winding-number back end: `"auto"` (default), `"bvh"` or
  `"bruteforce"`. See **Details**.

- threads:

  Number of threads for parallel computation. The default `NULL` applies
  the package thread policy (respecting `getOption("Ncpus")` and the
  `OMP_THREAD_LIMIT` environment variable, else 2). Set to 0 to use all
  available cores.

- accuracy:

  libigl accuracy-scale parameter for `method = "bvh"` (default 2);
  ignored by the brute-force back end.

## Value

A logical vector of length `nrow(points)` (`TRUE` = inside).

## Details

The mesh should be closed (watertight) and triangular; the result is
independent of face orientation (winding). It is intended as the
accelerated back end for `nat::pointsinside()`.

Two back ends are available, selected by `method`:

- `"bruteforce"`:

  A self-contained \\O(P \times F)\\ implementation (P points, F faces),
  parallelised over points with RcppThread. No setup cost, so it is
  fastest for small meshes.

- `"bvh"`:

  libigl's "Fast Winding Numbers for Soups and Clouds" (Barill et al.
  2018): a bounding-volume hierarchy is built once over the mesh and
  each query point is then evaluated in \\O(\log F)\\, so it scales to
  millions of points on meshes of tens of thousands of faces. `accuracy`
  tunes the multipole approximation.

`"auto"` (the default) picks `"bvh"` only for large meshes queried by
enough points to amortise building the hierarchy, and `"bruteforce"`
otherwise (so a simple mesh such as a cuboid always uses brute force).
The two back ends agree to within the winding-number tolerance.

## Examples

``` r
# tetrahedron
V <- rbind(c(0,0,0), c(1,0,0), c(0,1,0), c(0,0,1))
F <- rbind(c(1,3,2), c(1,2,4), c(1,4,3), c(2,3,4))
c_pointsinside(rbind(c(.2,.2,.2), c(2,2,2)), V, F)  # TRUE FALSE
#> [1]  TRUE FALSE
```
