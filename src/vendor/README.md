# Vendored third-party source

## libigl (`igl/`)

A minimal subset of [libigl](https://libigl.github.io/) providing
`igl::fast_winding_number` (the "Fast Winding Numbers for Soups and Clouds"
method of Barill et al. 2018), used as the accelerated back end for the
point-in-mesh test.

- Upstream: https://github.com/libigl/libigl
- Commit: `7100764` (branch `main`)
- Licence: MPL-2.0 (see `igl/LICENSE.MPL2`), which is compatible with this
  package's GPL (>= 3).

Only the transitive include closure of `igl/fast_winding_number.h` is vendored
(10 files). The bulk is `FastWindingNumberForSoups.h`, the self-contained HDK
solid-angle BVH implementation; it depends only on Eigen (provided here via
`RcppEigen`). No CGAL, Boost, GMP, MPFR or TBB is required or compiled (the
`tbb/*` includes in `parallel_for.h` are behind `IGL_PARALLEL_FOR_TBB`, which
is not defined).

To update: shallow-clone libigl, copy the closure of
`igl/fast_winding_number.h` (see `tools/`/commit history for the list) plus
`LICENSE.MPL2`, and record the new commit here.
