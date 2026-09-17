#' Fast point-in-mesh test via libigl fast winding number
#'
#' Accelerated back end for point-in-mesh classification, using libigl's "Fast
#' Winding Numbers for Soups and Clouds" (Barill et al. 2018): a bounding-volume
#' hierarchy is built once over the mesh and each query point is evaluated in
#' \eqn{O(\log F)}, so it scales to millions of points on meshes of tens of
#' thousands of faces. Results match the brute-force [c_mesh_winding_number()]
#' (which is retained as an \eqn{O(P \times F)} correctness reference).
#'
#' @inheritParams c_pointsinside
#' @param accuracy libigl accuracy-scale parameter (default 2).
#' @return For `c_fast_pointsinside`, a logical vector (`TRUE` = inside). For
#'   `c_fast_mesh_winding_number`, the approximate winding number per point.
#' @rdname c_fast_pointsinside
#' @keywords internal
#' @noRd
c_fast_pointsinside <- function(points, vertices, faces, threads = 4L,
                                accuracy = 2) {
  w <- c_fast_mesh_winding_number(as.matrix(points), as.matrix(vertices),
                                  matrix(as.integer(faces), ncol = 3L),
                                  threads = threads, accuracy = accuracy)
  abs(w) > 0.5
}
