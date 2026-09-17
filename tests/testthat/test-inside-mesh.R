# Analytic tetrahedron: interior winding number ~ +/-1, exterior ~ 0.
tetra <- function() {
  V <- rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
  F <- rbind(c(1, 3, 2), c(1, 2, 4), c(1, 4, 3), c(2, 3, 4))
  list(V = V, F = F)
}

test_that("winding number classifies an analytic tetrahedron", {
  m <- tetra()
  inside  <- rbind(c(0.2, 0.2, 0.2), c(0.1, 0.1, 0.1))
  outside <- rbind(c(2, 2, 2), c(-1, 0, 0), c(0.6, 0.6, 0.6))

  wi <- c_mesh_winding_number(inside,  m$V, m$F, threads = 2L)
  wo <- c_mesh_winding_number(outside, m$V, m$F, threads = 2L)
  expect_equal(abs(wi), c(1, 1), tolerance = 1e-6)
  expect_equal(wo, c(0, 0, 0), tolerance = 1e-6)

  expect_equal(c_pointsinside(inside,  m$V, m$F, threads = 2L), c(TRUE, TRUE))
  expect_equal(c_pointsinside(outside, m$V, m$F, threads = 2L),
               c(FALSE, FALSE, FALSE))
})

test_that("classification is independent of face orientation", {
  m <- tetra()
  Frev <- m$F[, c(1, 3, 2)]                 # uniform winding flip
  p <- rbind(c(0.2, 0.2, 0.2), c(2, 2, 2))
  w  <- c_mesh_winding_number(p, m$V, m$F,   threads = 2L)
  wr <- c_mesh_winding_number(p, m$V, Frev, threads = 2L)
  expect_equal(wr, -w, tolerance = 1e-9)    # sign flips, magnitude does not
  expect_equal(c_pointsinside(p, m$V, Frev, threads = 2L),
               c_pointsinside(p, m$V, m$F,  threads = 2L))
})

test_that("threads argument does not change the result", {
  m <- tetra()
  set.seed(1)
  p <- matrix(runif(300, -0.5, 1.5), ncol = 3)
  # keep <= 2 cores to respect CRAN's check-farm limit
  w1 <- c_mesh_winding_number(p, m$V, m$F, threads = 1L)
  w2 <- c_mesh_winding_number(p, m$V, m$F, threads = 2L)
  expect_identical(w1, w2)
})

test_that("input validation errors on malformed matrices", {
  m <- tetra()
  expect_error(c_mesh_winding_number(matrix(0, 2, 2), m$V, m$F), "Nx3")
  expect_error(c_mesh_winding_number(rbind(c(0, 0, 0)), matrix(0, 4, 2), m$F),
               "Nx3")
  bad <- rbind(c(1, 3, 2), c(1, 2, 99))     # vertex index out of range
  expect_error(c_mesh_winding_number(rbind(c(0, 0, 0)), m$V, bad),
               "vertex index")
})

test_that("real CA1 mesh: false positives outside, matches reference", {
  f <- test_path("testdata", "ca1_mesh.rds")
  skip_if_not(file.exists(f))
  d <- readRDS(f)

  # the four points a normal-based test wrongly called inside are all outside
  expect_false(any(c_pointsinside(d$false_positives, d$vertices, d$faces,
                                  threads = 2L)))
  expect_equal(c_mesh_winding_number(d$false_positives, d$vertices, d$faces,
                                     threads = 2L),
               rep(0, nrow(d$false_positives)), tolerance = 1e-3)

  # bbox-sampled points match an independent oracle (CGAL Side_of_triangle_mesh)
  set.seed(d$seed)
  bb <- apply(d$vertices, 2, range)
  P <- cbind(runif(d$n, bb[1, 1], bb[2, 1]),
             runif(d$n, bb[1, 2], bb[2, 2]),
             runif(d$n, bb[1, 3], bb[2, 3]))
  stopifnot(identical(dim(P), dim(d$points)))   # reproducible sample
  expect_equal(c_pointsinside(P, d$vertices, d$faces, threads = 2L),
               d$inside_ref)
})

test_that("fast (libigl) winding number agrees with brute force and oracle", {
  # analytic tetrahedron
  m <- tetra()
  p <- rbind(c(0.2, 0.2, 0.2), c(0.1, 0.1, 0.1), c(2, 2, 2), c(0.6, 0.6, 0.6))
  expect_equal(c_fast_pointsinside(p, m$V, m$F, threads = 2L),
               c_pointsinside(p, m$V, m$F, threads = 2L))

  f <- test_path("testdata", "ca1_mesh.rds")
  skip_if_not(file.exists(f))
  d <- readRDS(f)

  # the known false positives are outside by the fast method too
  expect_false(any(c_fast_pointsinside(d$false_positives, d$vertices, d$faces,
                                       threads = 2L)))

  # fast method matches the independent CGAL oracle on the bbox sample
  set.seed(d$seed)
  bb <- apply(d$vertices, 2, range)
  P <- cbind(runif(d$n, bb[1, 1], bb[2, 1]),
             runif(d$n, bb[1, 2], bb[2, 2]),
             runif(d$n, bb[1, 3], bb[2, 3]))
  expect_equal(c_fast_pointsinside(P, d$vertices, d$faces, threads = 2L),
               d$inside_ref)
})
