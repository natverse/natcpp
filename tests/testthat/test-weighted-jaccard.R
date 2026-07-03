test_that("c_weighted_jaccard_dense computes correct similarity", {
  m <- Matrix::sparseMatrix(
    i = c(1L, 2L, 1L, 2L, 3L, 3L),
    j = c(1L, 1L, 2L, 2L, 2L, 3L),
    x = c(4, 2, 1, 3, 3, 1),
    dims = c(3L, 3L)
  )

  # Naive reference: sim[a,b] = sum(pmin(col_a, col_b)) / sum(pmax(col_a, col_b))
  dm <- as.matrix(m)
  n <- ncol(dm)
  ref <- matrix(0, n, n)
  for (a in seq_len(n)) {
    for (b in seq_len(n)) {
      mins <- sum(pmin(dm[, a], dm[, b]))
      maxs <- sum(pmax(dm[, a], dm[, b]))
      ref[a, b] <- if (maxs > 0) mins / maxs else 0
    }
  }
  diag(ref) <- 1

  sim <- c_weighted_jaccard_dense(m, transpose = FALSE)
  expect_equal(sim, ref, tolerance = 1e-12)
})

test_that("c_weighted_jaccard_dense transpose works", {
  m <- Matrix::sparseMatrix(
    i = c(1L, 2L, 1L, 2L, 3L, 3L),
    j = c(1L, 1L, 2L, 2L, 2L, 3L),
    x = c(4, 2, 1, 3, 3, 1),
    dims = c(3L, 3L)
  )

  dm <- as.matrix(m)
  nr <- nrow(dm)
  ref_t <- matrix(0, nr, nr)
  for (a in seq_len(nr)) {
    for (b in seq_len(nr)) {
      mins <- sum(pmin(dm[a, ], dm[b, ]))
      maxs <- sum(pmax(dm[a, ], dm[b, ]))
      ref_t[a, b] <- if (maxs > 0) mins / maxs else 0
    }
  }
  diag(ref_t) <- 1

  sim_t <- c_weighted_jaccard_dense(m, transpose = TRUE)
  expect_equal(sim_t, ref_t, tolerance = 1e-12)
})

test_that("c_weighted_jaccard_sparse matches dense", {
  m <- Matrix::sparseMatrix(
    i = c(1L, 2L, 1L, 2L, 3L, 3L),
    j = c(1L, 1L, 2L, 2L, 2L, 3L),
    x = c(4, 2, 1, 3, 3, 1),
    dims = c(3L, 3L)
  )

  sim_dense <- c_weighted_jaccard_dense(m, transpose = FALSE)
  sim_sparse <- c_weighted_jaccard_sparse(m, transpose = FALSE, display_progress = FALSE)
  expect_equal(as.matrix(sim_sparse), sim_dense, tolerance = 1e-12)

  sim_dense_t <- c_weighted_jaccard_dense(m, transpose = TRUE)
  sim_sparse_t <- c_weighted_jaccard_sparse(m, transpose = TRUE, display_progress = FALSE)
  expect_equal(as.matrix(sim_sparse_t), sim_dense_t, tolerance = 1e-12)
})

test_that("c_weighted_jaccard_dense triangle returns dist object", {
  m <- Matrix::sparseMatrix(
    i = c(1L, 2L, 1L, 2L, 3L, 3L),
    j = c(1L, 1L, 2L, 2L, 2L, 3L),
    x = c(4, 2, 1, 3, 3, 1),
    dims = c(3L, 3L)
  )

  full <- c_weighted_jaccard_dense(m)

  # triangle=TRUE returns a dist object
  tri <- c_weighted_jaccard_dense(m, triangle = TRUE)
  expect_s3_class(tri, "dist")
  expect_equal(attr(tri, "Size"), ncol(m))

  # triangle similarity matches as.dist of full similarity
  # (as.dist extracts lower triangle; dist stores 1-sim by convention,
  #  but here we just compare the raw values)
  ref_lower <- full[lower.tri(full)]
  expect_equal(as.numeric(tri), ref_lower, tolerance = 1e-12)

  # triangle + distance matches as.dist(1 - full)
  tri_d <- c_weighted_jaccard_dense(m, triangle = TRUE, distance = TRUE)
  expect_s3_class(tri_d, "dist")
  expect_equal(as.numeric(tri_d), 1 - ref_lower, tolerance = 1e-12)
})

test_that("c_weighted_jaccard_dense distance works for full matrix", {
  m <- Matrix::sparseMatrix(
    i = c(1L, 2L, 1L, 2L, 3L, 3L),
    j = c(1L, 1L, 2L, 2L, 2L, 3L),
    x = c(4, 2, 1, 3, 3, 1),
    dims = c(3L, 3L)
  )

  full <- c_weighted_jaccard_dense(m)
  full_d <- c_weighted_jaccard_dense(m, distance = TRUE)
  expect_equal(full_d, 1 - full, tolerance = 1e-12)
})

test_that("c_weighted_jaccard_sparse triangle and distance args", {
  m <- Matrix::sparseMatrix(
    i = c(1L, 2L, 1L, 2L, 3L, 3L),
    j = c(1L, 1L, 2L, 2L, 2L, 3L),
    x = c(4, 2, 1, 3, 3, 1),
    dims = c(3L, 3L)
  )

  # Default: dgCMatrix
  sp <- c_weighted_jaccard_sparse(m, display_progress = FALSE)
  expect_true(inherits(sp, "dgCMatrix"))

  # triangle=TRUE: dsCMatrix
  sp_tri <- c_weighted_jaccard_sparse(m, display_progress = FALSE, triangle = TRUE)
  expect_true(inherits(sp_tri, "dsCMatrix"))
  expect_equal(as.matrix(sp_tri), as.matrix(sp), tolerance = 1e-12)

  # distance=TRUE warns
  expect_warning(
    sp_d <- c_weighted_jaccard_sparse(m, display_progress = FALSE, distance = TRUE),
    "distance=TRUE"
  )
  # Stored entries are correct (1-sim); structural zeros remain 0 instead of 1
  # — this is the inherent limitation the warning describes
  full <- c_weighted_jaccard_dense(m)
  ref_d <- 1 - full
  ref_d[full == 0 & row(full) != col(full)] <- 0  # zero out unstored positions
  expect_equal(as.matrix(sp_d), ref_d, tolerance = 1e-12)
})
