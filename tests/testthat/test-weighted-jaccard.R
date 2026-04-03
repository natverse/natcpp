test_that("weighted_jaccard_dense_cpp computes correct similarity", {
  skip_if_not_installed("Matrix")
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

  sim <- weighted_jaccard_dense_cpp(m, transpose = FALSE)
  expect_equal(sim, ref, tolerance = 1e-12)
})

test_that("weighted_jaccard_dense_cpp transpose works", {
  skip_if_not_installed("Matrix")
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

  sim_t <- weighted_jaccard_dense_cpp(m, transpose = TRUE)
  expect_equal(sim_t, ref_t, tolerance = 1e-12)
})

test_that("weighted_jaccard_sparse_fill_cpp matches dense", {
  skip_if_not_installed("Matrix")
  m <- Matrix::sparseMatrix(
    i = c(1L, 2L, 1L, 2L, 3L, 3L),
    j = c(1L, 1L, 2L, 2L, 2L, 3L),
    x = c(4, 2, 1, 3, 3, 1),
    dims = c(3L, 3L)
  )

  # Dense reference
  sim_dense <- weighted_jaccard_dense_cpp(m, transpose = FALSE)

  # Sparse fill: build pattern from binarised crossprod
  b <- m
  b@x <- rep(1, length(b@x))
  A <- methods::as(Matrix::crossprod(b), "generalMatrix")
  totals <- Matrix::colSums(m)

  A@x <- weighted_jaccard_sparse_fill_cpp(m, A, transpose = FALSE)

  # Convert min_sums to similarity
  col_idx <- rep(seq_along(diff(A@p)), diff(A@p))
  row_idx <- A@i + 1L
  denom <- totals[row_idx] + totals[col_idx] - A@x
  nonzero <- denom != 0
  A@x[nonzero] <- A@x[nonzero] / denom[nonzero]
  A@x[!nonzero] <- 0
  Matrix::diag(A) <- 1

  expect_equal(as.matrix(A), sim_dense, tolerance = 1e-12)
})
