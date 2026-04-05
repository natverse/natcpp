#' Sparse weighted Jaccard similarity via C++
#'
#' Compute the weighted Jaccard similarity matrix for a dgCMatrix, returning a
#' symmetric sparse result. Uses \code{weighted_jaccard_sparse_fill} to compute
#' min-sums only for column (or row) pairs that share at least one non-zero
#' feature, then normalises to similarity. Only the upper triangle is computed,
#' taking advantage of the symmetry of the Jaccard index.
#'
#' @param x A dgCMatrix (sparse column-compressed matrix)
#' @param transpose If \code{FALSE} (default), compare columns; if
#'   \code{TRUE}, compare rows.
#' @param display_progress Whether to show a text progress bar (default
#'   \code{TRUE}).
#' @param threads Number of threads for parallel computation (default 4).
#'   Set to 0 to use all available cores.
#' @return A symmetric sparse dsCMatrix similarity matrix
#' @export
#' @seealso \code{\link{c_weighted_jaccard_dense}} for the dense equivalent
#' @examples
#' \dontrun{
#' library(Matrix)
#' m <- sparseMatrix(i = c(1,2,1,2,3,3), j = c(1,1,2,2,2,3),
#'                   x = c(4,2,1,3,3,1), dims = c(3,3))
#' c_weighted_jaccard_sparse(m)
#' }
c_weighted_jaccard_sparse <- function(x, transpose = FALSE, display_progress = TRUE,
                                      threads = 4L) {
  crossfun <- if (transpose) Matrix::tcrossprod else Matrix::crossprod
  n <- if (transpose) nrow(x) else ncol(x)

  if (length(x@x) == 0L) {
    return(Matrix::sparseMatrix(i = seq_len(n), j = seq_len(n), x = 1,
                                dims = c(n, n), symmetric = TRUE))
  }

  # Sparsity pattern from binarised crossprod — dsCMatrix (upper triangle)
  b <- x
  b@x <- rep(1, length(b@x))
  A <- crossfun(b)

  # Column/row totals for denominator
  totals <- if (transpose) Matrix::rowSums(x) else Matrix::colSums(x)

  # Fill in min_sums using C++ (only upper triangle entries)
  A@x <- weighted_jaccard_sparse_fill(x, A, transpose = transpose,
                                      display_progress = display_progress,
                                      threads = threads)

  # Convert min_sums to similarity: sim = ms / (s_i + s_j - ms)
  col_idx <- rep(seq_along(diff(A@p)), diff(A@p))
  row_idx <- A@i + 1L
  denom <- totals[row_idx] + totals[col_idx] - A@x
  nonzero <- denom != 0
  A@x[nonzero] <- A@x[nonzero] / denom[nonzero]
  A@x[!nonzero] <- 0

  # Set diagonal to 1 — diagonal entries are in the upper triangle pattern
  diag_pos <- which(A@i + 1L == col_idx)
  A@x[diag_pos] <- 1

  A
}
