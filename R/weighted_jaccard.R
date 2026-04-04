#' Sparse weighted Jaccard similarity via C++
#'
#' Compute the weighted Jaccard similarity matrix for a dgCMatrix, returning a
#' sparse result. Uses \code{weighted_jaccard_sparse_fill} to compute
#' min-sums only for column (or row) pairs that share at least one non-zero
#' feature, then normalises to similarity.
#'
#' @param x A dgCMatrix (sparse column-compressed matrix)
#' @param transpose If \code{FALSE} (default), compare columns; if
#'   \code{TRUE}, compare rows.
#' @return A sparse dgCMatrix similarity matrix
#' @export
#' @seealso \code{\link{c_weighted_jaccard_dense}} for the dense equivalent
#' @examples
#' \dontrun{
#' library(Matrix)
#' m <- sparseMatrix(i = c(1,2,1,2,3,3), j = c(1,1,2,2,2,3),
#'                   x = c(4,2,1,3,3,1), dims = c(3,3))
#' c_weighted_jaccard_sparse(m)
#' }
c_weighted_jaccard_sparse <- function(x, transpose = FALSE) {
  crossfun <- if (transpose) Matrix::tcrossprod else Matrix::crossprod
  n <- if (transpose) nrow(x) else ncol(x)

  if (length(x@x) == 0L) {
    return(Matrix::sparseMatrix(i = seq_len(n), j = seq_len(n), x = 1,
                                dims = c(n, n)))
  }

  # Sparsity pattern from binarised crossprod
  b <- x
  b@x <- rep(1, length(b@x))
  A <- as(crossfun(b), "generalMatrix")

  # Column/row totals for denominator
  totals <- if (transpose) Matrix::rowSums(x) else Matrix::colSums(x)

  # Fill in min_sums using C++
  A@x <- weighted_jaccard_sparse_fill(x, A, transpose = transpose)

  # Convert min_sums to similarity: sim = ms / (s_i + s_j - ms)
  col_idx <- rep(seq_along(diff(A@p)), diff(A@p))
  row_idx <- A@i + 1L
  denom <- totals[row_idx] + totals[col_idx] - A@x
  nonzero <- denom != 0
  A@x[nonzero] <- A@x[nonzero] / denom[nonzero]
  A@x[!nonzero] <- 0
  Matrix::diag(A) <- 1

  A
}
