# Sparse weighted Jaccard similarity via C++

Compute the weighted Jaccard similarity matrix for a
[dgCMatrix](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html),
returning a sparse result. Uses `weighted_jaccard_sparse_fill` to
compute min-sums only for column (or row) pairs that share at least one
non-zero feature, then normalises to similarity. Only the upper triangle
is computed, taking advantage of the symmetry of the Jaccard index.

## Usage

``` r
c_weighted_jaccard_sparse(
  x,
  transpose = FALSE,
  display_progress = TRUE,
  threads = 4L,
  triangle = FALSE,
  distance = FALSE
)
```

## Arguments

- x:

  A [dgCMatrix](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html)
  (sparse column-compressed matrix)

- transpose:

  If `FALSE` (default), compare columns; if `TRUE`, compare rows.

- display_progress:

  Whether to show a text progress bar (default `TRUE`).

- threads:

  Number of threads for parallel computation (default 4). Set to 0 to
  use all available cores.

- triangle:

  If `TRUE`, return a symmetric `dsCMatrix` (upper triangle only). If
  `FALSE` (default), return a general
  [dgCMatrix](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html).

- distance:

  If `TRUE`, return distance (`1 - similarity`) instead of similarity.
  Default `FALSE`. A warning is issued since sparse distance matrices
  are typically dense.

## Value

A sparse similarity (or distance) matrix: `dsCMatrix` when
`triangle = TRUE`,
[dgCMatrix](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html)
otherwise.

## See also

[`c_weighted_jaccard_dense`](https://natverse.org/natcpp/reference/c_weighted_jaccard_dense.md)
for the dense equivalent

## Examples

``` r
if (FALSE) { # \dontrun{
library(Matrix)
m <- sparseMatrix(i = c(1,2,1,2,3,3), j = c(1,1,2,2,2,3),
                  x = c(4,2,1,3,3,1), dims = c(3,3))
c_weighted_jaccard_sparse(m)
} # }
```
