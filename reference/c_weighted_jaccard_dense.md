# Dense weighted Jaccard similarity via C++

Compute the full weighted Jaccard similarity matrix for a
[dgCMatrix](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html),
returning a dense matrix (or a
[`dist`](https://rdrr.io/r/stats/dist.html)-layout vector).

## Usage

``` r
c_weighted_jaccard_dense(
  x,
  transpose = FALSE,
  threads = NULL,
  triangle = FALSE,
  distance = FALSE
)
```

## Arguments

- x:

  A [dgCMatrix](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html)
  (sparse column-compressed matrix)

- transpose:

  If `FALSE`, compare columns; if `TRUE`, compare rows

- threads:

  Number of threads for parallel computation. The default `NULL` applies
  the package thread policy (respecting `getOption("Ncpus")` and the
  `OMP_THREAD_LIMIT` environment variable, else 2). Set to 0 to use all
  available cores.

- triangle:

  If `TRUE`, return only the lower triangle as a flat numeric vector in
  [`dist`](https://rdrr.io/r/stats/dist.html) layout. If `FALSE`
  (default), return a full square matrix.

- distance:

  If `TRUE`, return distance (`1 - similarity`) instead of similarity.
  Default `FALSE`.

## Value

A dense numeric similarity matrix, or a numeric vector in
[`dist`](https://rdrr.io/r/stats/dist.html) layout when
`triangle = TRUE`.

## Details

Uses an adaptive dense accumulation strategy: for small output matrices
a feature-oriented loop, switching to a column-oriented loop for larger
outputs for better cache performance.

## See also

[`c_weighted_jaccard_sparse`](https://natverse.org/natcpp/reference/c_weighted_jaccard_sparse.md)
for the sparse equivalent

## Examples

``` r
if (FALSE) { # \dontrun{
library(Matrix)
m <- sparseMatrix(i = c(1,2,1,2,3,3), j = c(1,1,2,2,2,3),
                  x = c(4,2,1,3,3,1), dims = c(3,3))
c_weighted_jaccard_dense(m)
} # }
```
