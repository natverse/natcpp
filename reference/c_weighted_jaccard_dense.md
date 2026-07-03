# Dense weighted Jaccard similarity via C++

Compute the full weighted Jaccard similarity matrix for a
[dgCMatrix](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html) using
an adaptive dense accumulation strategy. For small output matrices, uses
a feature-oriented loop; for large outputs, switches to a
column-oriented loop for better cache performance.

## Usage

``` r
c_weighted_jaccard_dense(
  x,
  transpose = FALSE,
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

  If `FALSE`, compare columns; if `TRUE`, compare rows

- threads:

  Number of threads (default 4). Set to 0 for all cores.

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
