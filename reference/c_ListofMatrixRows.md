# Convert a matrix into list of row vectors

Convert a matrix into list of row vectors

## Usage

``` r
c_ListofMatrixRows(object)
```

## Arguments

- object:

  An integer, numeric, character or logical matrix of N rows and M
  columns

## Value

a list containing N vectors of length M corresponding to the rows of
`object`.

## Details

Typically this will be for 3D coordinates but there are no limits on row
length.

## Examples

``` r
if (FALSE) { # \dontrun{
library(nat)
xyz=xyzmatrix(Cell07PNs)
mat2list = function(m) {
um=unname(m)
lapply(1:nrow(um), function(i) um[i,])
}
bench::mark(rcpp=c_ListofMatrixRows(xyz), r=mat2list(xyz))
} # }
```
