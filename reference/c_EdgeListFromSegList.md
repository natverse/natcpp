# Turn a segment list into an edgelist suitable for constructing an ngraph

Turn a segment list into an edgelist suitable for constructing an ngraph

## Usage

``` r
c_EdgeListFromSegList(L)
```

## Arguments

- L:

  a list containing integer vectors from `as.seglist`

## Value

An integer matrix of N rows and 2 columns

## Details

It is up to the caller to generate the `seglist`. Note that isolated
points will be dropped since they have no edges.

## Examples

``` r
if (FALSE) { # \dontrun{
library(nat)
# make a neuron with multiple subtrees
n=prune_vertices(Cell07PNs[[1]], 48L)
# Must use flatten=T if including all subtrees
sl=as.seglist(n, all = TRUE, flatten = TRUE)
c_EdgeListFromSegList(sl)
} # }
```
