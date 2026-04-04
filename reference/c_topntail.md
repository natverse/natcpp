# Find the first and last elements of all vectors in a list

`c_topntail` returns an 2xN matrix containing the start and end of each
of the vectors in the input list. Length 0 vectors are ignored, while
length 1 vectors are duplicated

For `c_topntail_list`, a list of the same length as `L` having the same
elements when their length is \<=2 or the first and last elements when
length\>2.

## Usage

``` r
c_topntail(L)

c_topntail_list(L)
```

## Arguments

- L:

  a list containing integer vectors, typically a `seglist`

## Value

For `c_topntail` an integer `matrix`. For `c_topntail_list` a `list`.
