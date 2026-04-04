# A simple function to compute the lengths of the elements of an R list

A simple function to compute the lengths of the elements of an R list

## Usage

``` r
c_listlengths(L)
```

## Arguments

- L:

  a list

## Value

An integer vector containing the length of each element of `L`

## Details

This is equivalent to the
[`base::lengths`](https://rdrr.io/r/base/lengths.html) however it it
much faster for long lists (and somewhat slower for short ones).
