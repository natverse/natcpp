# Compute summed segment lengths or total cable

`c_seglengths` computes the summed segment length equivalent to
`nat::seglengths(sumsegment = T)`

`c_total_cable` computes the summed total cable for a whole neuron. It's
intended use is the `nat::summary.neuron` function.

## Usage

``` r
c_seglengths(sl, x, y, z)

c_total_cable(sl, x, y, z)
```

## Arguments

- sl:

  A `seglist` with 1-indices into vectors x,y,z

- x, y, z:

  Numeric vectors with 3D coordinate data (which could be columns from a
  data frame)
