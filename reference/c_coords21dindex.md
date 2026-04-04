# Convert physical coordinates to 1d indices into image array

Convert physical coordinates to 1d indices into image array

## Usage

``` r
c_coords21dindex(xyz, origin, voxdims, dims, clamp = FALSE)
```

## Arguments

- xyz:

  Nx3 matrix or data.frame of physical coordinates

- origin:

  Numeric: 3d coordinates of the origin

- voxdims:

  Numeric: 3 numbers describing the voxel dimensions

- dims:

  Integer dimensions of the 3d image array

- clamp:

  Logical: whether or not to clamp values within the pixel boundaries of
  the image.

## Value

Nx3 integer matrix of pixel coordinates
