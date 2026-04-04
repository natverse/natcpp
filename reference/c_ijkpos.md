# Convert physical coordinates to pixel coordinates

Convert physical coordinates to pixel coordinates

## Usage

``` r
c_ijkpos(xyz, origin, voxdims, dims, clamp = FALSE)
```

## Arguments

- xyz:

  Nx3 matrix of physical coordinates

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
