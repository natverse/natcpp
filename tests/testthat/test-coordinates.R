test_that("ijkpos, c_sub2ind, c_coords21dindex", {
  origin  <- c(10, 20, 30)
  voxdims <- c(1, 2, 3)
  dims    <- c(20, 30, 40)

  xyz_mat <- matrix(c(10, 20, 30,
                      11, 22, 33,
                      40, 80, 160), ncol = 3, byrow = TRUE)
  xyz_df <- data.frame(x = c(10, 11),
                       y = c(20, 22),
                       z = c(30, 33))
  bl=structure(c(1L, 2L, 20L, 1L, 2L, 30L, 1L, 2L, 40L), dim = c(3L, 3L))

  expect_warning(expect_equal(
    c_ijkpos(xyz_mat, origin, voxdims, dims, clamp = TRUE),
    bl))

  expect_silent(ijk <- c_ijkpos(xyz_mat[1:2,], origin, voxdims, dims, clamp = TRUE))
  expect_equal(c_ijkpos(xyz_mat[1:2,], origin, voxdims, dims, clamp = TRUE),
               c_ijkpos(xyz_df, origin, voxdims, dims, clamp = TRUE))

  bl=c(1,622)
  expect_equal(c_sub2ind(dims, indices = ijk),bl)
  expect_equal(
    c_sub2ind(dims, indices = ijk),
    c_coords21dindex(xyz_mat[1:2,], origin = origin, voxdims = voxdims, dims = dims))

})
