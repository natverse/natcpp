test_that("multiplication works", {
  m=matrix(1:6, ncol=3, byrow = T)
  l=list(1:3, 4:6)
  colnames(m)=c("X","Y","Z")
  expect_equal(c_ListofMatrixRows(m), l)

  nm=matrix(rnorm(n = 300), ncol=3)

  expect_equal(c_ListofMatrixRows(nm),
               lapply(seq_len(nrow(nm)), function(i) nm[i,]))

  expect_equal(c_ListofMatrixRows(matrix(1L, nrow=2, ncol=0)),
               list(integer(), integer()))
  expect_equal(c_ListofMatrixRows(matrix(1L, nrow=0, ncol=2)),
               list())
  expect_error(c_ListofMatrixRows(l))
  # matrix of lists
  ml=matrix(rep(list(1), 4), ncol = 2)
  expect_error(c_ListofMatrixRows(ml))
})
