test_that("multiplication works", {
  m=matrix(1:6, ncol=3, byrow = T)
  l=list(1:3, 4:6)
  colnames(m)=c("X","Y","Z")
  expect_equal(c_ListofMatrixRows(m), l)

  nm=matrix(rnorm(n = 300), ncol=3)

  expect_equal(c_ListofMatrixRows(nm),
               lapply(seq_len(nrow(nm)), function(i) nm[i,]))
})
