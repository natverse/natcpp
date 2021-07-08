test_that("c_EdgeListFromSegList works", {
  n=readRDS(test_path('testdata/cell071.rds'))
  # el=nat:::EdgeListFromSegList(n$SegList)
  # mode(el)='integer'
  # dimnames(el)=NULL
  # expect_snapshot_value(el, style = 'json2')
  expect_snapshot_value(c_EdgeListFromSegList(n$SegList), style = 'json2')
})
