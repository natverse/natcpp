test_that("c_EdgeListFromSegList works", {
  n=readRDS(test_path('testdata/cell071.rds'))
  # el=nat:::EdgeListFromSegList(n$SegList)
  # mode(el)='integer'
  # dimnames(el)=NULL
  # expect_snapshot_value(el, style = 'json2')
  expect_snapshot_value(c_EdgeListFromSegList(n$SegList), style = 'json2')
})


test_that("topntail works", {
  n=readRDS(test_path('testdata/cell071.rds'))
  # slt=sapply(n$SegList, function(x) as.integer(c(x[1], x[length(x)])))
  # expect_snapshot_value(slt, style = 'json2')
  expect_snapshot_value(tntm <- c_topntail(n$SegList), style = 'json2')
  expect_equal(c_topntail_list(n$SegList), apply(tntm, 2, c, simplify = F))
})
