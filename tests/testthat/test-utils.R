test_that("natcpp_threads honours explicit values", {
  expect_identical(natcpp_threads(1L), 1L)
  expect_identical(natcpp_threads(8), 8L)
  # 0 (all cores) is passed through unchanged, not treated as "unset"
  expect_identical(natcpp_threads(0L), 0L)
})

test_that("natcpp_threads applies the default policy", {
  old_opt <- options(Ncpus = 3L)
  on.exit(options(old_opt), add = TRUE)
  expect_identical(natcpp_threads(NULL), 3L)

  options(Ncpus = NULL)
  old_env <- Sys.getenv("OMP_THREAD_LIMIT", unset = NA)
  on.exit(
    if (is.na(old_env)) Sys.unsetenv("OMP_THREAD_LIMIT")
    else Sys.setenv(OMP_THREAD_LIMIT = old_env),
    add = TRUE
  )

  Sys.setenv(OMP_THREAD_LIMIT = "5")
  expect_identical(natcpp_threads(NULL), 5L)

  Sys.unsetenv("OMP_THREAD_LIMIT")
  expect_identical(natcpp_threads(NULL), 2L)
})
