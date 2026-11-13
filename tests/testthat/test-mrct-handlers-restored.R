# mrct_region_sims() puts the caller's progressr.handlers option back on exit
# (TASK_cran_prep_policy_fixes_2026-11-12, fix 2). The call is mrct_rs_call()
# (helper-mrct-region-sims-record.R) on the S2 build (helper-mrct-analytic.R),
# 2 workers. The workers run the installed build: install before running.

test_that("mrct_region_sims() leaves an unset progressr.handlers unset", {
  skip_on_cran()
  mrct_skip_if_no_build()
  own <- getOption("progressr.handlers")
  on.exit(options(progressr.handlers = own), add = TRUE)
  options(progressr.handlers = NULL)
  expect_null(getOption("progressr.handlers"))
  res <- suppressMessages(mrct_rs_call(mrct_build_S2()))
  expect_s3_class(res, "data.table")
  expect_null(getOption("progressr.handlers"))
})

test_that("mrct_region_sims() leaves the caller's progressr handler in place", {
  skip_on_cran()
  mrct_skip_if_no_build()
  own <- getOption("progressr.handlers")
  on.exit(options(progressr.handlers = own), add = TRUE)
  progressr::handlers("void")
  before <- getOption("progressr.handlers")
  expect_false(is.null(before))
  res <- suppressMessages(mrct_rs_call(mrct_build_S2()))
  expect_s3_class(res, "data.table")
  expect_identical(getOption("progressr.handlers"), before)
})
