# Tests for .fs_grf_num_threads(): grf thread cap under _R_CHECK_LIMIT_CORES_.
# withr is not a declared dependency, so the variable is saved at the start
# of each test and restored on exit with base Sys.setenv()/Sys.unsetenv().

.restore_check_limit <- function(old) {
  if (is.na(old)) Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
  else Sys.setenv("_R_CHECK_LIMIT_CORES_" = old)
}

test_that(".fs_grf_num_threads() returns NULL when the variable is unset", {
  old <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA)
  on.exit(.restore_check_limit(old), add = TRUE)
  Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
  expect_null(.fs_grf_num_threads())
})

test_that(".fs_grf_num_threads() returns NULL for empty or false values", {
  old <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA)
  on.exit(.restore_check_limit(old), add = TRUE)
  for (v in c("", "false", "FALSE")) {
    Sys.setenv("_R_CHECK_LIMIT_CORES_" = v)
    expect_null(.fs_grf_num_threads())
  }
})

test_that(".fs_grf_num_threads() returns 2L when the variable is TRUE", {
  old <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA)
  on.exit(.restore_check_limit(old), add = TRUE)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")
  expect_identical(.fs_grf_num_threads(), 2L)
})

test_that(".fs_grf_num_threads() returns the user value when supplied", {
  old <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA)
  on.exit(.restore_check_limit(old), add = TRUE)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")
  expect_identical(.fs_grf_num_threads(4L), 4L)
  Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
  expect_identical(.fs_grf_num_threads(1L), 1L)
})
