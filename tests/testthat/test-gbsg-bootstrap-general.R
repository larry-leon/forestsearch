# ============================================================================
# test-gbsg-bootstrap-general.R
#
# generate_gbsg_bootstrap_general() takes its source data from
# survival::gbsg. With survival 3.8-12, data("gbsg", package = "survival")
# loads nothing (gbsg ships inside the 'cancer' data file), which made every
# call stop with "'data' must be a data frame" (audit finding A-004).
# ============================================================================

test_that("generate_gbsg_bootstrap_general() returns without error on a small call", {
  expect_no_error(
    out <- suppressMessages(
      generate_gbsg_bootstrap_general(n = 50, seed = 1, noise_level = 0.1)
    )
  )
  expect_s3_class(out, "data.frame")
  expect_identical(nrow(out), 50L)
  expect_identical(names(out), names(survival::gbsg))
  expect_equal(out$pid, 1:50)

  # Bounds of the continuous variables are preserved (preserve_bounds = TRUE)
  for (v in c("age", "size", "nodes", "pgr", "er", "rfstime")) {
    expect_true(all(out[[v]] >= min(survival::gbsg[[v]])), info = v)
    expect_true(all(out[[v]] <= max(survival::gbsg[[v]])), info = v)
  }
  # Binary variables stay 0/1
  for (v in c("meno", "hormon", "status")) {
    expect_true(all(out[[v]] %in% c(0, 1)), info = v)
  }
})

test_that("generate_gbsg_bootstrap_general() is reproducible for a fixed seed", {
  a <- suppressMessages(generate_gbsg_bootstrap_general(n = 40, seed = 7))
  b <- suppressMessages(generate_gbsg_bootstrap_general(n = 40, seed = 7))
  expect_identical(a, b)
})

test_that("generate_gbsg_bootstrap_general() default call returns 686 rows", {
  out <- suppressMessages(generate_gbsg_bootstrap_general())
  expect_identical(dim(out), c(686L, ncol(survival::gbsg)))
})
