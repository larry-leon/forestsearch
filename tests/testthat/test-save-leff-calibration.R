# Tests for save_leff_calibration(): output_dir is required (no default
# write to "_data" in the working directory).

test_that("save_leff_calibration() errors when output_dir is not supplied", {
  wd_files <- list.files(getwd(), all.files = TRUE)
  expect_error(
    save_leff_calibration(C = 0.22, alpha = 1.3, n_min = 60,
                          outcome_type = "binary"),
    "output_dir.*must be supplied.*\"_data\""
  )
  # nothing was created in the working directory
  expect_identical(list.files(getwd(), all.files = TRUE), wd_files)
})

test_that("save_leff_calibration() writes to the supplied output_dir", {
  out_dir <- tempfile("leff_cal_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  path <- suppressMessages(
    save_leff_calibration(C = 0.22, alpha = 1.3, n_min = 60,
                          outcome_type = "binary", output_dir = out_dir)
  )
  expect_identical(path, file.path(out_dir, "calibration_binary_leff.rds"))
  expect_true(file.exists(path))
  bundle <- readRDS(path)
  expect_identical(bundle$C, 0.22)
  expect_identical(bundle$outcome_type, "binary")
})
