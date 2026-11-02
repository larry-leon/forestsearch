# mrct_region_sims(record = "full") with fs_args mr_inference = TRUE
# (TASK_fs-exact-betaHhat-mr-record_2026-11-01_r2, sections 2 and 3.3): the
# fourteen MR columns are filled where declared and NA where not; record =
# "point" with the same fs_args returns today's columns; a "full" call without
# mr_inference has the columns all NA. The call is mrct_rs_call()
# (helper-mrct-region-sims-record.R) on the S2 build at 2 replicates, 2
# workers. The workers run the installed build: install before running.

mr_cols <- c("mr_ci_method", "mr_debiased_H", "mr_ij_lower_H", "mr_ij_upper_H",
             "mr_field_lower_1s_H", "mr_debiased_Hc", "mr_ij_lower_Hc", "mr_ij_upper_Hc",
             "mr_field_upper_1s_Hc", "mr_field_upper_1s_s_Hc", "mr_bonf_lower_H",
             "mr_bonf_upper_Hc", "mr_p_hat_H", "mr_selection_rate")
mr_point_cols <- c("sim", "n_itt", "hr_itt", "hr_ittX", "n_train", "hr_train", "n_test",
                   "hr_test", "any_found", "sg_found", "n_sg", "hr_sg", "POhr_sg", "prev_sg",
                   "n_sg_train", "hr_sg_train", "POhr_sg_train", "hr_sg_null")
mr_count_cols <- paste0(rep(c("TP", "FP", "FN", "TN"), 2),
                        rep(c("_nonregion", "_region"), each = 4))
mr_fs_extra <- list(mr_inference = TRUE,
                    mr_inference_args = list(ci_method = "field", include_complement = TRUE,
                                             field_complement = TRUE, draws = 200,
                                             field_R_out = 50, field_R_in = 25, seed = 1))

mr_memo <- new.env()
mr_get <- function(key, expr) {
  if (is.null(mr_memo[[key]])) mr_memo[[key]] <- expr
  mr_memo[[key]]
}
mr_full <- function() mr_get("full", mrct_rs_call(mrct_build_S2(), record = "full",
                                                  n_sims = 2L, fs_args_extra = mr_fs_extra))
mr_point <- function() mr_get("point", mrct_rs_call(mrct_build_S2(), record = "point",
                                                    n_sims = 2L, fs_args_extra = mr_fs_extra))
mr_full_nomr <- function() mr_get("nomr", mrct_rs_call(mrct_build_S2(), record = "full",
                                                       n_sims = 2L))

test_that("record = 'full' with MR: the columns come after the counts, filled where declared", {
  mrct_skip_if_no_build()
  pf <- mr_full()
  expect_identical(names(pf), c(mr_point_cols, mr_count_cols, mr_cols))
  expect_type(pf$mr_ci_method, "character")
  for (k in setdiff(mr_cols, "mr_ci_method")) expect_type(pf[[k]], "double")
  dec <- pf$any_found == 1
  expect_true(any(dec))
  d <- pf[dec]
  expect_identical(d$mr_ci_method, rep("field", nrow(d)))
  filled <- setdiff(mr_cols, "mr_ci_method")
  for (k in filled) expect_true(all(is.finite(d[[k]])), info = k)
  # forestsearch() passes return_reselection = TRUE unless told otherwise
  # (R/forestsearch_main.R), so p_hat is present: a share of the draws, at
  # most the selection rate (sum(p_hat) == selection_rate)
  expect_true(all(d$mr_p_hat_H > 0 & d$mr_p_hat_H <= d$mr_selection_rate))
  # ordering of the bounds on the HR scale
  expect_true(all(d$mr_ij_lower_H <= d$mr_debiased_H & d$mr_debiased_H <= d$mr_ij_upper_H))
  expect_true(all(d$mr_ij_lower_Hc <= d$mr_debiased_Hc & d$mr_debiased_Hc <= d$mr_ij_upper_Hc))
  expect_true(all(d$mr_bonf_lower_H <= d$mr_field_lower_1s_H))
  expect_true(all(d$mr_bonf_upper_Hc >= d$mr_field_upper_1s_Hc))
  expect_true(all(d$mr_selection_rate > 0 & d$mr_selection_rate <= 1))
  # undeclared trials, if any, are NA throughout
  if (any(!dec)) for (k in mr_cols) expect_true(all(is.na(pf[!dec][[k]])), info = k)
  # the point columns are those of the point call with the same fs_args
  pp <- mr_point()
  expect_identical(as.list(pf)[mr_point_cols], as.list(pp))
})

test_that("record = 'full' with p_hat: mr_p_hat_H is the selected candidate's re-selection share", {
  # mr_p_hat_H is read from reselection$p_hat[selected_index], present only
  # under return_reselection = TRUE; the helper is checked on a stand-in.
  mr <- list(ci_method = "ij", selected_index = 2L, selection_rate = 0.9,
             debiased = list(est = 1.5, lower = 0.8, upper = 2.8),
             reselection = list(p_hat = c(a = 0.1, b = 0.3, c = 0.5)))
  r <- .mrct_mr_record(mr)
  expect_identical(names(r), mr_cols)
  expect_identical(r$mr_p_hat_H, 0.3)
  expect_identical(r$mr_ci_method, "ij")
  expect_identical(r$mr_debiased_H, 1.5)
  # no complement, no field: NA
  for (k in c("mr_debiased_Hc", "mr_ij_lower_Hc", "mr_ij_upper_Hc", "mr_field_lower_1s_H",
              "mr_field_upper_1s_Hc", "mr_field_upper_1s_s_Hc", "mr_bonf_lower_H",
              "mr_bonf_upper_Hc")) expect_true(is.na(r[[k]]), info = k)
  # NULL and the short return are all NA
  for (m in list(NULL, list(selected_index = NA_integer_, note = "x", n_family = 3L))) {
    r0 <- .mrct_mr_record(m)
    expect_identical(names(r0), mr_cols)
    expect_true(all(vapply(r0, function(v) is.na(v), NA)))
  }
})

test_that("record = 'point' with MR fs_args returns exactly today's columns", {
  mrct_skip_if_no_build()
  pp <- mr_point()
  expect_identical(names(pp), mr_point_cols)
  expect_null(attr(pp, "estimators"))
})

test_that("record = 'full' without mr_inference has the MR columns all NA", {
  mrct_skip_if_no_build()
  pf <- mr_full_nomr()
  expect_identical(names(pf), c(mr_point_cols, mr_count_cols, mr_cols))
  for (k in mr_cols) expect_true(all(is.na(pf[[k]])), info = k)
  expect_type(pf$mr_ci_method, "character")
})
