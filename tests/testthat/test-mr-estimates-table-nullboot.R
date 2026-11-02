# mr_estimates_table(fs, boot_results = NULL): the table built from the
# forestsearch() result alone (TASK_fs-mrct-article_2026-11-01, section 1).
# Contracts checked on a small forestsearch(mr_inference = TRUE) fit with the
# field method at reduced field_R_out / field_R_in:
#   1. the rows (selected subgroup, complement) and the columns (descriptive,
#      the robust naive effect, MR de-biased with IJ, field bound, Bonferroni
#      pair) are present, with the robust header;
#   2. the numbers equal the elements of fs$mr_inference they come from
#      (1e-12), and the displayed cells are those numbers at `digits`;
#   3. without include_complement the complement cells are marked and the
#      complement-only columns do not appear when their quantities are absent;
#   4. supplying boot_results is still required to have FSsg_tab / SG_CIs.

CONF_MRT <- c("age", "stage", "sex", "noise")
.mrt_data <- .make_survival_data(N = 400L, HR_harm = 3.0)
.mrt_fit <- function(...) {
  args <- .fs_args_for("survival", confounders = CONF_MRT,
                       extra = list(use_grf = FALSE, use_lasso = FALSE,
                                    quiet = TRUE, mr_inference = TRUE,
                                    mr_inference_args = list(...)))
  .run_fs_capture(.mrt_data, args)$value
}

.mrt_field <- .mrt_fit(draws = 200L, seed = 11L, ci_method = "field",
                       include_complement = TRUE, field_complement = TRUE,
                       field_R_out = 60L, field_R_in = 30L)

.mrt_cells <- function(tbl) as.data.frame(tbl[["_data"]], stringsAsFactors = FALSE)
.f2 <- function(x) formatC(x, format = "f", digits = 2)

test_that("boot_results = NULL builds the table with the expected rows and columns", {
  fs <- .mrt_field
  skip_if(is.null(fs$sg.harm), "no subgroup identified on the synthetic data")
  g <- fs$mr_inference
  expect_false(is.null(g$field)); expect_false(is.null(g$field$complement))
  expect_false(is.null(g$field$joint))

  expect_identical(formals(mr_estimates_table)$boot_results, NULL)
  tbl <- mr_estimates_table(fs)
  expect_s3_class(tbl, "gt_tbl")
  d <- .mrt_cells(tbl)
  expect_identical(d$Subgroup, c("Qstnbl", "Recmnd"))
  expect_identical(colnames(d),
                   c("Subgroup", "n", "n1", "events", "m1", "m0", "RMST",
                     "HR (95% CI)", "MR_debiased", "MR_field", "MR_bonferroni"))
  # the naive column's header says robust
  lab <- tbl[["_boxhead"]]$column_label[tbl[["_boxhead"]]$var == "HR (95% CI)"][[1]]
  expect_match(as.character(lab), "robust")
  # sizes: the selected subgroup and its complement in the analysis data
  # (the n cells read "count (percent)")
  n_cell <- as.integer(sub("^\\s*([0-9]+).*$", "\\1", d$n))
  expect_equal(n_cell, c(sum(fs$df.est$treat.recommend == 0),
                         sum(fs$df.est$treat.recommend == 1)))
  expect_equal(n_cell[1], as.integer(g$n_selected))
  expect_equal(n_cell[2], as.integer(g$complement$n))
})

test_that("the table's numbers are the elements of fs$mr_inference (1e-12)", {
  fs <- .mrt_field
  skip_if(is.null(fs$sg.harm), "no subgroup identified on the synthetic data")
  g <- fs$mr_inference; cc <- g$complement; fl <- g$field
  tbl <- mr_estimates_table(fs, est.scale = "hr", digits = 2)
  v <- attr(tbl, "mr_values")
  expect_identical(v$row, c("H", "Hc"))
  ref <- rbind(
    H  = c(g$naive$est, g$naive$lower, g$naive$upper,
           g$debiased$est, g$debiased$lower, g$debiased$upper,
           fl$lower_1s, fl$joint$bonf_lower_H),
    Hc = c(cc$naive$est, cc$naive$lower, cc$naive$upper,
           cc$debiased$est, cc$debiased$lower, cc$debiased$upper,
           fl$complement$upper_1s, fl$joint$bonf_upper_Hc))
  got <- as.matrix(v[, c("naive_est", "naive_lower", "naive_upper",
                         "debiased_est", "ij_lower", "ij_upper",
                         "field_bound", "bonferroni_bound")])
  expect_true(all(is.finite(ref)))
  expect_lte(max(abs(got - ref)), 1e-12)

  # the displayed cells are those numbers at digits = 2
  d <- .mrt_cells(tbl)
  ci <- function(e, l, u) sprintf("%s (%s, %s)", .f2(e), .f2(l), .f2(u))
  expect_identical(d[["HR (95% CI)"]], ci(ref[, 1], ref[, 2], ref[, 3]))
  expect_identical(d$MR_debiased, ci(ref[, 4], ref[, 5], ref[, 6]))
  expect_identical(d$MR_field, c(paste("≥", .f2(ref[1, 7])), paste("≤", .f2(ref[2, 7]))))
  expect_identical(d$MR_bonferroni, c(paste("≥", .f2(ref[1, 8])), paste("≤", .f2(ref[2, 8]))))
  # the robust naive limits are those of the robust Cox fit on the subgroup
  H <- fs$df.est[fs$df.est$treat.recommend == 0, ]
  f <- summary(survival::coxph(survival::Surv(time, event) ~ treat, data = H, robust = TRUE))
  expect_equal(unname(f$conf.int[1, c(1, 3, 4)]), unname(ref["H", 1:3]), tolerance = 1e-6)

  # 1/hr inverts: estimates, swapped limits, bounds change side
  ti <- .mrt_cells(mr_estimates_table(fs, est.scale = "1/hr"))
  expect_identical(ti$MR_debiased[1], ci(1 / ref[1, 4], 1 / ref[1, 6], 1 / ref[1, 5]))
  expect_identical(ti$MR_field[1], paste("≤", .f2(1 / ref[1, 7])))
})

test_that("without include_complement the complement is marked and complement columns stay out", {
  fs <- .mrt_fit(draws = 200L, seed = 11L, ci_method = "field",
                 include_complement = FALSE, field_complement = FALSE,
                 field_R_out = 60L, field_R_in = 30L)
  skip_if(is.null(fs$sg.harm), "no subgroup identified on the synthetic data")
  g <- fs$mr_inference
  expect_null(g$complement); expect_null(g$field$joint)
  tbl <- mr_estimates_table(fs)
  d <- .mrt_cells(tbl)
  expect_identical(colnames(d)[-(1:8)], c("MR_debiased", "MR_field"))
  expect_identical(d$MR_debiased[2], "—")
  expect_identical(d$MR_field[2], "—")
  v <- attr(tbl, "mr_values")
  expect_lte(abs(v$debiased_est[1] - g$debiased$est), 1e-12)
  expect_lte(abs(v$field_bound[1] - g$field$lower_1s), 1e-12)
  expect_true(all(is.na(v[2, c("naive_est", "debiased_est", "field_bound")])))

  # ci_method = "ij": no field block, so no field column
  fs_ij <- .mrt_fit(draws = 200L, seed = 11L, ci_method = "ij")
  skip_if(is.null(fs_ij$sg.harm))
  d_ij <- .mrt_cells(mr_estimates_table(fs_ij))
  expect_false(any(c("MR_field", "MR_bonferroni") %in% colnames(d_ij)))
  expect_true("MR_debiased" %in% colnames(d_ij))
})

test_that("a supplied boot_results must still be a template", {
  fs <- .mrt_field
  skip_if(is.null(fs$sg.harm))
  expect_error(mr_estimates_table(fs, boot_results = list()), "FSsg_tab and SG_CIs")
  expect_error(mr_estimates_table(list(mr_inference = NULL)), "mr_inference is NULL")
})
