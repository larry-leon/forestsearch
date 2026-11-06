# test-fs-declaration-table.R
#
# fs_declaration_table() (TASK_wrapper_batch_2026-11-06, B4): on a small fit with the
# declaration field, the article's S17 entries are recovered exactly. s17_column() is
# transplanted from vignettes/articles/actg175-glm-outcomes.qmd.

.fdt_memo <- new.env()
.fdt_fit <- function() {
  if (is.null(.fdt_memo$fit)) {
    set.seed(1); n <- 500
    d <- data.frame(id = seq_len(n), treat = stats::rbinom(n, 1, 0.5), x1 = stats::rnorm(n), x2 = stats::rnorm(n))
    d$time <- stats::rexp(n, 0.1 * ifelse(d$treat == 1, ifelse(d$x1 < -0.6, 3, 0.45), 1))
    d$event <- as.integer(d$time < 30); d$time <- pmin(d$time, 30)
    .fdt_memo$fit <- forestsearch(
      df.analysis = d, outcome.name = "time", event.name = "event", treat.name = "treat", id.name = "id",
      confounders.name = c("x1", "x2"), use_grf = FALSE, use_lasso = FALSE, maxk = 1, fs.splits = 200,
      sg_focus = "hrMaxSG", mr_inference = TRUE,
      mr_inference_args = list(ci_method = "ij", draws = 500, keep_declaration_field = TRUE,
                               declaration_c0 = c(0.6, 0.7)),
      parallel_args = list(plan = "sequential"), quiet = TRUE, details = FALSE, plot.sg = FALSE,
      show_candidate_summary = FALSE)
  }
  .fdt_memo$fit
}
# ---- transplanted: vignettes/articles/actg175-glm-outcomes.qmd, declaration
s17_column <- function(fit, c0_grid, c0_primary, fmt_c0) {
  cal    <- fs_declaration_calibration(fit, alpha = 0.10, c0 = c0_grid)
  cal_c2 <- fs_declaration_calibration(fit, alpha = 0.10)   # c0 = NULL: the claim threshold
  stopifnot(identical(cal$c0$source, "capture"))            # read from the fit's own draws
  fld <- fit$mr_inference$declaration_field
  T_H <- cal$T_hat[[fld$family_id[fld$meta$selected_index]]]
  tb  <- cal$c0$table; i <- which(tb$c0 == c0_primary)
  admits_none <- all(lengths(cal$c0$admitted_calibrated) == 0) &&
    length(cal_c2$admitted_calibrated) == 0
  c(fit$mr_inference$n_selected,
    sprintf("%.3f / %.3f", T_H, max(cal$T_hat, na.rm = TRUE)),
    fmt_c0(c0_primary),
    sprintf("%.3f", tb$kappa_hat[i]),
    sprintf("%.3f", 2 * pnorm(tb$kappa_hat[i]) - 1),
    sprintf("%.3f", tb$fw_size[i]),
    sprintf("%.3f", cal_c2$kappa_hat),
    sprintf("%.3f", cal_c2$fw_size),
    if (admits_none) "admits nothing" else "admits a candidate")
}
# the same nine entries, read from the table
.fdt_entries <- function(tab, fmt_c0) {
  p <- tab[tab$primary, ]; cl <- tab[tab$row == "claim_threshold", ]
  c(p$n_region, sprintf("%.3f / %.3f", p$T_region, p$T_max), fmt_c0(p$c0), sprintf("%.3f", p$kappa_hat),
    sprintf("%.3f", p$implied_level), sprintf("%.3f", p$fw_size), sprintf("%.3f", cl$kappa_hat),
    sprintf("%.3f", cl$fw_size), if (p$admits_none) "admits nothing" else "admits a candidate")
}

test_that("the article's S17 entries are recovered exactly from the table", {
  fit <- .fdt_fit()
  fmt <- function(x) sprintf("HR %.2f", x)
  for (primary in c(0.6, 0.7)) {
    tab <- fs_declaration_table(fit, alpha = 0.10, c0 = c(0.6, 0.7), c0_primary = primary)
    expect_identical(.fdt_entries(tab, fmt), s17_column(fit, c(0.6, 0.7), primary, fmt))
  }
})

test_that("the rows, columns and direct values", {
  fit <- .fdt_fit()
  tab <- fs_declaration_table(fit, c0 = c(0.6, 0.7), c0_primary = 0.7)
  expect_identical(class(tab), "data.frame")
  expect_identical(names(tab), c("row", "c0", "primary", "alpha", "n_region", "T_region", "T_max", "kappa_hat",
                                 "implied_level", "fw_size", "n_admitted", "admits_none", "source"))
  expect_identical(tab$row, c("c0", "c0", "claim_threshold"))
  expect_identical(tab$c0, c(0.6, 0.7, NA))
  expect_identical(tab$primary, c(FALSE, TRUE, FALSE))
  expect_identical(tab$source, c("capture", "capture", NA))
  cal <- fs_declaration_calibration(fit, alpha = 0.10, c0 = c(0.6, 0.7))
  cl <- fs_declaration_calibration(fit, alpha = 0.10)
  expect_identical(tab$kappa_hat, c(cal$c0$table$kappa_hat, cl$kappa_hat))
  expect_identical(tab$fw_size, c(cal$c0$table$fw_size, cl$fw_size))
  expect_identical(tab$implied_level, 2 * pnorm(tab$kappa_hat) - 1)
  expect_identical(tab$n_admitted, c(unname(lengths(cal$c0$admitted_calibrated)), length(cl$admitted_calibrated)))
  expect_identical(tab$admits_none, rep(all(tab$n_admitted == 0L), 3))
  expect_true(all(tab$n_region == fit$mr_inference$n_selected))
  # no grid: the claim-threshold row alone
  one <- fs_declaration_table(fit)
  expect_identical(one$row, "claim_threshold")
  expect_identical(one$kappa_hat, cl$kappa_hat)
  expect_identical(names(one), names(tab))
  expect_false(any(grepl("significan|winner", names(tab), ignore.case = TRUE)))
})

test_that("informative errors", {
  fit <- .fdt_fit()
  bare <- fit; bare$mr_inference$declaration_field <- NULL
  expect_error(fs_declaration_table(bare), "keep_declaration_field = TRUE")
  no_mr <- structure(list(sg.harm = "x"), class = c("forestsearch", "list"))
  expect_error(fs_declaration_table(no_mr), "no declaration field")
  expect_error(fs_declaration_table(list(a = 1)), "no declaration field")
  expect_error(fs_declaration_table(fit, c0 = c(0.6, 0.7), c0_primary = 0.5), "one of the levels")
  expect_error(fs_declaration_table(fit, c0_primary = 0.7), "one of the levels")
})
