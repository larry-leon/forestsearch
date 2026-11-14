# test-mr_values.R
#
# mr_values(): the tidy accessor for fit$mr_inference (TASK_mr_values_2026-11-06).
# The substance is the identity contract (D6): for every fixture, every non-NA cell
# equals the value read directly from the raw mr_inference paths, every other cell
# is NA and corresponds to a component the fit genuinely lacks, and the quantities
# shared with mr_estimates_table()'s "mr_values" attribute agree.
#
# F1: the full-coverage fixture (standalone fs_mr_inference() on a hand-built family).
# F2: the NA-fill fixture (ci_method = "ij", no complement).
# The variants cover each remaining absence condition of the field map.

skip_on_cran()

.mrv_data <- function(n = 300L, seed = 11L) {
  set.seed(seed)
  data.frame(time  = stats::rexp(n, 0.2),
             event = stats::rbinom(n, 1L, 0.8),
             treat = stats::rbinom(n, 1L, 0.5),
             x1 = stats::rnorm(n), x2 = stats::rnorm(n))
}
.mrv_spec <- list(outcome_type = "survival", effect_measure = "HR",
                  treat.name = "treat", outcome.name = "time",
                  event.name = "event", offset.name = NULL,
                  adjust_covariates = NULL, adverse_outcome = TRUE)
.mrv_df <- .mrv_data()
.mrv_cands <- list(g1 = which(.mrv_df$x1 > 0), g2 = which(.mrv_df$x2 > 0),
                   g3 = which(.mrv_df$x1 + .mrv_df$x2 > 0.3))
# Selected: the candidate with the largest naive log-HR (the maxeff rule below).
.mrv_sel <- names(.mrv_cands)[which.max(vapply(.mrv_cands, function(ix) {
  unname(stats::coef(survival::coxph(survival::Surv(time, event) ~ treat,
                                     data = .mrv_df[ix, ])))
}, numeric(1)))]
.mrv_call <- function(..., df = .mrv_df, cands = .mrv_cands, spec = .mrv_spec,
                      sel = cands[[.mrv_sel]]) {
  fs_mr_inference(df = df, candidates = cands, spec = spec, selected_members = sel,
                  admission = list(effect_floor = NULL, consistency = NULL),
                  reselection = "maxeff", draws = 400L, multiplier = "poisson",
                  seed = 4242L, field_R_out = 300L, field_R_in = 100L, ...)
}

.mrv_F1 <- .mrv_call(ci_method = "field", include_complement = TRUE,
                     field_complement = TRUE, field_scale_complement = "selected",
                     return_reselection = TRUE)
.mrv_F2 <- .mrv_call(ci_method = "ij", include_complement = FALSE)

.mrv_cols <- c("row", "n", "naive_est", "naive_lower", "naive_upper", "debiased_est",
               "ij_lower", "ij_upper", "ij_se", "field_est2", "field_est2_s",
               "field_lower_1s", "field_upper_1s", "field_upper_1s_s", "bonf_lower_H",
               "bonf_upper_Hc", "bonf_lower_H_s", "bonf_upper_Hc_s", "selection_bias",
               "p_hat", "n_family", "n_selected", "timing_seconds",
               "field_timing_seconds", "field_complement_timing_seconds")

# The raw paths, written out here independently of the implementation.
.mrv_raw_H_base <- function(g) {
  list(n = g$n_selected, naive_est = g$naive$est, naive_lower = g$naive$lower,
       naive_upper = g$naive$upper, debiased_est = g$debiased$est,
       ij_lower = g$debiased$lower, ij_upper = g$debiased$upper,
       ij_se = g$debiased$se_ij_two_term, selection_bias = g$selection_bias,
       p_hat = g$reselection$p_hat[[g$selected_index]], n_family = g$n_family,
       n_selected = g$n_selected, timing_seconds = g$timing_seconds)
}
.mrv_raw_H_field <- function(g) {
  list(field_est2 = g$field$est2, field_lower_1s = g$field$lower_1s,
       field_timing_seconds = g$field$timing_seconds)
}
.mrv_raw_H_cfield <- function(g) {
  list(bonf_lower_H = g$field$joint$bonf_lower_H,
       field_complement_timing_seconds = g$field$complement$timing_seconds)
}
.mrv_raw_H_s <- function(g) list(bonf_lower_H_s = g$field$joint_s$bonf_lower_H)
.mrv_raw_Hc_base <- function(g) {
  cc <- g$complement
  list(n = cc$n, naive_est = cc$naive$est, naive_lower = cc$naive$lower,
       naive_upper = cc$naive$upper, debiased_est = cc$debiased$est,
       ij_lower = cc$debiased$lower, ij_upper = cc$debiased$upper,
       ij_se = cc$debiased$se_ij_two_term)
}
.mrv_raw_Hc_cfield <- function(g) {
  fc <- g$field$complement
  list(field_est2 = fc$est2, field_upper_1s = fc$upper_1s,
       bonf_upper_Hc = g$field$joint$bonf_upper_Hc)
}
.mrv_raw_Hc_s <- function(g) {
  fc <- g$field$complement
  list(field_est2_s = fc$est2_s, field_upper_1s_s = fc$upper_1s_s,
       bonf_upper_Hc_s = g$field$joint_s$bonf_upper_Hc)
}

# Every cell of one row: a listed column equals its raw value (and the raw value
# exists); every other column is NA.
.mrv_expect_row <- function(row, raw) {
  stopifnot(nrow(row) == 1L)
  for (cn in setdiff(.mrv_cols, "row")) {
    if (cn %in% names(raw)) {
      expect_length(raw[[cn]], 1L)
      expect_false(is.na(raw[[cn]]), label = paste("raw", cn))
      expect_equal(as.numeric(row[[cn]]), as.numeric(raw[[cn]]), tolerance = 1e-12,
                   label = paste(row$row, cn))
    } else {
      expect_true(is.na(row[[cn]]), label = paste(row$row, cn, "is NA"))
    }
  }
}
.mrv_row <- function(v, r) v[v$row == r, , drop = FALSE]

test_that("F1 (full coverage): every cell equals its raw mr_inference path", {
  g <- .mrv_F1
  v <- mr_values(g)
  expect_identical(class(v), "data.frame")
  expect_identical(names(v), .mrv_cols)
  expect_identical(v$row, c("H", "Hc"))
  .mrv_expect_row(.mrv_row(v, "H"), c(.mrv_raw_H_base(g), .mrv_raw_H_field(g),
                                      .mrv_raw_H_cfield(g), .mrv_raw_H_s(g)))
  .mrv_expect_row(.mrv_row(v, "Hc"), c(.mrv_raw_Hc_base(g), .mrv_raw_Hc_cfield(g),
                                       .mrv_raw_Hc_s(g)))
})

test_that("F2 (NA-fill): field, joint and complement are NA because they are absent", {
  g <- .mrv_F2
  expect_null(g$field)
  expect_null(g$complement)
  v <- mr_values(g)
  expect_identical(names(v), .mrv_cols)
  expect_identical(v$row, "H")
  .mrv_expect_row(v, .mrv_raw_H_base(g))
  expect_true(all(is.na(v[grep("^(field|bonf)_", names(v))])))
})

test_that("each remaining absence condition gives NA, the rest stays identical", {
  # field_complement = FALSE: the complement exists, its field and the joint do not
  g <- .mrv_call(ci_method = "field", include_complement = TRUE,
                 field_complement = FALSE)
  expect_null(g$field$complement); expect_null(g$field$joint); expect_null(g$field$joint_s)
  v <- mr_values(g)
  .mrv_expect_row(.mrv_row(v, "H"), c(.mrv_raw_H_base(g), .mrv_raw_H_field(g)))
  .mrv_expect_row(.mrv_row(v, "Hc"), .mrv_raw_Hc_base(g))

  # field_scale_complement = "none": no field-s companions and no joint_s
  g <- .mrv_call(ci_method = "field", include_complement = TRUE,
                 field_scale_complement = "none")
  expect_null(g$field$joint_s); expect_null(g$field$complement$upper_1s_s)
  v <- mr_values(g)
  .mrv_expect_row(.mrv_row(v, "H"), c(.mrv_raw_H_base(g), .mrv_raw_H_field(g),
                                      .mrv_raw_H_cfield(g)))
  .mrv_expect_row(.mrv_row(v, "Hc"), c(.mrv_raw_Hc_base(g), .mrv_raw_Hc_cfield(g)))

  # return_reselection = FALSE: no p_hat
  g <- .mrv_call(ci_method = "ij", include_complement = FALSE, return_reselection = FALSE)
  expect_null(g$reselection)
  raw <- .mrv_raw_H_base(.mrv_F2)
  v <- mr_values(g)
  expect_true(is.na(v$p_hat))
  expect_equal(v$debiased_est, raw$debiased_est, tolerance = 1e-12)
})

test_that("only the two-term IJ construction is exposed", {
  # ci_method = "wald": debiased$lower / $upper hold the robust-SE interval
  g <- .mrv_call(ci_method = "wald", include_complement = TRUE)
  v <- mr_values(g)
  expect_true(all(is.na(v$ij_lower)))
  expect_true(all(is.na(v$ij_upper)))
  expect_equal(v$ij_se, c(g$debiased$se_ij_two_term, g$complement$debiased$se_ij_two_term),
               tolerance = 1e-12)
  expect_equal(v$debiased_est, c(g$debiased$est, g$complement$debiased$est), tolerance = 1e-12)

  # ij_residual = "winner": the reported SE and limits are the winner variant
  g <- .mrv_call(ci_method = "ij", include_complement = TRUE, ij_residual = "winner")
  expect_false(isTRUE(all.equal(g$debiased$se_ij, g$debiased$se_ij_two_term)))
  v <- mr_values(g)
  expect_true(all(is.na(v$ij_lower)))
  expect_true(all(is.na(v$ij_upper)))
  expect_equal(v$ij_se, c(g$debiased$se_ij_two_term, g$complement$debiased$se_ij_two_term),
               tolerance = 1e-12)
  winner <- c(unlist(g$debiased[c("se_ij", "se_ij_winner", "se_ij_winner_floor", "lower",
                                  "upper", "lower_w", "upper_w", "lower_1s_w", "lower_wf",
                                  "upper_wf", "lower_1s_wf")]))
  cells <- as.numeric(unlist(v[setdiff(names(v), "row")]))
  expect_false(any(cells %in% winner, na.rm = TRUE))

  # On the default fit too: no winner-variant name and no winner-variant value
  v <- mr_values(.mrv_F1)
  expect_false(any(grepl("winner|_w$|_wf$", names(v))))
  winner <- unlist(.mrv_F1$debiased[c("se_ij_winner", "se_ij_winner_floor", "lower_w",
                                      "upper_w", "lower_1s_w", "lower_wf", "upper_wf",
                                      "lower_1s_wf")])
  cells <- as.numeric(unlist(v[setdiff(names(v), "row")]))
  expect_false(any(cells %in% winner, na.rm = TRUE))
})

test_that("a complement that could not be fit: an all-NA Hc row and no spurious values", {
  # The selected subgroup leaves 3 patients, too few to fit: complement is note-only,
  # and `$n` on it would partially match `note`.
  sel <- setdiff(seq_len(nrow(.mrv_df)), 1:3)
  g <- .mrv_call(ci_method = "field", include_complement = TRUE,
                 cands = c(.mrv_cands, list(big = sel)), sel = sel)
  expect_identical(names(g$complement), "note")
  expect_false(is.null(g$field$complement$note))
  expect_no_warning(v <- mr_values(g))
  expect_identical(v$row, c("H", "Hc"))
  expect_true(all(is.na(.mrv_row(v, "Hc")[setdiff(.mrv_cols, "row")])))
  .mrv_expect_row(.mrv_row(v, "H"), c(.mrv_raw_H_base(g), .mrv_raw_H_field(g)))
})

test_that("the short return gives an H row that is NA apart from n_family", {
  g <- .mrv_call(ci_method = "field", include_complement = TRUE, sel = 1:3)
  expect_true(is.na(g$selected_index))
  expect_false(is.null(g$note))
  v <- mr_values(g)
  expect_identical(names(v), .mrv_cols)
  expect_identical(v$row, "H")
  .mrv_expect_row(v, list(n_family = g$n_family))
  expect_identical(names(mr_values(g, scale = "log")), .mrv_cols)
})

test_that("scale = 'log' transforms the effect columns of ratio measures only", {
  eff <- c("naive_est", "naive_lower", "naive_upper", "debiased_est", "ij_lower",
           "ij_upper", "field_est2", "field_est2_s", "field_lower_1s", "field_upper_1s",
           "field_upper_1s_s", "bonf_lower_H", "bonf_upper_Hc", "bonf_lower_H_s",
           "bonf_upper_Hc_s")
  check_log <- function(g) {
    v <- mr_values(g); l <- mr_values(g, scale = "log")
    expect_identical(mr_values(g, scale = "effect"), v)
    expect_equal(as.matrix(l[eff]), log(as.matrix(v[eff])), tolerance = 1e-12)
    rest <- setdiff(names(v), eff)
    expect_identical(l[rest], v[rest])          # stored scale under both settings
  }
  expect_true(.mrv_F1$log_scale)
  check_log(.mrv_F1)

  # OR: a ratio measure on the GLM path
  set.seed(21)
  n <- 400L
  d <- data.frame(y = stats::rbinom(n, 1L, 0.4), trt = stats::rbinom(n, 1L, 0.5),
                  x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  fam <- list(a = which(d$x1 > 0), b = which(d$x2 > 0))
  spec_or <- list(outcome_type = "binary", effect_measure = "OR", treat.name = "trt",
                  outcome.name = "y", event.name = NULL, offset.name = NULL,
                  adjust_covariates = NULL, adverse_outcome = TRUE)
  g_or <- .mrv_call(ci_method = "field", include_complement = TRUE, df = d, cands = fam,
                    spec = spec_or, sel = fam$a)
  expect_identical(g_or$measure, "OR"); expect_true(g_or$log_scale)
  check_log(g_or)
})

test_that("MD: effect returns as stored, log is refused", {
  set.seed(31)
  n <- 400L
  d <- data.frame(y = stats::rnorm(n), trt = stats::rbinom(n, 1L, 0.5),
                  x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  fam <- list(a = which(d$x1 > 0), b = which(d$x2 > 0))
  spec_md <- list(outcome_type = "continuous", effect_measure = "MD", treat.name = "trt",
                  outcome.name = "y", event.name = NULL, offset.name = NULL,
                  adjust_covariates = NULL, adverse_outcome = TRUE)
  g <- .mrv_call(ci_method = "field", include_complement = TRUE, df = d, cands = fam,
                 spec = spec_md, sel = fam$a)
  expect_false(g$log_scale)
  v <- mr_values(g)
  .mrv_expect_row(.mrv_row(v, "H"), c(.mrv_raw_H_base(g), .mrv_raw_H_field(g),
                                      .mrv_raw_H_cfield(g), .mrv_raw_H_s(g)))
  .mrv_expect_row(.mrv_row(v, "Hc"), c(.mrv_raw_Hc_base(g), .mrv_raw_Hc_cfield(g),
                                       .mrv_raw_Hc_s(g)))
  expect_error(mr_values(g, scale = "log"), "ratio measures only")
})

test_that("input detection: a fit, a bare result, and informative errors otherwise", {
  fit <- structure(list(mr_inference = .mrv_F1), class = c("forestsearch", "list"))
  expect_identical(mr_values(fit), mr_values(.mrv_F1))
  no_mr <- structure(list(sg.harm = "x"), class = c("forestsearch", "list"))
  expect_error(mr_values(no_mr), "no mr_inference results")
  expect_error(mr_values(list(a = 1)), "forestsearch fit or a fs_mr_inference")
  expect_error(mr_values(data.frame(a = 1)), "forestsearch fit or a fs_mr_inference")
  expect_error(mr_values(NULL), "forestsearch fit or a fs_mr_inference")
  expect_error(mr_values(.mrv_F1, scale = "hr"))
})

test_that("the overlapping quantities equal mr_estimates_table()'s attribute", {
  skip_if_not_installed("gt")
  # The fit-shaped input mr_estimates_table() reads, around a standalone result.
  as_fs <- function(g) {
    d <- .mrv_df
    d$treat.recommend <- 1L
    d$treat.recommend[.mrv_cands[[.mrv_sel]]] <- 0L
    list(mr_inference = g, df.est = d, sg.harm = .mrv_sel,
         args_call_all = list(outcome_type = "survival", outcome.name = "time",
                              event.name = "event", treat.name = "treat"))
  }
  shared <- c("naive_est", "naive_lower", "naive_upper", "debiased_est", "ij_lower", "ij_upper")

  v <- mr_values(.mrv_F1)
  a <- attr(mr_estimates_table(as_fs(.mrv_F1), digits = 3), "mr_values")
  expect_identical(a$row, v$row)
  for (cn in shared) expect_equal(v[[cn]], a[[cn]], tolerance = 1e-12, label = cn)
  expect_equal(c(v$field_lower_1s[1], v$field_upper_1s[2]), a$field_bound, tolerance = 1e-12)
  expect_equal(c(v$bonf_lower_H[1], v$bonf_upper_Hc[2]), a$bonferroni_bound, tolerance = 1e-12)

  v <- mr_values(.mrv_F2)
  a <- attr(mr_estimates_table(as_fs(.mrv_F2), digits = 3), "mr_values")
  for (cn in shared) expect_equal(v[[cn]], a[[cn]][1], tolerance = 1e-12, label = cn)
  expect_true(all(is.na(a[a$row == "Hc", shared])))     # no Hc row in v: absent in both
  expect_false(any(c("field_bound", "bonferroni_bound") %in% names(a)))
})
