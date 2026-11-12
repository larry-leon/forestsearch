# ============================================================================
# test-lasso-intercept.R
#
# A-042 (F6, CRAN preparation part 4): LASSO selection, the intercept is not a
# selected variable.
#
# For the binary, continuous and count outcome types, lasso_selection() fits
# a glmnet model with an (unpenalised) intercept, and coef() of that fit has a
# first row named "(Intercept)".  `selected` and `omitted` describe the
# CONFOUNDERS that LASSO keeps or drops, so the intercept belongs in neither;
# it stays in `coefficients`.  With the intercept out of `selected`, the rule
# in get_FSdata() "if LASSO selects nothing, keep the full candidate pool"
# applies to every outcome type, not only to survival (where a Cox fit has no
# intercept row).  forestsearch(details = TRUE) labels its LASSO summary by
# outcome type, as get_FSdata() already does.
#
# Data come from helper-synthetic-dgm.R (fixed seeds).  All calls are
# sequential.
# ============================================================================

.f6_confs <- c("age", "biomarker", "biomarker_hi", "sex")

# The three non-survival data sets, with the arguments lasso_selection() and
# get_FSdata() need for each.  For continuous and count outcomes there is no
# censoring indicator, so a placeholder event column of 1s is added (this is
# what forestsearch() does before it calls get_FSdata()).
.f6_glm_cases <- function() {
  bin <- .make_binary_data(N = 200L, seed = 42L)
  con <- .make_continuous_data(N = 200L, seed = 42L)
  con$event <- 1L
  cnt <- .make_count_data(N = 200L, seed = 42L)
  cnt$event <- 1L
  list(
    binary     = list(df = bin, event.name = "y",     offset.name = NULL),
    continuous = list(df = con, event.name = "event", offset.name = NULL),
    count      = list(df = cnt, event.name = "event", offset.name = "ftime")
  )
}

# Variable named on the left of a candidate cut ("age <= 55" -> "age"; a
# bare categorical name is its own variable).
.f6_cut_vars <- function(confs) {
  unique(trimws(sub("<=.*$", "", confs)))
}


test_that("lasso_selection(): (Intercept) is in coefficients, not in selected or omitted [A-042]", {
  cases <- .f6_glm_cases()
  for (ot in names(cases)) {
    cs <- cases[[ot]]
    res <- lasso_selection(
      df               = cs$df,
      confounders.name = .f6_confs,
      outcome.name     = "y",
      event.name       = cs$event.name,
      outcome_type     = ot,
      offset.name      = cs$offset.name
    )

    # glmnet documentation (?glmnet, argument `intercept`, default TRUE):
    # gaussian / binomial / poisson fits include an intercept, and coef()
    # returns it as the first row, named "(Intercept)", followed by one row
    # per column of x.  So `coefficients` has 1 + 4 entries, the first being
    # the intercept, the rest the four confounders in the order supplied.
    expect_identical(names(res$coefficients), c("(Intercept)", .f6_confs),
                     info = ot)

    # The fix's specification: the intercept is not a confounder, so it is in
    # neither `selected` nor `omitted`.
    expect_false("(Intercept)" %in% res$selected, info = ot)
    expect_false("(Intercept)" %in% res$omitted,  info = ot)

    # From the roxygen ("selected, omitted variables") and the spec:
    # `selected` and `omitted` split the supplied confounders, with no
    # overlap, and nothing else.
    expect_length(intersect(res$selected, res$omitted), 0L)
    expect_setequal(c(res$selected, res$omitted), .f6_confs)

    # By definition of LASSO selection at lambda.min: a confounder is
    # selected exactly when its coefficient is non-zero.
    conf_coef <- res$coefficients[.f6_confs]
    expect_identical(res$selected, .f6_confs[conf_coef != 0], info = ot)
    expect_identical(res$omitted,  .f6_confs[conf_coef == 0], info = ot)
  }
})


test_that("get_FSdata(): lasso_selection() mocked to select nothing keeps the full pool [A-042]", {
  # The mock is what a non-survival LASSO fit that keeps no covariate
  # returns under the corrected contract: nothing selected, every confounder
  # omitted, the fitted intercept kept in `coefficients` (value arbitrary).
  local_mocked_bindings(
    lasso_selection = function(df, confounders.name, ...) {
      list(selected     = character(0),
           omitted      = confounders.name,
           coefficients = c("(Intercept)" = 0.3,
                            stats::setNames(rep(0, length(confounders.name)),
                                            confounders.name)),
           lambda_min   = NA_real_,
           cvfit        = NULL,
           fit          = NULL)
    }
  )
  cases <- .f6_glm_cases()
  for (ot in names(cases)) {
    cs <- cases[[ot]]
    common <- list(df.analysis = cs$df, confounders.name = .f6_confs,
                   outcome.name = "y", event.name = cs$event.name,
                   outcome_type = ot, offset.name = cs$offset.name,
                   use_grf = FALSE, cut_type = "default",
                   collapse_cuts = FALSE, details = FALSE)
    fs_lasso <- do.call(get_FSdata, c(common, list(use_lasso = TRUE)))
    fs_none  <- do.call(get_FSdata, c(common, list(use_lasso = FALSE)))

    # Nothing selected, so lassokeep is empty and every confounder omitted.
    expect_length(fs_lasso$lassokeep, 0L)
    expect_identical(fs_lasso$lassoomit, .f6_confs, info = ot)

    # Full candidate pool, worked out from get_FSdata()'s documented cut
    # construction: biomarker_hi and sex are factors with 2 levels
    # (cont.cutoff = 4), so each enters as itself; age and biomarker are
    # continuous, so each enters through default cuts (mean, median, qlow,
    # qhigh; cut_var()).  Every one of the four confounders is therefore
    # represented, and the two factors appear verbatim.
    expect_setequal(.f6_cut_vars(fs_lasso$confs), .f6_confs)
    expect_true(all(c("biomarker_hi", "sex") %in% fs_lasso$confs), info = ot)
    expect_true(any(grepl("^age <= ", fs_lasso$confs)), info = ot)
    expect_true(any(grepl("^biomarker <= ", fs_lasso$confs)), info = ot)

    # "If LASSO selects nothing, keep the full pool": with no confounder
    # screened out, the candidate cuts are those built with LASSO off
    # (cut_type = "default" builds default cuts for every continuous
    # confounder on both branches of get_FSdata()).
    expect_identical(fs_lasso$confs, fs_none$confs, info = ot)
  }
})


test_that("get_FSdata(): a real non-survival LASSO fit that keeps no covariate keeps the full pool [A-042]", {
  # Here lasso_selection() itself runs (cv.glmnet + glmnet).  Only coef() of
  # the glmnet fit is mocked so that every covariate row is 0 while the
  # fitted intercept is left as it is: exactly the case "LASSO at lambda.min
  # keeps no covariate".  No other cut source is supplied (no conf_force,
  # defaultcut_names, conf.cont_jcuts, conf.cont_medians_force, GRF or DINA
  # cuts), so the candidate pool can only come from the LASSO pool.
  real_coef <- stats::coef
  local_mocked_bindings(
    coef = function(object, ...) {
      m <- real_coef(object, ...)
      if (inherits(object, "glmnet")) {
        m[rownames(m) != "(Intercept)", ] <- 0
      }
      m
    }
  )
  cases <- .f6_glm_cases()
  for (ot in names(cases)) {
    cs <- cases[[ot]]
    res <- get_FSdata(df.analysis = cs$df, confounders.name = .f6_confs,
                      outcome.name = "y", event.name = cs$event.name,
                      outcome_type = ot, offset.name = cs$offset.name,
                      use_lasso = TRUE, use_grf = FALSE,
                      cut_type = "default", collapse_cuts = FALSE,
                      details = FALSE)

    # Every covariate coefficient is 0 by construction of the mock, so LASSO
    # selects nothing and omits all four confounders (in the order supplied,
    # which is the row order of coef()).
    expect_length(res$lassokeep, 0L)
    expect_identical(res$lassoomit, .f6_confs, info = ot)

    # Full pool, by the same hand derivation as in the previous test: all
    # four confounders represented, both factors verbatim.
    expect_setequal(.f6_cut_vars(res$confs), .f6_confs)
    expect_true(all(c("biomarker_hi", "sex") %in% res$confs), info = ot)
  }
})


test_that("forestsearch(details = TRUE): LASSO summary is labelled by outcome type [A-042]", {
  # Labels as get_FSdata() prints them (R/get_fsdata.R, switch on
  # outcome_type): Logistic- / Linear- / Poisson-LASSO, Cox-LASSO for
  # survival.  The count "of N candidate factors" is
  # length(lassokeep) + length(lassoomit) = the 4 confounders supplied
  # (the intercept is not a candidate factor).
  expect_label <- list(binary = "Logistic-LASSO", continuous = "Linear-LASSO",
                       count = "Poisson-LASSO", survival = "Cox-LASSO")
  data_for <- list(
    binary     = .make_binary_data(N = 200L, seed = 42L),
    continuous = .make_continuous_data(N = 200L, seed = 42L),
    count      = .make_count_data(N = 200L, seed = 42L),
    survival   = .make_survival_data(N = 200L, seed = 42L)
  )
  confs_for <- list(binary = .f6_confs, continuous = .f6_confs,
                    count = .f6_confs,
                    survival = c("age", "stage", "sex", "noise"))
  for (ot in names(expect_label)) {
    args <- .fs_args_for(ot, confounders = confs_for[[ot]],
                         extra = list(use_grf = FALSE, use_lasso = TRUE,
                                      details = TRUE, quiet = FALSE))
    out <- utils::capture.output(
      suppressWarnings(suppressMessages(
        do.call(forestsearch, c(list(df.analysis = data_for[[ot]]), args))
      ))
    )
    summary_line <- grep("LASSO selected: [0-9]+ of [0-9]+ candidate factors",
                         out, value = TRUE)
    expect_length(summary_line, 1L)
    expect_match(summary_line,
                 paste0("^", expect_label[[ot]],
                        " selected: [0-4] of 4 candidate factors$"),
                 info = ot)
  }
})
