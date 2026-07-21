# ============================================================================
# test-no-subgroup-bootstrap.R
#
# Precondition guard for the exported low-level bootstrap worker
# bootstrap_results() on a no-subgroup fit.
#
# Every bootstrap iteration splits df_boot on the treat.recommend column
# (subgroups H / H^c).  A no-subgroup fit -- in particular a DINA/GRF fit
# whose df.est is populated but carries NO treat.recommend column -- has
# nothing to bias-correct.  Called directly (bypassing the guarded entry
# point forestsearch_bootstrap_dofuture()), each iteration used to fail
# with the cryptic "object 'treat.recommend' not found", surfaced by
# .collect_bootstrap_results() only as the aggregate "All N bootstrap
# iterations failed".  bootstrap_results() now fails fast with the same
# clear precondition message the entry point uses.
#
# The entry point forestsearch_bootstrap_dofuture() already guards this
# case (its own no-subgroup stop at SECTION 1); that behaviour is not
# changed here and is asserted for completeness.
#
# A subgroup-found fit carries treat.recommend, so the guard never fires
# for the intended input -- the bootstrap path is byte-identical there.
# ============================================================================


.make_null_dg_fit_boot <- function(method) {
  forcing <- if (method == "grf") {
    list(hr.threshold = 10.0, pconsistency.threshold = 0.99,
         subgroup_method = "grf", dmin.grf = 1e6)
  } else {
    list(hr.threshold = 10.0, pconsistency.threshold = 0.99,
         subgroup_method = "dina")
  }
  df <- .make_survival_data(N = 160L, HR_harm = 1.0, seed = 7L)
  args <- .fs_args_for("survival", confounders = c("age", "stage", "sex"),
                       extra = forcing)
  fs <- suppressWarnings(do.call(forestsearch, c(list(df.analysis = df), args)))
  if (!is.null(fs$sg.harm)) {
    skip(sprintf("null DGM unexpectedly identified a subgroup (%s)", method))
  }
  # Unified no-subgroup contract (Option i): df.est = NULL for all identifiers.
  expect_null(fs$df.est)
  fs
}


for (method in c("dina", "grf")) {

  test_that(sprintf("bootstrap_results() fails fast with a clear message for no-subgroup %s fit",
                    method), {
    skip_on_cran()
    fs <- .make_null_dg_fit_boot(method)

    err <- tryCatch(
      bootstrap_results(
        fs.est           = fs,
        df_boot_analysis = fs$df.est,
        cox.formula.boot = survival::Surv(time, event) ~ treat,
        nb_boots         = 5L,
        show_three       = FALSE,
        H_obs            = log(1.5),
        Hc_obs           = log(0.9)),
      error = function(e) e)
    expect_s3_class(err, "error")
    msg <- conditionMessage(err)

    # Clear precondition message, NOT the cryptic per-iteration failure
    # surfaced as the aggregate "All N bootstrap iterations failed".
    expect_match(msg, "no identified subgroup to bias-correct", fixed = TRUE)
    expect_false(grepl("All 5 bootstrap iterations failed", msg, fixed = TRUE))
    expect_false(grepl("object 'treat.recommend' not found", msg, fixed = TRUE))
  })
}


test_that("forestsearch_bootstrap_dofuture() still fails fast on a no-subgroup fit", {
  skip_on_cran()
  fs <- .make_null_dg_fit_boot("dina")

  err <- tryCatch(
    suppressWarnings(forestsearch_bootstrap_dofuture(
      fs.est = fs, nb_boots = 20L,
      parallel_args = list(plan = "sequential", workers = 1L,
                           show_message = FALSE))),
    error = function(e) e)
  expect_s3_class(err, "error")
  # Under the unified no-subgroup contract (Option i) df.est = NULL, so the
  # entry point halts at its df.est non-empty-frame validator ("must be a
  # non-empty data frame") -- the same clean halt the consistency path always
  # produced -- before reaching the more specific "no subgroup was identified"
  # message.  Either is a correct fail-fast on a no-subgroup fit.
  expect_match(conditionMessage(err),
               "no subgroup was identified|must be a non-empty data frame")
})
