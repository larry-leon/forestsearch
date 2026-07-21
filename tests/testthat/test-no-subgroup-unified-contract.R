# ============================================================================
# test-no-subgroup-unified-contract.R
#
# Unified no-subgroup return contract (Option i): ALL THREE identifiers --
# consistency, DINA, and GRF -- return df.est = NULL, df.predict = NULL,
# df.test = NULL, and sg.harm = NULL when no subgroup is identified.
#
# Previously only the consistency path returned df.est = NULL; the DINA and
# GRF selection paths returned a populated df.est WITHOUT a treat.recommend
# column, which forced per-identifier guards in every downstream consumer.
# The divergence is now unified at the source (.forestsearch_dina_select /
# .forestsearch_grf_select no-detection returns), so detection keys uniformly
# on sg.harm and the no-subgroup df.est is uniformly NULL.
#
# Found fits are unchanged (only the no-detection return was touched).
# ============================================================================


# No-subgroup fit for a given identifier, forced with extreme thresholds.
.make_null_fit_contract <- function(method) {
  extra <- switch(method,
    consistency = list(hr.threshold = 10.0, pconsistency.threshold = 0.99),
    dina        = list(hr.threshold = 10.0, pconsistency.threshold = 0.99,
                       subgroup_method = "dina"),
    grf         = list(hr.threshold = 10.0, pconsistency.threshold = 0.99,
                       subgroup_method = "grf", dmin.grf = 1e6))
  df <- .make_survival_data(N = 160L, HR_harm = 1.0, seed = 7L)
  args <- .fs_args_for("survival", confounders = c("age", "stage", "sex"),
                       extra = extra)
  fs <- suppressWarnings(do.call(forestsearch, c(list(df.analysis = df), args)))
  if (!is.null(fs$sg.harm)) {
    skip(sprintf("null DGM unexpectedly identified a subgroup (%s)", method))
  }
  fs
}


for (method in c("consistency", "dina", "grf")) {

  test_that(sprintf("no-subgroup %s fit returns the unified NULL contract", method), {
    skip_on_cran()
    fs <- .make_null_fit_contract(method)

    # The unified contract: all frame slots NULL, sg.harm NULL.
    expect_null(fs$sg.harm)
    expect_null(fs$df.est)
    expect_null(fs$df.predict)
    expect_null(fs$df.test)

    # args_call_all must survive for downstream CV/bootstrap consumers.
    expect_false(is.null(fs$args_call_all))
  })
}


test_that("the three no-subgroup returns are contract-identical on the frame slots", {
  skip_on_cran()
  slots <- function(fs) list(sg.harm = fs$sg.harm, df.est = fs$df.est,
                             df.predict = fs$df.predict, df.test = fs$df.test)
  s_cons <- slots(.make_null_fit_contract("consistency"))
  s_dina <- slots(.make_null_fit_contract("dina"))
  s_grf  <- slots(.make_null_fit_contract("grf"))

  # All four slots are NULL in every identifier -> identical across the three.
  expect_identical(s_cons, s_dina)
  expect_identical(s_cons, s_grf)
})
