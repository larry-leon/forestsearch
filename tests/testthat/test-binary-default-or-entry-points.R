# ============================================================================
# test-binary-default-or-entry-points.R
#
# The rider of TASK_threshold_sync_2026-09-18.md (Step 4), following up
# Directive B's finding F2: the exported estimation-layer entry points still
# defaulted binary to "RD" while forestsearch() defaults it to "OR", so the
# same unset call meant a risk difference at one layer and an odds ratio at
# the other.
#
# No fit and no resample: make_effect_estimator() only builds a closure, and
# the resolution inside .consistency_glm_pieces() is evaluated as the single
# statement it is.
# ============================================================================

# The resolution site itself, evaluated as the single statement it is.  The
# returned closure only carries effect_measure on the branches that pass it to
# their maker (binary does; survival and continuous do not), so reading it back
# off the closure would not cover all four outcome types.
.bdep_resolution_stmt <- function(f) {
  stmts <- as.list(body(f))[-1L]
  txt   <- vapply(stmts, function(e) paste(deparse(e), collapse = " "),
                  character(1))
  ix <- which(grepl("is.null(effect_measure)", txt, fixed = TRUE))
  if (length(ix) != 1L)
    stop("expected exactly one effect_measure resolution site, found ",
         length(ix))
  stmts[[ix]]
}

.bdep_resolve <- function(f, outcome_type, effect_measure = NULL) {
  env <- new.env(parent = asNamespace("forestsearch"))
  env$outcome_type   <- outcome_type
  env$effect_measure <- effect_measure
  eval(.bdep_resolution_stmt(f), env)
  env$effect_measure
}


test_that("make_effect_estimator() resolves an unset binary measure to OR", {
  expect_identical(.bdep_resolve(make_effect_estimator, "binary"), "OR")

  # End to end for the branch that carries it: the closure the entry point
  # actually hands back is built for an odds ratio.
  est <- make_effect_estimator(outcome_type = "binary", treat.name = "Treat",
                               outcome.name = "Y")
  expect_identical(get0("effect_measure", envir = environment(est)), "OR")
})


test_that("make_effect_estimator() leaves an explicit measure alone", {
  for (em in c("RD", "RR", "OR")) {
    expect_identical(.bdep_resolve(make_effect_estimator, "binary", em), em)
    est <- make_effect_estimator(outcome_type = "binary", effect_measure = em,
                                 treat.name = "Treat", outcome.name = "Y")
    expect_identical(get0("effect_measure", envir = environment(est)), em)
  }
})


test_that("make_effect_estimator()'s other outcome-type defaults are unmoved", {
  expect_identical(.bdep_resolve(make_effect_estimator, "survival"),   "HR")
  expect_identical(.bdep_resolve(make_effect_estimator, "continuous"), "MD")
  expect_identical(.bdep_resolve(make_effect_estimator, "count"),      "IRR")
})


test_that("consistency_resample()'s GLM resolution defaults binary to OR", {
  # consistency_resample() keeps effect_measure = NULL and hands it to
  # .consistency_glm_pieces(), which holds the resolution.  Evaluate that one
  # statement rather than running the resample.
  pieces <- forestsearch:::.consistency_glm_pieces
  expect_null(eval(formals(consistency_resample)[["effect_measure"]]))

  expect_identical(.bdep_resolve(pieces, "binary"),     "OR")
  expect_identical(.bdep_resolve(pieces, "continuous"), "MD")
  expect_identical(.bdep_resolve(pieces, "count"),      "IRR")
})


test_that("an explicit measure survives .consistency_glm_pieces()' resolution", {
  pieces <- forestsearch:::.consistency_glm_pieces
  expect_identical(.bdep_resolve(pieces, "binary", "RD"), "RD")
  expect_identical(.bdep_resolve(pieces, "binary", "OR"), "OR")
})


test_that("consistency_resample_compare() carries no binary default at all", {
  # Directive B's F2 names three entry points, but this one has neither an
  # outcome_type nor an effect_measure formal: it is survival-only and calls
  # consistency_resample() without outcome_type, so it resolves through the
  # Cox branch.  There is nothing here to flip, and this test says so.
  fm <- names(formals(consistency_resample_compare))
  expect_false("effect_measure" %in% fm)
  expect_false("outcome_type" %in% fm)

  src <- paste(deparse(body(consistency_resample_compare)), collapse = "\n")
  expect_true(grepl("consistency_resample(", src, fixed = TRUE))
  expect_false(grepl("outcome_type", src, fixed = TRUE))
})


test_that("no \"RD\" binary default is left at the two entry points", {
  for (f in list(make_effect_estimator,
                 forestsearch:::.consistency_glm_pieces)) {
    src <- paste(deparse(body(f)), collapse = "\n")
    expect_match(src, "binary\\s*=\\s*\"OR\"")
    expect_false(grepl("binary\\s*=\\s*\"RD\"", src))
  }
})
