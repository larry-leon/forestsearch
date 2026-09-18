# =============================================================================
# Acceptance tests -- TASK_directive_C_2026-09-18.md
#
# Part 1  the DINA identity-scale guard, at the m_diff DERIVATION SITE
# Part 2  the frontier-key warning under subgroup_method = "dina"
# Part 3  dina_frontier() display caps -> Inf / Inf, warn when a value trims
# Part 4  the details-time frontier print retitle
# Part 5  the c2 / p* echo annotation under dina / grf
#
# Compute budget: at most two seeded micro-fits in this file, hard abort 5 min;
# the file itself aborts at 3 minutes.  Everything else is source inspection,
# direct helper calls, or string assertions.
# =============================================================================

.dc_t0 <- proc.time()[["elapsed"]]

# ---------------------------------------------------------------------------
# Source-inspection helpers.  A route is "covered" when the guard call appears
# in the same function body as the derivation AND textually precedes it.
# ---------------------------------------------------------------------------

.dc_body_text <- function(name) {
  ns <- asNamespace("forestsearch")
  paste(deparse(body(get(name, envir = ns))), collapse = "\n")
}

.dc_before <- function(txt, first, second) {
  i <- regexpr(first,  txt, fixed = TRUE)
  j <- regexpr(second, txt, fixed = TRUE)
  i > 0L && j > 0L && i < j
}

# The routes to the m_diff derivation, established from source (Part 1):
#   1. .forestsearch_dina_select()  -- subgroup_method = "dina"
#   2. forestsearch()               -- use_dina + dina_args$selected_only
.DC_ROUTES <- list(
  list(fn = ".forestsearch_dina_select",
       derivation = 'm_diff <- if (identical(da$fit$family, "gaussian"))'),
  list(fn = "forestsearch",
       derivation = 'm_diff_sel <- if (identical(da$fit$family, "gaussian"))')
)


# =============================================================================
# Part 1 -- the guard
# =============================================================================

test_that("the guard helper exists and refuses only identity-scale, non-gaussian", {
  guard <- forestsearch:::.dina_assert_ratio_estimand

  # Fires: RD / IRD on every non-gaussian family DINA supports.
  for (fam in c("binomial", "poisson", "cox")) {
    for (em in c("RD", "IRD")) {
      expect_error(guard(fam, em),
                   'does not support identity-scale estimands',
                   fixed = TRUE)
    }
  }

  # Cannot fire on gaussian -- mddina is a committed campaign.
  expect_true(guard("gaussian", "MD"))
  expect_true(guard("gaussian", "RD"))
  expect_true(guard("gaussian", "IRD"))

  # Cannot fire for the ratio estimands, nor on survival (effect_measure NULL).
  for (fam in c("binomial", "poisson", "cox", "gaussian")) {
    for (em in list("OR", "RR", "IRR", "HR", "MD", NULL, NA_character_)) {
      expect_true(guard(fam, em))
    }
  }
})

test_that("the guard message is the one the directive specifies", {
  guard <- forestsearch:::.dina_assert_ratio_estimand
  msg <- tryCatch(guard("binomial", "RD"), error = conditionMessage)

  expect_true(grepl(
    'subgroup_method = "dina" does not support identity-scale estimands (RD, IRD):',
    msg, fixed = TRUE))
  expect_true(grepl(
    "DINA's admission floor operates on the family link scale.", msg,
    fixed = TRUE))
  expect_true(grepl('effect_measure = "OR"', msg, fixed = TRUE))
  expect_true(grepl('subgroup_method = "consistency" or', msg, fixed = TRUE))
  expect_true(grepl('"grf" for RD.', msg, fixed = TRUE))

  # The named measure follows the request.
  msg_ird <- tryCatch(guard("poisson", "IRD"), error = conditionMessage)
  expect_true(grepl('"grf" for IRD.', msg_ird, fixed = TRUE))

  # ASCII only (CRAN hygiene).
  expect_false(grepl("[^\x01-\x7f]", msg))
})

test_that("every route to the m_diff derivation is covered by the guard", {
  for (rt in .DC_ROUTES) {
    txt <- .dc_body_text(rt$fn)
    expect_true(grepl(rt$derivation, txt, fixed = TRUE),
                info = paste("derivation still present in", rt$fn))
    expect_true(grepl(".dina_assert_ratio_estimand(", txt, fixed = TRUE),
                info = paste("guard called in", rt$fn))
    expect_true(.dc_before(txt, ".dina_assert_ratio_estimand(", rt$derivation),
                info = paste("guard precedes the derivation in", rt$fn))
  }

  # No OTHER function in the package derives m_diff from hr.threshold: the two
  # routes above are the whole list.  (A third site would need its own guard.)
  ns   <- asNamespace("forestsearch")
  objs <- ls(ns, all.names = TRUE)
  fns  <- objs[vapply(objs, function(o) is.function(get(o, envir = ns)),
                      logical(1))]
  pat <- 'm_diff(_sel)? <- if \\(identical\\(da\\$fit\\$family, "gaussian"\\)\\)'
  derivers <- Filter(function(f) {
    txt <- paste(deparse(body(get(f, envir = ns))), collapse = "\n")
    grepl(pat, txt)
  }, fns)
  expect_setequal(derivers, vapply(.DC_ROUTES, `[[`, character(1), "fn"))
})

test_that("the guard fires before DINA's model runs", {
  skip_if_not_installed("testthat")
  df <- .make_binary_data(N = 120L, seed = 11L)

  # dina() replaced by a tripwire: if the guard let execution through, the
  # error would be "MODEL RAN", not the refusal.
  testthat::local_mocked_bindings(
    dina = function(...) stop("MODEL RAN", call. = FALSE),
    .package = "forestsearch")

  err <- tryCatch(
    forestsearch:::.forestsearch_dina_select(
      df = df, df.predict = NULL, df.test = NULL,
      confounders.name = c("age", "biomarker"),
      outcome.name = "y", event.name = "y", treat.name = "treat",
      id.name = "id", outcome_type = "binary",
      hr.threshold = 0.07, n.min = 30L, sg_focus = "maxSG",
      selection_rule = "hr", effect_neighborhood = 0.05,
      dina_args = list(), dina_res = NULL, seedit = 1L, details = FALSE,
      effect_measure = "RD", adverse_outcome = TRUE),
    error = conditionMessage)

  expect_false(grepl("MODEL RAN", err, fixed = TRUE))
  expect_true(grepl("does not support identity-scale estimands", err,
                    fixed = TRUE))

  # The same call with a ratio estimand does reach the fit (tripwire trips),
  # which is what proves the guard is the thing that stopped the RD call.
  err_or <- tryCatch(
    forestsearch:::.forestsearch_dina_select(
      df = df, df.predict = NULL, df.test = NULL,
      confounders.name = c("age", "biomarker"),
      outcome.name = "y", event.name = "y", treat.name = "treat",
      id.name = "id", outcome_type = "binary",
      hr.threshold = 1.25, n.min = 30L, sg_focus = "maxSG",
      selection_rule = "hr", effect_neighborhood = 0.05,
      dina_args = list(), dina_res = NULL, seedit = 1L, details = FALSE,
      effect_measure = "OR", adverse_outcome = TRUE),
    error = conditionMessage)
  expect_true(grepl("MODEL RAN", err_or, fixed = TRUE))

  # And the mddina shape (gaussian family, MD) also reaches the fit: the guard
  # is provably unable to fire there.
  err_md <- tryCatch(
    forestsearch:::.forestsearch_dina_select(
      df = df, df.predict = NULL, df.test = NULL,
      confounders.name = c("age", "biomarker"),
      outcome.name = "y", event.name = "y", treat.name = "treat",
      id.name = "id", outcome_type = "continuous",
      hr.threshold = 0.30, n.min = 30L, sg_focus = "maxSG",
      selection_rule = "hr", effect_neighborhood = 0.05,
      dina_args = list(), dina_res = NULL, seedit = 1L, details = FALSE,
      effect_measure = "MD", adverse_outcome = TRUE),
    error = conditionMessage)
  expect_true(grepl("MODEL RAN", err_md, fixed = TRUE))
})

test_that("forestsearch(subgroup_method = 'dina') refuses RD and IRD end to end", {
  df <- .make_binary_data(N = 120L, seed = 11L)
  args <- .fs_args_for("binary",
                       confounders = c("age", "biomarker"),
                       extra = list(subgroup_method = "dina",
                                    effect_measure = "RD",
                                    effect.threshold = 0.07,
                                    consistency.threshold = 0.05,
                                    use_grf = FALSE, use_lasso = FALSE,
                                    n.min = 30L))
  expect_error(do.call(forestsearch, c(list(df.analysis = df), args)),
               "does not support identity-scale estimands", fixed = TRUE)

  dfc <- .make_count_data(N = 120L, seed = 11L)
  argsc <- .fs_args_for("count",
                        confounders = c("age", "biomarker"),
                        extra = list(subgroup_method = "dina",
                                     effect_measure = "IRD",
                                     effect.threshold = 0.02,
                                     consistency.threshold = 0.01,
                                     use_grf = FALSE, use_lasso = FALSE,
                                     n.min = 30L))
  expect_error(do.call(forestsearch, c(list(df.analysis = dfc), argsc)),
               "does not support identity-scale estimands", fixed = TRUE)
})

test_that("the guard cannot fire under consistency or grf", {
  df <- .make_binary_data(N = 120L, seed = 11L)
  for (m in c("consistency", "grf")) {
    args <- .fs_args_for("binary",
                         confounders = c("age", "biomarker"),
                         extra = list(subgroup_method = m,
                                      effect_measure = "RD",
                                      effect.threshold = 0.07,
                                      consistency.threshold = 0.05,
                                      use_grf = (m == "grf"),
                                      use_lasso = FALSE,
                                      use_dina = FALSE, n.min = 30L))
    out <- tryCatch(
      suppressWarnings(suppressMessages(
        do.call(forestsearch, c(list(df.analysis = df), args)))),
      error = conditionMessage)
    if (is.character(out)) {
      expect_false(grepl("identity-scale estimands", out, fixed = TRUE),
                   info = paste("subgroup_method =", m))
    } else {
      expect_true(TRUE)
    }
  }
})

test_that("the use_dina screening route reports the refusal and contributes no cuts", {
  # The guard sits inside the screening tryCatch, so on this OPTIONAL route it
  # surfaces as the existing "DINA analysis failed" warning and the run
  # continues with the consistency search -- which is correct for RD.
  df <- .make_binary_data(N = 120L, seed = 11L)
  args <- .fs_args_for("binary",
                       confounders = c("age", "biomarker"),
                       extra = list(subgroup_method = "consistency",
                                    effect_measure = "RD",
                                    effect.threshold = 0.07,
                                    consistency.threshold = 0.05,
                                    use_dina = TRUE, use_grf = FALSE,
                                    use_lasso = FALSE, n.min = 30L))
  cap <- .run_fs_capture(df, args)
  expect_true(any(grepl("does not support identity-scale estimands",
                        cap$warnings, fixed = TRUE)))
})


# =============================================================================
# Part 2 -- the frontier-key warning under subgroup_method = "dina"
# =============================================================================

.dc_dina_select_warnings <- function(dina_args, effect_measure = "OR",
                                     outcome_type = "binary") {
  df <- .make_binary_data(N = 120L, seed = 11L)
  seen <- character(0)
  testthat::local_mocked_bindings(
    dina = function(...) stop("MODEL RAN", call. = FALSE),
    .package = "forestsearch")
  withCallingHandlers(
    tryCatch(
      forestsearch:::.forestsearch_dina_select(
        df = df, df.predict = NULL, df.test = NULL,
        confounders.name = c("age", "biomarker"),
        outcome.name = "y", event.name = "y", treat.name = "treat",
        id.name = "id", outcome_type = outcome_type,
        hr.threshold = 1.25, n.min = 30L, sg_focus = "maxSG",
        selection_rule = "hr", effect_neighborhood = 0.05,
        dina_args = dina_args, dina_res = NULL, seedit = 1L, details = FALSE,
        effect_measure = effect_measure, adverse_outcome = TRUE),
      error = function(e) NULL),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  seen
}

test_that("the seven frontier keys are the documented set", {
  expect_setequal(forestsearch:::.DINA_FRONTIER_KEYS,
                  c("scope", "m_diff", "n_min", "direction",
                    "max_per_covariate", "max_subgroups", "digits"))
})

test_that("frontier keys under dina warn ONCE, naming every offending key", {
  keys <- forestsearch:::.DINA_FRONTIER_KEYS
  supplied <- list(scope = "wide", m_diff = 0.2, n_min = 30L,
                   direction = "both", max_per_covariate = 3L,
                   max_subgroups = 10L, digits = 3L)

  w <- .dc_dina_select_warnings(supplied)
  hits <- grep("frontier key", w, fixed = TRUE, value = TRUE)
  expect_length(hits, 1L)             # one warning per fit, never one per key
  for (k in keys) {
    expect_true(grepl(shQuote(k), hits[[1L]], fixed = TRUE), info = k)
  }
  expect_true(grepl('subgroup_method = "dina"', hits[[1L]], fixed = TRUE))

  # A subset names exactly that subset.
  w2 <- .dc_dina_select_warnings(list(max_subgroups = 5L, digits = 2L))
  hits2 <- grep("frontier key", w2, fixed = TRUE, value = TRUE)
  expect_length(hits2, 1L)
  expect_true(grepl(shQuote("max_subgroups"), hits2[[1L]], fixed = TRUE))
  expect_true(grepl(shQuote("digits"), hits2[[1L]], fixed = TRUE))
  expect_false(grepl(shQuote("scope"), hits2[[1L]], fixed = TRUE))
})

test_that("no frontier keys, no warning; fit / behaviour keys do not warn", {
  expect_length(grep("frontier key", .dc_dina_select_warnings(list()),
                     fixed = TRUE), 0L)
  expect_length(
    grep("frontier key",
         .dc_dina_select_warnings(list(family = "binomial", seed = 3L,
                                       max_depth = 1L, selected_only = TRUE)),
         fixed = TRUE),
    0L)
})

test_that("the same dina_args under consistency and grf warn nothing", {
  df <- .make_binary_data(N = 120L, seed = 11L)
  for (m in c("consistency", "grf")) {
    args <- .fs_args_for("binary",
                         confounders = c("age", "biomarker"),
                         extra = list(subgroup_method = m,
                                      use_grf = (m == "grf"),
                                      use_lasso = FALSE, use_dina = FALSE,
                                      n.min = 30L,
                                      dina_args = list(scope = "wide",
                                                       max_subgroups = 5L,
                                                       digits = 2L)))
    cap <- .run_fs_capture(df, args)
    expect_length(grep("frontier key", cap$warnings, fixed = TRUE), 0L)
  }
})


test_that("wall clock stays inside the file's abort budget", {
  elapsed <- proc.time()[["elapsed"]] - .dc_t0
  message(sprintf("[directive C] acceptance-test wall clock: %.1f s", elapsed))
  expect_lt(elapsed, 180)
})
