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


# =============================================================================
# Part 3 -- honest display caps
#
# The file's ONE micro-fit lives here and is reused by every cap assertion.
# =============================================================================

.dc_fit_env <- new.env(parent = emptyenv())

.dc_micro_fit <- function() {
  if (!is.null(.dc_fit_env$fit)) return(.dc_fit_env)
  set.seed(7L)
  n  <- 200L
  df <- data.frame(
    w  = stats::rbinom(n, 1L, 0.5),
    x1 = stats::runif(n, -1, 1),
    x2 = stats::runif(n, -1, 1),
    x3 = stats::runif(n, -1, 1))
  tau  <- 0.3 + 1.2 * df$x1 - 0.4 * df$x2
  df$y <- 0.5 * df$x1 + df$w * tau + stats::rnorm(n)
  t0 <- proc.time()[["elapsed"]]
  .dc_fit_env$fit  <- dina(df, outcome = "y", treatment = "w",
                           covariates = c("x1", "x2", "x3"),
                           family = "gaussian", seed = 1L)
  .dc_fit_env$secs <- proc.time()[["elapsed"]] - t0
  .dc_fit_env$df   <- df
  .dc_fit_env$cov  <- c("x1", "x2", "x3")
  .dc_fit_env
}

test_that("dina_frontier() caps default to Inf / Inf", {
  fm <- formals(dina_frontier)
  expect_identical(eval(fm$max_per_covariate), Inf)
  expect_identical(eval(fm$max_subgroups),     Inf)
})

test_that("the default frontier is untrimmed and silent", {
  e <- .dc_micro_fit()
  expect_silent(
    fr_inf <- dina_frontier(e$fit, e$df, covariates = e$cov, n_min = 40L))
  expect_gt(nrow(fr_inf), 0L)
  # Inf shows every non-dominated cut the extractor found.
  expect_identical(nrow(fr_inf), as.integer(attr(fr_inf, "n_frontier")))
  .dc_fit_env$fr_inf <- fr_inf
})

test_that("a finite cap that trims warns, naming kept and available", {
  e  <- .dc_micro_fit()
  n_full <- nrow(.dc_fit_env$fr_inf)
  skip_if(n_full < 3L, "frontier too small to exercise a trimming cap")

  cap <- n_full - 1L
  w <- tryCatch(dina_frontier(e$fit, e$df, covariates = e$cov, n_min = 40L,
                              max_subgroups = cap),
                warning = conditionMessage)
  expect_true(is.character(w))
  expect_true(grepl(paste0("max_subgroups = ", cap), w, fixed = TRUE))
  expect_true(grepl(paste0(cap, " of ", n_full), w, fixed = TRUE))
  expect_true(grepl("DISPLAY limit", w, fixed = TRUE))

  # ...and it really trims.
  fr <- suppressWarnings(
    dina_frontier(e$fit, e$df, covariates = e$cov, n_min = 40L,
                  max_subgroups = cap))
  expect_identical(nrow(fr), as.integer(cap))

  # The condition carries the class the screening path muffles.
  cond <- tryCatch(dina_frontier(e$fit, e$df, covariates = e$cov, n_min = 40L,
                                 max_subgroups = cap),
                   warning = function(w) w)
  expect_s3_class(cond, "dina_frontier_cap_trim")
})

test_that("a finite cap that does not trim is silent", {
  e <- .dc_micro_fit()
  n_full <- nrow(.dc_fit_env$fr_inf)
  expect_silent(dina_frontier(e$fit, e$df, covariates = e$cov, n_min = 40L,
                              max_subgroups = n_full,
                              max_per_covariate = n_full))
})

test_that("max_per_covariate trims and warns on its own terms", {
  e  <- .dc_micro_fit()
  per <- table(.dc_fit_env$fr_inf$covariate)
  skip_if(max(per) < 2L, "no covariate contributes 2+ cuts")
  cap <- max(per) - 1L
  w <- tryCatch(dina_frontier(e$fit, e$df, covariates = e$cov, n_min = 40L,
                              max_per_covariate = cap),
                warning = conditionMessage)
  expect_true(is.character(w))
  expect_true(grepl(paste0("max_per_covariate = ", cap), w, fixed = TRUE))
})

test_that("the discarded screening frontier muffles exactly the trim warning", {
  txt <- .dc_body_text("forestsearch")
  expect_true(grepl("dina_frontier_cap_trim = function(w)", txt, fixed = TRUE))
  expect_true(grepl("isTRUE(da$selected_only)", txt, fixed = TRUE))
  expect_true(grepl("muffleWarning", txt, fixed = TRUE))
  # The screening pool is untouched: .resolve_dina_args() still supplies the
  # finite caps, so Part 3 changes no candidate the search evaluates.
  da <- forestsearch:::.resolve_dina_args(list(), "binary",
                                          n_min_default = 60L, seed_default = 1L)
  expect_identical(da$frontier$max_per_covariate, 3L)
  expect_identical(da$frontier$max_subgroups, 10L)
})

test_that("the micro-fit stays inside the compute cap", {
  e <- .dc_micro_fit()
  message(sprintf("[directive C] micro-fit (1 x dina, gaussian, n = 200): %.2f s",
                  e$secs))
  expect_lt(e$secs, 300)
})


# =============================================================================
# Part 4 -- the details-time frontier print retitle
#
# Reuses the Part 3 micro-fit via `dina_res`, so no second model is fitted.
# =============================================================================

.dc_details_lines <- function() {
  e <- .dc_micro_fit()
  msgs <- character(0)
  withCallingHandlers(
    tryCatch(
      forestsearch:::.forestsearch_dina_select(
        df = e$df, df.predict = NULL, df.test = NULL,
        confounders.name = e$cov,
        outcome.name = "y", event.name = "y", treat.name = "w",
        id.name = "id", outcome_type = "continuous",
        hr.threshold = 0.10, n.min = 40L, sg_focus = "maxSG",
        selection_rule = "neighborhood", effect_neighborhood = 0.05,
        dina_args = list(family = "gaussian"), dina_res = e$fit,
        seedit = 1L, details = TRUE,
        effect_measure = "MD", adverse_outcome = TRUE),
      error = function(err) NULL),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) invokeRestart("muffleWarning"))
  paste(msgs, collapse = "")
}

test_that("the frontier print is titled as proposed single cuts", {
  out <- .dc_details_lines()

  expect_true(grepl("DINA frontier -- proposed single cuts", out, fixed = TRUE))
  expect_true(grepl("display only, not the searched family", out, fixed = TRUE))

  # The old title is gone.
  expect_false(grepl("DINA frontier candidates", out, fixed = TRUE))

  # Shown BESIDE the family counts, which are unchanged.
  expect_true(grepl("Candidates searched:", out, fixed = TRUE))
  expect_true(grepl("Candidates qualifying", out, fixed = TRUE))

  # Structure is unchanged: the same header block still precedes it.
  expect_true(grepl("[forestsearch] DINA selection", out, fixed = TRUE))
  expect_true(grepl("Harm floor:", out, fixed = TRUE))
})

test_that("neither old frontier title survives anywhere in the package", {
  ns  <- asNamespace("forestsearch")
  txt <- vapply(ls(ns, all.names = TRUE), function(o) {
    f <- get(o, envir = ns)
    if (!is.function(f)) return("")
    paste(deparse(body(f)), collapse = "\n")
  }, character(1))
  all_txt <- paste(txt, collapse = "\n")
  expect_false(grepl("DINA frontier candidates (per-covariate non-dominated)",
                     all_txt, fixed = TRUE))
  expect_false(grepl("DINA frontier: no candidates met the size constraint.",
                     all_txt, fixed = TRUE))
})


# =============================================================================
# Part 5 -- the c2 / p* echo annotation under dina / grf
#
# Display only.  No value changes, no warning, no error: the echo sites are
# found from source and each is exercised directly.
# =============================================================================

# The echo sites, established from source (Part 5).  Each is a place c2
# (consistency.threshold / hr.consistency) or p* (pconsistency.threshold) is
# printed back to the user.
.DC_ECHO_SITES <- c(
  "forestsearch",              # the two config banners (GLM + survival)
  "interpret_search_config",   # the Search Alignment Diagnostic thresholds
  "summary.forestsearch",      # Analysis Parameters
  "print_cv_params",           # ForestSearch parameters for CV folds
  "fs_family_report"           # the consistency-screen row
)

test_that("the annotation helper fires only on dina and grf", {
  note <- forestsearch:::.fs_c2_inert_note
  expect_identical(note("consistency"), "")
  expect_identical(note(NULL), "")
  expect_identical(note(NA_character_), "")
  expect_true(grepl('subgroup_method = "dina"', note("dina"), fixed = TRUE))
  expect_true(grepl('subgroup_method = "grf"', note("grf"), fixed = TRUE))
  expect_true(grepl("not used on this path", note("grf"), fixed = TRUE))
  expect_false(grepl("[^-]", note("grf")))
})

test_that("every echo site found in Part 5 carries the annotation", {
  for (fn in .DC_ECHO_SITES) {
    txt <- .dc_body_text(fn)
    hit <- grepl(".fs_c2_inert_note(", txt, fixed = TRUE) ||
           grepl("not used on this path", txt, fixed = TRUE)
    expect_true(hit, info = paste("echo site annotated:", fn))
  }
  # The bootstrap banner too (its else-branch; dina is already excluded there).
  expect_true(grepl(".fs_c2_inert_note(",
                    .dc_body_text("forestsearch_bootstrap_dofuture"),
                    fixed = TRUE))
})

test_that("the banner annotates c2 and p* under grf and not under consistency", {
  df <- .make_binary_data(N = 120L, seed = 11L)
  grab <- function(m) {
    args <- .fs_args_for("binary",
                         confounders = c("age", "biomarker"),
                         extra = list(subgroup_method = m,
                                      use_grf = (m == "grf"),
                                      use_lasso = FALSE, use_dina = FALSE,
                                      n.min = 30L, quiet = FALSE,
                                      details = FALSE))
    msgs <- character(0)
    withCallingHandlers(
      tryCatch(suppressWarnings(
        do.call(forestsearch, c(list(df.analysis = df), args))),
        error = function(e) NULL),
      message = function(mm) {
        msgs <<- c(msgs, conditionMessage(mm))
        invokeRestart("muffleMessage")
      })
    paste(msgs, collapse = "")
  }

  out_grf <- grab("grf")
  expect_true(grepl("Subgroup Identification Configuration", out_grf,
                    fixed = TRUE))
  expect_true(grepl('not used on this path: subgroup_method = "grf"',
                    out_grf, fixed = TRUE))

  out_cons <- grab("consistency")
  expect_true(grepl("Subgroup Identification Configuration", out_cons,
                    fixed = TRUE))
  expect_false(grepl("not used on this path", out_cons, fixed = TRUE))

  # The echoed VALUES are untouched -- annotation only.
  expect_true(grepl("Consistency rate threshold: 80%", out_grf, fixed = TRUE))
  expect_true(grepl("Consistency rate threshold: 80%", out_cons, fixed = TRUE))
})

test_that("interpret_search_config annotates only dina / grf", {
  cap <- function(m) {
    msgs <- character(0)
    withCallingHandlers(
      interpret_search_config(
        outcome_type = "binary", effect_measure = "OR",
        adverse_outcome = TRUE, effect_threshold = log(1.25),
        consistency_threshold = log(1.0),
        use_lasso = FALSE, use_grf = TRUE,
        outcome.name = "y", event.name = "y", treat.name = "treat",
        subgroup_method = m, quiet = FALSE),
      message = function(mm) {
        msgs <<- c(msgs, conditionMessage(mm))
        invokeRestart("muffleMessage")
      })
    paste(msgs, collapse = "")
  }
  expect_true(grepl('not used on this path: subgroup_method = "grf"',
                    cap("grf"), fixed = TRUE))
  expect_true(grepl('not used on this path: subgroup_method = "dina"',
                    cap("dina"), fixed = TRUE))
  expect_false(grepl("not used on this path", cap("consistency"),
                     fixed = TRUE))
  # Default keeps the pre-Part-5 output for every existing caller.
  expect_identical(eval(formals(interpret_search_config)$subgroup_method),
                   "consistency")
})

test_that("print_cv_params annotates only dina / grf", {
  base <- list(sg_focus = "maxSG", maxk = 2L, fs.splits = 100L,
               max_subgroups_search = Inf, hr.threshold = 1.25,
               hr.consistency = 1.0, pconsistency.threshold = 0.8,
               n.min = 40L, use_twostage = FALSE, use_lasso = FALSE,
               use_grf = TRUE, outcome_type = "survival")
  grab <- function(m) paste(utils::capture.output(
    forestsearch:::print_cv_params(utils::modifyList(
      base, list(subgroup_method = m)))), collapse = "\n")
  expect_true(grepl("not used on this path", grab("grf"), fixed = TRUE))
  expect_true(grepl("not used on this path", grab("dina"), fixed = TRUE))
  expect_false(grepl("not used on this path", grab("consistency"),
                     fixed = TRUE))
})

test_that("summary.forestsearch annotates only dina / grf", {
  stub <- function(m) structure(
    list(args_call_all = list(
      sg_focus = "maxSG", hr.threshold = 1.25, hr.consistency = 1.0,
      pconsistency.threshold = 0.8, n.min = 40, fs.splits = 100, maxk = 2,
      use_twostage = FALSE, use_lasso = FALSE, use_grf = TRUE,
      use_dina = FALSE, subgroup_method = m)),
    class = "forestsearch")
  grab <- function(m) paste(utils::capture.output(summary(stub(m))),
                            collapse = "\n")
  for (m in c("dina", "grf")) {
    out <- grab(m)
    expect_true(grepl(sprintf('not used on this path: subgroup_method = "%s"', m),
                      out, fixed = TRUE), info = m)
    # Both quantities annotated, and the values themselves unchanged.
    expect_true(grepl("hr.consistency: 1  [not used", out, fixed = TRUE))
    expect_true(grepl("pconsistency.threshold: 0.8  [not used", out,
                      fixed = TRUE))
  }
  out_c <- grab("consistency")
  expect_false(grepl("not used on this path", out_c, fixed = TRUE))
  expect_true(grepl("hr.consistency: 1", out_c, fixed = TRUE))
  expect_true(grepl("pconsistency.threshold: 0.8", out_c, fixed = TRUE))
})

test_that("Part 5 changes no resolved value on any path", {
  # Gate C's instrument: the threshold-resolution probe, all 175 cells,
  # compared to the baseline recorded before the first Directive C edit.
  base <- utils::read.csv("../../dev/reports/baseline_directive_C_2026-09-18.csv",
                          stringsAsFactors = FALSE, colClasses = "character")
  now  <- probe_threshold_sync(sync = TRUE, validate = TRUE, extended = TRUE)
  now  <- as.data.frame(lapply(now, as.character), stringsAsFactors = FALSE)
  base <- as.data.frame(lapply(base, as.character), stringsAsFactors = FALSE)
  base[is.na(base)] <- ""
  now[is.na(now)]   <- ""
  expect_identical(dim(now), dim(base))
  expect_identical(names(now), names(base))
  for (j in names(base)) {
    expect_identical(now[[j]], base[[j]], info = paste("probe column", j))
  }
})


# =============================================================================
# Gates A / B / C
# =============================================================================

test_that("Gate B: an mddina-shaped call passes silently", {
  # gaussian family, MD, subgroup_method = "dina" -- a committed campaign
  # shape.  The SECOND and last micro-fit in this file.
  df <- .make_continuous_data(N = 120L, seed = 5L)
  args <- .fs_args_for("continuous",
                       confounders = c("age", "biomarker"),
                       extra = list(subgroup_method = "dina",
                                    effect_measure = "MD",
                                    effect.threshold = 0.30,
                                    consistency.threshold = 0.20,
                                    use_grf = FALSE, use_lasso = FALSE,
                                    n.min = 30L))
  t0 <- proc.time()[["elapsed"]]
  cap <- .run_fs_capture(df, args)
  secs <- proc.time()[["elapsed"]] - t0
  message(sprintf("[directive C] Gate B micro-fit (forestsearch + dina, MD): %.2f s",
                  secs))
  expect_lt(secs, 300)

  # No Directive C condition fires on this shape.
  expect_length(grep("identity-scale estimands", cap$warnings, fixed = TRUE), 0L)
  expect_length(grep("frontier key", cap$warnings, fixed = TRUE), 0L)
  expect_length(grep("trimmed the frontier display", cap$warnings,
                     fixed = TRUE), 0L)
})

test_that("Gate A: the refusal is the only new failure mode", {
  guard <- forestsearch:::.dina_assert_ratio_estimand
  # The full (family, measure) grid: the guard fires on exactly six cells.
  fams <- c("cox", "binomial", "poisson", "gaussian")
  ems  <- c("HR", "OR", "RR", "IRR", "MD", "RD", "IRD")
  fires <- outer(fams, ems, Vectorize(function(f, m)
    inherits(tryCatch(guard(f, m), error = function(e) e), "error")))
  dimnames(fires) <- list(fams, ems)
  expect_identical(sum(fires), 6L)
  expect_true(all(fires[c("cox", "binomial", "poisson"), c("RD", "IRD")]))
  expect_false(any(fires["gaussian", ]))
  expect_false(any(fires[, c("HR", "OR", "RR", "IRR", "MD")]))
})

test_that("Gate A: neither new warning fires on a clean call", {
  # A plain consistency run touches none of the three conditions.
  df <- .make_binary_data(N = 120L, seed = 11L)
  args <- .fs_args_for("binary", confounders = c("age", "biomarker"),
                       extra = list(use_grf = FALSE, use_lasso = FALSE,
                                    use_dina = FALSE, n.min = 30L))
  cap <- .run_fs_capture(df, args)
  for (pat in c("identity-scale estimands", "frontier key",
                "trimmed the frontier display")) {
    expect_length(grep(pat, cap$warnings, fixed = TRUE), 0L)
  }
})


test_that("wall clock stays inside the file's abort budget", {
  elapsed <- proc.time()[["elapsed"]] - .dc_t0
  message(sprintf("[directive C] acceptance-test wall clock: %.1f s", elapsed))
  expect_lt(elapsed, 180)
})
