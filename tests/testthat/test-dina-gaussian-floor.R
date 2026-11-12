# =============================================================================
# The continuous-outcome (gaussian family) DINA floor is the resolved
# screening threshold (audit A-001).
#
# Two derivation sites feed dina_subgroup() its floor m_diff:
#   1. .forestsearch_dina_select()        -- subgroup_method = "dina"
#   2. forestsearch(), m_diff_sel         -- use_dina = TRUE, selected_only
# With effect.threshold and hr.threshold both unset, MD resolves the
# screening threshold to 0 and records 0 (threshold_config$screening,
# admission$effect_floor, args_call_all$effect.threshold); a bootstrap / CV
# replay of args_call_all therefore runs with 0.  The fit must use the same 0,
# not the raw hr.threshold default 1.25.
#
# The floor is observed by wrapping dina_subgroup() (the only consumer of
# m_diff at both sites) and recording the m_diff each call receives.
# Compute: a handful of small (N = 200) sequential micro-fits.
# =============================================================================

skip_on_cran()

.gf_df <- function() .make_continuous_data(N = 200L, MD_harm = 1.5, seed = 42L)

# Run `f()` with dina_subgroup() wrapped; return its value and the m_diff of
# every dina_subgroup() call made while it ran.
.gf_capture <- function(f) {
  floors <- numeric(0)
  orig <- get("dina_subgroup", envir = asNamespace("forestsearch"))
  testthat::local_mocked_bindings(
    dina_subgroup = function(...) {
      floors <<- c(floors, list(...)$m_diff)
      orig(...)
    },
    .package = "forestsearch")
  value <- suppressWarnings(suppressMessages(f()))
  list(value = value, floors = floors)
}

# forestsearch() arguments; `thresholds` is merged last, so NULL entries drop
# the argument (it is then genuinely missing() in the call).
.gf_args <- function(route, thresholds = list()) {
  route_args <- switch(route,
    dina = list(subgroup_method = "dina", use_dina = FALSE),
    use_dina = list(subgroup_method = "consistency", use_dina = TRUE))
  args <- .fs_args_for("continuous",
                       confounders = c("age", "biomarker"),
                       extra = c(route_args,
                                 list(use_grf = FALSE, use_lasso = FALSE,
                                      n.min = 40L, hr.threshold = NULL,
                                      hr.consistency = NULL)))
  modifyList(args, thresholds)
}

.gf_fit <- function(route, thresholds = list()) {
  df <- .gf_df()
  args <- .gf_args(route, thresholds)
  .gf_capture(function() do.call(forestsearch, c(list(df.analysis = df), args)))
}

# Replay the parent's args_call_all the way the bootstrap does
# (R/bootstrap_analysis_dofuture.R: args_FS_boot <- fs.est$args_call_all, the
# DINA / GRF inputs nulled, output suppressed, sequential plan).
.gf_replay <- function(fs) {
  a <- fs$args_call_all
  a$df.analysis <- .gf_df()
  a$grf_res <- NULL
  a$grf_cuts <- NULL
  a$dina_res <- NULL
  a$dina_cuts <- NULL
  a$details <- FALSE
  a$quiet <- TRUE
  a$plot.sg <- FALSE
  a$plot.grf <- FALSE
  a$parallel_args$plan <- "sequential"
  a$parallel_args$workers <- 1L
  a$parallel_args$show_message <- FALSE
  .gf_capture(function() do.call(forestsearch, a))
}


test_that("both thresholds unset, subgroup_method = 'dina': fit, record and replay floors agree", {
  fit <- .gf_fit("dina")
  fs <- fit$value

  expect_length(fit$floors, 1L)
  expect_identical(fs$threshold_config$screening, 0)
  expect_identical(fs$args_call_all$effect.threshold, 0)
  expect_identical(fs$admission$effect_floor, 0)
  expect_identical(fit$floors, fs$threshold_config$screening)

  rep <- .gf_replay(fs)
  expect_length(rep$floors, 1L)
  expect_identical(rep$floors, fit$floors)
})

test_that("both thresholds unset, use_dina selected-cut screening: fit, record and replay floors agree", {
  fit <- .gf_fit("use_dina")
  fs <- fit$value

  expect_length(fit$floors, 1L)
  expect_identical(fs$threshold_config$screening, 0)
  expect_identical(fs$args_call_all$effect.threshold, 0)
  expect_identical(fit$floors, fs$threshold_config$screening)

  rep <- .gf_replay(fs)
  expect_length(rep$floors, 1L)
  expect_identical(rep$floors, fit$floors)
})

test_that("a supplied threshold is the floor at both sites (either spelling)", {
  for (route in c("dina", "use_dina")) {
    for (spelling in c("effect.threshold", "hr.threshold")) {
      thr <- stats::setNames(list(0.4), spelling)
      fit <- .gf_fit(route, thr)
      info <- paste(route, spelling)
      expect_identical(fit$floors, 0.4, info = info)
      expect_identical(fit$value$threshold_config$screening, 0.4, info = info)
      expect_identical(fit$value$args_call_all$effect.threshold, 0.4,
                       info = info)
    }
  }
})

test_that("a supplied ratio-scale threshold that is remapped for RD is still the floor as given", {
  # Binary RD with the gaussian family forced through dina_args: a supplied
  # effect.threshold = 1.5 is remapped (with a warning) to the recorded 0.05,
  # but a supplied threshold leaves the floor unchanged from before A-001.
  set.seed(7); N <- 300L
  df <- data.frame(id = seq_len(N), treat = rbinom(N, 1L, 0.5),
                   x1 = round(rnorm(N), 3), x2 = round(runif(N, 0, 10), 2))
  df$y <- rbinom(N, 1L, plogis(-0.5 + 0.8 * df$treat * (df$x1 <= 0)))
  for (route in c("dina", "use_dina")) {
    route_args <- switch(route,
      dina = list(subgroup_method = "dina", use_dina = FALSE),
      use_dina = list(subgroup_method = "consistency", use_dina = TRUE))
    fit <- .gf_capture(function() do.call(forestsearch, c(list(
      df.analysis = df, confounders.name = c("x1", "x2"), outcome.name = "y",
      event.name = "y", treat.name = "treat", id.name = "id",
      outcome_type = "binary", effect_measure = "RD", adverse_outcome = TRUE,
      dina_args = list(family = "gaussian"), effect.threshold = 1.5,
      sg_focus = "maxSG", n.min = 40L, use_grf = FALSE, use_lasso = FALSE,
      details = FALSE, quiet = TRUE, plot.sg = FALSE, plot.grf = FALSE,
      parallel_args = list(plan = "sequential", workers = 1L,
                           show_message = FALSE)), route_args)))
    expect_identical(unique(fit$floors), 1.5, info = route)
    expect_identical(fit$value$args_call_all$effect.threshold, 0.05, info = route)
  }
})

test_that("the helper's gaussian floor defaults to hr.threshold for direct callers", {
  # .forestsearch_dina_select() called without effect.threshold (as
  # quarto/dina/smoke-tests/validate_screen3.R does) keeps hr.threshold.
  df <- .gf_df()
  sel <- function(...) {
    .gf_capture(function() .forestsearch_dina_select(
      df = df, df.predict = NULL, df.test = NULL,
      confounders.name = c("age", "biomarker"), outcome.name = "y",
      event.name = NULL, treat.name = "treat", id.name = "id",
      outcome_type = "continuous", hr.threshold = 0.3, n.min = 40L,
      sg_focus = "maxSG", selection_rule = "neighborhood",
      effect_neighborhood = 0.10, dina_args = list(), dina_res = NULL,
      seedit = 1L, details = FALSE, effect_measure = "MD",
      adverse_outcome = FALSE, ...))
  }
  expect_identical(sel()$floors, 0.3)
  expect_identical(sel(effect.threshold = 0)$floors, 0)
})
