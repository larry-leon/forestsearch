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
# DINA / GRF inputs nulled, output suppressed, sequential plan).  `df` is the
# parent's analysis data (default: the continuous fixture above).
.gf_replay <- function(fs, df = .gf_df()) {
  a <- fs$args_call_all
  a$df.analysis <- df
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

test_that("a supplied ratio-scale threshold that is remapped for RD is the floor as remapped", {
  # Binary RD with the gaussian family forced through dina_args: a supplied
  # effect.threshold = 1.5 is remapped (with a warning) to the recorded 0.05.
  # The gaussian floor is the resolved effect.threshold (?forestsearch,
  # dina_args: "identity for gaussian -- from the resolved effect.threshold"),
  # so the fit's floor is that 0.05 too (F14; before F14 it was the raw 1.5).
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
    # 0.05: the RD row of ?forestsearch "Resolved defaults, per estimand"
    # (c1 = 0.05), the value the remap warning names ("Remapping to
    # default: 0.05").
    expect_identical(unique(fit$floors), 0.05, info = route)
    expect_identical(fit$value$args_call_all$effect.threshold, 0.05, info = route)
  }
})

test_that("forced-gaussian corner (binary RD, count IRD): fit, record and replay floors agree", {
  # Setting: an identity-scale estimand (binary RD, count IRD) with
  # dina_args = list(family = "gaussian") forced and a supplied ratio-scale
  # threshold (effect.threshold = 1.5 > 1) that forestsearch() warns about and
  # remaps.  Expected recorded values, from ?forestsearch "Resolved defaults,
  # per estimand": RD -> c1 = 0.05, IRD -> c1 = 0.01 (identity scale, so
  # threshold_config$screening is that value itself, not its log).  The
  # gaussian DINA floor is "identity for gaussian -- from the resolved
  # effect.threshold" (?forestsearch, dina_args), so the floor the fit hands
  # dina_subgroup() must equal the recorded value, and a replay of
  # args_call_all (as the bootstrap / CV do) must use the same floor.
  set.seed(7); N <- 300L
  df_bin <- data.frame(id = seq_len(N), treat = rbinom(N, 1L, 0.5),
                       x1 = round(rnorm(N), 3), x2 = round(runif(N, 0, 10), 2))
  df_bin$y <- rbinom(N, 1L, plogis(-0.5 + 0.8 * df_bin$treat * (df_bin$x1 <= 0)))
  df_cnt <- .make_count_data(N = 300L, IRR_harm = 2.0, seed = 7L)

  settings <- list(
    RD = list(df = df_bin, expected = 0.05,
              args = list(confounders.name = c("x1", "x2"),
                          outcome.name = "y", event.name = "y",
                          treat.name = "treat", id.name = "id",
                          outcome_type = "binary", effect_measure = "RD",
                          adverse_outcome = TRUE)),
    IRD = list(df = df_cnt, expected = 0.01,
               args = list(confounders.name = c("age", "biomarker"),
                           outcome.name = "y", treat.name = "treat",
                           id.name = "id", offset.name = "ftime",
                           outcome_type = "count", effect_measure = "IRD",
                           adverse_outcome = TRUE)))

  for (em in names(settings)) {
    s <- settings[[em]]
    for (route in c("dina", "use_dina")) {
      route_args <- switch(route,
        dina = list(subgroup_method = "dina", use_dina = FALSE),
        use_dina = list(subgroup_method = "consistency", use_dina = TRUE))
      info <- paste(em, route)
      fit <- .gf_capture(function() do.call(forestsearch, c(
        list(df.analysis = s$df), s$args, list(
        dina_args = list(family = "gaussian"), effect.threshold = 1.5,
        sg_focus = "maxSG", n.min = 40L, use_grf = FALSE, use_lasso = FALSE,
        max_subgroups_search = 5L, fs.splits = 100L, seedit = 42L,
        details = FALSE, quiet = TRUE, plot.sg = FALSE, plot.grf = FALSE,
        parallel_args = list(plan = "sequential", workers = 1L,
                             show_message = FALSE)), route_args)))
      fs <- fit$value

      # One dina_subgroup() call per fit (one derivation site per route).
      expect_length(fit$floors, 1L)
      # The recorded values (documented resolved default for the estimand).
      expect_identical(fs$args_call_all$effect.threshold, s$expected, info = info)
      expect_identical(fs$threshold_config$screening, s$expected, info = info)
      # The fit's floor equals the recorded value.
      expect_identical(fit$floors, fs$args_call_all$effect.threshold, info = info)

      # The replay (bootstrap / CV) uses the same floor as the fit.
      rep <- .gf_replay(fs, df = s$df)
      expect_length(rep$floors, 1L)
      expect_identical(rep$floors, fit$floors, info = info)
    }
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
