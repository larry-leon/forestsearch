# =============================================================================
# Integration tests: GLM pipeline completeness
#
# These tests verify that forestsearch() with each GLM outcome type
# runs the FULL pipeline (GRF -> get_FSdata -> search -> consistency)
# without silent failures. They do NOT test statistical power or
# correctness of subgroup identification; they test infrastructure.
# =============================================================================

test_that("binary pipeline completes without internal failure", {
  set.seed(42)
  n <- 300
  df <- data.frame(
    id         = seq_len(n),
    treat      = rbinom(n, 1, 0.5),
    bm1        = as.factor(rbinom(n, 1, 0.70)),
    age        = round(rnorm(n, 55, 10)),
    progressed = rbinom(n, 1, 0.3)
  )
  # Strong planted signal so pipeline exercises all stages
  df$progressed <- ifelse(
    df$treat == 1 & df$bm1 == 0,
    rbinom(sum(df$treat == 1 & df$bm1 == 0), 1, 0.7),
    df$progressed
  )

  result <- forestsearch(
    df.analysis      = df,
    confounders.name = c("bm1", "age"),
    outcome.name     = "progressed",
    event.name       = "progressed",
    treat.name       = "treat",
    id.name          = "id",
    outcome_type     = "binary",
    effect_measure   = "OR",
    adverse_outcome  = TRUE,
    hr.threshold     = 1.10,
    hr.consistency   = 1.0,
    pconsistency.threshold = 0.30,
    fs.splits        = 20,
    n.min            = 20,
    maxk             = 1,
    use_lasso        = TRUE,
    use_grf          = TRUE,
    is.RCT           = TRUE,
    seedit           = 42,
    details          = FALSE,
    plot.sg          = FALSE,
    parallel_args    = list(plan = "sequential", workers = 1,
                            show_message = FALSE)
  )

  # Pipeline should complete and return full structure
  expect_true(is.list(result))
  expect_true(length(names(result)) > 1,
              info = "Result has only sg.harm -- pipeline failed silently")

  # If error_log exists, report it
  if (!is.null(result$error_log)) {
    fail(sprintf("Pipeline failed at %s: %s",
                 result$error_log$stage, result$error_log$message))
  }
})


test_that("continuous pipeline completes without internal failure", {
  set.seed(42)
  n <- 400
  df <- data.frame(
    id       = seq_len(n),
    treat    = rbinom(n, 1, 0.5),
    sero     = as.factor(rbinom(n, 1, 0.65)),
    age      = round(rnorm(n, 55, 12)),
    das28_chg = rnorm(n, -1.5, 1.5)
  )
  # Planted harm: sero=0 treated patients get worse
  idx_harm <- df$treat == 1 & df$sero == 0
  df$das28_chg[idx_harm] <- df$das28_chg[idx_harm] + 2.0

  result <- forestsearch(
    df.analysis      = df,
    confounders.name = c("sero", "age"),
    outcome.name     = "das28_chg",
    event.name       = "das28_chg",
    treat.name       = "treat",
    id.name          = "id",
    outcome_type     = "continuous",
    effect_measure   = "MD",
    adverse_outcome  = TRUE,
    hr.threshold     = 0.0,
    hr.consistency   = 0.0,
    pconsistency.threshold = 0.30,
    fs.splits        = 20,
    n.min            = 20,
    maxk             = 1,
    use_lasso        = TRUE,
    use_grf          = TRUE,
    is.RCT           = TRUE,
    seedit           = 42,
    details          = FALSE,
    plot.sg          = FALSE,
    parallel_args    = list(plan = "sequential", workers = 1,
                            show_message = FALSE)
  )

  expect_true(is.list(result))
  expect_true(length(names(result)) > 1,
              info = "Result has only sg.harm -- pipeline failed silently")

  if (!is.null(result$error_log)) {
    fail(sprintf("Pipeline failed at %s: %s",
                 result$error_log$stage, result$error_log$message))
  }
})


test_that("count pipeline completes without internal failure", {
  set.seed(42)
  n <- 500
  df <- data.frame(
    id           = seq_len(n),
    treat        = rbinom(n, 1, 0.5),
    eos_low      = as.factor(rbinom(n, 1, 0.65)),
    age          = round(rnorm(n, 60, 12)),
    follow_up    = round(runif(n, 0.5, 3.0), 2)
  )
  # Poisson outcome with planted harm in eos_low=0
  log_rate <- with(df, {
    -1.0 + 0.02 * (age - 60) +
      1.0  * treat * (eos_low == 0) +
      -0.5 * treat * (eos_low == 1) +
      log(follow_up)
  })
  df$exacerbations <- rpois(n, lambda = exp(log_rate))

  result <- forestsearch(
    df.analysis      = df,
    confounders.name = c("eos_low", "age"),
    outcome.name     = "exacerbations",
    event.name       = "exacerbations",
    treat.name       = "treat",
    id.name          = "id",
    outcome_type     = "count",
    effect_measure   = "IRR",
    offset.name      = "follow_up",
    adverse_outcome  = TRUE,
    hr.threshold     = 1.10,
    hr.consistency   = 1.0,
    pconsistency.threshold = 0.30,
    fs.splits        = 20,
    n.min            = 20,
    maxk             = 1,
    use_lasso        = TRUE,
    use_grf          = TRUE,
    is.RCT           = TRUE,
    seedit           = 42,
    details          = FALSE,
    plot.sg          = FALSE,
    parallel_args    = list(plan = "sequential", workers = 1,
                            show_message = FALSE)
  )

  expect_true(is.list(result))
  expect_true(length(names(result)) > 1,
              info = "Result has only sg.harm -- pipeline failed silently")

  if (!is.null(result$error_log)) {
    fail(sprintf("Pipeline failed at %s: %s",
                 result$error_log$stage, result$error_log$message))
  }
})


test_that("grf_cuts format is character vector after GRF step", {
  # Verify the normalization works for GLM GRF output
  set.seed(42)
  n <- 200
  df <- data.frame(
    id         = seq_len(n),
    treat      = rbinom(n, 1, 0.5),
    bm1        = as.factor(rbinom(n, 1, 0.70)),
    age        = round(rnorm(n, 55, 10)),
    progressed = rbinom(n, 1, 0.3)
  )

  grf_res <- grf.subg.harm.glm(
    data             = df,
    confounders.name = c("bm1", "age"),
    outcome.name     = "progressed",
    treat.name       = "treat",
    id.name          = "id",
    outcome_type     = "binary",
    effect_measure   = "log_OR",
    n.min            = 20,
    seedit           = 42,
    verbose          = FALSE
  )

  tree_cuts <- grf_res$tree.cuts

  # tree.cuts from GLM GRF is a named list -- verify normalization works
  if (!is.null(tree_cuts) && is.list(tree_cuts)) {
    # Simulate the normalization from forestsearch_main.R
    normalized <- unlist(lapply(names(tree_cuts), function(nm) {
      paste0(nm, " <= ", tree_cuts[[nm]])
    }))
    expect_true(is.character(normalized))
    expect_true(all(grepl(" <= ", normalized)))
  }

  # Verify get_FSdata can consume the normalized cuts
  if (!is.null(tree_cuts)) {
    normalized <- if (is.list(tree_cuts) && !is.null(names(tree_cuts))) {
      unlist(lapply(names(tree_cuts), function(nm) {
        paste0(nm, " <= ", tree_cuts[[nm]])
      }))
    } else {
      tree_cuts
    }

    fsdata <- get_FSdata(
      df.analysis      = df,
      confounders.name = c("bm1", "age"),
      outcome.name     = "progressed",
      event.name       = "progressed",
      outcome_type     = "binary",
      use_lasso        = TRUE,
      use_grf          = TRUE,
      grf_cuts         = normalized,
      cont.cutoff      = 4
    )
    expect_true(is.list(fsdata))
    expect_true("confs_names" %in% names(fsdata))
  }
})
