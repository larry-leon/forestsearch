# ============================================================================
# test-binary-default-or.R
#
# Directive B: forestsearch() resolves an unset effect_measure to "OR" for
# binary outcomes, not "RD".  These tests pin the resolved estimand and the
# thresholds and comparison scale that follow from it.
#
# Checks 1, 2, 7 and 8 fail against the pre-change source (the unset binary
# cell resolved to "RD" with 0.05 / 0.0 on the identity scale); checks 3-6
# and 9 pass both before and after and exist to catch collateral movement.
# ============================================================================

# -- shared fixtures ---------------------------------------------------------

.bd_data <- function(outcome_type) {
  switch(outcome_type,
    binary     = .make_binary_data(N = 150L, OR_harm = 3.0, seed = 7L),
    survival   = .make_survival_data(N = 150L, HR_harm = 2.0, seed = 7L),
    continuous = .make_continuous_data(N = 150L, MD_harm = 1.5, seed = 7L),
    count      = .make_count_data(N = 150L, IRR_harm = 2.0, seed = 7L))
}

.bd_conf <- function(outcome_type) {
  if (outcome_type == "survival") c("age", "stage", "sex")
  else c("age", "biomarker_hi", "sex")
}

# Thresholds are deliberately NOT passed unless `extra` supplies them: the
# defaults are what Directive B moves, and forestsearch() keys the remap on
# missing(hr.threshold).
.bd_fit <- function(outcome_type, effect_measure = NULL,
                    subgroup_method = "consistency", extra = list()) {
  args <- list(
    df.analysis            = .bd_data(outcome_type),
    confounders.name       = .bd_conf(outcome_type),
    outcome_type           = outcome_type,
    subgroup_method        = subgroup_method,
    treat.name             = "treat",
    id.name                = "id",
    outcome.name           = if (outcome_type == "survival") "time" else "y",
    sg_focus               = "maxSG",
    pconsistency.threshold = 0.80,
    max_subgroups_search   = 5L,
    use_grf                = TRUE,
    use_lasso              = FALSE,
    is.RCT                 = TRUE,
    maxk                   = 2L,
    n.min                  = 30L,
    d0.min                 = 5L,
    d1.min                 = 5L,
    fs.splits              = 20L,
    seedit                 = 42L,
    mr_inference           = FALSE,
    details                = FALSE,
    quiet                  = TRUE,
    plot.sg                = FALSE,
    plot.grf               = FALSE,
    parallel_args          = list(plan = "sequential", workers = 1L,
                                  show_message = FALSE))
  if (outcome_type == "survival") args$event.name  <- "event"
  if (outcome_type == "binary")   args$event.name  <- "y"
  if (outcome_type == "count")    args$offset.name <- "ftime"
  if (!is.null(effect_measure))   args$effect_measure <- effect_measure
  args <- utils::modifyList(args, extra)
  set.seed(20260918L)
  suppressWarnings(do.call(forestsearch, args))
}

.bd_tc <- function(fit)
  fit$threshold_config[c("effect_measure", "screening", "consistency",
                         "screening_natural", "consistency_natural", "scale")]


# -- check 1 -----------------------------------------------------------------

test_that("binary with effect_measure unset resolves to OR on ratio thresholds", {
  for (m in c("consistency", "dina", "grf")) {
    fit <- .bd_fit("binary", effect_measure = NULL, subgroup_method = m)
    expect_identical(fit$effect_measure, "OR", info = m)
    tc <- fit$threshold_config
    expect_identical(tc$effect_measure, "OR", info = m)
    expect_equal(tc$screening,   log(1.25), info = m)
    expect_equal(tc$consistency, log(1.00), info = m)
    expect_equal(tc$screening_natural,   1.25, info = m)
    expect_equal(tc$consistency_natural, 1.00, info = m)
    expect_identical(tc$scale, "log", info = m)
    # The identity-scale RD defaults must not survive anywhere in the config.
    expect_false(isTRUE(all.equal(tc$screening, 0.05)), info = m)
  }
})


# -- check 8 -----------------------------------------------------------------

test_that("the three identifiers screen a default binary call on one estimand", {
  tcs <- lapply(c("consistency", "dina", "grf"),
                function(m) .bd_tc(.bd_fit("binary", subgroup_method = m)))
  expect_identical(tcs[[2]], tcs[[1]])
  expect_identical(tcs[[3]], tcs[[1]])
})


# -- check 2 -----------------------------------------------------------------

test_that("binary unset and binary OR explicit resolve identically", {
  expect_identical(.bd_tc(.bd_fit("binary", effect_measure = NULL)),
                   .bd_tc(.bd_fit("binary", effect_measure = "OR")))
})


# -- check 3 -----------------------------------------------------------------

test_that("binary RD explicit keeps the identity scale and its own defaults", {
  tc <- .bd_tc(.bd_fit("binary", effect_measure = "RD"))
  expect_identical(tc$effect_measure, "RD")
  expect_equal(tc$screening,   0.05)
  expect_equal(tc$consistency, 0.00)
  expect_identical(tc$scale, "identity")
})


# -- check 4 -----------------------------------------------------------------

test_that("survival, MD and IRR resolution is untouched by the binary default", {
  sv <- .bd_fit("survival")
  expect_identical(sv$threshold_config$effect_measure, "HR")
  expect_equal(sv$threshold_config$screening_natural, 1.25)
  expect_identical(sv$threshold_config$scale, "log")

  for (em in list(NULL, "MD")) {
    tc <- .bd_tc(.bd_fit("continuous", effect_measure = em))
    expect_identical(tc$effect_measure, "MD")
    expect_equal(tc$screening,   0.0)
    expect_equal(tc$consistency, 0.0)
    expect_identical(tc$scale, "identity")
  }

  for (em in list(NULL, "IRR")) {
    tc <- .bd_tc(.bd_fit("count", effect_measure = em))
    expect_identical(tc$effect_measure, "IRR")
    expect_equal(tc$screening,   log(1.25))
    expect_equal(tc$consistency, log(1.00))
    expect_identical(tc$scale, "log")
  }
})


# -- check 5 -----------------------------------------------------------------

test_that("forestsearch() holds exactly one effect_measure resolution site", {
  src <- paste(deparse(body(forestsearch)), collapse = "\n")
  # The resolution is a switch on outcome_type whose arms are the three
  # non-survival types.  Count its occurrences in the deparsed body.
  hits <- gregexpr(
    "switch\\(outcome_type,\\s*binary\\s*=\\s*\"[A-Z]+\",\\s*continuous\\s*=",
    src)[[1]]
  n_sites <- if (identical(as.integer(hits[1]), -1L)) 0L else length(hits)
  expect_identical(n_sites, 1L)
  # And that one site yields OR for binary.
  expect_match(src, "binary\\s*=\\s*\"OR\"")
  expect_false(grepl("binary\\s*=\\s*\"RD\"", src))
})


# -- check 7 -----------------------------------------------------------------

test_that("a replayed args_call_all resolves the parent fit's estimand", {
  parent <- .bd_fit("binary", effect_measure = NULL)
  aca <- parent$args_call_all
  # The capture happens after resolution, so the replay carries the resolved
  # estimand rather than re-resolving it.
  expect_identical(aca$effect_measure, "OR")

  # Replay the way the bootstrap does: do.call(forestsearch, args_call_all).
  set.seed(20260918L)
  replay <- suppressWarnings(do.call(forestsearch, aca))
  expect_identical(replay$effect_measure, parent$effect_measure)
  expect_identical(.bd_tc(replay), .bd_tc(parent))
})


# -- check 6 -----------------------------------------------------------------

test_that("the campaign template's call re-resolves to OR on its own thresholds", {
  tpl <- testthat::test_path("..", "..", "quarto", "simulations", "actg175",
                             "binary_020", "sim_fs_mr_field_or_template.qmd")
  skip_if_not(file.exists(tpl), "campaign template not in this tree")
  txt <- readLines(tpl, warn = FALSE)

  # The template passes effect_measure explicitly, so Directive B cannot
  # reach it.  Assert that, rather than assuming it.
  expect_true(any(grepl('^\\s*effect_measure\\s*<-\\s*"OR"', txt)))
  expect_true(any(grepl('effect_measure\\s*=\\s*effect_measure', txt)))
  or_thr <- as.numeric(sub(".*<-\\s*([0-9.]+).*", "\\1",
                           grep("^or_threshold\\s*<-", txt, value = TRUE)[1]))
  or_con <- as.numeric(sub(".*<-\\s*([0-9.]+).*", "\\1",
                           grep("^or_consistency\\s*<-", txt, value = TRUE)[1]))
  expect_true(is.finite(or_thr) && is.finite(or_con))

  tc <- .bd_tc(.bd_fit("binary", effect_measure = "OR",
                       extra = list(effect.threshold = or_thr,
                                    consistency.threshold = or_con)))
  expect_identical(tc$effect_measure, "OR")
  expect_equal(tc$screening,   log(or_thr))
  expect_equal(tc$consistency, log(or_con))
  expect_identical(tc$scale, "log")
})


# -- check 9 -----------------------------------------------------------------

test_that("subgroup_search()'s hr.threshold roxygen no longer claims log for HR", {
  f <- testthat::test_path("..", "..", "R", "subgroup_search.R")
  skip_if_not(file.exists(f), "package source not in this tree")
  rox <- grep("^#'", readLines(f, warn = FALSE), value = TRUE)
  expect_false(any(grepl("log scale for ratio measures", rox)))
  expect_false(any(grepl("On the log scale for ratio measures \\(OR, HR\\)", rox)))
  expect_true(any(grepl("natural", rox)))
})
