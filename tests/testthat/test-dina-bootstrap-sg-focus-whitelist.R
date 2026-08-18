# =============================================================================
# dina_subgroup_bootstrap(): sg_focus whitelist parity with the other entries
# =============================================================================
# This whitelist had drifted to a five-element subset of the canonical seven,
# so a focus that forestsearch() and dina_subgroup() both accept was rejected
# here and the bootstrap stopped before its first iteration.

.dbs_demo_df <- function(n = 300L, seed = 1L) {
  set.seed(seed)
  df <- data.frame(
    w  = stats::rbinom(n, 1, 0.5),
    x1 = stats::runif(n, -1, 1),
    x2 = stats::runif(n, -1, 1)
  )
  tau_x  <- 0.4 + 1.2 * df$x1
  df$y   <- 0.5 * df$x1 + df$w * tau_x + stats::rnorm(n)
  df
}

test_that("the whitelist is the shared canonical constant, not a local literal", {
  # Assigned from .FS_SG_FOCUS_CANONICAL: a re-inlined literal would pass
  # every functional test below while reintroducing exactly this drift.
  src <- body(dina_subgroup_bootstrap)
  found <- FALSE
  walk <- function(e) {
    if (!is.call(e)) return(invisible(NULL))
    if (length(e) == 3L && is.symbol(e[[1L]]) &&
        as.character(e[[1L]]) %in% c("<-", "=") &&
        is.symbol(e[[2L]]) &&
        identical(as.character(e[[2L]]), "valid_sg_focus") &&
        is.symbol(e[[3L]]) &&
        identical(as.character(e[[3L]]), ".FS_SG_FOCUS_CANONICAL")) {
      found <<- TRUE
    }
    for (i in seq_along(e)) tryCatch(walk(e[[i]]), error = function(err) NULL)
    invisible(NULL)
  }
  walk(src)
  expect_true(found)
})

test_that("a 2-iteration bootstrap completes under maxeffCons", {
  skip_on_cran()
  df <- .dbs_demo_df()
  bs <- suppressWarnings(suppressMessages(dina_subgroup_bootstrap(
    df = df, outcome = "y", treatment = "w",
    covariates = c("x1", "x2"), family = "gaussian",
    m_diff = 0.5, n_boot = 2L, seed = 1L, n_min = 60L)))
  expect_true(is.list(bs))

  bs_mec <- suppressWarnings(suppressMessages(dina_subgroup_bootstrap(
    df = df, outcome = "y", treatment = "w",
    covariates = c("x1", "x2"), family = "gaussian",
    m_diff = 0.5, n_boot = 2L, seed = 1L, n_min = 60L,
    sg_focus = "maxeffCons")))
  expect_true(is.list(bs_mec))
})

test_that("maxeff and the eff/maxcons aliases are accepted too", {
  skip_on_cran()
  df <- .dbs_demo_df()
  for (f in c("maxeff", "eff", "maxcons", "hr")) {
    expect_no_error(
      suppressWarnings(suppressMessages(dina_subgroup_bootstrap(
        df = df, outcome = "y", treatment = "w",
        covariates = c("x1", "x2"), family = "gaussian",
        m_diff = 0.5, n_boot = 2L, seed = 1L, n_min = 60L,
        sg_focus = f))))
  }
})

test_that("the rejection message enumerates the seven spellings and maxcons", {
  df <- .dbs_demo_df(n = 100L)
  # suppressWarnings() as in the sibling tests above: n_boot = 2L trips the
  # "< 100, percentile CI unstable" advisory, which fires before the sg_focus
  # whitelist rejects and is incidental to the message being asserted here.
  err <- tryCatch(
    suppressWarnings(dina_subgroup_bootstrap(
      df = df, outcome = "y", treatment = "w",
      covariates = c("x1", "x2"), family = "gaussian",
      m_diff = 0.5, n_boot = 2L, seed = 1L, sg_focus = "not_a_focus")),
    error = conditionMessage)
  for (s in c("maxSG", "minSG", "eff", "effMaxSG", "effMinSG",
              "maxeff", "maxeffCons", "maxcons")) {
    expect_match(err, s, fixed = TRUE)
  }
})
