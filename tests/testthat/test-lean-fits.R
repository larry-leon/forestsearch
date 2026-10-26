# =============================================================================
# Lean search-stage fits (0.4.0.9000): the unadjusted survival fit calls
# survival::coxph.fit() directly (fit_cox_for_subgroup()), and the unadjusted
# binary OR fit runs through a tagged glm.fit() closure and vector path
# (.make_glm_or_estimator_fast(), fit_glm_for_subgroup_fast()).  Both must be
# bit-identical to the formula-interface fits they replace, warnings included.
# =============================================================================

# Collects the warning texts raised by expr (muffled) alongside its value.
.lf_with_warnings <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(expr, warning = function(cond) {
    w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w)
}

# CRAN 0.4.0 R/subgroup_search.R L834-L845, L873, L875-L881, L888-L894, verbatim
ref_fit_cox <- function(Y, E, Treat) {
  hr.cox <- try({
    fit.cox <- survival::coxph(survival::Surv(Y, E) ~ Treat, robust = FALSE)
    beta <- fit.cox$coefficients
    se_beta <- sqrt(diag(fit.cox$var))
    z_ci <- stats::qnorm((1 + 0.95) / 2, 0, 1)
    ci <- cbind(exp(beta), exp(-beta),
                exp(beta - z_ci * se_beta), exp(beta + z_ci * se_beta))
    dimnames(ci) <- list(names(beta),
                         c("exp(coef)", "exp(-coef)",
                           "lower .95", "upper .95"))
    ci
  }, silent = TRUE)
  if (inherits(hr.cox, "try-error")) return(NULL)
  if (is.null(dim(hr.cox))) {
    hr.cox <- matrix(hr.cox, nrow = 1L, dimnames = list("Treat", names(hr.cox)))
  }
  trow <- if ("Treat" %in% rownames(hr.cox)) hr.cox["Treat", ] else hr.cox[1, ]
  list(hr = trow[1], lower = trow[3], upper = trow[4], med0 = NA_real_, med1 = NA_real_)
}

.lf_expect_cox_same <- function(Y, E, Treat, label) {
  a <- .lf_with_warnings(
    forestsearch:::fit_cox_for_subgroup(Y, E, Treat, rep(1L, length(Y))))
  b <- .lf_with_warnings(ref_fit_cox(Y, E, Treat))
  expect_identical(a$value, b$value, label = paste(label, "value"))
  expect_identical(a$warnings, b$warnings, label = paste(label, "warnings"))
  invisible(b)
}


test_that("A: fit_cox_for_subgroup() is bit-identical to the coxph() formula fit", {
  g <- survival::gbsg
  Yg <- g$rfstime / 30.4375
  Eg <- g$status
  Tg <- g$hormon

  set.seed(1)
  for (i in 1:25) {
    idx <- sort(sample(nrow(g), sample(40:300, 1L)))
    .lf_expect_cox_same(Yg[idx], Eg[idx], Tg[idx], sprintf("random subset %d", i))
  }

  # All rows.
  .lf_expect_cox_same(Yg, Eg, Tg, "all rows")

  # Near ties: ten clusters of three times that aeqSurv() collapses to exact
  # ties (the tripwire for the timefix convention).
  Y_tie <- rep(seq(5, 50, by = 5), each = 3) + rep(c(0, 1e-9, 2e-9), 10)
  .lf_expect_cox_same(Y_tie, Eg[1:30], Tg[1:30], "near-tie set")

  # Monotone likelihood: every treated event precedes every control event.
  Y_mono <- c(21:40, 1:20)
  T_mono <- rep(c(0L, 1L), each = 20L)
  r_mono <- .lf_expect_cox_same(Y_mono, rep(1L, 40L), T_mono, "monotone likelihood")
  expect_true(any(grepl("Loglik converged before variable", r_mono$warnings)))

  # n = 8 with at least one event per arm.
  .lf_expect_cox_same(c(5, 8, 12, 3, 9, 15, 7, 11), c(1, 0, 1, 1, 1, 0, 1, 0),
                      c(0, 0, 0, 0, 1, 1, 1, 1), "n = 8")

  # Single arm: an all-NA row and no warning on both paths.
  r_one <- .lf_expect_cox_same(Yg[1:50], Eg[1:50], rep(0L, 50L), "single arm")
  expect_length(r_one$warnings, 0L)
  expect_true(all(is.na(unlist(r_one$value[c("hr", "lower", "upper")]))))
})


test_that("F: the fast OR closure and vector path are bit-identical to the slow chain", {
  df <- .make_binary_data(N = 300L, OR_harm = 2.5, seed = 7L)
  N <- nrow(df)

  # Slices as 0/1 membership vectors over df.
  memb <- function(rows) { x <- integer(N); x[rows] <- 1L; x }
  set.seed(11)
  slices <- lapply(1:25, function(i) memb(sample(N, sample(20:300, 1L))))
  names(slices) <- sprintf("random subset %d", 1:25)
  cell <- list(e0  = df$treat == 0 & df$y == 1, ne0 = df$treat == 0 & df$y == 0,
               e1  = df$treat == 1 & df$y == 1, ne1 = df$treat == 1 & df$y == 0)
  for (k in names(cell)) {
    slices[[paste("zero", k)]] <- memb(head(which(!cell[[k]]), 40L))
  }
  slices[["y == treat"]] <- memb(head(which(df$y == df$treat), 40L))
  # One event among many non-events, and one event per arm among many.
  slices[["one event"]] <- memb(c(head(which(df$y == 0), 120L), which(df$y == 1)[1L]))
  slices[["one event per arm"]] <- memb(c(head(which(df$y == 0), 120L),
                                          which(cell$e0)[1L], which(cell$e1)[1L]))
  slices[["five rows"]] <- memb(1:5)
  slices[["single arm"]] <- memb(head(which(df$treat == 1), 60L))

  for (adverse in c(TRUE, FALSE)) {
    slow <- forestsearch:::.make_glm_binary_estimator(
      "treat", "y", effect_measure = "OR", adverse_outcome = adverse)
    fast <- forestsearch:::.make_glm_or_estimator_fast(
      "treat", "y", adverse_outcome = adverse)
    for (nm in names(slices)) {
      id.x <- slices[[nm]]
      lab <- sprintf("[adverse = %s] %s", adverse, nm)
      sl <- df[id.x == 1, , drop = FALSE]
      a <- .lf_with_warnings(fast(sl))
      b <- .lf_with_warnings(slow(sl))
      expect_identical(a$value, b$value, label = paste(lab, "closure value"))
      expect_identical(a$warnings, b$warnings, label = paste(lab, "closure warnings"))
      a <- .lf_with_warnings(
        forestsearch:::fit_glm_for_subgroup_fast(df$y, df$treat, id.x, fast))
      b <- .lf_with_warnings(forestsearch:::fit_glm_for_subgroup(df, id.x, slow))
      expect_identical(a$value, b$value, label = paste(lab, "entry value"))
      expect_identical(a$warnings, b$warnings, label = paste(lab, "entry warnings"))
    }
    expect_null(forestsearch:::fit_glm_for_subgroup_fast(
      df$y, df$treat, slices[["five rows"]], fast))
    expect_true(isTRUE(forestsearch:::fit_glm_for_subgroup_fast(
      df$y, df$treat, slices[["single arm"]], fast)$nonestimable))
  }
})


test_that("routing: unadjusted binary OR takes the fast path; RD and adjusted OR do not", {
  skip_on_cran()
  count_fits <- function(expr) {
    cnt <- new.env(); cnt$slow <- 0L; cnt$fast <- 0L
    ns <- asNamespace("forestsearch")
    suppressMessages({
      trace("fit_glm_for_subgroup", where = ns, print = FALSE,
            tracer = bquote(assign("slow", get("slow", envir = .(cnt)) + 1L, envir = .(cnt))))
      trace("fit_glm_for_subgroup_fast", where = ns, print = FALSE,
            tracer = bquote(assign("fast", get("fast", envir = .(cnt)) + 1L, envir = .(cnt))))
    })
    on.exit(suppressMessages({
      untrace("fit_glm_for_subgroup", where = ns)
      untrace("fit_glm_for_subgroup_fast", where = ns)
    }), add = TRUE)
    force(expr)
    c(slow = cnt$slow, fast = cnt$fast)
  }
  df <- .make_binary_data(N = 300L, OR_harm = 2.5, seed = 7L)
  run <- function(...) {
    suppressWarnings(suppressMessages(forestsearch(
      df.analysis = df, confounders.name = c("age", "biomarker_hi", "sex"),
      outcome.name = "y", event.name = "y", treat.name = "treat", id.name = "id",
      outcome_type = "binary", adverse_outcome = TRUE,
      sg_focus = "maxSG", pconsistency.threshold = 0.80, is.RCT = TRUE,
      fs.splits = 20, n.min = 30, d0.min = 5L, d1.min = 5L, maxk = 2,
      use_lasso = FALSE, use_grf = FALSE, mr_inference = FALSE,
      seedit = 42L, details = FALSE, quiet = TRUE, plot.sg = FALSE,
      plot.grf = FALSE, parallel_args = list(plan = "sequential"), ...)))
  }

  or_un <- count_fits(run(effect_measure = "OR"))
  expect_identical(or_un[["slow"]], 0L)
  expect_gt(or_un[["fast"]], 0L)

  rd_un <- count_fits(run(effect_measure = "RD"))
  expect_gt(rd_un[["slow"]], 0L)
  expect_identical(rd_un[["fast"]], 0L)

  or_adj <- count_fits(run(effect_measure = "OR", adjust_covariates = "age"))
  expect_identical(or_adj[["fast"]], 0L)
})
