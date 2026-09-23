# MR's consistency floor uses the screen's rounded admission rule
# (TASK_mr_admission_alignment_2026-09-23).  The search admits on
# round(Pcons, digits) >= p_star, i.e. on .fs_pcons_eff(p_star, digits); MR
# used to admit on the exact cutoff qnorm((1 + p_star) / 2), which is a
# different selection map from the one that produced H-hat.
#
# .fs_mr_select() is stubbed to record, per draw, the admitted set MR passes
# it together with the draw's consistency-standardized statistic, so the
# admission rule is checked on MR's own draws without restating its internals.

.mr_rounded_fixture <- function() {
  set.seed(20260923)
  n <- 400
  df <- data.frame(treat = rep(0:1, n / 2), x1 = rbinom(n, 1, .5),
                   x2 = rbinom(n, 1, .5), x3 = rbinom(n, 1, .5))
  df$y <- rbinom(n, 1, plogis(-0.4 + 0.5 * df$treat * df$x1))
  cands <- list(a = which(df$x1 == 1), b = which(df$x2 == 1),
                c = which(df$x3 == 1), d = which(df$x1 == 1 & df$x2 == 1),
                e = which(df$x1 == 1 | df$x3 == 1), f = which(df$x2 == 0),
                g = which(df$x1 == 1 & df$x3 == 0), h = which(df$x3 == 0))
  spec <- list(outcome_type = "binary", effect_measure = "OR",
               treat.name = "treat", outcome.name = "y", event.name = "y",
               offset.name = NULL, adjust_covariates = NULL,
               adverse_outcome = TRUE)
  list(df = df, cands = cands, spec = spec)
}

.mr_record_admission <- function(fx, admission, digits) {
  ns <- asNamespace("forestsearch")
  orig <- get(".fs_mr_select", envir = ns)
  rec <- new.env(); rec$draws <- list()
  unlockBinding(".fs_mr_select", ns)
  on.exit({
    assign(".fs_mr_select", orig, envir = ns)
    lockBinding(".fs_mr_select", ns)
  }, add = TRUE)
  assign(".fs_mr_select", function(bs, zc, sz, pass, ...) {
    rec$draws[[length(rec$draws) + 1L]] <- list(bs = bs, zc = zc, pass = pass)
    orig(bs, zc, sz, pass, ...)
  }, envir = ns)
  res <- forestsearch:::fs_mr_inference(
    df = fx$df, candidates = fx$cands, spec = fx$spec,
    selected_members = fx$cands$a, admission = admission,
    reselection = "maxcons", draws = 300L, ci_method = "ij", seed = 11L,
    include_complement = FALSE, pconsistency.digits = digits)
  list(res = res, draws = rec$draws)
}

test_that("MR admits on the rounded rule's effective threshold", {
  fx <- .mr_rounded_fixture()
  adm <- list(effect_floor = 0, consistency = list(c_cons = 0, p_star = 0.90))
  z_exact <- stats::qnorm((1 + 0.90) / 2)
  z_eff <- stats::qnorm((1 + forestsearch:::.fs_pcons_eff(0.90, 2L)) / 2)
  expect_lt(z_eff, z_exact)

  r2 <- .mr_record_admission(fx, adm, 2L)
  expect_gt(length(r2$draws), 0L)
  ok <- vapply(r2$draws, function(d)
    setequal(d$pass, which(d$bs >= 0 & d$zc >= z_eff - 1e-12)), logical(1))
  expect_true(all(ok))

  # The band [z_eff, z_exact) is populated on this fixture: some draw admits
  # a candidate the exact cutoff would have refused.  Without this the test
  # above could pass with the old rule.
  in_band <- vapply(r2$draws, function(d)
    any(d$zc[d$pass] < z_exact), logical(1))
  expect_true(any(in_band))
})

test_that("NULL digits falls back to the subgroup.consistency() default", {
  fx <- .mr_rounded_fixture()
  adm <- list(effect_floor = 0, consistency = list(c_cons = 0, p_star = 0.90))
  d0 <- eval(formals(forestsearch::subgroup.consistency)$pconsistency.digits)
  a <- .mr_record_admission(fx, adm, NULL)
  b <- .mr_record_admission(fx, adm, as.integer(d0))
  a$res$timing_seconds <- b$res$timing_seconds <- NULL
  expect_identical(a$res, b$res)
})

test_that("at 6 digits the floor is the exact cutoff to within rounding", {
  fx <- .mr_rounded_fixture()
  adm <- list(effect_floor = 0, consistency = list(c_cons = 0, p_star = 0.90))
  z_exact <- stats::qnorm((1 + 0.90) / 2)
  z6 <- stats::qnorm((1 + forestsearch:::.fs_pcons_eff(0.90, 6L)) / 2)
  r6 <- .mr_record_admission(fx, adm, 6L)
  ok <- vapply(r6$draws, function(d)
    setequal(d$pass, which(d$bs >= 0 & d$zc >= z6 - 1e-12)), logical(1))
  expect_true(all(ok))
  expect_lt(z_exact - z6, 1e-5)
})

test_that("no consistency floor: digits does not enter admission", {
  # GRF / DINA carry an effect floor only; p_star and digits must not reach
  # their admission set.
  fx <- .mr_rounded_fixture()
  adm <- list(effect_floor = 0, consistency = NULL)
  run <- function(dg) {
    r <- forestsearch:::fs_mr_inference(
      df = fx$df, candidates = fx$cands, spec = fx$spec,
      selected_members = fx$cands$a, admission = adm,
      reselection = "maxeff", draws = 300L, ci_method = "ij", seed = 11L,
      include_complement = FALSE, pconsistency.digits = dg)
    r$timing_seconds <- NULL
    r
  }
  expect_identical(run(2L), run(6L))
  expect_identical(run(2L), run(NULL))
})
