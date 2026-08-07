# Contract tests for the consolidated betaHhat implementation (R/betaHhat_truth.R).
#
# These are the standing exception to the project's "no testthat scaffolding"
# convention: they assert STATIC properties of membership resolution and the
# record contract.  Everything integration-shaped stays in Quarto.
#
# Every test pairs its assertion with a plausible-but-wrong form asserted to
# FAIL.  An assertion that only ever passes cannot distinguish a correct
# implementation from one that happens to be close.
#
# Spec: dev/betaHhat-consolidation/SPEC_betaHhat_package_function.md

.resolve <- forestsearch:::.fs_resolve_membership
.rulecols <- forestsearch:::.fs_rule_columns

# A frame carrying the columns every family needs, so the SAME rule can be
# scored under each outcome_type and membership compared (T8).
.mk_frame <- function(n = 400L, seed = 11L) {
  set.seed(seed)
  x1 <- runif(n, 0, 100)
  x2 <- runif(n, 0, 100)
  treat_sim <- rbinom(n, 1L, 0.5)
  mu0 <- 10 + 0.10 * x1
  mu1 <- mu0 + ifelse(x1 > 50 & x2 <= 50, 8, 2)
  y_cont <- ifelse(treat_sim == 1L, mu1, mu0) + rnorm(n, sd = 2)
  y_bin  <- rbinom(n, 1L, plogis(-0.5 + 0.01 * x1 + 0.8 * treat_sim))
  t_time <- rexp(n, rate = 0.05 * exp(0.4 * treat_sim))
  event_sim <- rbinom(n, 1L, 0.85)
  data.frame(id = seq_len(n), x1 = x1, x2 = x2, treat_sim = treat_sim,
             mu0 = mu0, mu1 = mu1, y_cont = y_cont, y_bin = y_bin,
             t_time = t_time, event_sim = event_sim)
}

.RULES <- list(
  conjunction = c("{x1 > 50}", "{x2 <= 50}"),
  single_cut  = c("{x1 > 50}"),
  negation    = c("!{x1 <= 50}"),
  disjunction = "(x1 > 50 & x2 <= 50) | (x1 <= 20)",
  empty       = c("{x1 > 1000}")
)

.FAMILIES <- list(
  survival   = list(outcome_type = "survival",   effect_measure = NULL,
                    outcome.name = "t_time", event.name = "event_sim"),
  binary     = list(outcome_type = "binary",     effect_measure = NULL,
                    outcome.name = "y_bin",  event.name = "event_sim"),
  continuous = list(outcome_type = "continuous", effect_measure = "MD",
                    outcome.name = "y_cont", event.name = "event_sim")
)

.one <- function(rule, frame, fam, focus = "harm", ...) {
  f <- .FAMILIES[[fam]]
  fs_betaHhat_one(rule, frame, focus = focus,
                  outcome_type = f$outcome_type,
                  effect_measure = f$effect_measure,
                  outcome.name = f$outcome.name,
                  event.name = f$event.name,
                  treat.name = "treat_sim", ...)
}


# --- T1: the partition invariant -------------------------------------------
test_that("T1: nH_eval + nHc_eval == N for every rule shape and family", {
  fr <- .mk_frame()
  N  <- nrow(fr)
  for (fam in names(.FAMILIES)) {
    for (rn in names(.RULES)) {
      for (fc in c("harm", "benefit")) {
        r <- .one(.RULES[[rn]], fr, fam, focus = fc)
        expect_identical(r$nH_eval + r$nHc_eval, N,
                         info = paste(fam, rn, fc))
      }
    }
  }
})

test_that("T1 negative control: partial-NA membership FAILS the sum", {
  fr <- .mk_frame()
  N  <- nrow(fr)
  # What get_dfpred() returns for a rule naming a column the frame lacks.
  pred <- suppressWarnings(get_dfpred(fr, c("{x1 > 50}", "{nosuch > 1}")))
  inH  <- pred$treat.recommend == 0L
  # This is the shape the resolver refuses to return.
  expect_true(anyNA(inH))
  expect_lt(sum(inH, na.rm = TRUE) + sum(!inH, na.rm = TRUE), N)
  # And the resolver does refuse it.
  expect_identical(.resolve(fr, c("{x1 > 50}", "{nosuch > 1}"))$status,
                   "unresolved")
})


# --- T2: unresolvable rules -------------------------------------------------
test_that("T2: a rule naming a missing column gives an all-NA unresolved record", {
  fr <- .mk_frame()
  for (fam in names(.FAMILIES)) {
    r <- .one(c("{nosuch > 1}"), fr, fam)
    expect_identical(r$status, "unresolved")
    expect_true(is.na(r$betaHhat_H))
    expect_true(is.na(r$betaHhat_Hc))
    expect_true(is.na(r$nH_eval))
    expect_true(is.na(r$nHc_eval))
    expect_identical(r$missing_cols, "nosuch")
  }
})

test_that("T2 negative control: the unguarded path yields NA membership", {
  fr <- .mk_frame()
  # The guard is doing work: without it, membership carries NA.
  pred <- suppressWarnings(get_dfpred(fr, c("{nosuch > 1}")))
  expect_true(anyNA(pred$treat.recommend))
  expect_identical(sum(is.na(pred$treat.recommend)), nrow(fr))
  # A partial rule is worse: it silently keeps the resolvable half.
  pred2 <- suppressWarnings(get_dfpred(fr, c("{x1 > 50}", "{nosuch > 1}")))
  expect_true(anyNA(pred2$treat.recommend))
  expect_gt(sum(!is.na(pred2$treat.recommend)), 0L)
})


# --- T3: GRF disjunctions ---------------------------------------------------
test_that("T3: a realized disjunction and its structured sg_def agree", {
  fx <- readRDS(test_path("fixtures", "grf_disjunction.rds"))
  fr <- fx$frame
  expect_true(grepl("|", fx$definition, fixed = TRUE))
  expect_true(fx$sg_def$is_disjunction)

  by_string <- .resolve(fr, fx$definition)
  by_struct <- .resolve(fr, rule = NULL, sg_def_struct = fx$sg_def)

  expect_identical(by_string$status, "ok")
  expect_identical(by_struct$status, "ok")
  expect_identical(by_string$in_region, by_struct$in_region)
  expect_gt(sum(by_struct$in_region), 0L)
})

test_that("T3 negative control: split-first FAILS on the same string", {
  fx <- readRDS(test_path("fixtures", "grf_disjunction.rds"))
  fr <- fx$frame
  ok <- .resolve(fr, fx$definition)$in_region

  # The pre-fix dispatch: split on " & " before testing for "|".
  shredded <- strsplit(fx$definition, " & ", fixed = TRUE)[[1L]]
  expect_gt(length(shredded), 1L)
  pred <- suppressWarnings(
    tryCatch(get_dfpred(fr, shredded), error = function(e) NULL))
  bad <- if (is.null(pred)) NULL else pred$treat.recommend
  # It does not reproduce the correct membership: either it errors, or it
  # returns NA membership.  It must not silently agree.
  expect_false(identical(bad, as.integer(!ok)))
  expect_true(is.null(bad) || anyNA(bad))
})


# --- T4: negation -----------------------------------------------------------
test_that("T4: negation is equivalent to the flipped comparison, and partitions", {
  fr <- .mk_frame()
  a <- .resolve(fr, c("!{x1 <= 50}"))$in_region
  b <- .resolve(fr, c("{x1 > 50}"))$in_region
  expect_identical(a, b)

  # A rule and its negation partition the frame.
  neg <- .resolve(fr, c("{x1 <= 50}"))$in_region
  expect_identical(a, !neg)
  expect_identical(sum(a) + sum(neg), nrow(fr))
})


# --- T5: the no-subgroup case ----------------------------------------------
test_that("T5: no subgroup gives nHc_eval == N and the ITT effect", {
  fr <- .mk_frame()
  N  <- nrow(fr)
  for (fam in names(.FAMILIES)) {
    r <- .one(NULL, fr, fam)
    expect_identical(r$nH_eval, 0L)
    expect_identical(r$nHc_eval, N)
    expect_true(is.na(r$betaHhat_H))          # empty region has no target
    f <- .FAMILIES[[fam]]
    itt <- forestsearch:::.fs_region_effect(
      fr, rep(TRUE, N), f$outcome_type, f$effect_measure,
      f$outcome.name, f$event.name, "treat_sim")
    expect_identical(r$betaHhat_Hc, itt)      # complement IS the ITT population
  }
})

test_that("T5 negative control: a real subgroup must NOT give nHc_eval == N", {
  fr <- .mk_frame()
  r <- .one(c("{x1 > 50}"), fr, "continuous")
  expect_gt(r$nH_eval, 0L)
  expect_lt(r$nHc_eval, nrow(fr))
})


# --- T8: membership is family-independent -----------------------------------
test_that("T8: the same rule gives identical membership across families", {
  fr <- .mk_frame()
  for (rn in names(.RULES)) {
    recs <- lapply(names(.FAMILIES), function(fam) .one(.RULES[[rn]], fr, fam))
    nH  <- vapply(recs, function(x) x$nH_eval,  integer(1))
    nHc <- vapply(recs, function(x) x$nHc_eval, integer(1))
    expect_identical(length(unique(nH)),  1L, info = rn)
    expect_identical(length(unique(nHc)), 1L, info = rn)
  }
})


# --- T9: exactness ----------------------------------------------------------
test_that("T9: repeated evaluation of a resolved target is bit-identical", {
  fr <- .mk_frame()
  for (fam in names(.FAMILIES)) {
    v <- vapply(1:15, function(i) .one(c("{x1 > 50}"), fr, fam)$betaHhat_H,
                numeric(1))
    expect_identical(max(v) - min(v), 0)      # not a tolerance
    expect_identical(length(unique(v)), 1L)
  }
})


# --- focus ------------------------------------------------------------------
test_that("focus is required, validated, and inverts the region", {
  fr <- .mk_frame()
  expect_error(fs_betaHhat_one(c("{x1 > 50}"), fr,
                               outcome_type = "continuous",
                               effect_measure = "MD"),
               "focus.*required")
  expect_error(.one(c("{x1 > 50}"), fr, "continuous", focus = "either"),
               "exactly one of")
  h <- .one(c("{x1 > 50}"), fr, "continuous", focus = "harm")
  b <- .one(c("{x1 > 50}"), fr, "continuous", focus = "benefit")
  expect_identical(h$nH_eval + b$nH_eval, nrow(fr))
  expect_identical(h$betaHhat_H, b$betaHhat_Hc)
})


# --- rule column extraction -------------------------------------------------
test_that(".fs_rule_columns finds the referenced columns in every rule shape", {
  expect_setequal(.rulecols(c("{x1 > 50}", "{x2 <= 50}")), c("x1", "x2"))
  expect_setequal(.rulecols(c("!{x1 <= 50}")), "x1")
  expect_setequal(.rulecols("(a > 1 & b <= 2) | (c > 3)"), c("a", "b", "c"))
  expect_setequal(.rulecols("{x1 > 50} & {x2 <= 50}"), c("x1", "x2"))
})


# --- T6: resolution accounting ----------------------------------------------
test_that("T6: counters match the unresolvable rules in a synthetic bundle", {
  fr <- .mk_frame()
  # 10 replicates: 4 resolvable (2 distinct), 4 unresolvable (2 distinct),
  # 2 undetected.
  sg <- c("{x1 > 50}", "{x1 > 50}", "{x2 <= 50}", "{x2 <= 50}",
          "{nosuch > 1}", "{nosuch > 1}", "{alsomissing <= 0}",
          "{alsomissing <= 0}", NA_character_, NA_character_)
  tb <- fs_betaHhat_table(sg, fr, focus = "harm",
                          outcome_type = "continuous", effect_measure = "MD")
  cnt <- fs_betaHhat_counts(tb)

  expect_identical(cnt$n_rules_total,      4L)
  expect_identical(cnt$n_rules_resolved,   2L)
  expect_identical(cnt$n_rules_unresolved, 2L)
  expect_identical(cnt$n_reps_total,      10L)
  expect_identical(cnt$n_reps_resolved,    4L)
  expect_identical(cnt$n_reps_unresolved,  4L)
  expect_identical(cnt$n_reps_undetected,  2L)
  # the accounting closes
  expect_identical(cnt$n_reps_resolved + cnt$n_reps_unresolved +
                     cnt$n_reps_undetected, cnt$n_reps_total)

  # and it survives the attach
  res <- data.frame(sim_id = seq_along(sg), sg_def = sg,
                    stringsAsFactors = FALSE)
  at <- fs_attach_betaHhat(res, fr, focus = "harm",
                           outcome_type = "continuous", effect_measure = "MD")
  expect_identical(fs_betaHhat_counts(at), cnt)
  expect_identical(sum(at$betaHhat_status == "unresolved", na.rm = TRUE), 4L)
  expect_identical(sum(is.na(at$betaHhat_H)), 6L)   # 4 unresolved + 2 undetected
})

test_that("T6 negative control: a clean bundle reports zero unresolved", {
  fr <- .mk_frame()
  sg <- rep(c("{x1 > 50}", "{x2 <= 50}"), 5)
  cnt <- fs_betaHhat_counts(
    fs_betaHhat_table(sg, fr, focus = "harm",
                      outcome_type = "continuous", effect_measure = "MD"))
  expect_identical(cnt$n_rules_unresolved, 0L)
  expect_identical(cnt$n_reps_unresolved,  0L)
  expect_identical(cnt$n_reps_resolved,   10L)
})


# --- T7: the parity guard ---------------------------------------------------
.mk_cov <- function() {
  base <- expand.grid(block = c("H", "Hc"),
                      estimator = c("oracle", "naive", "MR"),
                      target = c("C_dagger", "C_betaHhat"),
                      stringsAsFactors = FALSE)
  base$src <- "cell.rds"
  base$subgroup_method <- "consistency"
  base$n_sample <- 1000L
  base$n_eff <- 95L
  base$coverage <- 0.95
  base
}

test_that("T7: the parity guard passes on a clean table", {
  expect_silent(fs_betaHhat_neff_parity(.mk_cov()))
  expect_identical(fs_betaHhat_neff_parity(.mk_cov()), .mk_cov())
})

test_that("T7: the guard FIRES on an injected dropped target", {
  cv <- .mk_cov()
  i <- which(cv$target == "C_betaHhat")[1]
  cv$n_eff[i] <- cv$n_eff[i] - 1L          # one replicate silently dropped
  expect_error(fs_betaHhat_neff_parity(cv), "n_eff parity FAILED")
  expect_error(fs_betaHhat_neff_parity(cv), "1 of 6 cells")
})

test_that("T7: the guard FIRES when a C_betaHhat row is missing entirely", {
  cv <- .mk_cov()
  cv <- cv[-which(cv$target == "C_betaHhat")[1], , drop = FALSE]
  expect_error(fs_betaHhat_neff_parity(cv), "n_eff parity FAILED")
})

test_that("T7: strict = FALSE warns without stopping", {
  cv <- .mk_cov()
  cv$n_eff[which(cv$target == "C_betaHhat")[1]] <- 0L
  expect_warning(out <- fs_betaHhat_neff_parity(cv, strict = FALSE),
                 "n_eff parity FAILED")
  expect_s3_class(out, "data.frame")        # it returned rather than stopped
})

test_that("T7: a table without the required columns is a no-op, not an error", {
  expect_silent(fs_betaHhat_neff_parity(data.frame(a = 1)))
})

test_that("T7 companion: the report separates unresolvable from degenerate", {
  fr <- .mk_frame()
  sg <- c("{x1 > 50}", "{nosuch > 1}", "{x1 > 1000}")
  res <- data.frame(sim_id = 1:3, sg_def = sg, stringsAsFactors = FALSE)
  at <- fs_attach_betaHhat(res, fr, focus = "harm",
                           outcome_type = "continuous", effect_measure = "MD")
  rep <- fs_betaHhat_neff_report(at, "H")
  expect_setequal(rep$sg_def, c("{nosuch > 1}", "{x1 > 1000}"))
  # the two causes are distinguishable by status
  expect_identical(rep$status[rep$sg_def == "{nosuch > 1}"], "unresolved")
  expect_identical(rep$status[rep$sg_def == "{x1 > 1000}"],  "ok")
  # and a clean bundle reports nothing
  at2 <- fs_attach_betaHhat(data.frame(sg_def = "{x1 > 50}"), fr,
                            focus = "harm", outcome_type = "continuous",
                            effect_measure = "MD")
  expect_message(r2 <- fs_betaHhat_neff_report(at2, "H"), "no NA")
  expect_null(r2)
})
