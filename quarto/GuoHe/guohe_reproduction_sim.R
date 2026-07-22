# guohe_reproduction_sim.R
#
# Reproduction of the Guo & He (2021) simulation study, Section 5 -- the
# published operating characteristics (coverage, distance, bias) that the
# estimator-level validation in `guohe_algorithm3.R` does NOT cover.
#
# PURPOSE. `quarto/GuoHe/guohe_algorithm3.R` is validated as an ESTIMATOR: to
# machine precision against the authors' `synthetic code.R`, and to ~1e-3
# against published Table 2. Neither check touches coverage or bias. Tables 3-6
# of the paper report exactly those, in the authors' own regime (small, fixed
# candidate family). Reproducing them decouples two explanations of the
# large-family tuning sensitivity observed on GBSG at k = 1744:
#
#     (a) the method is tuning-sensitive at scale;   (b) our code is wrong.
#
# If the published tables reproduce in the authors' regime, (b) is excluded.
#
# SCOPE. This script does not modify, wrap or reimplement `guohe_algorithm3()`
# or `guohe_adaptive_r()`. It supplies a data-generating mechanism, an
# independent naive comparator, and scoring against the published targets.
#
# DESIGN, extracted verbatim from Guo & He (2021) Section 5.1 (see
# `guohe_simulation_dgm_extraction.md` for the extraction and its provenance):
#
#   lambda(t) = lambda_0(t) exp(beta_i D),  i = 1, ..., k
#   lambda_0  : Weibull(1, 1) baseline, i.e. the unit exponential
#   membership: equal probability 1/k over the k disjoint subgroups
#   D         : Bernoulli(0.5)
#   censoring : log(C) ~ Unif(-1.25, 1.00); implied rate ~41% at beta = 0
#   Tables 3-5: k = 2, n = 400, beta_1 = 0, beta_2 in {0, .1, .2, .3, .4, .5}
#   Table 6   : k in {2, 6, 10, 12}, n = 200k, all beta_i = 0
#   replicates: 2000;  r grid: 1/3, 1/12, 1/21, 1/30;  adaptive: v = 5
#
# ORIENTATION. Section 2.1 sets larger beta = better treatment effect, and the
# authors' code negates the Cox coefficient. The inferential parameter is
# therefore the NEGATED log hazard ratio, so subgroup i has true oriented effect
# -beta_i, and `orient = -1` in `guohe_algorithm3()`. With beta_1 = 0 and
# beta_2 >= 0, subgroup 1 is the best; at beta_2 = 0 the two coincide and
# |H| = 2, which is the maximal-bias configuration.
#
# Run:  Rscript quarto/GuoHe/guohe_reproduction_sim.R          # verify + time
#       Rscript quarto/GuoHe/guohe_reproduction_sim.R --full   # full study

suppressMessages(library(survival))

# guohe_algorithm3() / guohe_adaptive_r() are package API (exported by
# forestsearch); the sourced copies this script used to load were promoted
# to R/ and no longer live in this directory.
suppressMessages(library(forestsearch))

# Published grid, Tables 3-6. NOT the package defaults of `guohe_adaptive_r()`.
GH_R_GRID <- c(1 / 3, 1 / 12, 1 / 21, 1 / 30)

# ---- data-generating mechanism --------------------------------------------

#' Generate one replication under the Guo & He Section 5 design
#'
#' @param beta Numeric vector of length `k`; the DGM log hazard ratios. The
#'   oriented (inferential) effects are `-beta`.
#' @param n Sample size.
#' @return A data frame with `time`, `event`, `treat`, subgroup label `z`, and
#'   membership columns `S1, ..., Sk`.
gh_sim_data <- function(beta, n) {
  k <- length(beta)
  z <- sample.int(k, n, replace = TRUE)
  d <- stats::rbinom(n, 1L, 0.5)
  t_event <- stats::rexp(n, rate = exp(beta[z] * d))
  t_cens <- exp(stats::runif(n, -1.25, 1.00))
  out <- data.frame(
    time = pmin(t_event, t_cens),
    event = as.integer(t_event <= t_cens),
    treat = d,
    z = z
  )
  for (j in seq_len(k)) out[[paste0("S", j)]] <- as.integer(z == j)
  out
}

gh_candidates <- function(k) paste0("S", seq_len(k))

# ---- naive comparator ------------------------------------------------------
# "select the better subgroup from beta_hat_i and proceed as if the subgroup
# were selected independent of the data" (Section 5.1). Computed here rather
# than taken from `guohe_algorithm3()`, which does not return standard errors,
# so the naive column is an INDEPENDENT calculation.

#' Per-subgroup oriented effect and standard error on the full sample
gh_subgroup_fits <- function(df, beta, orient = -1) {
  k <- length(beta)
  est <- se <- rep(NA_real_, k)
  for (j in seq_len(k)) {
    sub <- df[df[[paste0("S", j)]] == 1L, , drop = FALSE]
    fit <- try(suppressWarnings(
      survival::coxph(survival::Surv(time, event) ~ treat, data = sub)
    ), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      est[j] <- orient * unname(fit$coefficients[[1L]])
      se[j] <- unname(sqrt(diag(stats::vcov(fit)))[1L])
    }
  }
  list(est = est, se = se)
}

#' Naive selection, lower bound and bias
gh_naive <- function(fits, beta, orient = -1, level = 0.05) {
  truth <- orient * beta
  ok <- is.finite(fits$est)
  sel <- which.max(replace(fits$est, !ok, -Inf))
  beta_s <- truth[sel]
  gamma_max <- fits$est[sel]
  lower <- gamma_max - stats::qnorm(1 - level) * fits$se[sel]
  list(
    sel = sel, beta_s = beta_s, point = gamma_max, lower = lower,
    cover = as.integer(lower <= beta_s), dist = beta_s - lower,
    bias = gamma_max - beta_s
  )
}

# ---- back-transform: report scale -> oriented score scale ------------------
# Exact inverse of `guohe_algorithm3()`'s `to_report()`.
gh_to_score <- function(report, orient = -1, outcome = "survival") {
  if (outcome == "continuous") orient * report else orient * log(report)
}

# ---- one replication -------------------------------------------------------

#' One replication: naive plus one Guo-He fit per value of `r`
#'
#' @param adaptive Logical; if `TRUE`, additionally run Algorithm 2. This is the
#'   dominant cost -- it fits `v * length(r_grid)` bootstraps per replication.
gh_one_rep <- function(beta, n, r_grid = GH_R_GRID, B = 2000L, v = 5L,
                       level = 0.05, adaptive = FALSE, min_events = 5L,
                       seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  k <- length(beta)
  cand <- gh_candidates(k)
  truth <- -beta

  df <- gh_sim_data(beta, n)
  fits <- gh_subgroup_fits(df, beta)
  nv <- gh_naive(fits, beta, level = level)

  row <- data.frame(
    k = k, n = n, beta2 = if (k == 2L) beta[2L] else NA_real_,
    cens_rate = mean(df$event == 0L),
    naive_sel = nv$sel, naive_beta_s = nv$beta_s,
    naive_cover = nv$cover, naive_dist = nv$dist, naive_bias = nv$bias
  )

  # COMMON RANDOM NUMBERS across the r grid. `guohe_algorithm3()` draws its
  # resample indices under `seed` BEFORE `r` enters (see its lines 352-353), so
  # passing one `boot_seed` to every r gives all columns identical bootstrap
  # draws. Verified by execution: reconstructing the bound for each r from a
  # single score matrix agrees with the engine to 0.000e+00. This removes
  # between-column Monte Carlo noise, which matters because Table 3's r columns
  # differ by as little as 0.000-0.002.
  boot_seed <- if (is.null(seed)) NULL else as.integer(seed) + 500000L

  for (i in seq_along(r_grid)) {
    fit <- try(suppressWarnings(guohe_algorithm3(
      data = df, outcome = "survival", treatment = "treat",
      candidates = cand, time = "time", event = "event",
      orient = -1, B = B, r = r_grid[i], level = level, seed = boot_seed,
      min_events = min_events, diagnostics = FALSE
    )), silent = TRUE)
    tag <- sprintf("r%d", i)
    if (inherits(fit, "try-error")) {
      row[[paste0(tag, "_cover")]] <- NA_integer_
      row[[paste0(tag, "_dist")]] <- NA_real_
      row[[paste0(tag, "_bias")]] <- NA_real_
      next
    }
    sel_i <- match(fit$selected, cand)
    beta_s <- truth[sel_i]
    lower <- gh_to_score(fit$bound_one_sided)
    debiased <- gh_to_score(fit$debiased)
    row[[paste0(tag, "_cover")]] <- as.integer(lower <= beta_s)
    row[[paste0(tag, "_dist")]] <- beta_s - lower
    row[[paste0(tag, "_bias")]] <- debiased - beta_s
  }

  if (isTRUE(adaptive)) {
    ar <- try(suppressWarnings(guohe_adaptive_r(
      data = df, outcome = "survival", treatment = "treat",
      candidates = cand, time = "time", event = "event",
      orient = -1, r_grid = r_grid, v = v, B = B, level = level,
      min_events = min_events, refit = TRUE
    )), silent = TRUE)
    if (inherits(ar, "try-error")) {
      row$ad_r <- NA_real_
      row$ad_cover <- NA_integer_
      row$ad_dist <- row$ad_bias <- NA_real_
    } else {
      sel_i <- match(ar$fit$selected, cand)
      beta_s <- truth[sel_i]
      lower <- gh_to_score(ar$fit$bound_one_sided)
      debiased <- gh_to_score(ar$fit$debiased)
      row$ad_r <- ar$r_hat
      row$ad_cover <- as.integer(lower <= beta_s)
      row$ad_dist <- beta_s - lower
      row$ad_bias <- debiased - beta_s
    }
  }
  row
}

# ---- scenario driver -------------------------------------------------------

#' Run one scenario to `n_rep` replications and summarise against the tables
gh_run_scenario <- function(beta, n, n_rep, r_grid = GH_R_GRID, B = 2000L,
                            v = 5L, adaptive = FALSE, seed = 1L,
                            progress_every = 0L) {
  rows <- vector("list", n_rep)
  for (m in seq_len(n_rep)) {
    rows[[m]] <- gh_one_rep(beta, n, r_grid = r_grid, B = B, v = v,
                            adaptive = adaptive, seed = seed + m)
    if (progress_every > 0L && m %% progress_every == 0L) {
      cat(sprintf("    rep %d / %d\n", m, n_rep))
      utils::flush.console()
    }
  }
  do.call(rbind, rows)
}

#' Collapse replications into the published table layout
gh_summarise <- function(res, r_grid = GH_R_GRID) {
  mk <- function(tag, label) {
    data.frame(
      column = label,
      cover = mean(res[[paste0(tag, "_cover")]], na.rm = TRUE),
      dist = mean(res[[paste0(tag, "_dist")]], na.rm = TRUE),
      bias = mean(res[[paste0(tag, "_bias")]], na.rm = TRUE)
    )
  }
  out <- do.call(rbind, lapply(seq_along(r_grid), function(i) {
    mk(sprintf("r%d", i), sprintf("r = %.4f", r_grid[i]))
  }))
  out <- rbind(out, mk("naive", "Naive"))
  if ("ad_cover" %in% names(res)) out <- rbind(out, mk("ad", "Adaptive"))
  out$mcse_cover <- sqrt(out$cover * (1 - out$cover) / nrow(res))
  out
}

# ---- published targets -----------------------------------------------------

GH_TABLE_3 <- data.frame(
  beta2 = c(0, 1, 2, 3, 4, 5) / 10,
  r_1_3 = c(0.933, 0.926, 0.928, 0.941, 0.939, 0.952),
  r_1_12 = c(0.950, 0.945, 0.949, 0.957, 0.955, 0.965),
  r_1_21 = c(0.952, 0.947, 0.951, 0.959, 0.956, 0.965),
  r_1_30 = c(0.952, 0.947, 0.951, 0.959, 0.957, 0.966),
  naive = c(0.896, 0.912, 0.910, 0.919, 0.927, 0.934),
  adaptive = c(0.943, 0.936, 0.939, 0.947, 0.945, 0.953)
)

GH_TABLE_5 <- data.frame(
  beta2 = c(0, 1, 2, 3, 4, 5) / 10,
  r_1_3 = c(0.028, 0.024, 0.005, -0.003, -0.018, -0.027),
  r_1_12 = c(0.008, 0.002, -0.021, -0.045, -0.063, -0.067),
  r_1_21 = c(0.007, 0.000, -0.022, -0.036, -0.065, -0.070),
  r_1_30 = c(0.006, -0.001, -0.023, -0.037, -0.066, -0.071),
  naive = c(0.107, 0.100, 0.077, 0.061, 0.029, 0.022),
  adaptive = c(0.018, 0.012, -0.008, -0.018, -0.042, -0.040)
)

GH_TABLE_6 <- data.frame(
  k = c(2, 6, 10, 12),
  cover_r_1_3 = c(0.929, 0.891, 0.866, 0.860),
  cover_r_1_12 = c(0.952, 0.941, 0.944, 0.946),
  cover_r_1_21 = c(0.953, 0.943, 0.949, 0.950),
  cover_r_1_30 = c(0.953, 0.945, 0.950, 0.950),
  cover_naive = c(0.900, 0.739, 0.594, 0.543),
  cover_adaptive = c(0.939, 0.930, 0.927, 0.925),
  bias_r_1_3 = c(0.029, 0.060, 0.066, 0.062),
  bias_r_1_12 = c(0.006, 0.011, 0.009, 0.003),
  bias_r_1_21 = c(0.004, 0.009, 0.006, 0.001),
  bias_r_1_30 = c(0.004, 0.008, 0.005, 0.001),
  bias_naive = c(0.105, 0.240, 0.290, 0.302),
  bias_adaptive = c(0.014, 0.029, 0.031, 0.026)
)

# ---- verification and timing ----------------------------------------------
# Runs when the script is executed rather than sourced. Fail-loud: every block
# prints one PASS/FAIL line and the script exits non-zero on any failure.

if (sys.nframe() == 0L) {
  PF <- function(ok) if (isTRUE(ok)) "PASS" else "FAIL"
  fails <- character(0)
  note <- function(id, ok, msg) {
    cat(sprintf("%s %s -- %s\n", id, PF(ok), msg))
    if (!isTRUE(ok)) fails <<- c(fails, id)
  }
  args <- commandArgs(trailingOnly = TRUE)

  cat("=== V1: DGM structure, k = 2, beta_2 = 0, large sample ===\n")
  set.seed(20260719)
  big <- gh_sim_data(c(0, 0), 400000L)
  cr <- mean(big$event == 0L)
  pr_s1 <- mean(big$S1)
  pr_tr <- mean(big$treat)
  # analytic censoring probability under beta = 0:
  #   P(cens) = E_U[exp(-e^U)], U ~ Unif(-1.25, 1)
  cr_theory <- stats::integrate(function(u) exp(-exp(u)), -1.25, 1.00)$value / 2.25
  note("V1a", abs(cr - cr_theory) < 0.005,
       sprintf("censoring %.4f vs analytic %.4f (paper: 'about 40%%')",
               cr, cr_theory))
  note("V1b", abs(pr_s1 - 0.5) < 0.005 && abs(pr_tr - 0.5) < 0.005,
       sprintf("P(S1) = %.4f, P(treat) = %.4f, both target 0.5", pr_s1, pr_tr))

  cat("\n=== V2: DGM recovers the specified log hazard ratios ===\n")
  set.seed(20260720)
  big2 <- gh_sim_data(c(0, 0.5), 400000L)
  b_hat <- vapply(1:2, function(j) {
    s <- big2[big2[[paste0("S", j)]] == 1L, ]
    unname(survival::coxph(survival::Surv(time, event) ~ treat, data = s)$coefficients[[1L]])
  }, numeric(1))
  note("V2", max(abs(b_hat - c(0, 0.5))) < 0.02,
       sprintf("fitted log-HR = (%.4f, %.4f) vs specified (0.000, 0.500)",
               b_hat[1], b_hat[2]))
  cr2 <- mean(big2$event == 0L)
  cat(sprintf("     censoring at beta_2 = 0.5: %.4f (analytic 0.336)\n", cr2))

  cat("\n=== V3: offsets and orientation, one replication ===\n")
  set.seed(11)
  df1 <- gh_sim_data(c(0, 0), 400L)
  fits1 <- gh_subgroup_fits(df1, c(0, 0))
  fit1 <- suppressWarnings(guohe_algorithm3(
    data = df1, outcome = "survival", treatment = "treat",
    candidates = gh_candidates(2L), time = "time", event = "event",
    orient = -1, B = 200L, r = 1 / 30, level = 0.05, seed = 99L,
    diagnostics = FALSE
  ))
  gmax_hand <- max(fits1$est)
  off_hand <- (1 - 400^(1 / 30 - 0.5)) * (gmax_hand - fits1$est)
  note("V3a", abs(gh_to_score(fit1$naive) - gmax_hand) < 1e-10,
       sprintf("gamma_max: engine %.10f vs independent coxph %.10f",
               gh_to_score(fit1$naive), gmax_hand))
  note("V3b", identical(match(fit1$selected, gh_candidates(2L)),
                        which.max(fits1$est)),
       sprintf("selected %s; independent argmax S%d",
               fit1$selected, which.max(fits1$est)))
  cat(sprintf("     hand offsets d = (%.6f, %.6f)\n", off_hand[1], off_hand[2]))

  cat("\n=== V4: report/score round trip is exact ===\n")
  rt <- gh_to_score(fit1$bound_one_sided)
  note("V4", abs(exp(-rt) - fit1$bound_one_sided) < 1e-12,
       sprintf("bound: report %.8f <-> score %.8f", fit1$bound_one_sided, rt))

  cat("\n=== V5: de-biasing moves the estimate toward the null ===\n")
  # At beta_1 = beta_2 = 0 the naive estimate is optimistic by construction
  # (|H| = 2). The de-biased score must be no larger, and the certified lower
  # bound must be below both. This is a direction check, not a coverage claim.
  s_naive <- gh_to_score(fit1$naive)
  s_deb <- gh_to_score(fit1$debiased)
  s_low <- gh_to_score(fit1$bound_one_sided)
  note("V5", s_deb <= s_naive && s_low < s_deb,
       sprintf("naive %.4f >= debiased %.4f > lower %.4f",
               s_naive, s_deb, s_low))

  cat("\n=== V6: timing ===\n")
  time_one <- function(lab, beta, n, B, adaptive) {
    t0 <- proc.time()[["elapsed"]]
    invisible(gh_one_rep(beta, n, B = B, adaptive = adaptive, seed = 7L))
    el <- proc.time()[["elapsed"]] - t0
    cat(sprintf("   %-46s %8.2f s\n", lab, el))
    el
  }
  t_fix_b200 <- time_one("k=2  n=400   B=200   fixed r only", c(0, 0), 400L, 200L, FALSE)
  t_fix_b2k <- time_one("k=2  n=400   B=2000  fixed r only", c(0, 0), 400L, 2000L, FALSE)
  t_ad_b200 <- time_one("k=2  n=400   B=200   fixed + adaptive", c(0, 0), 400L, 200L, TRUE)
  t_k12_b2k <- time_one("k=12 n=2400  B=2000  fixed r only", rep(0, 12), 2400L, 2000L, FALSE)

  cat("\n   --- extrapolation to the full study (single core) ---\n")
  hrs <- function(s) s / 3600
  cat(sprintf("   Tables 3-5 fixed r, 6 beta_2 x 2000 reps, B=2000 : %8.1f h\n",
              hrs(6 * 2000 * t_fix_b2k)))
  cat(sprintf("   Tables 3-5 adaptive, 6 beta_2 x 2000 reps, B=200 : %8.1f h\n",
              hrs(6 * 2000 * (t_ad_b200 - t_fix_b200))))
  cat(sprintf("   Table 6 fixed r, k=12 slice only, 2000 reps      : %8.1f h\n",
              hrs(2000 * t_k12_b2k)))

  cat("\n")
  if (length(fails)) {
    cat("FAILURES: ", paste(fails, collapse = ", "), "\n")
    quit(status = 1L)
  }
  cat("All verification blocks passed.\n")

  if ("--full" %in% args) {
    stop("--full is intended for the compute host; see guohe_reproduction_RUN.md")
  }
}
