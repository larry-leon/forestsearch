# guohe_sec52_sim.R
#
# Reproduction of Guo & He (2021) Section 5.2, Table 7 -- empirical coverage of
# the 95% lower bound of gamma_s (the best SELECTED subgroup effect) in the
# post-hoc identified, NESTED-family setting: S(c) = {W <= c}, c in [30, 60].
#
# PURPOSE. The Section 5.1 reproduction (Tables 3-6, disjoint families) is
# complete and committed. Section 5.2 is the nested regime closest to the
# `forestsearch` maxeff family on GBSG and is the remaining validation target
# for `guohe_algorithm3()` as an operating-characteristics engine.
#
# SCOPE. This script does not modify, wrap or reimplement `guohe_algorithm3()`.
# `guohe_adaptive_r()` is NOT used -- Table 7 has no Adaptive column. The DGM
# and the analytic censoring mixture live in `guohe_sec52_truth.R` (sourced
# below) and are NOT re-implemented here; this file supplies the candidate
# family, an independent naive comparator, the per-replicate harness, and the
# published Table 7 targets.
#
# DESIGN, stated in Guo & He (2021) Section 5.2 and verified against the PDF
# text layer on 2026-07-21 (see guohe_sec52_truth.R header for the full DGM):
#
#   n = 400;  hazard lambda_0(t) e^{b(W) D}, lambda_0 unit exponential;
#   D ~ Bernoulli(0.5) indep. of W ~ Unif[0, 80];
#   log(C) ~ Unif(-1.25, 1.00)  (inherited verbatim from Sec 5.1, ~40% cens);
#   b(w) = beta_2 on {w <= 30}, 0 on {w > 30};  beta_2 in {0,...,5}/10;
#   family S(c) = {W <= c}, c in [30, 60];  2000 Monte Carlo samples;
#   r grid {1/3, 1/12, 1/21, 1/30} plus Naive;  NO Adaptive column.
#
# ORIENTATION (resolved 2026-07-21; supersedes any negated-scale framing).
# Section 5.2 pins the direction for the nested family: "the best subgroup,
# S(30), is more distinctive from the others" as beta_2 grows. Under the
# printed DGM the within-S(c) Cox coefficient attains its maximum
# beta(30) = beta_2 at c = 30 and dilutes monotonically toward zero as c
# grows, so selection is the argmax of the RAW coefficient: `orient = +1` in
# every `guohe_algorithm3()` call, the truth curve, the selection, the bounds
# and the coverage indicator all live on the raw coefficient scale, and NO
# negation appears anywhere in this file. Do NOT copy `orient = -1` from
# `gh_one_rep()` in guohe_reproduction_sim.R: Section 5.1's disjoint families
# leave the sign unidentified by the published tables (relabeling symmetry);
# the nesting breaks the symmetry, and under negation selection drifts to
# S(60), inverting the design. Verification V4 below checks the direction.
#
# CANDIDATE FAMILY -- data-determined, never a fixed grid. Per replicate the
# distinct subgroups are indexed by cutpoints {30} union {W_i in (30, 60]}:
# the indicator 1{W_i <= c} changes only as c crosses an observed W_i, so the
# supremum over c in [30, 60] is an EXACT maximum over these order-statistic
# cutpoints. Expected family size 1 + 400 * 3/8 = 151, RANDOM across
# replicates. A fixed ~31-point integer grid is a ~5x smaller family and
# understates selection bias, most visibly in the Naive column (Table 6:
# Naive 0.900 at k = 2 -> 0.543 at k = 12); verification V2 fails on it.
#
# ESTIMAND. gamma_s = beta(c_hat): the truth curve evaluated at the
# replicate's selected cutpoint via `gh52_truth_at()` (default "smooth", the
# isotonic projection). Coverage per column: lower <= gamma_s. The engine's
# selection is r-invariant (it happens on full-sample scores before r
# enters); this is asserted per replicate. The naive comparator is an
# INDEPENDENT calculation (own coxph loop, own argmax, own lower bound),
# scored against gamma_s at its OWN selected cutpoint.
#
# Run:  Rscript quarto/GuoHe/guohe_sec52_sim.R     # verification block + timing
# Full study: guohe_sec52_run.R (requires the truth caches; see there).

suppressMessages(library(survival))

# guohe_algorithm3() is package API, exported by forestsearch at d516a92.
suppressMessages(library(forestsearch))

# DGM, analytic censoring, truth curve and gh52_truth_at() -- path-robust
# whether this file is sourced or executed.
.gh52_sim_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f[1])) else normalizePath(getwd())
  # sourced from a working directory that IS quarto/GuoHe (e.g. the qmd):
  if (!file.exists(file.path(d, "guohe_sec52_truth.R"))) d <- normalizePath(getwd())
  d
})
if (!file.exists(file.path(.gh52_sim_dir, "guohe_sec52_truth.R"))) {
  stop("guohe_sec52_truth.R not found beside guohe_sec52_sim.R nor in the ",
       "working directory (", .gh52_sim_dir, ").", call. = FALSE)
}
source(file.path(.gh52_sim_dir, "guohe_sec52_truth.R"))

# Published grid, Table 7. Same grid as Tables 3-6.
GH52_R_GRID <- c(1 / 3, 1 / 12, 1 / 21, 1 / 30)

# ---- published targets: Table 7 (transcription-verified 2026-07-21) --------
# Empirical coverage of the 95% lower bound of gamma_s. The qmd sources this
# object -- targets are never retyped. Note the proposed columns run
# CONSERVATIVE (0.947-0.973 vs nominal 0.95): adjacent nested candidates
# differ by a single subject, so the effective number of independent
# comparisons is far below the family size. A reproduction landing
# near-nominal is a red flag for the family construction (most likely a
# hard-coded grid).

GH52_TABLE_7 <- data.frame(
  beta2 = c(0, 1, 2, 3, 4, 5) / 10,
  r_1_3 = c(0.947, 0.960, 0.958, 0.959, 0.962, 0.964),
  r_1_12 = c(0.961, 0.972, 0.966, 0.969, 0.968, 0.972),
  r_1_21 = c(0.962, 0.972, 0.967, 0.970, 0.968, 0.973),
  r_1_30 = c(0.962, 0.972, 0.967, 0.970, 0.968, 0.973),
  naive = c(0.872, 0.879, 0.890, 0.895, 0.906, 0.901)
)

# ---- candidate family ------------------------------------------------------

#' Data-determined nested candidate family for one replicate
#'
#' Cutpoints {30} union {W_i : W_i in (30, 60]}, sorted; one 0/1 membership
#' column per cutpoint appended to the data frame. The family size is random
#' across replicates (expected 151 at n = 400) -- by design, never a grid.
#'
#' @param df One replicate from `gh52_sim_data()`.
#' @return List: `df` (with membership columns appended), `cuts` (cutpoint
#'   vector aligned to `names`), `names` (membership column names).
gh52_candidates <- function(df) {
  cuts <- c(GH52_C_LO, sort(df$w[df$w > GH52_C_LO & df$w <= GH52_C_HI]))
  nm <- sprintf("cand_%03d", seq_along(cuts))
  memb <- vapply(cuts, function(cc) as.integer(df$w <= cc),
                 integer(nrow(df)))
  colnames(memb) <- nm
  list(df = cbind(df, as.data.frame(memb)), cuts = cuts, names = nm)
}

# ---- independent per-candidate fits and the naive comparator ---------------
# "select the better subgroup ... and proceed as if the subgroup were selected
# independent of the data". Computed here rather than taken from
# `guohe_algorithm3()` (which does not return standard errors), so the naive
# column is an INDEPENDENT calculation -- the Section 5.1 discipline.

#' Per-candidate raw Cox coefficient and standard error on the full sample
gh52_subgroup_fits <- function(df, cand) {
  k <- length(cand$names)
  est <- se <- rep(NA_real_, k)
  for (j in seq_len(k)) {
    sub <- df[df[[cand$names[j]]] == 1L, , drop = FALSE]
    fit <- try(suppressWarnings(
      survival::coxph(survival::Surv(time, event) ~ treat, data = sub)
    ), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      est[j] <- unname(fit$coefficients[[1L]])
      se[j] <- unname(sqrt(diag(stats::vcov(fit)))[1L])
    }
  }
  list(est = est, se = se)
}

#' Naive selection, one-sided lower bound, and scoring at its own c_hat
gh52_naive <- function(fits, cand, truth, level = 0.05) {
  ok <- is.finite(fits$est)
  sel <- which.max(replace(fits$est, !ok, -Inf))
  c_hat <- cand$cuts[sel]
  gamma_s <- gh52_truth_at(truth, c_hat)
  point <- fits$est[sel]
  lower <- point - stats::qnorm(1 - level) * fits$se[sel]
  list(
    sel = sel, c_hat = c_hat, gamma_s = gamma_s,
    point = point, lower = lower,
    cover = as.integer(lower <= gamma_s),
    dist = gamma_s - lower, bias = point - gamma_s
  )
}

# ---- back-transform: report scale -> raw score scale -----------------------
# Exact inverse of `guohe_algorithm3()`'s `to_report()` at orient = +1 for a
# survival outcome: report = exp(score). No negation (see ORIENTATION above);
# do not reuse `gh_to_score()` from guohe_reproduction_sim.R at its default.
gh52_to_score <- function(report) log(report)

# ---- one replication -------------------------------------------------------

#' One replication: independent fits, naive, and one Guo-He fit per r
#'
#' COMMON RANDOM NUMBERS across the r grid: `guohe_algorithm3()` draws all
#' resample indices under `seed` BEFORE `r` enters, so passing one
#' `boot_seed = seed + 500000L` to every r gives all columns identical
#' bootstrap draws (property verified by execution in the Section 5.1
#' reproduction; it removes between-column Monte Carlo noise, which matters
#' because adjacent Table 7 r columns differ by as little as 0.000-0.001).
#'
#' @param beta2 Step height of `b(.)` on `W <= 30`.
#' @param n Sample size (published: 400).
#' @param truth A `gh52_truth` object for THIS beta2 (production cache).
#' @param r_grid Shrinkage grid (published: 1/3, 1/12, 1/21, 1/30).
#' @param B Bootstrap resamples per engine call (an INFERENCE; paper silent).
#' @param level One-sided level (0.05 -> 95% lower bound).
#' @param min_events Engine estimability guard, passed through.
#' @param seed Replicate seed; also derives the CRN bootstrap seed.
gh52_one_rep <- function(beta2, n = 400L, truth, r_grid = GH52_R_GRID,
                         B = 2000L, level = 0.05, min_events = 5L,
                         seed = NULL) {
  stopifnot(inherits(truth, "gh52_truth"), isTRUE(truth$beta2 == beta2))
  if (!is.null(seed)) set.seed(seed)
  df0 <- gh52_sim_data(beta2, n)
  cand <- gh52_candidates(df0)
  fits <- gh52_subgroup_fits(cand$df, cand)
  nv <- gh52_naive(fits, cand, truth, level = level)

  boot_seed <- if (is.null(seed)) NULL else as.integer(seed) + 500000L

  sel_name <- NA_character_
  c_hat_gh <- gamma_s <- NA_real_
  n_sel <- NA_integer_
  eng <- list()
  for (i in seq_along(r_grid)) {
    fit <- try(suppressWarnings(guohe_algorithm3(
      data = cand$df, outcome = "survival", treatment = "treat",
      candidates = cand$names, time = "time", event = "event",
      orient = +1, B = B, r = r_grid[i], level = level, seed = boot_seed,
      min_events = min_events, diagnostics = FALSE
    )), silent = TRUE)
    tag <- sprintf("r%d", i)
    if (inherits(fit, "try-error")) {
      eng[[paste0(tag, "_cover")]] <- NA_integer_
      eng[[paste0(tag, "_dist")]] <- NA_real_
      eng[[paste0(tag, "_bias")]] <- NA_real_
      next
    }
    if (is.na(sel_name)) {
      # selection is r-invariant (full-sample scores; r enters only through
      # the offsets) -- record once, assert on every later r
      sel_name <- fit$selected
      sel_i <- match(sel_name, cand$names)
      c_hat_gh <- cand$cuts[sel_i]
      n_sel <- sum(cand$df[[sel_name]])
      gamma_s <- gh52_truth_at(truth, c_hat_gh)
    } else if (!identical(fit$selected, sel_name)) {
      stop("Selection not r-invariant within a replicate: ", sel_name,
           " vs ", fit$selected, " at r = ", r_grid[i])
    }
    lower <- gh52_to_score(fit$bound_one_sided)
    debiased <- gh52_to_score(fit$debiased)
    eng[[paste0(tag, "_cover")]] <- as.integer(lower <= gamma_s)
    eng[[paste0(tag, "_dist")]] <- gamma_s - lower
    eng[[paste0(tag, "_bias")]] <- debiased - gamma_s
  }

  cbind(
    data.frame(
      beta2 = beta2, n = n, n_cand = length(cand$cuts),
      cens_rate = mean(df0$event == 0L),
      c_hat_gh = c_hat_gh, c_hat_naive = nv$c_hat,
      sel_agree = as.integer(!is.na(c_hat_gh) && c_hat_gh == nv$c_hat),
      n_sel = n_sel, gamma_s = gamma_s,
      naive_point = nv$point, naive_lower = nv$lower,
      naive_cover = nv$cover, naive_dist = nv$dist, naive_bias = nv$bias,
      gamma_s_naive = nv$gamma_s
    ),
    as.data.frame(eng)
  )
}

# ---- scenario driver -------------------------------------------------------

#' Run one scenario to `n_rep` replications
gh52_run_scenario <- function(beta2, truth, n_rep, n = 400L,
                              r_grid = GH52_R_GRID, B = 2000L, level = 0.05,
                              min_events = 5L, seed = 1L,
                              progress_every = 0L) {
  rows <- vector("list", n_rep)
  for (m in seq_len(n_rep)) {
    rows[[m]] <- gh52_one_rep(beta2, n = n, truth = truth, r_grid = r_grid,
                              B = B, level = level, min_events = min_events,
                              seed = seed + m)
    if (progress_every > 0L && m %% progress_every == 0L) {
      cat(sprintf("    rep %d / %d\n", m, n_rep))
      utils::flush.console()
    }
  }
  do.call(rbind, rows)
}

#' Collapse replications into the published table layout (no Adaptive column)
gh52_summarise <- function(res, r_grid = GH52_R_GRID) {
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
  out$mcse_cover <- sqrt(out$cover * (1 - out$cover) / nrow(res))
  out
}

# ---- verification and timing ----------------------------------------------
# Runs when the script is executed rather than sourced. Fail-loud: every block
# prints PASS/FAIL lines and the script exits non-zero on any failure.

if (sys.nframe() == 0L) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  PF <- function(ok) if (isTRUE(ok)) "PASS" else "FAIL"
  fails <- character(0)
  note <- function(id, ok, msg) {
    cat(sprintf("%s %s -- %s\n", id, PF(ok), msg))
    if (!isTRUE(ok)) fails <<- c(fails, id)
  }

  # A truth object is needed for scoring in V3/V5. Use the production cache if
  # present; otherwise compute a smoke-scale curve (its gates still apply --
  # a gate failure aborts this block, which is the intended fail-loud path).
  get_truth <- function(beta2) {
    f <- file.path(.gh52_sim_dir,
                   sprintf("guohe_sec52_truth_beta2_%02d.rds", round(beta2 * 10)))
    if (file.exists(f)) {
      cat(sprintf("     [truth: production cache %s]\n", basename(f)))
      readRDS(f)
    } else {
      cat(sprintf("     [truth: no cache; computing smoke curve for beta2 = %.1f]\n",
                  beta2))
      gh52_truth_curve(beta2, n_big = 2e5, c_step = 1, seed = 20260721L +
                         as.integer(round(100 * beta2)), verbose = FALSE)
    }
  }

  cat("=== V1: DGM censoring vs analytic mixture ===\n")
  set.seed(20260721)
  for (b2 in c(0, 0.5)) {
    big <- gh52_sim_data(b2, 400000L)
    cr <- mean(big$event == 0L)
    ca <- gh52_censoring_analytic(b2)
    note(sprintf("V1(beta2=%.1f)", b2), abs(cr - ca) < 0.005,
         sprintf("censoring %.4f vs analytic %.4f (paper: 'about 40%%')", cr, ca))
  }

  cat("\n=== V2: family size is data-determined, ~151, and random ===\n")
  set.seed(20260722)
  nc <- vapply(seq_len(500), function(m) {
    length(gh52_candidates(gh52_sim_data(0, 400L))$cuts)
  }, numeric(1))
  note("V2a", abs(mean(nc) - 151) < 3,
       sprintf("mean family size %.2f over 500 draws (target 151 = 1 + 400*3/8)",
               mean(nc)))
  note("V2b", stats::sd(nc) > 1,
       sprintf("family size varies across replicates (sd %.2f); a constant or ~31 means a hard-coded grid",
               stats::sd(nc)))

  cat("\n=== V3: engine agreement and round trip, one replicate ===\n")
  tr0 <- get_truth(0)
  set.seed(11)
  df1 <- gh52_sim_data(0, 400L)
  cand1 <- gh52_candidates(df1)
  fits1 <- gh52_subgroup_fits(cand1$df, cand1)
  fit1 <- suppressWarnings(guohe_algorithm3(
    data = cand1$df, outcome = "survival", treatment = "treat",
    candidates = cand1$names, time = "time", event = "event",
    orient = +1, B = 200L, r = 1 / 30, level = 0.05, seed = 99L,
    diagnostics = FALSE
  ))
  gmax_hand <- max(fits1$est, na.rm = TRUE)
  note("V3a", abs(gh52_to_score(fit1$naive) - gmax_hand) < 1e-10,
       sprintf("gamma_max: engine %.10f vs independent coxph %.10f",
               gh52_to_score(fit1$naive), gmax_hand))
  note("V3b", identical(match(fit1$selected, cand1$names),
                        which.max(replace(fits1$est, !is.finite(fits1$est), -Inf))),
       sprintf("selected %s (c_hat = %.3f); independent argmax cand_%03d",
               fit1$selected,
               cand1$cuts[match(fit1$selected, cand1$names)],
               which.max(replace(fits1$est, !is.finite(fits1$est), -Inf))))
  rt <- gh52_to_score(fit1$bound_one_sided)
  note("V3c", abs(exp(rt) - fit1$bound_one_sided) < 1e-12,
       sprintf("report/score round trip at orient = +1: report %.8f <-> score %.8f",
               fit1$bound_one_sided, rt))

  cat("\n=== V4: orientation direction -- selection concentrates near c = 30 ===\n")
  # Selection needs only the full-sample argmax (r-invariant), so the engine is
  # not required here. Under orient = +1 the selected cutpoint moves TOWARD 30
  # as beta_2 grows ("the best subgroup, S(30)"); concentration near 60 at
  # beta_2 = 0.5 means the Section 5.1 negated convention leaked in.
  sel_c <- function(beta2, n_rep, seed0) {
    vapply(seq_len(n_rep), function(m) {
      set.seed(seed0 + m)
      d <- gh52_sim_data(beta2, 400L)
      cd <- gh52_candidates(d)
      ft <- gh52_subgroup_fits(cd$df, cd)
      cd$cuts[which.max(replace(ft$est, !is.finite(ft$est), -Inf))]
    }, numeric(1))
  }
  c0 <- sel_c(0, 40L, 30000L)
  c5 <- sel_c(0.5, 40L, 40000L)
  cat(sprintf("     c_hat quantiles at beta2 = 0.0: %s\n",
              paste(sprintf("%.1f", stats::quantile(c0, c(.1, .25, .5, .75, .9))),
                    collapse = " / ")))
  cat(sprintf("     c_hat quantiles at beta2 = 0.5: %s\n",
              paste(sprintf("%.1f", stats::quantile(c5, c(.1, .25, .5, .75, .9))),
                    collapse = " / ")))
  note("V4", mean(c5) < mean(c0),
       sprintf("mean c_hat %.2f at beta2 = 0.5 < %.2f at beta2 = 0",
               mean(c5), mean(c0)))

  cat("\n=== V5: timing ===\n")
  time_one <- function(lab, B) {
    t0 <- proc.time()[["elapsed"]]
    invisible(gh52_one_rep(0, n = 400L, truth = tr0, B = B, seed = 7L))
    el <- proc.time()[["elapsed"]] - t0
    cat(sprintf("   %-42s %8.2f s\n", lab, el))
    el
  }
  t_b500 <- time_one("beta2=0  n=400  B=500   all four r", 500L)
  t_b2k <- time_one("beta2=0  n=400  B=2000  all four r", 2000L)
  cat("\n   --- extrapolation to the full study (single core) ---\n")
  hrs <- function(s) s / 3600
  cat(sprintf("   Table 7, 6 beta_2 x 2000 reps, B=500  : %8.1f h\n",
              hrs(6 * 2000 * t_b500)))
  cat(sprintf("   Table 7, 6 beta_2 x 2000 reps, B=2000 : %8.1f h\n",
              hrs(6 * 2000 * t_b2k)))
  cat("   (Production choice of B is an owner decision; run the driver's\n")
  cat("    --pilot mode on the compute host and report the projection.)\n")

  cat("\n")
  if (length(fails)) {
    cat("FAILURES: ", paste(fails, collapse = ", "), "\n")
    quit(status = 1L)
  }
  cat("All verification blocks passed.\n")
}
