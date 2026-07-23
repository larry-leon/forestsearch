# guohe_sec52_truth.R
#
# Truth curve beta(c) for the Guo & He (2021) Section 5.2 reproduction
# (Table 7, post-hoc identified case). This file carries the pieces that do
# NOT depend on `forestsearch`: the Section 5.2 data-generating mechanism, the
# analytic censoring probability, the large-sample computation of the truth
# curve beta(c) with its hard verification gates, and the interpolator that
# scores gamma_s = beta(c_hat) at the realized selected cutpoint.
#
# DESIGN, stated in Guo & He (2021) Section 5.2 and verified against the PDF
# text layer on 2026-07-21:
#
#   hazard    : lambda_0(t) e^{b(W) D}, lambda_0 = Weibull(1, 1) (unit expo)
#   D         : Bernoulli(0.5), independent of W;  W ~ Unif[0, 80]
#   censoring : log(C) ~ Unif(-1.25, 1.00), inherited verbatim from Sec 5.1;
#               "about 40% across different choices of b(.)"
#   effect    : b(w) = beta_1 = 0 for w > 30;  b(w) = beta_2 for w <= 30
#   family    : S(c) = {W <= c}, c in [30, 60];  beta(c) = subgroup effect,
#               well defined under misspecification by Lin & Wei (1989),
#               "a weighted average of b(.) in the range [0, c]"
#   scenarios : beta_2 in {0, 1, 2, 3, 4, 5} / 10;  2000 Monte Carlo samples
#
# ORIENTATION (resolved 2026-07-21; supersedes any negated-scale framing).
# Section 2.1 sets larger beta = better effect without loss of generality.
# Section 5.2 pins the direction for the nested family: "the best subgroup,
# S(30), is more distinctive from the others" as beta_2 grows. Under the
# printed DGM the within-S(c) Cox coefficient attains its maximum
# beta(30) = beta_2 at c = 30 and dilutes monotonically toward zero as c
# grows. Selection is therefore the argmax of the RAW coefficient:
# `orient = +1` in guohe_algorithm3(), the truth curve lives on the raw
# coefficient scale, and NO negation appears anywhere in the Section 5.2
# reproduction. (Section 5.1 used orient = -1. For disjoint families the sign
# is unidentified by the published tables -- relabeling symmetry -- so that
# convention reproduced Tables 3-6. The nesting breaks the symmetry and the
# paper's sentence adjudicates. Do not copy `orient = -1` from gh_one_rep().)
#
# STATED vs INFERRED. The DGM, censoring, family, scenario grid and Table 7
# targets are stated in the paper. Inferred / chosen here, recorded as such:
# the numerical method for beta(c) (one large sample per scenario, prefix
# sweep over a fine c grid, Cox fit per grid point), `n_big`, `c_step`, the
# isotonic non-increasing projection used as the interpolation basis where no
# closed form exists, and the gate tolerances (calibration notes at each gate
# below). NOT inferred: at beta2 = 0 the curve beta(c) = 0 is exact, not
# estimated (see `beta_exact` below), so that scenario carries no truth-curve
# Monte Carlo error at all.
#
# Run (production defaults; per-scenario .rds caches, skip-if-exists):
#   Rscript quarto/GuoHe/guohe_sec52_truth.R
# Reduced smoke (container-scale):
#   Rscript quarto/GuoHe/guohe_sec52_truth.R --nbig=200000 --cstep=1

suppressMessages(library(survival))

# ---- constants -------------------------------------------------------------

GH52_W_MAX  <- 80
GH52_C_LO   <- 30
GH52_C_HI   <- 60
# P(W <= 30, D = 1): the elevated-hazard cell of the DGM
GH52_P_ELEV <- (GH52_C_LO / GH52_W_MAX) * 0.5

# ---- data-generating mechanism ---------------------------------------------

#' Generate one sample under the Guo & He Section 5.2 design
#'
#' @param beta2 Scalar; the step height of `b(.)` on `W <= 30`.
#' @param n Sample size.
#' @return Data frame with `time`, `event`, `treat`, and the post-hoc
#'   covariate `w`. Membership columns are NOT materialized here; the
#'   candidate family is data-determined per replicate and belongs to the
#'   simulation harness.
gh52_sim_data <- function(beta2, n) {
  w <- stats::runif(n, 0, GH52_W_MAX)
  d <- stats::rbinom(n, 1L, 0.5)
  b <- ifelse(w <= GH52_C_LO, beta2, 0)
  t_event <- stats::rexp(n, rate = exp(b * d))
  t_cens <- exp(stats::runif(n, -1.25, 1.00))
  data.frame(
    time = pmin(t_event, t_cens),
    event = as.integer(t_event <= t_cens),
    treat = d,
    w = w
  )
}

# ---- analytic censoring -----------------------------------------------------

#' Analytic censoring probability under the Section 5.2 DGM
#'
#' P(cens | rate rho) = E_U[exp(-rho e^U)], U ~ Unif(-1.25, 1). The rate is
#' e^{beta2} on the cell {W <= 30, D = 1} (probability 3/16) and 1 elsewhere.
#' At beta2 = 0 this reduces to the Section 5.1 analytic value (~0.41,
#' the paper's "about 40%").
gh52_censoring_analytic <- function(beta2) {
  i_rho <- function(rho) {
    stats::integrate(function(u) exp(-rho * exp(u)), -1.25, 1.00)$value / 2.25
  }
  (1 - GH52_P_ELEV) * i_rho(1) + GH52_P_ELEV * i_rho(exp(beta2))
}

# ---- truth curve ------------------------------------------------------------

#' Compute the truth curve beta(c) by large-sample simulation, with hard gates
#'
#' beta(c) is the probability limit of the treatment coefficient of the Cox
#' partial-likelihood fit (treatment as sole covariate) on the subpopulation
#' S(c) = {W <= c} -- the Lin & Wei (1989) estimand under misspecification.
#' It is computed here by one large sample from the DGM (INCLUDING the
#' censoring, on which the misspecified estimand depends), a single sort by
#' `w`, and a prefix sweep over the `c` grid.
#'
#' Two exact anchors gate the computation, per the hand-off:
#'   beta(30) = beta2 exactly (b(.) is constant on [0, 30]);
#'   beta(c) is monotone non-increasing for c > 30 (dilution by the null
#'   region only).
#' All gates are evaluated, every failure is reported, and the function then
#' hard-stops -- a failed truth curve must never reach a coverage run.
#'
#' @param beta2 Scalar step height.
#' @param n_big Large-sample size (production default 2e6).
#' @param c_step Grid step over [30, 60] (production default 0.25).
#' @param seed Required; one fixed seed per scenario.
#' @param anchor_k Gate G2 multiple: |beta_hat(30) - beta2| <= anchor_k * SE(30).
#' @param null_k Gate G5 multiple (beta2 = 0 only): max |beta_hat(c)| <=
#'   null_k * max SE(c).
#' @param tol_mono Gate G3 ceiling on any single positive increment of the raw
#'   curve. Structural-error detector: adjacent grid fits share almost all
#'   rows, so genuine Monte Carlo increments are a few 1e-3 at n_big = 2e5 and
#'   shrink like 1/sqrt(n_big); a local jump beyond 0.02 indicates a
#'   subsetting or scale error, not noise.
#' @param tol_total_rise Gate G4 floor component: the end-to-end rise
#'   beta_hat(60) - beta_hat(30) must not exceed max(tol_total_rise,
#'   2 * SE(30)). A gently INCREASING curve is the signature of the negated
#'   (Section 5.1) orientation applied here by mistake; per-step checks alone
#'   cannot catch it because the wrong-direction drift per step is small.
#' @param drop_frac Gate G6 floor (beta2 > 0 only): total decrease
#'   beta_hat(30) - beta_hat(60) >= drop_frac * beta2. Dilution by the null
#'   region makes the true decrease a substantial fraction of beta2 (~0.49 in
#'   the smoke runs, linear in beta2); 0.15 is a wide-margin floor that a sign error (decrease
#'   ~ -0.45 * beta2) or a no-op subset bug (decrease ~ 0) cannot pass.
#' @param cens_tol Gate G1: |empirical - analytic| censoring tolerance.
#' @param verbose Print progress and the gate table.
#' @return Object of class `gh52_truth`: the grids, raw and isotonic-smoothed
#'   curves, per-point SEs and event counts, censoring diagnostics, the gate
#'   table, and bookkeeping (`beta2`, `n_big`, `c_step`, `seed`,
#'   `elapsed_sec`, `r_version`).
gh52_truth_curve <- function(beta2, n_big = 2e6, c_step = 0.25, seed,
                             anchor_k = 4, null_k = 4,
                             tol_mono = 0.02, tol_total_rise = 0.01,
                             drop_frac = 0.15, cens_tol = 0.005,
                             verbose = TRUE) {
  stopifnot(length(beta2) == 1L, beta2 >= 0, n_big >= 1e4, !missing(seed))
  t0 <- proc.time()[["elapsed"]]
  set.seed(seed)
  df <- gh52_sim_data(beta2, n_big)
  cens_emp <- mean(df$event == 0L)
  cens_ana <- gh52_censoring_analytic(beta2)

  # prefix sweep: sort by w once; S(c) is rows 1:findInterval(c, w_sorted)
  df <- df[order(df$w), , drop = FALSE]
  c_grid <- seq(GH52_C_LO, GH52_C_HI, by = c_step)
  n_sub <- findInterval(c_grid, df$w)
  m <- length(c_grid)
  beta_raw <- se <- rep(NA_real_, m)
  n_ev <- integer(m)
  for (j in seq_len(m)) {
    sub <- df[seq_len(n_sub[j]), , drop = FALSE]
    fit <- survival::coxph(survival::Surv(time, event) ~ treat, data = sub,
                           model = FALSE, x = FALSE, y = FALSE)
    beta_raw[j] <- unname(fit$coefficients[[1L]])
    se[j] <- unname(sqrt(diag(stats::vcov(fit)))[1L])
    n_ev[j] <- sum(sub$event)
    if (verbose && (j %% 20L == 0L || j == m)) {
      cat(sprintf("    beta2 = %.1f : grid point %d / %d\r", beta2, j, m))
      utils::flush.console()
    }
  }
  if (verbose) cat("\n")

  # isotonic non-increasing projection: beta(c) is monotone non-increasing on
  # the truth, so the projection strictly reduces Monte Carlo error and is the
  # default interpolation basis where no closed form exists (an INFERRED
  # choice; `gh52_truth_at()` exposes the raw curve as well).
  beta_smooth <- -stats::isoreg(c_grid, -beta_raw)$yf

  # EXACT truth where a closed form exists. At beta2 = 0 the DGM sets
  # b(w) = 0 for every w, so within EVERY S(c) the treatment-only Cox model is
  # correctly specified with true coefficient exactly 0: beta(c) = 0 on the
  # whole grid. Scoring against a simulated estimate of a known constant would
  # inject a FIXED per-scenario offset -- one draw's Monte Carlo error, shared
  # by all coverage replicates and therefore NOT averaging out -- so the exact
  # curve becomes the scoring basis while the simulated curve is retained in
  # full for the gates and for the offset diagnostics below. No closed form
  # exists for beta2 > 0 (the misspecified Lin-Wei estimand is a weighted
  # average of b(.) over [0, c] with censoring-dependent weights), where
  # `beta_exact` is NULL and scoring falls back to `beta_smooth`.
  beta_exact <- if (beta2 == 0) rep(0, m) else NULL

  # Residual truth-curve offset, for the manuscript's inferred ledger. The
  # anchor deviation is observable for every scenario because beta(30) = beta2
  # exactly; where `beta_exact` exists the curve-wide mean deviation is too.
  offset_anchor <- beta_raw[1] - beta2
  offset_mean <- if (is.null(beta_exact)) NA_real_ else mean(beta_raw - beta_exact)

  # ---- gates: evaluate all, report all, then hard-stop on any failure ------
  rise_tol <- max(tol_total_rise, 2 * se[1])
  gates <- data.frame(
    gate = c("G1 censoring", "G2 anchor beta(30)", "G3 max single rise",
             "G4 end-to-end rise",
             if (beta2 == 0) "G5 null band" else "G6 dilution drop"),
    value = c(abs(cens_emp - cens_ana),
              abs(beta_raw[1] - beta2),
              max(diff(beta_raw)),
              beta_raw[m] - beta_raw[1],
              if (beta2 == 0) max(abs(beta_raw))
              else beta_raw[1] - beta_raw[m]),
    threshold = c(cens_tol,
                  anchor_k * se[1],
                  tol_mono,
                  rise_tol,
                  if (beta2 == 0) null_k * max(se) else drop_frac * beta2),
    direction = c("<=", "<=", "<=", "<=", if (beta2 == 0) "<=" else ">="),
    stringsAsFactors = FALSE
  )
  gates$pass <- ifelse(gates$direction == "<=",
                       gates$value <= gates$threshold,
                       gates$value >= gates$threshold)

  if (verbose) {
    cat(sprintf(
      "    beta2 = %.1f : beta_hat(30) = %+.4f (SE %.4f)  drop to c=60 = %+.4f  cens %.4f vs analytic %.4f\n",
      beta2, beta_raw[1], se[1], beta_raw[1] - beta_raw[m], cens_emp, cens_ana))
    for (i in seq_len(nrow(gates))) {
      cat(sprintf("    %-20s %s  value %+.4f  %s %+.4f\n", gates$gate[i],
                  if (gates$pass[i]) "PASS" else "FAIL", gates$value[i],
                  gates$direction[i], gates$threshold[i]))
    }
  }
  if (!all(gates$pass)) {
    stop("Truth curve gates FAILED for beta2 = ", beta2, ": ",
         paste(gates$gate[!gates$pass], collapse = ", "),
         ". A failed truth curve must not reach a coverage run; fix before ",
         "proceeding (see guohe_sec52_truth.R gate documentation).",
         call. = FALSE)
  }

  structure(
    list(
      beta2 = beta2, n_big = n_big, c_step = c_step, seed = seed,
      c_grid = c_grid, beta_raw = beta_raw, beta_smooth = beta_smooth,
      beta_exact = beta_exact,
      scoring_basis = if (is.null(beta_exact)) "smooth" else "exact",
      offset_anchor = offset_anchor, offset_mean = offset_mean,
      se = se, n_sub = n_sub, n_events = n_ev,
      cens_emp = cens_emp, cens_analytic = cens_ana,
      gates = gates,
      elapsed_sec = proc.time()[["elapsed"]] - t0,
      r_version = R.version.string
    ),
    class = "gh52_truth"
  )
}

#' @export
print.gh52_truth <- function(x, ...) {
  cat(sprintf("gh52_truth: beta2 = %.1f, n_big = %g, c in [%d, %d] by %.2f\n",
              x$beta2, x$n_big, GH52_C_LO, GH52_C_HI, x$c_step))
  cat(sprintf("  beta_hat(30) = %+.4f (SE %.4f, target %.4f)\n",
              x$beta_raw[1], x$se[1], x$beta2))
  cat(sprintf("  beta_hat(60) = %+.4f   end-to-end drop %+.4f\n",
              x$beta_raw[length(x$beta_raw)],
              x$beta_raw[1] - x$beta_raw[length(x$beta_raw)]))
  cat(sprintf("  censoring %.4f (analytic %.4f)   gates: %s\n",
              x$cens_emp, x$cens_analytic,
              if (all(x$gates$pass)) "all PASS" else "FAILURES PRESENT"))
  cat(sprintf("  scoring basis: %s%s\n",
              if (is.null(x$beta_exact)) "simulated (isotonic)" else
                "EXACT beta(c) = 0",
              if (is.null(x$beta_exact))
                sprintf("   residual anchor offset %+.4f", x$offset_anchor)
              else ""))
  invisible(x)
}

#' Score the truth curve at the realized selected cutpoint
#'
#' gamma_s = beta(c_hat) for the replicate's selected c_hat, by linear
#' interpolation on the fine grid. Selected cutpoints are order statistics in
#' [30, 60], so `rule = 2` only guards floating-point edge contact.
#'
#' @param truth A `gh52_truth` object.
#' @param c_hat Numeric vector of selected cutpoints.
#' @param use Scoring basis. `"auto"` (default) uses the exact curve when one
#'   exists (beta2 = 0) and the isotonic projection otherwise; this is the
#'   basis the coverage harness should use. `"exact"` demands the closed form
#'   and errors where none exists; `"smooth"` and `"raw"` force the simulated
#'   curve and exist for diagnostics -- including quantifying what the exact
#'   basis removes. Objects written before `beta_exact` existed fall back to
#'   `"smooth"` under `"auto"`, so old caches remain readable.
gh52_truth_at <- function(truth, c_hat, use = c("auto", "exact", "smooth", "raw")) {
  use <- match.arg(use)
  has_exact <- !is.null(truth$beta_exact)
  if (use == "exact" && !has_exact) {
    stop("No exact truth curve exists for beta2 = ", truth$beta2,
         "; a closed form is available only at beta2 = 0.", call. = FALSE)
  }
  y <- switch(
    use,
    auto = if (has_exact) truth$beta_exact else truth$beta_smooth,
    exact = truth$beta_exact,
    smooth = truth$beta_smooth,
    raw = truth$beta_raw
  )
  stats::approx(truth$c_grid, y, xout = c_hat, rule = 2)$y
}

# ---- executable block: gates for all scenarios, caching, timing -------------
# Runs when executed with Rscript, not when sourced. Per-scenario caches
# `guohe_sec52_truth_beta2_XX.rds` are written next to this script (or --out),
# skipped if present unless --force. Exits non-zero on any gate failure.

if (sys.nframe() == 0L) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")

  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grep("^--file=", .a)])
  .dir <- if (length(.f)) dirname(normalizePath(.f[1])) else normalizePath(getwd())
  args <- commandArgs(trailingOnly = TRUE)
  opt <- function(nm, default) {
    hit <- grep(paste0("^--", nm, "="), args, value = TRUE)
    if (!length(hit)) default else sub(paste0("^--", nm, "="), "", hit[1])
  }
  n_big <- as.numeric(opt("nbig", "2e6"))
  c_step <- as.numeric(opt("cstep", "0.25"))
  b2_grid <- as.numeric(strsplit(opt("beta2", "0,0.1,0.2,0.3,0.4,0.5"), ",")[[1]])
  seed_base <- as.integer(opt("seed-base", "20260721"))
  out_dir <- opt("out", .dir)
  force <- any(args == "--force")

  cat(sprintf("Guo & He Sec 5.2 truth curves: n_big = %g, c_step = %.2f, seeds %d + 100*beta2\n\n",
              n_big, c_step, seed_base))
  rows <- list()
  n_fail <- 0L
  for (b2 in b2_grid) {
    id <- sprintf("%02d", round(b2 * 10))
    f_out <- file.path(out_dir, paste0("guohe_sec52_truth_beta2_", id, ".rds"))
    if (file.exists(f_out) && !force) {
      cat(sprintf("[skip] beta2 = %.1f (cache exists: %s)\n", b2, basename(f_out)))
      tr <- readRDS(f_out)
    } else {
      tr <- tryCatch(
        gh52_truth_curve(b2, n_big = n_big, c_step = c_step,
                         seed = seed_base + as.integer(round(100 * b2))),
        error = function(e) e
      )
      if (inherits(tr, "error")) {
        cat(sprintf("[FAIL] beta2 = %.1f : %s\n", b2, conditionMessage(tr)))
        n_fail <- n_fail + 1L
        next
      }
      saveRDS(tr, f_out)
      cat(sprintf("[done] beta2 = %.1f in %.1f s -> %s\n",
                  b2, tr$elapsed_sec, basename(f_out)))
    }
    rows[[length(rows) + 1L]] <- data.frame(
      beta2 = b2, beta30_hat = tr$beta_raw[1], se30 = tr$se[1],
      anchor_z = (tr$beta_raw[1] - tr$beta2) / tr$se[1],
      scoring = tr$scoring_basis,
      offset_anchor = tr$offset_anchor, offset_mean = tr$offset_mean,
      drop_30_60 = tr$beta_raw[1] - tr$beta_raw[length(tr$beta_raw)],
      max_rise = max(diff(tr$beta_raw)),
      cens_emp = tr$cens_emp, cens_analytic = tr$cens_analytic,
      elapsed_sec = tr$elapsed_sec
    )
  }

  if (length(rows)) {
    cat("\n==== summary ====\n")
    print(do.call(rbind, rows), row.names = FALSE, digits = 4)
    tot <- do.call(rbind, rows)
    per_fit <- sum(tot$elapsed_sec) / (nrow(tot) * length(seq(GH52_C_LO, GH52_C_HI, by = c_step)))
    n_prod <- length(seq(GH52_C_LO, GH52_C_HI, by = 0.25))
    cat(sprintf(
      "\n~%.2f s per Cox fit at n_big = %g. Production (n_big = 2e6, c_step = 0.25):\n~%.0f min per scenario serial, x6 scenarios; embarrassingly parallel over\nscenarios and over grid points (fits are independent given the dataset).\n",
      per_fit, n_big, per_fit * (2e6 / n_big) * n_prod / 60))
  }
  if (n_fail > 0L) {
    cat(sprintf("\n%d scenario(s) FAILED truth-curve gates.\n", n_fail))
    quit(status = 1L)
  }
  cat("\nAll requested truth curves computed and gated.\n")
}
