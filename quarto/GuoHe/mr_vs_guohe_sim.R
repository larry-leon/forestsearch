# mr_vs_guohe_sim.R
#
# Stage 1 adapter + identity harness for the MR-vs-Guo&He comparison
# (dev/tasks/TASK_mr_vs_guohe_2026-09-04.md + addendum A).
#
# PURPOSE. Pair the multiplier-resampling (MR) de-biased estimate
# (forestsearch:::fs_mr_inference(), internal) with the committed Guo & He
# replication (Sec 5.1 Tables 3-6, Sec 5.2 Table 7) on IDENTICAL simulated
# replicates, regenerated bit-for-bit from the bundles' stored seed bases.
# This file supplies:
#   * the adapter (replicate regeneration, candidate families as row-index
#     lists, the Sec 5.1 orientation flip, the MR call with pure argmax and
#     no admission screen);
#   * the addendum-A quantities (Sigma-hat, A6 retained mass, M_eff, implied
#     tie residual);
#   * the Stage 1 verification block (identities 1b + A1) and the 1c smoke
#     test, executable and fail-loud.
# The comparison qmd (Stage 3) sources this file; per the Section 5.1/5.2
# convention nothing here is retyped there.
#
# SCOPE. No replication file is modified and no replication compute re-run:
# Guo & He per-replicate values are read from the committed bundles (Sec 5.2
# directly; Sec 5.1 via the selection-recovery arithmetic in Stage 2). The
# only package change this task makes is the D6 add-only argument
# `return_reselection` on fs_mr_inference() (default FALSE, prior return
# object reproduced exactly).
#
# ORIENTATION ADAPTER (Sec 5.1 only). fs_mr_inference() has no `orient`
# argument; its `maxeff` rule maximizes the raw fitted coefficient. Sec 5.2
# selects the argmax of the raw coefficient (orient = +1) and needs no
# adapter. Sec 5.1's inferential effect is the NEGATED log-HR (orient = -1);
# swapping the treatment labels (`treat_flip = 1 - treat`) negates the
# treatment-only Cox coefficient and its dfbeta exactly (the partial
# likelihood is symmetric under the flip), so on the flipped data the raw
# coefficient IS the oriented effect and maxeff is the replication's argmax.
# Verification V6 measures this to 1e-10; the gate's arithmetic is untouched.
#
# MR SETTINGS (task D2/D3 defaults): reselection = "maxeff";
# admission = list(effect_floor = NULL, consistency = NULL) (the engine's
# unrestricted path -- every estimable candidate admissible); draws = 5000;
# multiplier = "poisson" (centred Poisson); ci_method = "ij".
# MR seed convention: seed_mr = seed_base + m + 700000L (disjoint from the
# data seed base+m and the replication's bootstrap seed base+m+500000).
#
# Run:  Rscript quarto/GuoHe/mr_vs_guohe_sim.R    # verification + smoke + timing

suppressMessages(library(survival))
suppressMessages(library(forestsearch))

.mv_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f[1])) else normalizePath(getwd())
  if (!file.exists(file.path(d, "guohe_sec52_sim.R"))) d <- normalizePath(getwd())
  d
})
source(file.path(.mv_dir, "guohe_reproduction_sim.R"))  # gh_sim_data etc. (Sec 5.1)
source(file.path(.mv_dir, "guohe_sec52_sim.R"))         # gh52_* (Sec 5.2; sources truth)

MV_DRAWS      <- 5000L        # D3 default
MV_MULTIPLIER <- "poisson"    # D3 default (centred Poisson)
MV_SEED_MR    <- 700000L      # MR seed offset from the replicate data seed

# ---- replicate regeneration from the bundles' seed scheme ------------------

#' Sec 5.1 scenario parameters from a bundle id ("t35_beta2_03", "t6_k12")
mv_gh51_scenario <- function(id) {
  if (grepl("^t35_beta2_", id)) {
    b2 <- as.integer(sub("^t35_beta2_", "", id)) / 10
    list(beta = c(0, b2), n = 400L)
  } else if (grepl("^t6_k", id)) {
    k <- as.integer(sub("^t6_k", "", id))
    list(beta = rep(0, k), n = 200L * k)
  } else stop("unknown Sec 5.1 scenario id: ", id)
}

#' Seed bases, exactly as in the two run drivers (verified against bundles)
mv_gh51_base <- function(id) 1000000L + as.integer(sum(utf8ToInt(id)) * 997L)
mv_gh52_base <- function(id) 1000000L + as.integer(sum(utf8ToInt(id))) * 100003L

#' Regenerate Sec 5.1 replicate m of a scenario (bit-for-bit: the one-rep
#' function calls set.seed(seed) then draws the data first). Adds the
#' orientation flip column.
mv_gh51_regen <- function(id, m) {
  sc <- mv_gh51_scenario(id)
  set.seed(mv_gh51_base(id) + m)
  df <- gh_sim_data(sc$beta, sc$n)
  df$treat_flip <- 1L - df$treat
  df
}

#' Regenerate Sec 5.2 replicate m (data + order-statistic candidate family)
mv_gh52_regen <- function(id, m) {
  b2 <- as.integer(sub("^t7_beta2_", "", id)) / 10
  set.seed(mv_gh52_base(id) + m)
  df0 <- gh52_sim_data(b2, 400L)
  gh52_candidates(df0)   # list(df, cuts, names)
}

# ---- candidate families as row-index lists (the gate's input format) -------

mv_cand_idx_51 <- function(df, k) {
  stats::setNames(lapply(seq_len(k), function(j) which(df[[paste0("S", j)]] == 1L)),
                  paste0("S", seq_len(k)))
}

mv_cand_idx_52 <- function(cand) {
  stats::setNames(lapply(cand$names, function(nm) which(cand$df[[nm]] == 1L)),
                  cand$names)
}

# ---- specs -----------------------------------------------------------------

mv_spec52 <- list(outcome_type = "survival", effect_measure = "HR",
                  treat.name = "treat", outcome.name = "time",
                  event.name = "event", offset.name = NULL,
                  adjust_covariates = NULL, adverse_outcome = NULL)

# Sec 5.1: identical except the flipped treatment (orientation adapter above)
mv_spec51 <- utils::modifyList(mv_spec52, list(treat.name = "treat_flip"))

# ---- the MR call (pure argmax, no screen; D2/D3 defaults) ------------------

mv_mr <- function(df, cands, sel_label, spec, draws = MV_DRAWS,
                  multiplier = MV_MULTIPLIER, seed = NULL,
                  ci_method = "ij", field_R_out = 1000L, field_R_in = 500L,
                  field_uniform = FALSE) {
  # ci_method = "field" (TASK_mr_field_vs_guohe_2026-09-05, E2 defaults) adds
  # the field element; the debiased element is identical to the "ij" path, so
  # one call yields both MR rows.  Default "ij" keeps the 2026-09-04 behavior.
  args <- list(
    df = df, candidates = cands, spec = spec,
    selected_members = cands[[sel_label]],
    admission = list(effect_floor = NULL, consistency = NULL),
    reselection = "maxeff",
    draws = draws, multiplier = multiplier,
    ci_method = ci_method, seed = seed, return_reselection = TRUE)
  if (identical(ci_method, "field"))
    args <- c(args, list(field_R_out = field_R_out, field_R_in = field_R_in,
                         # kappa(Sigma-hat) sweep (TASK_mr_field_uniform_2026-09-05);
                         # FALSE keeps the 2026-09-05 field output byte-identical.
                         field_uniform = field_uniform))
  do.call(forestsearch:::fs_mr_inference, args)
}

#' Engine-identical influence assembly (read-only use of the internal)
mv_assemble <- function(df, cands, spec) forestsearch:::.fs_mr_assemble(df, cands, spec)

# ---- addendum A quantities -------------------------------------------------

#' Cross-covariance of the candidate effects: Sigma-hat = B_eff' B_eff
mv_sigma <- function(asm) crossprod(asm$B)

#' (A6) retained mass and its standardized form (addendum A1.3)
mv_a6 <- function(Sigma, p_hat, sel) {
  mass <- 2 * sum(p_hat * Sigma[sel, ]) + Sigma[sel, sel]
  list(mass = mass, std = mass / Sigma[sel, sel])
}

#' c(M): expected maximum of M independent standard normals, real M >= 1
mv_c_of_M <- function(M) {
  M * stats::integrate(function(z) z * stats::dnorm(z) * stats::pnorm(z)^(M - 1),
                       -Inf, Inf, rel.tol = 1e-10)$value
}

#' Monte Carlo E[max ζ], ζ ~ N(0, R), chunked; fixed seed; tiny ridge for the
#' Cholesky (the nested Sec 5.2 family is near-singular: adjacent candidates
#' differ by one subject). E[max] is insensitive to the 1e-8 ridge.
mv_m0_hat <- function(R, ndraw = 2e5, seed = 20260904L, ridge = 1e-8,
                      chunk = 5e4) {
  K <- nrow(R)
  L <- chol(R + diag(ridge, K))
  set.seed(seed)
  tot <- 0; tot2 <- 0; done <- 0
  while (done < ndraw) {
    nb <- as.integer(min(chunk, ndraw - done))
    Z <- matrix(stats::rnorm(nb * K), nb, K) %*% L
    mx <- Z[, 1]
    if (K > 1L) for (j in 2:K) mx <- pmax(mx, Z[, j])
    tot <- tot + sum(mx); tot2 <- tot2 + sum(mx^2); done <- done + nb
  }
  m0 <- tot / ndraw
  list(m0 = m0, se = sqrt(max(tot2 / ndraw - m0^2, 0) / ndraw))
}

#' Solve c(M) = m0 for M (effective competition, addendum A1.4)
mv_M_eff <- function(m0, K) {
  if (m0 <= 0) return(1)
  stats::uniroot(function(M) mv_c_of_M(M) - m0, c(1, K + 1),
                 extendInt = "upX", tol = 1e-6)$root
}

mv_R_from_sigma <- function(S) {
  d <- sqrt(diag(S))
  S / tcrossprod(d)
}

#' All addendum-A record columns for one MR result (Stage 2 schema A3)
mv_a_record <- function(asm, mr, ndraw_m0 = 2e5, seed_m0 = 20260904L) {
  Sigma <- mv_sigma(asm)
  sel   <- mr$selected_index
  p_hat <- mr$reselection$p_hat
  a6    <- mv_a6(Sigma, p_hat, sel)
  top3  <- order(p_hat, decreasing = TRUE)[seq_len(min(3L, length(p_hat)))]
  m0    <- mv_m0_hat(mv_R_from_sigma(Sigma), ndraw = ndraw_m0, seed = seed_m0)
  Meff  <- mv_M_eff(m0$m0, nrow(Sigma))
  list(p_hat_H = unname(p_hat[sel]),
       p_hat_top3_labels = names(p_hat)[top3],
       p_hat_top3 = unname(p_hat[top3]),
       Sigma_HH = Sigma[sel, sel],
       A6_mass = a6$mass, A6_mass_std = a6$std,
       m0_hat = m0$m0, m0_mc_se = m0$se, M_eff = Meff,
       tie_resid_implied = (1 - 2^(-0.5)) * m0$m0 * sqrt(Sigma[sel, sel]))
}

# ---- Stage 1 verification, smoke test, timing ------------------------------
# Fail-loud: PASS/FAIL per block, non-zero exit on any failure. Tolerances
# stated inline and printed with the measured values.

if (sys.nframe() == 0L) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  PF <- function(ok) if (isTRUE(ok)) "PASS" else "FAIL"
  fails <- character(0)
  note <- function(id, ok, msg) {
    cat(sprintf("%s %s -- %s\n", id, PF(ok), msg))
    if (!isTRUE(ok)) fails <<- c(fails, id)
  }

  cat("=== V1: influence identities, every candidate (tol 1e-6 relative) ===\n")
  # Sum-to-zero: |sum(dfbeta)| / sigma_D < 1e-6.  Variance: sigma_D^2 =
  # sum(dfbeta^2) vs the Lin-Wei robust (sandwich) variance from an
  # independent coxph(robust = TRUE), relative tolerance 1e-6.
  v1_check <- function(lab, df, cands, spec) {
    asm <- mv_assemble(df, cands, spec)
    stopifnot(length(asm$names) == length(cands))  # nothing dropped
    worst_sum <- 0; worst_var <- 0
    for (g in seq_along(asm$names)) {
      idx <- cands[[asm$names[g]]]
      db  <- asm$B[idx, g]
      worst_sum <- max(worst_sum, abs(sum(db)) / asm$sigma_D[g])
      rob <- survival::coxph(
        stats::reformulate(spec$treat.name,
          response = sprintf("survival::Surv(%s, %s)", spec$outcome.name, spec$event.name)),
        data = df[idx, , drop = FALSE], robust = TRUE)
      vr <- unname(rob$var[1, 1])
      worst_var <- max(worst_var, abs(asm$sigma_D[g]^2 - vr) / vr)
    }
    note(paste0("V1sum(", lab, ")"), worst_sum < 1e-6,
         sprintf("%d candidates, worst |sum(dfbeta)|/sigma_D = %.2e", length(cands), worst_sum))
    note(paste0("V1var(", lab, ")"), worst_var < 1e-6,
         sprintf("worst |sigma_D^2 - robust var| / robust var = %.2e", worst_var))
  }
  c52 <- mv_gh52_regen("t7_beta2_00", 1L)
  v1_check("Sec5.2 t7_b2_00 m=1", c52$df, mv_cand_idx_52(c52), mv_spec52)
  d51 <- mv_gh51_regen("t6_k12", 1L)
  v1_check("Sec5.1 t6_k12 m=1", d51, mv_cand_idx_51(d51, 12L), mv_spec51)

  cat("\n=== V2: K = 1 -- de-biased equals naive; V/sigma_D^2 in [3.6, 4.4] ===\n")
  # With one candidate the winner is that candidate on every draw, so
  # debiased - naive = -2 * mean_b D(b), sd = 2 sigma_D / sqrt(B): tolerance
  # 4 sd (fixed seed, deterministic).  The IJ variance in the same regime is
  # the separated-regime limit 4 sigma_D^2 minus the finite-B correction
  # (N/B) * 4 sigma_D^2; the task band is [3.6, 4.4].
  df1 <- mv_gh51_regen("t35_beta2_00", 1L)
  k1  <- list(S1 = which(df1$S1 == 1L))
  mr1 <- mv_mr(df1, k1, "S1", mv_spec51, draws = 5000L, seed = 12345L)
  d_log <- log(mr1$debiased$est) - log(mr1$naive$est)
  tolK1 <- 4 * 2 * mr1$debiased$se_wald / sqrt(5000)
  note("V2a", abs(d_log) < tolK1,
       sprintf("debiased - naive (log) = %+.5f, tol (4 sd) = %.5f", d_log, tolK1))
  ratio <- mr1$debiased$var_ij / mr1$debiased$se_wald^2
  note("V2b", ratio >= 3.6 && ratio <= 4.4,
       sprintf("V_ij / sigma_D^2 = %.3f (band [3.6, 4.4]; ij_source = %s)",
               ratio, mr1$debiased$ij_source))

  cat("\n=== V3: exchangeable candidates -> exchangeable selection (chi-sq p > 0.01) ===\n")
  # 'All K candidates identical' read as identical IN DISTRIBUTION (K = 5
  # disjoint, equal probability, all beta = 0): within one replicate the
  # re-selection legitimately concentrates on the realized argmax, so
  # exchangeability is tested ACROSS replicates -- the modal MR winner over
  # 200 independent replicates must be uniform on the 5 labels.
  # (Literally identical membership vectors give identical perturbed effects
  # on every draw and which.max resolves ties to the first index
  # deterministically; that tests tie-breaking, not exchangeability.)
  nrep3 <- 200L
  modal <- integer(nrep3)
  for (m in seq_len(nrep3)) {
    set.seed(610000L + m)
    dfe <- gh_sim_data(rep(0, 5), 500L)
    dfe$treat_flip <- 1L - dfe$treat
    ce  <- mv_cand_idx_51(dfe, 5L)
    est <- vapply(ce, function(ix) {
      f <- survival::coxph(survival::Surv(time, event) ~ treat_flip, data = dfe[ix, ])
      unname(f$coefficients[[1L]])
    }, numeric(1))
    mre <- mv_mr(dfe, ce, names(ce)[which.max(est)], mv_spec51,
                 draws = 200L, seed = 620000L + m)
    modal[m] <- which.max(mre$reselection$p_hat)
  }
  cnt <- tabulate(modal, nbins = 5L)
  p3  <- stats::chisq.test(cnt, p = rep(1 / 5, 5))$p.value
  note("V3", p3 > 0.01,
       sprintf("modal-winner counts %s over %d replicates; chi-sq p = %.3f",
               paste(cnt, collapse = "/"), nrep3, p3))

  cat("\n=== V4: M_eff exchangeable targets (addendum A1.4; tol 0.1) ===\n")
  for (tv in list(c(0.3, 6.22), c(0.6, 3.65), c(0.9, 1.80))) {
    Rho <- matrix(tv[1], 10, 10); diag(Rho) <- 1
    m0  <- mv_m0_hat(Rho, ndraw = 2e5, seed = 20260904L)
    m0a <- sqrt(1 - tv[1]) * mv_c_of_M(10)   # exchangeable closed form
    Me  <- mv_M_eff(m0$m0, 10L)
    note(sprintf("V4(rho=%.1f)", tv[1]),
         abs(Me - tv[2]) < 0.1 && abs(m0$m0 - m0a) < 4 * m0$se,
         sprintf("m0 = %.4f (MC se %.4f; analytic %.4f), M_eff = %.2f (target %.2f)",
                 m0$m0, m0$se, m0a, Me, tv[2]))
  }

  cat("\n=== V5: disjoint family (R = I, K = 10) -> M_eff = K (tol 0.2) ===\n")
  m0i <- mv_m0_hat(diag(10), ndraw = 2e5, seed = 20260905L)
  Mei <- mv_M_eff(m0i$m0, 10L)
  note("V5", abs(Mei - 10) < 0.2,
       sprintf("m0 = %.4f (c(10) = %.4f), M_eff = %.2f", m0i$m0, mv_c_of_M(10), Mei))

  cat("\n=== V6: treatment-flip orientation adapter (rel tol 1e-10) ===\n")
  # Flip negates every beta-hat(g) and every column of B_eff; Sigma-hat is
  # invariant.
  dfo <- mv_gh51_regen("t35_beta2_03", 2L)
  co  <- mv_cand_idx_51(dfo, 2L)
  asm_raw  <- mv_assemble(dfo, co, utils::modifyList(mv_spec51, list(treat.name = "treat")))
  asm_flip <- mv_assemble(dfo, co, mv_spec51)
  rel_b <- max(abs(asm_flip$beta_hat + asm_raw$beta_hat)) / max(abs(asm_raw$beta_hat))
  rel_B <- max(abs(asm_flip$B + asm_raw$B)) / max(abs(asm_raw$B))
  rel_S <- max(abs(mv_sigma(asm_flip) - mv_sigma(asm_raw))) / max(abs(mv_sigma(asm_raw)))
  note("V6", rel_b < 1e-10 && rel_B < 1e-10 && rel_S < 1e-10,
       sprintf("rel err: beta %.1e, B_eff columns %.1e, Sigma %.1e", rel_b, rel_B, rel_S))

  cat("\n=== V7: (A6) retained mass on both family types ===\n")
  # Sec 5.1 disjoint: off-diagonal Sigma-hat is EXACTLY zero (disjoint
  # dfbeta supports), so A6_mass_std = 1 + 2 p_hat_H in [1, 3].
  asm71 <- mv_assemble(d51, mv_cand_idx_51(d51, 12L), mv_spec51)
  S71   <- mv_sigma(asm71)
  offmax <- max(abs(S71[upper.tri(S71)]))
  est71 <- asm71$beta_hat
  mr71  <- mv_mr(d51, mv_cand_idx_51(d51, 12L),
                 asm71$names[which.max(est71)], mv_spec51,
                 draws = 2000L, seed = 710001L)
  a71 <- mv_a6(S71, mr71$reselection$p_hat, mr71$selected_index)
  note("V7a", offmax == 0 &&
         abs(a71$std - (1 + 2 * mr71$reselection$p_hat[mr71$selected_index])) < 1e-12,
       sprintf("Sec 5.1 disjoint: max |offdiag| = %g; A6_std = %.3f = 1 + 2 p_hat_H", offmax, a71$std))
  # Sec 5.2 nested: structural claim is non-negative cross-terms with the
  # winner, hence A6_mass_std >= 1.  Measured, not assumed.
  asm72 <- mv_assemble(c52$df, mv_cand_idx_52(c52), mv_spec52)
  S72   <- mv_sigma(asm72)
  mr72  <- mv_mr(c52$df, mv_cand_idx_52(c52),
                 asm72$names[which.max(asm72$beta_hat)], mv_spec52,
                 draws = 2000L, seed = 720001L)
  a72 <- mv_a6(S72, mr72$reselection$p_hat, mr72$selected_index)
  note("V7b", a72$std >= 1,
       sprintf("Sec 5.2 nested: A6_std = %.3f (>= 1); min cross-term with winner = %+.3e",
               a72$std, min(S72[mr72$selected_index, ])))

  cat("\n=== V8: D6 return contract ===\n")
  mr_off <- forestsearch:::fs_mr_inference(
    df = df1, candidates = k1, spec = mv_spec51, selected_members = k1$S1,
    admission = list(effect_floor = NULL, consistency = NULL),
    reselection = "maxeff", draws = 50L, multiplier = "poisson",
    ci_method = "ij", seed = 5L)
  exp_names <- c("selected_index", "selected_label", "measure", "log_scale",
                 "ci_method", "naive", "debiased", "selection_bias",
                 "fixed_bias", "selection_rate", "mean_r", "mean_r_c",
                 "complement", "settings", "harm_flag", "n_family",
                 "n_selected", "timing_seconds")
  note("V8a", identical(names(mr_off), exp_names),
       "default return names unchanged (no 'reselection' element)")
  note("V8b", length(mr72$reselection$winner) == 2000L &&
         abs(sum(mr72$reselection$p_hat) - mr72$selection_rate) < 1e-12,
       sprintf("winner length = draws; sum(p_hat) = %.6f = selection_rate = %.6f",
               sum(mr72$reselection$p_hat), mr72$selection_rate))

  cat("\n=== S1: smoke test -- t7_beta2_00, replicates 1..5, paired by seed ===\n")
  bun <- readRDS(file.path(.mv_dir, "guohe_repro_t7_beta2_00.rds"))
  tru <- readRDS(file.path(.mv_dir, "guohe_sec52_truth_beta2_00.rds"))
  stopifnot(bun$seed_base == mv_gh52_base("t7_beta2_00"))
  t_mr <- numeric(5)
  ok_naive <- ok_sel <- ok_mr <- logical(5)
  smoke <- vector("list", 5)
  for (m in 1:5) {
    t0 <- proc.time()[["elapsed"]]
    cand <- mv_gh52_regen("t7_beta2_00", m)
    fits <- gh52_subgroup_fits(cand$df, cand)
    nv   <- gh52_naive(fits, cand, tru)
    row  <- bun$results[m, ]
    ok_naive[m] <- identical(nv$point, row$naive_point) &&
      identical(nv$lower, row$naive_lower) && identical(nv$c_hat, row$c_hat_naive) &&
      identical(nv$gamma_s, row$gamma_s_naive)
    sel_lab <- cand$names[which.max(replace(fits$est, !is.finite(fits$est), -Inf))]
    ok_sel[m] <- identical(cand$cuts[match(sel_lab, cand$names)], row$c_hat_gh)
    ci  <- mv_cand_idx_52(cand)
    mr  <- mv_mr(cand$df, ci, sel_lab, mv_spec52,
                 seed = bun$seed_base + m + MV_SEED_MR)
    t_mr[m] <- proc.time()[["elapsed"]] - t0
    ok_mr[m] <- is.finite(mr$debiased$est) && is.finite(mr$debiased$se) &&
      identical(mr$selected_label, sel_lab)
    # Guo & He per-replicate values, from the STORED columns (no re-run):
    # debiased = r*_bias + gamma_s ; lower = gamma_s - r*_dist  (r4 = 1/30)
    smoke[[m]] <- data.frame(
      m = m, gamma_s = row$gamma_s,
      naive = log(mr$naive$est),
      gh_debiased_r4 = row$r4_bias + row$gamma_s,
      gh_lower_r4 = row$gamma_s - row$r4_dist,
      mr_debiased = log(mr$debiased$est),
      mr_lower_1s = log(mr$debiased$lower_1s),
      mr_se_ij = mr$debiased$se_ij, mr_bias_sel = mr$selection_bias,
      mr_bias_fix = mr$fixed_bias, mr_ij_source = mr$debiased$ij_source,
      mr_p_hat_H = unname(mr$reselection$p_hat[mr$selected_index]),
      t_mr_s = t_mr[m])
  }
  note("S1a", all(ok_naive),
       "recomputed naive point/lower/c_hat/gamma_s identical() to stored, all 5 replicates")
  note("S1b", all(ok_sel), "recomputed argmax cutpoint equals stored c_hat_gh, all 5")
  note("S1c", all(ok_mr), "MR ran with finite estimate and SE on all 5")
  sm <- do.call(rbind, smoke)
  print(sm, digits = 4)
  cat("\n  Addendum-A record columns, replicate 1 (schema viability):\n")
  cand1 <- mv_gh52_regen("t7_beta2_00", 1L)
  asm1  <- mv_assemble(cand1$df, mv_cand_idx_52(cand1), mv_spec52)
  fits1 <- gh52_subgroup_fits(cand1$df, cand1)
  sel1  <- cand1$names[which.max(replace(fits1$est, !is.finite(fits1$est), -Inf))]
  mrA   <- mv_mr(cand1$df, mv_cand_idx_52(cand1), sel1, mv_spec52,
                 seed = bun$seed_base + 1L + MV_SEED_MR)
  arec  <- mv_a_record(asm1, mrA)
  cat(sprintf("  p_hat_H %.3f | top3 %s (%s) | Sigma_HH %.5f | A6_mass %.5f (std %.3f)\n",
              arec$p_hat_H, paste(sprintf("%.3f", arec$p_hat_top3), collapse = "/"),
              paste(arec$p_hat_top3_labels, collapse = "/"),
              arec$Sigma_HH, arec$A6_mass, arec$A6_mass_std))
  cat(sprintf("  m0_hat %.4f (mc se %.5f) | M_eff %.2f of K = %d | tie_resid_implied %.4f\n",
              arec$m0_hat, arec$m0_mc_se, arec$M_eff, nrow(mv_sigma(asm1)),
              arec$tie_resid_implied))

  cat("\n=== T: timing and full-grid projection ===\n")
  # Per-replicate Stage 2 cost = regen + naive + MR + the addendum-A record
  # (whose M_eff Monte Carlo at 2e5 draws dominates on the ~151-candidate
  # Sec 5.2 family).  Both components measured and printed.
  time_51 <- function(id, m) {
    t0 <- proc.time()[["elapsed"]]
    df <- mv_gh51_regen(id, m)
    k  <- length(mv_gh51_scenario(id)$beta)
    ci <- mv_cand_idx_51(df, k)
    asm <- mv_assemble(df, ci, mv_spec51)
    mr <- mv_mr(df, ci, asm$names[which.max(asm$beta_hat)], mv_spec51,
                seed = mv_gh51_base(id) + m + MV_SEED_MR)
    stopifnot(is.finite(mr$debiased$est))
    t1 <- proc.time()[["elapsed"]]
    invisible(mv_a_record(asm, mr))
    c(mr = t1 - t0, arec = proc.time()[["elapsed"]] - t1)
  }
  t35v  <- rowMeans(vapply(1:3, function(m) time_51("t35_beta2_00", m), numeric(2)))
  tk12v <- rowMeans(vapply(1:2, function(m) time_51("t6_k12", m), numeric(2)))
  # Sec 5.2: MR pass measured in the smoke block; A-record measured here
  t0 <- proc.time()[["elapsed"]]
  invisible(mv_a_record(asm1, mrA))
  t52_arec <- proc.time()[["elapsed"]] - t0
  t35 <- sum(t35v); tk12 <- sum(tk12v); t52 <- mean(t_mr) + t52_arec
  cat(sprintf("   per-replicate (regen + naive + MR at B = %d + A-record), measured:\n", MV_DRAWS))
  cat(sprintf("     t35 (k=2, n=400)     %6.2f s  (MR %5.2f + A %5.2f; 3 reps)\n",
              t35, t35v[1], t35v[2]))
  cat(sprintf("     t6_k12 (k=12,n=2400) %6.2f s  (MR %5.2f + A %5.2f; 2 reps)\n",
              tk12, tk12v[1], tk12v[2]))
  cat(sprintf("     t7 (~151 cands,n=400)%6.2f s  (MR %5.2f + A %5.2f; 5 + 1 reps)\n",
              t52, mean(t_mr), t52_arec))
  # t6 k = 6, 10 interpolated linearly in k between the measured k = 2 and 12
  interp <- function(k) t35 + (tk12 - t35) * (k - 2) / 10
  core_h <- (7 * 2000 * t35 + 2000 * (interp(6) + interp(10) + tk12) +
             6 * 2000 * t52) / 3600
  workers <- 120
  cat(sprintf("   full grid, 16 cells x 2000 reps: %.1f core-h -> ~%.1f min wall at %d workers\n",
              core_h, core_h / workers * 60, workers))

  cat("\n")
  if (length(fails)) {
    cat("FAILURES: ", paste(fails, collapse = ", "), "\n")
    quit(status = 1L)
  }
  cat("All Stage 1 verification blocks passed.\n")
}
