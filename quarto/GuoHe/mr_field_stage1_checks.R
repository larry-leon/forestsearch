# mr_field_stage1_checks.R
#
# Stage 1 verification for ci_method = "field" on fs_mr_inference
# (dev/tasks/TASK_mr_field_vs_guohe_2026-09-05.md; PoC provenance:
# dev/tasks/POC_mr_interval_alternatives_2026-09-05.md documents the first
# proof of concept, dev/tasks/poc2_results_2026-09-05.csv is the
# authoritative reference for the values quoted in F5).
#
# Blocks: F0 default-path byte-identity (three fixed-seed cases captured with
# the pre-change installed package); F1 K = 1 identities; F2 perturbation
# covariance vs B'B; F3 seed reproducibility; F4 fast-path equivalence to
# .fs_mr_select; F5 tie sign identity (K = 10 disjoint, 200 replicates, with
# a K = 2 reference row); F6 smoke -- 5 paired replicates each of t35/t6/t7
# with naive identical() to the replication bundles and MR-current identical()
# to the 2026-09-04 comparison bundles; F7 loaded-machine projection at the
# intended worker count (the 2026-09-04 single-core projection missed by 9x).
#
# Fail-loud: PASS/FAIL per block, non-zero exit on failure.
# Run:  Rscript quarto/GuoHe/mr_field_stage1_checks.R [--ref=<path to prechange rds>]

suppressMessages(library(survival))
suppressMessages(library(forestsearch))

.mf_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else normalizePath(getwd())
})
source(file.path(.mf_dir, "mr_vs_guohe_sim.R"))

args <- commandArgs(trailingOnly = TRUE)
ref_opt <- grep("^--ref=", args, value = TRUE)
REF_PATH <- if (length(ref_opt)) sub("^--ref=", "", ref_opt[1]) else
  "/tmp/claude-1000/-home-larryleon-Downloads/d52862e8-776b-4a4d-a050-eeeba8daf6a5/scratchpad/prechange_reference_objects.rds"

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
PF <- function(ok) if (isTRUE(ok)) "PASS" else "FAIL"
fails <- character(0)
note <- function(id, ok, msg) {
  cat(sprintf("%s %s -- %s\n", id, PF(ok), msg))
  if (!isTRUE(ok)) fails <<- c(fails, id)
}

cat("=== F0: default-path byte-identity on three fixed-seed cases ===\n")
# The reference objects were captured with the PRE-change installed package
# (736088e3) using calls identical to the ones below; identical() on the full
# return object is the byte-identity check the task mandates.
stopifnot(file.exists(REF_PATH))
ref <- readRDS(REF_PATH)
dfA <- mv_gh51_regen("t35_beta2_00", 1L); ciA <- mv_cand_idx_51(dfA, 2L)
fA <- gh_subgroup_fits(dfA, c(0, 0))
selA <- which.max(replace(fA$est, !is.finite(fA$est), -Inf))
A2 <- mv_mr(dfA, ciA, names(ciA)[selA], mv_spec51,
            seed = mv_gh51_base("t35_beta2_00") + 1L + MV_SEED_MR)
cB <- mv_gh52_regen("t7_beta2_00", 1L); ciB <- mv_cand_idx_52(cB)
fB <- gh52_subgroup_fits(cB$df, cB)
selB <- which.max(replace(fB$est, !is.finite(fB$est), -Inf))
B2 <- mv_mr(cB$df, ciB, cB$names[selB], mv_spec52,
            seed = mv_gh52_base("t7_beta2_00") + 1L + MV_SEED_MR)
k1 <- list(S1 = which(dfA$S1 == 1L))
C2 <- forestsearch:::fs_mr_inference(
  df = dfA, candidates = k1, spec = mv_spec51, selected_members = k1$S1,
  admission = list(effect_floor = NULL, consistency = NULL),
  reselection = "maxeff", draws = 5000L, multiplier = "poisson",
  ci_method = "ij", seed = 12345L)
# timing_seconds is the wall-clock of the selection loop (proc.time output):
# it differs between ANY two runs, pre-change included, so it is excluded --
# measured on 2026-09-05: the sole differing component was timing_seconds
# (0.086 vs 0.085 s), every substantive field byte-identical.  The exclusion
# is asserted to be the only one: names must match and every other component
# must be identical().
strip_t <- function(x) { x$timing_seconds <- NULL; x }
f0_cmp <- function(new, old) {
  identical(names(new), names(old)) &&
    identical(strip_t(new), strip_t(old)) &&
    is.numeric(new$timing_seconds) && length(new$timing_seconds) == 1L
}
note("F0a", f0_cmp(A2, ref$A),
     "Sec 5.1 t35_beta2_00 m=1: all components identical (timing_seconds excluded, wall-clock)")
note("F0b", f0_cmp(B2, ref$B),
     "Sec 5.2 t7_beta2_00 m=1: all components identical (timing_seconds excluded)")
note("F0c", f0_cmp(C2, ref$C),
     "K = 1: all components identical (timing_seconds excluded)")

cat("\n=== F1: K = 1 field identities ===\n")
# At K = 1 re-selection is degenerate: Lambda*_r = zeta*_r - const, so
# lambda_mean ~ N(0, sigma_D^2 (1/R_out + 1/R_in)) (tol 4 sd), lambda_sd
# estimates sigma_D (band [0.95, 1.05], task-stated), est2 = beta_deb -
# lambda_mean (same 4-sd tol on the log difference), and q95 -> 1.645 sigma_D
# (tol 5%, task-stated).  Gated at R_out = 20,000 so the quantile Monte Carlo
# error (~1.5%) sits inside the 5% tolerance; the production configuration
# (1000/500) is reported alongside, ungated.
k1_field <- function(R_out, R_in) forestsearch:::fs_mr_inference(
  df = dfA, candidates = k1, spec = mv_spec51, selected_members = k1$S1,
  admission = list(effect_floor = NULL, consistency = NULL),
  reselection = "maxeff", draws = 5000L, multiplier = "poisson",
  ci_method = "field", seed = 12345L,
  field_R_out = R_out, field_R_in = R_in)
mrF <- k1_field(20000L, 500L)
sD <- mrF$debiased$se_wald
tol_lm <- 4 * sD * sqrt(1 / 20000 + 1 / 500)
note("F1a", abs(mrF$field$lambda_mean) < tol_lm,
     sprintf("lambda_mean = %+.5f (4-sd tol %.5f)", mrF$field$lambda_mean, tol_lm))
note("F1b", mrF$field$lambda_sd / sD >= 0.95 && mrF$field$lambda_sd / sD <= 1.05,
     sprintf("lambda_sd / sigma_D = %.4f (band [0.95, 1.05])", mrF$field$lambda_sd / sD))
note("F1c", abs(log(mrF$field$est2) - log(mrF$debiased$est)) < tol_lm,
     sprintf("log est2 - log beta_deb = %+.5f (tol %.5f)",
             log(mrF$field$est2) - log(mrF$debiased$est), tol_lm))
note("F1d", abs(mrF$field$q95 / (qnorm(0.95) * sD) - 1) < 0.05,
     sprintf("q95 / (1.645 sigma_D) = %.4f (tol 5%%)", mrF$field$q95 / (qnorm(0.95) * sD)))
mrFp <- k1_field(1000L, 500L)
cat(sprintf("     production config (1000/500): lambda_mean %+.4f, lambda_sd/sigma_D %.3f, q95 ratio %.3f\n",
            mrFp$field$lambda_mean, mrFp$field$lambda_sd / sD,
            mrFp$field$q95 / (qnorm(0.95) * sD)))

cat("\n=== F2: perturbation covariance equals B'B (nested family, R = 20,000) ===\n")
# zeta = B' xi with xi ~ N(0, I) has cov B'B by construction; measured here on
# the Sec 5.2 m=1 family (K = 154, all entries populated).  Deviations are
# scaled by sqrt(Sigma_gg Sigma_hh) (the correlation scale, so exact zeros in
# a disjoint family would not divide by zero); MC sd per entry <= sqrt(2/R) =
# 0.010, max over ~12k entries expected ~ 0.043; tolerance 0.06.
asmB <- mv_assemble(cB$df, ciB, mv_spec52)
SigB <- crossprod(asmB$B)
set.seed(424242)
R2 <- 20000L
Z2 <- crossprod(asmB$B, matrix(rnorm(nrow(asmB$B) * R2), nrow(asmB$B), R2))
Chat <- tcrossprod(Z2) / R2
dscale <- sqrt(diag(SigB))
dev <- abs(Chat - SigB) / tcrossprod(dscale)
note("F2", max(dev) < 0.06,
     sprintf("K = %d, max correlation-scale deviation = %.4f (tol 0.06, MC sd 0.010/entry)",
             nrow(SigB), max(dev)))

cat("\n=== F3: seed reproducibility of the field output ===\n")
B3 <- mv_mr(cB$df, ciB, cB$names[selB], mv_spec52,
            seed = mv_gh52_base("t7_beta2_00") + 1L + MV_SEED_MR, ci_method = "field")
B3b <- mv_mr(cB$df, ciB, cB$names[selB], mv_spec52,
             seed = mv_gh52_base("t7_beta2_00") + 1L + MV_SEED_MR, ci_method = "field")
f3 <- B3$field; f3b <- B3b$field
f3$timing_seconds <- f3b$timing_seconds <- NULL
note("F3", identical(f3, f3b), "two same-seed calls: field element identical (timing excluded)")

cat("\n=== F4: fast-path equivalence to .fs_mr_select (1,000 random vectors) ===\n")
# The field block vectorizes S as which.max / max.col(ties = "first") when
# reselection = "maxeff" with no admission floor.  Equivalence to the gate's
# .fs_mr_select on arbitrary vectors is asserted directly.
set.seed(77)
K4 <- 17L; sz4 <- sample(50:300, K4, replace = TRUE)
ok4 <- TRUE
for (i in 1:1000) {
  v <- rnorm(K4)
  s_ref <- forestsearch:::.fs_mr_select(v, NULL, sz4, seq_len(K4), "maxeff",
                                        0.10, "neighborhood", TRUE)
  if (!identical(as.integer(s_ref), which.max(v))) { ok4 <- FALSE; break }
}
V4 <- matrix(rnorm(K4 * 500), K4, 500)
ok4b <- identical(max.col(t(V4), ties.method = "first"),
                  apply(V4, 2, which.max))
note("F4", ok4 && ok4b,
     "which.max == .fs_mr_select(maxeff, unrestricted); max.col('first') == column-wise which.max")

cat("\n=== F5: tie sign identity (K = 10 disjoint, 200 replicates) ===\n")
# All-null exchangeable family, theta = 0: lambda_mean > 0 on average and
# est2 closer to theta than beta_deb.  PoC-2 reference
# (dev/tasks/poc2_results_2026-09-05.csv): retained bias 0.0634 -> 0.0249 at
# K = 10 (S5), 0.0177 -> 0.0041 at K = 2 (S4), normal-means model.
tie_run <- function(K, ng, n_rep, seed0) {
  rows <- parallel::mclapply(seq_len(n_rep), function(m) {
    set.seed(seed0 + m)
    d <- gh_sim_data(rep(0, K), K * ng)
    d$treat_flip <- 1L - d$treat
    ci <- mv_cand_idx_51(d, K)
    fits <- vapply(ci, function(ix) {
      f <- survival::coxph(survival::Surv(time, event) ~ treat_flip, data = d[ix, ])
      unname(f$coefficients[[1L]])
    }, numeric(1))
    mr <- mv_mr(d, ci, names(ci)[which.max(fits)], mv_spec51,
                seed = seed0 + m + 500000L, ci_method = "field")
    c(bt = log(mr$debiased$est), e2 = log(mr$field$est2),
      lm = mr$field$lambda_mean)
  }, mc.cores = 40L, mc.preschedule = FALSE)
  do.call(rbind, rows)
}
t10 <- tie_run(10L, 200L, 200L, 8100000L)
lm_t <- mean(t10[, "lm"]) / (sd(t10[, "lm"]) / sqrt(nrow(t10)))
note("F5a", mean(t10[, "lm"]) > 0 && lm_t > 3,
     sprintf("K=10: mean lambda_mean = %.4f (t = %.1f > 3)", mean(t10[, "lm"]), lm_t))
note("F5b", abs(mean(t10[, "e2"])) < abs(mean(t10[, "bt"])),
     sprintf("K=10 retained bias: beta_deb %.4f -> est2 %.4f (MCSE %.4f; PoC-2 normal-means ref 0.063 -> 0.025)",
             mean(t10[, "bt"]), mean(t10[, "e2"]), sd(t10[, "e2"]) / sqrt(nrow(t10))))
t2 <- tie_run(2L, 200L, 200L, 8200000L)
cat(sprintf("     K=2 reference row: beta_deb %.4f -> est2 %.4f, mean lambda_mean %.4f (PoC-2 ref 0.018 -> 0.004)\n",
            mean(t2[, "bt"]), mean(t2[, "e2"]), mean(t2[, "lm"])))

cat("\n=== F6: smoke -- 5 paired replicates each of t35_beta2_00, t6_k12, t7_beta2_00 ===\n")
smoke_one <- function(id, m) {
  sec52 <- grepl("^t7_", id)
  rep_bun <- readRDS(file.path(.mf_dir, paste0("guohe_repro_", id, ".rds")))
  cmp_bun <- readRDS(file.path(.mf_dir, paste0("mr_vs_guohe_", id, ".rds")))
  row_r <- rep_bun$results[m, ]; row_c <- cmp_bun$results[m, ]
  base <- if (sec52) mv_gh52_base(id) else mv_gh51_base(id)
  if (sec52) {
    tru <- readRDS(file.path(.mf_dir, sprintf("guohe_sec52_truth_beta2_%02d.rds",
                                              as.integer(sub("^t7_beta2_", "", id)))))
    cand <- mv_gh52_regen(id, m)
    fits <- gh52_subgroup_fits(cand$df, cand)
    nv <- gh52_naive(fits, cand, tru)
    naive_ok <- identical(nv$point, row_r$naive_point) &&
      identical(nv$lower, row_r$naive_lower) && identical(nv$c_hat, row_r$c_hat_naive)
    sel <- which.max(replace(fits$est, !is.finite(fits$est), -Inf))
    mr <- mv_mr(cand$df, mv_cand_idx_52(cand), cand$names[sel], mv_spec52,
                seed = base + m + MV_SEED_MR, ci_method = "field")
  } else {
    sc <- mv_gh51_scenario(id)
    df <- mv_gh51_regen(id, m)
    fits <- gh_subgroup_fits(df, sc$beta)
    nv <- gh_naive(fits, sc$beta)
    naive_ok <- identical(nv$sel, row_r$naive_sel) &&
      identical(nv$beta_s, row_r$naive_beta_s) && identical(nv$dist, row_r$naive_dist)
    mr <- mv_mr(df, mv_cand_idx_51(df, length(sc$beta)), paste0("S", nv$sel),
                mv_spec51, seed = base + m + MV_SEED_MR, ci_method = "field")
  }
  cur_ok <- identical(log(mr$debiased$est), row_c$mr_est) &&
    identical(mr$debiased$se_ij, row_c$mr_se_ij) &&
    identical(mr$selection_bias, row_c$mr_bias_sel) &&
    identical(mr$fixed_bias, row_c$mr_bias_fix)
  fld <- mr$field
  fld_ok <- is.finite(log(fld$est2)) && is.finite(log(fld$lower_1s)) &&
    is.finite(fld$lambda_sd) && fld$n_out_used == fld$R_out
  c(naive_ok = naive_ok, cur_ok = cur_ok, fld_ok = fld_ok,
    t_field = fld$timing_seconds)
}
for (id in c("t35_beta2_00", "t6_k12", "t7_beta2_00")) {
  sm <- t(vapply(1:5, function(m) smoke_one(id, m), numeric(4)))
  note(sprintf("F6(%s)", id),
       all(sm[, "naive_ok"] == 1) && all(sm[, "cur_ok"] == 1) && all(sm[, "fld_ok"] == 1),
       sprintf("naive identical 5/5; MR-current est/se_ij/bias terms identical to 2026-09-04 bundle 5/5; field finite 5/5 (mean field pass %.2f s)",
               mean(sm[, "t_field"])))
}

cat("\n=== F7: loaded-machine projection (intended worker count = 120) ===\n")
# Per-replicate work of the future Stage 2 driver (regeneration + naive
# recompute + one ci_method = "field" gate call + pairing asserts; NO M_eff
# pass, E5), timed under real 120-way mclapply load -- the 1d requirement
# after the 2026-09-04 single-core projection missed by 9x.
proj_cell <- function(id, n_probe) {
  t0 <- proc.time()[["elapsed"]]
  ok <- parallel::mclapply(seq_len(n_probe), function(m) {
    r <- smoke_one(id, m)
    all(r[c("naive_ok", "cur_ok", "fld_ok")] == 1)
  }, mc.cores = 120L, mc.preschedule = FALSE)
  stopifnot(all(unlist(ok)))
  (proc.time()[["elapsed"]] - t0) / n_probe * 120  # amortized core-s per rep
}
cs_t35 <- proj_cell("t35_beta2_00", 360L)
cs_k12 <- proj_cell("t6_k12", 240L)
cs_t7  <- proj_cell("t7_beta2_00", 360L)
interp <- function(k) cs_t35 + (cs_k12 - cs_t35) * (k - 2) / 10
core_h <- (7 * 2000 * cs_t35 + 2000 * (interp(6) + interp(10) + cs_k12) +
           6 * 2000 * cs_t7) / 3600
cat(sprintf("   loaded core-s/rep: t35 %.1f | t6_k12 %.1f | t7 %.1f (probes 360/240/360 at 120 workers)\n",
            cs_t35, cs_k12, cs_t7))
cat(sprintf("   FULL GRID 16 cells x 2000 reps: %.0f core-h -> ~%.0f min wall at 120 workers\n",
            core_h, core_h / 120 * 60))
cat("   (probes reuse Stage 2 seeds; smoke asserts ran inside every probe replicate)\n")

cat("\n")
if (length(fails)) {
  cat("FAILURES: ", paste(fails, collapse = ", "), "\n")
  quit(status = 1L)
}
cat("All Stage 1 field checks passed.\n")
