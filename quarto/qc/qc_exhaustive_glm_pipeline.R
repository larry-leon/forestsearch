#!/usr/bin/env Rscript
# =============================================================================
# EXHAUSTIVE QC: GLM Simulation Pipeline Pressure Testing
#
# Self-contained: requires only base R + survival.
# Sources generate_glm_dgm.R, simulate_from_glm_dgm.R,
# calibrate_glm_interaction.R, and glm_effect_estimators.R directly.
#
# 12 Sections, 80+ checks covering:
#   1. Estimator sign conventions (all types, both directions)
#   2. DGM potential outcomes correctness
#   3. Calibration accuracy & monotonicity
#   4. Simulation fidelity (DGM → simulated data agreement)
#   5. Edge cases: subgroup sizes (tiny, large, 1-var, 3-var)
#   6. k_inter extremes (zero, huge positive, huge negative)
#   7. Offset edge cases (zero, negative, constant, tiny)
#   8. Input validation (wrong types, missing args)
#   9. Reproducibility (seed stability)
#  10. print.glm_dgm method
#  11. Cross-outcome type isolation (no leaking columns)
#  12. Null model across all types (k_inter forced to 0)
# =============================================================================

library(survival)

cat("=============================================================\n")
cat("  EXHAUSTIVE QC: GLM DGM / Calibration / Simulation\n")
cat("  R", R.version.string, "\n")
cat("=============================================================\n\n")

# ─── Source pipeline ─────────────────────────────────────────────────────
# Adjust these paths to your environment. With devtools::load_all() in
# RStudio, these are not needed — the package namespace handles it.
# For standalone execution:
# src <- function(f) {
#   for (dir in c(".", "R", "latest_R/R")) {
#     fp <- file.path(dir, f)
#     if (file.exists(fp)) { source(fp, local = FALSE); return(invisible()) }
#   }
#   stop("Cannot find: ", f)
# }
# src("globals.R")
# src("generate_aft_dgm_helpers.R")
# src("glm_effect_estimators.R")
# src("generate_glm_dgm.R")
# src("simulate_from_glm_dgm.R")
# src("calibrate_glm_interaction.R")


library(forestsearch)

pass <- 0L; fail <- 0L; warn_count <- 0L

check <- function(label, cond, msg = "") {
  if (isTRUE(cond)) {
    pass <<- pass + 1L; cat(sprintf("  [PASS] %s\n", label))
  } else {
    fail <<- fail + 1L; cat(sprintf("  [FAIL] %s %s\n", label, msg))
  }
}

# ─── Shared synthetic data ───────────────────────────────────────────────
set.seed(8316951L)
N <- 1000
syn <- data.frame(
  id = seq_len(N), treat = rep(0:1, each = N / 2),
  z1 = as.factor(sample(0:1, N, TRUE)),
  z2 = as.factor(sample(0:1, N, TRUE)),
  z3 = as.factor(sample(0:1, N, TRUE)),
  z4 = as.factor(sample(0:1, N, TRUE))
)
p_base <- plogis(-0.5 + 0.3 * syn$treat)
syn$y_bin   <- rbinom(N, 1, p_base)
syn$y_cont  <- 50 + 5 * syn$treat + rnorm(N, sd = 15)
syn$y_count <- rpois(N, exp(1.0 + 0.3 * syn$treat))
syn$follow_up <- runif(N, 0.5, 5.0)
syn$y_events  <- rpois(N, exp(-0.5 + 0.3 * syn$treat) * syn$follow_up)

cat(sprintf("Synthetic: N=%d, Q(z1&z2)=%.1f%%\n\n",
            N, 100 * mean(syn$z1 == 1 & syn$z2 == 1)))


# ═══════════════════════════════════════════════════════════════════════════
cat("─── 1. Estimator Sign Conventions ─────────────────────────\n")
# ═══════════════════════════════════════════════════════════════════════════

set.seed(42)
df_p <- data.frame(treat = rep(0:1, each = 300),
                   y = c(rbinom(300, 1, .3), rbinom(300, 1, .7)))
df_n <- data.frame(treat = rep(0:1, each = 300),
                   y = c(rbinom(300, 1, .7), rbinom(300, 1, .3)))

est_or <- make_effect_estimator("binary", "treat", "y", effect_measure = "OR")
check("OR: treat↑ → log(OR)>0", est_or(df_p)$estimate > 0)
check("OR: treat↓ → log(OR)<0", est_or(df_n)$estimate < 0)
check("OR: SE > 0", est_or(df_p)$se > 0)

est_rd <- make_effect_estimator("binary", "treat", "y", effect_measure = "RD")
check("RD: treat↑ → RD>0", est_rd(df_p)$estimate > 0)

df_cp <- data.frame(treat = rep(0:1, each = 300),
                    y = c(rnorm(300, 50, 10), rnorm(300, 80, 10)))
df_cn <- data.frame(treat = rep(0:1, each = 300),
                    y = c(rnorm(300, 80, 10), rnorm(300, 50, 10)))
est_md <- make_effect_estimator("continuous", "treat", "y", effect_measure = "MD")
check("MD: treat↑ → MD>0", est_md(df_cp)$estimate > 20)
check("MD: treat↓ → MD<0", est_md(df_cn)$estimate < -20)
check("MD: SE > 0", est_md(df_cp)$se > 0)

df_ct1 <- data.frame(treat = rep(0:1, each = 300),
                     y = c(rpois(300, 2), rpois(300, 5)), fu = 1)
df_ct2 <- data.frame(treat = rep(0:1, each = 300),
                     y = c(rpois(300, 5), rpois(300, 2)), fu = 1)
est_irr <- make_effect_estimator("count", "treat", "y",
                                  offset.name = "fu", effect_measure = "IRR")
check("IRR: treat↑ → log(IRR)>0", est_irr(df_ct1)$estimate > 0)
check("IRR: treat↓ → log(IRR)<0", est_irr(df_ct2)$estimate < 0)

# IRR with large offset — same IRR regardless
df_off5 <- data.frame(treat = rep(0:1, each = 300),
                      y = c(rpois(300, 2 * 5), rpois(300, 5 * 5)), fu = 5)
est_off5 <- make_effect_estimator("count", "treat", "y",
                                   offset.name = "fu", effect_measure = "IRR")
check("IRR+offset=5: IRR≈2.5",
      abs(exp(est_off5(df_off5)$estimate) - 2.5) < 1,
      sprintf("got %.2f", exp(est_off5(df_off5)$estimate)))

cat("\n")


# ═══════════════════════════════════════════════════════════════════════════
cat("─── 2. DGM Potential Outcomes ─────────────────────────────\n")
# ═══════════════════════════════════════════════════════════════════════════

# Binary null
dgm_bn <- generate_glm_dgm(
  syn, c("z1", "z2", "z3"), "y_bin", "treat", "binary",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "null", n_super = 3000L, seed = 42)
check("Bin null: p0,p1 exist",
      all(c("p0", "p1") %in% names(dgm_bn$df_super)))
check("Bin null: OR(Q)≈OR(Qc)",
      abs(dgm_bn$hazard_ratios$harm_subgroup -
          dgm_bn$hazard_ratios$no_harm_subgroup) < 0.1)

# Binary alt — positive k_inter
dgm_ba <- generate_glm_dgm(
  syn, c("z1", "z2", "z3"), "y_bin", "treat", "binary",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "alt", k_inter = 1.5, n_super = 3000L, seed = 42)
qi <- dgm_ba$df_super$flag_harm == 1
check("Bin alt: mean(p1)>mean(p0) in Q",
      mean(dgm_ba$df_super$p1[qi]) > mean(dgm_ba$df_super$p0[qi]))
check("Bin alt: OR(Q)>OR(Qc)",
      dgm_ba$hazard_ratios$harm_subgroup > dgm_ba$hazard_ratios$no_harm_subgroup)

# Continuous alt
dgm_ca <- generate_glm_dgm(
  syn, c("z1", "z2", "z3"), "y_cont", "treat", "continuous", "MD",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "alt", k_inter = 20, n_super = 3000L, seed = 42)
qc <- dgm_ca$df_super$flag_harm == 1
check("Cont alt: mu0,mu1 exist",
      all(c("mu0", "mu1") %in% names(dgm_ca$df_super)))
check("Cont alt: mu1>mu0 in Q",
      mean(dgm_ca$df_super$mu1[qc]) > mean(dgm_ca$df_super$mu0[qc]))
check("Cont alt: MD(Q)>MD(Qc)",
      dgm_ca$hazard_ratios$harm_subgroup > dgm_ca$hazard_ratios$no_harm_subgroup)
check("Cont alt: sigma stored",
      !is.null(dgm_ca$model_params$sigma) && dgm_ca$model_params$sigma > 0)

# Count alt
dgm_ct <- generate_glm_dgm(
  syn, c("z1", "z2", "z3"), "y_count", "treat", "count", "IRR",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "alt", k_inter = 0.5, n_super = 3000L, seed = 42)
check("Count alt: all mu>0",
      all(dgm_ct$df_super$mu0 > 0) && all(dgm_ct$df_super$mu1 > 0))
check("Count alt: IRR(Q)>IRR(Qc)",
      dgm_ct$hazard_ratios$harm_subgroup > dgm_ct$hazard_ratios$no_harm_subgroup)

# Count with offset
dgm_off <- generate_glm_dgm(
  syn, c("z1", "z2", "z3"), "y_events", "treat", "count", "IRR",
  offset_var = "follow_up",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "alt", k_inter = 0.5, n_super = 3000L, seed = 42)
check("Count+offset: offset_var stored",
      dgm_off$model_params$offset_var == "follow_up")
check("Count+offset: mu0 scales with FU",
      cor(dgm_off$df_super$mu0, dgm_off$df_super$follow_up) > 0.3)

# Negative k_inter reverses direction
dgm_neg <- generate_glm_dgm(
  syn, c("z1", "z2", "z3"), "y_bin", "treat", "binary",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "alt", k_inter = -1.5, n_super = 3000L, seed = 42)
check("Neg k_inter: OR(Q)<OR(Qc)",
      dgm_neg$hazard_ratios$harm_subgroup < dgm_neg$hazard_ratios$no_harm_subgroup)

cat("\n")


# ═══════════════════════════════════════════════════════════════════════════
cat("─── 3. Calibration Accuracy & Monotonicity ────────────────\n")
# ═══════════════════════════════════════════════════════════════════════════

dgm_cb <- calibrate_glm_interaction(
  syn, c("z1", "z2", "z3"), "y_bin", "treat", 2.0,
  "binary", "OR", subgroup_vars = c("z1", "z2"),
  subgroup_cuts = list(z1 = 1L, z2 = 1L),
  k_inter_range = c(0, 5), grid_step = 0.1, n_super = 3000L, seed = 42)
check("Cal binary: |OR(Q)-2|<0.3",
      abs(dgm_cb$hazard_ratios$harm_subgroup - 2.0) < 0.3,
      sprintf("%.3f", dgm_cb$hazard_ratios$harm_subgroup))

dgm_cc <- calibrate_glm_interaction(
  syn, c("z1", "z2", "z3"), "y_cont", "treat", 25,
  "continuous", "MD", subgroup_vars = c("z1", "z2"),
  subgroup_cuts = list(z1 = 1L, z2 = 1L),
  k_inter_range = c(0, 50), grid_step = 1, n_super = 3000L, seed = 42)
check("Cal continuous: |MD(Q)-25|<3",
      abs(dgm_cc$hazard_ratios$harm_subgroup - 25) < 3,
      sprintf("%.1f", dgm_cc$hazard_ratios$harm_subgroup))

dgm_cct <- calibrate_glm_interaction(
  syn, c("z1", "z2", "z3"), "y_count", "treat", 2.0,
  "count", "IRR", subgroup_vars = c("z1", "z2"),
  subgroup_cuts = list(z1 = 1L, z2 = 1L),
  k_inter_range = c(0, 3), grid_step = 0.1, n_super = 3000L, seed = 42)
check("Cal count: |IRR(Q)-2|<0.3",
      abs(dgm_cct$hazard_ratios$harm_subgroup - 2.0) < 0.3,
      sprintf("%.3f", dgm_cct$hazard_ratios$harm_subgroup))

dgm_coff <- calibrate_glm_interaction(
  syn, c("z1", "z2", "z3"), "y_events", "treat", 2.0,
  "count", "IRR", offset_var = "follow_up",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  k_inter_range = c(0, 3), grid_step = 0.1, n_super = 3000L, seed = 42)
check("Cal count+offset: |IRR(Q)-2|<0.5",
      abs(dgm_coff$hazard_ratios$harm_subgroup - 2.0) < 0.5,
      sprintf("%.3f", dgm_coff$hazard_ratios$harm_subgroup))

# Monotonicity
or_mono <- vapply(seq(0, 2, by = 0.4), function(k)
  generate_glm_dgm(syn, c("z1", "z2", "z3"), "y_bin", "treat", "binary",
    subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
    model = "alt", k_inter = k, n_super = 2000L, seed = 42
  )$hazard_ratios$harm_subgroup, numeric(1))
check("Mono binary: OR(Q) ↑ with k_inter", all(diff(or_mono) > 0),
      paste(round(or_mono, 3), collapse = " → "))

md_mono <- vapply(seq(0, 20, by = 5), function(k)
  generate_glm_dgm(syn, c("z1", "z2", "z3"), "y_cont", "treat",
    "continuous", "MD",
    subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
    model = "alt", k_inter = k, n_super = 2000L, seed = 42
  )$hazard_ratios$harm_subgroup, numeric(1))
check("Mono continuous: MD(Q) ↑ with k_inter", all(diff(md_mono) > 0),
      paste(round(md_mono, 1), collapse = " → "))

irr_mono <- vapply(seq(0, 1, by = 0.25), function(k)
  generate_glm_dgm(syn, c("z1", "z2", "z3"), "y_count", "treat",
    "count", "IRR",
    subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
    model = "alt", k_inter = k, n_super = 2000L, seed = 42
  )$hazard_ratios$harm_subgroup, numeric(1))
check("Mono count: IRR(Q) ↑ with k_inter", all(diff(irr_mono) > 0),
      paste(round(irr_mono, 3), collapse = " → "))

cat("\n")


# ═══════════════════════════════════════════════════════════════════════════
cat("─── 4. Simulation Fidelity ────────────────────────────────\n")
# ═══════════════════════════════════════════════════════════════════════════

# Binary — simulated OR matches DGM over replicates
or_reps <- replicate(50, {
  df <- simulate_from_glm_dgm(dgm_ba, n = 500, seed = NULL)
  fit <- glm(y_sim ~ treat_sim, data = df, family = binomial)
  exp(coef(fit)["treat_sim"])
})
check("Bin sim: mean(OR) within 30% of DGM",
      abs(mean(or_reps) - dgm_ba$hazard_ratios$overall) /
        dgm_ba$hazard_ratios$overall < 0.30,
      sprintf("sim=%.3f dgm=%.3f", mean(or_reps), dgm_ba$hazard_ratios$overall))

# Continuous — simulated MD matches DGM
md_reps <- replicate(50, {
  df <- simulate_from_glm_dgm(dgm_ca, n = 500, seed = NULL)
  mean(df$y_sim[df$treat_sim == 1]) - mean(df$y_sim[df$treat_sim == 0])
})
check("Cont sim: |mean(MD)−DGM|<5",
      abs(mean(md_reps) - dgm_ca$hazard_ratios$overall) < 5,
      sprintf("sim=%.1f dgm=%.1f", mean(md_reps), dgm_ca$hazard_ratios$overall))

# Binary y_sim ∈ {0,1}
df_bs <- simulate_from_glm_dgm(dgm_ba, n = 500, seed = 1)
check("Bin sim: y ∈ {0,1}", all(df_bs$y_sim %in% c(0L, 1L)))

# Count y_sim ≥ 0 and integer
df_cts <- simulate_from_glm_dgm(dgm_ct, n = 500, seed = 1)
check("Count sim: y≥0", all(df_cts$y_sim >= 0))
check("Count sim: integers", all(df_cts$y_sim == floor(df_cts$y_sim)))

# Continuous SD ≈ sigma
df_cs <- simulate_from_glm_dgm(dgm_ca, n = 2000, seed = 1)
check("Cont sim: SD within 30% of sigma",
      abs(sd(df_cs$y_sim) - dgm_ca$model_params$sigma) /
        dgm_ca$model_params$sigma < 0.30,
      sprintf("sd=%.1f sig=%.1f", sd(df_cs$y_sim), dgm_ca$model_params$sigma))

# flag_harm preserved
check("Sim: flag_harm exists", "flag_harm" %in% names(df_bs))
check("Sim: flag_harm has {0,1}", all(c(0, 1) %in% df_bs$flag_harm))

# n=NULL → full superpop
df_full <- simulate_from_glm_dgm(dgm_ba, n = NULL, seed = 1)
check("Sim n=NULL: nrow=superpop", nrow(df_full) == nrow(dgm_ba$df_super))

# Offset preserved
df_os <- simulate_from_glm_dgm(dgm_off, n = 500, seed = 1)
check("Sim: follow_up preserved", "follow_up" %in% names(df_os))

# Subgroup-specific effect in simulated data
md_Q <- replicate(30, {
  df <- simulate_from_glm_dgm(dgm_ca, n = 800, seed = NULL)
  mean(df$y_sim[df$flag_harm == 1 & df$treat_sim == 1]) -
    mean(df$y_sim[df$flag_harm == 1 & df$treat_sim == 0])
})
md_Qc <- replicate(30, {
  df <- simulate_from_glm_dgm(dgm_ca, n = 800, seed = NULL)
  mean(df$y_sim[df$flag_harm == 0 & df$treat_sim == 1]) -
    mean(df$y_sim[df$flag_harm == 0 & df$treat_sim == 0])
})
check("Cont sim: MD(Q)>MD(Qc) on avg",
      mean(md_Q) > mean(md_Qc),
      sprintf("Q=%.1f Qc=%.1f", mean(md_Q), mean(md_Qc)))
check("Cont sim: MD(Q) near DGM",
      abs(mean(md_Q) - dgm_ca$hazard_ratios$harm_subgroup) < 5,
      sprintf("sim=%.1f dgm=%.1f", mean(md_Q), dgm_ca$hazard_ratios$harm_subgroup))

cat("\n")


# ═══════════════════════════════════════════════════════════════════════════
cat("─── 5. Edge Cases: Subgroup Sizes ─────────────────────────\n")
# ═══════════════════════════════════════════════════════════════════════════

syn$z_rare <- as.factor(ifelse(runif(N) < 0.08, 1L, 0L))
dgm_tiny <- generate_glm_dgm(
  syn, c("z1", "z2", "z_rare"), "y_bin", "treat", "binary",
  subgroup_vars = "z_rare", subgroup_cuts = list(z_rare = 1L),
  model = "alt", k_inter = 2.0, n_super = 5000L, seed = 42)
check("Tiny SG (~8%): valid", inherits(dgm_tiny, "glm_dgm"))
check("Tiny SG: OR(Q)>OR(Qc)",
      dgm_tiny$hazard_ratios$harm_subgroup > dgm_tiny$hazard_ratios$no_harm_subgroup)

syn$z_big <- as.factor(ifelse(runif(N) < 0.75, 1L, 0L))
dgm_big <- generate_glm_dgm(
  syn, c("z1", "z2", "z_big"), "y_bin", "treat", "binary",
  subgroup_vars = "z_big", subgroup_cuts = list(z_big = 1L),
  model = "alt", k_inter = 1.0, n_super = 5000L, seed = 42)
check("Large SG (~75%): prev>65%", dgm_big$subgroup_info$proportion > 0.65)

dgm_3v <- generate_glm_dgm(
  syn, c("z1", "z2", "z3", "z4"), "y_bin", "treat", "binary",
  subgroup_vars = c("z1", "z2", "z3"),
  subgroup_cuts = list(z1 = 1L, z2 = 1L, z3 = 1L),
  model = "alt", k_inter = 1.5, n_super = 5000L, seed = 42)
check("3-var SG: prev≈12.5%",
      abs(dgm_3v$subgroup_info$proportion - 0.125) < 0.05,
      sprintf("%.3f", dgm_3v$subgroup_info$proportion))

dgm_1v <- generate_glm_dgm(
  syn, c("z1", "z2", "z3"), "y_bin", "treat", "binary",
  subgroup_vars = "z1", subgroup_cuts = list(z1 = 1L),
  model = "alt", k_inter = 1.0, n_super = 3000L, seed = 42)
check("1-var SG: prev≈50%",
      abs(dgm_1v$subgroup_info$proportion - 0.50) < 0.10)

cat("\n")


# ═══════════════════════════════════════════════════════════════════════════
cat("─── 6. k_inter Extremes ───────────────────────────────────\n")
# ═══════════════════════════════════════════════════════════════════════════

dgm_k0 <- generate_glm_dgm(
  syn, c("z1", "z2", "z3"), "y_bin", "treat", "binary",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "alt", k_inter = 0, n_super = 3000L, seed = 42)
check("k=0: OR(Q)≈OR(Qc)",
      abs(dgm_k0$hazard_ratios$harm_subgroup -
          dgm_k0$hazard_ratios$no_harm_subgroup) < 0.05)

dgm_huge <- tryCatch(
  generate_glm_dgm(syn, c("z1", "z2", "z3"), "y_bin", "treat", "binary",
    subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
    model = "alt", k_inter = 10, n_super = 2000L, seed = 42),
  error = function(e) e)
check("k=10 binary: no crash", inherits(dgm_huge, "glm_dgm"))
if (inherits(dgm_huge, "glm_dgm")) {
  check("k=10: no NaN", !any(is.nan(unlist(dgm_huge$hazard_ratios))))
  check("k=10: p1∈[0,1]",
        all(dgm_huge$df_super$p1 >= 0 & dgm_huge$df_super$p1 <= 1))
  check("k=10: OR(Q)>10",
        dgm_huge$hazard_ratios$harm_subgroup > 10,
        sprintf("%.1f", dgm_huge$hazard_ratios$harm_subgroup))
}

dgm_neg_big <- generate_glm_dgm(
  syn, c("z1", "z2", "z3"), "y_bin", "treat", "binary",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "alt", k_inter = -10, n_super = 2000L, seed = 42)
check("k=-10: OR(Q)<0.01",
      dgm_neg_big$hazard_ratios$harm_subgroup < 0.01,
      sprintf("%.6f", dgm_neg_big$hazard_ratios$harm_subgroup))
check("k=-10: p1∈[0,1]",
      all(dgm_neg_big$df_super$p1 >= 0 & dgm_neg_big$df_super$p1 <= 1))

dgm_c_huge <- generate_glm_dgm(
  syn, c("z1", "z2", "z3"), "y_cont", "treat", "continuous", "MD",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "alt", k_inter = 500, n_super = 2000L, seed = 42)
check("Cont k=500: MD(Q)≈505",
      abs(dgm_c_huge$hazard_ratios$harm_subgroup - 505) < 10,
      sprintf("%.1f", dgm_c_huge$hazard_ratios$harm_subgroup))

dgm_ct_huge <- generate_glm_dgm(
  syn, c("z1", "z2", "z3"), "y_count", "treat", "count", "IRR",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "alt", k_inter = 5, n_super = 2000L, seed = 42)
check("Count k=5: IRR(Q)>100",
      dgm_ct_huge$hazard_ratios$harm_subgroup > 100,
      sprintf("%.1f", dgm_ct_huge$hazard_ratios$harm_subgroup))
check("Count k=5: mu1 finite", all(is.finite(dgm_ct_huge$df_super$mu1)))

cat("\n")


# ═══════════════════════════════════════════════════════════════════════════
cat("─── 7. Offset Edge Cases ──────────────────────────────────\n")
# ═══════════════════════════════════════════════════════════════════════════

err_z <- tryCatch(
  generate_glm_dgm(transform(syn, bad = 0), c("z1"), "y_events", "treat",
    "count", "IRR", offset_var = "bad",
    subgroup_vars = "z1", subgroup_cuts = list(z1 = 1L),
    model = "null", n_super = 100L),
  error = function(e) e)
check("Offset zero → error", inherits(err_z, "error"))

err_neg <- tryCatch(
  generate_glm_dgm(transform(syn, bad = -1), c("z1"), "y_events", "treat",
    "count", "IRR", offset_var = "bad",
    subgroup_vars = "z1", subgroup_cuts = list(z1 = 1L),
    model = "null", n_super = 100L),
  error = function(e) e)
check("Offset negative → error", inherits(err_neg, "error"))

syn$fu_one <- rep(1, N)
dgm_off1 <- generate_glm_dgm(
  syn, c("z1", "z2"), "y_count", "treat", "count", "IRR",
  offset_var = "fu_one",
  subgroup_vars = "z1", subgroup_cuts = list(z1 = 1L),
  model = "alt", k_inter = 0.5, n_super = 3000L, seed = 42)
dgm_nooff <- generate_glm_dgm(
  syn, c("z1", "z2"), "y_count", "treat", "count", "IRR",
  subgroup_vars = "z1", subgroup_cuts = list(z1 = 1L),
  model = "alt", k_inter = 0.5, n_super = 3000L, seed = 42)
check("Offset=1 ≈ no offset",
      abs(dgm_off1$hazard_ratios$harm_subgroup -
          dgm_nooff$hazard_ratios$harm_subgroup) < 0.3,
      sprintf("%.3f vs %.3f", dgm_off1$hazard_ratios$harm_subgroup,
              dgm_nooff$hazard_ratios$harm_subgroup))

# Variable offset: mu scales with FU
dgm_voff <- generate_glm_dgm(
  syn, c("z1", "z2"), "y_events", "treat", "count", "IRR",
  offset_var = "follow_up",
  subgroup_vars = "z1", subgroup_cuts = list(z1 = 1L),
  model = "alt", k_inter = 0.5, n_super = 3000L, seed = 42)
check("Var offset: mu0~FU corr>0.5",
      cor(dgm_voff$df_super$mu0, dgm_voff$df_super$follow_up) > 0.5)

# Tiny offset (near zero, positive)
syn$fu_tiny <- rep(0.001, N)
syn$y_tiny <- rpois(N, exp(-0.5) * 0.001)
dgm_to <- tryCatch(
  generate_glm_dgm(syn, c("z1"), "y_tiny", "treat", "count", "IRR",
    offset_var = "fu_tiny",
    subgroup_vars = "z1", subgroup_cuts = list(z1 = 1L),
    model = "alt", k_inter = 0.5, n_super = 1000L, seed = 42),
  error = function(e) e)
check("Offset tiny(0.001): no crash", inherits(dgm_to, "glm_dgm"))

cat("\n")


# ═══════════════════════════════════════════════════════════════════════════
cat("─── 8. Input Validation ───────────────────────────────────\n")
# ═══════════════════════════════════════════════════════════════════════════

check("Non-binary for binary → error", inherits(tryCatch(
  generate_glm_dgm(transform(syn, yb = rnorm(N)), c("z1"), "yb", "treat",
    "binary", subgroup_vars = "z1", subgroup_cuts = list(z1 = 1L),
    model = "null", n_super = 100L), error = function(e) e), "error"))

check("Neg counts → error", inherits(tryCatch(
  generate_glm_dgm(transform(syn, yn = -rpois(N, 3)), c("z1"), "yn", "treat",
    "count", subgroup_vars = "z1", subgroup_cuts = list(z1 = 1L),
    model = "null", n_super = 100L), error = function(e) e), "error"))

check("Non-numeric continuous → error", inherits(tryCatch(
  generate_glm_dgm(transform(syn, yc = as.character(y_bin)), c("z1"), "yc",
    "treat", "continuous", subgroup_vars = "z1",
    subgroup_cuts = list(z1 = 1L), model = "null", n_super = 100L),
  error = function(e) e), "error"))

check("Non-glm_dgm to simulate → error", inherits(tryCatch(
  simulate_from_glm_dgm(list(a = 1), n = 100),
  error = function(e) e), "error"))

check("Missing SG vars in calibrate → error", inherits(tryCatch(
  calibrate_glm_interaction(syn, c("z1"), "y_bin", "treat", 2.0,
    "binary", subgroup_vars = NULL, subgroup_cuts = NULL,
    k_inter_range = c(0, 3), grid_step = 1),
  error = function(e) e), "error"))

check("Non-integer counts → error", inherits(tryCatch(
  generate_glm_dgm(transform(syn, yf = runif(N, 0, 5)), c("z1"), "yf",
    "treat", "count", subgroup_vars = "z1",
    subgroup_cuts = list(z1 = 1L), model = "null", n_super = 100L),
  error = function(e) e), "error"))

cat("\n")


# ═══════════════════════════════════════════════════════════════════════════
cat("─── 9. Reproducibility ────────────────────────────────────\n")
# ═══════════════════════════════════════════════════════════════════════════

dgm_r1 <- generate_glm_dgm(syn, c("z1", "z2"), "y_bin", "treat", "binary",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "alt", k_inter = 1, n_super = 2000L, seed = 12345)
dgm_r2 <- generate_glm_dgm(syn, c("z1", "z2"), "y_bin", "treat", "binary",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "alt", k_inter = 1, n_super = 2000L, seed = 12345)
check("Same seed → identical DGM",
      identical(dgm_r1$hazard_ratios, dgm_r2$hazard_ratios))

df_s1 <- simulate_from_glm_dgm(dgm_r1, n = 200, seed = 99)
df_s2 <- simulate_from_glm_dgm(dgm_r1, n = 200, seed = 99)
check("Same seed → identical sim", identical(df_s1$y_sim, df_s2$y_sim))

df_s3 <- simulate_from_glm_dgm(dgm_r1, n = 200, seed = 100)
check("Diff seed → diff sim", !identical(df_s1$y_sim, df_s3$y_sim))

cat("\n")


# ═══════════════════════════════════════════════════════════════════════════
cat("─── 10. print.glm_dgm ────────────────────────────────────\n")
# ═══════════════════════════════════════════════════════════════════════════

check("print binary: has OR",
      any(grepl("OR", capture.output(print(dgm_ba)))))
check("print continuous: has MD",
      any(grepl("MD", capture.output(print(dgm_ca)))))
check("print count: has IRR",
      any(grepl("IRR", capture.output(print(dgm_ct)))))

cat("\n")


# ═══════════════════════════════════════════════════════════════════════════
cat("─── 11. Cross-Outcome Type Isolation ──────────────────────\n")
# ═══════════════════════════════════════════════════════════════════════════

check("Binary: no mu0/mu1", !("mu0" %in% names(dgm_ba$df_super)))
check("Continuous: no p0/p1", !("p0" %in% names(dgm_ca$df_super)))
check("Count: no p0/p1", !("p0" %in% names(dgm_ct$df_super)))
check("Count: has mu0/mu1",
      all(c("mu0", "mu1") %in% names(dgm_ct$df_super)))
check("Binary: has p0/p1",
      all(c("p0", "p1") %in% names(dgm_ba$df_super)))

cat("\n")


# ═══════════════════════════════════════════════════════════════════════════
cat("─── 12. Null Model Across All Types ──────────────────────\n")
# ═══════════════════════════════════════════════════════════════════════════

for (info in list(
  list(ot = "binary",     ov = "y_bin",   em = "OR",  tol = 0.15),
  list(ot = "continuous", ov = "y_cont",  em = "MD",  tol = 2),
  list(ot = "count",      ov = "y_count", em = "IRR", tol = 0.15)
)) {
  dgm_n <- generate_glm_dgm(syn, c("z1", "z2", "z3"), info$ov, "treat",
    info$ot, info$em,
    subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
    model = "null", k_inter = 999, n_super = 3000L, seed = 42)
  d <- abs(dgm_n$hazard_ratios$harm_subgroup - dgm_n$hazard_ratios$no_harm_subgroup)
  check(sprintf("Null %s: Q≈Qc", info$ot), d < info$tol,
        sprintf("|%.3f−%.3f|=%.3f", dgm_n$hazard_ratios$harm_subgroup,
                dgm_n$hazard_ratios$no_harm_subgroup, d))
}

cat("\n")


# ═══════════════════════════════════════════════════════════════════════════
cat("=============================================================\n")
cat(sprintf("  RESULTS: %d passed, %d failed, %d warnings\n",
            pass, fail, warn_count))
cat("=============================================================\n")
if (fail > 0) {
  cat("\n*** FAILURES — review above ***\n")
  quit(status = 1)
} else {
  cat("\n=== All checks passed ===\n")
}
