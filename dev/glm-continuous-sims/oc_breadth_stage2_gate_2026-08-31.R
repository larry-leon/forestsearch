# =============================================================================
# oc_breadth_stage2_gate_2026-08-31.R -- Stage 2 §2(a) forecast read-out, §3 DGM
# and family gates, pairing check, and the §2(c) single-replicate equivalence
# check (winner's fitted effect == recorder's nv_H_est; threshold acts as a
# post-hoc filter on it).  Task: dev/tasks/cc_task_oc_breadth_stage2_2026-08-31.md
# Reads only.  Writes dev/glm-continuous-sims/oc_breadth_stage2_2026-08-31_gate.rds
# =============================================================================
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
stopifnot(requireNamespace("speff2trial", quietly = TRUE))
DIR  <- file.path("dev", "glm-continuous-sims")
say <- function(...) cat(sprintf(...), "\n", sep = "")
tic <- function() proc.time()[["elapsed"]]
f6 <- function(x) sprintf("%.6f", x)

# ---- §2(a): the locked forecast, read from c9cb0ca2's artefacts --------------
FC <- readRDS(file.path(DIR, "oc_breadth_ladder_2026-08-30_forecast120.rds"))
G0 <- readRDS(file.path(DIR, "oc_breadth_ladder_2026-08-30_gate.rds"))
CORR <- readRDS(file.path(DIR, "oc_wrapper_grid_corrected_2026-08-30.rds"))
stopifnot(FC$q == 120, all(FC$gates), all(FC$null_gates))
k_inter <- FC$k_inter; c1_star <- FC$c1_star; c1_05 <- FC$c1_05
say("[2a] forecast rung q = %g ; k_inter = %s (= k40 %s + s(%d)*(120-40)) ; c1* = %s ; c1_05 = %s ; c1_10 = %s",
    FC$q, f6(k_inter), f6(G0$k40), G0$s, f6(c1_star), f6(c1_05), f6(FC$c1_10))
stopifnot(abs(k_inter - (G0$k40 + G0$s * (120 - 40))) < 1e-12)
qs <- c("det_rate", "det_rate_se", "EnH", "Esens", "Espec", "Eppv", "Enpv", "EbetaH", "Enaive_bias", "mass_below")
fc_tab <- data.frame(quantity = qs,
                     resample = sapply(qs, function(q) FC$star$resample[[q]]),
                     split    = sapply(qs, function(q) FC$star$split[[q]]),
                     null_resample = sapply(qs, function(q) FC$null_at_star[[q]]))
say("[2a] forecast table at (n = 500, c1* = %s, c2 = 10, pconsistency 0.9), both gates, and the null at c1*:", f6(c1_star))
print(format(fc_tab, digits = 7), row.names = FALSE)
say("[2a] sel_c mass on PQg >= 0.95: resample %s split %s ; mass on Q: %s / %s",
    f6(FC$sel_star$resample$mass_PQg95), f6(FC$sel_star$split$mass_PQg95), f6(FC$sel_star$resample$mass_Q), f6(FC$sel_star$split$mass_Q))
say("[2a] at the driver's c1 = 30 (resample): det %s EnH %s Esens %s Espec %s Eppv %s EbetaH %s naive %s",
    f6(FC$named$resample$driver$det_rate), f6(FC$named$resample$driver$EnH), f6(FC$named$resample$driver$Esens),
    f6(FC$named$resample$driver$Espec), f6(FC$named$resample$driver$Eppv), f6(FC$named$resample$driver$EbetaH), f6(FC$named$resample$driver$Enaive_bias))
grid_alt <- FC$grid$resample$table
say("[2a] forecast rung resample grid: %d rows, c1 %d..%d ; det_rate at 30 %s, 133 %s, 136 %s", nrow(grid_alt), min(grid_alt$c1), max(grid_alt$c1),
    f6(grid_alt$det_rate[grid_alt$c1 == 30]), f6(grid_alt$det_rate[grid_alt$c1 == 133]), f6(grid_alt$det_rate[grid_alt$c1 == 136]))
null_grid <- CORR$null$sweep$table[CORR$null$sweep$table$consistency_method == "resample", ]
null_inv  <- CORR$null$invert$table[CORR$null$invert$table$consistency_method == "resample", ]
say("[2a] null curve artefacts: corrected null sweep (resample) c1 %s ; null inversions (resample) at targets %s -> c1 %s",
    paste(range(null_grid$c1), collapse = ".."), paste(null_inv$target, collapse = "/"), paste(f6(null_inv$value), collapse = "/"))

# ---- the driver frame and the direct build (Stage 1 script, verbatim) --------
driver_frame <- function() {
  actg_arms <- c(1L, 3L); actg_treat_arm <- 1L; actg_age_cut <- 34; actg_preanti_cut <- 744.5
  analysis_binary_vars <- c("hemo","homo","drugs","race","gender","symptom", "str2")
  actg_df <- subset(speff2trial::ACTG175, arms %in% actg_arms)
  actg_df$id <- seq_len(nrow(actg_df))
  actg_df$treat_orig <- ifelse(actg_df$arms == actg_treat_arm, 1L, 0L)
  actg_df$treat      <- 1L - actg_df$treat_orig
  actg_df$cd4_change <- actg_df$cd420 - actg_df$cd40
  actg_df <- actg_df[!is.na(actg_df$cd420), ]
  actg_df$z1  <- as.factor(ifelse(actg_df$age > actg_age_cut, 1L, 0L))
  actg_df$z2  <- as.factor(ifelse(actg_df$preanti <= actg_preanti_cut, 1L, 0L))
  actg_df$z3  <- as.factor(ifelse(actg_df$wtkg <= 75, 1L, 0L))
  actg_df$z4  <- as.factor(ifelse(actg_df$karnof <= median(actg_df$karnof), 1L, 0L))
  actg_df$z5  <- as.factor(ifelse(actg_df$cd40 <= median(actg_df$cd40), 1L, 0L))
  actg_df$z6  <- as.factor(ifelse(actg_df$cd80 <= median(actg_df$cd80), 1L, 0L))
  actg_df$z7  <- as.factor(actg_df$hemo);  actg_df$z8  <- as.factor(actg_df$homo)
  actg_df$z9  <- as.factor(actg_df$drugs); actg_df$z10 <- as.factor(actg_df$race)
  actg_df$z11 <- as.factor(actg_df$gender);actg_df$z12 <- as.factor(actg_df$symptom)
  for (v in analysis_binary_vars) actg_df[[v]] <- as.factor(actg_df[[v]])
  actg_df
}
build_direct <- function(data, k_inter) generate_glm_dgm(
  data = data, factor_vars = paste0("z", 1:12),
  outcome_var = "cd4_change", treatment_var = "treat",
  outcome_type = "continuous", effect_measure = "MD",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "alt", k_treat = 1, k_inter = k_inter, adverse_outcome = FALSE,
  n_super = 5000L, seed = 8316951L, verbose = FALSE)
fs_args <- CORR$fs_args
stopifnot(fs_args$effect.threshold == 30, fs_args$consistency.threshold == 10, fs_args$pconsistency.threshold == 0.9)

# ---- §3 DGM gates ---------------------------------------------------------------
t0 <- tic()
df_drv <- driver_frame()
dgm120 <- build_direct(df_drv, k_inter)
dgm40  <- build_direct(df_drv, G0$k40)
sc <- fs_dgm_scale(dgm120)$regions; mt <- setNames(sc$m_tau, sc$region)
g_Q  <- abs(abs(mt[["Q"]]) - 120) < 1e-9
g_Qc <- abs(abs(mt[["Qc"]]) - abs(G0$m_tau[["Qc"]])) < 1e-9
say("[3] fs_dgm_scale: |m_tau[Q]| = %.12f (gate 120: %s) ; |m_tau[Qc]| = %.12f vs MD40 %.12f (gate: %s) ; effect_Q = %.10f",
    abs(mt[["Q"]]), g_Q, abs(mt[["Qc"]]), abs(G0$m_tau[["Qc"]]), g_Qc, dgm120$hazard_ratios$harm_subgroup)
stopifnot(g_Q, g_Qc)
t1 <- tic()
fam <- fs_oc_family_enumerate(dgm120, fs_args, n = 500, max_M = 5000L, verbose = FALSE)
t_fam <- tic() - t1
F0 <- G0$family
cmp <- sapply(c("lab", "Pg", "PQg", "se_g", "sens_g", "spec_g", "M"), function(f) identical(fam[[f]], F0[[f]]))
cmp <- c(cmp, beta_g_forecast = identical(fam$beta_g, FC$beta_g))
say("[3] family M = %d (%.0f s); identical() to the stored forecast-rung family (gate.rds fields + forecast120 beta_g): %s",
    fam$M, t_fam, paste(names(cmp), cmp, collapse = " "))
stopifnot(all(cmp))
# df_super: identical to MD40's except mu1 (the shift of Q's treated outcomes)
ds40 <- dgm40$df_super; ds120 <- dgm120$df_super
same_cols <- sapply(names(ds40), function(v) identical(ds40[[v]], ds120[[v]]))
say("[3] df_super columns differing between MD40 and MD120: %s", paste(names(same_cols)[!same_cols], collapse = ", "))
inQ <- ds120$flag_harm == 1L
dmu1 <- ds120$mu1 - ds40$mu1
say("[3] mu1 shift: on Q range [%.10f, %.10f] ; off Q max|.| = %.3e ; mu0 identical %s",
    min(dmu1[inQ]), max(dmu1[inQ]), max(abs(dmu1[!inQ])), identical(ds40$mu0, ds120$mu0))

# ---- pairing: the same seed draws the same rows, treatment and residuals -------
seed_base <- 8316951L
pair <- t(sapply(1:5, function(i) {
  sd_i <- seed_base + i
  RNGkind("L'Ecuyer-CMRG"); set.seed(sd_i); a <- simulate_from_glm_dgm(dgm40,  n = 500, seed = sd_i)
  RNGkind("L'Ecuyer-CMRG"); set.seed(sd_i); b <- simulate_from_glm_dgm(dgm120, n = 500, seed = sd_i)
  covs <- setdiff(names(a), c("y_sim", "mu1"))
  # residual stream: y differs from MD40's y only by the mu1 shift (exactly -80 on Q-treated, 0 elsewhere)
  # -- the RNG stream (rows, treatment, rnorm) is consumed identically
  c(sim_id = i, rows_treat_identical = identical(a[covs], b[covs]),
    max_abs_dy_offQ_or_control = max(abs((b$y_sim - a$y_sim)[!(a$flag_harm == 1L & a$treat_sim == 1L)])),
    dy_Q_treated_min = min((b$y_sim - a$y_sim)[a$flag_harm == 1L & a$treat_sim == 1L]),
    dy_Q_treated_max = max((b$y_sim - a$y_sim)[a$flag_harm == 1L & a$treat_sim == 1L]))
}))
say("[3] pairing check on replicates 1..5 (seed_base + sim_id, RNGkind L'Ecuyer-CMRG as the driver):"); print(pair, digits = 12)

# ---- §2(c) single-replicate check: winner's fitted effect and the threshold ----
# the driver's forestsearch() call, verbatim (sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_1000.qmd L560-590)
confs <- c("age","preanti","wtkg","karnof","cd40","cd80","hemo","homo","drugs","race","gender","symptom","str2")
run_fs <- function(df, c1, sd_i) {
  t <- tic()
  msgs <- utils::capture.output(fs <- suppressWarnings(forestsearch(
    df.analysis = df, confounders.name = confs,
    outcome.name = "y_sim", treat.name = "treat_sim", id.name = "id",
    outcome_type = "continuous", effect_measure = "MD",
    effect.threshold = c1, consistency.threshold = 10,
    pconsistency.threshold = 0.90, fs.splits = 400L,
    n.min = 60L, d0.min = 12L, d1.min = 12L, maxk = 2L,
    vi.grf.min = -0.2, sg_focus = "maxeffCons",
    selection_rule = "neighborhood", effect_neighborhood = 0.10, stop_threshold = NULL,
    consistency_method = "resample", conf.cont_jcuts = list(age = 10, preanti = 10),
    use_lasso = FALSE, use_dina = FALSE, use_grf = FALSE, use_twostage = TRUE, is.RCT = TRUE,
    adverse_outcome = FALSE, details = FALSE, quiet = FALSE, seedit = sd_i,
    parallel_args = list(plan = "sequential"), mr_inference = TRUE,
    mr_inference_args = list(ci_method = "ij", draws = 5000L, include_complement = TRUE))), type = "message")
  list(fs = fs, secs = tic() - t)
}
one <- list()
for (i in 1:2) {
  sd_i <- seed_base + i
  RNGkind("L'Ecuyer-CMRG"); set.seed(sd_i); df <- simulate_from_glm_dgm(dgm120, n = 500, seed = sd_i)
  df$id <- seq_len(nrow(df))
  r30 <- run_fs(df, 30, sd_i)
  fs <- r30$fs
  win <- if (!is.null(fs$sg.harm)) paste(fs$sg.harm, collapse = " & ") else NA_character_
  # the winner's fitted effect as the search/consistency table carries it, and as the recorder stores it (mr naive est)
  osg <- fs$grp.consistency$out_sg
  nv <- fs$mr_inference$naive$est
  say("[2c] rep %d at c1 = 30: %.2f s ; winner = %s ; n_harm = %d ; mr naive est = %.10f", i, r30$secs, win,
      sum(fs$grp.consistency$sg.harm.id == 1L), nv)
  r1 <- as.data.frame(osg$result)[1, c("Pcons", "hr", "N", "E", "K", "M.1", "M.2")]
  cat("      consistency table row 1 (the winner):\n"); print(r1, digits = 12)
  say("      recorder's nv_H_est minus the table's hr = %.3e ; candidates entering consistency: %d", nv - r1$hr, fs$grp.consistency$n_candidates_total)
  # the same replicate at c1 = c1*: declares iff nv >= c1*, same winner
  rst <- run_fs(df, c1_star, sd_i)
  fs2 <- rst$fs
  win2 <- if (!is.null(fs2$sg.harm)) paste(fs2$sg.harm, collapse = " & ") else NA_character_
  nv2 <- if (!is.null(fs2$mr_inference)) fs2$mr_inference$naive$est else NA_real_
  say("[2c] rep %d at c1 = c1* = %s: %.2f s ; declared %s (post-hoc rule says %s) ; winner = %s ; naive est = %s ; membership identical %s",
      i, f6(c1_star), rst$secs, !is.null(fs2$sg.harm), nv >= c1_star, win2, if (is.na(nv2)) "NA" else f6(nv2),
      if (!is.null(fs2$sg.harm)) identical(fs$grp.consistency$sg.harm.id, fs2$grp.consistency$sg.harm.id) else NA)
  one[[i]] <- list(sim_id = i, secs30 = r30$secs, secs_star = rst$secs, winner30 = win, nv30 = nv, table_row1 = r1,
                   declared_star = !is.null(fs2$sg.harm), winner_star = win2, nv_star = nv2,
                   memb_identical = if (!is.null(fs2$sg.harm)) identical(fs$grp.consistency$sg.harm.id, fs2$grp.consistency$sg.harm.id) else NA,
                   n_candidates = list(c1_30 = fs$grp.consistency$n_candidates_total, c1_star = fs2$grp.consistency$n_candidates_total))
}
saveRDS(list(forecast = fc_tab, k_inter = k_inter, c1_star = c1_star, c1_05 = c1_05, c1_10 = FC$c1_10,
             mass_PQg95 = c(resample = FC$sel_star$resample$mass_PQg95, split = FC$sel_star$split$mass_PQg95),
             m_tau = mt, gates = c(Q = g_Q, Qc = g_Qc, cmp), fam_secs = t_fam, M = fam$M,
             df_super_diff_cols = names(same_cols)[!same_cols], pairing = pair, one_rep = one,
             forecast_commit = "c9cb0ca2", pkg_version = as.character(packageVersion("forestsearch")), built_at = Sys.time()),
        file.path(DIR, "oc_breadth_stage2_2026-08-31_gate.rds"))
say("[gate] written (total %.0f s)", tic() - t0)
