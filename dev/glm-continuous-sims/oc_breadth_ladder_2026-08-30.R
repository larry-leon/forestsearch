# =============================================================================
# oc_breadth_ladder_2026-08-30.R -- the interaction ladder on the MD40 DGM
# Task: dev/tasks/cc_task_oc_breadth_ladder_2026-08-30.md
#
# Parts (each a separate process; logs in *_<part>.log):
#   gate            §2  rebuild MD40 as the driver did; identity gates; one
#                       fs_oc_predict() identical to the corrected run
#   rung <q>        §3+§4  one ladder rung: DGM gates, family gates, the c1 grid
#                       0..300, the inversions, the se_g diagnostic
#   forecast <q>    §5+§6  the forecast rung under both gates, c1*, the null rate
#                       at c1*, the selection distribution, the se_direct
#                       sensitivity
# All settings are read from oc_wrapper_grid_corrected_2026-08-30.rds, unchanged.
# =============================================================================
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
stopifnot(requireNamespace("speff2trial", quietly = TRUE))
args <- commandArgs(trailingOnly = TRUE)
part <- args[1]; stopifnot(part %in% c("gate", "rung", "forecast"))
DIR  <- file.path("dev", "glm-continuous-sims")
STEM <- "oc_breadth_ladder_2026-08-30"
outp <- function(s) file.path(DIR, sprintf("%s_%s.rds", STEM, s))
CORR <- readRDS(file.path(DIR, "oc_wrapper_grid_corrected_2026-08-30.rds"))
fs_args <- CORR$fs_args
SEED <- CORR$seed; DRAWS <- CORR$draws
stopifnot(SEED == 20260825L, DRAWS == 2e5, CORR$alt$sweep$block == Inf)
C1_DRV <- fs_args$effect.threshold; C2_DRV <- fs_args$consistency.threshold
PCONS  <- fs_args$pconsistency.threshold
stopifnot(C1_DRV == 30, C2_DRV == 10, PCONS == 0.9)
N <- 500   # double, as the corrected script called run_alt(500)
say <- function(...) cat(sprintf(...), "\n", sep = "")
tic <- function() proc.time()[["elapsed"]]

# ---- the ACTG175 frame, verbatim from the driver (build-dgm chunk, L282-302) --
# sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_1000.qmd
driver_frame <- function() {
  actg_arms <- c(1L, 3L); actg_treat_arm <- 1L; actg_age_cut <- 34; actg_preanti_cut <- 744.5
  analysis_binary_vars <- c("hemo","homo","drugs","race","gender","symptom", "str2")  # include_str2 = TRUE
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
# the corrected script's frame (oc_wrapper_grid_corrected_2026-08-30.R, actg_frame)
corrected_frame <- function() {
  actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
  actg_df$id <- seq_len(nrow(actg_df))
  actg_df$treat <- 1L - ifelse(actg_df$arms == 1L, 1L, 0L)
  actg_df$cd4_change <- actg_df$cd420 - actg_df$cd40
  actg_df <- actg_df[!is.na(actg_df$cd420), ]
  actg_df$z1  <- as.factor(ifelse(actg_df$age > 34, 1L, 0L))
  actg_df$z2  <- as.factor(ifelse(actg_df$preanti <= 744.5, 1L, 0L))
  actg_df$z3  <- as.factor(ifelse(actg_df$wtkg <= 75, 1L, 0L))
  actg_df$z4  <- as.factor(ifelse(actg_df$karnof <= median(actg_df$karnof), 1L, 0L))
  actg_df$z5  <- as.factor(ifelse(actg_df$cd40 <= median(actg_df$cd40), 1L, 0L))
  actg_df$z6  <- as.factor(ifelse(actg_df$cd80 <= median(actg_df$cd80), 1L, 0L))
  actg_df$z7  <- as.factor(actg_df$hemo);  actg_df$z8  <- as.factor(actg_df$homo)
  actg_df$z9  <- as.factor(actg_df$drugs); actg_df$z10 <- as.factor(actg_df$race)
  actg_df$z11 <- as.factor(actg_df$gender);actg_df$z12 <- as.factor(actg_df$symptom)
  actg_df
}
# generate_glm_dgm() with every driver argument fixed; only k_inter (and data) vary
build_direct <- function(data, k_inter) generate_glm_dgm(
  data = data, factor_vars = paste0("z", 1:12),
  outcome_var = "cd4_change", treatment_var = "treat",
  outcome_type = "continuous", effect_measure = "MD",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "alt", k_treat = 1, k_inter = k_inter, adverse_outcome = FALSE,
  n_super = 5000L, seed = 8316951L, verbose = FALSE)

family_fields <- c("lab", "Pg", "PQg", "se_g", "sens_g", "spec_g", "M")
enumerate <- function(dgm) fs_oc_family_enumerate(dgm, fs_args, n = N, max_M = 5000L, verbose = TRUE)
mtau <- function(dgm) { r <- fs_dgm_scale(dgm)$regions; setNames(r$m_tau, r$region) }

# ---- the se_g diagnostic (§4): se_direct from fs_dgm_scale on each candidate ---
se_diag <- function(dgm, fam) {
  regs <- lapply(seq_len(fam$M), function(j) fam$memb[, j])
  names(regs) <- sprintf("g%04d", seq_len(fam$M))
  sc <- fs_dgm_scale(dgm, regions = regs)$regions
  stopifnot(isTRUE(all.equal(sc$P_g, fam$Pg)))
  se_direct <- sqrt(sc$V_eff / (N * fam$Pg))
  ratio <- fam$se_g / se_direct
  band <- cut(fam$PQg, c(-Inf, 0.25, 0.95, Inf), right = FALSE,
              labels = c("PQg<0.25", "0.25<=PQg<0.95", "PQg>=0.95"))
  qs <- function(v) c(min = min(v), q05 = unname(quantile(v, .05)), q25 = unname(quantile(v, .25)),
                      median = median(v), q75 = unname(quantile(v, .75)), q95 = unname(quantile(v, .95)),
                      max = max(v), n = length(v))
  list(se_direct = se_direct, V_eff = sc$V_eff, ratio = ratio,
       overall = qs(ratio),
       by_band = do.call(rbind, lapply(split(ratio, band), qs)),
       cor_Pg = cor(ratio, fam$Pg), cor_PQg = cor(ratio, fam$PQg),
       outside_2pct = sum(ratio < 0.98 | ratio > 1.02),
       outside_md40_band = sum(ratio < 0.992 | ratio > 1.015))
}
fam_record <- function(fam) fam[c("M", "counts", "lab", "Pg", "PQg", "beta_g", "se_g", "sens_g",
                                  "spec_g", "PQ", "cuts", "args_used", "null")]
grid_point <- function(sw, c1) {
  i <- which(sw$table$c1 == c1); stopifnot(length(i) == 1L)
  r <- sw$results[[i]]; r$family <- NULL; r
}
row_of <- function(p) data.frame(
  det_rate = p$det_rate, det_rate_se = p$det_rate_se, EnH = p$EnH, Esens = p$Esens,
  Espec = p$Espec, Eppv = p$Eppv, Enpv = p$Enpv, EbetaH = p$EbetaH,
  Enaive_bias = p$Enaive_bias, mass_below = p$mass_below)

# =============================================================================
# §2 GATE
# =============================================================================
run_gate <- function() {
  t0 <- tic()
  df_drv <- driver_frame()
  say("[gate] driver frame: %d rows, %d cols", nrow(df_drv), ncol(df_drv))
  # 2a. rebuild EXACTLY as the driver did: calibrate_glm_interaction()
  dgm_cal <- calibrate_glm_interaction(
    data = df_drv, factor_vars = paste0("z", 1:12),
    outcome_var = "cd4_change", treatment_var = "treat",
    target_effect = -40, outcome_type = "continuous",
    effect_measure = "MD", subgroup_vars = c("z1", "z2"),
    subgroup_cuts = list(z1 = 1L, z2 = 1L), k_inter_range = c(0, 120),
    grid_step = 2, n_super = 5000L, seed = 8316951L, verbose = FALSE)
  k40 <- dgm_cal$model_params$k_inter
  say("[gate] calibrated: k_inter = %.10f  beta_inter = %.10f  (%.0f s)",
      k40, dgm_cal$model_params$beta_inter, tic() - t0)
  mt <- mtau(dgm_cal); print(mt, digits = 12)
  # against the payload truth the corrected run gated on
  pl <- readRDS(CORR$alt$payloads[["500"]])
  say("[gate] payload truth beta_inter = %.10f ; |k40 - truth| = %.3e", pl$truth$beta_inter, abs(k40 - pl$truth$beta_inter))
  gate_k <- abs(k40 - pl$truth$beta_inter) < 1e-9
  say("[gate] k40 == payload truth$beta_inter (1e-9): %s", gate_k)
  # 2b. direct route on the SAME frame
  dgm_dir <- build_direct(df_drv, k40)
  g_df <- identical(dgm_cal$df_super, dgm_dir$df_super)
  say("[gate] identical(df_super calibrate, df_super direct): %s", g_df)
  stopifnot(gate_k, g_df)
  # 2c. family at n = 500 against the stored corrected family
  t1 <- tic(); fam <- enumerate(dgm_dir); t_fam <- tic() - t1
  S <- CORR$alt$families[["500"]]
  cmp <- c(lab = identical(fam$lab, S$lab), Pg = identical(fam$Pg, S$Pg), PQg = identical(fam$PQg, S$PQg),
           beta_g = identical(fam$beta_g, S$beta_g), se_g = identical(fam$se_g, S$se_g),
           sens_g = identical(fam$sens_g, S$sens_g), spec_g = identical(fam$spec_g, S$spec_g),
           M = identical(fam$M, S$M))
  say("[gate] family M = %d (stored %d), %.1f s", fam$M, S$M, t_fam); print(cmp)
  say("[gate] NOTE: the stored family record (family_record() in the corrected script) holds no ovl/memb; those cannot be gated against the store")
  if (!all(cmp)) {
    # diagnose: is it the frame?  (driver converts the analysis binaries to factor; corrected script did not)
    dgm_c <- build_direct(corrected_frame(), k40)
    fam_c <- enumerate(dgm_c)
    say("[gate] DIAGNOSIS on the corrected script's frame (integer binaries): lab %s Pg %s se_g %s M %d",
        identical(fam_c$lab, S$lab), identical(fam_c$Pg, S$Pg), identical(fam_c$se_g, S$se_g), fam_c$M)
    say("[gate] df_super identical across the two frames (columns in common)? %s",
        identical(dgm_c$df_super[names(dgm_c$df_super)], dgm_dir$df_super[names(dgm_c$df_super)]))
  }
  stopifnot(all(cmp), fam$M == 1696L)
  # 2d. draw-level: the exact call the corrected script made (fs_oc_predict, one block)
  t2 <- tic()
  pred <- fs_oc_predict(family = fam, n = N, c1 = C1_DRV, c2 = C2_DRV,
                        consistency_method = "resample", pconsistency = PCONS,
                        draws = DRAWS, seed = SEED)
  t_pred <- tic() - t2
  pred$family <- NULL
  S_pred <- CORR$alt$runs$n500_resample$pred
  g_pred <- identical(pred, S_pred)
  say("[gate] fs_oc_predict(n=500, c1=%g, c2=%g, resample, pcons=%g, draws=%g, seed=%d) identical() to stored: %s  (%.0f s; det_rate %.5f vs %.5f)",
      C1_DRV, C2_DRV, PCONS, DRAWS, SEED, g_pred, t_pred, pred$det_rate, S_pred$det_rate)
  if (!g_pred) { fld <- names(pred); print(sapply(fld, function(f) identical(pred[[f]], S_pred[[f]]))) }
  stopifnot(g_pred)
  # scale table for Q, Qc, plus the se_g diagnostic at MD40
  sc <- fs_dgm_scale(dgm_dir)
  diag <- se_diag(dgm_dir, fam)
  say("[gate] se_g diagnostic MD40: ratio range [%.5f, %.5f] median %.5f; outside [0.992,1.015]: %d", diag$overall["min"], diag$overall["max"], diag$overall["median"], diag$outside_md40_band)
  print(diag$by_band, digits = 5)
  saveRDS(list(k40 = k40, s = sign(mt[["Q"]]), m_tau = mt, scale = sc$regions,
               calibrate_call = "calibrate_glm_interaction(data = driver frame, factor_vars = z1..z12, outcome_var = cd4_change, treatment_var = treat, target_effect = -40, outcome_type = continuous, effect_measure = MD, subgroup_vars = c(z1,z2), subgroup_cuts = list(z1=1L,z2=1L), k_inter_range = c(0,120), grid_step = 2, n_super = 5000L, seed = 8316951L)",
               gates = c(k40_vs_truth = gate_k, df_super_identical = g_df, cmp, predict_identical = g_pred),
               family = fam_record(fam), pred = pred, pred_secs = t_pred, fam_secs = t_fam,
               se_diag = diag[c("overall", "by_band", "cor_Pg", "cor_PQg", "outside_2pct", "outside_md40_band")],
               se_direct = diag$se_direct, V_eff_g = diag$V_eff,
               pkg_version = as.character(packageVersion("forestsearch")), built_at = Sys.time()),
          outp("gate"))
  say("[gate] written %s  (total %.0f s)", outp("gate"), tic() - t0)
}

# =============================================================================
# §3 + §4  one rung
# =============================================================================
build_rung <- function(q, G) {
  k_q <- G$k40 + G$s * (q - 40)
  dgm <- build_direct(driver_frame(), k_q)
  mt <- mtau(dgm)
  gq  <- abs(abs(mt[["Q"]]) - q) < 1e-9
  gqc <- abs(abs(mt[["Qc"]]) - abs(G$m_tau[["Qc"]])) < 1e-9
  say("[rung %g] k_inter = %.10f ; |m_tau[Q]| = %.12f (gate %s) ; |m_tau[Qc]| = %.12f vs MD40 %.12f (gate %s)",
      q, k_q, abs(mt[["Q"]]), gq, abs(mt[["Qc"]]), abs(G$m_tau[["Qc"]]), gqc)
  stopifnot(gq, gqc)
  fam <- enumerate(dgm)
  F0 <- G$family
  cmp <- sapply(family_fields, function(f) identical(fam[[f]], F0[[f]]))
  db <- fam$beta_g - F0$beta_g - (q - 40) * fam$PQg
  cmp <- c(cmp, beta_shift = max(abs(db)) < 1e-9)
  say("[rung %g] family gates (M = %d): %s ; max|beta_g - beta_g_MD40 - (q-40)PQg| = %.3e",
      q, fam$M, paste(names(cmp), cmp, collapse = " "), max(abs(db)))
  stopifnot(all(cmp))
  bint <- abs(mt[["Q"]]) - abs(mt[["Qc"]])
  list(dgm = dgm, fam = fam, k_q = k_q, m_tau = mt, bint = bint, gates = cmp)
}

run_rung <- function(q) {
  t0 <- tic(); G <- readRDS(outp("gate"))
  R <- build_rung(q, G); fam <- R$fam
  C1_GRID <- 0:300
  t1 <- tic()
  sw <- fs_oc_grid(family = fam, n = N, c1 = C1_GRID, c2 = C2_DRV,
                   consistency_method = "resample", pconsistency = PCONS,
                   draws = DRAWS, block = Inf, seed = SEED, verbose = TRUE)
  say("[rung %g] grid: %d rows, %.0f s", q, nrow(sw$table), tic() - t1); print(sw$timing)
  c1_05 <- CORR$null$invert$table$value[CORR$null$invert$table$target == 0.05 & CORR$null$invert$table$consistency_method == "resample"]
  c1_10 <- CORR$null$invert$table$value[CORR$null$invert$table$target == 0.10 & CORR$null$invert$table$consistency_method == "resample"]
  named <- list(driver = grid_point(sw, C1_DRV), c1_05 = grid_point(sw, round(c1_05)), c1_10 = grid_point(sw, round(c1_10)))
  t2 <- tic()
  iv <- fs_oc_invert(family = fam, n = N, target = c(0.80, 0.90, 0.95), solve_for = "c1",
                     c2 = C2_DRV, consistency_method = "resample", pconsistency = PCONS,
                     draws = DRAWS, seed = SEED)
  tb <- attr(iv, "table"); tb$secs <- tic() - t2
  tb$grid_c1_nearest <- round(tb$value)
  tb$grid_det_rate_nearest <- sapply(tb$grid_c1_nearest, function(v) if (is.na(v)) NA_real_ else sw$table$det_rate[sw$table$c1 == v])
  print(tb, digits = 5)
  diag <- se_diag(R$dgm, fam)
  say("[rung %g] se_g diagnostic: ratio range [%.5f, %.5f] median %.5f; cor(ratio,Pg) = %.4f; outside [0.992,1.015]: %d ; outside +/-2%%: %d",
      q, diag$overall["min"], diag$overall["max"], diag$overall["median"], diag$cor_Pg, diag$outside_md40_band, diag$outside_2pct)
  print(diag$by_band, digits = 5)
  d <- named$driver
  say("[rung %g] at driver c1 = %g: det %.4f (se %.4f) EnH %.2f Esens %.4f Espec %.4f Eppv %.4f Enpv %.4f EbetaH %.2f naive %.2f mass_below %.4f",
      q, C1_DRV, d$det_rate, d$det_rate_se, d$EnH, d$Esens, d$Espec, d$Eppv, d$Enpv, d$EbetaH, d$Enaive_bias, d$mass_below)
  say("[rung %g] power at c1_05 = %.3f (grid %d): %.4f ; at c1_10 = %.3f (grid %d): %.4f",
      q, c1_05, round(c1_05), named$c1_05$det_rate, c1_10, round(c1_10), named$c1_10$det_rate)
  out <- list(q = q, k_inter = R$k_q, bint = R$bint, m_tau = R$m_tau, gates = R$gates,
              beta_g = fam$beta_g, M = fam$M,
              grid = list(table = sw$table, c1 = C1_GRID, c2 = C2_DRV, block = Inf, seed = SEED, draws = DRAWS, timing = sw$timing),
              named = named, c1_05 = c1_05, c1_10 = c1_10,
              invert = list(table = tb, targets = c(0.80, 0.90, 0.95), objects = iv),
              se_diag = diag[c("overall", "by_band", "cor_Pg", "cor_PQg", "outside_2pct", "outside_md40_band")],
              se_direct = diag$se_direct, ratio = diag$ratio,
              secs = tic() - t0, built_at = Sys.time())
  saveRDS(out, outp(sprintf("rung%g", q)))
  say("[rung %g] written %s (%.0f s)", q, outp(sprintf("rung%g", q)), tic() - t0)
}

# =============================================================================
# §5 + §6  the forecast rung
# =============================================================================
run_forecast <- function(q) {
  t0 <- tic(); G <- readRDS(outp("gate"))
  R <- build_rung(q, G); fam <- R$fam
  nf <- CORR$null$family
  null_fam <- structure(list(M = nf$M, lab = nf$lab, Pg = nf$Pg, PQg = nf$PQg, beta_g = nf$beta_g,
                             se_g = nf$se_g, sens_g = nf$sens_g, spec_g = nf$spec_g, PQ = nf$PQ,
                             null = nf$null), class = "fs_oc_family")
  # the stored null record carries no ovl; re-enumerate the null family and gate it identical
  dgm_null <- generate_glm_dgm(
    data = driver_frame(), factor_vars = paste0("z", 1:12),
    outcome_var = "cd4_change", treatment_var = "treat",
    outcome_type = "continuous", effect_measure = "MD",
    subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
    model = "null", n_super = 5000L, seed = 8316951L, verbose = FALSE)
  nfam <- enumerate(dgm_null)
  ncmp <- sapply(c("lab", "Pg", "PQg", "beta_g", "se_g", "sens_g", "spec_g", "M"), function(f) identical(nfam[[f]], nf[[f]]))
  say("[forecast %g] re-enumerated null family gates vs stored: %s", q, paste(names(ncmp), ncmp, collapse = " ")); print(ncmp)
  stopifnot(all(ncmp))

  C1_GRID <- 0:300
  sw <- list(); named <- list(); inv <- list()
  for (g in c("resample", "split")) {
    t1 <- tic()
    sw[[g]] <- fs_oc_grid(family = fam, n = N, c1 = C1_GRID, c2 = C2_DRV,
                          consistency_method = g, pconsistency = PCONS,
                          draws = DRAWS, block = Inf, seed = SEED, verbose = TRUE)
    say("[forecast %g] grid %s: %.0f s", q, g, tic() - t1)
  }
  c1_05 <- CORR$null$invert$table$value[CORR$null$invert$table$target == 0.05 & CORR$null$invert$table$consistency_method == "resample"]
  c1_10 <- CORR$null$invert$table$value[CORR$null$invert$table$target == 0.10 & CORR$null$invert$table$consistency_method == "resample"]
  for (g in c("resample", "split"))
    named[[g]] <- list(driver = grid_point(sw[[g]], C1_DRV), c1_05 = grid_point(sw[[g]], round(c1_05)), c1_10 = grid_point(sw[[g]], round(c1_10)))
  # c1*: fs_oc_invert at 0.80 under resample (plus 0.90/0.95 for the ladder table)
  t2 <- tic()
  iv <- fs_oc_invert(family = fam, n = N, target = c(0.80, 0.90, 0.95), solve_for = "c1",
                     c2 = C2_DRV, consistency_method = "resample", pconsistency = PCONS,
                     draws = DRAWS, seed = SEED)
  tb <- attr(iv, "table"); tb$secs <- tic() - t2
  tb$grid_c1_nearest <- round(tb$value)
  tb$grid_det_rate_nearest <- sapply(tb$grid_c1_nearest, function(v) if (is.na(v)) NA_real_ else sw$resample$table$det_rate[sw$resample$table$c1 == v])
  print(tb, digits = 5)
  c1_star <- tb$value[tb$target == 0.80]
  say("[forecast %g] c1* = %.4f attainable %s", q, c1_star, tb$attainable[tb$target == 0.80])
  stopifnot(!is.na(c1_star))
  # split-gate inversion too, for the record
  iv_s <- fs_oc_invert(family = fam, n = N, target = c(0.80, 0.90, 0.95), solve_for = "c1",
                       c2 = C2_DRV, consistency_method = "split", pconsistency = PCONS, draws = DRAWS, seed = SEED)
  tb_s <- attr(iv_s, "table"); tb_s$grid_c1_nearest <- round(tb_s$value)
  tb_s$grid_det_rate_nearest <- sapply(tb_s$grid_c1_nearest, function(v) if (is.na(v)) NA_real_ else sw$split$table$det_rate[sw$split$table$c1 == v])
  print(tb_s, digits = 5)
  # the full rows at (500, c1*, c2) under both gates
  star <- list()
  for (g in c("resample", "split")) {
    t3 <- tic()
    p <- fs_oc_predict(family = fam, n = N, c1 = c1_star, c2 = C2_DRV, consistency_method = g,
                       pconsistency = PCONS, draws = DRAWS, seed = SEED)
    p$family <- NULL; star[[g]] <- p
    say("[forecast %g] at c1* %s: det %.4f (se %.4f) EnH %.2f Esens %.4f Espec %.4f Eppv %.4f Enpv %.4f EbetaH %.2f naive %.2f mass_below %.4f (%.0f s)",
        q, g, p$det_rate, p$det_rate_se, p$EnH, p$Esens, p$Espec, p$Eppv, p$Enpv, p$EbetaH, p$Enaive_bias, p$mass_below, tic() - t3)
  }
  # the null false-declaration rate at c1*: the stored null family, resample
  t4 <- tic()
  pn <- fs_oc_predict(family = nfam, n = N, c1 = c1_star, c2 = C2_DRV, consistency_method = "resample",
                      pconsistency = PCONS, draws = DRAWS, seed = SEED)
  pn$family <- NULL
  say("[forecast %g] null at c1*: false declaration %.5f (se %.5f) EnH %.2f (%.0f s)", q, pn$det_rate, pn$det_rate_se, pn$EnH, tic() - t4)
  # null at the driver c1 must reproduce the stored null point (guard)
  # (the stored null runs$resample$pred was fs_oc_predict at c1 = 30; a fresh call here would cost 400 s -- the family identity above is the guard)
  # selection distribution given declaration
  sel <- function(p) {
    o <- order(p$sel_c, decreasing = TRUE)[1:15]
    list(mass_Q = sum(p$sel_c[fam$PQg == 1 & fam$sens_g == 1]),
         mass_PQg95 = sum(p$sel_c[fam$PQg >= 0.95]),
         top15 = data.frame(lab = fam$lab[o], sel_c = p$sel_c[o], Pg = fam$Pg[o], PQg = fam$PQg[o], beta_g = fam$beta_g[o]))
  }
  sel_star <- lapply(star, sel)
  for (g in names(sel_star)) { say("[forecast %g] %s: sel_c mass on Q %.4f ; PQg>=0.95 %.4f", q, g, sel_star[[g]]$mass_Q, sel_star[[g]]$mass_PQg95); print(sel_star[[g]]$top15, digits = 4) }
  # §6 sensitivity: se_direct family
  diag <- se_diag(R$dgm, fam)
  fam_dir <- fam; fam_dir$se_g <- diag$se_direct
  t5 <- tic()
  p6 <- fs_oc_predict(family = fam_dir, n = N, c1 = c1_star, c2 = C2_DRV, consistency_method = "resample",
                      pconsistency = PCONS, draws = DRAWS, seed = SEED)
  p6$family <- NULL
  say("[forecast %g] §6 se_direct family at c1* resample: det %.4f (se %.4f) EnH %.2f Eppv %.4f Esens %.4f EbetaH %.2f mass_Q %.4f (%.0f s)",
      q, p6$det_rate, p6$det_rate_se, p6$EnH, p6$Eppv, p6$Esens, p6$EbetaH, sel(p6)$mass_Q, tic() - t5)
  out <- list(q = q, k_inter = R$k_q, bint = R$bint, m_tau = R$m_tau, gates = R$gates, null_gates = ncmp,
              beta_g = fam$beta_g, M = fam$M,
              grid = lapply(sw, function(s) list(table = s$table, timing = s$timing)),
              grid_settings = list(c1 = C1_GRID, c2 = C2_DRV, block = Inf, seed = SEED, draws = DRAWS),
              named = named, c1_05 = c1_05, c1_10 = c1_10,
              invert = list(resample = tb, split = tb_s, objects = list(resample = iv, split = iv_s)),
              c1_star = c1_star, star = star, null_at_star = pn, sel_star = sel_star,
              sens6 = list(pred = p6, sel = sel(p6), se_direct = diag$se_direct),
              se_diag = diag[c("overall", "by_band", "cor_Pg", "cor_PQg", "outside_2pct", "outside_md40_band")],
              ratio = diag$ratio, secs = tic() - t0, built_at = Sys.time())
  saveRDS(out, outp(sprintf("forecast%g", q)))
  say("[forecast %g] written %s (%.0f s)", q, outp(sprintf("forecast%g", q)), tic() - t0)
}

switch(part,
       gate     = run_gate(),
       rung     = run_rung(as.numeric(args[2])),
       forecast = run_forecast(as.numeric(args[2])))
