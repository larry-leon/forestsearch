# =============================================================================
# oc_wrapper_null_2026-08-29.R -- the null-cell comparison (§5 of the null task)
# -----------------------------------------------------------------------------
# Rebuilds the null DGM deterministically as the driver does under
# null_cell = TRUE (generate_glm_dgm(model = "null"), same ACTG175 frame, same
# factor/subgroup spec, seed 8316951), gated against the tracked null payload's
# committed truth; enumerates the family at n = 500 under the null branch;
# evaluates fs_oc_predict() at the driver's (c1, c2) on both gates at 2e5 draws
# (seed 20260825); then the c1 sweep and the inversions through the
# order-statistic reduction.  The measured column is read from the null
# payload's $oc / $results (computed, never typed; the payload's third element
# is an NA-named NULL -- everything is indexed by name).
# Writes a SIBLING .rds: dev/glm-continuous-sims/oc_wrapper_null_2026-08-29.rds
# =============================================================================
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
stopifnot(requireNamespace("speff2trial", quietly = TRUE))
PAY <- file.path("quarto", "simulations", "actg175", "continuous", "mr_md_harm",
                 "fs_maxeffCons_mr_mdnull_knoise0_n500_s1000_d5000",
                 "fs_maxeffCons_mr_mdnull_knoise0_n500_res_1_1000.rds")
OUT <- file.path("dev", "glm-continuous-sims", "oc_wrapper_null_2026-08-29.rds")
G_ALT <- readRDS(file.path("dev", "glm-continuous-sims", "oc_wrapper_grid_2026-08-29.rds"))
SEED <- G_ALT$seed; DRAWS <- G_ALT$draws
SWEEP_C1 <- G_ALT$sweep$c1                 # the same c1 grid as the alternative
TARGETS  <- c(0.05, 0.10, 0.20, 0.50, 0.80, 0.90, 0.95)

# ---- the null payload, indexed by name (third element is an NA-named NULL) --
pl <- readRDS(PAY)
nm <- names(pl); nm[is.na(nm)] <- ""
pl <- setNames(pl, nm)
truth <- pl[["truth"]]; oc <- pl[["oc"]]; res <- pl[["results"]]; meta <- pl[["meta"]]
stopifnot(!is.null(truth), !is.null(oc), !is.null(res), !is.null(meta))
stopifnot(is.na(truth$effect_Q), truth$prevalence_Q == 0, truth$beta_inter == 0,
          meta$consistency_method == "resample", meta$n_sample == 500)

# ---- the null DGM, as the driver builds it ----------------------------------
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
dgm <- generate_glm_dgm(
  data = actg_df, factor_vars = paste0("z", 1:12), outcome_var = "cd4_change",
  treatment_var = "treat", outcome_type = "continuous", effect_measure = "MD",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "null", n_super = 5000L, seed = 8316951L, verbose = FALSE)
gate <- c(Q_empty      = sum(dgm$df_super$flag_harm) == 0L,
          effect_Q_NA  = is.na(dgm$hazard_ratios$harm_subgroup),
          effect_Qc    = abs(dgm$hazard_ratios$no_harm_subgroup - truth$effect_Qc) < 1e-9,
          beta_inter   = isTRUE(all.equal(dgm$model_params$beta_inter, truth$beta_inter)),
          effect_ITT   = abs(dgm$hazard_ratios$overall - truth$effect_ITT) < 1e-9,
          model_null   = identical(dgm$model, "null"))
cat("null rebuild gate:\n"); print(gate); stopifnot(all(gate))

fs_args <- G_ALT$fs_args
stopifnot(meta$effect_threshold == fs_args$effect.threshold,
          meta$consistency_threshold == fs_args$consistency.threshold)
c1 <- fs_args$effect.threshold; c2 <- fs_args$consistency.threshold
n_trial <- 500

# ---- the family under the null floor ------------------------------------------
t0 <- proc.time()[["elapsed"]]
fam <- fs_oc_family_enumerate(dgm, fs_args, n = n_trial, max_M = 5000L, verbose = TRUE)
t_enum <- proc.time()[["elapsed"]] - t0
print(fam)
stopifnot(isTRUE(fam$null), all(fam$PQg == 0), all(is.na(fam$sens_g)))
cat(sprintf("null family: M = %d; common effect %.6f; S-row V_eff %.4f; alt M at n = 500 was %d\n",
            fam$M, fam$beta_g[1], fam$scale$regions$V_eff, G_ALT$families[["500"]]$M))

# ---- measured column from the payload ----------------------------------------
idc <- oc$identification[oc$identification$convention == "conditional", ]
est <- oc$estimation
nvH <- est[est$block == "H" & est$estimator == "naive", ]
det_rows <- res[res$detected %in% TRUE, ]
measured <- list(
  det_rate = idc$detection, det_rate_se = sqrt(idc$detection * (1 - idc$detection) / meta$n_sims),
  EnH = idc$mean_size_H, EnH_se = stats::sd(det_rows$n_harm, na.rm = TRUE) / sqrt(sum(!is.na(det_rows$n_harm))),
  Esens = idc$sens,                       # NA in the payload: 0/0
  Espec = idc$spec, Espec_se = stats::sd(det_rows$spec, na.rm = TRUE) / sqrt(sum(!is.na(det_rows$spec))),
  Eppv = idc$ppv, Enpv = idc$npv,
  # E[beta(Hhat)]: under the null every rule has the common effect; from results
  EbetaH = mean(oc$targets$orient * det_rows$betaHhat_H, na.rm = TRUE),
  EbetaH_se = stats::sd(det_rows$betaHhat_H, na.rm = TRUE) / sqrt(sum(!is.na(det_rows$betaHhat_H))),
  Enaive_bias = nvH$bias_beta,
  Enaive_bias_se = stats::sd(det_rows$nv_H_est - oc$targets$orient * det_rows$betaHhat_H, na.rm = TRUE) /
                   sqrt(sum(!is.na(det_rows$nv_H_est))),
  n_used = idc$n_used, n_sims = meta$n_sims,
  source = "null payload $oc$identification[conditional], $oc$estimation[H, naive], $results")
cat("measured (null payload): det", measured$det_rate, " EnH", measured$EnH, " spec", measured$Espec,
    " ppv", measured$Eppv, " npv", measured$Enpv, " EbetaH", measured$EbetaH, " naive bias", measured$Enaive_bias, "\n")

# ---- the two single-point evaluations at the driver's (c1, c2) ----------------
runs <- list()
for (g in c("resample", "split")) {
  t0 <- proc.time()[["elapsed"]]
  p <- fs_oc_predict(family = fam, n = n_trial, c1 = c1, c2 = c2,
                     consistency_method = g, pconsistency = fs_args$pconsistency.threshold,
                     draws = DRAWS, seed = SEED)
  secs <- proc.time()[["elapsed"]] - t0
  p1 <- p$P1
  L_eff <- log(1 - p$det_rate) / log(1 - max(p1))    # the document's own definition
  p$family <- NULL
  runs[[g]] <- list(pred = p, secs = secs, p1_range = range(p1), L_eff = L_eff)
  cat(sprintf("%-8s false declaration %.4f (SE %.4f); P1 range %.4f..%.4f; L_eff %.2f; EnH %.1f; EbetaH %.3f; naive %.2f; mass_below %.3f; %.0f s\n",
              g, p$det_rate, p$det_rate_se, min(p1), max(p1), L_eff, p$EnH, p$EbetaH, p$Enaive_bias, p$mass_below, secs))
}

# ---- the c1 sweep (type-I error curve) and the inversions ----------------------
t0 <- proc.time()[["elapsed"]]
sw <- fs_oc_grid(family = fam, n = n_trial, c1 = SWEEP_C1, c2 = c2,
                 consistency_method = c("resample", "split"),
                 pconsistency = fs_args$pconsistency.threshold,
                 draws = DRAWS, block = Inf, seed = SEED, verbose = TRUE)
cat(sprintf("sweep: %d rows in %.0f s\n", nrow(sw$table), proc.time()[["elapsed"]] - t0))
for (g in c("resample", "split")) {
  i <- which(sw$table$consistency_method == g & sw$table$c1 == c1)
  new <- sw$results[[i]]; new$family <- NULL
  cat(sprintf("%-8s sweep point at c1 = %g identical() to fs_oc_predict(): %s\n", g, c1, identical(new, runs[[g]]$pred)))
  stopifnot(identical(new, runs[[g]]$pred))
}
print(sw$timing)

inv <- list(); inv_rows <- list()
for (g in c("resample", "split")) {
  t0 <- proc.time()[["elapsed"]]
  iv <- fs_oc_invert(family = fam, n = n_trial, target = TARGETS, solve_for = "c1", c2 = c2,
                     consistency_method = g, pconsistency = fs_args$pconsistency.threshold,
                     draws = DRAWS, seed = SEED)
  secs <- proc.time()[["elapsed"]] - t0
  tb <- attr(iv, "table"); tb$consistency_method <- g; tb$secs <- secs
  inv[[g]] <- iv; inv_rows[[g]] <- tb
  cat(sprintf("invert %-8s: %d targets in %.0f s (one draw set)\n", g, length(TARGETS), secs))
}
inv_tab <- do.call(rbind, inv_rows); rownames(inv_tab) <- NULL
print(inv_tab[, c("consistency_method", "target", "value", "achieved", "achieved_se", "ceiling", "attainable")], digits = 5)

# ---- do the gates differ? --------------------------------------------------------
cmp <- merge(sw$table[sw$table$consistency_method == "resample", c("c1", "det_rate", "det_rate_se")],
             sw$table[sw$table$consistency_method == "split", c("c1", "det_rate")],
             by = "c1", suffixes = c("_resample", "_split"))
cmp$diff <- cmp$det_rate_split - cmp$det_rate_resample
cat("\nsplit - resample false-declaration rate along c1:\n"); print(cmp, digits = 4, row.names = FALSE)

out <- list(
  runs = runs, family = list(M = fam$M, counts = fam$counts, floor = fam$args_used$n.min / n_trial,
                             Pg = fam$Pg, se_g = fam$se_g, beta_common = fam$beta_g[1],
                             scale = fam$scale, lab = fam$lab, secs = t_enum),
  measured = measured, sweep = list(table = sw$table, c1 = SWEEP_C1, c2 = c2, timing = sw$timing),
  invert = list(table = inv_tab, targets = TARGETS, objects = inv),
  gate_compare = cmp, truth = truth, meta = meta, fs_args = fs_args, n = n_trial,
  seed = SEED, draws = DRAWS, payload = PAY, built_at = Sys.time(),
  pkg_version = as.character(utils::packageVersion("forestsearch")))
saveRDS(out, OUT)
cat("written:", OUT, "\n")
