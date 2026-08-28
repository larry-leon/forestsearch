# =============================================================================
# oc_wrapper_grid_2026-08-29.R -- four fs_oc_predict() evaluations to one .rds
# -----------------------------------------------------------------------------
# Extends fixture_run_fs_oc_predict_2026-08-28.R (left untouched).  For each of
# n in {500, 700} x gate in {"resample", "split"}:
#   - the fs_oc_predict() object (with its MC SEs), the family's M and stage
#     counts, the floor n.min/n, settings and seed, wall-clock;
# plus, per n, the MEASURED column read from the tracked payload's $oc (and
# $results where $oc has no counterpart) -- computed, never typed.
# Output: dev/glm-continuous-sims/oc_wrapper_grid_2026-08-29.rds
#
# Run from the repository root:
#   Rscript dev/glm-continuous-sims/oc_wrapper_grid_2026-08-29.R
# =============================================================================
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
stopifnot(requireNamespace("speff2trial", quietly = TRUE))

PAY_DIR <- file.path("quarto", "simulations", "actg175", "continuous", "mr_md_harm")
pay_path <- function(n) file.path(
  PAY_DIR, sprintf("fs_maxeffCons_mr_md40_knoise0_n%d_s1000_d5000", n),
  sprintf("fs_maxeffCons_mr_md40_knoise0_n%d_res_1_1000.rds", n))
OUT <- file.path("dev", "glm-continuous-sims", "oc_wrapper_grid_2026-08-29.rds")

SEED  <- 20260825L      # the prediction document's worked-predictions seed
DRAWS <- 2e5

# ---- 1. the fixture (deterministic MD40 rebuild, gated) ---------------------
pl500 <- readRDS(pay_path(500)); pl700 <- readRDS(pay_path(700))
truth <- pl500$truth
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
  model = "alt", k_treat = 1, k_inter = truth$beta_inter,
  n_super = 5000L, seed = 8316951L, verbose = FALSE)
gate <- c(abs(dgm$hazard_ratios$harm_subgroup - truth$effect_Q),
          abs(dgm$subgroup_info$proportion - truth$prevalence_Q),
          abs(dgm$model_params$beta_inter - truth$beta_inter))
stopifnot(all(gate < 1e-9))
sc <- fs_dgm_scale(dgm)
stopifnot(isTRUE(all.equal(sc$regions, pl500$scale$regions)),
          isTRUE(all.equal(sc$regions, pl700$scale$regions)))
cat("fixture rebuilt and gated; scale identical to both payloads\n")

# ---- 2. the driver's arguments ---------------------------------------------
fs_args <- list(
  confounders.name = c("age", "preanti", "wtkg", "karnof", "cd40", "cd80",
                       "hemo", "homo", "drugs", "race", "gender", "symptom"),
  conf.cont_jcuts = list(age = 10, preanti = 10),
  n.min = 60L, maxk = 2L, sg_focus = "maxeffCons",
  effect.threshold = 30, consistency.threshold = 10,
  pconsistency.threshold = 0.90)
for (n in c(500, 700)) {
  pm <- (if (n == 500) pl500 else pl700)$meta
  stopifnot(pm$effect_threshold == fs_args$effect.threshold,
            pm$consistency_threshold == fs_args$consistency.threshold,
            pm$sg_focus == fs_args$sg_focus, pm$n_sample == n,
            pm$consistency_method == "resample")
}

# ---- 3. the measured columns, from the payloads -----------------------------
measured_from_payload <- function(pl) {
  oc <- pl$oc; res <- pl$results
  idc <- oc$identification[oc$identification$convention == "conditional", ]
  est <- oc$estimation
  orH <- est[est$block == "H" & est$estimator == "oracle", ]
  nvH <- est[est$block == "H" & est$estimator == "naive",  ]
  det_rows <- res[res$detected %in% TRUE, ]
  list(
    det_rate    = idc$detection,
    EnH         = idc$mean_size_H,
    Esens       = idc$sens,
    Espec       = idc$spec,
    Eppv        = idc$ppv,
    Enpv        = idc$npv,
    # E[beta(Hhat)] oriented: the oracle row's avg minus its bias against
    # beta(Hhat) -- i.e. the mean of the per-replicate truth on Hhat.  Also
    # taken directly from results$betaHhat_H (oriented by targets$orient).
    EbetaH      = orH$avg - orH$bias_beta,
    EbetaH_results = mean(oc$targets$orient * det_rows$betaHhat_H, na.rm = TRUE),
    Enaive_bias = nvH$bias_beta,
    n_used      = idc$n_used,
    n_sims      = pl$meta$n_sims,
    # MC SEs on the measured side, from the per-replicate records
    det_rate_se = sqrt(idc$detection * (1 - idc$detection) / pl$meta$n_sims),
    EnH_se      = stats::sd(det_rows$n_harm, na.rm = TRUE) / sqrt(sum(!is.na(det_rows$n_harm))),
    Esens_se    = stats::sd(det_rows$sens, na.rm = TRUE) / sqrt(sum(!is.na(det_rows$sens))),
    Espec_se    = stats::sd(det_rows$spec, na.rm = TRUE) / sqrt(sum(!is.na(det_rows$spec))),
    Eppv_se     = stats::sd(det_rows$ppv,  na.rm = TRUE) / sqrt(sum(!is.na(det_rows$ppv))),
    Enpv_se     = stats::sd(det_rows$npv,  na.rm = TRUE) / sqrt(sum(!is.na(det_rows$npv))),
    EbetaH_se   = stats::sd(det_rows$betaHhat_H, na.rm = TRUE) / sqrt(sum(!is.na(det_rows$betaHhat_H))),
    Enaive_bias_se = stats::sd(det_rows$nv_H_est - oc$targets$orient * det_rows$betaHhat_H, na.rm = TRUE) /
                     sqrt(sum(!is.na(det_rows$nv_H_est))),
    source = "payload$oc$identification[conditional], payload$oc$estimation[H: oracle, naive], payload$results",
    M = NA_integer_)
}
measured <- list(`500` = measured_from_payload(pl500), `700` = measured_from_payload(pl700))
cat("measured columns read from payload$oc / $results\n")

# ---- 4. the four evaluations ------------------------------------------------
runs <- list(); families <- list()
for (n in c(500, 700)) {
  t0 <- proc.time()[3]
  fam <- fs_oc_family_enumerate(dgm, fs_args, n = n, max_M = 5000L, verbose = TRUE)
  t_fam <- proc.time()[3] - t0
  families[[as.character(n)]] <- list(
    M = fam$M, counts = fam$counts, n.min = fam$args_used$n.min,
    floor_Pg = fam$args_used$n.min / n, Pg_range = range(fam$Pg),
    PQ = fam$PQ, cuts = fam$cuts, lab = fam$lab, Pg = fam$Pg, PQg = fam$PQg,
    beta_g = fam$beta_g, se_g = fam$se_g, secs = t_fam)
  for (gate in c("resample", "split")) {
    gc()
    t0 <- proc.time()[3]
    mem0 <- sum(gc()[, 2])
    pred <- fs_oc_predict(family = fam, n = n,
                          c1 = fs_args$effect.threshold,
                          c2 = fs_args$consistency.threshold,
                          consistency_method = gate,
                          pconsistency = fs_args$pconsistency.threshold,
                          draws = DRAWS, seed = SEED)
    secs <- proc.time()[3] - t0
    peak <- sum(gc()[, 6])          # max used since last gc reset (MB)
    pred$family <- NULL             # the family is stored once, above
    runs[[sprintf("n%d_%s", n, gate)]] <- list(
      n = n, gate = gate, pred = pred, secs = secs, peak_MB = peak,
      settings = pred$settings, seed = SEED, draws = DRAWS)
    cat(sprintf("n = %d  gate = %-8s  M = %d  det_rate = %.4f  EnH = %.1f  %.0f s  peak %.0f MB\n",
                n, gate, fam$M, pred$det_rate, pred$EnH, secs, peak))
  }
}

out <- list(
  runs = runs, families = families, measured = measured,
  fs_args = fs_args, seed = SEED, draws = DRAWS,
  scale = sc, payloads = c(`500` = pay_path(500), `700` = pay_path(700)),
  built_at = Sys.time(),
  pkg_version = as.character(utils::packageVersion("forestsearch")))
saveRDS(out, OUT)
cat("written:", OUT, "\n")
