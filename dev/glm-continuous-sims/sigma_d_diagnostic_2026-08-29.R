# =============================================================================
# sigma_d_diagnostic_2026-08-29.R -- §2.3 of the gate-and-n700 task
# -----------------------------------------------------------------------------
# ONE simulated trial (n = 500, seed stated) from the deterministic MD40 DGM.
# For ~10 candidate subgroups spanning the prevalence range, compare
#   sigma_D   : what consistency_resample()'s own machinery computes
#   se_model  : classical SE of the treatment coefficient, subgroup lm()
#   se_hc0    : sqrt(sum(dfbeta[, treat]^2)) computed directly from
#               stats::lm.influence()$coefficients (not the package's helper)
#   se_g      : the wrapper's population SE for that candidate at n = 500
# This is an identity check, not a result.  Nothing here is a reported figure.
# =============================================================================
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))

PAY <- file.path("quarto", "simulations", "actg175", "continuous", "mr_md_harm",
                 "fs_maxeffCons_mr_md40_knoise0_n500_s1000_d5000",
                 "fs_maxeffCons_mr_md40_knoise0_n500_res_1_1000.rds")
truth <- readRDS(PAY)$truth
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
stopifnot(abs(dgm$hazard_ratios$harm_subgroup - truth$effect_Q) < 1e-9)

fs_args <- list(
  confounders.name = c("age", "preanti", "wtkg", "karnof", "cd40", "cd80",
                       "hemo", "homo", "drugs", "race", "gender", "symptom"),
  conf.cont_jcuts = list(age = 10, preanti = 10),
  n.min = 60L, maxk = 2L, sg_focus = "maxeffCons",
  effect.threshold = 30, consistency.threshold = 10)
fam <- fs_oc_family_enumerate(dgm, fs_args, n = 500, max_M = 5000L)

# ---- one trial ---------------------------------------------------------------
SEED_TRIAL <- 20260829L
trial <- simulate_from_glm_dgm(dgm, n = 500, seed = SEED_TRIAL)
cat("one trial: n =", nrow(trial), " seed =", SEED_TRIAL,
    " treat table:", paste(table(trial$treat_sim), collapse = "/"), "\n")

# ---- ~10 candidates spanning the prevalence range ----------------------------
ord <- order(fam$Pg)
pick <- unique(round(seq(1, fam$M, length.out = 10)))
pick <- ord[pick]
memb_trial <- fam$memb[trial$id, pick, drop = FALSE]   # sample membership via id

c2 <- fs_args$consistency.threshold
rows <- lapply(seq_along(pick), function(j) {
  m <- pick[j]
  sub <- trial[memb_trial[, j], ]
  # (1) the package's own path, exactly as the GLM resample call site invokes it
  rr <- consistency_resample(
    sub, method = "closed", outcome_type = "continuous", effect_measure = "MD",
    comparison_threshold = c2, treat.name = "treat_sim", outcome.name = "y_sim",
    adverse_outcome = FALSE)
  # (2) classical SE, oriented as the package orients (adverse_outcome = FALSE
  #     negates Y); the SE is sign-invariant
  fit <- stats::lm(I(-y_sim) ~ treat_sim, data = sub)
  se_model <- summary(fit)$coefficients["treat_sim", "Std. Error"]
  # (3) direct dfbeta from base R's influence (change in coef when obs i dropped)
  dfb <- stats::lm.influence(fit)$coefficients[, "treat_sim"]
  se_hc0 <- sqrt(sum(dfb^2))
  # (3b) plain HC0 sandwich for reference (no 1/(1-h) factor)
  X <- stats::model.matrix(fit); A <- solve(crossprod(X))
  inf0 <- (X %*% A)[, "treat_sim"] * stats::residuals(fit)
  se_sand0 <- sqrt(sum(inf0^2))
  data.frame(candidate = fam$lab[m], Pg_pop = fam$Pg[m], n_g = nrow(sub),
             beta_hat = rr$beta_hat, sigma_D = rr$sigma_D,
             se_model = se_model, se_hc0 = se_hc0, se_hc0_plain = se_sand0,
             se_g_pop = fam$se_g[m],
             r_model = rr$sigma_D / se_model, r_hc0 = rr$sigma_D / se_hc0,
             r_se_g = rr$sigma_D / fam$se_g[m],
             rate_closed = rr$rate_closed, stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows)
tab <- tab[order(tab$Pg_pop), ]
cat("\n=== sigma_D identification on one trial (n = 500) ===\n")
print(tab, digits = 5, row.names = FALSE)
cat(sprintf("\nsigma_D / se_hc0   : range %.10f .. %.10f  (identical to 1e-10: %s)\n",
            min(tab$r_hc0), max(tab$r_hc0), all(abs(tab$r_hc0 - 1) < 1e-10)))
cat(sprintf("sigma_D / se_model : range %.4f .. %.4f;  Spearman with Pg: %.3f\n",
            min(tab$r_model), max(tab$r_model), stats::cor(tab$r_model, tab$Pg_pop, method = "spearman")))
cat(sprintf("sigma_D / se_g_pop : range %.4f .. %.4f;  Spearman with Pg: %.3f\n",
            min(tab$r_se_g), max(tab$r_se_g), stats::cor(tab$r_se_g, tab$Pg_pop, method = "spearman")))
cat(sprintf("qnorm((1 + 0.90)/2) = %.6f\n", stats::qnorm(0.95)))
cat(sprintf("rate_closed vs 2*pnorm((beta_hat - c2)/sigma_D) - 1 (clamped at 0): max |diff| = %.2e\n",
            max(abs(tab$rate_closed - pmax(0, 2 * stats::pnorm((tab$beta_hat - c2) / tab$sigma_D) - 1)))))
