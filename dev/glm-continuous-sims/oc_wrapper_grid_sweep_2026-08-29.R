# =============================================================================
# oc_wrapper_grid_sweep_2026-08-29.R -- the c1 sweep and the inversions on the
# MD40 fixture, appended to oc_wrapper_grid_2026-08-29.rds
# -----------------------------------------------------------------------------
# Reads the 08-29 precompute (families, measured columns, settings), rebuilds
# the fixture exactly as oc_wrapper_grid_2026-08-29.R does, and adds
#   $sweep   : fs_oc_grid() over c1 = SWEEP_C1 at the driver's c2, both gates,
#              n in {500, 700}, one draw set per (n, gate), block = Inf
#              (exact), seed = the document's;
#   $invert  : fs_oc_invert() for declaration targets 0.80 / 0.90 / 0.95,
#              per n and gate, solving for c1 at the driver's c2;
#   $sweep_timing : enumerate / draw / sweep seconds per (n, gate).
# The four single-point evaluations of the 08-29 grid are left as they are.
# Writes the SAME .rds (extended); the verification document reads it.
# =============================================================================
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
stopifnot(requireNamespace("speff2trial", quietly = TRUE))
RDS <- file.path("dev", "glm-continuous-sims", "oc_wrapper_grid_2026-08-29.rds")
G <- readRDS(RDS)
stopifnot(!is.null(G$runs$n500_resample), !is.null(G$measured[["700"]]))

SEED  <- G$seed;  DRAWS <- G$draws
SWEEP_C1 <- seq(20, 120, by = 5)
TARGETS  <- c(0.80, 0.90, 0.95)

# ---- fixture (identical construction to oc_wrapper_grid_2026-08-29.R) --------
pl500 <- readRDS(G$payloads[["500"]]); truth <- pl500$truth
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
stopifnot(abs(dgm$hazard_ratios$harm_subgroup - truth$effect_Q) < 1e-9,
          isTRUE(all.equal(fs_dgm_scale(dgm)$regions, G$scale$regions)))
fs_args <- G$fs_args
c2 <- fs_args$consistency.threshold

# ---- the sweep: one draw set per (n, gate), thresholds swept for free -------
t0 <- proc.time()[["elapsed"]]
sw <- fs_oc_grid(dgm, fs_args, n = c(500, 700), c1 = SWEEP_C1, c2 = c2,
                 consistency_method = c("resample", "split"),
                 pconsistency = fs_args$pconsistency.threshold,
                 draws = DRAWS, block = Inf, seed = SEED, verbose = TRUE,
                 max_M = 5000L)
cat(sprintf("sweep: %d rows in %.0f s\n", nrow(sw$table), proc.time()[["elapsed"]] - t0))
for (k in names(sw$families)) stopifnot(sw$families[[k]]$M == G$families[[k]]$M)
# the sweep's point at the driver's c1 must be identical() to the 08-29 grid's
for (n in c(500, 700)) for (g in c("resample", "split")) {
  i <- which(sw$table$n == n & sw$table$consistency_method == g &
             sw$table$c1 == fs_args$effect.threshold)
  old <- G$runs[[sprintf("n%d_%s", n, g)]]$pred; new <- sw$results[[i]]; new$family <- NULL
  cat(sprintf("n = %d %-8s sweep point at c1 = 30 identical() to 08-29 fs_oc_predict(): %s\n",
              n, g, identical(old, new)))
  stopifnot(identical(old, new))
}
print(sw$timing)

# ---- the inversions ----------------------------------------------------------
inv_rows <- list(); inv_objs <- list()
for (n in c(500, 700)) {
  fam <- sw$families[[as.character(n)]]$family
  for (g in c("resample", "split")) for (tg in TARGETS) {
    t0 <- proc.time()[["elapsed"]]
    iv <- fs_oc_invert(family = fam, n = n, target = tg, solve_for = "c1", c2 = c2,
                       consistency_method = g,
                       pconsistency = fs_args$pconsistency.threshold,
                       draws = DRAWS, seed = SEED)
    secs <- proc.time()[["elapsed"]] - t0
    inv_objs[[sprintf("n%d_%s_%.2f", n, g, tg)]] <- iv
    inv_rows[[length(inv_rows) + 1L]] <- data.frame(
      n = n, consistency_method = g, target = tg, c1 = iv$value,
      achieved = iv$achieved, achieved_se = iv$achieved_se,
      next_step_rate = if (iv$attainable) iv$next_step_rate else NA_real_,
      ceiling = iv$ceiling, ceiling_se = iv$ceiling_se,
      binding = iv$binding, attainable = iv$attainable,
      iterations = iv$iterations, secs = secs, stringsAsFactors = FALSE)
    cat(sprintf("invert n = %d %-8s target %.2f -> c1 = %s achieved %.4f (ceiling %.4f) %.0f s\n",
                n, g, tg, if (is.na(iv$value)) "NA" else sprintf("%.3f", iv$value),
                iv$achieved, iv$ceiling, secs))
  }
}
inv_tab <- do.call(rbind, inv_rows); rownames(inv_tab) <- NULL
print(inv_tab[, c("n", "consistency_method", "target", "c1", "achieved", "achieved_se", "ceiling", "attainable")])

G$sweep <- list(table = sw$table, c1 = SWEEP_C1, c2 = c2, block = Inf,
                seed = SEED, draws = DRAWS)
G$sweep_timing <- sw$timing
G$invert <- list(table = inv_tab, objects = inv_objs, targets = TARGETS)
G$extended_at <- Sys.time()
G$pkg_version_sweep <- as.character(utils::packageVersion("forestsearch"))
saveRDS(G, RDS)
cat("extended:", RDS, "\n")
