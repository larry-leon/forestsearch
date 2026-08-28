# =============================================================================
# invert_reduction_guard_2026-08-29.R -- §3 GATE of the null-path task
# -----------------------------------------------------------------------------
# Re-runs the twelve inversions stored in oc_wrapper_grid_2026-08-29.rds
# (slow path: bisection, one draw set per inversion, ~2.4 h) under the
# order-statistic reduction (one draw set per (n, gate), vector targets), same
# seed and draw count, and compares the inverted c1 values to the stored ones
# at the resolution the draw count supports: the local order-statistic step,
# i.e. the gap between the k-th largest T and its neighbours, and the old
# bisection tolerance (1e-3).  Also records the new wall-clock.
# =============================================================================
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
RDS <- "dev/glm-continuous-sims/oc_wrapper_grid_2026-08-29.rds"
G <- readRDS(RDS); old <- G$invert$table
stopifnot(nrow(old) == 12L)
SEED <- G$seed; DRAWS <- G$draws; fs_args <- G$fs_args
c2 <- fs_args$consistency.threshold

# fixture, as before
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
stopifnot(abs(dgm$hazard_ratios$harm_subgroup - truth$effect_Q) < 1e-9)

t_all <- proc.time()[["elapsed"]]
rows <- list()
for (n in c(500, 700)) {
  fam <- fs_oc_family_enumerate(dgm, fs_args, n = n, max_M = 5000L)
  stopifnot(fam$M == G$families[[as.character(n)]]$M)
  for (g in c("resample", "split")) {
    t0 <- proc.time()[["elapsed"]]
    iv <- fs_oc_invert(family = fam, n = n, target = G$invert$targets,
                       solve_for = "c1", c2 = c2, consistency_method = g,
                       pconsistency = fs_args$pconsistency.threshold,
                       draws = DRAWS, seed = SEED)
    secs <- proc.time()[["elapsed"]] - t0
    # the local order-statistic step around k, from the same draw set
    Sg <- (fam$ovl / sqrt(outer(fam$Pg, fam$Pg))) * outer(fam$se_g, fam$se_g)
    set.seed(SEED); dr <- forestsearch:::.fs_oc_draw(Sg, fam$beta_g, fam$M, DRAWS, g)
    red <- forestsearch:::.fs_oc_reduce(dr, fam, c2, if (g == "resample") fs_args$pconsistency.threshold else NA_real_, g)
    Ts <- sort(red$T, decreasing = TRUE)
    for (o in iv) {
      k <- o$k
      step <- max(Ts[k - 1L] - Ts[k], Ts[k] - Ts[k + 1L])
      so <- old[old$n == n & old$consistency_method == g & abs(old$target - o$target) < 1e-9, ]
      rows[[length(rows) + 1L]] <- data.frame(
        n = n, gate = g, target = o$target, c1_old = so$c1, c1_new = o$value,
        diff = o$value - so$c1, step_local = step, tol_old = 1e-3,
        achieved_old = so$achieved, achieved_new = o$achieved,
        within_one_step = abs(o$value - so$c1) <= max(step, 1e-3),
        secs_new_per_gate = secs, secs_old = so$secs)
    }
    rm(dr, red, Ts); gc()
  }
}
tab <- do.call(rbind, rows); rownames(tab) <- NULL
t_all <- proc.time()[["elapsed"]] - t_all
print(tab, digits = 6)
cat(sprintf("\nlargest |c1_new - c1_old| = %.3e;  largest local step = %.3e;  all within one step or tol: %s\n",
            max(abs(tab$diff)), max(tab$step_local), all(tab$within_one_step)))
cat(sprintf("achieved rates identical to stored: %s (max |diff| %.2e)\n",
            isTRUE(all.equal(tab$achieved_old, tab$achieved_new)), max(abs(tab$achieved_old - tab$achieved_new))))
cat(sprintf("wall-clock: reduction, all twelve inversions incl. 4 draw sets and 2 enumerations: %.0f s;  old slow path: %.0f s\n",
            t_all, sum(old$secs)))
saveRDS(list(table = tab, secs_new_total = t_all, secs_old_total = sum(old$secs)),
        "dev/glm-continuous-sims/invert_reduction_guard_2026-08-29.rds")
cat("REDUCTION GUARD:", if (all(tab$within_one_step)) "PASS" else "FAIL", "\n")
quit(status = if (all(tab$within_one_step)) 0L else 1L)
