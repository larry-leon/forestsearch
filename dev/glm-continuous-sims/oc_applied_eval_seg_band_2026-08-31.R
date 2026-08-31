# =============================================================================
# oc_applied_eval_seg_band_2026-08-31.R
# Applied OC evaluation §5: the se_g band at q in {0.01, 20, T_obs}.
# For each of the three rungs: fs_dgm_scale(dgm, regions = the 4508 memb
# columns) -> se_direct = sqrt(V_eff[g] / (n * Pg)), ratio to the family's
# prevalence-scaled se_g; range, median, and the three purity bands
# (PQg >= 0.95, mid, < 0.25).  Reported, not acted on (the task adds a
# sensitivity run only if the top rung's band exceeds ~5%).
# Outputs: oc_applied_eval_seg_band_2026-08-31.{log,rds}.
# =============================================================================
suppressPackageStartupMessages({
  library(forestsearch)
  library(speff2trial)
})

out_dir <- "dev/glm-continuous-sims"
scratch <- Sys.getenv("FS_OC_SCRATCH")
stopifnot(dir.exists(out_dir), nzchar(scratch), dir.exists(scratch))

build <- readRDS(file.path(out_dir, "oc_applied_eval_build_2026-08-31.rds"))
q_rungs    <- build$q_rungs
beta_treat <- build$context$beta_treat
n          <- build$context$N

# Data prep, identical to the build script
cont_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
bin_vars  <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
actg_df$id <- seq_len(nrow(actg_df))
actg_df$treat <- ifelse(actg_df$arms == 1L, 1L, 0L)
actg_df <- actg_df[!is.na(actg_df$cd420), ]
actg_df$y_decline <- actg_df$cd40 - actg_df$cd420
for (v in bin_vars) actg_df[[v]] <- as.factor(actg_df[[v]])
stopifnot(nrow(actg_df) == n)

build_dgm <- function(q) {
  generate_glm_dgm(
    data            = actg_df,
    factor_vars     = bin_vars,
    continuous_vars = cont_vars,
    outcome_var     = "y_decline",
    treatment_var   = "treat",
    outcome_type    = "continuous",
    effect_measure  = "MD",
    subgroup_vars   = build$subgroup$vars,
    subgroup_cuts   = build$subgroup$cuts,
    model           = "alt",
    k_treat         = 1,
    k_inter         = q - beta_treat,
    adverse_outcome = TRUE,
    n_super         = 5000L,
    seed            = 8316951L,
    verbose         = FALSE
  )
}

band_rungs <- c(1L, 7L, 11L)      # q = 0.01, 20, T_obs
stopifnot(isTRUE(all.equal(q_rungs[7L], 20)))

res <- list()
for (i in band_rungs) {
  q <- q_rungs[i]
  cat(sprintf("[§5] rung %d (q = %.6f): scale over the %s memb columns...\n",
              i, q, "4508"))
  fam <- readRDS(file.path(scratch, sprintf("fam_rung_%02d.rds", i)))
  stopifnot(fam$M == 4508L)
  dgm <- build_dgm(q)
  regions <- lapply(seq_len(fam$M), function(j) fam$memb[, j])
  names(regions) <- sprintf("g%04d", seq_len(fam$M))
  t0 <- proc.time()[["elapsed"]]
  sc <- fs_dgm_scale(dgm, regions = regions)
  el <- proc.time()[["elapsed"]] - t0
  stopifnot(nrow(sc$regions) == fam$M)
  se_direct <- sqrt(sc$regions$V_eff / (n * fam$Pg))
  ratio <- se_direct / fam$se_g

  band <- cut(fam$PQg, breaks = c(-Inf, 0.25, 0.95, Inf), right = FALSE,
              labels = c("PQg < 0.25", "mid", "PQg >= 0.95"))
  by_band <- do.call(rbind, lapply(levels(band), function(b) {
    r <- ratio[band == b]
    data.frame(band = b, n_members = length(r),
               min = min(r), median = median(r), max = max(r))
  }))
  cat(sprintf("[§5] rung %d: ratio range [%.6f, %.6f], median %.6f (%.1f s)\n",
              i, min(ratio), max(ratio), median(ratio), el))
  print(by_band, digits = 6)
  res[[as.character(i)]] <- list(
    rung = i, q = q, ratio = ratio, PQg = fam$PQg, Pg = fam$Pg,
    se_g = fam$se_g, se_direct = se_direct,
    range = range(ratio), median = median(ratio), by_band = by_band,
    secs = el)
}

top <- res[[as.character(11L)]]
max_dev <- max(abs(top$ratio - 1))
cat(sprintf("\n[§5] top rung max |ratio - 1| = %.4f (%.2f%%) -> sensitivity %s\n",
            max_dev, 100 * max_dev,
            if (max_dev > 0.05) "TRIGGERED (run it)" else "not triggered"))
res$max_dev_top <- max_dev
res$sensitivity_triggered <- max_dev > 0.05

saveRDS(res, file.path(out_dir, "oc_applied_eval_seg_band_2026-08-31.rds"))
cat("[done] se_g band saved.\n")
