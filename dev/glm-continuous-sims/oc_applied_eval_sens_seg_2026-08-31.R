# =============================================================================
# oc_applied_eval_sens_seg_2026-08-31.R
# Applied OC evaluation §5, the triggered sensitivity: the top rung's band of
# se_direct/se_g exceeded ~5% (max 9.14%), so the top rung (q = T_obs) is
# re-evaluated once with the family's se_g replaced by the member-direct
# se_direct = sqrt(V_eff[g]/(n*Pg)) -- one draw set, the same §4 grid and
# named points, saved beside the primary.  Labeled exactly as the breadth
# ladder labeled it: a sensitivity, not a constructor, not adopted.
# =============================================================================
suppressPackageStartupMessages(library(forestsearch))

out_dir <- "dev/glm-continuous-sims"
scratch <- Sys.getenv("FS_OC_SCRATCH")
stopifnot(dir.exists(out_dir), nzchar(scratch), dir.exists(scratch))

build <- readRDS(file.path(out_dir, "oc_applied_eval_build_2026-08-31.rds"))
band  <- readRDS(file.path(out_dir, "oc_applied_eval_seg_band_2026-08-31.rds"))
T_obs <- build$context$T_obs
n     <- build$context$N
q     <- build$q_rungs[11L]

fam <- readRDS(file.path(scratch, "fam_rung_11.rds"))
stopifnot(fam$M == 4508L, fam$orientation$s == 1)
se_direct <- band[["11"]]$se_direct
stopifnot(length(se_direct) == fam$M,
          isTRUE(all.equal(band[["11"]]$se_g, fam$se_g)))
fam$se_g <- se_direct
cat(sprintf("[sens] q = %.9f: se_g replaced by se_direct (max |ratio-1| = %.4f)\n",
            q, band$max_dev_top))

c1_grid <- sort(c(0:200, T_obs))
t0 <- proc.time()[["elapsed"]]
g <- fs_oc_grid(family = fam, n = n, c1 = c1_grid, c2 = 5,
                consistency_method = "resample", pconsistency = 0.90,
                draws = 2e5, block = 5e4, seed = 8316951L, verbose = TRUE)
t_grid <- proc.time()[["elapsed"]] - t0
cat(sprintf("[sens] grid wall-clock: %.1f s\n", t_grid))

tb <- g$table
get_pred <- function(c1v) {
  j <- which(abs(tb$c1 - c1v) < 1e-9)
  stopifnot(length(j) == 1L)
  g$results[[j]]
}
strip <- function(p) {
  p$family$memb  <- NULL
  p$family$ovl   <- NULL
  p$family$scale <- NULL
  p$family$cuts  <- NULL
  p
}
named <- list(analyst = strip(get_pred(10)), T_obs = strip(get_pred(T_obs)))

np_file <- file.path(scratch, "named_points.rds")
t_wait <- 0
while (!file.exists(np_file) && t_wait < 10800) {
  Sys.sleep(30); t_wait <- t_wait + 30
}
np <- if (file.exists(np_file)) readRDS(np_file) else NULL
missing_named <- is.null(np) || is.na(np$c1_05) || is.na(np$c1_10)
if (!missing_named) {
  named$c1_05 <- strip(get_pred(np$c1_05))
  named$c1_10 <- strip(get_pred(np$c1_10))
}

out <- list(
  label = "sensitivity: se_g -> se_direct at the top rung; not adopted",
  rung = 11L, q = q, n = n, T_obs = T_obs,
  table = tb, named = named, named_points = np,
  missing_named = missing_named,
  timing = g$timing, t_grid = t_grid, settings = g$settings)
saveRDS(out, file.path(out_dir, "oc_applied_eval_sens_seg_2026-08-31.rds"))
cat(sprintf("[sens] DONE (%.1f s grid; checkpoint saved)\n", t_grid))
