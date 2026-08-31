# =============================================================================
# oc_applied_eval_rung_2026-08-31.R  <rung index 1..11>
# Applied OC evaluation, one rung (task cc_task_oc_applied_evaluation_
# 2026-08-31.md §4): one 2e5-draw fs_oc_grid() over the full c1 sweep at the
# analyst's gate, on the family stashed by oc_applied_eval_build_2026-08-31.R.
# Checkpoints the grid table plus the full fs_oc_predict objects at the named
# points only (the analyst's c1 = 10, c1 = T_obs, and c1_05 / c1_10 once the
# zero-plus rung publishes them); per-threshold results lists are not saved,
# and the heavy family carriers (memb, ovl, scale, cuts) are dropped from the
# saved predict objects -- every numeric OC quantity is kept.
# Rung 1 (q = 0.01, the zero-plus structural null) derives c1_05 / c1_10 from
# its own table and publishes them for the other rungs.
# =============================================================================
suppressPackageStartupMessages(library(forestsearch))

args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1L])
stopifnot(length(i) == 1L, !is.na(i), i >= 1L, i <= 11L)

out_dir <- "dev/glm-continuous-sims"
scratch <- Sys.getenv("FS_OC_SCRATCH")
stopifnot(dir.exists(out_dir), nzchar(scratch), dir.exists(scratch))

build <- readRDS(file.path(out_dir, "oc_applied_eval_build_2026-08-31.rds"))
q_rungs <- build$q_rungs
q     <- q_rungs[i]
T_obs <- build$context$T_obs
n     <- build$context$N

fam <- readRDS(file.path(scratch, sprintf("fam_rung_%02d.rds", i)))
stopifnot(fam$M == 4508L, fam$orientation$s == 1)
cat(sprintf("[rung %d] q = %.9f  M = %d  s = +1  n = %d\n", i, q, fam$M, n))

c1_grid <- sort(c(0:200, T_obs))

t0 <- proc.time()[["elapsed"]]
g <- fs_oc_grid(family = fam, n = n, c1 = c1_grid, c2 = 5,
                consistency_method = "resample", pconsistency = 0.90,
                draws = 2e5, block = 5e4, seed = 8316951L, verbose = TRUE)
t_grid <- proc.time()[["elapsed"]] - t0
cat(sprintf("[rung %d] grid wall-clock: %.1f s\n", i, t_grid))

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
smallest_int_c1 <- function(rate_max) {
  ok <- tb$c1 %% 1 == 0 & tb$det_rate <= rate_max
  if (!any(ok)) NA_real_ else min(tb$c1[ok])
}
if (i == 1L) {
  # the zero-plus rung defines the fixed-type-I thresholds (1-unit resolution)
  np <- list(c1_05 = smallest_int_c1(0.05), c1_10 = smallest_int_c1(0.10))
  tmp <- paste0(np_file, ".tmp")
  saveRDS(np, tmp)
  file.rename(tmp, np_file)
  cat(sprintf("[rung 1] published named points: c1_05 = %s, c1_10 = %s\n",
              format(np$c1_05), format(np$c1_10)))
} else {
  # wait for the zero-plus rung (launched concurrently; finishes ~together)
  t_wait <- 0
  while (!file.exists(np_file) && t_wait < 10800) {
    Sys.sleep(30); t_wait <- t_wait + 30
  }
  np <- if (file.exists(np_file)) readRDS(np_file) else NULL
  cat(sprintf("[rung %d] named points %s after %d s wait\n",
              i, if (is.null(np)) "MISSING" else "received", t_wait))
}

missing_named <- is.null(np) || is.na(np$c1_05) || is.na(np$c1_10)
if (!missing_named) {
  named$c1_05 <- strip(get_pred(np$c1_05))
  named$c1_10 <- strip(get_pred(np$c1_10))
}

out <- list(
  rung = i, q = q, n = n, T_obs = T_obs,
  table = tb,
  named = named,
  named_points = np,
  missing_named = missing_named,
  timing = g$timing, t_grid = t_grid,
  settings = g$settings
)
saveRDS(out, file.path(out_dir,
                       sprintf("oc_applied_eval_rung%02d_2026-08-31.rds", i)))
cat(sprintf("[rung %d] DONE  (%.1f s grid; checkpoint saved)\n", i, t_grid))
