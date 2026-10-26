# timing_common.R -- shared timing protocol for stageA.R / stageC.R
# (task dev/tasks/cc_task_lean_fits_2026-10-26.md, sections 2 and 6).
# Sources ../vi_deprecation_2026-10-26/vi_common.R unchanged.  Run from this
# directory.
source("../vi_deprecation_2026-10-26/vi_common.R")

# Per fixture (survival, continuous, binary), "ii" args: one untimed warm-up,
# three timed sequential forestsearch() calls, one nb = 20 sequential
# bootstrap on the warm-up fit.  Saves stage<S>_<fixture>_boot.rds and
# returns the timing table.
run_timing <- function(stage) {
  rows <- list()
  for (fx_name in names(FIXTURES)) {
    fx <- FIXTURES[[fx_name]]()
    args <- variant_args(fx$args, "ii")
    warm <- run_collect(args)
    cat(sprintf("[%s] warm-up done (warnings = %d)\n", fx_name, length(warm$warnings)))
    calls <- vapply(1:3, function(i) {
      system.time(run_collect(args))[["elapsed"]]
    }, numeric(1))
    bt_time <- system.time(bt <- run_boot(warm$fs))[["elapsed"]]
    rep_mean <- mean(bt$results$tmins_iteration) * 60
    saveRDS(extract_boot(bt), file.path(OUT, sprintf("stage%s_%s_boot.rds", stage, fx_name)))
    rows[[fx_name]] <- data.frame(
      fixture = fx_name, call_s_1 = calls[1], call_s_2 = calls[2], call_s_3 = calls[3],
      call_median_s = stats::median(calls), boot_outer_s = bt_time, rep_mean_s = rep_mean)
    cat(sprintf("[%s] calls %.2f %.2f %.2f (median %.2f); boot outer %.1f s, rep mean %.2f s\n",
                fx_name, calls[1], calls[2], calls[3], stats::median(calls), bt_time, rep_mean))
  }
  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  write.csv(tab, file.path(OUT, sprintf("stage%s_timing.csv", stage)), row.names = FALSE)
  tab
}
