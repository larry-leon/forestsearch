# Pre-change baseline of fs_sim_bias_coverage() on the seven s7/map1 bundles, BOTH blocks,
# all estimators the function accepts, both sides -- saved for a byte-identity comparison
# after the add-only scale argument lands.
suppressPackageStartupMessages(library(forestsearch))
files <- sort(unlist(lapply(c("results/fs_*_s7_combined_1_2000.rds", "results/fs_*_map1_combined_1_2000.rds"), Sys.glob)))
out <- list()
for (f in files) {
  r <- readRDS(f)$results
  for (blk in c("H", "Hc")) for (sd in c("lower", "upper")) {
    key <- paste(basename(f), blk, sd, sep = "|")
    out[[key]] <- tryCatch(suppressMessages(fs_sim_bias_coverage(r, block = blk,
                     estimators = c("naive", "mr", "fld", "mr_w", "mr_wf"), side = sd)),
                     error = function(e) conditionMessage(e))
  }
}
saveRDS(out, file.path(Sys.getenv("SP"), "bias_coverage_baseline_prechange.rds"))
cat("saved", length(out), "baseline tables from", length(files), "bundles\n")
