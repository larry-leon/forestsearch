suppressPackageStartupMessages(library(forestsearch))
base <- readRDS(file.path(Sys.getenv("SP"), "bias_coverage_baseline_prechange.rds"))
files <- sort(unlist(lapply(c("results/fs_*_s7_combined_1_2000.rds", "results/fs_*_map1_combined_1_2000.rds"), Sys.glob)))
n_id <- 0L; bad <- character(0)
for (f in files) {
  r <- readRDS(f)$results
  for (blk in c("H", "Hc")) for (sd in c("lower", "upper")) {
    key <- paste(basename(f), blk, sd, sep = "|")
    now <- tryCatch(suppressMessages(fs_sim_bias_coverage(r, block = blk,
                     estimators = c("naive", "mr", "fld", "mr_w", "mr_wf"), side = sd)),
                    error = function(e) conditionMessage(e))
    if (identical(now, base[[key]])) n_id <- n_id + 1L else bad <- c(bad, key)
  }
}
cat(sprintf("BYTE-IDENTITY (default scale) vs pre-change baseline: %d of %d tables identical()%s\n",
            n_id, length(base), if (length(bad)) paste0("; DIFFER: ", paste(bad, collapse = ", ")) else ""))
