# 14-point identity fixture for fs_sim_bias_coverage() (default scale), reproducing
# summary_bias_coverage.qmd's read chunk and comparing to dev/tasks/bias_coverage_points.csv.
suppressPackageStartupMessages(library(forestsearch))
globs <- c("results/fs_*_s7_combined_1_2000.rds", "results/fs_*_map1_combined_1_2000.rds")
cell_label <- function(meta) {
  hr <- meta$target_hr_harm
  base <- if (hr < 1) sprintf("protective %.2f", hr) else if (hr == 1) "null 1.0" else sprintf("harm %.4g", hr)
  kn <- if (!is.null(meta$k_random_noise) && meta$k_random_noise > 0) sprintf(", %d noise", meta$k_random_noise) else ""
  nn <- if (is.null(meta$n_sample) || meta$n_sample == 500) if (nzchar(kn) || hr < 1) "" else ", n=500" else sprintf(", n=%d", meta$n_sample)
  paste0(base, nn, kn)
}
files <- sort(unlist(lapply(globs, Sys.glob)))
tab <- do.call(rbind, lapply(files, function(f) {
  b <- readRDS(f); m <- b$meta
  if (is.null(m$target_hr_harm)) m$target_hr_harm <- as.numeric(sub("^.*_m1_h(\\d{3})_.*$", "\\1", basename(f))) / 100
  if (is.null(m$k_random_noise)) m$k_random_noise <- as.integer(sub("^.*_knoise(\\d+)_.*$", "\\1", basename(f)))
  s <- fs_sim_bias_coverage(b$results, block = "H", estimators = c("mr", "fld"))
  s$cell <- cell_label(m); s
}))
tab$method <- ifelse(tab$estimator == "mr", "IJ", "field")
fx <- read.csv("../../../dev/tasks/bias_coverage_points.csv", stringsAsFactors = FALSE)
m <- merge(fx, tab[, c("cell", "method", "b", "r", "cov1", "cov2", "cov1_ref", "cov2_ref")], by = c("cell", "method"), suffixes = c("_fx", ""))
stopifnot(nrow(m) == 14L)
d <- data.frame(cell = m$cell, method = m$method,
                db = abs(m$b - m$b_fx), dr = abs(m$r - m$r_fx), dc1 = abs(m$cov1 - m$c1), dc2 = abs(m$cov2 - m$c2),
                dref1 = abs(m$cov1_ref - m$c1_pred), dref2 = abs(m$cov2_ref - m$c2_pred))
print(d, digits = 4)
cat(sprintf("\nFIXTURE (%d rows): max |b| %.4f, |r| %.4f, |cov1| %.4f, |cov2| %.4f (tol 0.005); max |ref1| %.2e, |ref2| %.2e (tol 0.001) -> %s\n",
            nrow(m), max(d$db), max(d$dr), max(d$dc1), max(d$dc2), max(d$dref1), max(d$dref2),
            if (max(d$db, d$dr, d$dc1, d$dc2) <= 0.005 && max(d$dref1, d$dref2) <= 0.001) "PASS" else "FAIL"))
