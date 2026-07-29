# =============================================================================
# Phase 0 -- cross-check: does the evaluated prefix contain the global winner?
# =============================================================================
# The divergence measurement in phase0_summarise.R rests on one claim:
#
#   under the as-run early-stopping configuration, the qualifier with the
#   smallest candidate index m among the EVALUATED candidates is also the
#   smallest-m qualifier over the FULL candidate family.
#
# It follows from the preview ordering (setorder(-HR, K), so m ascending is
# effect descending) plus the fact that every candidate before the stop WAS
# evaluated.  This script checks it against execution rather than argument, by
# re-running the same seeds with stop_threshold = NULL -- which forces the
# entire family to be evaluated -- and comparing the smallest-m qualifier.
#
# It also quantifies what the early stop never looked at: the number of
# qualifiers in the full family versus the number the as-run pass saw.
#
# Usage (from the repository root):
#   Rscript dev/sg-focus-work/R/phase0_crosscheck.R

source("dev/sg-focus-work/R/phase0_cells.R")

load_pool <- function(pattern) {
  files <- Sys.glob(pattern)
  if (!length(files)) return(NULL)
  bundles <- lapply(files, readRDS)
  cells <- vapply(bundles, `[[`, character(1), "cell")
  out <- list()
  for (cl in unique(cells)) {
    r <- do.call(rbind, lapply(bundles[cells == cl], `[[`, "results"))
    r <- r[order(r$sim_id), ]
    out[[cl]] <- r[!duplicated(r$sim_id), ]
  }
  out
}

asrun <- load_pool("dev/sg-focus-work/out/run_*.rds")
full  <- load_pool("dev/sg-focus-work/out/xchk_[a-z]*.rds")
if (is.null(full)) stop("no cross-check runs found")

rows <- list()
for (cl in intersect(names(asrun), names(full))) {
  a <- asrun[[cl]]; f <- full[[cl]]
  ids <- intersect(a$sim_id, f$sim_id)
  a <- a[match(ids, a$sim_id), ]; f <- f[match(ids, f$sim_id), ]

  # Sanity: the cross-check really did evaluate everything.
  full_eval <- all(f$n_cand_evaluated == f$n_cand_total, na.rm = TRUE)

  # Compare only where the as-run pass found at least one qualifier.
  k <- !is.na(a$first_qual_m) & !is.na(f$first_qual_m)

  rows[[cl]] <- data.frame(
    cell = cl,
    n_compared = sum(k),
    full_eval_ok = full_eval,
    # THE claim under test
    first_qual_m_agrees = sum(a$first_qual_m[k] == f$first_qual_m[k]),
    first_qual_hr_agrees = sum(abs(a$first_qual_hr[k] - f$first_qual_hr[k]) < 1e-9),
    # what early stopping never saw
    med_qual_seen  = stats::median(a$n_qual[k]),
    med_qual_total = stats::median(f$n_qual[k]),
    med_qual_unseen = stats::median(f$n_qual[k] - a$n_qual[k]),
    max_qual_total = max(f$n_qual[k]),
    med_cand_eval_asrun = stats::median(a$n_cand_evaluated[k]),
    med_cand_total = stats::median(f$n_cand_total[k]),
    stringsAsFactors = FALSE)
}

tab <- do.call(rbind, rows)
ord <- match(names(PHASE0_CELLS), tab$cell); tab <- tab[ord[!is.na(ord)], ]

cat("\n=========== CROSS-CHECK: evaluated prefix vs full candidate family ===========\n")
print(tab, row.names = FALSE, digits = 4)
cat("\n'first_qual_m_agrees' == 'n_compared' confirms that the smallest-m\n",
    "qualifier seen under early stopping is the smallest-m qualifier of the\n",
    "whole family -- i.e. the divergence measurement needs no second run.\n", sep = "")

saveRDS(tab, "dev/sg-focus-work/out/summary_crosscheck.rds")
cat("\nsaved -> dev/sg-focus-work/out/summary_crosscheck.rds\n")
