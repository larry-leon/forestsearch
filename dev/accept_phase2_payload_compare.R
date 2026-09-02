#!/usr/bin/env Rscript
# =========================================================================
# PHASE 2 PAYLOAD COMPARISON GATE
#
# The ported vignettes (run_subgroup_sims() + summary() + plot()) must
# reproduce the committed 10,000-sim results exactly.  For each design,
# this script compares the payload written by the ported render against
# the committed payload (produced by the pre-port vignettes at commit
# 002f4f37):
#
#   matrices sim_hrs / sim_ubs / sim_ns .... identical()
#   design, n_sims, subgroups ............. identical()
#   cens_adjust, k_treat, hr_true ......... identical()
#
# and additionally checks that summary() applied to the *committed*
# payload reproduces the ported document's results table -- the property
# that lets the Phase 3 memo consolidation swap .col_stats() for
# summary() with no numeric change.
#
# Usage:  Rscript dev/accept_phase2_payload_compare.R [results_dir]
#         default results_dir = quarto/extreme_subgroups/fixed_random/results
# =========================================================================
suppressMessages(library(forestsearch))

args <- commandArgs(trailingOnly = TRUE)
res_dir <- if (length(args) >= 1) args[1] else
  "quarto/extreme_subgroups/fixed_random/results"

.fails <- 0L
ok <- function(cond, msg) {
  cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg))
  if (!isTRUE(cond)) .fails <<- .fails + 1L
  invisible(cond)
}

for (design in c("resample", "fixed")) {
  f_old <- file.path(res_dir,
                     sprintf("extreme_sims_%s_10000_payload.rds", design))
  f_new <- file.path(res_dir,
                     sprintf("extreme_sims_%s_10000_ported.rds", design))
  for (f in c(f_old, f_new)) {
    if (!file.exists(f)) {
      stop("missing payload: ", f,
           if (f == f_new) "\n  (render the ported vignette first: Gate 3)"
           else "\n  (committed payload expected from the 002f4f37 renders)")
    }
  }
  old <- readRDS(f_old)
  new <- readRDS(f_new)

  cat(sprintf("\n-- design: %-8s  old: %s | new: %s --\n",
              design, basename(f_old), basename(f_new)))
  ok(identical(old$sim_hrs, new$sim_hrs),
     sprintf("%s: sim_hrs identical (%d x %d)", design,
             nrow(new$sim_hrs), ncol(new$sim_hrs)))
  ok(identical(old$sim_ubs, new$sim_ubs), paste0(design, ": sim_ubs identical"))
  ok(identical(old$sim_ns,  new$sim_ns),  paste0(design, ": sim_ns  identical"))
  ok(identical(old$design, new$design) && identical(old$n_sims, new$n_sims),
     paste0(design, ": design label and n_sims identical"))
  ok(identical(old$subgroups, new$subgroups),
     paste0(design, ": subgroup definitions identical"))
  ok(identical(old$cens_adjust, new$cens_adjust),
     sprintf("%s: cens_adjust identical (%.6f)", design, new$cens_adjust))
  ok(identical(old$k_treat, new$k_treat) && identical(old$hr_true, new$hr_true),
     paste0(design, ": k_treat and hr_true identical"))
  ok(inherits(new, "subgroup_sims"),
     paste0(design, ": ported payload carries the subgroup_sims class"))

  # summary() reproduces on the committed (legacy, class-less) payload:
  # the memo-consolidation property.
  S_old <- summary(structure(old, class = c("subgroup_sims", "list")),
                   hr_true = old$hr_true)
  S_new <- summary(new, hr_true = new$hr_true)
  ok(identical(S_old$results_tbl, S_new$results_tbl),
     paste0(design, ": summary() results table identical on old vs new payload"))
}

cat("\n")
if (.fails == 0L) {
  cat("==== PHASE 2 PAYLOAD COMPARISON: ALL GATES PASSED ====\n")
} else {
  cat(sprintf("==== PHASE 2 PAYLOAD COMPARISON: %d GATE(S) FAILED ====\n",
              .fails))
  quit(status = 1L)
}
