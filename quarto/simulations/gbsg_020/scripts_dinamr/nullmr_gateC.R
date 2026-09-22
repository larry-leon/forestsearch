#!/usr/bin/env Rscript
# GATE C -- nullmr (TASK_null_gbsg_mr_2026-09-21), one per cell: the three
# identifiers ran on IDENTICAL draws.  nullthr_gateC.R as is; the ONE change is
# the bundle path, which with MR on carries no _nomr token (template :481).
#
# The fingerprint is itt_est / itt_se, the whole-trial Cox fit of treatment
# alone, recorded on EVERY replicate before the search runs.  It is computed
# from the simulated data only, so it is engine-independent, and unlike
# nullid's or_Hc_* it does not depend on any engine declaring -- the check is
# on all rows, at any declaration rate.  Invariants only.
#
# usage: nullthr_gateC.R <cell> <hr> <n> <c1> <c2> <campaign> <reps> [quickrun TRUE|FALSE]
args <- commandArgs(trailingOnly = TRUE)
CELL <- args[1]; HR <- as.numeric(args[2]); N <- as.integer(args[3])
C1 <- as.numeric(args[4]); C2 <- as.numeric(args[5]); CAMP <- args[6]; REPS <- as.integer(args[7])
QR <- identical(if (length(args) >= 8) args[8] else "FALSE", "TRUE")
qd <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
setwd(qd)

fails <- character(0); nchk <- 0L
chk <- function(ok, what, got = "") {
  nchk <<- nchk + 1L
  if (isTRUE(ok)) cat(sprintf("[ ok ] %s\n", what))
  else { cat(sprintf("[FAIL] %s %s\n", what, got)); fails <<- c(fails, what) }
}
thr_tag <- if (abs(C1 - 0.90) < 1e-12 && abs(C2 - 0.80) < 1e-12) "" else
  sprintf("_c%03dc%03d", round(100 * C1), round(100 * C2))
path <- function(tag)
  sprintf("results/%s_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d_null%03d_nb20%s_%s%s_res_1_%d.rds",
          tag, round(100 * HR), N, round(1000 * HR), thr_tag, CAMP, if (QR) "_quickrun" else "", REPS)

cat(sprintf("GATE C -- cell %s (hr %.3f, n %d, c1 %.2f, c2 %.2f, %d reps)\n", CELL, HR, N, C1, C2, REPS))
tags <- c(consistency = "fs", dina = "dina", grf = "grf")
bs <- list()
for (e in names(tags)) {
  p <- path(tags[[e]])
  chk(file.exists(p), sprintf("%s bundle present", e), sprintf("(%s)", p))
  if (file.exists(p)) bs[[e]] <- readRDS(p)
}
if (length(bs) == 3L) {
  rs <- lapply(bs, function(b) b$results[order(b$results$sim_id), , drop = FALSE])
  ref <- rs[["consistency"]]
  for (e in names(rs))
    chk(nrow(rs[[e]]) == REPS && all(is.finite(rs[[e]]$itt_est)),
        sprintf("%s: itt_est finite on all %d rows", e, REPS))
  for (e in c("dina", "grf")) {
    chk(identical(rs[[e]]$sim_id, ref$sim_id), sprintf("%s: same sim_id vector", e))
    for (k in c("itt_est", "itt_se")) {
      a <- unname(rs[[e]][[k]]); c0 <- unname(ref[[k]])
      chk(identical(a, c0), sprintf("%s: %s identical to consistency on all %d rows", e, k, length(c0)),
          sprintf("(max |diff| %s)", format(suppressWarnings(max(abs(a - c0), na.rm = TRUE)))))
    }
    chk(isTRUE(all.equal(bs[[e]]$meta$k_treat, bs[["consistency"]]$meta$k_treat)),
        sprintf("%s: same calibrated k_treat", e))
    chk(identical(bs[[e]]$truth, bs[["consistency"]]$truth), sprintf("%s: same truth object", e))
  }
  cat(sprintf("\n  realized ITT HR: median %.4f ; Q1 %.4f ; Q3 %.4f\n",
              stats::median(ref$itt_est), stats::quantile(ref$itt_est, 0.25),
              stats::quantile(ref$itt_est, 0.75)))
  cat("  declarations per identifier:\n")
  for (e in names(rs))
    cat(sprintf("    %-11s %4d / %d (%.4f)\n", e, sum(rs[[e]]$detected == 1L), REPS,
                mean(rs[[e]]$detected == 1L)))
}
cat(sprintf("GATE C %s (%d checks, %d failed)\n",
            if (length(fails)) "FAILED" else "PASS", nchk, length(fails)))
quit(status = if (length(fails)) 1 else 0)
