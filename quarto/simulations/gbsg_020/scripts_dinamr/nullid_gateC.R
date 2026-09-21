#!/usr/bin/env Rscript
# GATE C -- nullid, one per cell: the three identifiers ran on IDENTICAL draws.
#
# Under the structural null flag_harm == 0 on every subject, so the template's
# oracle complement refit or_Hc_* is a Cox fit of treatment alone on the WHOLE
# trial.  It is computed from the simulated data only and never touches the
# engine, so it is an engine-independent fingerprint of the draw.  Equality of
# or_Hc_est / or_Hc_se across the three bundles, row for row, is therefore the
# same-draws check.  n_true is 0 everywhere here and cannot discriminate, which
# is why it is not used (it is the check the alt-design campaigns use).
#
# usage: nullid_gateC.R <cell> <hr> <n>
args <- commandArgs(trailingOnly = TRUE)
CELL <- args[1]; HR <- as.numeric(args[2]); N <- as.integer(args[3])
qd <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
setwd(qd)

fails <- character(0); nchk <- 0L
chk <- function(ok, what, got = "") {
  nchk <<- nchk + 1L
  if (isTRUE(ok)) cat(sprintf("[ ok ] %s\n", what))
  else { cat(sprintf("[FAIL] %s %s\n", what, got)); fails <<- c(fails, what) }
}
path <- function(tag)
  sprintf("results/%s_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d_null%03d_nb20_nomr_nullid_res_1_2000.rds",
          tag, round(100 * HR), N, round(1000 * HR))

cat(sprintf("GATE C -- cell %s (hr %.3f, n %d)\n", CELL, HR, N))
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
  for (e in c("dina", "grf")) {
    chk(identical(rs[[e]]$sim_id, ref$sim_id), sprintf("%s: same sim_id vector", e))
    chk(identical(rs[[e]]$n_true, ref$n_true), sprintf("%s: same n_true vector", e))
    for (k in c("or_Hc_est", "or_Hc_se")) {
      a <- rs[[e]][[k]]; c0 <- ref[[k]]
      ok <- identical(a, c0) ||
            isTRUE(all(abs(a - c0) <= 1e-12 * pmax(1, abs(c0)), na.rm = TRUE) &&
                   identical(is.na(a), is.na(c0)))
      chk(ok, sprintf("%s: same draws via %s", e, k),
          sprintf("(max |diff| %s)", format(suppressWarnings(max(abs(a - c0), na.rm = TRUE)))))
    }
  }
  # the design point is the same object in all three
  for (e in c("dina", "grf")) {
    chk(isTRUE(all.equal(bs[[e]]$meta$k_treat, bs[["consistency"]]$meta$k_treat)),
        sprintf("%s: same calibrated k_treat", e))
    chk(identical(bs[[e]]$truth, bs[["consistency"]]$truth), sprintf("%s: same truth object", e))
  }
  cat("\n  declarations per identifier:\n")
  for (e in names(rs))
    cat(sprintf("    %-11s %4d / 2000 (%.4f)\n", e, sum(rs[[e]]$detected == 1L),
                mean(rs[[e]]$detected == 1L)))
}
cat(sprintf("GATE C %s (%d checks, %d failed)\n",
            if (length(fails)) "FAILED" else "PASS", nchk, length(fails)))
quit(status = if (length(fails)) 1 else 0)
