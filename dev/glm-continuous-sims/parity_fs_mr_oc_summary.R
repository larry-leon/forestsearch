# =============================================================================
# parity_fs_mr_oc_summary.R -- acceptance test for the fs_mr_oc_summary()
#                              extraction, against committed payloads
# =============================================================================
#
# fs_mr_oc_summary() was transplanted out of the `.build_block()` machinery
# duplicated across the MD40 batch documents.  This script re-implements that
# machinery verbatim, runs both against the same committed payloads, and
# requires all thirteen numeric columns to agree EXACTLY.
#
# It cannot live in tests/testthat/ because it needs the payloads, which are
# large and are not fixtures.
#
# Expect: max |difference| of 0.000e+00 on every payload.
#
# NOTE ON THE ORACLE ROWS.  Payloads written before commit 7c517531 carry
# swapped oracle CI bounds, so their oracle coverage reads 0 and their oracle
# width is negative.  Both implementations reproduce that faithfully, and must:
# this script tests that the extraction is EQUIVALENT, not that the data are
# correct.  The oracle rows become meaningful only once the payloads are
# regenerated.
#
# Usage: Rscript parity_fs_mr_oc_summary.R <payload.rds> [<payload.rds> ...]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  stop("usage: Rscript parity_fs_mr_oc_summary.R <payload.rds> ...", call. = FALSE)
}
# Load the WORKING TREE, not an installed build.  This is an acceptance test
# for the source: an installed forestsearch lacking the function fails
# confusingly, and one carrying an OLDER copy of it would silently test the
# wrong code and appear to pass.  devtools::load_all() is correct here because
# this script is single-process; parallel harnesses still require
# devtools::install(), since workers only see the installed package.
if (file.exists("DESCRIPTION") && requireNamespace("devtools", quietly = TRUE)) {
  suppressMessages(devtools::load_all(".", quiet = TRUE))
} else {
  suppressPackageStartupMessages(library(forestsearch))
}
if (!exists("fs_mr_oc_summary")) {
  stop("fs_mr_oc_summary() not found -- run from the package root, or ",
       "install a build that contains it.", call. = FALSE)
}


# -- The batch documents' inline machinery, transcribed verbatim --------------

doc_build <- function(pl) {
  results <- pl$results
  truth   <- pl$truth
  meta    <- pl$meta
  adverse_outcome <- meta$adverse_outcome
  nb_boots <- meta$nb_boots
  run_fb <- !is.null(nb_boots) && nb_boots >= 1L

  rd <- results[results$status %in% "DETECTED", , drop = FALSE]
  .orient   <- if (isTRUE(adverse_outcome)) 1 else -1
  rd$bH_or  <- .orient * rd$betaHhat_H
  rd$bHc_or <- .orient * rd$betaHhat_Hc

  .mfin  <- function(x) { x <- x[is.finite(x)]; if (length(x)) mean(x) else NA_real_ }
  .sdfin <- function(x) { x <- x[is.finite(x)]; if (length(x) > 1L) sd(x) else NA_real_ }
  .cover <- function(target, lo, hi) {
    target <- rep_len(target, length(lo))
    ok <- is.finite(target) & is.finite(lo) & is.finite(hi)
    if (!any(ok)) return(NA_real_)
    mean(target[ok] >= lo[ok] & target[ok] <= hi[ok])
  }

  est_keys <- c(oracle = "or", naive = "nv", FB = "fb", MR = "mr")
  if (!run_fb) est_keys <- est_keys[names(est_keys) != "FB"]

  .blk <- list(
    H  = list(beta = "bH_or",  struct = .orient * truth$effect_Q,
              marg = .orient * truth$marg_H,  or_col = "or_H_est"),
    Hc = list(beta = "bHc_or", struct = .orient * truth$effect_Qc,
              marg = .orient * truth$marg_Hc, or_col = "or_Hc_est")
  )
  .se_col_blk <- function(k, sfx) {
    if (k == "mr") paste0("mr_", sfx, "_se_ij") else paste0(k, "_", sfx, "_se")
  }

  build <- function(sfx) {
    b <- .blk[[sfx]]
    beta_ref <- rd[[b$beta]]
    or_ref   <- rd[[b$or_col]]
    do.call(rbind, lapply(names(est_keys), function(lab) {
      k  <- est_keys[[lab]]
      e  <- rd[[paste0(k, "_", sfx, "_est")]]
      if (is.null(e)) return(NULL)
      lo <- rd[[paste0(k, "_", sfx, "_lo")]]
      hi <- rd[[paste0(k, "_", sfx, "_hi")]]
      se <- rd[[.se_col_blk(k, sfx)]]
      if (is.null(se)) se <- rep(NA_real_, nrow(rd))
      data.frame(
        Block = sfx, Estimator = lab, n = sum(is.finite(e)),
        Avg = .mfin(e), SD_emp = .sdfin(e), SE_hat = .mfin(se),
        b_beta = .mfin(e - beta_ref), b_or = .mfin(e - or_ref),
        b_str  = .mfin(e - b$struct), b_marg = .mfin(e - b$marg),
        C_beta = .cover(beta_ref, lo, hi), C_or = .cover(or_ref, lo, hi),
        C_str  = .cover(b$struct, lo, hi), C_marg = .cover(b$marg, lo, hi),
        Width  = .mfin(hi - lo), stringsAsFactors = FALSE
      )
    }))
  }
  rbind(build("H"), build("Hc"))
}


# -- Compare ------------------------------------------------------------------

COLS <- c(n = "n", Avg = "avg", SD_emp = "sd_emp", SE_hat = "se_hat",
          b_beta = "bias_beta", b_or = "bias_oracle", b_str = "bias_struct",
          b_marg = "bias_marg", C_beta = "cov_beta", C_or = "cov_oracle",
          C_str = "cov_struct", C_marg = "cov_marg", Width = "width")

ok_all <- TRUE

for (f in args) {
  pl  <- readRDS(f)
  ref <- doc_build(pl)
  new <- fs_mr_oc_summary(pl)$estimation

  cat(sprintf("payload: %s  (n = %s, sims = %s)\n",
              basename(f), pl$meta$n_sample, pl$meta$n_sims))

  if (nrow(ref) != nrow(new) ||
      !identical(ref$Block, new$block) ||
      !identical(ref$Estimator, new$estimator)) {
    cat("  *** row structure differs ***\n\n")
    ok_all <- FALSE
    next
  }

  worst <- 0
  for (a in names(COLS)) {
    d <- suppressWarnings(max(abs(ref[[a]] - new[[COLS[[a]]]]), na.rm = TRUE))
    if (!is.finite(d)) d <- 0
    na_ok <- identical(is.na(ref[[a]]), is.na(new[[COLS[[a]]]]))
    worst <- max(worst, d)
    if (d > 0 || !na_ok) {
      cat(sprintf("  MISMATCH %-8s max|diff| %.3e   NA-pattern ok %s\n",
                  a, d, na_ok))
      ok_all <- FALSE
    }
  }
  cat(sprintf("  %d rows, 13 columns, max |difference| %.3e   %s\n\n",
              nrow(ref), worst, if (worst == 0) "EXACT" else "NOT EXACT"))
  if (worst != 0) ok_all <- FALSE
}

cat(if (ok_all) "PARITY EXACT ON ALL PAYLOADS\n" else "*** PARITY FAILURE ***\n")
if (!ok_all) quit(status = 1L)
