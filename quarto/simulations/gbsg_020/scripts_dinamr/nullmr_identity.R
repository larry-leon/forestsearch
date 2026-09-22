#!/usr/bin/env Rscript
# nullmr Step 3.2 -- IDENTITY GATE, same machine: "MR cannot change
# identification", machine-checked.  Per identifier, the MR-on smoke of
# null0657_n500 (sim_id 1-20) against nullc125's committed MR-off render of the
# same replicates at the same (unset = 0.90 / 0.80) screen:
#   results/*_nomr_nullc125inertunset_quickrun_res_1_20.rds
# Identical (identical()) on every shared results column except the classes
# below, which are fixed here before either smoke output is read, and on the
# truth object.  Exit status 1 on any difference.
#
# Excluded, by class:
#   timing       fit_mr_secs, fb_secs, fld_H_secs, fld_Hc_secs, fld_H_uniform_secs
#                (nullthr_inertness.R:14)
#   MR / field / IJ products: every mr_* and fld_* column (mr_ok and
#                mr_harm_flag included), ij_source, n_family, p_hat_*
#   the naive (unadjusted) block nv_*: with MR on it is MR's own g$naive
#                (template :1455, :1466-1469); with MR off it is a template
#                Cox refit (:1403-1413).  Two code paths for one quantity, so
#                it is an MR-path product here; its agreement is PRINTED
#                (max |diff|), not gated.
#   MR-off-branch-only fields: p_sel, p_max_qual (:1391-1395) and nv_H_lo1s
#                (:1408-1409), NA by construction with MR on.
# Everything else -- status, detected, label, sg_def, covs, n_sel, n_harm,
# betaHhat_*, classification, band / qualifying counts, family counts, maxT,
# resolved thresholds, itt_*, oracle, admitted_n, fb_* -- is gated.
#
# usage: nullmr_identity.R <smoke-campaign-tag>
args <- commandArgs(trailingOnly = TRUE)
CAMP <- if (length(args)) args[1] else "nullmrsmoke"
qd <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
setwd(qd)
mr_on  <- function(tag) sprintf("results/%s_effMaxSG_fb_mr_field_m1_h066_knoise0_n500_null657_nb20_%s_quickrun_res_1_20.rds", tag, CAMP)
mr_off <- function(tag) sprintf("results/%s_effMaxSG_fb_mr_field_m1_h066_knoise0_n500_null657_nb20_nomr_nullc125inertunset_quickrun_res_1_20.rds", tag)
tags <- c(consistency = "fs", dina = "dina", grf = "grf")
timing <- c("fit_mr_secs", "fb_secs", "fld_H_secs", "fld_Hc_secs", "fld_H_uniform_secs")
excl_class <- function(nm) {
  if (nm %in% timing) return("timing")
  if (grepl("^(mr_|fld_|p_hat_)", nm) || nm %in% c("ij_source", "n_family")) return("MR/field/IJ")
  if (nm %in% c("p_sel", "p_max_qual", "nv_H_lo1s")) return("MR-off-only")
  if (grepl("^nv_", nm)) return("naive (MR-path)")
  NA_character_
}
fails <- 0L
cat("nullmr Step 3.2 -- identity gate: MR on vs nullc125's MR-off render, same replicates\n")
cat(sprintf("  host %s ; R %s ; forestsearch %s ; %s\n", Sys.info()[["nodename"]],
            getRversion(), utils::packageVersion("forestsearch"), format(Sys.time())))
cat(sprintf("  cell null0657_n500, sim_id 1-20 ; MR-on campaign %s ; MR-off campaign nullc125inertunset\n\n", CAMP))
for (e in names(tags)) {
  a <- readRDS(mr_on(tags[[e]])); b <- readRDS(mr_off(tags[[e]]))
  ra <- a$results[order(a$results$sim_id), , drop = FALSE]
  rb <- b$results[order(b$results$sim_id), , drop = FALSE]
  shared <- intersect(names(ra), names(rb))
  cls <- vapply(shared, excl_class, character(1))
  gated <- shared[is.na(cls)]
  diff <- gated[!vapply(gated, function(cn) identical(unname(ra[[cn]]), unname(rb[[cn]])), logical(1))]
  tr_ok <- identical(a$truth, b$truth)
  cat(sprintf("== %s ==\n", e))
  cat(sprintf("  meta mr_inference on/off: %s / %s ; rows %d / %d\n",
              a$meta$mr_inference, b$meta$mr_inference, nrow(ra), nrow(rb)))
  cat(sprintf("  declarations on/off: %d / %d ; mr_ok on declaring (MR on): %d\n",
              sum(ra$detected == 1L), sum(rb$detected == 1L), sum(ra$mr_ok == 1L & ra$detected == 1L)))
  cat(sprintf("  shared columns %d ; excluded %d (%s) ; gated %d\n", length(shared), sum(!is.na(cls)),
              paste(sprintf("%s %d", names(table(cls)), as.integer(table(cls))), collapse = ", "),
              length(gated)))
  cat(sprintf("  identical on %d of %d gated columns%s ; truth identical %s\n",
              length(gated) - length(diff), length(gated),
              if (length(diff)) paste0(" -- DIFFER: ", paste(diff, collapse = ", ")) else "", tr_ok))
  for (cn in c("label", "n_sel", "sg_def")) if (cn %in% diff) {
    w <- which(!mapply(identical, ra[[cn]], rb[[cn]]))
    cat(sprintf("    %s differs at sim_id %s\n", cn, paste(ra$sim_id[w], collapse = ",")))
  }
  nvc <- shared[cls %in% "naive (MR-path)"]
  for (cn in nvc) {
    x <- unname(ra[[cn]]); y <- unname(rb[[cn]])
    both <- is.finite(x) & is.finite(y)
    cat(sprintf("    [printed] %-10s finite on/off %2d / %2d ; same NA pattern %s ; max |diff| %s\n",
                cn, sum(is.finite(x)), sum(is.finite(y)), identical(is.finite(x), is.finite(y)),
                if (any(both)) format(max(abs(x[both] - y[both])), digits = 3) else "-"))
  }
  cat("\n")
  if (length(diff) || !tr_ok) fails <- fails + 1L
}
cat(sprintf("IDENTITY GATE %s (%d failing identifiers)\n", if (fails) "FAILED" else "PASS", fails))
quit(status = if (fails) 1 else 0)
