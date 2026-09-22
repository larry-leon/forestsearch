#!/usr/bin/env Rscript
# nullc125 Step 1.5 -- inertness of the threshold knobs, same machine.
# Compares, per identifier, the pre-edit render (HEAD template) against the
# post-edit render with the knobs unset and with them set to 0.90 / 0.80.
# All three must be identical on every shared non-timing results column and on
# the truth object.  Exit status 1 on any difference.
qd <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
setwd(qd)
stem <- function(tag, camp)
  sprintf("results/%s_effMaxSG_fb_mr_field_m1_h066_knoise0_n500_null657_nb20_nomr_%s_quickrun_res_1_20.rds",
          tag, camp)
camps  <- c(pre = "nullc125inertpre", unset = "nullc125inertunset", explicit = "nullc125inertexpl")
tags   <- c(consistency = "fs", dina = "dina", grf = "grf")
timing <- c("fit_mr_secs", "fb_secs", "fld_H_secs", "fld_Hc_secs", "fld_H_uniform_secs")
fails  <- 0L
cat("nullc125 Step 1.5 -- threshold-knob inertness\n")
cat(sprintf("  host %s ; R %s ; forestsearch %s ; %s\n", Sys.info()[["nodename"]],
            getRversion(), utils::packageVersion("forestsearch"), format(Sys.time())))
cat("  cell null0657_n500, sim_id 1-20, MR off, effMaxSG eps 0.20, FS_S7_QUICKRUN=TRUE\n")
cat("  pre      = HEAD template (before the edit)\n  unset    = edited template, FS_S7_C1/C2 unset\n")
cat("  explicit = edited template, FS_S7_C1=0.90 FS_S7_C2=0.80\n\n")
for (e in names(tags)) {
  bs <- lapply(camps, function(cp) readRDS(stem(tags[[e]], cp)))
  rs <- lapply(bs, function(b) b$results[order(b$results$sim_id), , drop = FALSE])
  shared <- Reduce(intersect, lapply(rs, names))
  shared <- setdiff(shared, timing)
  added  <- setdiff(names(rs$unset), names(rs$pre))
  cat(sprintf("== %s ==\n", e))
  cat(sprintf("  declarations pre/unset/explicit: %s\n",
              paste(vapply(rs, function(r) sum(r$detected == 1L), integer(1)), collapse = " / ")))
  cat(sprintf("  columns added by the edit: %s\n", paste(added, collapse = ", ")))
  for (k in c("unset", "explicit")) {
    diff <- shared[!vapply(shared, function(cn)
      identical(unname(rs$pre[[cn]]), unname(rs[[k]][[cn]])), logical(1))]
    tr_ok <- identical(bs$pre$truth, bs[[k]]$truth)
    cat(sprintf("  pre vs %-8s: identical on %d of %d shared non-timing columns%s ; truth identical %s\n",
                k, length(shared) - length(diff), length(shared),
                if (length(diff)) paste0(" -- DIFFER: ", paste(diff, collapse = ", ")) else "",
                tr_ok))
    if (length(diff) || !tr_ok) fails <- fails + 1L
  }
  # the added columns agree between the two post-edit renders
  diff2 <- added[!vapply(added, function(cn)
    identical(unname(rs$unset[[cn]]), unname(rs$explicit[[cn]])), logical(1))]
  cat(sprintf("  unset vs explicit on the added columns: %s\n",
              if (length(diff2)) paste("DIFFER:", paste(diff2, collapse = ", ")) else "identical"))
  if (length(diff2)) fails <- fails + 1L
  m <- bs$explicit$meta
  cat(sprintf("  explicit meta c1/c2/pstar/dmin_grf/dina_m_diff: %s / %s / %s / %s / %.6f\n",
              m$c1, m$c2, m$pstar, m$dmin_grf, m$dina_m_diff))
  cat(sprintf("  resolved c1 / c2 (explicit): %s ; NA on %d of 20\n",
              paste(unique(rs$explicit$c1_resolved), unique(rs$explicit$c2_resolved), sep = " / "),
              sum(is.na(rs$explicit$c1_resolved) | is.na(rs$explicit$c2_resolved))))
  cat(sprintf("  itt_est finite on %d of 20\n\n", sum(is.finite(rs$explicit$itt_est))))
}
cat(sprintf("INERTNESS %s (%d failing comparisons)\n", if (fails) "FAILED" else "PASS", fails))
quit(status = if (fails) 1 else 0)
