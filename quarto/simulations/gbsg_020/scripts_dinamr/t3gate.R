# Gate T3 comparison (TASK_partB_enabling_2026-09-12, Part 1), stop-on-failure.
#
# HALF 1 -- unset.  t3pre (template before the edit) against t3post (after, with
#   FS_S7_MR unset), ON THIS MACHINE at the standing identity cell: the same
#   columns in the same order, every non-timing column identical(), truth
#   identical().  The change is required to be default-inert.
# HALF 2 -- MR off.  t3post against t3mroff (FS_S7_MR=FALSE), same seeds: the
#   identification and classification columns present, populated on detected
#   rows, and identical() to the MR-on run.  MR cannot change the identified
#   subgroup (R/forestsearch_main.R:1041-1044); this half tests that, and tests
#   that the recorder's MR-off routes reproduce what the MR object recorded.
#   n_family is EXEMPT and must be NA: it is MR's own fitted family size and no
#   result field carries it (template comment in record_replicate(); report).
# SUPPLEMENTARY (not the gate as specified; run if the bundles exist): the same
#   half-2 comparison for DINA and GRF, whose label route differs.
#
# usage (from scripts_dinamr/): Rscript t3gate.R
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")
RES <- file.path(QMD_DIR, "results")
f <- function(eng, tag, nomr = FALSE) file.path(RES, sprintf(
  "%s_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_nb20%s_%s_res_1_5.rds",
  eng, if (nomr) "_nomr" else "", tag))
TIMING <- c("fit_mr_secs", "fb_secs", "fld_H_secs", "fld_Hc_secs",
            "fld_H_uniform_secs", "fld_recov_secs", "fld_joint_secs")
# Identification and classification: what Part B records.  Everything else in
# the record is an MR (nv_*, mr_*, fld_*, p_hat_*), FB (fb_*) or timing column.
ID <- c("sim_id", "detected", "status", "n_harm", "n_true", "label", "sg_def",
        "covs", "err_msg", "n_sel", "n_cons_qual", "band_n", "admitted_n",
        "sens", "spec", "ppv", "npv", "betaHhat_H", "betaHhat_Hc",
        "betaHhat_status", "nH_eval", "nHc_eval",
        "or_H_est", "or_H_lo", "or_H_hi", "or_H_se",
        "or_Hc_est", "or_Hc_lo", "or_Hc_hi", "or_Hc_se")
POPULATED <- c("n_sel", "label", "detected", "sens", "spec", "ppv", "npv")
fail <- 0L; pass <- 0L
say <- function(ok, msg) { if (ok) pass <<- pass + 1L else fail <<- fail + 1L
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", msg)) }

cat("=== HALF 1: unset, pre-change vs post-change (consistency) ===\n")
pre <- readRDS(f("fs", "t3pre")); post <- readRDS(f("fs", "t3post"))
rp <- pre$results; rq <- post$results
say(nrow(rp) == 5L && nrow(rq) == 5L, sprintf("5 rows each (pre %d, post %d)", nrow(rp), nrow(rq)))
say(identical(names(rp), names(rq)), sprintf("identical column names and order (%d columns)", ncol(rq)))
shared <- setdiff(intersect(names(rp), names(rq)), TIMING)
bad <- shared[!vapply(shared, function(cc) identical(rp[[cc]], rq[[cc]]), logical(1))]
say(!length(bad), sprintf("every non-timing column identical() across %d columns (mismatched: %s)",
                          length(shared), if (length(bad)) paste(bad, collapse = ", ") else "<none>"))
say(identical(pre$truth, post$truth), "truth identical()")
say(is.null(pre$meta$mr_inference) && isTRUE(post$meta$mr_inference),
    "meta$mr_inference absent pre-change, TRUE post-change (record only)")

half2 <- function(eng, on_tag, off_tag, label) {
  cat(sprintf("\n=== %s: MR on (unset) vs MR off, %s ===\n", label, eng))
  a <- readRDS(f(eng, on_tag)); b <- readRDS(f(eng, off_tag, nomr = TRUE))
  ra <- a$results; rb <- b$results
  say(nrow(ra) == 5L && nrow(rb) == 5L, sprintf("5 rows each (on %d, off %d)", nrow(ra), nrow(rb)))
  say(identical(names(ra), names(rb)), "identical column set")
  say(isTRUE(a$meta$mr_inference) && identical(b$meta$mr_inference, FALSE),
      "meta$mr_inference TRUE (on) / FALSE (off)")
  d <- rb$detected %in% 1L
  cat(sprintf("detected: on %d/5, off %d/5\n", sum(ra$detected %in% 1L), sum(d)))
  miss <- setdiff(POPULATED, names(rb))
  say(!length(miss), sprintf("columns present (missing: %s)", if (length(miss)) paste(miss, collapse = ", ") else "<none>"))
  unpop <- POPULATED[vapply(POPULATED, function(cc) any(is.na(rb[[cc]][d])), logical(1))]
  say(any(d) && !length(unpop),
      sprintf("populated on every detected MR-off row (NA in: %s)", if (length(unpop)) paste(unpop, collapse = ", ") else "<none>"))
  if (eng == "fs") {
    unpop2 <- c("n_cons_qual", "band_n")[vapply(c("n_cons_qual", "band_n"), function(cc) any(is.na(rb[[cc]][d])), logical(1))]
    say(!length(unpop2), sprintf("n_cons_qual, band_n populated on detected rows, consistency (NA in: %s)",
                                 if (length(unpop2)) paste(unpop2, collapse = ", ") else "<none>"))
  }
  idc <- intersect(ID, names(ra))
  bad <- idc[!vapply(idc, function(cc) identical(ra[[cc]], rb[[cc]]), logical(1))]
  say(!length(bad), sprintf("identification + classification columns identical() on/off across %d columns (mismatched: %s)",
                            length(idc), if (length(bad)) paste(bad, collapse = ", ") else "<none>"))
  for (cc in bad) { cat("  ", cc, "\n"); print(data.frame(sim_id = ra$sim_id, on = ra[[cc]], off = rb[[cc]])) }
  say(identical(a$truth, b$truth), "truth identical()")
  say(all(is.na(rb$n_family)), "n_family NA on every MR-off row (EXEMPT: MR-only quantity, documented)")
  cat(sprintf("  n_family with MR on: %s\n", paste(ra$n_family, collapse = ", ")))
  mrcols <- grep("^(nv_|mr_H|mr_Hc|fld_|p_hat_)", names(rb), value = TRUE)
  say(all(rb$mr_ok == 0L) && all(vapply(mrcols, function(cc) all(is.na(rb[[cc]])), logical(1))),
      sprintf("MR off: mr_ok 0 and all %d MR-derived columns NA (by design)", length(mrcols)))
  other <- setdiff(names(ra), c(idc, TIMING, mrcols, "n_family", "mr_ok"))
  diff_other <- other[!vapply(other, function(cc) identical(ra[[cc]], rb[[cc]]), logical(1))]
  cat(sprintf("  columns outside the ID/MR/timing sets that differ (informational): %s\n",
              if (length(diff_other)) paste(diff_other, collapse = ", ") else "<none>"))
  cat("  selected rules (off):\n"); print(rb[, c("sim_id", "label", "sg_def", "n_sel", "n_cons_qual", "band_n")])
}
half2("fs", "t3post", "t3mroff", "HALF 2")
gate_fail <- fail
for (e in c("dina", "grf")) {
  if (file.exists(f(e, paste0("t3", e))) && file.exists(f(e, paste0("t3", e, "off"), nomr = TRUE)))
    half2(e, paste0("t3", e), paste0("t3", e, "off"), "SUPPLEMENTARY")
}
cat(sprintf("\nGATE T3 (halves 1 and 2): %d failures -- %s\n", gate_fail,
            if (gate_fail) "FAIL: REVERT AND STOP" else "PASS"))
cat(sprintf("supplementary DINA/GRF: %d failures\n", fail - gate_fail))
if (fail) quit(status = 1L)
