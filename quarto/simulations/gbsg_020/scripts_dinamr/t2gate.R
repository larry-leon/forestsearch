# Gate T2 comparison (TASK_grfmr_campaign_2026-09-11, Part T2), stop-on-failure.
#
# Pre-change against post-change ON THIS MACHINE at the standing identity cell,
# FS_S7_METHOD unset (the consistency path).  The requirement is that the
# recorder change is ADD-ONLY: every column that existed before must be
# identical(), and `truth` must be identical().  The new column `admitted_n` is
# expected to be present post-change and absent pre-change, and is the ONLY
# permitted difference in the column set.
#
# Timing columns are excluded by name: they measure the machine, not the record.
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")
RES <- file.path(QMD_DIR, "results")
f <- function(tag) file.path(RES, sprintf(
  "fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_nb20_%s_res_1_5.rds", tag))
pre <- readRDS(f("t2pre")); post <- readRDS(f("t2post"))
rp <- pre$results; rq <- post$results
TIMING <- c("fit_mr_secs", "fb_secs", "fld_H_secs", "fld_Hc_secs",
            "fld_recov_secs", "fld_joint_secs")
fail <- 0L; pass <- 0L
say <- function(ok, msg) { if (ok) pass <<- pass + 1L else fail <<- fail + 1L
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", msg)) }

say(nrow(rp) == 5L && nrow(rq) == 5L,
    sprintf("5 rows each (pre %d, post %d)", nrow(rp), nrow(rq)))
new_cols <- setdiff(names(rq), names(rp))
gone     <- setdiff(names(rp), names(rq))
say(identical(new_cols, "admitted_n"),
    sprintf("the ONLY new column is admitted_n (new: %s)",
            if (length(new_cols)) paste(new_cols, collapse = ", ") else "<none>"))
say(!length(gone),
    sprintf("no column removed (removed: %s)",
            if (length(gone)) paste(gone, collapse = ", ") else "<none>"))
# Column ORDER of the pre-existing columns must be preserved as a block-relative
# order; admitted_n is inserted, so compare the subsequence.
say(identical(setdiff(names(rq), "admitted_n"), names(rp)),
    "pre-existing columns keep their names AND their order")

shared <- setdiff(names(rp), TIMING)
bad <- shared[!vapply(shared, function(cc) identical(rp[[cc]], rq[[cc]]), logical(1))]
say(!length(bad),
    sprintf("every non-timing column identical() across %d columns (mismatched: %s)",
            length(shared), if (length(bad)) paste(bad, collapse = ", ") else "<none>"))
say(identical(pre$truth, post$truth), "truth identical()")
say(all(is.na(rq$admitted_n)),
    "admitted_n is NA on every row of a NON-GRF (consistency) run")
tim <- intersect(TIMING, names(rp))
cat(sprintf("\ntiming columns excluded by name (%d present): %s\n",
            length(tim), paste(tim, collapse = ", ")))
cat(sprintf("\nGATE T2: %d passes, %d failures -- %s\n", pass, fail,
            if (fail) "FAIL: REVERT AND STOP" else "PASS"))
if (fail) quit(status = 1L)
