# Section 1 of CC_BRIEF_sim_integrity.md: bit-identical comparison of the
# pre-edit and post-edit per-replicate records, matched on sim_id, tolerance 0.
# Usage: Rscript compare_rds.R <pre.rds> <post.rds>

args <- commandArgs(TRUE)
pre  <- readRDS(args[1])
post <- readRDS(args[2])

rp <- pre$results
rq <- post$results

cat("== meta ==\n")
cat("pre  meta: ", paste(sprintf("%s=%s", names(pre$meta),
      vapply(pre$meta, function(x) paste(x, collapse=","), "")), collapse="; "), "\n")
cat("post meta: ", paste(sprintf("%s=%s", names(post$meta),
      vapply(post$meta, function(x) paste(x, collapse=","), "")), collapse="; "), "\n\n")

cat("== truth targets ==\n")
cat("identical(pre$truth, post$truth): ", identical(pre$truth, post$truth), "\n")
if (!identical(pre$truth, post$truth)) print(all.equal(pre$truth, post$truth))
cat("\n")

## ---- name mapping: t2_secs first (specific), then generic prefixes ----------
map_names <- function(nm) {
  nm[nm == "t2_secs"] <- "fit_mr_secs"
  nm <- sub("^t1_", "fb_", nm)
  nm <- sub("^t2_", "mr_", nm)
  nm
}
np_raw <- names(rp)
np <- map_names(np_raw)
nq <- names(rq)

cat("== column sets after mapping (t1_->fb_, t2_->mr_, t2_secs->fit_mr_secs) ==\n")
cat("pre  cols: ", length(np), "   post cols: ", length(nq), "\n")
only_pre  <- setdiff(np, nq)
only_post <- setdiff(nq, np)
cat("in mapped-PRE but not POST : ",
    if (length(only_pre))  paste(only_pre,  collapse=", ") else "(none)", "\n")
cat("in POST but not mapped-PRE : ",
    if (length(only_post)) paste(only_post, collapse=", ") else "(none)", "\n")
cat("column ORDER identical      : ", identical(np, nq), "\n\n")

## ---- match on sim_id --------------------------------------------------------
names(rp) <- np
cat("== sim_id alignment ==\n")
cat("pre  sim_id: ", paste(rp$sim_id, collapse=", "), "\n")
cat("post sim_id: ", paste(rq$sim_id, collapse=", "), "\n")
common <- intersect(rp$sim_id, rq$sim_id)
cat("common: ", length(common), "\n\n")
rp <- rp[match(common, rp$sim_id), , drop = FALSE]
rq <- rq[match(common, rq$sim_id), , drop = FALSE]

## ---- per-column comparison at tolerance 0 -----------------------------------
WALL <- c("fb_secs", "fit_mr_secs")
cols <- intersect(np, nq)
res <- data.frame(column = character(0), class_pre = character(0),
                  class_post = character(0), status = character(0),
                  detail = character(0), stringsAsFactors = FALSE)

for (cl in cols) {
  a <- rp[[cl]]; b <- rq[[cl]]
  exempt <- cl %in% WALL
  if (!identical(class(a), class(b))) {
    st <- "CLASS-DIFF"; det <- sprintf("%s vs %s", paste(class(a), collapse="/"),
                                       paste(class(b), collapse="/"))
  } else if (identical(a, b)) {
    st <- if (exempt) "IDENTICAL (wall-clock, exempt)" else "IDENTICAL"; det <- ""
  } else {
    d <- which(!( (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b) ))
    det <- paste(sprintf("sim_id %s: %s vs %s", common[d],
                         format(a[d], digits = 12), format(b[d], digits = 12)),
                 collapse = " | ")
    st <- if (exempt) "DIFFERS (wall-clock, exempt)" else "DIFFERS"
  }
  res <- rbind(res, data.frame(column = cl, class_pre = paste(class(a), collapse="/"),
                               class_post = paste(class(b), collapse="/"),
                               status = st, detail = det,
                               stringsAsFactors = FALSE))
}

cat("== per-column comparison, tolerance 0 ==\n")
cat(sprintf("columns compared        : %d\n", nrow(res)))
cat(sprintf("identical               : %d\n", sum(grepl("^IDENTICAL", res$status))))
cat(sprintf("differing (wall-clock)  : %d\n", sum(res$status == "DIFFERS (wall-clock, exempt)")))
cat(sprintf("differing (SUBSTANTIVE) : %d\n", sum(res$status == "DIFFERS")))
cat(sprintf("class mismatches        : %d\n\n", sum(res$status == "CLASS-DIFF")))

bad <- res[res$status %in% c("DIFFERS", "CLASS-DIFF"), , drop = FALSE]
if (nrow(bad)) {
  cat("!! SUBSTANTIVE DIFFERENCES !!\n")
  for (i in seq_len(nrow(bad)))
    cat(sprintf("  %-14s [%s] %s\n", bad$column[i], bad$status[i], bad$detail[i]))
} else {
  cat("No substantive differences: every non-wall-clock column is bit-identical.\n")
}

cat("\n== wall-clock columns (exempt) ==\n")
for (cl in intersect(WALL, cols))
  cat(sprintf("  %-12s pre: %s\n  %-12s post: %s\n", cl,
              paste(round(rp[[cl]], 2), collapse=", "), "",
              paste(round(rq[[cl]], 2), collapse=", ")))

cat("\n== full status table ==\n")
print(res[, c("column", "status")], row.names = FALSE)
