# Regenerates quarto/simulations/gbsg_020/STATUS_p12x20.md (task §8) from the
# directory and git -- never from chat records.  No pin claim: the HEAD recorded
# is informational; the pin is claimed once, at the merge, in a separate task.
# usage: Rscript status_p12x20.R   (run by run_p12x20.sh at halt and at closeout)
args_all <- commandArgs(FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])))
qmd  <- normalizePath(file.path(here, ".."))
repo <- system2("git", c("-C", qmd, "rev-parse", "--show-toplevel"), stdout = TRUE)
rq   <- sub(paste0("^", repo, "/"), "", qmd)
git  <- function(...) suppressWarnings(system2("git", c("-C", repo, ...), stdout = TRUE, stderr = FALSE))
TAG <- "p12x20"; DATE <- "2026-09-12"; BASE <- "0ab5d1c560cfd0c536efd914f3a6266aa4511bab"
paydir <- paste0(TAG, "_", DATE); pay <- file.path(qmd, paydir)
cells <- data.frame(cell = paste0("A", 1:9), hr = c(1.50,1.50,1.50,1.75,1.75,1.75,1.00,1.00,1.00),
                    n = rep(c(500L, 1000L, 1500L), 3), stringsAsFactors = FALSE)
hbf <- file.path(qmd, "LOG_p12x20_progress.txt")
hb  <- if (file.exists(hbf)) readLines(hbf) else character(0)
haltf <- file.path(qmd, "HALT_p12x20.md")
tracked <- function(p) length(git("ls-files", "--error-unmatch", "--", p)) > 0L

L <- c("# STATUS — p12x20 (Part A: FS effMaxSG eps 0.20, the nine 12.4% cells)", "",
       sprintf("- Generated (UTC): %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
       sprintf("- Branch: `%s`; HEAD at generation: `%s` (informational; no pin claim on this file).",
               git("rev-parse", "--abbrev-ref", "HEAD"), git("rev-parse", "--short", "HEAD")),
       sprintf("- Base SHA: `%s`.", BASE),
       "- Task: `dev/tasks/TASK_p12x20_partA_2026-09-12_v2.md`; addendum `dev/tasks/ADDENDUM_p12x20_stage2_unattended_2026-09-12_v2.md`.",
       sprintf("- Stage 0 record: `%s/REPORT_p12x20_stage0_2026-09-12.md`.", rq),
       sprintf("- Gate 1 record: `%s/REPORT_p12x20_gate1_2026-09-12.md`.", rq),
       sprintf("- Runner: `%s/scripts_p12x20/run_p12x20.sh`; payload directory `%s/%s/`.", rq, rq, paydir),
       sprintf("- HALT file: %s.", if (file.exists(haltf)) sprintf("PRESENT (`%s/HALT_p12x20.md`)", rq) else "absent"),
       sprintf("- Open items: %s.", if (file.exists(file.path(qmd, "OPEN_ITEMS_p12x20.md"))) sprintf("`%s/OPEN_ITEMS_p12x20.md`", rq) else "none"),
       "")

s0 <- file.path(qmd, "REPORT_p12x20_stage0_2026-09-12.md")
if (file.exists(s0)) {
  x <- readLines(s0); i0 <- grep("^## 0a", x); i1 <- grep("^## 0g", x)
  L <- c(L, "## Stage 0 install record (verbatim from the Stage 0 report, 0a-0f)", "",
         if (length(i0) && length(i1)) sub("^## ", "### ", x[i0:(i1 - 1)]) else "- section markers not found", "")
}

L <- c(L, "## Cells", "",
       "| cell | HR | n | status | cell wall (s) | replicates | gate counts | cell commit | gate record |",
       "|---|---|---|---|---|---|---|---|---|")
tot <- c(run = 0L, passed = 0L, failed = 0L); inv <- character(0)
for (i in seq_len(nrow(cells))) {
  ce <- cells$cell[i]; h3 <- sprintf("%03d", round(100 * cells$hr[i]))
  stem <- sprintf("fs_effMaxSG_fb_mr_field_m1_h%s_knoise0_n%d_nb20_%s", h3, cells$n[i], TAG)
  comb <- file.path(rq, paydir, paste0(stem, "_combined_1_2000.rds"))
  gname <- sprintf("GATE2_%s_%s_%s.txt", TAG, ce, DATE); grec <- file.path(pay, gname)
  rname <- sprintf("REPORT_%s_%s_%s.md", TAG, ce, DATE); rep <- file.path(pay, rname)
  done <- tracked(comb) && tracked(file.path(rq, paydir, gname))
  hbc <- grep(sprintf("\t%s\t", ce), hb, value = TRUE)
  last <- if (length(hbc)) strsplit(tail(hbc, 1), "\t")[[1]][3] else ""
  status <- if (done) "done" else if (identical(last, "halt")) "halted" else if (identical(last, "start")) "in flight / interrupted" else "pending"
  cnt <- if (file.exists(grec)) sub("^GATE_COUNTS ", "", grep("^GATE_COUNTS", readLines(grec), value = TRUE)) else ""
  if (done && length(cnt) && nzchar(cnt)) {
    v <- as.integer(sub(".*=", "", strsplit(cnt, " ")[[1]])); tot <- tot + v
  }
  wall <- if (file.exists(rep)) sub(".*: ([0-9]+) s\\.$", "\\1", grep("^- Cell wall-clock", readLines(rep), value = TRUE)) else ""
  reps <- if (file.exists(rep)) sub("^- Replicates combined: ([0-9]+)\\.$", "\\1", grep("^- Replicates combined", readLines(rep), value = TRUE)) else ""
  dsha <- if (length(grep("\tdone\t", hbc))) sub(".*commit=", "", tail(grep("\tdone\t", hbc, value = TRUE), 1)) else ""
  L <- c(L, sprintf("| %s | %.2f | %d | %s | %s | %s | %s | %s | %s |", ce, cells$hr[i], cells$n[i], status,
                    paste(wall, collapse = ""), paste(reps, collapse = ""), paste(cnt, collapse = ""), dsha,
                    if (file.exists(grec)) sprintf("`%s`", gname) else ""))
  if (dir.exists(pay)) for (f in sort(list.files(pay, pattern = paste0("^(", stem, "_|p12x20_", ce, "_|GATE2_", TAG, "_", ce, "_|REPORT_", TAG, "_", ce, "_)"))))
    inv <- c(inv, sprintf("- `%s/%s/%s`: %d B, %s", rq, paydir, f, file.size(file.path(pay, f)),
                          if (tracked(file.path(rq, paydir, f))) "tracked" else "untracked"))
}
L <- c(L, "", sprintf("- Gate counts over done cells: run=%d passed=%d failed=%d.", tot[["run"]], tot[["passed"]], tot[["failed"]]), "",
       "## Payload inventory (from the directory)", "", if (length(inv)) inv else "- none", "",
       "## Heartbeat (`LOG_p12x20_progress.txt`, verbatim)", "", "```", hb, "```", "")
writeLines(L, file.path(qmd, "STATUS_p12x20.md"))
cat("wrote", file.path(qmd, "STATUS_p12x20.md"), "\n")
