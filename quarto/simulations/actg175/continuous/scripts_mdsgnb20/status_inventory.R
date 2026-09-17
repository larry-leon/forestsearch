# current_status.md section 3.5 -- the payload inventory of quarto/simulations/actg175/continuous,
# regenerated FROM THE DIRECTORY, never from chat records.
# Transplant of quarto/simulations/gbsg_020/scripts_dinamr/status_inventory.R (TASK_md_field_rerun_2026-09-15 §3.4):
# the rule table names this directory's layout (bundles under mr_md_harm/<stem>_d5000/, renders in the root,
# scripts_mdf1/ and scripts_mdsgnb20/); the matching, sizing and tracked logic are unchanged.
# Every file is assigned to exactly one row, first match wins, so the rows sum to the total.  Sizes are
# APPARENT size (st_size).  Tracked status is `git ls-files`.
# usage (from scripts_mdsgnb20/): Rscript status_inventory.R   -> markdown on stdout
QMD_DIR <- normalizePath(Sys.getenv("MDSG_QMD_DIR", unset = ".."))
old <- setwd(QMD_DIR); on.exit(setwd(old))
f <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
f <- f[basename(f) != ".DS_Store"]
tracked <- system2("git", c("ls-files", "--", "."), stdout = TRUE)
sz <- file.info(f)$size
g2r <- function(g) { g <- gsub(".", "[.]", g, fixed = TRUE)
  g <- gsub("**", "\001", g, fixed = TRUE); g <- gsub("*", "[^/]*", g, fixed = TRUE)
  paste0("^", gsub("\001", ".*", g, fixed = TRUE), "$") }
rules <- list(
  c("mr_md_harm/*_mdsgnb20_d5000/*_combined_*.rds", "`mdsgnb20` combined bundles (2,000 replicates per cell)"),
  c("mr_md_harm/*_mdsgnb20_d5000/*_res_*.rds", "`mdsgnb20` per-replicate batch bundles + metas"),
  # TASK_md_dina_campaign_2026-09-17 (carried fix): mdgrf and mddina rows, copied from mdsgnb20's.
  c("mr_md_harm/*_mdgrf_d5000/*_combined_*.rds", "`mdgrf` combined bundles (2,000 replicates per cell)"),
  c("mr_md_harm/*_mdgrf_d5000/*_res_*.rds", "`mdgrf` per-replicate batch bundles + metas"),
  c("mr_md_harm/*_mddina_d5000/*_combined_*.rds", "`mddina` combined bundles (2,000 replicates per cell)"),
  c("mr_md_harm/*_mddina_d5000/*_res_*.rds", "`mddina` per-replicate batch bundles + metas"),
  c("mr_md_harm/*_mdf1_d5000/*_combined_*.rds", "`mdf1` combined bundles (2,000 replicates per cell)"),
  c("mr_md_harm/*_mdf1_d5000/*_res_*.rds", "`mdf1` per-replicate batch bundles + metas"),
  c("mr_md_harm/*_mdf1_d5000/gate2_flips.txt", "`mdf1` Gate 2 enumerated rows per cell"),
  c("mr_md_harm/*_mdsmoke_d5000/**", "`mdsmoke` smoke bundles (Stage 1; deleted at closeout)"),
  c("mr_md_harm/*_mdcal*_d5000/**", "`mdcal` calibration bundles (Stage 1; deleted at closeout)"),
  c("mr_md_harm/*_s1000_d5000/**", "the continuous twin's committed MR (IJ) bundles (`sim_fs_maxeffCons_mr_md*_batch_1_1000.qmd`)"),
  c("mr_md_harm/*_s100_d5000/**", "the twin's 100-replicate bundle"),
  c("mr_md_harm/**", "other bundles under `mr_md_harm/`"),
  c("fb_mr_md_harm/**", "FB bundles (the twin's FB runs; joined by `mdf1` on md40 n500 sim_id 1-100, never re-run)"),
  c("scripts_mdsgnb20/*", "`mdsgnb20` runner, Gate 2 checker, smoke identity checker, memory sampler, catalog generator and checkers"),
  c("scripts_mdf1/*", "`mdf1` session scripts (Mac; records of what ran)"),
  c("logs_mdsgnb20/**", "`mdsgnb20` raw render, smoke and calibration logs (untracked)"),
  c("fs_*_mdsgnb20_combine_1_2000.html", "`mdsgnb20` combine renders, one per cell"),
  c("scripts_mdgrf/*", "`mdgrf` runner, Gate 2 checker, smoke checker, memory sampler"),
  c("logs_mdgrf/**", "`mdgrf` raw render, smoke and calibration logs (untracked)"),
  c("grf_*_mdgrf_combine_1_2000.html", "`mdgrf` combine renders, one per cell"),
  c("scripts_mddina/*", "`mddina` runner, Gate 2 checker, smoke checker, memory sampler"),
  c("logs_mddina/**", "`mddina` raw render, smoke and calibration logs (untracked)"),
  c("dina_*_mddina_combine_1_2000.html", "`mddina` combine renders, one per cell"),
  c("fs_*_mdf1_combine_1_2000.html", "`mdf1` combine renders, one per cell"),
  c("summary_continuous_field_*.qmd", "cross-cell summary sources (`mdf1`, `mdsgnb20`)"),
  c("summary_continuous_field_*.html", "cross-cell summary renders"),
  c("sim_fs_maxeffCons_mr_field_md_template.qmd", "**the** MD template (`mdf1`, `mdsgnb20`)"),
  c("sim_*.qmd", "the continuous twin's batch documents (earlier campaigns)"),
  c("sim_*.html", "renders of those batch documents"),
  c("md_field_metrics.csv", "the `mdsgnb20` long-format extract"),
  c("COLUMNS_md_field.md", "column and metric definitions of the extract"),
  c("md_grf_metrics.csv", "the `mdgrf` long-format extract (with FS comparator rows)"),
  c("COLUMNS_md_grf.md", "column and metric definitions of the `mdgrf` extract"),
  c("md_dina_metrics.csv", "the `mddina` long-format extract (with FS and GRF comparator rows)"),
  c("COLUMNS_md_dina.md", "column and metric definitions of the `mddina` extract"),
  c("fig_*.png", "figures embedded in records"),
  c("REPORT_md_field_rerun_*.md", "`mdsgnb20` records (Stage 0, Stage 1, Gate 2, final)"),
  c("REPORT_md_grf_*.md", "`mdgrf` records (Stage 1, Stage 1 resumed, Gate 2, final)"),
  c("REPORT_md_dina_stage1_*.md", "`mddina` Gate 1 record"),
  c("REPORT_md_dina_gate2_*.md", "`mddina` Gate 2 record"),
  c("REPORT_md_dina_2026-*.md", "`mddina` final record"),
  c("REPORT_continuous_field_*.md", "`mdf1` records (Stage 0, Stage 1, Gate 2, Stage 3)"),
  c("REPORT_*.md", "other REPORT documents"),
  c("LOG_mdsgnb20_progress.txt", "`mdsgnb20` per-render heartbeat"),
  c("HALT_mdsgnb20.md", "`mdsgnb20` halt file (absent unless the campaign halted)"),
  c("LOG_mdgrf_progress.txt", "`mdgrf` per-render heartbeat"),
  c("HALT_mdgrf.md", "`mdgrf` halt file (absent unless the campaign halted)"),
  c("LOG_mddina_progress.txt", "`mddina` per-render heartbeat"),
  c("HALT_mddina.md", "`mddina` halt file (absent unless the campaign halted)"),
  c("current_status.md", "this file -- the directory's catalog at a pin"),
  c("status_curated.md", "the hand-maintained sections of `current_status.md`"),
  c("*.qmd", "other documents (analytic derivations, OC wrapper verification, coverage sweeps, the simulations overview)"),
  c("*.html", "renders of those documents"),
  c("*.R", "top-level R scripts (truth, readouts, batch drivers)"),
  c("*.md", "other notes in the directory"),
  c("**", "everything else"))
assigned <- rep(NA_integer_, length(f))
for (i in seq_along(rules)) { hit <- is.na(assigned) & grepl(g2r(rules[[i]][1]), f); assigned[hit] <- i }
stopifnot(!anyNA(assigned))
hs <- function(b) if (b >= 1024^2) sprintf("%.2f MB", b / 1024^2) else sprintf("%.0f KB", b / 1024)
hs1 <- function(b) if (b >= 1024^2) sprintf("%.2f MB", b / 1024^2) else sprintf("%d KB", as.integer(round(b / 1024)))
cat("| path pattern (first match wins) | tracked/disk | total | largest single file | what it is |\n")
cat("|---|---|---|---|---|\n")
for (i in seq_along(rules)) {
  k <- which(assigned == i); if (!length(k)) next
  j <- k[which.max(sz[k])]
  cat(sprintf("| `%s` | %d/%d | %s | `%s` %s | %s |\n", rules[[i]][1], sum(f[k] %in% tracked),
      length(k), hs(sum(sz[k])), basename(f[j]), hs1(sz[j]), rules[[i]][2])) }
j <- which.max(sz)
cat(sprintf("| **all** | **%d/%d** | **%s** | `%s` %s | |\n", sum(f %in% tracked), length(f), hs(sum(sz)), basename(f[j]), hs1(sz[j])))
