# current_status.md section 3.5 -- the payload inventory of quarto/simulations/actg175/binary_020,
# regenerated FROM THE DIRECTORY, never from chat records.
# Transplant of ../../continuous/scripts_mdsgnb20/status_inventory.R
# (TASK_actg175_binary_campaign_2026-09-17 §3.4): the rule table names THIS directory's layout
# (the committed study driver and its mr_sweep/ payloads; the new campaigns' bundles under
# mr_or_harm/<stem>_d5000/, their combine renders in the root, scripts_or/ and logs_or/); the
# matching, sizing and tracked logic are unchanged.
# Every file is assigned to exactly one row, first match wins, so the rows sum to the total.
# Sizes are APPARENT size (st_size).  Tracked status is `git ls-files`.
# usage (from scripts_or/): Rscript status_inventory.R   -> markdown on stdout
QMD_DIR <- normalizePath(Sys.getenv("ORSG_QMD_DIR", unset = ".."))
old <- setwd(QMD_DIR); on.exit(setwd(old))
f <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
f <- f[basename(f) != ".DS_Store"]
tracked <- system2("git", c("ls-files", "--", "."), stdout = TRUE)
sz <- file.info(f)$size
g2r <- function(g) { g <- gsub(".", "[.]", g, fixed = TRUE)
  g <- gsub("**", "\001", g, fixed = TRUE); g <- gsub("*", "[^/]*", g, fixed = TRUE)
  paste0("^", gsub("\001", ".*", g, fixed = TRUE), "$") }
rules <- list(
  c("mr_or_harm/*_orfs_d5000/*_combined_*.rds",   "`orfs` combined bundles (2,000 replicates per cell)"),
  c("mr_or_harm/*_orfs_d5000/*_res_*.rds",        "`orfs` per-replicate batch bundles + metas"),
  c("mr_or_harm/*_orgrf_d5000/*_combined_*.rds",  "`orgrf` combined bundles (2,000 replicates per cell)"),
  c("mr_or_harm/*_orgrf_d5000/*_res_*.rds",       "`orgrf` per-replicate batch bundles + metas"),
  c("mr_or_harm/*_ordina_d5000/*_combined_*.rds", "`ordina` combined bundles (2,000 replicates per cell)"),
  c("mr_or_harm/*_ordina_d5000/*_res_*.rds",      "`ordina` per-replicate batch bundles + metas"),
  c("mr_or_harm/*_orsmoke_d5000/**",              "`orsmoke` smoke bundles (Stage 1; deleted at closeout)"),
  c("mr_or_harm/*_orcal*_d5000/**",               "`orcal` calibration bundles (Stage 1; deleted at closeout)"),
  c("mr_or_harm/**",                              "other bundles under `mr_or_harm/`"),
  c("mr_sweep/*/mr_coverage_grid_*.rds",          "the committed study's coverage grid (21 cells; read, never recomputed)"),
  c("mr_sweep/*/*.rds",                           "the committed study's 21 per-cell payloads (`<method>_mr_n<n>_res.rds`)"),
  c("mr_sweep/**",                                "other files under the committed study's `mr_sweep/`"),
  c("maxeffCons_mr_coverage_sweep_or075.qmd",     "**the committed study driver** (the supplement's Figures S9/S10; not edited)"),
  c("maxeffCons_mr_coverage_sweep_or075.html",    "the committed study driver's render"),
  c("_sim_mr_coverage_or075.qmd",                 "the study driver's figure fragment"),
  c("sim_fs_mr_field_or_template.qmd",            "**the** OR template (`orfs`, `orgrf`, `ordina`)"),
  c("scripts_or/*",                               "`or` runner, Gate 2 checker, smoke identity checker, memory sampler, catalog generator and checkers"),
  c("logs_or/**",                                 "raw render, smoke and calibration logs (untracked)"),
  c("fs_*_orfs_combine_1_2000.html",              "`orfs` combine renders, one per cell"),
  c("grf_*_orgrf_combine_1_2000.html",            "`orgrf` combine renders, one per cell"),
  c("dina_*_ordina_combine_1_2000.html",          "`ordina` combine renders, one per cell"),
  c("summary_actg175_or.qmd",                     "the cross-cell summary source (all three campaigns, six cells)"),
  c("summary_actg175_or.html",                    "the cross-cell summary render"),
  c("or_metrics.csv",                             "the long-format extract (all three identifiers, with the quoted study rows)"),
  c("COLUMNS_or.md",                              "column, metric, scale and reading definitions of the extract"),
  c("fig_or_*.png",                               "figures written by the summary and embedded in the record"),
  c("REPORT_actg175_or_stage1_*.md",              "Gate 1 record (template and script transplants, smoke, calibration, projection)"),
  c("REPORT_actg175_or_gate2_*.md",               "Gate 2 record (per cell)"),
  c("REPORT_actg175_or_2026-*.md",                "the campaigns' final record"),
  c("REPORT_*.md",                                "other REPORT documents"),
  c("LOG_or_progress.txt",                        "per-render heartbeat"),
  c("HALT_or.md",                                 "halt file (absent unless a campaign halted)"),
  c("current_status.md",                          "this file -- the directory's catalog at a pin"),
  c("status_curated.md",                          "the hand-maintained sections of `current_status.md`"),
  c("*.qmd",                                      "other documents in the directory"),
  c("*.html",                                     "renders of those documents"),
  c("*.R",                                        "top-level R scripts"),
  c("*.csv",                                      "other extracts"),
  c("*.md",                                       "other notes in the directory"),
  c("**",                                         "everything else"))
assigned <- rep(NA_integer_, length(f))
for (i in seq_along(rules)) { hit <- is.na(assigned) & grepl(g2r(rules[[i]][1]), f); assigned[hit] <- i }
stopifnot(!anyNA(assigned))
hs  <- function(b) if (b >= 1024^2) sprintf("%.2f MB", b / 1024^2) else sprintf("%.0f KB", b / 1024)
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
