# current_status.md section 3 -- the payload inventory, regenerated FROM THE
# DIRECTORY (standing rule 2, scripts_dinamr/README.md), never from chat records.
#
# Every file under quarto/simulations/gbsg_020 is assigned to exactly one row,
# first match wins, so the rows sum to the total.  Sizes are APPARENT size
# (st_size).  Tracked status is `git ls-files`.  Excluded throughout: the
# gitignored _gateT_pre_template_files/ and .DS_Store.
#
# usage (from scripts_dinamr/): Rscript status_inventory.R   -> markdown on stdout
QMD_DIR <- normalizePath(Sys.getenv("DINAMR_QMD_DIR", unset = ".."))
old <- setwd(QMD_DIR); on.exit(setwd(old))
f <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
f <- f[!grepl("^_gateT_pre_template_files/", f) & basename(f) != ".DS_Store"]
tracked <- system2("git", c("ls-files", "--", "."), stdout = TRUE)
sz <- file.info(f)$size
# glob -> anchored regex, '*' not crossing '/', '**' crossing it
g2r <- function(g) { g <- gsub(".", "[.]", g, fixed = TRUE)
  g <- gsub("**", "\001", g, fixed = TRUE); g <- gsub("*", "[^/]*", g, fixed = TRUE)
  paste0("^", gsub("\001", ".*", g, fixed = TRUE), "$") }
rules <- list(
  c("results/*dinamr*.rds", "`dinamr` per-replicate bundles + metas (batch and combined)"),
  c("results/*grfmr*.rds", "`grfmr` per-replicate bundles + metas (batch and combined), and the `grfmrsmk` Stage 1 smoke"),
  c("results/*grfprobe*.rds", "`grfprobe` cost probes, 36 replicates each"),
  c("results/*_nomr_idsweep_*.rds", "`idsweep` per-replicate bundles (MR off, 500 replicates, one batch per run; all three engines)"),
  c("results/fs_*.rds", "FS bundles -- the comparator grid plus every earlier FS campaign"),
  c("results/*.rds", "other bundles in `results/`"),
  c("p12x20_2026-09-12/*_res_*.rds", "`p12x20` per-replicate batch bundles + metas (campaign directory, not `results/`)"),
  c("p12x20_2026-09-12/*_combined_*.rds", "`p12x20` combined bundles (campaign directory)"),
  c("p12x20_2026-09-12/*.html", "`p12x20` batch and combine renders (campaign directory, not the root)"),
  c("p12x20_2026-09-12/GATE2_*.txt", "`p12x20` per-cell Gate 2 records"),
  c("p12x20_2026-09-12/REPORT_*.md", "`p12x20` per-cell reports"),
  c("p12x20_2026-09-12/**", "`p12x20` campaign directory, anything else"),
  c("scripts_p12x20/logs/*", "`p12x20` runner, driver and per-render logs -- `WALL_SECONDS` / `CELL DONE wall=`"),
  c("scripts_p12x20/*", "`p12x20` runner, driver, Gate 2 checker, STATUS generator"),
  c("STATUS_p12x20.md", "`p12x20` campaign-scoped record (no pin claim)"),
  c("LOG_p12x20_progress.txt", "`p12x20` per-cell heartbeat"),
  c("mr_sweep/**", "`mr_sweep/` -- an earlier seed-table sweep, superseded, kept for provenance"),
  c("scripts_dinamr/logs/*", "per-render, driver and gate logs -- `WALL_SECONDS` / `CELL DONE wall=`"),
  c("scripts_dinamr/*.R", "drivers, checkers, projections, extractions (R)"),
  c("scripts_dinamr/*.sh", "render/campaign drivers and the closeout checker (shell)"),
  c("scripts_dinamr/*.py", "transplant / chunk-diff helpers (Python)"),
  c("scripts_dinamr/*.cells", "cell lists, one line per cell"),
  c("scripts_dinamr/*.rds", "saved derived objects (projections, walls, extracted tables)"),
  c("scripts_dinamr/*.md", "`README.md` -- the standing rules, and what each script is"),
  c("scripts_dinamr/*.txt", "captured script output (GRF mechanism probe)"),
  c("dinamr_*.html", "`dinamr` batch and combine renders"),
  c("grfmr_*.html", "`grfmr` batch and combine renders"),
  c("grfprobe_*.html", "`grfprobe` renders"),
  c("idsweep_*.html", "`idsweep` batch renders, one per cell-run"),
  c("probe_*.html", "Gate 1 cost-probe renders (`dinamr` era)"),
  c("summary_*.html", "summary rendered outputs, all campaigns"),
  c("summary_*.qmd", "summary sources, all campaigns"),
  c("sim_fs_maxeffCons_fb_mr_field_m1_template.qmd", "**the** template `dinamr` / `grfmr` / the FS grid all render"),
  c("sim_*.qmd", "other simulation templates (earlier campaigns and variants)"),
  c("sim_*.html", "renders of those other templates"),
  c("fs_*.html", "FS campaign batch/combine renders (`p12ext`, `tier2`, `e1stud`, `cert20`, earlier)"),
  c("*.html", "remaining renders (smoke, gate, dflt, compare)"),
  c("REPORT_*.md", "REPORT documents"),
  c("TABLES_*.md", "TABLES documents"),
  c("REVIEW_*.md", "REVIEW documents"),
  c("current_status.md", "this file -- the directory's catalog at a pin"),
  c("*.md", "other notes in the directory"),
  c("*.qmd", "remaining `.qmd`"),
  c("*.R", "top-level ad-hoc R scripts"),
  c("**", "everything else"))
assigned <- rep(NA_integer_, length(f))
for (i in seq_along(rules)) {
  hit <- is.na(assigned) & grepl(g2r(rules[[i]][1]), f)
  assigned[hit] <- i }
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
cat(sprintf("| **total** | **%d/%d** | **%.0f MB** | `%s` %s | every file, each counted once |\n",
    sum(f %in% tracked), length(f), sum(sz) / 1024^2, basename(f[j]), hs1(sz[j])))
cat(sprintf("\nFILES OVER 50 MB: %d ; OVER 100 MB: %d\n", sum(sz > 50 * 1024^2), sum(sz > 100 * 1024^2)))
if (any(!f %in% tracked)) cat("UNTRACKED:\n", paste0("  ", f[!f %in% tracked], "\n"), sep = "")
