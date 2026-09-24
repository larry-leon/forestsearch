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
  c("results/*_nomr_nullid_*.rds", "`nullid` per-replicate bundles (structural null, MR off, 2,000 replicates, one batch per run; all three engines)"),
  c("results/*_nomr_c125c100_nullc125_*.rds", "`nullc125` per-replicate bundles (structural null at c1 1.25 / c2 1.00, MR off, 2,000 replicates, one batch per run; all three engines)"),
  c("results/*_nullc125*_quickrun_*.rds", "`nullc125` Step 3.1 smoke (`nullc125smoke`, at 1.25 / 1.00) and the Step 1.5 threshold-knob inertness renders (`nullc125inert*`, at 0.90 / 0.80), 20 replicates each"),
  c("results/*_nb20_nullmr_*.rds", "`nullmr` per-replicate bundles (structural null at c1 0.90 / c2 0.80, MR on, 2,000 replicates, one batch per run; all three engines)"),
  c("results/*_nullmrsmoke_quickrun_*.rds", "`nullmr` Step 3.1 smoke (`nullmrsmoke`, MR on), 20 replicates each"),
  c("results/declcalc0_smoke_*.rds", "`declcalc0` Stage 1 smoke (cell B5, 5 replicates)"),
  c("results/declcalc0_*.rds", "`declcalc0` per-replicate payloads (the `declcal` B and C blocks re-run with `declaration_c0` = 0.70 / 0.75 / 0.80 / 0.85; `declcal` schema + per-c0 columns + aux + meta; 2,000 replicates per cell)"),
  c("results/declcal_c0approx_*.rds", "`declcal` c0 approximate table (120 field captures, B 2000; plug-in fixed-cutoff rates)"),
  c("results/declcal_pilot_*.rds", "`declcal` Stage 1 pilot (cell A2, 200 replicates, B 2000 sub-sampled)"),
  c("results/declcal_*.rds", "`declcal` per-replicate payloads (section-5 schema + aux + meta; FS only, c1 = c2 = 1.0, 2,000 replicates per cell; blocks bnull / inull / power)"),
  c("results/*null*_quickrun_*.rds", "`nullid` Step 2 smoke (`nullsmk`) and the alt-path inertness check (`nullinert`), 20 replicates each"),
  c("results/*_mrs5probe_*.rds", "`mrs5` probe: A7 after arm (200 replicates) and its Gate T (10)"),
  c("results/*_mrs5sweep_*.rds", "`mrs5` sweep: after arm, 17 Section 5 cells, 200 replicates each"),
  c("results/*_mrs5pre_*.rds", "`mrs5` parent-build (`ba595f4b`) checks, 10 replicates: A7, cert20 HR 1.00 n 500, e1stud HR 1.50 n 500"),
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
  c("mrs5probe/**", "`mrs5` probe runner, paired read-out, logs (`WALL_SECONDS`, worker build checks)"),
  c("mrs5sweep/**", "`mrs5` sweep driver, runner, paired read-out, cell list, logs (`WALL_SECONDS`, worker build checks)"),
  c("mrs5pre/**", "`mrs5` parent-build check logs"),
  c("mr_sweep/**", "`mr_sweep/` -- an earlier seed-table sweep, superseded, kept for provenance"),
  c("scripts_dinamr/logs/*", "per-render, driver and gate logs -- `WALL_SECONDS` / `CELL DONE wall=`"),
  c("scripts_dinamr/*.R", "drivers, checkers, projections, extractions (R)"),
  c("scripts_dinamr/*.sh", "render/campaign drivers and the closeout checker (shell)"),
  c("scripts_dinamr/*.py", "transplant / chunk-diff helpers (Python)"),
  c("scripts_dinamr/*.cells", "cell lists, one line per cell"),
  c("scripts_dinamr/resume/*.cells", "`declcalc0` Stage 2 resume cell lists (same basenames as the block files, so the campaign tags are unchanged; the cells not completed by the stopped run)"),
  c("scripts_dinamr/*.rds", "saved derived objects (projections, walls, extracted tables)"),
  c("scripts_dinamr/*.md", "`README.md` -- the standing rules, and what each script is"),
  c("scripts_dinamr/*.txt", "captured script output (GRF mechanism probe)"),
  c("mrs5*_*.html", "`mrs5` probe, sweep and parent-build check renders"),
  c("dinamr_*.html", "`dinamr` batch and combine renders"),
  c("grfmr_*.html", "`grfmr` batch and combine renders"),
  c("grfprobe_*.html", "`grfprobe` renders"),
  c("idsweep_*.html", "`idsweep` batch renders, one per cell-run"),
  c("nullid_*.html", "`nullid` batch renders, one per cell-run"),
  c("nullsmk_*.html", "`nullid` Step 2 smoke renders"),
  c("nullinert_*.html", "the `nullid` alt-path inertness-check render"),
  c("nullc125_*.html", "`nullc125` batch renders, one per cell-run"),
  c("nullc125smoke_*.html", "`nullc125` Step 3.1 smoke renders"),
  c("nullmr_*.html", "`nullmr` batch renders, one per cell-run"),
  c("nullmrsmoke_*.html", "`nullmr` Step 3.1 smoke renders"),
  c("nullc125inert*_*.html", "the `nullc125` threshold-knob inertness renders (pre-edit, unset, explicit)"),
  c("probe_*.html", "Gate 1 cost-probe renders (`dinamr` era)"),
  c("summary_*.html", "summary rendered outputs, all campaigns"),
  c("summary_*.qmd", "summary sources, all campaigns"),
  c("sim_fs_maxeffCons_fb_mr_field_m1_template.qmd", "**the** template `dinamr` / `grfmr` / the FS grid all render"),
  c("sim_*.qmd", "other simulation templates (earlier campaigns and variants)"),
  c("sim_*.html", "renders of those other templates"),
  c("fs_*.html", "FS campaign batch/combine renders (`p12ext`, `tier2`, `e1stud`, `cert20`, earlier)"),
  c("*.html", "remaining renders (smoke, gate, dflt, compare)"),
  c("HALT_*.txt", "halt records (a cell abandoned mid-campaign)"),
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
