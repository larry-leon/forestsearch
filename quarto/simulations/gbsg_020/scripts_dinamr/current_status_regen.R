# Regenerates quarto/simulations/gbsg_020/current_status.md.
#
# Two sources, never mixed:
#   * GENERATED -- the pin, the campaign layout, cells, git and the payload
#     inventory, from bundle metas, file listings and git (TASK_p12x20_merge_2026-09-13
#     §5); never from a prior status file, a report or a chat record.
#   * CURATED -- status_curated.md, hand-maintained and included VERBATIM
#     (TASK_status_curated_restore_2026-09-13).  Its four marked blocks are placed
#     where they sat relative to the inventory at 0ab5d1c5: the preamble bullets
#     after the pin, §1-§2 before §3, the table commentary straight after the
#     inventory table, §4-§7 after §3.  A missing file or marker is a hard error:
#     current_status.md is never written without the curated sections.
#
# The pin is HEAD at generation.  Commit the file ALONE, as HEAD's child, then run
#   bash scripts_dinamr/check_current_status.sh --commit
# usage (from scripts_dinamr/): Rscript current_status_regen.R
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
qmd  <- normalizePath(file.path(here, ".."))
repo <- system2("git", c("-C", qmd, "rev-parse", "--show-toplevel"), stdout = TRUE)
rq   <- sub(paste0("^", repo, "/"), "", qmd)
# system2() does not quote: shQuote every argument, or a format string with spaces splits.
git  <- function(...) suppressWarnings(system2("git", shQuote(c("-C", repo, ...)), stdout = TRUE, stderr = FALSE))
`%||%` <- function(a, b) if (is.null(a) || !length(a) || all(is.na(a))) b else a
old <- setwd(qmd); on.exit(setwd(old))

# ---- curated source: fail loudly --------------------------------------------
cur_f <- file.path(qmd, "status_curated.md")
if (!file.exists(cur_f))
  stop("status_curated.md not found at ", cur_f,
       " -- refusing to write current_status.md without the curated sections.", call. = FALSE)
cur   <- readLines(cur_f, encoding = "UTF-8", warn = TRUE)
marks <- c("preamble", "before-inventory", "after-inventory-table", "after-inventory")
pos   <- match(sprintf("<!-- curated:%s -->", marks), cur)
if (anyNA(pos) || is.unsorted(pos))
  stop("status_curated.md: markers missing or out of order (",
       paste(marks[is.na(pos)], collapse = ", "), ") -- refusing to write current_status.md.", call. = FALSE)
C <- lapply(seq_along(pos), function(i) {
  a <- pos[i] + 1L; b <- if (i < length(pos)) pos[i + 1L] - 1L else length(cur)
  if (b < a) character(0) else cur[a:b] })
names(C) <- marks
if (any(lengths(C) == 0L))
  stop("status_curated.md: empty block(s): ", paste(marks[lengths(C) == 0L], collapse = ", "), call. = FALSE)

pin    <- git("rev-parse", "--short", "HEAD")
branch <- git("rev-parse", "--abbrev-ref", "HEAD")
trk_all <- git("ls-files", "--", rq)
trk <- function(p) file.path(rq, p) %in% trk_all
MB  <- function(b) sprintf("%.2f MB", b / 1024^2)
tick <- function(x) paste0("`", x, "`")
PAYDIR <- "p12x20_2026-09-12"

camps <- list(
  dinamr = list(bdir = "results", bre = "^dina_.*_dinamr_(res|combined)_[0-9]+_[0-9]+[.]rds$",
                rdir = ".", rre = "^dinamr_.*_(batch|combine)_[0-9]+[.]html$",
                scripts = "scripts_dinamr", sre = "dinamr|campaign[.]sh|gate2[.]R|block|probe[.]sh|project[.]R|checkpoint",
                docre = "^(REPORT|REVIEW|TABLES|summary)_dinamr[_.]"),
  grfmr  = list(bdir = "results", bre = "^grf_.*_grfmr_(res|combined)_[0-9]+_[0-9]+[.]rds$",
                rdir = ".", rre = "^grfmr_.*_(batch|combine)_[0-9]+[.]html$",
                scripts = "scripts_dinamr", sre = "grfmr|gate2G|gate3|stage1G|projectG|wallsGC|t2gate|t3gate",
                docre = "^(REPORT|REVIEW|TABLES|summary)_grfmr[_.]"),
  p12x20 = list(bdir = PAYDIR, bre = "^fs_.*_p12x20_(res|combined)_[0-9]+_[0-9]+[.]rds$",
                rdir = PAYDIR, rre = "^p12x20_A[0-9]_.*_(batch|combine)_[0-9]+[.]html$",
                scripts = "scripts_p12x20", sre = ".",
                docre = "^(REPORT|STATUS|LOG)_p12x20[_.]"),
  idsweep = list(bdir = "results", bre = "^(fs|dina|grf)_.*_nomr_idsweep_(res|combined)_[0-9]+_[0-9]+[.]rds$",
                rdir = ".", rre = "^idsweep_.*[.]html$",
                scripts = "scripts_dinamr", sre = "idsweep",
                docre = "^(REPORT|summary)_idsweep[_.]"),
  # nullid (TASK_null_gbsg_identification_2026-09-21): the structural-null cell.
  # Batch-only like idsweep -- one 2,000-replicate batch per cell-run, MR off,
  # no combine -- so its engine line and its cells come from the batch metas.
  # The smoke (`nullsmk`) and inertness (`nullinert`) bundles are quickrun-tagged
  # and are deliberately OUTSIDE bre, so they never enter the cell table.
  nullid = list(bdir = "results", bre = "^(fs|dina|grf)_.*_nomr_nullid_(res|combined)_[0-9]+_[0-9]+[.]rds$",
                rdir = ".", rre = "^nullid_.*[.]html$",
                scripts = "scripts_dinamr", sre = "nullid",
                # Tightened for nullc125: "_null_gbsg[_.]" also matched that
                # campaign's REPORT_null_gbsg_thresholds_*.
                docre = "^(REPORT|summary)_null_gbsg_identification[_.]"),
  # nullc125 (TASK_null_gbsg_thresholds_2026-09-21_v2): the structural null at
  # the stricter screen c1 1.25 / c2 1.00.  Batch-only, as nullid.  Its smoke
  # (`nullc125smoke`) and inertness (`nullc125inert*`) bundles are
  # quickrun-tagged and outside bre; `^nullc125_` does not match their renders.
  nullc125 = list(bdir = "results", bre = "^(fs|dina|grf)_.*_nomr_c125c100_nullc125_(res|combined)_[0-9]+_[0-9]+[.]rds$",
                rdir = ".", rre = "^nullc125_.*[.]html$",
                scripts = "scripts_dinamr", sre = "nullthr|nullc125",
                docre = "^(REPORT|summary)_null_gbsg_thresholds[_.]"))

meta_rows <- list(); nav <- character(0)
for (cn in names(camps)) {
  cc <- camps[[cn]]
  bf <- list.files(cc$bdir, pattern = cc$bre); comb <- grep("_combined_", bf, value = TRUE); bat <- grep("_res_", bf, value = TRUE)
  rf <- list.files(cc$rdir, pattern = cc$rre); rsz <- sum(file.size(file.path(cc$rdir, rf)))
  sf <- list.files(cc$scripts); sf <- sf[grepl(cc$sre, sf) & !dir.exists(file.path(cc$scripts, sf))]
  docs <- list.files(".", pattern = cc$docre)
  ms <- lapply(comb, function(f) readRDS(file.path(cc$bdir, f))$meta)
  # A batch-only campaign (idsweep: one 500-replicate batch per cell-run, no
  # combine) has no combined metas; its engine line and its cells come from the
  # batch metas instead.
  batch_only <- !length(comb) && length(bat) > 0L
  bms <- if (batch_only) lapply(bat, function(f) readRDS(file.path(cc$bdir, f))) else list()
  eng <- if (batch_only) {
    bm_ <- lapply(bms, `[[`, "meta")
    sprintf("%s / %s / eps %s (batch metas; MR %s)",
            paste(unique(vapply(bm_, function(m) m$subgroup_method, "")), collapse = ", "),
            paste(unique(vapply(bm_, function(m) m$sg_focus, "")), collapse = ", "),
            paste(sort(unique(vapply(bm_, function(m) format(m$effect_neighborhood), ""))), collapse = ", "),
            paste(unique(vapply(bm_, function(m) if (isTRUE(m$mr_inference)) "on" else "off", "")), collapse = "/"))
  } else paste(unique(vapply(ms, function(m) sprintf("%s / %s / eps %s", m$subgroup_method, m$sg_focus, format(m$effect_neighborhood)), "")), collapse = "; ")
  bpath <- if (cc$bdir == "results") sprintf("`results/` (in `gbsg_020/`)") else sprintf("`%s/`", cc$bdir)
  rpath <- if (cc$rdir == ".") "`gbsg_020/` root" else sprintf("`%s/`", cc$rdir)
  nav <- c(nav, sprintf("| `%s` | %s | %d | %d | %s | %d, %s, in %s | `%s/` (%d files named for it) | %s |",
    cn, eng, length(comb), length(bat), bpath, length(rf), MB(rsz), rpath, cc$scripts, length(sf),
    paste(tick(docs), collapse = ", ")))
  for (k in seq_along(comb)) {
    f <- comb[k]; m <- ms[[k]]; b <- readRDS(file.path(cc$bdir, f))
    stem <- sub("_combined_[0-9]+_[0-9]+[.]rds$", "", f)
    bm <- lapply(list.files(cc$bdir, pattern = paste0("^", stem, "_res_")), function(x) readRDS(file.path(cc$bdir, x))$meta)
    h3 <- sprintf("h%03d", round(100 * m$target_hr_harm))
    cell <- if (cn == "p12x20") sub("^p12x20_(A[0-9])_.*$", "\\1",
               grep(sprintf("_%s_n%d_combine_", h3, m$n_sample), rf, value = TRUE)[1] %||% NA_character_) else "—"
    meta_rows[[length(meta_rows) + 1]] <- data.frame(
      campaign = cn, cell = cell, prev = m$harm_prevalence_super, hr = m$target_hr_harm, n = m$n_sample,
      rows = nrow(b$results), sim = paste(range(b$results$sim_id), collapse = "-"), batches = length(bm),
      host = paste(unique(vapply(bm, function(x) x$hostname %||% NA_character_, "")), collapse = "/"),
      workers = paste(unique(vapply(bm, function(x) as.character(x$n_workers %||% NA), "")), collapse = "/"),
      version = paste(unique(vapply(bm, function(x) x$forestsearch_version %||% NA_character_, "")), collapse = "/"),
      path = file.path(cc$bdir, f), size = file.size(file.path(cc$bdir, f)), tracked = trk(file.path(cc$bdir, f)),
      stringsAsFactors = FALSE)
  }
  if (batch_only) {
    key <- vapply(bms, function(b) with(b$meta, sprintf("%.5f|%.6f|%d", harm_prevalence_super, target_hr_harm, n_sample)), "")
    for (kk in unique(key)) {
      ii <- which(key == kk); m <- bms[[ii[1]]]$meta; bm_ <- lapply(bms[ii], `[[`, "meta")
      meta_rows[[length(meta_rows) + 1]] <- data.frame(
        campaign = cn, cell = "—", prev = m$harm_prevalence_super, hr = m$target_hr_harm, n = m$n_sample,
        rows = nrow(bms[[ii[1]]]$results),
        sim = paste(range(unlist(lapply(bms[ii], function(b) b$results$sim_id))), collapse = "-"), batches = length(ii),
        host = paste(unique(vapply(bm_, function(x) x$hostname %||% NA_character_, "")), collapse = "/"),
        workers = paste(unique(vapply(bm_, function(x) as.character(x$n_workers %||% NA), "")), collapse = "/"),
        version = paste(unique(vapply(bm_, function(x) x$forestsearch_version %||% NA_character_, "")), collapse = "/"),
        path = file.path(cc$bdir, bat[ii[1]]), size = sum(file.size(file.path(cc$bdir, bat[ii]))),
        tracked = all(trk(file.path(cc$bdir, bat[ii]))), stringsAsFactors = FALSE)
    }
  }
}
cells <- do.call(rbind, meta_rows)
cells <- cells[order(cells$campaign, cells$prev, cells$hr, cells$n), ]
# Two decimals is right for the alt grid (1.00 / 1.50 / 1.75) and WRONG for the
# structural-null design points (0.657 / 0.721), which it rounds to 0.66 / 0.72.
# Print the third decimal only when it carries information.
.hrfmt <- function(x) if (isTRUE(all.equal(x, round(x, 2)))) sprintf("%.2f", x) else sprintf("%.3f", x)

# ---- preamble, then curated §1-§2 --------------------------------------------
L <- c(sprintf("# current_status — `%s`", rq), "",
  sprintf("- **Pin:** `%s` on `%s` — HEAD at the time this file was committed; the closeout commit that adds this file is its child.", pin, branch),
  "- **Updated:** 2026-09-21 (closeout of `TASK_null_gbsg_thresholds_2026-09-21_v2` — the structural null at the stricter screen c1 1.25 / c2 1.00, `nullc125`, \u00a72.9: 6 cells \u00d7 3 identifiers \u00d7 2,000 replicates, identification and classification only, no MR anywhere; curated sections from `status_curated.md`).",
  "- **Generated by** `scripts_dinamr/current_status_regen.R`: the pin and §3 from the directory and from git (§3's table by `scripts_dinamr/status_inventory.R`). Everything else is included verbatim from the hand-maintained `status_curated.md`.",
  C[["preamble"]],
  "", "---", "",
  C[["before-inventory"]])

# ---- §3, generated -----------------------------------------------------------
L <- c(L, "## 3. Payload inventory", "",
  "Regenerated from the directory and from git by `scripts_dinamr/current_status_regen.R`; the table in 3.5 by `scripts_dinamr/status_inventory.R`. The commentary after that table, and §1, §2 and §4–§7, are included verbatim from `status_curated.md`.", "",
  "### 3.1 Where each campaign lives", "",
  "**The layouts differ.** `dinamr`, `grfmr` and `idsweep` keep bundles in `results/` and renders in the `gbsg_020/` root. `idsweep` has no combined bundles: it is one 500-replicate batch per cell-run, MR off. `p12x20` keeps bundles, renders, per-cell gate records and per-cell reports together in `p12x20_2026-09-12/`, and its campaign-scoped record is `STATUS_p12x20.md`.", "",
  "| campaign | engine / focus / ε (combined metas) | combined bundles | batch bundles | bundles in | renders (batch + combine) | scripts | campaign documents in `gbsg_020/` |",
  "|---|---|---|---|---|---|---|---|", nav, "",
  sprintf("- `p12x20` per-cell gate records and reports: %d `GATE2_p12x20_*` and %d `REPORT_p12x20_A*` files in `%s/`.",
          length(list.files(PAYDIR, "^GATE2_p12x20_")), length(list.files(PAYDIR, "^REPORT_p12x20_A")), PAYDIR),
  "- `p12x20` campaign-scoped record: `STATUS_p12x20.md` (cells, gate counts, payload inventory, Stage 0 install record); per-cell heartbeat `LOG_p12x20_progress.txt`; runner logs `scripts_p12x20/logs/`.",
  "- `dinamr`, `grfmr` and `idsweep` share `scripts_dinamr/`; `scripts_dinamr/README.md` states each script's role. The script counts above are files whose names match the campaign. `idsweep`'s record is `REPORT_idsweep_2026-09-12.md` with `summary_idsweep.html`; its logs are `scripts_dinamr/logs/idsweep*`.",
  "", "### 3.2 Cells, from the combined bundle metas (from the batch metas for `idsweep`, which has none)", "")
for (cn in names(camps)) {
  d <- cells[cells$campaign == cn, ]
  bo <- nrow(d) > 0L && !any(grepl("_combined_", d$path))
  L <- c(L, sprintf("#### `%s` — %d cells%s", cn, nrow(d), if (bo) " (no combined bundles: one batch bundle per cell-run)" else ""), "",
    if (bo) "| cell | prevalence (super-population) | HR | n | rows per run | sim_id | runs (batch bundles) | host | workers | forestsearch | first batch bundle | total size | tracked (all) |" else
    "| cell | prevalence (super-population) | HR | n | rows | sim_id | batches | host | workers | forestsearch | combined bundle | size | tracked |",
    "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    sprintf("| %s | %.5f | %s | %d | %d | %s | %d | %s | %s | %s | `%s` | %s | %s |", d$cell, d$prev, .hrfmt(d$hr), d$n, d$rows, d$sim,
            d$batches, d$host, d$workers, d$version, d$path, MB(d$size), ifelse(d$tracked, "yes", "no")), "")
}
seeds <- unique(unlist(lapply(cells$path, function(p) readRDS(p)$meta$seed_base)))
L <- c(L, sprintf("- `seed_base` in every meta above (combined metas; the first batch meta per `idsweep` cell): %s (per replicate `seed_base + sim_id`).", paste(seeds, collapse = ", ")), "")

L <- c(L, "### 3.3 Git", "")
for (cn in names(camps)) {
  specs <- switch(cn,
    dinamr = c(file.path(rq, "results/*_dinamr_*"), file.path(rq, "dinamr_*.html")),
    grfmr  = c(file.path(rq, "results/*_grfmr_*"),  file.path(rq, "grfmr_*.html")),
    p12x20 = c(file.path(rq, PAYDIR), file.path(rq, "scripts_p12x20"), file.path(rq, "STATUS_p12x20.md")),
    idsweep = c(file.path(rq, "results/*_nomr_idsweep_*"), file.path(rq, "idsweep_*.html")),
    nullid  = c(file.path(rq, "results/*_nomr_nullid_*"),  file.path(rq, "nullid_*.html")),
    nullc125 = c(file.path(rq, "results/*_nomr_c125c100_nullc125_*"), file.path(rq, "nullc125_*.html")))
  # A missing branch here makes switch() return NULL, and `git log --` with an
  # empty pathspec logs the ENTIRE REPOSITORY.  Guard it: a campaign with no
  # pathspec is a bug in this registry, not a campaign with thousands of commits.
  if (is.null(specs) || !length(specs))
    stop("current_status_regen.R: no git pathspec for campaign '", cn,
         "' -- add one to `specs` above.", call. = FALSE)
  lg <- git("log", "--format=%h %ad %s", "--date=short", "--", specs)
  L <- c(L, sprintf("- `%s` bundles/renders: %d commits; first `%s`; last `%s`.", cn, length(lg),
                    sub(" .*", "", tail(lg, 1) %||% "none"), sub(" .*", "", lg[1] %||% "none")))
}
# The campaign branch is local on the machine that ran it and remote-tracking elsewhere.
cref <- "campaign/p12x20"; cp <- git("rev-parse", "--verify", "--quiet", "--short", cref)
if (!length(cp)) { cref <- "origin/campaign/p12x20"; cp <- git("rev-parse", "--verify", "--quiet", "--short", cref) }
anc <- length(cp) && identical(system2("git", c("-C", repo, "merge-base", "--is-ancestor", cref, "HEAD")), 0L)
L <- c(L, sprintf("- `%s` tip `%s`; ancestor of HEAD: %s.", cref, cp %||% "absent", if (anc) "yes" else "no"), "")

rn <- list.files("results", pattern = "_(res|combined)_[0-9]+_[0-9]+[.]rds$")
tg <- sub("^.*_([^_]+)_(res|combined)_[0-9]+_[0-9]+[.]rds$", "\\1", rn); kd <- sub("^.*_(res|combined)_[0-9]+_[0-9]+[.]rds$", "\\1", rn)
tt <- table(tg, kd); tags <- rownames(tt)
L <- c(L, "### 3.4 Campaign tags in `results/` (last token of each bundle stem)", "",
  "| tag | combined | batch |", "|---|---|---|",
  sprintf("| `%s` | %d | %d |", tags, if ("combined" %in% colnames(tt)) tt[, "combined"] else 0L, if ("res" %in% colnames(tt)) tt[, "res"] else 0L), "",
  sprintf("- `%s/`: %d combined, %d batch (tag `p12x20`).", PAYDIR,
          length(list.files(PAYDIR, "_p12x20_combined_")), length(list.files(PAYDIR, "_p12x20_res_"))), "")

inv <- local({ o <- setwd(here); on.exit(setwd(o)); system2("Rscript", "status_inventory.R", stdout = TRUE) })
L <- c(L, "### 3.5 Payload inventory table", "", inv, "",
  C[["after-inventory-table"]], "",
  C[["after-inventory"]])
writeLines(L, "current_status.md")
cat("wrote", file.path(qmd, "current_status.md"), "pin", pin, "\n")
