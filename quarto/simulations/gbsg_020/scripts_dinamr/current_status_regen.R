# Regenerates quarto/simulations/gbsg_020/current_status.md FROM THE DIRECTORY AND
# FROM GIT (TASK_p12x20_merge_2026-09-13 §5) -- never from a prior status file, a
# report, or a chat record.  Numbers come from bundle metas, file listings and git;
# §4 is status_inventory.R's output.
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
                docre = "^(REPORT|STATUS|LOG)_p12x20[_.]"))

meta_rows <- list(); nav <- character(0)
for (cn in names(camps)) {
  cc <- camps[[cn]]
  bf <- list.files(cc$bdir, pattern = cc$bre); comb <- grep("_combined_", bf, value = TRUE); bat <- grep("_res_", bf, value = TRUE)
  rf <- list.files(cc$rdir, pattern = cc$rre); rsz <- sum(file.size(file.path(cc$rdir, rf)))
  sf <- list.files(cc$scripts); sf <- sf[grepl(cc$sre, sf) & !dir.exists(file.path(cc$scripts, sf))]
  docs <- list.files(".", pattern = cc$docre)
  ms <- lapply(comb, function(f) readRDS(file.path(cc$bdir, f))$meta)
  eng <- paste(unique(vapply(ms, function(m) sprintf("%s / %s / eps %s", m$subgroup_method, m$sg_focus, format(m$effect_neighborhood)), "")), collapse = "; ")
  bpath <- if (cc$bdir == "results") sprintf("`results/` (in `gbsg_020/`)") else sprintf("`%s/`", cc$bdir)
  rpath <- if (cc$rdir == ".") "`gbsg_020/` root" else sprintf("`%s/`", cc$rdir)
  nav <- c(nav, sprintf("| `%s` | %s | %d | %d | %s | %d, %s, in %s | `%s/` (%d files named for it) | %s |",
    cn, eng, length(comb), length(bat), bpath, length(rf), MB(rsz), rpath, cc$scripts, length(sf),
    paste(tick(docs), collapse = ", ")))
  # per-cell rows from the combined metas and their batch metas
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
}
cells <- do.call(rbind, meta_rows)
cells <- cells[order(cells$campaign, cells$prev, cells$hr, cells$n), ]

L <- c(sprintf("# current_status — `%s`", rq), "",
  sprintf("- **Pin:** `%s` on `%s` — HEAD at the time this file was committed; the closeout commit that adds this file is its child.", pin, branch),
  "- **Updated:** 2026-09-13 (`TASK_p12x20_merge_2026-09-13`: `campaign/p12x20` merged into `feature/glm-extension`).",
  "- **Generated by** `scripts_dinamr/current_status_regen.R` from the directory and from git; §4 by `scripts_dinamr/status_inventory.R`. Nothing here is carried from a prior status file, a report or a chat record.",
  "- **Maintenance:** regenerated as the closeout step of a task that touches this directory; the pin must equal HEAD at commit time (`scripts_dinamr/check_current_status.sh`). Campaign-scoped records such as `STATUS_p12x20.md` make no pin claim.",
  "", "---", "",
  "## 1. Where each campaign lives", "",
  "**The layouts differ.** `dinamr` and `grfmr` keep bundles in `results/` and renders in the `gbsg_020/` root. `p12x20` keeps bundles, renders, per-cell gate records and per-cell reports together in `p12x20_2026-09-12/`. Its campaign-scoped record is `STATUS_p12x20.md`.", "",
  "| campaign | engine / focus / ε (combined metas) | combined bundles | batch bundles | bundles in | renders (batch + combine) | scripts | campaign documents in `gbsg_020/` |",
  "|---|---|---|---|---|---|---|---|", nav, "",
  sprintf("- `p12x20` per-cell gate records and reports: %d `GATE2_p12x20_*` and %d `REPORT_p12x20_A*` files in `%s/`.",
          length(list.files(PAYDIR, "^GATE2_p12x20_")), length(list.files(PAYDIR, "^REPORT_p12x20_A")), PAYDIR),
  "- `p12x20` campaign-scoped record: `STATUS_p12x20.md` (cells, gate counts, payload inventory, Stage 0 install record); per-cell heartbeat `LOG_p12x20_progress.txt`; runner logs `scripts_p12x20/logs/`.",
  "- `dinamr` and `grfmr` share `scripts_dinamr/`; `scripts_dinamr/README.md` states each script's role. The script counts above are files whose names match the campaign.",
  "", "## 2. Cells, from the combined bundle metas", "")
for (cn in names(camps)) {
  d <- cells[cells$campaign == cn, ]
  L <- c(L, sprintf("### `%s` — %d cells", cn, nrow(d)), "",
    "| cell | prevalence (super-population) | HR | n | rows | sim_id | batches | host | workers | forestsearch | combined bundle | size | tracked |",
    "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    sprintf("| %s | %.5f | %.2f | %d | %d | %s | %d | %s | %s | %s | `%s` | %s | %s |", d$cell, d$prev, d$hr, d$n, d$rows, d$sim,
            d$batches, d$host, d$workers, d$version, d$path, MB(d$size), ifelse(d$tracked, "yes", "no")), "")
}
seeds <- unique(unlist(lapply(cells$path, function(p) readRDS(p)$meta$seed_base)))
L <- c(L, sprintf("- `seed_base` in every combined meta above: %s (per replicate `seed_base + sim_id`).", paste(seeds, collapse = ", ")), "")

L <- c(L, "## 3. Git", "")
for (cn in names(camps)) {
  cc <- camps[[cn]]
  specs <- switch(cn,
    dinamr = c(file.path(rq, "results/*_dinamr_*"), file.path(rq, "dinamr_*.html")),
    grfmr  = c(file.path(rq, "results/*_grfmr_*"),  file.path(rq, "grfmr_*.html")),
    p12x20 = c(file.path(rq, PAYDIR), file.path(rq, "scripts_p12x20"), file.path(rq, "STATUS_p12x20.md")))
  lg <- git("log", "--format=%h %ad %s", "--date=short", "--", specs)
  L <- c(L, sprintf("- `%s` bundles/renders: %d commits; first `%s`; last `%s`.", cn, length(lg),
                    sub(" .*", "", tail(lg, 1) %||% "none"), sub(" .*", "", lg[1] %||% "none")))
}
cp <- git("rev-parse", "--short", "campaign/p12x20")
anc <- length(cp) && identical(system2("git", c("-C", repo, "merge-base", "--is-ancestor", "campaign/p12x20", "HEAD")), 0L)
L <- c(L, sprintf("- `campaign/p12x20` tip `%s`; ancestor of HEAD: %s.", cp %||% "absent", if (anc) "yes" else "no"), "")

inv <- local({ o <- setwd(here); on.exit(setwd(o)); system2("Rscript", "status_inventory.R", stdout = TRUE) })
L <- c(L, "## 4. Payload inventory (from the directory, `status_inventory.R`)", "", inv, "")

rn <- list.files("results", pattern = "_(res|combined)_[0-9]+_[0-9]+[.]rds$")
tg <- sub("^.*_([^_]+)_(res|combined)_[0-9]+_[0-9]+[.]rds$", "\\1", rn); kd <- sub("^.*_(res|combined)_[0-9]+_[0-9]+[.]rds$", "\\1", rn)
tt <- table(tg, kd)
tags <- rownames(tt)
L <- c(L, "## 5. Campaign tags in `results/` (last token of each bundle stem)", "",
  "| tag | combined | batch |", "|---|---|---|",
  sprintf("| `%s` | %d | %d |", tags, if ("combined" %in% colnames(tt)) tt[, "combined"] else 0L, if ("res" %in% colnames(tt)) tt[, "res"] else 0L), "",
  sprintf("- `%s/`: %d combined, %d batch (tag `p12x20`).", PAYDIR,
          length(list.files(PAYDIR, "_p12x20_combined_")), length(list.files(PAYDIR, "_p12x20_res_"))), "")

prev <- git("log", "-1", "--format=%h", "--", file.path(rq, "current_status.md"))
L <- c(L, "## 6. Not regenerated", "",
  sprintf("- The previous `current_status.md` (last committed at `%s`) carried curated prose that is not derivable from the directory or from git. It is not carried forward here; read it with `git show %s:%s/current_status.md`.",
          prev %||% "none", prev %||% "none", rq), "")
writeLines(L, "current_status.md")
cat("wrote", file.path(qmd, "current_status.md"), "pin", pin, "\n")
