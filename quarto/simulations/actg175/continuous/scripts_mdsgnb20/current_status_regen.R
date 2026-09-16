# Regenerates quarto/simulations/actg175/continuous/current_status.md.
# Transplant of quarto/simulations/gbsg_020/scripts_dinamr/current_status_regen.R
# (TASK_md_field_rerun_2026-09-15 §3.4).  Same two sources, never mixed:
#   * GENERATED -- the pin, the campaign layout, cells, git and the payload inventory,
#     from bundle metas, file listings and git; never from a report or a chat record.
#   * CURATED -- status_curated.md, hand-maintained and included VERBATIM; its four
#     marked blocks are placed exactly as gbsg_020's: preamble after the pin, §1-§2
#     before §3, the table commentary after the inventory table, §4-§7 after §3.
#     A missing file or marker is a hard error.
# Named changes: the campaigns are mdf1 and mdsgnb20; a cell is (MD target or null, n)
# rather than (prevalence, HR, n); bundles live in mr_md_harm/<stem>_d5000/; the
# inventory table comes from scripts_mdsgnb20/status_inventory.R.
# The pin is HEAD at generation.  Commit the file ALONE, as HEAD's child, then run
#   bash scripts_mdsgnb20/check_current_status.sh --commit
# usage (from scripts_mdsgnb20/): Rscript current_status_regen.R
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
qmd  <- normalizePath(file.path(here, ".."))
repo <- system2("git", c("-C", qmd, "rev-parse", "--show-toplevel"), stdout = TRUE)
rq   <- sub(paste0("^", repo, "/"), "", qmd)
git  <- function(...) suppressWarnings(system2("git", shQuote(c("-C", repo, ...)), stdout = TRUE, stderr = FALSE))
`%||%` <- function(a, b) if (is.null(a) || !length(a) || all(is.na(a))) b else a
old <- setwd(qmd); on.exit(setwd(old))

# ---- curated source: fail loudly --------------------------------------------
cur_f <- file.path(qmd, "status_curated.md")
if (!file.exists(cur_f)) stop("status_curated.md not found at ", cur_f, " -- refusing to write current_status.md without the curated sections.", call. = FALSE)
cur   <- readLines(cur_f, encoding = "UTF-8", warn = TRUE)
marks <- c("preamble", "before-inventory", "after-inventory-table", "after-inventory")
pos   <- match(sprintf("<!-- curated:%s -->", marks), cur)
if (anyNA(pos) || is.unsorted(pos)) stop("status_curated.md: markers missing or out of order (", paste(marks[is.na(pos)], collapse = ", "), ")", call. = FALSE)
C <- lapply(seq_along(pos), function(i) { a <- pos[i] + 1L; b <- if (i < length(pos)) pos[i + 1L] - 1L else length(cur); if (b < a) character(0) else cur[a:b] })
names(C) <- marks
if (any(lengths(C) == 0L)) stop("status_curated.md: empty block(s): ", paste(marks[lengths(C) == 0L], collapse = ", "), call. = FALSE)

pin    <- git("rev-parse", "--short", "HEAD")
branch <- git("rev-parse", "--abbrev-ref", "HEAD")
trk_all <- git("ls-files", "--", rq)
trk <- function(p) file.path(rq, p) %in% trk_all
MB  <- function(b) sprintf("%.2f MB", b / 1024^2)
tick <- function(x) paste0("`", x, "`")

camps <- list(
  mdf1     = list(dre = "_mdf1_d5000$",     rre = "^fs_.*_mdf1_combine_[0-9]+_[0-9]+[.]html$",     scripts = "scripts_mdf1",     docre = "^(REPORT_continuous_field_|summary_continuous_field_mdf1)"),
  mdsgnb20 = list(dre = "_mdsgnb20_d5000$", rre = "^fs_.*_mdsgnb20_combine_[0-9]+_[0-9]+[.]html$", scripts = "scripts_mdsgnb20", docre = "^(REPORT_md_field_rerun_|summary_continuous_field_mdsgnb20|LOG_mdsgnb20|HALT_mdsgnb20|md_field_metrics|COLUMNS_md_field)"))
meta_rows <- list(); nav <- character(0)
for (cn in names(camps)) {
  cc <- camps[[cn]]
  dirs <- list.files("mr_md_harm", pattern = cc$dre, full.names = TRUE)
  comb <- unlist(lapply(dirs, list.files, pattern = "_combined_[0-9]+_[0-9]+[.]rds$", full.names = TRUE))
  bat  <- unlist(lapply(dirs, list.files, pattern = "_res_[0-9]+_[0-9]+[.]rds$", full.names = TRUE))
  rf <- list.files(".", pattern = cc$rre); rsz <- sum(file.size(rf))
  sf <- list.files(cc$scripts); sf <- sf[!dir.exists(file.path(cc$scripts, sf))]
  docs <- list.files(".", pattern = cc$docre)
  ms <- lapply(comb, function(f) readRDS(f)$meta)
  eng <- paste(unique(vapply(ms, function(m) sprintf("%s / %s / eps %s / rule %s / scale-complement %s", m$subgroup_method, m$sg_focus, format(m$effect_neighborhood %||% NA), m$selection_rule %||% NA, m$field_scale_complement %||% "not recorded"), "")), collapse = "; ")
  nav <- c(nav, sprintf("| `%s` | %s | %d | %d | `mr_md_harm/<stem>%s/` | %d, %s, in the directory root | `%s/` (%d files) | %s |",
    cn, eng, length(comb), length(bat), sub("[$]$", "", cc$dre), length(rf), MB(rsz), cc$scripts, length(sf), paste(tick(docs), collapse = ", ")))
  for (f in comb) {
    b <- readRDS(f); m <- b$meta
    stem <- sub("_combined_[0-9]+_[0-9]+[.]rds$", "", basename(f))
    bm <- lapply(list.files(dirname(f), pattern = paste0("^", stem, "_res_"), full.names = TRUE), function(x) readRDS(x)$meta)
    meta_rows[[length(meta_rows) + 1]] <- data.frame(
      campaign = cn, cell = if (isTRUE(m$null_cell)) "null" else sprintf("MD %g", abs(m$target_md_harm)), n = m$n_sample,
      rows = nrow(b$results), sim = paste(range(b$results$sim_id), collapse = "-"), batches = length(bm),
      declared = sum(b$results$detected %in% 1L),
      host = paste(unique(vapply(bm, function(x) x$hostname %||% NA_character_, "")), collapse = "/"),
      workers = paste(unique(vapply(bm, function(x) as.character(x$n_workers %||% NA), "")), collapse = "/"),
      version = paste(unique(vapply(bm, function(x) x$pkg_version %||% NA_character_, "")), collapse = "/"),
      path = f, size = file.size(f), tracked = trk(f), stringsAsFactors = FALSE)
  }
}
cells <- do.call(rbind, meta_rows)

L <- c(sprintf("# current_status — `%s`", rq), "",
  sprintf("- **Pin:** `%s` on `%s` — HEAD at the time this file was committed; the closeout commit that adds this file is its child.", pin, branch),
  "- **Generated by** `scripts_mdsgnb20/current_status_regen.R`: the pin and §3 from the directory and from git (§3.5's table by `scripts_mdsgnb20/status_inventory.R`). Everything else is included verbatim from `status_curated.md`.",
  C[["preamble"]], "", "---", "", C[["before-inventory"]])
L <- c(L, "## 3. Payload inventory", "",
  "Regenerated from the directory and from git by `scripts_mdsgnb20/current_status_regen.R`; the table in 3.5 by `scripts_mdsgnb20/status_inventory.R`. The commentary after that table, and §1, §2, §4–§7, are curated (`status_curated.md`).", "",
  "### 3.1 Where each campaign lives", "",
  "Both campaigns keep their bundles under `mr_md_harm/<stem>_d5000/` (one directory per cell: two batch bundles and the combined bundle) and their combine renders in the directory root.", "",
  "| campaign | engine / focus / ε / rule / complement scale (combined metas) | combined bundles | batch bundles | bundles in | renders (combine) | scripts | campaign documents |",
  "|---|---|---|---|---|---|---|---|", nav, "",
  "### 3.2 Cells, from the combined bundle metas", "")
for (cn in names(camps)) {
  d <- cells[cells$campaign == cn, ]
  L <- c(L, sprintf("#### `%s` — %d cells", cn, nrow(d)), "",
    "| cell | n | rows | sim_id | batches | declared | host | workers | forestsearch | combined bundle | size | tracked |", "|---|---|---|---|---|---|---|---|---|---|---|---|",
    sprintf("| %s | %d | %d | %s | %d | %d | %s | %s | %s | `%s` | %s | %s |", d$cell, d$n, d$rows, d$sim, d$batches, d$declared, d$host, d$workers, d$version, d$path, MB(d$size), ifelse(d$tracked, "yes", "no")), "")
}
seeds <- unique(unlist(lapply(cells$path, function(p) readRDS(p)$meta$seed_base)))
L <- c(L, sprintf("- `seed_base` in every combined meta above: %s (per replicate `seed_base + sim_id`; the two campaigns share every draw cell for cell).", paste(seeds, collapse = ", ")), "")
L <- c(L, "### 3.3 Git", "")
for (cn in names(camps)) {
  specs <- c(file.path(rq, sprintf("mr_md_harm/*_%s_d5000", cn)), file.path(rq, sprintf("fs_*_%s_combine_*.html", cn)))
  lg <- git("log", "--format=%h %ad %s", "--date=short", "--", specs)
  L <- c(L, sprintf("- `%s` bundles/renders: %d commits; first `%s`; last `%s`.", cn, length(lg), sub(" .*", "", tail(lg, 1) %||% "none"), sub(" .*", "", lg[1] %||% "none")))
}
L <- c(L, "")
dn <- list.files("mr_md_harm"); tg <- sub("_d[0-9]+$", "", dn); tg <- sub("^.*_", "", tg); tt <- table(tg)
L <- c(L, "### 3.4 Campaign tags under `mr_md_harm/` (last token of each bundle directory's stem)", "", "| tag | cell directories |", "|---|---|",
  sprintf("| `%s` | %d |", names(tt), as.integer(tt)), "")
inv <- local({ o <- setwd(here); on.exit(setwd(o)); system2("Rscript", "status_inventory.R", stdout = TRUE) })
L <- c(L, "### 3.5 Payload inventory table", "", inv, "", C[["after-inventory-table"]], "", C[["after-inventory"]])
writeLines(L, "current_status.md")
cat("wrote", file.path(qmd, "current_status.md"), "pin", pin, "\n")
