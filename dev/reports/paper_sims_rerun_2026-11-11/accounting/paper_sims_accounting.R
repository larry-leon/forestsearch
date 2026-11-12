# paper_sims_accounting.R -- Part R step 6 of TASK_classification_outside_package_2026-11-11_v2.
#
# Compares every reported number of two pkgdown builds of vignettes/articles/paper_simulations.qmd:
# the 2026-10-30 record and the 2026-11-11 re-render. "Reported" = the article's visible output:
# prose, console output (<pre> inside a Quarto cell-output div) and tables. Code (code-fold
# <details>, sourceCode blocks) and section numbers in headings are not counted.
#
# Lines of the two pages are aligned by GNU diff on their text with every number masked as "#",
# so a line whose numbers changed still pairs with its counterpart. Table rows carry their
# table's column labels (header grid with colspan/rowspan expanded).
#
# Usage: Rscript paper_sims_accounting.R <old.html> <new.html> <out_dir>
# Writes <out_dir>/accounting_numbers.csv (one row per number; read by paper_sims_groups.R, which
#        writes accounting_grouped.csv with every column of it plus the groups) and
#        <out_dir>/accounting_unpaired_lines.csv (lines with no counterpart).
suppressMessages(library(xml2))
a <- commandArgs(TRUE)
old_f <- a[1]; new_f <- a[2]; out <- a[3]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

squish <- function(x) trimws(gsub("[[:space:] ]+", " ", x))
NUM_RE <- "(?<![A-Za-z0-9_.\\x{0302}\\x{0303}])[-−]?\\d+(?:,\\d{3})*(?:\\.\\d+)?(?:[eE][-+]?\\d+)?"
tokens <- function(s) {
  m <- gregexpr(NUM_RE, s, perl = TRUE)[[1]]
  if (m[1] == -1) return(data.frame(tok = character(0), start = integer(0)))
  data.frame(tok = regmatches(s, list(m))[[1]], start = as.integer(m))
}
mask <- function(s) gsub(NUM_RE, "#", s, perl = TRUE)
as_num <- function(tok) as.numeric(gsub(",", "", sub("−", "-", tok)))

has_anc <- function(n, xp) length(xml_find_first(n, xp)) > 0 && !is.na(xml_find_first(n, xp))

header_labels <- function(tab, ncol) {
  hr <- xml_find_all(tab, "./thead/tr")
  lab <- rep("", ncol); title <- character(0)
  if (!length(hr)) return(list(lab = lab, title = title))
  grid <- matrix(NA_character_, nrow = length(hr), ncol = max(ncol, 1))
  for (r in seq_along(hr)) {
    cells <- xml_find_all(hr[[r]], "./th|./td")
    spans <- as.integer(xml_attr(cells, "colspan")); spans[is.na(spans)] <- 1L
    if (length(cells) == 1L && spans[1] >= ncol && ncol > 1) { title <- c(title, squish(xml_text(cells[[1]]))); next }
    c0 <- 1L
    for (k in seq_along(cells)) {
      while (c0 <= ncol && !is.na(grid[r, c0])) c0 <- c0 + 1L
      rs <- as.integer(xml_attr(cells[[k]], "rowspan")); if (is.na(rs)) rs <- 1L
      rows <- r:min(length(hr), r + rs - 1L); cols <- c0:min(ncol, c0 + spans[k] - 1L)
      if (c0 <= ncol) grid[rows, cols] <- squish(xml_text(cells[[k]]))
      c0 <- c0 + spans[k]
    }
  }
  for (j in seq_len(ncol)) { v <- grid[, j]; v <- v[!is.na(v) & nzchar(v)]; if (length(v)) lab[j] <- v[length(v)] }
  list(lab = lab, title = title)
}

extract <- function(f) {
  h <- read_html(f, encoding = "UTF-8")
  main <- xml_find_first(h, "//main")
  if (is.na(main)) main <- xml_find_first(h, "//body")
  nodes <- xml_find_all(main, paste0(".//*[self::h1 or self::h2 or self::h3 or self::h4 or self::h5 or self::p",
                                     " or self::li or self::table or self::pre or self::dt or self::dd]"))
  L <- list(); block <- "(title block)"; ti <- 0L
  add <- function(kind, text, cells = NULL, rowlab = "", collab = NULL, tab = NA) {
    L[[length(L) + 1L]] <<- list(block = block, kind = kind, text = text, cells = cells,
                                 rowlab = rowlab, collab = collab, tab = tab)
  }
  for (n in nodes) {
    nm <- xml_name(n)
    if (has_anc(n, "ancestor::details") || has_anc(n, "ancestor::*[contains(@class,'sourceCode')]")) next
    if (nm %in% c("h1", "h2", "h3", "h4", "h5")) {
      if (has_anc(n, "ancestor::table")) next
      block <- squish(xml_text(n)); add("heading", block); next
    }
    if (nm == "pre") {
      if (!has_anc(n, "ancestor::div[contains(@class,'cell-output')]")) next
      for (ln in strsplit(xml_text(n), "\n", fixed = TRUE)[[1]]) { ln <- squish(ln); if (nzchar(ln)) add("output", ln) }
      next
    }
    if (nm == "table") {
      if (has_anc(n, "ancestor::table")) next
      ti <- ti + 1L
      body <- xml_find_all(n, "./tbody/tr|./tfoot/tr|./tr")
      ncol <- max(c(0L, vapply(body, function(r) sum(pmax(1L, as.integer(xml_attr(xml_find_all(r, "./th|./td"), "colspan")), na.rm = TRUE)), 1L)))
      hl <- header_labels(n, ncol)
      cap <- xml_find_first(n, "./caption")
      if (!is.na(cap)) add("table-title", squish(xml_text(cap)), tab = ti)
      for (t in hl$title) add("table-title", t, tab = ti)
      for (r in body) {
        cells <- vapply(xml_find_all(r, "./th|./td"), function(cc) squish(xml_text(cc)), "")
        if (!length(cells)) next
        if (length(cells) == 1L) { add("table-span", cells, tab = ti); next }
        add("table-row", paste(cells, collapse = " | "), cells = cells, rowlab = cells[1],
            collab = hl$lab[seq_along(cells)], tab = ti)
      }
      next
    }
    if (has_anc(n, "ancestor::table") || has_anc(n, "ancestor::pre")) next
    if (nm == "p" && has_anc(n, "ancestor::li")) next
    if (nm == "li") {
      parts <- xml_find_all(n, "./text()|./*[not(self::ul or self::ol)]")
      txt <- squish(paste(xml_text(parts), collapse = " "))
    } else txt <- squish(xml_text(n))
    if (nzchar(txt)) add("prose", txt)
  }
  L
}

O <- extract(old_f); N <- extract(new_f)
mk <- function(L) vapply(L, function(x) paste(x$kind, mask(x$text), sep = "\t"), "")
fo <- tempfile(); fn <- tempfile(); writeLines(mk(O), fo); writeLines(mk(N), fn)
ops <- suppressWarnings(system2("diff", c("--unchanged-line-format='=\n'", "--old-line-format='-\n'",
                                          "--new-line-format='+\n'", fo, fn), stdout = TRUE))
# pair: '=' pairs; within a hunk of k '-' and k '+' pair in order (text changed); else unpaired
pairs <- list(); unp <- list(); i <- 0L; j <- 0L; k <- 1L
while (k <= length(ops)) {
  if (ops[k] == "=") { i <- i + 1L; j <- j + 1L; pairs[[length(pairs) + 1L]] <- c(i, j, 0L); k <- k + 1L; next }
  del <- integer(0); ins <- integer(0)
  while (k <= length(ops) && ops[k] != "=") {
    if (ops[k] == "-") { i <- i + 1L; del <- c(del, i) } else { j <- j + 1L; ins <- c(ins, j) }
    k <- k + 1L
  }
  if (length(del) == length(ins)) for (q in seq_along(del)) pairs[[length(pairs) + 1L]] <- c(del[q], ins[q], 1L)
  else {
    for (d in del) unp[[length(unp) + 1L]] <- data.frame(side = "2026-10-30 only", block = O[[d]]$block, kind = O[[d]]$kind, text = O[[d]]$text)
    for (s in ins) unp[[length(unp) + 1L]] <- data.frame(side = "2026-11-11 only", block = N[[s]]$block, kind = N[[s]]$kind, text = N[[s]]$text)
  }
}
stopifnot(i == length(O), j == length(N))

ctx <- function(s, st) squish(substr(s, max(1L, st - 40L), st - 1L))
aft <- function(s, st, tok) substr(s, st + nchar(tok), st + nchar(tok) + 39L)
rows <- list()
for (p in pairs) {
  o <- O[[p[1]]]; nn <- N[[p[2]]]
  if (o$kind == "heading") next
  if (o$kind == "table-row" && length(o$cells) == length(nn$cells)) {
    for (cc in seq_along(o$cells)) {
      to <- tokens(o$cells[cc]); tn <- tokens(nn$cells[cc])
      if (!nrow(to) && !nrow(tn)) next
      same_n <- nrow(to) == nrow(tn)
      for (q in seq_len(max(nrow(to), nrow(tn)))) rows[[length(rows) + 1L]] <- data.frame(
        block = o$block, kind = o$kind, table = o$tab, row_label = o$rowlab,
        column_label = if (cc <= length(o$collab)) o$collab[cc] else "",
        context = if (cc == 1L) "row label" else "", after = "",
        old = if (q <= nrow(to)) to$tok[q] else NA, new = if (q <= nrow(tn)) tn$tok[q] else NA,
        text_changed = p[3] == 1L || !same_n, old_line = o$text, new_line = nn$text)
    }
  } else {
    to <- tokens(o$text); tn <- tokens(nn$text)
    for (q in seq_len(max(nrow(to), nrow(tn)))) rows[[length(rows) + 1L]] <- data.frame(
      block = o$block, kind = o$kind, table = o$tab, row_label = o$rowlab, column_label = "",
      context = if (q <= nrow(to)) ctx(o$text, to$start[q]) else ctx(nn$text, tn$start[q]),
      after = if (q <= nrow(to)) aft(o$text, to$start[q], to$tok[q]) else aft(nn$text, tn$start[q], tn$tok[q]),
      old = if (q <= nrow(to)) to$tok[q] else NA, new = if (q <= nrow(tn)) tn$tok[q] else NA,
      text_changed = p[3] == 1L || nrow(to) != nrow(tn), old_line = o$text, new_line = nn$text)
  }
}
R <- do.call(rbind, rows)
R$old_value <- as_num(R$old); R$new_value <- as_num(R$new)
R$difference <- R$new_value - R$old_value
R$changed <- is.na(R$old) | is.na(R$new) | R$old != R$new
write.csv(R, file.path(out, "accounting_numbers.csv"), row.names = FALSE)
U <- if (length(unp)) do.call(rbind, unp) else data.frame(side = character(0), block = character(0), kind = character(0), text = character(0))
write.csv(U, file.path(out, "accounting_unpaired_lines.csv"), row.names = FALSE)
cat(sprintf("lines: old %d, new %d; paired %d (%d with text changed); unpaired old %d, new %d\n",
            length(O), length(N), length(pairs), sum(vapply(pairs, `[`, 1L, 3) == 1L),
            sum(U$side == "2026-10-30 only"), sum(U$side == "2026-11-11 only")))
cat(sprintf("numbers in paired lines: %d; changed: %d\n", nrow(R), sum(R$changed)))
