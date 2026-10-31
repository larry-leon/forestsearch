# actg175_article_gate.R -- output-only gate on a rendered ACTG 175 article (part B, 2026-10-30).
#
#   Rscript actg175_article_gate.R <rendered.html> <out.txt>
#
# Reads the seven tables the article prints (in document order: Table S16 run and published,
# Table S17 run and published, Table S3 run and published, the cost table) and a set of fixed
# strings from the article's prose and printed output, compares them with the published
# values, writes one line per check to <out.txt>, and exits 1 if any check fails.
# Nothing here re-computes an analysis; it reads what the article printed.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, file.exists(args[1]))
doc <- xml2::read_html(args[1])
xml2::xml_remove(xml2::xml_find_all(doc, "//script | //style"))

lines <- character()
ok <- TRUE
chk <- function(name, cond, detail = "") {
  pass <- isTRUE(cond)
  lines <<- c(lines, sprintf("%-58s %s  %s", name, if (pass) "PASS" else "FAIL", detail))
  if (!pass) ok <<- FALSE
  invisible(pass)
}

# ---- the tables ---------------------------------------------------------------------------
as_matrix <- function(tab) {
  rows <- xml2::xml_find_all(tab, ".//tr")
  cells <- lapply(rows, function(r) {
    trimws(gsub("\\s+", " ", xml2::xml_text(xml2::xml_find_all(r, "./th | ./td"))))
  })
  widths <- lengths(cells)
  if (length(unique(widths)) != 1L) return(NULL)
  do.call(rbind, cells)
}
tabs <- lapply(xml2::xml_find_all(doc, "//table"), as_matrix)
chk("seven tables, each rectangular", length(tabs) == 7L && !any(vapply(tabs, is.null, TRUE)),
    sprintf("found %d", length(tabs)))
if (length(tabs) != 7L || any(vapply(tabs, is.null, TRUE))) {
  writeLines(lines, args[2]); quit(status = 1)
}
s16 <- tabs[[1]]; s16p <- tabs[[2]]; s17 <- tabs[[3]]; s17p <- tabs[[4]]
s3 <- tabs[[5]]; s3p <- tabs[[6]]; cost <- tabs[[7]]

show_diff <- function(a, b, keep) {
  d <- which(a != b, arr.ind = TRUE)
  d <- d[keep(d[, "row"]), , drop = FALSE]
  list(n = nrow(d),
       txt = paste(apply(d, 1, function(i) sprintf("[%s | col %d] run '%s' vs published '%s'",
                                                   a[i[1], 1], i[2], a[i[1], i[2]], b[i[1], i[2]])),
                   collapse = "; "))
}

# Table S16: 1 header row + 14 rows x 3 columns; every cell but the Region row equal.
chk("S16 shape 15 x 3 (run and published)",
    identical(dim(s16), c(15L, 3L)) && identical(dim(s16p), c(15L, 3L)))
if (identical(dim(s16), dim(s16p))) {
  ri <- which(s16[, 1] == "Region")
  d16 <- show_diff(s16, s16p, function(r) r != ri)
  chk("S16 every cell but the Region row", length(ri) == 1L && d16$n == 0L, d16$txt)
  chk("S16 Region, binary: the package's label",
      grepl("cd40 <= 380", s16[ri, 2], fixed = TRUE) && grepl("wtkg <= 86", s16[ri, 2], fixed = TRUE),
      s16[ri, 2])
  chk("S16 Region, continuous: the package's label",
      grepl("cd40 <= 507", s16[ri, 3], fixed = TRUE) && grepl("gender", s16[ri, 3], fixed = TRUE),
      s16[ri, 3])
}

# Table S17: 1 header row + 9 rows x 3 columns; identical.
chk("S17 shape 10 x 3 (run and published)",
    identical(dim(s17), c(10L, 3L)) && identical(dim(s17p), c(10L, 3L)))
if (identical(dim(s17), dim(s17p))) {
  d17 <- show_diff(s17, s17p, function(r) rep(TRUE, length(r)))
  chk("S17 every cell", d17$n == 0L, d17$txt)
}

# Table S3: 1 header row + 6 rows x 6 columns.  Body only (the headers differ by design).
chk("S3 shape 7 x 6 (run and published)",
    identical(dim(s3), c(7L, 6L)) && identical(dim(s3p), c(7L, 6L)))
if (identical(dim(s3), c(7L, 6L)) && identical(dim(s3p), c(7L, 6L))) {
  b <- s3[-1, , drop = FALSE]; bp <- s3p[-1, , drop = FALSE]
  chk("S3 Method, SG, n (%), Naive", identical(b[, 1:4], bp[, 1:4]),
      paste(which(b[, 1:4] != bp[, 1:4]), collapse = ","))
  fs_h <- which(b[, 1] == "FS" & b[, 2] == "H")
  other <- setdiff(seq_len(nrow(b)), fs_h)
  chk("S3 Multiplier, every row but FS H", length(fs_h) == 1L && identical(b[other, 6], bp[other, 6]),
      paste(sprintf("%s %s: run '%s' vs '%s'", b[other, 1], b[other, 2], b[other, 6], bp[other, 6])[
        b[other, 6] != bp[other, 6]], collapse = "; "))
  chk("S3 Multiplier, FS H: 1.79 (0.63, 5.09) under 0.4.1 (Gate C, 10-30)",
      length(fs_h) == 1L && identical(b[fs_h, 6], "1.79 (0.63, 5.09)") &&
        identical(bp[fs_h, 6], "1.79 (0.63, 5.08)"),
      sprintf("run '%s'", b[fs_h, 6]))
  bs <- suppressWarnings(as.numeric(b[, 5]))
  chk("S3 bootstrap column: six finite positive numbers", length(bs) == 6L && all(is.finite(bs) & bs > 0),
      paste(b[, 5], collapse = " "))
}

# Cost table: 1 header row + 3 rows x 6 columns, numeric body.
chk("cost shape 4 x 6", identical(dim(cost), c(4L, 6L)))
if (identical(dim(cost), c(4L, 6L))) {
  num <- suppressWarnings(matrix(as.numeric(gsub(",", "", cost[-1, -1])), nrow = 3L))
  chk("cost table: methods FS, DINA, GRF; numbers finite",
      identical(cost[-1, 1], c("FS", "DINA", "GRF")) && all(is.finite(num)),
      paste(apply(cost[-1, , drop = FALSE], 1, paste, collapse = " "), collapse = " | "))
}

# ---- fixed strings in prose and printed output ----------------------------------------------
txt <- gsub("\\s+", " ", xml2::xml_text(doc))
count <- function(s) {
  m <- gregexpr(s, txt, fixed = TRUE)[[1]]
  if (m[1] == -1L) 0L else length(m)
}
want <- list(
  c("FS and GRF classify 1060 of 1083 patients alike", 1),
  c("brings it to 1.773", 1),
  c("field lower bound of 0.686", 1),
  c("field-s upper bound of 0.786", 1),
  c("leaves 22.1 on the adverse side", 1),
  c("field lower bound of -38.3", 1),
  c("field-s upper bound of -17.6 cells", 1),
  c("harm subgroup in 0.838 (binary) and 0.909 (continuous)", 1),
  c("The implied levels, 0.996 and 0.997", 1),
  c("interval ends at 5.09 where the table prints 5.08", 1),
  c("1.79 by multiplier resampling", 1),
  c("about half the naive 3.58", 1),
  c("NaN", 0),
  c("?@", 0))
for (w in want) {
  n <- count(w[1])
  chk(sprintf("count \"%s\" == %s", w[1], w[2]), n == as.integer(w[2]), sprintf("found %d", n))
}
lines <- c(lines, sprintf("%-58s %s", "recorded (not gated): count \"Warning\"", count("Warning")))

lines <- c(lines, sprintf("RESULT: %s", if (ok) "OK" else "FAIL"))
writeLines(lines, args[2])
cat(lines, sep = "\n")
quit(status = if (ok) 0L else 1L)
