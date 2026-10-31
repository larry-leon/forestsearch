# vignette_gate.R -- output-only gate on the built vignette's "Where to go next" (part B, 2026-10-30).
#
#   Rscript vignette_gate.R <forestsearch.html from the built tarball's inst/doc> <out.txt>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, file.exists(args[1]))
doc <- xml2::read_html(args[1])
xml2::xml_remove(xml2::xml_find_all(doc, "//script | //style"))
txt <- gsub("\\s+", " ", xml2::xml_text(doc))
count <- function(s) { m <- gregexpr(s, txt, fixed = TRUE)[[1]]; if (m[1] == -1L) 0L else length(m) }
want <- list(
  c("The Guo and He article sets the bounds beside", 1),
  c("The ACTG 175 analyses are an article on the package website", 1),
  c("Table S3)", 1),
  c("with examples on this trial and on the binary and continuous ACTG 175 endpoints", 1),
  c("examples for this analysis", 0))
hrefs <- xml2::xml_attr(xml2::xml_find_all(doc, "//a"), "href")
ok <- TRUE; out <- character()
for (w in want) {
  n <- count(w[1]); pass <- n == as.integer(w[2]); ok <- ok && pass
  out <- c(out, sprintf("%-82s %s  found %d", sprintf("count \"%s\" == %s", w[1], w[2]),
                        if (pass) "PASS" else "FAIL", n))
}
for (h in c("https://larry-leon.github.io/forestsearch/articles/guohe-nested-continuum.html",
            "https://larry-leon.github.io/forestsearch/articles/")) {
  n <- sum(hrefs == h, na.rm = TRUE); pass <- n == 1L; ok <- ok && pass
  out <- c(out, sprintf("%-82s %s  found %d", sprintf("link %s == 1", h), if (pass) "PASS" else "FAIL", n))
}
out <- c(out, sprintf("RESULT: %s", if (ok) "OK" else "FAIL"))
writeLines(out, args[2]); cat(out, sep = "\n")
quit(status = if (ok) 0L else 1L)
