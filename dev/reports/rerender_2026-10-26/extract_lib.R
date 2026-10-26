# Shared extraction for compare.R: code-output lines from quarto / rmarkdown HTML.
suppressPackageStartupMessages(library(xml2))
noise_re <- paste0("fit time|Computation time|elapsed|seconds|minutes|Using .* cores|R version|Platform:|",
                   "Running under|BLAS|LAPACK|time zone|tzcode|attached base|other attached|loaded via|",
                   "^\\[ *[0-9]+\\] [a-zA-Z0-9._]+_[0-9]",
                   # CC addition (sessionInfo cruft not covered by the task's list; recorded as a deviation)
                   "|Matrix products|^locale:|^\\[1\\] [A-Z]+/|^\\[1\\] [a-z][a-z]_[A-Z]+|^\\[1\\] stats +graphics",
                   # CC addition: gc() memory tables, per-iteration timings, machine core counts, absolute payload paths
                   "|^(Ncells|Vcells) |used +\\(Mb\\)|Average per iteration|System max cores|payload to: /|LC_[A-Z]+=")
iso_re <- "[0-9]{4}-[0-9]{2}-[0-9]{2}"
tbl_lines <- function(t) {
  rows <- xml_find_all(t, ".//tr")
  vapply(rows, function(r) paste(trimws(xml_text(xml_find_all(r, "./th|./td"))), collapse = " | "), "")
}
extract_outputs <- function(path) {
  h <- read_html(path)
  quarto <- length(xml_find_all(h, "//div[contains(concat(' ', normalize-space(@class), ' '), ' cell-output')]")) > 0
  if (quarto) {
    nodes <- xml_find_all(h, paste0("//div[contains(concat(' ', normalize-space(@class), ' '), ' cell-output-stdout ') or ",
                                    "contains(concat(' ', normalize-space(@class), ' '), ' cell-output-stderr ')]//pre | ",
                                    "//div[contains(concat(' ', normalize-space(@class), ' '), ' cell-output-display ')]//table[not(ancestor::table)]"))
  } else {
    # V2 uses knitr's default "## " output prefix, not "#>" (CC deviation: both accepted)
    nodes <- xml_find_all(h, "//pre[starts-with(normalize-space(.), '#>') or starts-with(normalize-space(.), '##')] | //div[@class='section level1' or contains(@class,'section')]//table[not(ancestor::table)]")
  }
  out <- character()
  for (n in nodes) {
    if (xml_name(n) == "table") {
      tl <- tbl_lines(n)
      # CC addition: Computational Timing tables are wall-clock noise
      if (length(tl) && grepl("Computational Timing", tl[1])) next
      out <- c(out, "<<table>>", tl)
    }
    else {
      l <- strsplit(xml_text(n), "\n", fixed = TRUE)[[1]]
      out <- c(out, "<<block>>", l)
    }
  }
  out <- trimws(sub("^(#>|##) ?", "", trimws(out)))
  out <- out[nzchar(out)]
  out <- out[!grepl(noise_re, out) & !grepl(iso_re, out)]
  out
}
