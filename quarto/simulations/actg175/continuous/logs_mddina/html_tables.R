# Extract every table (with its caption) from the rendered summary into markdown pipe tables, in order.
suppressPackageStartupMessages(library(xml2))
h <- read_html("summary_continuous_field_mddina.html")
tabs <- xml_find_all(h, "//table")
for (i in seq_along(tabs)) {
  t <- tabs[[i]]
  cap <- xml_text(xml_find_first(t, ".//caption"))
  sec <- xml_attr(xml_find_first(t, "ancestor::section[1]"), "id")
  hd <- xml_text(xml_find_all(t, ".//thead//th"))
  rows <- xml_find_all(t, ".//tbody/tr")
  drop1 <- length(hd) && !nzchar(trimws(hd[1])); if (drop1) hd <- hd[-1]
  cat(sprintf("\n=== TABLE %d [section %s] %s\n", i, sec, gsub("\\s+", " ", cap)))
  cat("|", paste(hd, collapse = "|"), "|\n", sep = "")
  cat("|", paste(rep("---", length(hd)), collapse = "|"), "|\n", sep = "")
  for (r in rows) { v <- gsub("\\|", "/", trimws(xml_text(xml_find_all(r, "./td")))); if (drop1) v <- v[-1]; cat("|", paste(v, collapse = "|"), "|\n", sep = "") }
}
cat("\nEXTRACT LINE:", grep("extract:", xml_text(xml_find_all(h, "//pre|//code")), value = TRUE)[1], "\n")
