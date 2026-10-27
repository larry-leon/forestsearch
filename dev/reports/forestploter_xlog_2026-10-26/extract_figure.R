# usage: Rscript extract_figure.R <html> <out.png>
# Extracts the base64 <img> whose figure caption is the forest-plot chunk's fig.cap.
suppressPackageStartupMessages(library(xml2))
a <- commandArgs(TRUE)
cap <- "Subgroup forest plot including identified subgroups"
h <- read_html(a[1])
img <- xml_find_all(h, sprintf("//*[contains(normalize-space(.), '%s') and (@class='caption' or self::figcaption)]/preceding-sibling::img | //img[@alt='%s']", cap, cap))
stopifnot(length(img) >= 1)
src <- xml_attr(img[[1]], "src")
stopifnot(startsWith(src, "data:image/png;base64,"))
writeBin(jsonlite::base64_dec(sub("^data:image/png;base64,", "", src)), a[2])
cat(a[2], file.size(a[2]), "bytes; matches:", length(img), "\n")
