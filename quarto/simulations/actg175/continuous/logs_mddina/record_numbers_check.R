# Every number in the record's tables must be found in md_dina_metrics.csv at printed precision
# (value, mc_se, wilson_lo, wilson_hi, n; sprintf formatting as the tables print).
x <- read.csv("md_dina_metrics.csv")
pool <- unique(c(x$value, x$mc_se, x$wilson_lo, x$wilson_hi, x$n)); pool <- pool[is.finite(pool)]
md <- readLines("REPORT_md_dina_2026-09-17.md")
i0 <- grep("^## Tables", md); i1 <- grep("^## The confound", md)
tl <- md[i0:i1]; tl <- tl[grepl("^\\|", tl) & !grepl("^\\|---", tl) & !grepl("^\\|cell\\|", tl)]
cells <- unlist(lapply(strsplit(tl, "|", fixed = TRUE), function(r) r[-(1:2)]))   # drop the cell label column
cells <- sub("^(Bonferroni|separate).*$", "", cells)                               # pair labels carry 95 / 0.025
nums <- unique(unlist(regmatches(cells, gregexpr("-?[0-9]+\\.?[0-9]*", cells))))
fmt <- function(s) { d <- if (grepl("\\.", s)) nchar(sub(".*\\.", "", s)) else 0; s %in% c(sprintf(paste0("%.", d, "f"), pool), sprintf(paste0("%.", d, "f"), -pool)) || as.numeric(s) %in% round(pool, d) }
found <- vapply(nums, fmt, logical(1))
cat("distinct numbers in the record's tables:", length(nums), "| found in the CSV at printed precision:", sum(found), "\n")
if (any(!found)) { cat("NOT FOUND:", nums[!found], "\n"); quit(status = 1) }
