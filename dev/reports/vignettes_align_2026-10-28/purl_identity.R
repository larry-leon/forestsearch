# §1.2 gate: purled code of HEAD's vignettes/forestsearch.qmd vs the v4 file
args <- commandArgs(TRUE)  # head.qmd new.qmd outdir
a <- knitr::purl(args[1], output = file.path(args[3], "purl_head.R"), quiet = TRUE)
b <- knitr::purl(args[2], output = file.path(args[3], "purl_new.R"), quiet = TRUE)
same <- identical(readBin(a, "raw", file.size(a)), readBin(b, "raw", file.size(b)))
cat("purl_head md5:", tools::md5sum(a), "\npurl_new  md5:", tools::md5sum(b), "\nbyte-identical:", same, "\n")
if (!same) quit(status = 1)
