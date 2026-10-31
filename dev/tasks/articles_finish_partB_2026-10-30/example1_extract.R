# example1_extract.R -- write ?forestsearch Example 1, as man/forestsearch.Rd ships it, to a script
# that runs it and records what mr_estimates_table() returns (part B, 2026-10-30).
#
#   Rscript example1_extract.R man/forestsearch.Rd <record dir>
#
# Writes <record dir>/forestsearch_examples.R (every example, \dontrun{} uncommented) and
# <record dir>/example1.R (Example 1 alone, plus the recording lines below).
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, file.exists(args[1]), dir.exists(args[2]))
all_ex <- file.path(args[2], "forestsearch_examples.R")
tools::Rd2ex(args[1], all_ex, commentDontrun = FALSE)
ex <- readLines(all_ex, encoding = "UTF-8")
i1 <- grep("^## Example 1 --", ex)
i2 <- grep("^## Example 2 --", ex)
stopifnot(length(i1) == 1L, length(i2) == 1L, i2 > i1)
body <- c("library(forestsearch)   # R CMD check attaches the package; a script must", ex[i1:(i2 - 1L)])
stopifnot(sum(grepl("mr_estimates_table(fs_surv, boot_results = fs_boot", body, fixed = TRUE)) == 1L,
          !any(grepl("mr_estimates_table(fs_surv)", body, fixed = TRUE)))
record <- c(
  "",
  "## ---- recording (appended by example1_extract.R; not part of the example) ----",
  "cat(\"REGION\", paste(fs_surv$sg.harm, collapse = \" & \"), \"\\n\")",
  "g <- fs_surv$mr_inference",
  "cat(sprintf(\"MR corrected %.3f (%.3f, %.3f); complement %.3f (%.3f, %.3f)\\n\",",
  "            g$debiased$est, g$debiased$lower, g$debiased$upper,",
  "            g$complement$debiased$est, g$complement$debiased$lower, g$complement$debiased$upper))",
  "cat(sprintf(\"BOOT nb_boots %d; H2 %.3f\\n\", 20L, fs_boot$H_estimates$H2))",
  "tab <- mr_estimates_table(fs_surv, boot_results = fs_boot, est.scale = \"hr\")",
  "stopifnot(inherits(tab, \"gt_tbl\"))",
  sprintf("gt::gtsave(tab, %s)", deparse(file.path(normalizePath(args[2]), "example1_mr_table.html"))),
  "cat(\"EXAMPLE1_OK\\n\")")
writeLines(c(body, record), file.path(args[2], "example1.R"), useBytes = TRUE)
cat("wrote", file.path(args[2], "example1.R"), "(", length(body), "example lines )\n")
