# Spelling and link checks -- TASK_cran_prep_release_2026-11-12, Part 3 step 8.
# Usage: Rscript --vanilla spell_links.R <package tree> <out_prefix>
# Spelling: the `spelling` package is not installed, so base R's utils::aspell()
# family is used with the installed GNU Aspell (en_US) on the Rd files, the
# DESCRIPTION and the CRAN vignette.  Links: urlchecker::url_check() (installed,
# 2.0.0) on the package tree; it reads each URL over the network.
a <- commandArgs(TRUE); tree <- a[1]; out <- a[2]
sink_txt <- paste0(out, "_spelling.txt")
words_of <- function(x) if (is.null(x) || !nrow(x)) character(0) else
  sort(unique(x$Original))
res <- list()
res$Rd <- tryCatch(utils::aspell_package_Rd_files(tree, dictionaries = "en_stats",
                                                  program = "aspell"),
                   error = function(e) e)
res$DESCRIPTION <- tryCatch(utils:::aspell_package_description(tree, program = "aspell"),
                            error = function(e) e)
vig <- file.path(tree, "vignettes", "forestsearch.qmd")
res$vignette <- tryCatch(utils::aspell(vig, filter = "md", program = "aspell"),
                         error = function(e) e)
res$NEWS <- tryCatch(utils::aspell(file.path(tree, "NEWS.md"), filter = "md",
                                   program = "aspell"), error = function(e) e)
con <- file(sink_txt, "w")
for (nm in names(res)) {
  r <- res[[nm]]
  writeLines(sprintf("== %s", nm), con)
  if (inherits(r, "error")) { writeLines(paste("ERROR:", conditionMessage(r)), con); next }
  writeLines(sprintf("flagged occurrences: %d; distinct words: %d", nrow(r), length(words_of(r))), con)
  if (nrow(r)) {
    agg <- stats::aggregate(File ~ Original, data = transform(r, File = basename(File)),
                            FUN = function(f) paste(unique(f), collapse = ", "))
    writeLines(sprintf("  %-28s %s", agg$Original, agg$File), con)
  }
}
close(con)
cat(readLines(sink_txt), sep = "\n")

if (requireNamespace("urlchecker", quietly = TRUE)) {
  t0 <- proc.time()[["elapsed"]]
  u <- tryCatch(urlchecker::url_check(tree, progress = FALSE), error = function(e) e)
  cat("\n== urlchecker::url_check (", round(proc.time()[["elapsed"]] - t0, 1), "s )\n")
  if (inherits(u, "error")) cat("ERROR:", conditionMessage(u), "\n") else {
    print(u)
    saveRDS(u, paste0(out, "_urls.rds"))
    cat("problems:", nrow(u), "\n")
  }
}
