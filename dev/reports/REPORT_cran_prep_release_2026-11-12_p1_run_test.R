args <- commandArgs(TRUE); tree <- args[1]; f <- args[2]
suppressMessages(pkgload::load_all(tree, quiet = TRUE))
future::plan("sequential")
t0 <- proc.time()[["elapsed"]]
res <- testthat::test_file(file.path(tree, "tests/testthat", f), reporter = "summary", package = "forestsearch", load_package = "none", stop_on_failure = FALSE)
df <- as.data.frame(res)
print(df[, c("test", "nb", "failed", "skipped", "error", "warning")])
cat("TOTAL expectations", sum(df$nb), "failed", sum(df$failed), "errors", sum(df$error), "\n")
cat("elapsed", round(proc.time()[["elapsed"]] - t0, 2), "s\n")
