# §2.5 (recorded, not gated): FS bootstrap at B = 25, sequential vs multisession(2).
# usage: Rscript plan_invariance.R <purled vignette script>
ex <- parse(commandArgs(TRUE)[1])
env <- globalenv()
for (e in ex) { eval(e, env); if (exists("fit_fs", envir = env, inherits = FALSE)) break }
stopifnot(exists("fit_fs"), !exists("fit_dina"))
cat("forestsearch", as.character(packageVersion("forestsearch")), "\n")
run <- function(pa) {
  t0 <- proc.time()[["elapsed"]]
  b <- forestsearch_bootstrap_dofuture(fs.est = fit_fs, nb_boots = 25L, seed = 8316951L,
                                       mr_in_replicates = FALSE, parallel_args = pa)
  list(b = b, secs = proc.time()[["elapsed"]] - t0)
}
seq_r <- run(list(plan = "sequential"))
par_r <- run(list(plan = "multisession", workers = 2L))
cat(sprintf("wall-clock: sequential %.1f s; multisession(2) %.1f s\n", seq_r$secs, par_r$secs))
for (k in c("H_estimates", "Hc_estimates"))
  cat(sprintf("%s$H2: sequential %.10f; multisession %.10f; identical %s\n", k,
              seq_r$b[[k]]$H2, par_r$b[[k]]$H2, identical(seq_r$b[[k]]$H2, par_r$b[[k]]$H2)))
rs <- seq_r$b$results; rp <- par_r$b$results
cat("results: class", class(rs)[1], "dims", paste(dim(rs), collapse = "x"), "vs",
    paste(dim(rp), collapse = "x"), "\n")
cat("results identical():", identical(rs, rp), "\n")
cat("results all.equal():", paste(isTRUE(all.equal(as.data.frame(rs), as.data.frame(rp))),
    collapse = " "), "\n")
if (!isTRUE(all.equal(as.data.frame(rs), as.data.frame(rp))))
  print(all.equal(as.data.frame(rs), as.data.frame(rp)))
# Timing-free comparison: drop tmins_* columns and attributes.
strip <- function(x) { x <- as.data.frame(x); x <- x[, !grepl("^tmins", names(x)), drop = FALSE]
  attributes(x) <- attributes(x)[c("names", "row.names", "class")]; x }
cat("results without tmins_* columns and attributes: identical", identical(strip(rs), strip(rp)),
    "; columns compared", ncol(strip(rs)), "\n")
