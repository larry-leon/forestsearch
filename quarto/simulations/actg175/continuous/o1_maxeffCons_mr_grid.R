# O1 -- maxeffCons MR-only grid over n, sibling run of the committed FB/MR twin.
#
# Document-level only; no R/ changes, no DGM changes.  The twin is the single
# source of configuration: this script PURLS it and evaluates everything up to
# (not including) its run-batch chunk, then overrides exactly four knobs:
#
#   nb_boots      <- 0L        MR-only.  run_fb is !is.null(nb_boots) &&
#                              nb_boots >= 1L, so the bootstrap is never
#                              called -- the skip is structural, not a flag.
#   parallel_mode <- "sims"    replicates fan out; inner_parallel unused.
#   n_sims        <- 1000L     per cell
#   n_sample      <- cell n    500 / 1000 / 2000 / 4000
#
# The J = 10 grids from Step 1 are DELIBERATELY NOT APPLIED here: O1's
# pre-registered question is the focus comparison against the hr grid holding
# everything else equal, so the cut spec must stay as committed.
setwd("/home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/actg175/continuous")
library(forestsearch); library(doFuture); library(foreach)
PURL <- tempfile(fileext=".R")
knitr::purl("sim_fs_maxeffCons_fb_mr_md40_knoise0_n500_batch_1_50.qmd", output=PURL, quiet=TRUE)
src <- readLines(PURL)
i <- grep("^sim_ids <- seq", src)[1]
eval(parse(text=paste(src[1:(i-1)], collapse="\n")), envir=globalenv())

nb_boots      <<- 0L
run_fb        <<- !is.null(nb_boots) && nb_boots >= 1L
parallel_mode <<- "sims"
n_sims        <<- 1000L
stopifnot(identical(run_fb, FALSE))
cat(sprintf("O1: run_fb=%s  parallel_mode=%s  n_sims=%d  mr_draws=%d  workers=%d\n",
            run_fb, parallel_mode, n_sims, mr_draws, n_workers))

outdir <- "o1_maxeffCons_mr_grid"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
n_grid <- c(500L, 1000L, 2000L, 4000L)

for (nn in n_grid) {
  n_sample <<- nn
  cat(sprintf("\n=============== CELL n = %d ===============\n", nn)); flush.console()
  t0 <- proc.time()[3]
  registerDoFuture(); future::plan(future::multisession, workers=n_workers)
  ids <- seq(1L, length.out=n_sims)
  results <- foreach(s = ids, .combine=rbind,
                     .options.future=list(packages=c("forestsearch"), seed=TRUE)) %dofuture%
    tryCatch(record_replicate(s), error=function(e) NULL)
  future::plan(future::sequential); gc(verbose=FALSE)
  el <- proc.time()[3] - t0
  results <- fs_attach_betaHhat(results, eval_df, focus="harm",
                                outcome_type="continuous", effect_measure="MD")
  out <- file.path(outdir, sprintf("o1_maxeffCons_mr_n%d_s%d.rds", nn, n_sims))
  saveRDS(list(results=results, truth=truth,
               meta=list(cell_n=nn, n_sims=n_sims, nb_boots=nb_boots,
                         mr_draws=mr_draws, sg_focus=sg_focus,
                         consistency_method=consistency_method,
                         effect_threshold=md_threshold,
                         consistency_threshold=md_consistency,
                         target_md_harm=target_md_harm,
                         adverse_outcome=adverse_outcome, seed_base=seed_base,
                         parallel_mode=parallel_mode, n_workers=n_workers,
                         elapsed_sec=el,
                         pkg_version=as.character(utils::packageVersion("forestsearch")),
                         built_at=Sys.time())), out)
  cat(sprintf("cell n=%d done in %.1f s (%.2f h); saved %s\n", nn, el, el/3600, out))
  cat(sprintf("  detected %d/%d   mr_ok %d/%d\n",
              sum(results$detected %in% 1L), nrow(results),
              sum(results$mr_ok %in% 1L), nrow(results))); flush.console()
}
cat("\nO1 COMPLETE\n")
