# O2 -- FB/MR batches 51 onward on the anchor cell.
# Document-level only; no R/ changes.  Purls the committed twin (now
# parallel_mode = "sims") and overrides only sim_id_start per batch.
# Usage: Rscript o2_fb_mr_batches.R <sim_id_start>
setwd("/home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/actg175/continuous")
args <- commandArgs(trailingOnly=TRUE)
START <- as.integer(args[1])
library(forestsearch); library(doFuture); library(foreach)
PURL <- tempfile(fileext=".R")
knitr::purl("sim_fs_maxeffCons_fb_mr_md40_knoise0_n500_batch_1_50.qmd", output=PURL, quiet=TRUE)
src <- readLines(PURL)
i <- grep("^sim_ids <- seq", src)[1]
eval(parse(text=paste(src[1:(i-1)], collapse="\n")), envir=globalenv())

stopifnot(identical(parallel_mode, "sims"), identical(inner_parallel$plan, "sequential"),
          identical(run_fb, TRUE), nb_boots == 300L, n_sims == 50L)
sim_id_start <<- START
sim_ids <- seq(sim_id_start, length.out = n_sims)
cat(sprintf("O2 batch %d..%d | FB nb=%d | MR draws=%d | mode=%s inner=%s | workers=%d\n",
            min(sim_ids), max(sim_ids), nb_boots, mr_draws, parallel_mode,
            inner_parallel$plan, n_workers)); flush.console()

t0 <- proc.time()[3]
registerDoFuture(); future::plan(future::multisession, workers=n_workers)
results <- foreach(s = sim_ids, .combine=rbind,
                   .options.future=list(packages=c("forestsearch"), seed=TRUE)) %dofuture%
  tryCatch(record_replicate(s), error=function(e) NULL)
future::plan(future::sequential); gc(verbose=FALSE)
el <- proc.time()[3] - t0
results <- fs_attach_betaHhat(results, eval_df, focus="harm",
                              outcome_type="continuous", effect_measure="MD")
out <- file.path(results_dir, sprintf("%s_res_%d_%d.rds", rds_stem, min(sim_ids), max(sim_ids)))
saveRDS(list(results=results, truth=truth,
             meta=list(n_sample=n_sample, n_sims=length(sim_ids), nb_boots=nb_boots,
                       mr_draws=mr_draws, sg_focus=sg_focus, selection_rule=selection_rule,
                       consistency_method=consistency_method, subgroup_method=subgroup_method,
                       target_md_harm=target_md_harm, effect_threshold=md_threshold,
                       consistency_threshold=md_consistency, adverse_outcome=adverse_outcome,
                       seed_base=seed_base, sim_id_start=sim_id_start,
                       parallel_mode=parallel_mode, inner_plan=inner_parallel$plan,
                       n_workers=n_workers, elapsed_sec=el,
                       pkg_version=as.character(utils::packageVersion("forestsearch")),
                       built_at=Sys.time())), out)
cat(sprintf("\nSAVED %s  (%.1f s = %.2f h)\n", out, el, el/3600))
d <- results[results$status=="DETECTED",]
cat(sprintf("detected %d/%d | mr_ok %d/%d | fb_err %d\n",
            sum(results$detected %in% 1L), nrow(results),
            sum(results$mr_ok %in% 1L), nrow(results), sum(!is.na(results$fb_err))))
cat("\nbatch means (detected):\n")
print(colMeans(d[,c("nv_H_est","or_H_est","fb_H_est","mr_H_est","betaHhat_H")], na.rm=TRUE), digits=10)
cat("\nfb_src decomposition:\n")
print(summary(d[,c("fb_src1","fb_src2","fb_nres")]))
# algebra check: fb_src1 + fb_src2 must reconcile with H0 - H2 (raw, oriented back)
chk <- (d$nv_H_est - d$fb_H_est) - (-(d$fb_src1 + d$fb_src2))
cat(sprintf("\nALGEBRA CHECK max|(nv-fb) + (src1+src2)| = %.3e\n", max(abs(chk), na.rm=TRUE)))
cat("\nfb_secs:\n"); print(summary(results$fb_secs))
cat("\nfit_mr_secs:\n"); print(summary(results$fit_mr_secs))
