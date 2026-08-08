src <- readLines("dev/glm-continuous-sims/verification/mr_mechanism_A1.R")
cut <- grep("^## --- run ---", src)
eval(parse(text = paste(src[1:(cut-1)], collapse="\n")))
library(doFuture); library(foreach); library(future)
probe <- function(i, n_sample) {
  sd_i <- SEED_TABLE[i]; RNGkind("L'Ecuyer-CMRG"); set.seed(sd_i)
  df <- simulate_from_glm_dgm(dgm, n=n_sample, seed=sd_i); df$id <- seq_len(nrow(df))
  fs <- tryCatch(suppressWarnings(forestsearch(df.analysis=df,
    confounders.name=confounders_analysis, outcome.name=outcome_name,
    treat.name=treat_name, id.name=id_name, outcome_type="continuous",
    effect_measure="MD", effect.threshold=md_threshold,
    consistency.threshold=md_consistency, pconsistency.threshold=pconsistency,
    fs.splits=fs_splits, n.min=n_min, d0.min=d0_min, d1.min=d1_min, maxk=maxk,
    vi.grf.min=vi_grf_min, consistency_method="resample",
    use_twostage=use_twostage, is.RCT=is_rct, adverse_outcome=adverse_outcome,
    details=FALSE, quiet=TRUE, seedit=sd_i,
    parallel_args=list(plan="sequential"), mr_inference=FALSE)), error=function(e) NULL)
  if (is.null(fs) || is.null(fs$sg.harm)) return(NULL)
  g <- fs$grp.consistency
  data.frame(sim_id=i, n=n_sample,
    early_stop = isTRUE(g$early_stop_triggered),
    n_eval = if (!is.null(g$n_candidates_evaluated)) g$n_candidates_evaluated else NA_integer_,
    n_total = if (!is.null(g$n_candidates_total)) g$n_candidates_total else NA_integer_,
    n_passed = if (!is.null(g$n_passed)) g$n_passed else NA_integer_,
    stringsAsFactors=FALSE)
}
plan("multisession", workers=40)
for (cl in list(list(n=1000L,k=40L), list(n=4000L,k=20L))) {
  res <- foreach(i=seq_len(cl$k), .options.future=list(packages="forestsearch", seed=TRUE)) %dofuture% probe(i, cl$n)
  r <- do.call(rbind, Filter(Negate(is.null), res))
  cat(sprintf("\n=== n = %d, %d replicates ===\n", cl$n, nrow(r)))
  cat(sprintf("  early_stop_triggered : %d/%d = %.3f\n", sum(r$early_stop), nrow(r), mean(r$early_stop)))
  cat(sprintf("  candidates evaluated : median %.0f of total median %.0f  (frac %.3f)\n",
      median(r$n_eval, na.rm=TRUE), median(r$n_total, na.rm=TRUE),
      median(r$n_eval/r$n_total, na.rm=TRUE)))
  cat(sprintf("  n_passed             : median %.0f\n", median(r$n_passed, na.rm=TRUE)))
}
plan("sequential")
