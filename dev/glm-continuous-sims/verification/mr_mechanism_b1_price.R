src <- readLines("dev/glm-continuous-sims/verification/mr_mechanism_A1.R")
cut <- grep("^## --- run ---", src); eval(parse(text=paste(src[1:(cut-1)],collapse="\n")))
NB <- 500L      # documented recommendation is 500-1000; <100 warns
sd_i <- SEED_TABLE[1]; RNGkind("L'Ecuyer-CMRG"); set.seed(sd_i)
df <- simulate_from_glm_dgm(dgm, n=1000L, seed=sd_i); df$id <- seq_len(nrow(df))
fs <- suppressWarnings(forestsearch(df.analysis=df, confounders.name=confounders_analysis,
  outcome.name=outcome_name, treat.name=treat_name, id.name=id_name, outcome_type="continuous",
  effect_measure="MD", effect.threshold=md_threshold, consistency.threshold=md_consistency,
  pconsistency.threshold=pconsistency, fs.splits=fs_splits, n.min=n_min, d0.min=d0_min,
  d1.min=d1_min, maxk=maxk, vi.grf.min=vi_grf_min, consistency_method="resample",
  use_twostage=use_twostage, is.RCT=is_rct, adverse_outcome=adverse_outcome,
  details=FALSE, quiet=TRUE, seedit=sd_i, parallel_args=list(plan="sequential"), mr_inference=FALSE))
cat("base fit ok; sg.harm:", paste(fs$sg.harm, collapse=" & "), "\n")
cat("consistency_method in args_call_all:", paste(format(fs$args_call_all$consistency_method), collapse=","), "\n")
t0 <- proc.time()[3]
fb <- tryCatch(forestsearch_bootstrap_dofuture(fs.est=fs, nb_boots=NB, seed=sd_i,
        details=FALSE, parallel_args=list(plan="multisession", workers=100L)),
       error=function(e) structure(list(msg=conditionMessage(e)), class="b1err"))
el <- proc.time()[3]-t0
if (inherits(fb,"b1err")) { cat("\nFB ERRORED:\n"); cat(fb$msg,"\n") } else {
  cat(sprintf("\nFB OK: nb_boots=%d, wall %.1f s\n", NB, el))
  cat("returned names:", paste(head(names(fb),25), collapse=", "), "\n")
}
