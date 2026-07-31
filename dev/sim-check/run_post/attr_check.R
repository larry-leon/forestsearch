suppressPackageStartupMessages({library(forestsearch);library(survival)
  library(data.table);library(foreach);library(doFuture);library(future)})
source("betaHhat_truth.R")
seed_base <- 8316951L; sim_id <- 1L; n_sample <- 500L
n_workers <- ceiling(0.90*max(1L, min(parallel::detectCores(logical=FALSE)-1L)))
inner_parallel <- list(plan="multisession", workers=n_workers)
k_inter <- calibrate_k_inter(target_hr_harm=1.0, model="alt", use_ahr=FALSE)
dgm <- compute_dgm_cde(setup_gbsg_dgm(model="alt", k_inter=k_inter,
                                      n_super=100000L, seed=seed_base))
df <- simulate_from_dgm(dgm, n=n_sample, analysis_time=84,
                        cens_adjust=log(1.5), seed=seed_base+sim_id)
df$id <- seq_len(nrow(df))
confs <- intersect(c("er","age","meno","pgr","nodes","size","grade"), names(df))
fs.est <- forestsearch(df.analysis=df, outcome.name="y_sim", event.name="event_sim",
  treat.name="treat_sim", id.name="id", flag_harm.name="flag_harm",
  confounders.name=confs, is.RCT=TRUE, seedit=seed_base+sim_id, quiet=TRUE,
  sg_focus="maxcons", subgroup_method="consistency", hr.threshold=0.90,
  hr.consistency=0.80, pconsistency.threshold=0.90, n.min=NULL,
  selection_rule="neighborhood", effect_neighborhood=0.10, stop_threshold=NULL,
  parallel_args=inner_parallel, mr_inference=TRUE,
  mr_inference_args=list(ci_method="ij", draws=500L, include_complement=TRUE),
  consistency_method="resample", use_lasso=FALSE, use_grf=FALSE,
  use_twostage=TRUE, use_dina=FALSE,
  conf_force=c("meno == 0","er <= 0","pgr <= 0"), conf.cont_jcuts=list(er=10),
  fs.splits=400L, maxk=2L, d0.min=10L, d1.min=10L)
bF <- suppressWarnings(forestsearch_bootstrap_dofuture(fs.est, nb_boots=20L,
        seed=seed_base+sim_id, parallel_args=inner_parallel, mr_in_replicates=FALSE))
bT <- suppressWarnings(forestsearch_bootstrap_dofuture(fs.est, nb_boots=20L,
        seed=seed_base+sim_id, parallel_args=inner_parallel, mr_in_replicates=TRUE))
cat("\n### PER-REPLICATE TABLE COMPARISON\n")
x <- bF$results; y <- bT$results
cat("class:", class(x)[1], " dim:", paste(dim(x), collapse="x"), "\n")
cat("attr(FALSE):", paste(names(attributes(x)), collapse=", "), "\n")
cat("attr(TRUE) :", paste(names(attributes(y)), collapse=", "), "\n")
cat("extra attr on TRUE:", paste(setdiff(names(attributes(y)), names(attributes(x))), collapse=", "), "\n")
attributes(y) <- attributes(y)[names(attributes(x))]
cat("per-replicate table identical after dropping that attribute:", identical(x, y), "\n")
cat("columns:", paste(names(x), collapse=", "), "\n")
for (nm in c("SG_CIs","FSsg_tab","Ystar_mat","original_sg","est.scale"))
  cat(sprintf("  %-12s identical: %s\n", nm, identical(bF[[nm]], bT[[nm]])))
plan("sequential")
cat("### DONE2\n")

cat("\n### COLUMN-LEVEL DIFF OF THE PER-REPLICATE TABLE\n")
x <- as.data.frame(bF$results); y <- as.data.frame(bT$results)
for (nm in names(x)) {
  a <- x[[nm]]; b <- y[[nm]]
  if (!identical(a, b)) {
    n <- sum(!((is.na(a)&is.na(b)) | (!is.na(a)&!is.na(b)&a==b)))
    cat(sprintf("  DIFFERS %-16s (%d/%d rows)  max|d| = %s\n", nm, n, length(a),
                if (is.numeric(a)) format(max(abs(a-b), na.rm=TRUE), digits=4) else "n/a"))
  }
}
cat("  -- all other columns identical --\n")
cat("attr 'timing' identical:", identical(attr(bF$results,"timing"), attr(bT$results,"timing")), "\n")
cat("### DONE3\n")
