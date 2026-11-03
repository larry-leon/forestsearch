## =============================================================================
## dev/article/make_rep5_trial.R
##
## Writes inst/extdata/mrct-consistency/rep5_trial.csv: replicate 5 of
## cell S2 of the mrct-consistency campaign, as drawn, one row per patient.
## (TASK_fs-mrct-article-r3_2026-11-02, section 3.) The article reads the CSV;
## it never reads this script.
##
## The draw is the one mrct_case_study() makes (R/mrct_case_study.R): under
## L'Ecuyer-CMRG, simulate_from_dgm() with the run's sim_args overridden by
## n = 1000, analysis_time = 40, cens_adjust = 50 and seed = seed_S2 + 5, on
## the S2 build of record (latent_R.csv supplied, so no latentcor fit).
## The script checks the draw against mrct_case_study(mr_inference = FALSE)
## and checks that the CSV reads back to the same values exactly.
##
## Run from the package root:  Rscript dev/article/make_rep5_trial.R
## =============================================================================

suppressPackageStartupMessages(library(forestsearch))
dd <- file.path("inst", "extdata", "mrct-consistency")
stopifnot(file.exists(file.path(dd, "PROVENANCE.md")))
rd <- function(f, ...) read.csv(file.path(dd, f), stringsAsFactors = FALSE, check.names = FALSE, ...)

## ---- the S2 build of record (the article's build_cell(h = 2)) --------------
actg <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
actg$treat_mrct <- as.integer(actg$arms == 1L)
actg$time_months <- actg$days / 30.4375
rownames(actg) <- NULL
latent_R <- as.matrix(rd("latent_R.csv", row.names = 1))
dgm_S2 <- suppressWarnings(inject_mrct_structure(
  actg,
  outcome_var = "time_months", event_var = "cens", treatment_var = "treat_mrct",
  continuous_vars = c("age", "wtkg", "karnof", "preanti", "cd40", "cd80"),
  factor_vars = c("hemo", "homo", "drugs", "race", "gender", "symptom"),
  x_pred = "cd40",
  spline_spec = list(knot = 300, zeta = 800, log_hrs = rep(log(0.55), 3)),
  region = list(prevalence = 0.20, or_pred = 1, x1_vars = "symptom", or_x1 = 3,
                loghr = -log(4), band = list(cut = 230, or = 1)),
  x3 = list(vars = "cd40", cuts = list(cd40 = 230), loghr = log(2 / 0.55)),
  region_treat_loghr = NULL,
  n_super = 20000, expand = "copula", seed = 7L,
  continuous_vars_cens = c("age", "cd40"), factor_vars_cens = "symptom",
  latent_R = latent_R))
dgm_S2 <- scale_dgm_hazard(dgm_S2, m = 3)

tf <- rd("timing_full.csv")
seed_S2 <- tf$seed[tf$setting == "S2_rmlog4"]
stopifnot(length(seed_S2) == 1)
k <- 5L
call_seed <- seed_S2 + k

## ---- the draw, as mrct_case_study() makes it --------------------------------
sim_args <- list(rand_ratio = 1, draw_treatment = TRUE, max_entry = 11, cens_adjust = 50)
sim_call <- utils::modifyList(sim_args, list(
  dgm = dgm_S2, n = 1000,
  entry_var = if ("entrytime" %in% names(dgm_S2$df_super)) "entrytime" else NULL,
  analysis_time = 40, cens_adjust = 50, seed = call_seed))
old_kind <- RNGkind()
RNGkind("L'Ecuyer-CMRG")
trial <- as.data.frame(data.table::as.data.table(do.call(simulate_from_dgm, sim_call)))
do.call(RNGkind, as.list(old_kind))

## against mrct_case_study() (identifier only; no MR)
confounders <- analysis_covariates(dgm_S2, "observed")
cs <- mrct_case_study(
  dgm_S2, seed = seed_S2, k = k, n = 1000, region_var = "z_region",
  fs_args = list(sg_focus = "hrMaxSG", maxk = 2,
                 hr.threshold = 0.90, hr.consistency = 0.80, pconsistency.threshold = 0.90,
                 confounders.name = confounders, conf_force = character(0),
                 conf.cont_jcuts = list(z_cd40 = 10, z_wtkg = 10, z_preanti = 10),
                 effect_neighborhood = 0.20, use_grf = FALSE, use_lasso = FALSE),
  sim_args = sim_args, analysis_time = 40, cens_adjust = 50, mr_inference = FALSE)
stopifnot(identical(trial, cs$trial))

## ---- the shipped columns ----------------------------------------------------
cols <- c("id", "y_sim", "event_sim", "treat_sim", "z_region", confounders, "flag_harm")
stopifnot(all(cols %in% names(trial)), !anyDuplicated(cols), !anyDuplicated(trial$id))
out <- trial[, cols]

## shortest decimal that reads back to the same double
fmt <- function(x) {
  if (!is.numeric(x)) return(x)
  s <- sprintf("%.15g", x)
  bad <- as.numeric(s) != x
  s[bad] <- sprintf("%.17g", x[bad])
  s
}
txt <- as.data.frame(lapply(out, fmt), stringsAsFactors = FALSE)
f_out <- file.path(dd, "rep5_trial.csv")
write.csv(txt, f_out, row.names = FALSE, quote = FALSE)

## ---- read-back check --------------------------------------------------------
back <- read.csv(f_out)
stopifnot(identical(names(back), cols), nrow(back) == nrow(out))
for (v in cols) stopifnot(all(as.numeric(back[[v]]) == as.numeric(out[[v]])))

cat(sprintf("wrote %s: %d rows, %d columns; call seed %s (seed_S2 %s + %d); H-hat (identifier, no MR): %s\n",
            f_out, nrow(out), ncol(out), format(call_seed, scientific = FALSE),
            format(seed_S2, scientific = FALSE), k, cs$Hhat))
cat("columns:", paste(cols, collapse = ", "), "\n")
cat("sha256 (12):", substr(digest::digest(file = f_out, algo = "sha256"), 1, 12), "\n")
