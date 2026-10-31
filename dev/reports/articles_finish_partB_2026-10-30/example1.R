library(forestsearch)   # R CMD check attaches the package; a script must
## Example 1 -- Survival (GBSG, hormone therapy): the analysis of
## León and Anderson (2026a, Sec. 6) and of the survival-post-selection
## article on the package website. Forest search over the cut grid (no
## model-based front ends), effMaxSG selection in a 10% effect
## neighbourhood, multiplier-resampling inference for the reported
## subgroup computed inside the fit.
library(survival)
df_gbsg <- within(gbsg, {
  id          <- seq_len(nrow(gbsg))
  time_months <- rfstime / 30.4375
  grade3      <- ifelse(grade == "3", 1, 0)
})
fs_surv <- forestsearch(
  df.analysis            = df_gbsg,
  outcome.name           = "time_months",
  event.name             = "status",
  treat.name             = "hormon",
  id.name                = "id",
  confounders.name       = c("age", "meno", "size", "grade3", "nodes",
                             "pgr", "er"),
  is.RCT                 = TRUE,
  use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
  conf.cont_jcuts        = list(er = 10, pgr = 10),
  conf_force             = c("er <= 0", "pgr <= 0"),
  subgroup_method        = "consistency",
  sg_focus               = "effMaxSG",
  selection_rule         = "neighborhood",
  effect_neighborhood    = 0.10,
  hr.threshold           = 1.0,
  hr.consistency         = 1.0,
  pconsistency.threshold = 0.90,
  consistency_method     = "resample",
  fs.splits              = 1000,
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
  mr_inference           = TRUE,
  mr_inference_args      = list(draws = 5000L),
  parallel_args          = list(plan = "multisession", workers = 4),
  seedit                 = 8316951
)
print(fs_surv)

## mr_estimates_table() lays the multiplier-resampling estimates out in the
## format of the bootstrap tables, so it takes a bootstrap result as its
## template (a few replicates are enough for the layout; the
## survival-post-selection article compares the two at full size).
fs_boot <- forestsearch_bootstrap_dofuture(
  fs.est        = fs_surv,
  nb_boots      = 20,
  parallel_args = list(plan = "multisession", workers = 4)
)
mr_estimates_table(fs_surv, boot_results = fs_boot, est.scale = "hr")


## ---- recording (appended by example1_extract.R; not part of the example) ----
cat("REGION", paste(fs_surv$sg.harm, collapse = " & "), "\n")
g <- fs_surv$mr_inference
cat(sprintf("MR corrected %.3f (%.3f, %.3f); complement %.3f (%.3f, %.3f)\n",
            g$debiased$est, g$debiased$lower, g$debiased$upper,
            g$complement$debiased$est, g$complement$debiased$lower, g$complement$debiased$upper))
cat(sprintf("BOOT nb_boots %d; H2 %.3f\n", 20L, fs_boot$H_estimates$H2))
tab <- mr_estimates_table(fs_surv, boot_results = fs_boot, est.scale = "hr")
stopifnot(inherits(tab, "gt_tbl"))
gt::gtsave(tab, "/home/larryleon/Documents/GitHub/forestsearch/dev/reports/articles_finish_partB_2026-10-30/example1_mr_table.html")
cat("EXAMPLE1_OK\n")
