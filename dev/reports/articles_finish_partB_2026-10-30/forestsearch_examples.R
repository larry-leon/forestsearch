### Name: forestsearch
### Title: ForestSearch: Exploratory Subgroup Identification
### Aliases: forestsearch

### ** Examples


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

## Example 2 -- Binary outcome, odds ratio (ACTG175, ZDV+ddI vs ddI):
## the forest-search analysis of León and Anderson (2026a, Table S3); the
## binary analysis of León and Anderson (2026b) selects the same region
## with effect_neighborhood = 0.20. Adverse-scale outcome (1 = no CD4
## improvement at week 20), OR thresholds at the null, same selection rule
## and multiplier-resampling gate as Example 1.
if (requireNamespace("speff2trial", quietly = TRUE)) {
  actg <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
  actg <- actg[!is.na(actg$cd420), ]
  actg$id           <- seq_len(nrow(actg))
  actg$treat        <- as.integer(actg$arms == 1L)
  actg$y_binary_neg <- as.integer(!(actg$cd420 > actg$cd40))
  bin_vars <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
  for (v in bin_vars) actg[[v]] <- as.numeric(actg[[v]])
  fs_bin <- forestsearch(
    df.analysis            = actg,
    outcome.name           = "y_binary_neg",
    event.name             = "y_binary_neg",
    treat.name             = "treat",
    id.name                = "id",
    confounders.name       = c("age", "preanti", "wtkg", "karnof", "cd40",
                               "cd80", bin_vars),
    outcome_type           = "binary",
    effect_measure         = "OR",
    adverse_outcome        = TRUE,
    use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
    conf.cont_jcuts        = list(cd40 = 10, wtkg = 10),
    subgroup_method        = "consistency",
    sg_focus               = "effMaxSG",
    selection_rule         = "neighborhood",
    effect_neighborhood    = 0.10,
    effect.threshold       = 1.0,
    consistency.threshold  = 1.0,
    pconsistency.threshold = 0.90,
    consistency_method     = "resample",
    fs.splits              = 500,
    n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
    mr_inference           = TRUE,
    mr_inference_args      = list(draws = 5000L),
    parallel_args          = list(plan = "multisession", workers = 4),
    seedit                 = 8316951
  )
  print(fs_bin)
}

## Example 3 -- Continuous outcome, mean difference (ACTG175): the
## continuous application of León and Anderson (2026b), with the field
## (calibrated one-sided) bounds for the selected subgroup and its
## complement. Outcome is CD4 decline from baseline to week 20 (adverse
## scale); thresholds are on the mean-difference scale.
if (requireNamespace("speff2trial", quietly = TRUE)) {
  actg <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
  actg <- actg[!is.na(actg$cd420), ]
  actg$id        <- seq_len(nrow(actg))
  actg$treat     <- as.integer(actg$arms == 1L)
  actg$y_decline <- actg$cd40 - actg$cd420
  bin_vars <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
  for (v in bin_vars) actg[[v]] <- as.factor(actg[[v]])
  fs_cont <- forestsearch(
    df.analysis            = actg,
    outcome.name           = "y_decline",
    treat.name             = "treat",
    id.name                = "id",
    confounders.name       = c("age", "preanti", "wtkg", "karnof", "cd40",
                               "cd80", bin_vars),
    outcome_type           = "continuous",
    effect_measure         = "MD",
    adverse_outcome        = TRUE,
    use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
    conf.cont_jcuts        = list(age = 10, preanti = 10, wtkg = 10,
                                  karnof = 10, cd40 = 10, cd80 = 10),
    sg_focus               = "effMaxSG",
    selection_rule         = "neighborhood",
    effect_neighborhood    = 0.20,
    effect.threshold       = 10,
    consistency.threshold  = 5,
    pconsistency.threshold = 0.90,
    consistency_method     = "resample",
    fs.splits              = 500,
    n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
    mr_inference           = TRUE,
    mr_inference_args      = list(ci_method = "field", draws = 5000L,
                                  include_complement = TRUE,
                                  field_complement = TRUE,
                                  field_scale_complement = "selected",
                                  return_reselection = TRUE),
    parallel_args          = list(plan = "sequential"),
    seedit                 = 8316951
  )
  print(fs_cont)
}





