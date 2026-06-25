# =============================================================================
# build_actg175_glm_dgm.R
# -----------------------------------------------------------------------------
# Builds the ACTG175-anchored binary / odds-ratio DGM for the GLM coverage
# study (post-selection manuscript, review item A7).  This is the OR-scale
# analogue of the survival GBSG super-population: a known harm subgroup with a
# known effect, on the same trial, outcome coding, and covariate set as the
# manuscript's Section 6.2 ACTG175 application.
#
# Anchor (mirrors analysis_actg175_binary_multimethod_psi_v3a.qmd exactly):
#   - arms 1 (ZDV+ddI, experimental) vs 3 (ddI, control); treat = (arms == 1)
#   - outcome y_neg = 1 - (cd420 > cd40): the week-20 "no CD4 improvement"
#     (adverse) indicator; OR > 1 = harm
#   - 6 continuous + 6 binary baseline covariates
#
# Planted harm subgroup (choice (i), mirroring the FS/GRF discovery):
#   H = { wtkg > q70 AND cd40 > q70 }  -- the upper-tail (high weight, high CD4)
#   conjunction, in-family for the consistency / GRF search.
#
# Verified truth targets at the defaults below (n_super = 25000, seed 8316951):
#   prevalence(H) = 9.88%
#   marginal  OR(H)  = 1.505   OR(Hc)  = 0.656   OR(ITT)  = 0.717   (theta-dagger)
#   CDE       CDE(H) = 1.553   CDE(Hc) = 0.631   CDE(ITT) = 0.766   (theta-double-dagger)
#
# NOTE: forestsearch must be INSTALLED (devtools::install()), not load_all(),
# if this DGM is later consumed by doFuture parallel workers.
# =============================================================================

library(forestsearch)
library(speff2trial)

# Null-coalescing helper (local; independent of forestsearch's internal one).
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a


# ── Knobs ────────────────────────────────────────────────────────────────────
analysis_seed <- 8316951L
n_super       <- 25000L     # super-population size (n-invariant truth)
target_or_h   <- 1.5        # planted marginal OR in the harm subgroup H
sg_quantile   <- 0.70       # q-cut on each subgroup variable:
                            #   0.70 -> ~10% prevalence (matches discovered 7-9%)
                            #   0.65 -> ~13% prevalence (matches GBSG H's 12.5%)

out_path      <- "dgm_actg175_or15.rds"   # where to save the built DGM (NULL = don't save)

# ── Variable roles (identical to the Section 6.2 application) ────────────────
cont_vars   <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
bin_vars    <- c("hemo", "homo", "drugs", "race", "gender", "symptom")

# ── Planted harm subgroup H (choice (i): high weight AND high CD4) ────────────
# Upper-tail cuts are encoded with type = "greater" (a bare numeric or
# type = "quantile" would encode the LOWER tail, '<=', as in the GBSG H).
subgroup_vars <- c("wtkg", "cd40")
subgroup_cuts <- list(
  wtkg = list(type = "greater", quantile = sg_quantile),
  cd40 = list(type = "greater", quantile = sg_quantile)
)

# -- Alternative definition (ii), noted for reference (a simpler single cut):
#      subgroup_vars <- "cd40"
#      subgroup_cuts <- list(cd40 = list(type = "greater", quantile = 0.75))
#    -> H = { cd40 > q75 }, ~25% prevalence on the ACTG175 super-population.

# ── Data: ACTG175, prepped exactly as the application driver ──────────────────
set.seed(analysis_seed)
actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
actg_df$id    <- seq_len(nrow(actg_df))
actg_df$treat <- ifelse(actg_df$arms == 1L, 1L, 0L)   # ZDV+ddI = 1, ddI = 0
actg_df <- actg_df[!is.na(actg_df$cd420), ]           # drop missing week-20 CD4
actg_df$y_pos <- as.integer(actg_df$cd420 > actg_df$cd40)  # 1 = improved (better)
actg_df$y_neg <- 1L - actg_df$y_pos                        # 1 = no improvement (worse)

cat(sprintf("ACTG175: n = %d, adverse rate (y_neg) = %.1f%%\n",
            nrow(actg_df), 100 * mean(actg_df$y_neg)))

# ── Build + calibrate the DGM to the target harm OR ──────────────────────────
# adverse_outcome = FALSE here is INTENTIONAL: with the outcome already coded on
# the adverse y_neg scale, a positive interaction shift raises P(no improvement)
# for H under treatment, planting OR(H) > 1 (harm) directly on the analysis
# scale.  The COMPLEMENT inherits the fitted ACTG175 treatment effect, ~0.66
# (ZDV+ddI is protective outside H) -- the theta-dagger(Hc) ~ 0.65 parallel.
#
# The downstream forestsearch() / oracle calls in the sweep run on this same
# y_neg scale with adverse_outcome = TRUE (no Y-flip), matching the application.
dgm <- calibrate_glm_interaction(
  data            = actg_df,
  factor_vars     = bin_vars,
  continuous_vars = cont_vars,
  outcome_var     = "y_neg",
  treatment_var   = "treat",
  target_effect   = target_or_h,
  outcome_type    = "binary",
  effect_measure  = "OR",
  subgroup_vars   = subgroup_vars,
  subgroup_cuts   = subgroup_cuts,
  k_treat         = 1,
  adverse_outcome = FALSE,
  k_inter_range   = c(0.3, 1.5),
  grid_step       = 0.025,
  n_super         = n_super,
  seed            = analysis_seed,
  verbose         = FALSE
)

# ── Report the verified truth targets ────────────────────────────────────────
hr <- dgm$hazard_ratios
cat("\n===== Calibrated ACTG175 GLM DGM (OR scale, y_neg; >1 = harm) =====\n")
cat(sprintf("  calibrated k_inter : %.4f\n", dgm$model_params$k_inter))
cat(sprintf("  prevalence(H)      : %.2f%%  (%d / %d)\n",
            100 * dgm$subgroup_info$proportion, dgm$subgroup_info$size, dgm$n_super))
cat("  -- marginal OR (theta-dagger; AHR_*) --\n")
cat(sprintf("     OR(H)  = %.4f   OR(Hc)  = %.4f   OR(ITT)  = %.4f\n",
            hr$AHR_harm, hr$AHR_no_harm, hr$AHR))
cat("  -- CDE (theta-double-dagger; Jensen-corrected mean-odds ratio) --\n")
cat(sprintf("     CDE(H) = %.4f   CDE(Hc) = %.4f   CDE(ITT) = %.4f\n",
            hr$CDE_harm, hr$CDE_no_harm, hr$CDE))

# Truth-target list in the shape the coverage sweep consumes
# (marg_* = theta-dagger, cde_* = theta-double-dagger).
truth <- list(
  or_causal = hr$AHR          %||% NA_real_,   # overall marginal OR
  marg_H    = hr$AHR_harm     %||% NA_real_,   # theta-dagger (H)
  marg_Hc   = hr$AHR_no_harm  %||% NA_real_,   # theta-dagger (Hc)
  cde_H     = hr$CDE_harm     %||% NA_real_,   # theta-double-dagger (H)
  cde_Hc    = hr$CDE_no_harm  %||% NA_real_    # theta-double-dagger (Hc)
)

if (!is.null(out_path)) {
  saveRDS(list(dgm = dgm, truth = truth, actg_df = actg_df,
               cont_vars = cont_vars, bin_vars = bin_vars,
               subgroup_vars = subgroup_vars, subgroup_cuts = subgroup_cuts,
               sg_quantile = sg_quantile, target_or_h = target_or_h,
               seed = analysis_seed, n_super = n_super),
          out_path)
  cat(sprintf("\nDGM bundle saved to %s\n", out_path))
}
