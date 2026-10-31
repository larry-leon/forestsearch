# ---- transplant (verbatim chunk bodies; fences removed) ----
# source: forestsearch/quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd
# forestsearch commit 43b051b6 (2026-08-16 12:07; the render that wrote the payload at 2026-08-16 16:50:40 UTC ran this text — its HTML is committed in the same commit; identical to HEAD; differs from fd1ce74f, the last commit before the build, in fs_selection_rule both->neighborhood)
#   chunk 'setup': fence lines 39-46
#   chunk 'parallel-setup': fence lines 53-243
#   chunk 'load-packages': fence lines 317-336
#   chunk 'analysis-roles': fence lines 340-360
#   chunk 'data-treat': fence lines 367-374
#   chunk 'data-outcomes': fence lines 378-393
#   chunk 'data-covariates': fence lines 397-411
# ----------------------------------------------------------

## ==== chunk 'setup' ====
#| include: false
knitr::opts_chunk$set(comment = "#>", fig.align = "center")

# Initialize timing
timings <- list()
t_vignette_start <- proc.time()

## ==== chunk 'parallel-setup' ====
# Detect available cores
n_cores_total <- parallel::detectCores()
n_cores <- max(1L, floor(0.80 * parallel::detectCores(logical = FALSE)))
cat("Using", n_cores, "of", n_cores_total, "total cores for parallel processing")

# Number of bootstrap iterations
# Use 500-2000 for production; reduced here for demo / pressure test
NB <- 1000L

# ---------------------------------------------------------------------------
# Validation toggles (the slowest parts of the document).  Independent:
# run_cv  controls K-Fold cross-validation; run_loo controls leave-one-out
# (N-fold / OOB).  Either can be turned off without affecting the other.
# ---------------------------------------------------------------------------
run_cv  <- FALSE    # K-Fold cross-validation (kfold-cv + diagnostics)
run_loo <- FALSE    # Leave-one-out / N-fold OOB (oob-cv)
# 10-fold CV with multiple iterations.  Use Ksims >= 50 for production.
Ksims <- 50L

# Pre-seed CV/LOO timings so the Computational Timing table is well-formed
# whether or not these run (skipped analyses report NA).
timings$kfold <- NA_real_
timings$oob   <- NA_real_

# ---------------------------------------------------------------------------
# Per-analysis identification parameters (tune each analysis independently).
#   sg_focus       : which candidate subgroup the selection step chooses
#   or_threshold   : harm floor on the OR scale (adverse scale: OR > c1 = harm)
#   or_consistency : minimum per-split consistency OR (consistency search only)
# ---------------------------------------------------------------------------
# Main ForestSearch (consistency search)
fs_sg_focus              <- "effMaxSG"   # GLM-natural alias of "hrMaxSG"
fs_selection_rule        <- "neighborhood"       # neighborhood ∩ Pareto (effMaxSG/effMinSG)
fs_or_threshold          <- 1.0
fs_or_consistency        <- 1.0
fs_pconsistency          <- 0.90
# IMPORTANT: the Tier-2 de-biased gate on the GLM path runs ONLY under the
# resample-consistency engine (consistency_method = "resample").  Under the
# default "split" engine the gate is silently skipped for GLM outcomes.
fs_consistency_method    <- "resample"
# Tier-2 de-biased gate: a fast multiplier-bootstrap approximation of the
# bootstrap bias-correction, computed inside the fit and compared to the full
# bootstrap below.  It runs for whichever selection engine fs_subgroup_method
# names, de-biasing against that engine's own candidate family (it re-derives
# the per-candidate GLM dfbeta itself).
run_mr <- TRUE     # add mr_inference = TRUE to the forestsearch call
mr_draws      <- 5000L    # multiplier-bootstrap draws for the gate

# ---------------------------------------------------------------------------
# Output location for this document's results:  <base>/<dirout>/
#   results_dir = NULL -> base is <qmd dir>/_payloads
#   dirout      = NULL -> the .qmd stem
# Setting results_dir redirects the base only; each run still gets its own
# subdirectory named by dirout.
#
# Examples:
#   results_dir <- NULL;  dirout <- NULL
#     -> <qmd dir>/_payloads/<stem>/<stem>_payload.rds
#   results_dir <- NULL;  dirout <- "gbsg_survival_multimethod"
#     -> <qmd dir>/_payloads/gbsg_survival_multimethod/gbsg_survival_multimethod_payload.rds
#   results_dir <- "~/GitHub/fs-post-selection-manuscript_jasa/_payloads"; dirout <- NULL
#     -> ~/GitHub/fs-post-selection-manuscript_jasa/_payloads/<stem>/<stem>_payload.rds
# ---------------------------------------------------------------------------
results_dir <- NULL
dirout      <- NULL

# ---------------------------------------------------------------------------
# FIXED-FAMILY CONFIGURATION (use_lasso = FALSE, use_grf = FALSE)
# ---------------------------------------------------------------------------
# Fully enumerated candidate family: no LASSO or GRF screening, and as of
# 2026-08 no max_subgroups_search truncation.
#
# This is the fixed-family configuration the method's Theorem 1 certifies:
# the family is enumerated (default cuts, with J = 10 cuts for wtkg and
# cd40) rather than model-screened, so the bias correction captures all
# selection optimism.
#
# The published supplement Table S11 ({wtkg>86} & {cd40>380}, n=72,
# naive OR 3.58, FB 2.43, MR 1.80) came from an earlier run of this document
# that capped candidate evaluation at 50 (max_subgroups_search = 50L). That
# cap has been removed so the analysis follows the recommended specification
# rather than the earlier run. Results may therefore differ from the
# published table; the difference is intentional.
#
# Prior run (superseded): payload built 2026-07-08 on the 128-core machine,
# 5000 MR draws, max_subgroups_search = 50L.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Selection engine for the MAIN analysis.  The headline fit, its bootstrap,
# the de-biased gate, CV/LOO, and the subgroup outcome plot all flow from this
# single fit, so flipping subgroup_method re-runs the entire main analysis
# under that engine (each downstream step inherits the method from the fitted
# object):
#   "consistency" - the splitting-consistency search (default);
#   "dina"        - DINA performs selection (consistency search bypassed);
#   "grf"         - GRF performs selection (policy tree or DR-frontier).
# The Tier-2 gate de-biases against whichever engine's OWN candidate family,
# with re-selection matching that engine's selection rule.
# ---------------------------------------------------------------------------
fs_subgroup_method <- "consistency"  # "consistency" | "dina" | "grf"
fs_grf_selection   <- "frontier"     # only if subgroup_method == "grf": "tree" | "frontier"
fs_grf_select_statistic <- "effect"  # only if subgroup_method=="grf" & grf_selection=="frontier": "dr" | "effect".
                                     #   "dr"     = native DR-score frontier winner;
                                     #   "effect" = re-rank the DR frontier on the inferential GLM OR (the
                                     #              statistic the Tier-2 gate de-biases).  No-op for the policy tree.
fs_grf_depth       <- 2L             # GRF depth (screening + GRF selection)
fs_dmin_grf        <- 0              # GRF DR-score screening floor (RD scale for binary; 0 = permissive)
fs_dina_args       <- list()         # extra args forwarded to DINA selection
fs_dina_select_statistic <- "effect" # only if subgroup_method=="dina": "dina" | "effect".
                                     #   "dina"   = native DINA tau-hat winner;
                                     #   "effect" = re-rank the DINA family on the inferential GLM OR (the
                                     #              statistic the Tier-2 gate de-biases).
# Tier-2 gate (main analysis): the gate's re-selection family for "dina"/"grf" is
# restricted to candidates competitive on the method's NATIVE statistic (DINA
# tau-hat / GRF DR-score).  NULL = use the recommended default, which matches the
# band to effect_neighborhood automatically (the band the engine selects within);
# give a number in (0, 1) to set the band explicitly, or >= 1 to de-bias against
# the full family.  Ignored for "consistency" (its maxcons re-selection already
# matches the search).
fs_family_native_nbhd <- NULL

# (A) ForestSearch + DINA screening (consistency search selects; DINA only
#     contributes candidate cuts)
dina_screen_sg_focus             <- "effMaxSG"
dina_screen_or_threshold         <- 1.0
dina_screen_or_consistency       <- 0.9
# Screening cut set: FALSE => full DINA frontier; TRUE => only the single
# dina_subgroup()-selected cut (uses dina_screen_sg_focus / selection_rule).
dina_screen_selected_only        <- TRUE

# (C) ForestSearch subgroup_method = "dina" (DINA performs selection; the
#     consistency search is bypassed, so or_consistency does not apply here).
#     The standalone DINA analysis (B) reuses these so (B) and (C) match.
fs_dina_sg_focus     <- "effMaxSG"
fs_dina_or_threshold <- 1.0
# Tier-2 gate: the DINA re-selection family is restricted to candidates competitive
# on DINA's native tau-hat.  NULL = recommended default (matches effect_neighborhood
# automatically); set a number in (0, 1) for an explicit band, or >= 1 for the full
# family.  Note: this only tunes candidate competition -- the gate holds the forest
# fixed, so the full bootstrap stays the reference for DINA regardless.
fs_dina_family_native_nbhd <- NULL

# (D) ForestSearch subgroup_method = "grf" (GRF performs selection via a causal
#     forest policy tree or DR-score frontier; the consistency search, LASSO and
#     DINA are bypassed).  Dedicated demo: always runs GRF regardless of
#     fs_subgroup_method, mirroring the standalone DINA selection (C).
fs_grf_demo_selection     <- "frontier"  # "frontier" | "tree"
fs_grf_demo_depth         <- 2L          # GRF tree / DR-pair depth
fs_grf_demo_dmin          <- 0.0         # DR-score harm floor (RD scale; 0 = permissive default)
fs_grf_sg_focus           <- "effMaxSG"  # -> GRF frontier/tree selection rule
# Tier-2 gate: the GRF re-selection family is restricted to candidates competitive
# on GRF's native (RD-based) DR-score.  NULL = recommended default (matches
# effect_neighborhood automatically); set a number in (0, 1) for an explicit band,
# or >= 1 for the full family.  Note: this only tunes candidate competition -- the
# gate holds the forest fixed, so the full bootstrap stays the reference for GRF
# regardless (a dominant DR-score winner can collapse the family to 1).
fs_grf_family_native_nbhd <- NULL

# ---------------------------------------------------------------------------
# Numeric precision (two-stage control).
#   construction_digits : upstream sprintf precision in the bootstrap / CV
#     table builders.  Acts as a CEILING for downstream display.
#   display_digits      : precision rendered in user-facing tables.
# ---------------------------------------------------------------------------
construction_digits <- 2L
display_digits      <- 2L

# Display-only tidy of a subgroup label for printing: rounds the post-operator
# thresholds to `d` decimals via the package's tidy_cut_display().  The
# membership-exact sg.harm (re-parsed by get_dfpred and carried through the
# bootstrap) is never altered -- this only affects what is printed.  Used at
# every site in this document that renders sg.harm.
fmt_sg <- function(x, d = display_digits) {
  if (is.null(x) || length(x) == 0L) return(x)
  paste(forestsearch:::tidy_cut_display(x, digits = as.integer(d)), collapse = " & ")
}

# ---------------------------------------------------------------------------
# Pareto-frontier diagnostic parameters (drives compute_frontier_cis(),
# plot_pareto_frontier(), pareto_frontier_table(), explain_pareto_selection()).
# For OR (ratio-scale), the effect axis is auto-rendered on log scale.
# ---------------------------------------------------------------------------
frontier_params <- list(
  n_splits      = 1000L,
  seed          = 1L,
  xlim          = c(0.35, 10),   # OR auto-detected as ratio scale
  digits        = 2L
)

## ==== chunk 'load-packages' ====
#| code-fold: false
library(forestsearch)
library(speff2trial)
library(survival)
library(data.table)
library(ggplot2)
library(gt)
library(grf)
library(policytree)
library(doFuture)
library(future)

# Optional packages for enhanced output
library(patchwork)
library(forestploter)

# Set ggplot theme
theme_set(theme_minimal(base_size = 12))

## ==== chunk 'analysis-roles' ====
#| code-fold: false

# ─── Reproducibility ────────────────────────────────────────────────────
analysis_seed <- 8316951L
set.seed(analysis_seed)

# ─── Confounders ──────────────────────────────────────────────────────
cont_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
bin_vars  <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
confounders.name <- c(cont_vars, bin_vars)
id.name   <- "id"

# ─── Variables (single adverse scale) ──────────────────────────────────
treat.name   <- "treat"
outcome.name <- "y_binary_neg"   # FS adverse-scale outcome (1 = no improvement)

# ─── ForestSearch GLM parameters (binary, OR scale) ────────────────────
fs_outcome_type   <- "binary"
fs_effect_measure <- "OR"

## ==== chunk 'data-treat' ====
# ACTG175: arm 1 = ZDV+ddI (experimental), arm 3 = ddI (control)
actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
actg_df$id <- seq_len(nrow(actg_df))

# Treatment: ZDV+ddI = 1, ddI = 0 (no switching)
actg_df$treat <- ifelse(actg_df$arms == 1L, 1L, 0L)

## ==== chunk 'data-outcomes' ====
# Drop subjects missing the week-20 CD4 measurement
actg_df <- actg_df[!is.na(actg_df$cd420), ]

# Beneficial scale: CD4 improvement indicator (1 = improved, higher = better)
actg_df$y_binary_pos <- as.integer(actg_df$cd420 > actg_df$cd40)

# Adverse scale: NO improvement indicator (1 = did not improve, higher = worse).
# This is the single analysis scale for the entire multimethod workflow.
actg_df$y_binary_neg <- 1L - actg_df$y_binary_pos

cat("Sample size:", nrow(actg_df), "\n")
cat("Adverse events (no improvement):", sum(actg_df$y_binary_neg),
    sprintf("(%.1f%%)\n", 100 * mean(actg_df$y_binary_neg)))
cat("Baseline factors:", paste(confounders.name, collapse = ", "), "\n")

## ==== chunk 'data-covariates' ====
# Keep ALL candidate covariates numeric.  forestsearch() builds and evaluates
# "var <= value" candidate cuts, and its DINA paths require numeric columns, so
# binary indicators left as factors trigger hundreds of "'<=' not meaningful for
# factors" warnings during the consistency search and "covariate must be numeric"
# failures in DINA (both subgroup_method = "dina" selection and use_dina
# screening).  With 0/1 numeric coding, forestsearch's cont.cutoff (= 4 unique
# values) treats the binary candidates as categorical automatically -- no factor
# conversion is needed, and every engine (consistency, DINA, GRF, standalone
# DINA) receives numeric covariates.
for (v in bin_vars) actg_df[[v]] <- as.numeric(actg_df[[v]])

cat("Continuous candidates:", paste(cont_vars, collapse = ", "), "\n")
cat("Binary candidates:    ", paste(bin_vars, collapse = ", "), "\n")

