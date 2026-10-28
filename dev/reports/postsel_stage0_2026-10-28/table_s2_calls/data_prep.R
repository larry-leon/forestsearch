# ---------------------------------------------------------------------------
# Table S2 of León & Anderson (2026a), data preparation and configuration -- verbatim transplant.
# Source: forestsearch quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd
# Commit: cf4d64328b6288ce2ec009f2c2fb15689de6a4cb (2026-08-15 21:12 -0700), the last
#   commit touching the .qmd before the Table S2 payload (built 2026-08-15 22:50:30 PDT,
#   Mac-Studio-3.local, forestsearch 0.2.0); render html committed in 43b051b6.
#   At HEAD the .qmd differs from cf4d6432 only at lines 144-145 (run_cv / run_loo
#   FALSE -> TRUE, 6ad5a724); no transplanted chunk below is affected.
# Chunks (name:fence lines, inclusive): setup:92-99 parallel-setup:106-270 load-packages:314-331 data-setup:498-522 dina-config:1233-1243
# Chunk bodies are copied verbatim; only these header comments and the
# '## ---- chunk: <name> ----' separators were added. Transplanted 2026-10-28.
# ---------------------------------------------------------------------------

## ---- chunk: setup (lines 93-98) ----
#| include: false
knitr::opts_chunk$set(comment = "#>", fig.align = "center")

# Initialize timing
timings <- list()
t_vignette_start <- proc.time()

## ---- chunk: parallel-setup (lines 107-269) ----
# PARALLELISM NOTES
#   * FS search -- forestsearch() parallelises the consistency search via its
#     own parallel_args; the plan is set per call and restored on exit.
#   * FB  -- forestsearch_bootstrap_dofuture() uses %dofuture% with
#     .options.future = list(seed = TRUE), so each replicate holds its own
#     L'Ecuyer-CMRG substream.  It inherits the plan set here.
#   * MR  -- vectorised matrix algebra (a single crossprod over all draws), NOT
#     a parallel loop.  It needs no plan, only a threaded BLAS to benefit.  Do
#     not wrap it in a future.
#   * CV / LOO -- forestsearch_Kfold() and forestsearch_tenfold() pin each
#     fold's inner refit to plan = "sequential", workers = 1, so the fold loop
#     workers are the only ones in play.  No nested parallelism to guard.
# Every parallel block below is followed by plan("sequential").

# Worker count for every parallel block below: 80% of PHYSICAL cores (at least
# one).  This is a render-time setting, not the package default -- forestsearch's
# own .default_parallel_workers() is what applies the CRAN 2-core cap under
# R CMD check, and it is not in play here.
n_cores_total <- parallel::detectCores()
n_cores <- max(1L, floor(0.80 * parallel::detectCores(logical = FALSE)))
cat("Using", n_cores, "of", n_cores_total, "total cores for parallel processing")

# Number of bootstrap iterations
# Use 500-2000 for production; reduced here for vignette
NB <- 1000

# ---------------------------------------------------------------------------
# Validation toggles (the slowest parts of the document).  Independent:
# run_cv  controls K-Fold cross-validation; run_loo controls leave-one-out
# (N-fold / OOB).  Either can be turned off without affecting the other.
# ---------------------------------------------------------------------------
#   run_cv  : TRUE runs repeated K-fold CV (kfold-cv plus the diagnostic
#             tabulations and the identification plot).  FALSE skips all of
#             them; the cv-skipped-note chunk states so in the rendered output
#             and timings$kfold stays NA.
#   run_loo : TRUE runs leave-one-out / N-fold OOB.  FALSE skips it, the
#             loo-skipped-note chunk says so, and timings$oob stays NA.
run_cv  <- FALSE    # K-Fold cross-validation (kfold-cv + diagnostics)
run_loo <- FALSE    # Leave-one-out / N-fold OOB (oob-cv)
# K-fold CV with multiple iterations
# Use Ksims >= 50 for production
Ksims    <- 50
Kfolds_cv <- 10     # folds per iteration; feeds BOTH the call and the narrative,
                    # so the "Ksims x Kfolds" count in the prose cannot drift

# Pre-seed CV/LOO timings so the Computational Timing table is well-formed
# whether or not these run (skipped analyses report NA).
timings$kfold <- NA_real_
timings$oob   <- NA_real_

# ---------------------------------------------------------------------------
# Per-analysis identification parameters (tune each analysis independently).
#   sg_focus      : which candidate subgroup the selection step chooses
#   hr_threshold  : harm/benefit floor on the HR scale
#   hr_consistency: minimum consistency rate (consistency search only)
# ---------------------------------------------------------------------------
# Main ForestSearch (consistency search)
fs_sg_focus              <- "effMaxSG"
fs_hr_threshold          <- 1.0
fs_hr_consistency        <- 1.0

fs_conf_force <- c("er <= 0", "pgr <= 0")

fs_consistency_method    <- "resample"
# MR: a fast multiplier-bootstrap approximation of the
# bootstrap bias-correction, computed inside the fit and compared to the full
# bootstrap below.  It runs for whichever selection engine fs_subgroup_method
# names, de-biasing against that engine's own candidate family (it re-derives
# the per-candidate Cox dfbeta itself).
# run_mr : TRUE passes mr_inference = TRUE to every forestsearch() call, so the
#          MR de-biased estimate is computed inside each fit.  FALSE passes
#          mr_inference = FALSE and no MR compute runs at all; the MR rows and
#          tables then report that it was skipped rather than showing blanks.
run_mr <- TRUE     # add mr_inference = TRUE to the forestsearch call
mr_draws      <- 5000L    # multiplier-bootstrap draws for MR

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
# Selection engine for the MAIN analysis.  The headline fit, its bootstrap,
# the MR, CV/LOO, and the KM plot all flow from this single fit,
# so flipping subgroup_method re-runs the entire main analysis under that
# engine (each downstream step inherits the method from the fitted object):
#   "consistency" - the splitting-consistency search (default);
#   "dina"        - DINA performs selection (consistency search bypassed);
#   "grf"         - GRF performs selection (policy tree or DR-frontier).
# The MR de-biases against whichever engine's OWN candidate family,
# with re-selection matching that engine's selection rule.
# ---------------------------------------------------------------------------
fs_subgroup_method <- "consistency"  # "consistency" | "dina" | "grf"
fs_grf_selection   <- "frontier"      # only if fs_subgroup_method == "grf": "tree" | "frontier"
fs_grf_select_statistic <- "effect"       # only if subgroup_method=="grf" & grf_selection=="frontier": "dr" | "effect".
                                      #   "dr"     = native DR-score frontier winner;
                                      #   "effect" = re-rank the DR frontier on the inferential Cox HR (the
                                      #              statistic the MR de-biases).  No-op for the policy tree.
fs_grf_depth       <- 2L              # GRF depth (screening + GRF selection)
fs_dmin_grf        <- 4               # GRF DR-score screening floor
fs_dina_args       <- list()          # extra args forwarded to DINA selection
fs_dina_select_statistic <- "effect"    # only if subgroup_method=="dina": "dina" | "effect".
                                      #   "dina"   = native DINA tau-hat winner;
                                      #   "effect" = re-rank the DINA family on the inferential Cox HR (the
                                      #              statistic the MR de-biases).
# MR (main analysis): MR's de-biasing family for "dina"/"grf" is
# restricted to candidates competitive on the method's NATIVE statistic (DINA
# tau-hat / GRF DR-score).  NULL = use the recommended default, which matches the
# band to effect_neighborhood automatically (the band the engine selects within);
# give a number in (0, 1) to set the band explicitly, or >= 1 to de-bias against
# the full family.  Ignored for "consistency" (its maxcons rule already
# matches the search).
fs_family_native_nbhd <- NULL

# (A) ForestSearch + DINA screening (consistency search selects; DINA only
#     contributes candidate cuts)
dina_screen_sg_focus             <- "effMaxSG"
dina_screen_hr_threshold         <- 0.9
dina_screen_hr_consistency       <- 0.8
# Screening cut set: FALSE => full DINA frontier; TRUE => only the single
# dina_subgroup()-selected cut (uses dina_screen_sg_focus / selection_rule).
dina_screen_selected_only        <- TRUE

# (C) ForestSearch subgroup_method = "dina" (DINA performs selection; the
#     consistency search is bypassed, so hr_consistency does not apply here).
#     The standalone DINA analysis (B) reuses these so (B) and (C) match.
fs_dina_sg_focus     <- "effMaxSG"
fs_dina_hr_threshold <- 0.90
# MR: the DINA de-biasing family is restricted to candidates competitive
# on DINA's native tau-hat.  NULL = recommended default (matches effect_neighborhood
# automatically); set a number in (0, 1) for an explicit band, or >= 1 for the full
# family.  Note: this only tunes candidate competition -- MR holds the forest
# fixed, so the full bootstrap stays the reference for DINA regardless.
fs_dina_family_native_nbhd <- NULL

# (D) ForestSearch subgroup_method = "grf" (GRF performs selection via a causal
#     forest policy tree or DR-score frontier; the consistency search, LASSO and
#     DINA are bypassed).  Dedicated demo: always runs GRF regardless of
#     fs_subgroup_method, mirroring the standalone DINA selection (C).
fs_grf_demo_selection     <- "frontier"  # "frontier" | "tree"
fs_grf_demo_depth         <- 2L          # GRF tree / DR-pair depth
fs_grf_demo_dmin          <- 0.0         # DR-score harm floor (0 = permissive default)
fs_grf_sg_focus           <- "effMaxSG"  # -> GRF frontier/tree selection rule
# MR: the GRF de-biasing family is restricted to candidates competitive
# on GRF's native (RMST-based) DR-score.  NULL = recommended default (matches
# effect_neighborhood automatically); set a number in (0, 1) for an explicit band,
# or >= 1 for the full family.  Note: this only tunes candidate competition -- the
# MR holds the forest fixed, so the full bootstrap stays the reference for GRF
# regardless (a dominant DR-score winner can collapse the family to 1).
fs_grf_family_native_nbhd <- NULL

## ---- chunk: load-packages (lines 315-330) ----
#| code-fold: false
library(forestsearch)
library(survival)
library(data.table)
library(ggplot2)
library(gt)
library(grf)
library(policytree)
library(doFuture)

# Optional packages for enhanced output
library(patchwork)
library(weightedsurv)

# Set ggplot theme
theme_set(theme_minimal(base_size = 12))

## ---- chunk: data-setup (lines 499-521) ----
# Load GBSG data (included in forestsearch package)
df.analysis <- gbsg

# Prepare analysis variables
df.analysis <- within(df.analysis, {
  id <- seq_len(nrow(df.analysis))
  time_months <- rfstime / 30.4375
  grade3 <- ifelse(grade == "3", 1, 0)
  treat <- hormon
})

# Define variable roles
confounders.name <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")
outcome.name <- "time_months"
event.name <- "status"
id.name <- "id"
treat.name <- "hormon"

# Display data structure
cat("Sample size:", nrow(df.analysis), "\n")
cat("Events:", sum(df.analysis[[event.name]]),
    sprintf("(%.1f%%)\n", 100 * mean(df.analysis[[event.name]])))
cat("Baseline factors:", paste(confounders.name, collapse = ", "), "\n")

## ---- chunk: dina-config (lines 1234-1242) ----
#| code-fold: false

# Shared DINA *fit* settings (Cox / survival outcome => family mapped to "cox").
# Selection parameters (sg_focus, hr_threshold) live in the setup chunk under
# the fs_dina_* knobs so they can be tuned independently; (B) and (C) read them.
dina_seed   <- 8316951
dina_nfolds <- 2
# Harm/benefit floor on the natural-parameter (log-HR) scale for ratio measures
m_diff_dina <- log(fs_dina_hr_threshold)

