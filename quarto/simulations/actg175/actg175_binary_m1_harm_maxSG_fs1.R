# ============================================================================
# actg175_binary_m1_harm_maxSG_fs1.R
#
# Background simulation runner for configuration:
#   study     : ACTG175
#   outcome   : binary
#   DGM       : m1   — calibrated GLM interaction (target OR(Q) = 2.0)
#   objective : harm
#   sg_focus  : maxSG
#   FS bundle : fs1
#
# Performs:
#   1. Data preparation (ACTG175 arms 1 & 3, binary adverse outcome)
#   2. DGM calibration (H1: target OR(Q) = 2.0; H0: null)
#   3. Parallel simulation runs via run_sim_actg175_binary_m1_harm_maxSG_fs1()
#   4. saveRDS() of the resulting bundle
#
# ----------------------------------------------------------------------------
# Usage
# ----------------------------------------------------------------------------
# From the forestsearch package root (working directory):
#
#   # Interactive / RStudio Background Jobs:
#   source("quarto/simulations/actg175/actg175_binary_m1_harm_maxSG_fs1.R")
#
#   # Terminal:
#   Rscript quarto/simulations/actg175/actg175_binary_m1_harm_maxSG_fs1.R \
#     > quarto/simulations/actg175/logs/run_$(date +%Y%m%d_%H%M%S).log 2>&1
#
# Output: quarto/simulations/actg175/actg175_binary_m1_harm_maxSG_fs1.rds
#
# ----------------------------------------------------------------------------
# Companion files (same stem):
#   actg175_binary_m1_harm_maxSG_fs1.qmd      — pkgdown article (vignette)
#   sim_actg175_binary_m1_harm_maxSG_fs1.R    — simulation core (sourced below)
# ============================================================================

# ── 1. Libraries ────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(forestsearch)
  library(speff2trial)
  library(data.table)
  library(doFuture)
  library(foreach)
  options(warn = -1)
})

# ── 2. Configuration ────────────────────────────────────────────────────────
CONFIG_ID  <- "actg175_binary_m1_harm_maxSG_fs1"
CONFIG_DIR <- file.path("quarto", "simulations", "actg175")
HELPER     <- file.path(CONFIG_DIR, paste0("sim_", CONFIG_ID, ".R"))
OUT_PATH   <- file.path(CONFIG_DIR, paste0(CONFIG_ID, ".rds"))

# Fail fast if the working directory is not the package root
if (!file.exists(HELPER)) {
  stop(sprintf(
    "Helper not found at '%s'.\n  Expected working directory: forestsearch package root.\n  Current wd: %s",
    HELPER, getwd()
  ), call. = FALSE)
}
source(HELPER)

# Ensure output directory exists
dir.create(CONFIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ── 3. Run-scale parameters ─────────────────────────────────────────────────
# These are the production sample sizes for this configuration.
# Edit here to scale up/down for a particular run.
nsims_alt  <- 20L
nsims_null <- 20L
sim_n      <- 1000L

n_workers  <- max(1L, floor(0.95 * parallel::detectCores(logical = FALSE)))

# ── 4. Seeds ────────────────────────────────────────────────────────────────
master_seed <- 8316951L
set.seed(master_seed)

# ── 5. FS parameters (fs1) ──────────────────────────────────────────────────
c1 <- 1.25                       # hr.threshold (OR scale)
c2 <- 1.0                        # hr.consistency
fs_pconsistency  <- 0.80
max_look         <- 10
fs_splits        <- 1000L
fs_n_min         <- 60L
fs_d0_min        <- 10L
fs_d1_min        <- 0L
fs_maxk          <- 2L
fs_dmin_grf      <- 0.0
fs_seedit        <- master_seed
fs_sg_focus      <- "maxSG"
fs_use_lasso     <- FALSE
fs_use_grf       <- TRUE
fs_return_selected_cuts_only <- FALSE
fs_is_rct        <- TRUE

fs_params <- list(
  outcome.name              = "y_sim",
  event.name                = "y_sim",
  treat.name                = "treat_sim",
  id.name                   = "id",
  outcome_type              = "binary",
  effect_measure            = "OR",
  use_lasso                 = fs_use_lasso,
  use_grf                   = fs_use_grf,
  return_selected_cuts_only = fs_return_selected_cuts_only,
  max_subgroups_search      = max_look,
  use_twostage              = TRUE,
  hr.threshold              = c1,
  hr.consistency            = c2,
  pconsistency.threshold    = fs_pconsistency,
  stop_threshold            = 0.90,
  sg_focus                  = fs_sg_focus,
  fs.splits                 = fs_splits,
  n.min                     = fs_n_min,
  d0.min                    = fs_d0_min,
  d1.min                    = fs_d1_min,
  maxk                      = fs_maxk,
  dmin.grf                  = fs_dmin_grf,
  is.RCT                    = fs_is_rct,
  seedit                    = fs_seedit,
  showten_subgroups         = FALSE,
  details                   = FALSE,
  quiet                     = TRUE,
  parallel_args             = list(plan = "sequential", workers = 1,
                                   show_message = FALSE)
)

# ── 6. GRF parameters ───────────────────────────────────────────────────────
grf_params <- list(
  dmin.grf     = 0.10,
  maxdepth     = 2L,
  sg.criterion = "mDiff"
)

# ── 7. Diagnostics & calibration target ─────────────────────────────────────
diag_keep         <- c("fs_full", "candidate_table")
diag_keep_first_n <- 5L
cal_target_or     <- 2.0

# ── 8. Data preparation ─────────────────────────────────────────────────────
# ACTG175: arm 1 = ZDV+ddI (experimental), arm 3 = ddI (control), no switching.
actg_df    <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
actg_df$id <- seq_len(nrow(actg_df))
actg_df$treat <- ifelse(actg_df$arms == 1L, 1L, 0L)

# Binary adverse outcome: failure to improve (cd420 <= cd40)
actg_df$y_binary <- as.integer(actg_df$cd420 <= actg_df$cd40)
actg_df <- actg_df[!is.na(actg_df$cd420), ]

# DGM factor variables
actg_df$z1  <- as.factor(ifelse(actg_df$cd40 > 375, 1L, 0L))
actg_df$z2  <- as.factor(ifelse(actg_df$wtkg > 75,  1L, 0L))
actg_df$z7  <- as.factor(actg_df$hemo)
actg_df$z8  <- as.factor(actg_df$homo)
actg_df$z9  <- as.factor(actg_df$drugs)
actg_df$z10 <- as.factor(actg_df$race)
actg_df$z11 <- as.factor(actg_df$gender)
actg_df$z12 <- as.factor(actg_df$symptom)

actg_df$ar_naive   <- as.factor(ifelse(actg_df$str2 == 0,            1L, 0L))
actg_df$prior_6mo  <- as.factor(ifelse(actg_df$preanti <= 6 * 30.4375, 1L, 0L))

actg_df$flag_harm_ref <- as.integer(actg_df$z1 == 1 & actg_df$z2 == 1)

# Analysis covariates (what the analyst sees)
cont_vars   <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
bin_vars    <- c("hemo", "homo", "drugs", "race", "gender", "symptom",
                 "ar_naive", "prior_6mo")
for (v in bin_vars) actg_df[[v]] <- as.factor(actg_df[[v]])
confounders <- c(cont_vars, bin_vars)

# conf_force depends on data-prep factor variables; inject after they exist
fs_params$conf_force <- c("ar_naive == 1", "prior_6mo == 1")

# DGM variable sets
dgm_factor_vars   <- c("z1", "z2", paste0("z", 7:12))
dgm_cont_vars     <- cont_vars
dgm_subgroup_vars <- c("z1", "z2")
dgm_subgroup_cuts <- list(z1 = 1L, z2 = 1L)

# ── 9. DGM construction ─────────────────────────────────────────────────────
cat(sprintf("[%s] Calibrating H1 DGM (target OR(Q) = %.2f)...\n",
            CONFIG_ID, cal_target_or))
t_cal <- proc.time()[3]
dgm_calibrated <- calibrate_glm_interaction(
  data            = actg_df,
  factor_vars     = dgm_factor_vars,
  continuous_vars = dgm_cont_vars,
  outcome_var     = "y_binary",
  treatment_var   = "treat",
  target_effect   = cal_target_or,
  outcome_type    = "binary",
  effect_measure  = "OR",
  subgroup_vars   = dgm_subgroup_vars,
  subgroup_cuts   = dgm_subgroup_cuts,
  k_inter_range   = c(0, 3),
  grid_step       = 0.05,
  n_super         = 5000L,
  seed            = master_seed,
  verbose         = FALSE
)
cat(sprintf("  OR(Q):   %.3f (target: %.2f)\n",
            dgm_calibrated$hazard_ratios$harm_subgroup, cal_target_or))
cat(sprintf("  OR(Qc):  %.3f\n",
            dgm_calibrated$hazard_ratios$no_harm_subgroup))
cat(sprintf("  OR(ITT): %.3f\n",
            dgm_calibrated$hazard_ratios$overall))
cat(sprintf("  done in %.1f sec\n\n", proc.time()[3] - t_cal))

cat(sprintf("[%s] Building H0 (null) DGM...\n", CONFIG_ID))
dgm_null <- generate_glm_dgm(
  data            = actg_df,
  factor_vars     = dgm_factor_vars,
  continuous_vars = dgm_cont_vars,
  outcome_var     = "y_binary",
  treatment_var   = "treat",
  outcome_type    = "binary",
  subgroup_vars   = dgm_subgroup_vars,
  subgroup_cuts   = dgm_subgroup_cuts,
  model           = "null",
  n_super         = 5000L,
  seed            = master_seed,
  verbose         = FALSE
)

# ── 10. Parallel plan ───────────────────────────────────────────────────────
# multisession works cross-platform (Linux + macOS).
plan(multisession, workers = n_workers)
on.exit(plan(sequential), add = TRUE)
cat(sprintf("[%s] Using %d parallel workers\n\n", CONFIG_ID, n_workers))

# ── 11. Run simulations ─────────────────────────────────────────────────────
t_total <- proc.time()[3]
bundle <- run_sim_actg175_binary_m1_harm_maxSG_fs1(
  dgm_calibrated    = dgm_calibrated,
  dgm_null          = dgm_null,
  n_alt             = nsims_alt,
  n_null            = nsims_null,
  sim_n             = sim_n,
  confounders       = confounders,
  fs_params         = fs_params,
  grf_params        = grf_params,
  diag_keep         = diag_keep,
  diag_keep_first_n = diag_keep_first_n,
  run_alt           = TRUE,
  run_null          = TRUE,
  verbose           = TRUE,
  config_id         = CONFIG_ID
)
t_total <- proc.time()[3] - t_total

# ── 12. Augment bundle metadata ─────────────────────────────────────────────
git_sha <- tryCatch({
  out <- suppressWarnings(system2("git", c("rev-parse", "HEAD"),
                                  stdout = TRUE, stderr = FALSE))
  if (length(out) >= 1L && nzchar(out[1])) out[1] else NA_character_
}, error = function(e) NA_character_, warning = function(w) NA_character_)

bundle$meta$git_sha       <- git_sha
bundle$meta$run_host      <- Sys.info()[["nodename"]]
bundle$meta$script_path   <- normalizePath(
  file.path(CONFIG_DIR, paste0(CONFIG_ID, ".R")), mustWork = FALSE
)
bundle$meta$total_elapsed_sec <- unname(t_total)
bundle$meta$master_seed   <- master_seed
bundle$meta$cal_target_or <- cal_target_or

# ── 13. Persist ─────────────────────────────────────────────────────────────
saveRDS(bundle, OUT_PATH)
cat(sprintf("\n[%s] Bundle saved to: %s\n", CONFIG_ID, OUT_PATH))
cat(sprintf("[%s] Total elapsed: %.1f min\n", CONFIG_ID, t_total / 60))
cat(sprintf("[%s] Failures: H1 = %s, H0 = %s\n",
            CONFIG_ID,
            format(bundle$meta$n_failed_alt),
            format(bundle$meta$n_failed_null)))
