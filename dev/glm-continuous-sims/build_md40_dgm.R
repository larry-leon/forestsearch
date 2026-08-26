# =============================================================================
# build_md40_dgm.R -- deterministic reconstruction of the maxeffCons MD40 DGM
# =============================================================================
#
# WHY THIS EXISTS
#   The DGM for the continuous/MD workstream was never persisted.  The batch
#   documents build it inline at render time:
#     quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_md40_knoise0_n700_batch_1_1000.qmd
#     chunk `build-dgm`, lines 282-336.
#   A full screen of 3,199 .rds files under ~/ and the repo on 2026-08-26 found
#   zero objects of class glm_dgm.  There is nothing to recover, so the fixture
#   is rebuilt instead.
#
# WHY REBUILDING IS SAFE HERE
#   The construction is a deterministic function of (ACTG175, seed, parameters):
#     - generate_glm_dgm() calls set.seed(seed) before its only sample() call
#       (R/generate_glm_dgm.R:395-396), so df_super is fixed.
#     - No other RNG is touched on this path.
#   The calibrator is BYPASSED.  calibrate_glm_interaction() solves for k_inter
#   by uniroot and then makes one final generate_glm_dgm() call with the root
#   (R/calibrate_glm_interaction.R:140-190).  That root is committed in both
#   payloads as truth$beta_inter = -13.744764123964107 at full double precision
#   (beta_inter == k_inter for continuous outcomes, R/generate_glm_dgm.R:368).
#   Passing it directly reproduces the identical final call, with no root-find
#   in the path.
#
#   Note further that the quantity this fixture is being rebuilt FOR does not
#   depend on k_inter at all: sigma is stats::sigma(fit_base) (L489-493) and
#   mu0 is predict(fit_base, treat = 0) (L454), both taken from the base fit,
#   which is computed at L344 before any interaction shift is applied.  So
#   B(Q) = sigma^2 + V_Q[mu_0] is invariant to the calibration entirely.
#
# GATE
#   The rebuilt object must reproduce all four committed truth constants to
#   1e-9.  These come from the payloads, i.e. from committed data, not from a
#   re-derivation.  If any fails, nothing is written.
#
# COMPUTE
#   One GLM fit on ~1,083 rows plus a 5,000-row prediction.  Sub-second.  No
#   uniroot, no simulation, no parallelism.
#
# USAGE
#   Rscript dev/glm-continuous-sims/build_md40_dgm.R <path-to-a-payload.rds> [outpath]
#
#   The payload argument is the authority for beta_inter and for the four gate
#   constants.  Either the n=500 or the n=700 payload works; both carry
#   identical truth blocks.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("usage: Rscript build_md40_dgm.R <payload.rds> [outpath]", call. = FALSE)
}
payload_path <- args[[1L]]
out_path <- if (length(args) >= 2L) args[[2L]] else
  path.expand("~/actg175/continuous/mr_md_harm/dgm_md40_maxeffCons_seed8316951.rds")

if (!file.exists(payload_path)) stop("no such payload: ", payload_path, call. = FALSE)
if (!requireNamespace("speff2trial", quietly = TRUE))
  stop("package 'speff2trial' is required (provides ACTG175).", call. = FALSE)
suppressPackageStartupMessages(library(forestsearch))

cat("=====================================================================\n")
cat("Deterministic rebuild of the maxeffCons MD40 DGM\n")
cat("=====================================================================\n")
cat("payload : ", normalizePath(payload_path), "\n", sep = "")
cat("out     : ", out_path, "\n", sep = "")
cat("pkg ver : ", as.character(utils::packageVersion("forestsearch")), "\n\n", sep = "")


# -- 1. Authority: the committed truth block ---------------------------------

pl    <- readRDS(payload_path)
truth <- pl$truth

cat("--- 1. committed truth (from payload) ---\n")
for (nm in c("effect_Q", "effect_Qc", "beta_inter", "prevalence_Q", "effect_ITT"))
  cat(sprintf("  %-13s %s\n", nm, format(truth[[nm]], digits = 17)))
cat(sprintf("  meta: n_sample %d  sg_focus %s  seed_base %d  target_md_harm %s\n\n",
            pl$meta$n_sample, pl$meta$sg_focus, pl$meta$seed_base,
            format(pl$meta$target_md_harm)))


# -- 2. Rebuild actg_df, transcribed verbatim --------------------------------
#    Source: sim_fs_maxeffCons_mr_md40_knoise0_n700_batch_1_1000.qmd:282-303
#    Order is load-bearing: the cd420 NA filter runs BEFORE z1..z12, so the
#    medians defining z4/z5/z6 are taken on the filtered frame.

actg_arms        <- c(1L, 3L)
actg_treat_arm   <- 1L
actg_age_cut     <- 34
actg_preanti_cut <- 744.5

actg_df <- subset(speff2trial::ACTG175, arms %in% actg_arms)
actg_df$id <- seq_len(nrow(actg_df))
# TREATMENT IS SWITCHED (ddI = 1).  The -40 calibration is defined against THIS
# coding; its negation would calibrate a different cell under the same name.
actg_df$treat_orig <- ifelse(actg_df$arms == actg_treat_arm, 1L, 0L)
actg_df$treat      <- 1L - actg_df$treat_orig
actg_df$cd4_change <- actg_df$cd420 - actg_df$cd40
actg_df <- actg_df[!is.na(actg_df$cd420), ]
actg_df$z1  <- as.factor(ifelse(actg_df$age > actg_age_cut, 1L, 0L))
actg_df$z2  <- as.factor(ifelse(actg_df$preanti <= actg_preanti_cut, 1L, 0L))
actg_df$z3  <- as.factor(ifelse(actg_df$wtkg <= 75, 1L, 0L))
actg_df$z4  <- as.factor(ifelse(actg_df$karnof <= median(actg_df$karnof), 1L, 0L))
actg_df$z5  <- as.factor(ifelse(actg_df$cd40 <= median(actg_df$cd40), 1L, 0L))
actg_df$z6  <- as.factor(ifelse(actg_df$cd80 <= median(actg_df$cd80), 1L, 0L))
actg_df$z7  <- as.factor(actg_df$hemo);  actg_df$z8  <- as.factor(actg_df$homo)
actg_df$z9  <- as.factor(actg_df$drugs); actg_df$z10 <- as.factor(actg_df$race)
actg_df$z11 <- as.factor(actg_df$gender);actg_df$z12 <- as.factor(actg_df$symptom)

cat("--- 2. source frame ---\n")
cat(sprintf("  rows after arms %s and cd420 NA filter : %d\n",
            paste(actg_arms, collapse = "/"), nrow(actg_df)))
cat(sprintf("  treat coding: ddI = 1, table = %s\n",
            paste(names(table(actg_df$treat)), table(actg_df$treat),
                  sep = ":", collapse = "  ")))
cat(sprintf("  median karnof %.1f  cd40 %.1f  cd80 %.1f\n\n",
            median(actg_df$karnof), median(actg_df$cd40), median(actg_df$cd80)))


# -- 3. Build, with the calibrator bypassed ----------------------------------
#    Reproduces calibrate_glm_interaction()'s final call using the committed
#    root.  Every other argument is transcribed from the .qmd (L321-336).

dgm <- generate_glm_dgm(
  data           = actg_df,
  factor_vars    = paste0("z", 1:12),
  outcome_var    = "cd4_change",
  treatment_var  = "treat",
  outcome_type   = "continuous",
  effect_measure = "MD",
  subgroup_vars  = c("z1", "z2"),
  subgroup_cuts  = list(z1 = 1L, z2 = 1L),
  model          = "alt",
  k_treat        = 1,
  k_inter        = truth$beta_inter,     # the committed uniroot root
  n_super        = 5000L,
  seed           = 8316951L,
  verbose        = FALSE
)


# -- 4. Gate against committed truth -----------------------------------------

cat("--- 3. gate: rebuilt vs committed truth ---\n")
gate <- list(
  effect_Q     = c(dgm$hazard_ratios$harm_subgroup,    truth$effect_Q),
  effect_Qc    = c(dgm$hazard_ratios$no_harm_subgroup, truth$effect_Qc),
  effect_ITT   = c(dgm$hazard_ratios$overall,          truth$effect_ITT),
  prevalence_Q = c(dgm$subgroup_info$proportion,       truth$prevalence_Q),
  beta_inter   = c(dgm$model_params$beta_inter,        truth$beta_inter)
)

ok <- TRUE
for (nm in names(gate)) {
  got <- gate[[nm]][1]; want <- gate[[nm]][2]
  d <- abs(got - want)
  pass <- d < 1e-9
  ok <- ok && pass
  cat(sprintf("  [%s] %-13s rebuilt %20.12f   committed %20.12f   |diff| %.3e\n",
              if (pass) "PASS" else "FAIL", nm, got, want, d))
}
cat(sprintf("\n  |Q| rebuilt %d   (committed prevalence x 5000 = %.1f)\n",
            dgm$subgroup_info$size, truth$prevalence_Q * 5000))
cat(sprintf("  sigma (residual SD) : %.9f\n", dgm$model_params$sigma))
cat(sprintf("  noise_scheme        : %s\n\n", dgm$noise_scheme))

if (!ok) {
  cat("GATE FAILED -- nothing written.\n")
  cat("The rebuild does not reproduce the committed fixture.  Report verbatim\n")
  cat("and stop; do not loosen the tolerance and do not write the object.\n")
  quit(status = 1L)
}


# -- 5. Write -----------------------------------------------------------------

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(dgm, out_path)

cat("--- 4. written ---\n")
cat(sprintf("  path  : %s\n", out_path))
cat(sprintf("  bytes : %d\n", file.info(out_path)$size))
cat("\n  NOTE: this path is OUTSIDE the repository by design.  Do not git add\n")
cat("  the .rds; the reproducible artifact is this script, not the object.\n\n")
cat("=====================================================================\n")
cat("end\n")
cat("=====================================================================\n")
