# =============================================================================
# diagnose_md_harm_readout_c_and_mr_skip.R
#
# Follow-up to diagnose_md_harm_pilot_zero_detection.R, correcting one defect
# in that script and capturing one finding it surfaced.
#
# 1. READOUT (c) WAS WRONG in the first script. It searched for "age"/"preanti"
#    in the `grp` column of find.grps$out.found$hr.subgroups. `grp` is the
#    integer combination index, not a label -- format_search_results()
#    (R/subgroup_search.R:877-879) names the columns
#      c("grp","K","n","E","d1","m1","m0","HR","L(HR)","U(HR)", colnames(Z))
#    so candidate membership lives in the trailing q1..qL indicator columns.
#    The first script's "candidates naming age or preanti : 0" was therefore an
#    ARTEFACT OF THE EXTRACTION, not a finding. Corrected here.
#
# 2. RUN B returned with `mr_inference` NULL even though the MR guard passed.
#    This captures the verbatim skip reason by re-running with quiet = FALSE.
#
# INVESTIGATION ONLY. No harness edit, no configuration fix applied, no R/
# change. Flags-off probe exactly as authorized: use_lasso = FALSE,
# use_grf = FALSE, use_dina = FALSE, everything else byte-identical to the
# harness -- except `quiet`, which is a verbosity flag and is set FALSE here
# solely to make the skip reason audible.
#
# Run from the package root, against the INSTALLED package:
#   Rscript dev/glm-continuous-sims/verification/diagnose_md_harm_readout_c_and_mr_skip.R
# =============================================================================

library(forestsearch)
library(speff2trial)

`%||%` <- function(a, b) if (is.null(a)) b else a
.rule <- function(txt) cat("\n", strrep("=", 78), "\n", txt, "\n",
                           strrep("=", 78), "\n", sep = "")

OUT_RDS <- file.path(Sys.getenv("SCRATCH", tempdir()), "runB_replicate1.rds")

## --- harness constants, transcribed from the qmd ---------------------------
actg_arms        <- c(1L, 3L)
actg_treat_arm   <- 1L
actg_age_cut     <- 34
actg_preanti_cut <- 744.5
dgm_factor_vars   <- paste0("z", 1:12)
dgm_subgroup_vars <- c("z1", "z2")
dgm_subgroup_cuts <- list(z1 = 1L, z2 = 1L)
dgm_n_super       <- 5000L
cal_target_md     <- -40
cal_k_grid_range  <- c(0, 120)
cal_grid_step     <- 2

md_threshold    <- 30
md_consistency  <- 10
adverse_outcome <- FALSE
pconsistency    <- 0.90
fs_splits <- 400L; maxk <- 2L; n_min <- 60L; d0_min <- 12L; d1_min <- 12L
vi_grf_min <- -0.2
use_twostage <- TRUE; is_rct <- TRUE
mr_draws <- 5000L

analysis_continuous_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
analysis_binary_vars     <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
confounders_analysis     <- c(analysis_continuous_vars, analysis_binary_vars)

outcome_name <- "y_sim"; treat_name <- "treat_sim"
id_name <- "id"; harm_col <- "flag_harm"

seed_base <- 8316951L
MAX_SIMS  <- 5000L
set.seed(seed_base)
SEED_TABLE <- sample.int(.Machine$integer.max - 1L, MAX_SIMS)
sd_1 <- SEED_TABLE[1L]

## --- DGM, exactly as the qmd builds it -------------------------------------
actg_df <- subset(speff2trial::ACTG175, arms %in% actg_arms)
actg_df$id <- seq_len(nrow(actg_df))
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
actg_df$z7  <- as.factor(actg_df$hemo);   actg_df$z8  <- as.factor(actg_df$homo)
actg_df$z9  <- as.factor(actg_df$drugs);  actg_df$z10 <- as.factor(actg_df$race)
actg_df$z11 <- as.factor(actg_df$gender); actg_df$z12 <- as.factor(actg_df$symptom)
for (v in analysis_binary_vars) actg_df[[v]] <- as.factor(actg_df[[v]])

dgm <- calibrate_glm_interaction(
  data = actg_df, factor_vars = dgm_factor_vars,
  outcome_var = "cd4_change", treatment_var = "treat",
  target_effect = cal_target_md, outcome_type = "continuous",
  effect_measure = "MD", subgroup_vars = dgm_subgroup_vars,
  subgroup_cuts = dgm_subgroup_cuts, k_inter_range = cal_k_grid_range,
  grid_step = cal_grid_step, n_super = dgm_n_super, seed = seed_base,
  verbose = FALSE)

df1 <- simulate_from_glm_dgm(dgm, n = 1000L, seed = sd_1)
df1[[id_name]] <- seq_len(nrow(df1))

# =============================================================================
# RUN B', quiet = FALSE, so the MR skip reason is audible
# =============================================================================
.rule("RUN B' -- flags off, quiet = FALSE (to capture the MR skip reason)")

runB <- forestsearch(
  df.analysis = df1, confounders.name = confounders_analysis,
  outcome.name = outcome_name, treat.name = treat_name, id.name = id_name,
  outcome_type = "continuous", effect_measure = "MD",
  effect.threshold = md_threshold, consistency.threshold = md_consistency,
  pconsistency.threshold = pconsistency, fs.splits = fs_splits,
  n.min = n_min, d0.min = d0_min, d1.min = d1_min, maxk = maxk,
  vi.grf.min = vi_grf_min,
  use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
  use_twostage = use_twostage, is.RCT = is_rct,
  adverse_outcome = adverse_outcome,
  details = FALSE, quiet = FALSE, seedit = sd_1,
  parallel_args = list(plan = "sequential"),
  mr_inference = TRUE,
  mr_inference_args = list(ci_method = "ij", draws = mr_draws,
                           include_complement = TRUE))

saveRDS(runB, OUT_RDS)
cat("\nrunB saved to: ", OUT_RDS, "\n", sep = "")

.rule("MR eligibility, resolved")
cat("mr_inference on object      : ", !is.null(runB$mr_inference), "\n", sep = "")
cat("mr_harm_confirmed           : ", format(runB$mr_harm_confirmed), "\n", sep = "")
cat("consistency_method resolved : ",
    paste(format(runB$args_call_all$consistency_method), collapse = ", "),
    "\n", sep = "")
cat("consistency_method FORMAL default (match.arg takes the first) : ",
    paste(format(formals(forestsearch)$consistency_method), collapse = ", "),
    "\n", sep = "")
cat("\nEligibility test at R/forestsearch_main.R:3142-3146:\n")
cat("  .mr_glm_ok <- consistency_method == \"resample\" && !is.null(estimator_fn)\n")
cat("  The harness never passes consistency_method, so it resolves to\n")
cat("  \"split\" and .mr_glm_ok is FALSE. MR is skipped via .mr_skip()\n")
cat("  at R/forestsearch_main.R:3152-3168 -- silently, because the harness\n")
cat("  passes quiet = TRUE.\n")

# =============================================================================
# READOUT (c), CORRECTED
# =============================================================================
.rule("READOUT (c) CORRECTED -- truth-adjacent candidates")

fg  <- runB$find.grps
tbl <- as.data.frame(fg$out.found$hr.subgroups)
labs <- runB$confounders.evaluated
qs   <- runB$confounders.candidate

cat("qualifying candidates : ", nrow(tbl), "\n", sep = "")
cat("indicator columns     : ", length(qs), "\n", sep = "")

# Map each row's indicator columns to the human-readable cut labels.
qcols <- intersect(qs, names(tbl))
stopifnot(length(qcols) == length(qs))
lab_of <- stats::setNames(labs, qs)

row_label <- function(i) {
  on <- qcols[which(as.numeric(tbl[i, qcols]) == 1)]
  if (!length(on)) return(NA_character_)
  paste(lab_of[on], collapse = " & ")
}
tbl$label <- vapply(seq_len(nrow(tbl)), row_label, character(1))

cat("\ntop 15 qualifying candidates by ORIENTED effect (positive = harm):\n")
o <- order(-tbl$HR)
print(tbl[o[1:15], c("label", "K", "n", "HR", "L(HR)", "U(HR)")],
      row.names = FALSE)

## --- candidates naming the true-boundary variables -------------------------
age_labs <- grep("^age ", labs, value = TRUE)
pre_labs <- grep("^preanti ", labs, value = TRUE)
cat("\nage cut labels     : ", paste(age_labs, collapse = " | "), "\n", sep = "")
cat("preanti cut labels : ", paste(pre_labs, collapse = " | "), "\n", sep = "")

has_var <- function(v) grepl(paste0("(^|& )", v, " "), tbl$label)
sel_age <- has_var("age"); sel_pre <- has_var("preanti")
cat("\nqualifying candidates naming age     : ", sum(sel_age), "\n", sep = "")
cat("qualifying candidates naming preanti : ", sum(sel_pre), "\n", sep = "")
cat("qualifying candidates naming both    : ", sum(sel_age & sel_pre), "\n", sep = "")

## --- the truth-adjacent candidates -----------------------------------------
# True region Q is {age > 34} AND {preanti <= 744.5}.  On the enumerated grid
# the nearest expressions are the NEGATION of "age <= 34" and the cut
# "preanti <= 770" (the nearest cut above 744.5).  The search's dummy()
# expansion carries both directions of each cut, so both are representable.
cat("\n--- truth-adjacent candidates (naming age and/or preanti) ---\n")
ta <- tbl[sel_age | sel_pre, , drop = FALSE]
if (!nrow(ta)) {
  cat("NONE qualified.\n")
} else {
  ta <- ta[order(-ta$HR), ]
  cat("count : ", nrow(ta), "   oriented-effect range : [",
      sprintf("%.4f, %.4f", min(ta$HR), max(ta$HR)), "]\n", sep = "")
  print(utils::head(ta[, c("label", "K", "n", "HR", "L(HR)", "U(HR)")], 25L),
        row.names = FALSE)
}

cat("\nthreshold compared against (oriented scale, MD >= x) : ",
    md_threshold, "\n", sep = "")
cat("selected subgroup (sg.harm) : ",
    paste(runB$sg.harm, collapse = " & "), "\n", sep = "")

.rule("END -- investigation only; no configuration change was applied")
