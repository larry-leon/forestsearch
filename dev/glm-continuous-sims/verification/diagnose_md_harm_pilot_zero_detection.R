# =============================================================================
# diagnose_md_harm_pilot_zero_detection.R
#
# Blocker 3 of dev/glm-continuous-sims/HANDOFF_md_mr_harness_session.md (v2):
# diagnose "detected 0 of 50" in the pilot bundle
#   quarto/simulations/actg175/continuous/mr_sweep_md_harm/
#     md_harm_s50_pilot/fs_md_harm_n1000_res.rds
#
# INVESTIGATION ONLY.  This script changes nothing: it does not edit the
# harness, does not apply a configuration fix, and does not touch R/.  It
# reproduces ONE replicate (global replicate id 1, under the harness's own
# pre-generated seed) and prints the readouts blocker 3 asks for.
#
# It runs the replicate TWICE, deliberately:
#
#   RUN A -- the harness configuration verbatim, including
#            use_lasso = TRUE, use_grf = TRUE, mr_inference = TRUE.
#            This is the configuration that produced the 0/50 bundle.
#
#   RUN B -- the same replicate with use_lasso = FALSE, use_grf = FALSE
#            (the MR-COMPATIBLE configuration the handoff adopts in
#            principle).  RUN B exists only so that the candidate-family,
#            cut-grid, threshold and funnel readouts have something to be
#            read off; it is a diagnostic probe, NOT a proposed harness edit
#            and NOT a re-pilot.
#
# Requires the INSTALLED package (devtools::install()), matching the harness.
#
# Run from the package root:
#   Rscript dev/glm-continuous-sims/verification/diagnose_md_harm_pilot_zero_detection.R
# =============================================================================

library(forestsearch)
library(speff2trial)

`%||%` <- function(a, b) if (is.null(a)) b else a

.rule <- function(txt) cat("\n", strrep("=", 78), "\n", txt, "\n",
                           strrep("=", 78), "\n", sep = "")

BUNDLE <- file.path("quarto", "simulations", "actg175", "continuous",
                    "mr_sweep_md_harm", "md_harm_s50_pilot",
                    "fs_md_harm_n1000_res.rds")

# =============================================================================
# STEP 0 -- the object being diagnosed
# =============================================================================
.rule("STEP 0 -- the pilot bundle")

stopifnot(file.exists(BUNDLE))
cat("absolute path : ", normalizePath(BUNDLE), "\n", sep = "")
cat("size (bytes)  : ", file.size(BUNDLE), "\n", sep = "")
cat("mtime         : ", format(file.mtime(BUNDLE), "%Y-%m-%d %H:%M:%OS3 %z"),
    "\n", sep = "")

bnd <- readRDS(BUNDLE)
res <- bnd$results

cat("\n-- replicate counts --\n")
cat("n replicates  : ", nrow(res), "   (expected 50)\n", sep = "")
cat("sum(detected) : ", sum(res$detected), "   (expected 0)\n", sep = "")
cat("sum(mr_ok)    : ", sum(res$mr_ok), "\n", sep = "")

cat("\n-- per-replicate status tabulation --\n")
cat("detected:\n");        print(table(res$detected,        useNA = "ifany"))
cat("mr_ok:\n");           print(table(res$mr_ok,           useNA = "ifany"))
cat("betaHhat_status:\n"); print(table(res$betaHhat_status, useNA = "ifany"))
cat("nH_eval:\n");         print(table(res$nH_eval,         useNA = "ifany"))
cat("nHc_eval:\n");        print(table(res$nHc_eval,        useNA = "ifany"))
cat("label all NA  : ", all(is.na(res$label)), "\n", sep = "")
cat("sg_def all NA : ", all(is.na(res$sg_def)), "\n", sep = "")

cat("\n-- betaHhat_counts attribute --\n")
print(attr(res, "betaHhat_counts"))

cat("\n-- per-replicate wall clock inside forestsearch() (t2_secs) --\n")
print(summary(res$t2_secs))
cat("NOTE: t2_secs is measured around the forestsearch() call itself.\n")

cat("\n-- meta signature --\n")
for (k in names(bnd$meta)) {
  v <- bnd$meta[[k]]
  cat(sprintf("  %-22s : %s\n", k,
              if (is.null(v)) "NULL" else paste(format(v), collapse = ", ")))
}

cat("\n-- truth block --\n")
for (k in names(bnd$truth))
  cat(sprintf("  %-14s : %.10f\n", k, bnd$truth[[k]]))

# =============================================================================
# STEP 1 -- reproduce replicate 1 under the harness's own configuration
#           and the harness's own seed for GLOBAL replicate id 1
# =============================================================================
.rule("STEP 1 -- rebuild the harness's DGM and replicate 1")

## --- the harness's constants, transcribed from the qmd (do not re-derive) ---
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
use_lasso <- TRUE; use_grf <- TRUE; use_twostage <- TRUE; is_rct <- TRUE

analysis_continuous_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
analysis_binary_vars     <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
confounders_analysis     <- c(analysis_continuous_vars, analysis_binary_vars)

outcome_name <- "y_sim"; treat_name <- "treat_sim"
id_name <- "id"; harm_col <- "flag_harm"

mr_draws <- 5000L

## --- the harness's seed scheme: pre-generated, indexed by GLOBAL sim_id -----
seed_base <- 8316951L
MAX_SIMS  <- 5000L
set.seed(seed_base)
SEED_TABLE <- sample.int(.Machine$integer.max - 1L, MAX_SIMS)
seed_for <- function(sim_id) {
  if (sim_id < 1L || sim_id > MAX_SIMS)
    stop("sim_id ", sim_id, " outside the pre-generated seed table.")
  SEED_TABLE[sim_id]
}

SIM_ID <- 1L
sd_1 <- seed_for(SIM_ID)
cat("seed_for(1)                  : ", sd_1, "\n", sep = "")
cat("bundle results$seed[1]       : ", res$seed[1], "\n", sep = "")
cat("seed reproduction matches    : ", identical(sd_1, res$seed[1]), "\n", sep = "")

## --- the DGM, exactly as the qmd builds it ---------------------------------
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

eval_df <- fs_build_eval_frame(dgm, outcome_type = "continuous")
cat("df_super n                   : ", nrow(eval_df), "\n", sep = "")
cat("P(Q) on df_super             : ",
    sprintf("%.6f", mean(eval_df[[harm_col]] == 1L)), "\n", sep = "")

## --- replicate 1's trial data ----------------------------------------------
df1 <- simulate_from_glm_dgm(dgm, n = 1000L, seed = sd_1)
df1[[id_name]] <- seq_len(nrow(df1))
cat("replicate 1 n                : ", nrow(df1), "\n", sep = "")
cat("replicate 1 n_true (flag_harm): ", sum(df1[[harm_col]] == 1L), "\n", sep = "")
cat("bundle results$n_true[1]     : ", res$n_true[1], "\n", sep = "")
cat("n_true reproduction matches  : ",
    identical(as.integer(sum(df1[[harm_col]] == 1L)),
              as.integer(res$n_true[1])), "\n", sep = "")

# =============================================================================
# RUN A -- the harness configuration VERBATIM
# =============================================================================
.rule("RUN A -- harness configuration verbatim (use_lasso = TRUE, use_grf = TRUE)")

fs_args_common <- list(
  df.analysis = df1, confounders.name = confounders_analysis,
  outcome.name = outcome_name, treat.name = treat_name, id.name = id_name,
  outcome_type = "continuous", effect_measure = "MD",
  effect.threshold = md_threshold, consistency.threshold = md_consistency,
  pconsistency.threshold = pconsistency, fs.splits = fs_splits,
  n.min = n_min, d0.min = d0_min, d1.min = d1_min, maxk = maxk,
  vi.grf.min = vi_grf_min,
  use_twostage = use_twostage, is.RCT = is_rct,
  adverse_outcome = adverse_outcome,
  details = FALSE, quiet = TRUE, seedit = sd_1,
  parallel_args = list(plan = "sequential"),
  mr_inference = TRUE,
  mr_inference_args = list(ci_method = "ij", draws = mr_draws,
                           include_complement = TRUE))

tA <- proc.time()[3]
runA <- tryCatch(
  do.call(forestsearch, c(fs_args_common,
                          list(use_lasso = TRUE, use_grf = TRUE))),
  error = function(e) structure(list(msg = conditionMessage(e)),
                                class = "diag_error"))
secsA <- proc.time()[3] - tA

cat("elapsed (s)                  : ", sprintf("%.4f", secsA), "\n", sep = "")
if (inherits(runA, "diag_error")) {
  cat("RESULT                       : ERROR (forestsearch() stopped)\n")
  cat("--- verbatim error message ---\n")
  cat(runA$msg, "\n")
  cat("--- end error message ---\n\n")
  cat("The harness wraps this call in tryCatch(error = function(e) NULL)\n",
      "(qmd line 290) inside suppressWarnings(), so the stop() is swallowed,\n",
      "fs.est becomes NULL, record_replicate() returns at qmd line 292, and\n",
      "the replicate is recorded with detected = 0.\n", sep = "")
} else {
  cat("RESULT                       : forestsearch() RETURNED\n")
  cat("sg.harm                      : ",
      paste(runA$sg.harm %||% "NULL", collapse = " & "), "\n", sep = "")
  cat("family_status                : ", runA$family_status %||% "NULL", "\n", sep = "")
}

# =============================================================================
# RUN B -- MR-compatible probe (use_lasso = FALSE, use_grf = FALSE)
#          Diagnostic only.  Not a proposed harness edit.
# =============================================================================
.rule("RUN B -- MR-compatible probe (use_lasso = FALSE, use_grf = FALSE)")

tB <- proc.time()[3]
runB <- tryCatch(
  suppressWarnings(do.call(forestsearch,
                           c(fs_args_common,
                             list(use_lasso = FALSE, use_grf = FALSE,
                                  use_dina = FALSE)))),
  error = function(e) structure(list(msg = conditionMessage(e)),
                                class = "diag_error"))
secsB <- proc.time()[3] - tB
cat("elapsed (s)                  : ", sprintf("%.2f", secsB), "\n", sep = "")

if (inherits(runB, "diag_error")) {
  cat("RESULT                       : ERROR\n")
  cat(runB$msg, "\n")
  quit(save = "no", status = 0)
}

cat("RESULT                       : forestsearch() RETURNED\n")
cat("family_status                : ", runB$family_status %||% "NULL", "\n", sep = "")
cat("sg.harm                      : ",
    if (is.null(runB$sg.harm)) "NULL"
    else paste(runB$sg.harm, collapse = " & "), "\n", sep = "")
cat("mr_inference present         : ", !is.null(runB$mr_inference), "\n", sep = "")

# =============================================================================
# STEP 2 -- required readouts (a) - (e), read off RUN B
# =============================================================================

## ---------------------------------------------------------------------------
## (a) candidate family size and the cut grids for age / preanti
## ---------------------------------------------------------------------------
.rule("STEP 2 (a) -- candidate family and the cut grids")

cat("harness arguments reaching get_FSdata():\n")
cat("  cut_type   : ",
    format(runB$args_call_all$cut_type %||% "<not supplied; formal default>"),
    "\n", sep = "")
cat("  conf.cont  : ",
    if (is.null(runB$args_call_all$conf.cont)) "NULL (not supplied)"
    else paste(format(runB$args_call_all$conf.cont), collapse = ", "),
    "\n", sep = "")
cat("  confounders.name : ", paste(confounders_analysis, collapse = ", "),
    "\n", sep = "")
cat("  maxk = ", maxk, "   n.min = ", n_min,
    "   d0.min = ", d0_min, "   d1.min = ", d1_min, "\n", sep = "")

cuts_all <- runB$confounders.evaluated
cat("\nconfounders.evaluated (the enumerated cut labels), n = ",
    length(cuts_all), ":\n", sep = "")
print(cuts_all)

cat("\nconfounders.candidate (internal names), n = ",
    length(runB$confounders.candidate), ":\n", sep = "")
print(runB$confounders.candidate)

.grid_for <- function(var, labels) {
  hit <- grep(paste0("^", var, "\\b|^", var, "[^a-zA-Z0-9_]"), labels, value = TRUE)
  if (!length(hit)) hit <- grep(var, labels, value = TRUE, fixed = TRUE)
  hit
}
for (v in c("age", "preanti")) {
  cat("\ncut grid for '", v, "':\n", sep = "")
  g <- .grid_for(v, cuts_all)
  if (!length(g)) cat("  <no cut on this variable in the evaluated family>\n")
  else print(g)
}

## ---------------------------------------------------------------------------
## (b) bracketing of the true boundaries
## ---------------------------------------------------------------------------
.rule("STEP 2 (b) -- bracketing of age > 34 and preanti <= 744.5")

.numeric_cuts <- function(labels) {
  m <- regmatches(labels, regexpr("-?[0-9]+\\.?[0-9]*$", labels))
  sort(unique(as.numeric(m[nzchar(m)])))
}
for (spec in list(list(v = "age",     b = actg_age_cut),
                  list(v = "preanti", b = actg_preanti_cut))) {
  g  <- .grid_for(spec$v, cuts_all)
  cv <- .numeric_cuts(g)
  cat("\nvariable      : ", spec$v, "\n", sep = "")
  cat("true boundary : ", spec$b, "\n", sep = "")
  cat("cut values    : ", if (length(cv)) paste(cv, collapse = ", ") else "<none>",
      "\n", sep = "")
  below <- cv[cv <= spec$b]; above <- cv[cv > spec$b]
  cat("nearest at or below : ",
      if (length(below)) max(below) else "<none>", "\n", sep = "")
  cat("nearest above       : ",
      if (length(above)) min(above) else "<none>", "\n", sep = "")
  cat("brackets the boundary: ", length(below) > 0 && length(above) > 0,
      "\n", sep = "")
  cat("observed range on replicate 1 : [",
      paste(range(df1[[spec$v]], na.rm = TRUE), collapse = ", "), "]\n", sep = "")
}

## ---------------------------------------------------------------------------
## (c) truth-adjacent candidates and their oriented effects vs the threshold
## ---------------------------------------------------------------------------
.rule("STEP 2 (c) -- truth-adjacent candidates, oriented effects vs the threshold")

fg <- runB$find.grps
cat("find.grps class : ", paste(class(fg), collapse = ","), "\n", sep = "")
if (is.list(fg)) cat("find.grps names : ", paste(names(fg), collapse = ", "),
                     "\n", sep = "")

# subgroup.search() returns out.found = list(hr.subgroups = <data.table>), or
# out.found = NULL when nothing qualified (format_search_results(), :858-868).
.tbl <- if (is.list(fg) && is.list(fg$out.found)) fg$out.found$hr.subgroups else NULL

if (is.null(.tbl) || !nrow(.tbl)) {
  cat("\nout.found is NULL: NO candidate qualified. The surviving set is EMPTY.\n")
} else {
  .tbl <- as.data.frame(.tbl)
  cat("\nqualifying candidates : ", nrow(.tbl), "\n", sep = "")
  cat("columns (first 10) : ",
      paste(utils::head(names(.tbl), 10L), collapse = ", "), "\n", sep = "")
  keep <- intersect(c("grp", "K", "n", "E", "d1", "m1", "m0",
                      "HR", "L(HR)", "U(HR)"), names(.tbl))
  ord <- order(-.tbl$HR)
  cat("\ntop 20 by HR (the ORIENTED effect; MD on the adverse_outcome=FALSE\n",
      "negated scale, so positive = harm):\n", sep = "")
  print(utils::head(.tbl[ord, keep, drop = FALSE], 20L), row.names = FALSE)
  tr <- grepl("age", .tbl$grp) | grepl("preanti", .tbl$grp)
  cat("\ncandidates naming age or preanti : ", sum(tr), "\n", sep = "")
  if (any(tr)) print(.tbl[tr, keep, drop = FALSE], row.names = FALSE)
  cat("\nrange of oriented HR column : [",
      sprintf("%.4f, %.4f", min(.tbl$HR, na.rm = TRUE),
              max(.tbl$HR, na.rm = TRUE)), "]\n", sep = "")
}

cat("\nthreshold compared against (oriented scale) : ", md_threshold, "\n", sep = "")

## ---------------------------------------------------------------------------
## (d) threshold value, sign and orientation on THIS executing path
## ---------------------------------------------------------------------------
.rule("STEP 2 (d) -- threshold value, sign, orientation on this path")

cat("user_set_threshold inputs:\n")
cat("  effect.threshold supplied      : ", md_threshold, "\n", sep = "")
cat("  hr.threshold supplied          : <missing>\n")
cat("  consistency.threshold supplied : ", md_consistency, "\n", sep = "")

cat("\nthreshold_config as returned by forestsearch():\n")
tc <- runB$threshold_config
if (is.null(tc)) cat("  NULL\n") else {
  for (k in names(tc))
    cat(sprintf("  %-24s : %s\n", k,
                if (is.null(tc[[k]])) "NULL"
                else paste(format(tc[[k]]), collapse = ", ")))
}

cat("\nadverse_outcome as executed     : ", runB$args_call_all$adverse_outcome, "\n", sep = "")
cat("effect_measure                  : ", runB$effect_measure, "\n", sep = "")
cat("outcome_type                    : ", runB$outcome_type, "\n", sep = "")
cat("is_identity (no log transform)  : ",
    runB$effect_measure %in% c("RD", "IRD", "MD"), "\n", sep = "")

cat("\nadmission set MR re-selected over:\n")
print(runB$admission)

## ---------------------------------------------------------------------------
## (e) the failure funnel
## ---------------------------------------------------------------------------
.rule("STEP 2 (e) -- the failure funnel")

## Status-code verification: which gate a status code corresponds to, read
## from the source rather than assumed.  The status-6 rejection at
## R/subgroup_search.R:626-630 is the EFFECT FLOOR, and it is conditional on
## disable_effect_floor:
##
##   # Status 6: Check effect threshold.  Skipped when the effect floor is
##   # disabled (sg_focus = "maxeff" retains the full estimable family so the
##   # argmax is unconditional -- see forestsearch() search_overrides).
##   if (!disable_effect_floor && glm_result$hr <= hr.threshold) {
##     return(list(status = 6L, result = NULL))
##   }
##
## So n_passed_hr counts candidates with oriented effect STRICTLY ABOVE
## hr.threshold, and a candidate exactly AT the threshold is rejected.
cat("effect floor comparison in force on this run:\n")
cat("  site            : R/subgroup_search.R:626-630 (GLM path)\n")
cat("  test            : !disable_effect_floor && glm_result$hr <= hr.threshold -> status 6\n")
cat("  hr.threshold    : ", md_threshold, "\n", sep = "")
cat("  sg_focus        : ", runB$sg_focus %||% "<NULL>", "\n", sep = "")
cat("  admission$effect_floor : ",
    if (is.null(runB$admission$effect_floor)) "NULL (floor NOT applied)"
    else format(runB$admission$effect_floor), "\n", sep = "")
cat("  => effect floor applied on this path : ",
    !is.null(runB$admission$effect_floor), "\n\n", sep = "")

fc <- if (is.list(fg)) fg$filter_counts else NULL
if (is.null(fc)) {
  cat("filter_counts not carried on the returned object.\n")
} else {
  # Gate names and their status codes, from evaluate_combination_with_status()
  # (R/subgroup_search.R:553-673) and the accumulator at :276-299.
  gates <- c(
    n_evaluated          = "combinations evaluated              (status >= 0)",
    n_passed_variance    = "passed variance check               (status >= 1)",
    n_passed_prevalence  = "passed prevalence / minp            (status >= 2)",
    n_passed_redundancy  = "passed redundancy pruning           (status >= 3)",
    n_passed_events      = "passed per-arm counts d0.min/d1.min (status >= 4)",
    n_passed_sample_size = "passed sample size n.min            (status >= 5)",
    n_passed_cox         = "model fit succeeded                 (status >= 6)",
    n_passed_hr          = "passed EFFECT FLOOR effect_threshold(status >= 7)")
  cat("filter_counts funnel (subgroup.search / evaluate_combination_with_status):\n")
  prev <- NA_integer_
  for (k in names(gates)) {
    v <- fc[[k]]
    drop <- if (is.na(prev)) "" else sprintf("   (-%d)", prev - v)
    cat(sprintf("  %-52s : %6s%s\n", gates[[k]], format(v), drop))
    prev <- v
  }
  cat("\nfirst gate that empties the qualifying set:\n")
  emptied <- names(gates)[vapply(names(gates), function(k) identical(fc[[k]], 0L),
                                 logical(1))]
  if (!length(emptied)) {
    cat("  none -- the qualifying set is non-empty after every search gate.\n")
  } else {
    k <- emptied[1]
    src <- c(n_passed_variance    = "R/subgroup_search.R:568",
             n_passed_prevalence  = "R/subgroup_search.R:573",
             n_passed_redundancy  = "R/subgroup_search.R:579",
             n_passed_events      = "R/subgroup_search.R:599/606/647",
             n_passed_sample_size = "R/subgroup_search.R:614/653",
             n_passed_cox         = "R/subgroup_search.R:620/661",
             n_passed_hr          = "R/subgroup_search.R:626/667 (effect floor)")
    cat("  ", k, "  at ", src[[k]] %||% "<n/a>", "\n", sep = "")
  }
}

cat("\ngrp.consistency:\n")
gc_ <- runB$grp.consistency
if (is.null(gc_)) cat("  NULL\n") else {
  cat("  names : ", paste(names(gc_), collapse = ", "), "\n", sep = "")
  cat("  sg.harm : ",
      if (is.null(gc_$sg.harm)) "NULL"
      else paste(gc_$sg.harm, collapse = " & "), "\n", sep = "")
  if (!is.null(gc_$sg.harm.id))
    cat("  |sg.harm.id == 1| : ", sum(gc_$sg.harm.id == 1L), "\n", sep = "")
  if (!is.null(gc_$max.consistency))
    cat("  max.consistency : ", gc_$max.consistency, "\n", sep = "")
}

.rule("END -- investigation only; no configuration change was applied")
