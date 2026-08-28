# =============================================================================
# fixture_run_fs_oc_predict_2026-08-28.R -- fs_oc_predict() on the MD40 fixture
# -----------------------------------------------------------------------------
# Read-only apart from its printed output.  Rebuilds the maxeffCons MD40 DGM
# exactly as build_md40_dgm.R does (gated against the tracked n = 500 payload's
# committed truth), enumerates the family under the driver's forestsearch
# arguments (as in REPORT_oc_wrapper_build_2026-08-28.md, M = 1601), counts
# what the rmin and n.min floors remove on their own, runs fs_oc_predict() at
# n = 500 with the driver's c1/c2, draws = 2e5, seed = 20260825, and prints a
# three-column comparison whose two comparison columns are READ FROM THE
# PREDICTION DOCUMENT by evaluating its own chunks (worked-predictions and
# appendix-n500), not transcribed.
#
# Run from the repository root:
#   Rscript dev/glm-continuous-sims/fixture_run_fs_oc_predict_2026-08-28.R
# =============================================================================

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library(mvtnorm)
})
stopifnot(requireNamespace("speff2trial", quietly = TRUE))

PAY <- file.path("quarto", "simulations", "actg175", "continuous", "mr_md_harm",
                 "fs_maxeffCons_mr_md40_knoise0_n500_s1000_d5000",
                 "fs_maxeffCons_mr_md40_knoise0_n500_res_1_1000.rds")
DOC <- file.path("quarto", "simulations", "actg175", "continuous",
                 "analytic_verification_and_prediction_md_harm.qmd")
stopifnot(file.exists(PAY), file.exists(DOC))

# ---- 1. the fixture, rebuilt as build_md40_dgm.R builds it ------------------
pl    <- readRDS(PAY)
truth <- pl$truth
actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
actg_df$id <- seq_len(nrow(actg_df))
actg_df$treat_orig <- ifelse(actg_df$arms == 1L, 1L, 0L)
actg_df$treat      <- 1L - actg_df$treat_orig
actg_df$cd4_change <- actg_df$cd420 - actg_df$cd40
actg_df <- actg_df[!is.na(actg_df$cd420), ]
actg_df$z1  <- as.factor(ifelse(actg_df$age > 34, 1L, 0L))
actg_df$z2  <- as.factor(ifelse(actg_df$preanti <= 744.5, 1L, 0L))
actg_df$z3  <- as.factor(ifelse(actg_df$wtkg <= 75, 1L, 0L))
actg_df$z4  <- as.factor(ifelse(actg_df$karnof <= median(actg_df$karnof), 1L, 0L))
actg_df$z5  <- as.factor(ifelse(actg_df$cd40 <= median(actg_df$cd40), 1L, 0L))
actg_df$z6  <- as.factor(ifelse(actg_df$cd80 <= median(actg_df$cd80), 1L, 0L))
actg_df$z7  <- as.factor(actg_df$hemo);  actg_df$z8  <- as.factor(actg_df$homo)
actg_df$z9  <- as.factor(actg_df$drugs); actg_df$z10 <- as.factor(actg_df$race)
actg_df$z11 <- as.factor(actg_df$gender);actg_df$z12 <- as.factor(actg_df$symptom)
dgm <- generate_glm_dgm(
  data = actg_df, factor_vars = paste0("z", 1:12), outcome_var = "cd4_change",
  treatment_var = "treat", outcome_type = "continuous", effect_measure = "MD",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "alt", k_treat = 1, k_inter = truth$beta_inter,
  n_super = 5000L, seed = 8316951L, verbose = FALSE)
gate <- c(effect_Q     = abs(dgm$hazard_ratios$harm_subgroup    - truth$effect_Q),
          effect_Qc    = abs(dgm$hazard_ratios$no_harm_subgroup - truth$effect_Qc),
          prevalence_Q = abs(dgm$subgroup_info$proportion       - truth$prevalence_Q),
          beta_inter   = abs(dgm$model_params$beta_inter        - truth$beta_inter))
cat("fixture rebuild gate (|diff| vs committed truth):\n"); print(gate)
stopifnot(all(gate < 1e-9))
sc <- fs_dgm_scale(dgm)
cat("fs_dgm_scale regions identical to payload $scale$regions:",
    isTRUE(all.equal(sc$regions, pl$scale$regions)), "\n\n")

# ---- 2. the driver's forestsearch arguments ---------------------------------
fs_args <- list(
  confounders.name = c("age", "preanti", "wtkg", "karnof", "cd40", "cd80",
                       "hemo", "homo", "drugs", "race", "gender", "symptom"),
  conf.cont_jcuts = list(age = 10, preanti = 10),
  n.min = 60L, maxk = 2L, sg_focus = "maxeffCons",
  effect.threshold = 30, consistency.threshold = 10)
n_trial <- 500

fam <- fs_oc_family_enumerate(dgm, fs_args, n = n_trial, max_M = 5000L,
                              verbose = TRUE)
print(fam)

# ---- 3. rmin alone / n.min alone over the 2628 combinations -----------------
# Re-walk the enumeration with the wrapper's own pieces so that each floor is
# evaluated on EVERY combination, independent of the others.
df_cut <- dgm$df_super; df_cut$.fs_oc_y <- 0; df_cut$.fs_oc_event <- 1
FSdata <- get_FSdata(df.analysis = df_cut, use_lasso = FALSE, use_grf = FALSE,
                     confounders.name = fs_args$confounders.name,
                     conf.cont_jcuts = fs_args$conf.cont_jcuts,
                     outcome.name = ".fs_oc_y", event.name = ".fs_oc_event",
                     details = FALSE, outcome_type = "continuous")
Zdf <- dummy(FSdata$df[, FSdata$confs_names, drop = FALSE])
Z <- as.matrix(Zdf); storage.mode(Z) <- "numeric"
L <- ncol(Z); maxk <- 2L
combo <- generate_combination_indices(L, maxk)
tot   <- calculate_max_combinations(L, maxk)
rmin <- eval(formals(subgroup.search)$rmin); n.min <- fs_args$n.min
flag <- data.frame(kk = seq_len(tot), k = 0L, empty = FALSE, minp = FALSE,
                   rmin_prop = FALSE, rmin_pop = FALSE, size = FALSE)
for (kk in seq_len(tot)) {
  covs.in <- get_covs_in(kk, maxk, L, combo$counts_1, combo$indices_1,
                         combo$counts_2, combo$indices_2,
                         combo$counts_3, combo$indices_3)
  sel <- which(covs.in == 1); x <- Z[, sel, drop = FALSE]
  flag$k[kk]         <- length(sel)
  flag$empty[kk]     <- !has_positive_variance(x)
  flag$minp[kk]      <- !meets_prevalence_threshold(x, 0.025)
  flag$rmin_prop[kk] <- forestsearch:::.fs_oc_redundant(x, rmin / n_trial)  # wrapper's rule
  flag$rmin_pop[kk]  <- extract_idx_flagredundancy(x, rmin)$flag.redundant # raw count on df_super
  flag$size[kk]      <- mean(get_subgroup_membership(Z, covs.in)) < n.min / n_trial
}
stopifnot(nrow(flag) == 2628L)
cat("\n--- floors evaluated independently on all", tot, "combinations ---\n")
cat(sprintf("rmin alone (wrapper: shrink <= rmin/n = %g in proportion) : %d  [singles %d, pairs %d]\n",
            rmin / n_trial, sum(flag$rmin_prop),
            sum(flag$rmin_prop & flag$k == 1), sum(flag$rmin_prop & flag$k == 2)))
cat(sprintf("  of which also empty-intersection                       : %d\n",
            sum(flag$rmin_prop & flag$empty)))
cat(sprintf("  of which also below n.min/n                            : %d\n",
            sum(flag$rmin_prop & flag$size)))
cat(sprintf("rmin alone (raw subject count on df_super, rmin = %g)      : %d\n",
            rmin, sum(flag$rmin_pop)))
cat(sprintf("n.min alone (Pg < %g/%g = %.3f)                         : %d  [singles %d, pairs %d]\n",
            n.min, n_trial, n.min / n_trial, sum(flag$size),
            sum(flag$size & flag$k == 1), sum(flag$size & flag$k == 2)))
cat(sprintf("  of which also empty-intersection                       : %d\n",
            sum(flag$size & flag$empty)))
cat(sprintf("empty-intersection alone                                 : %d\n", sum(flag$empty)))
cat(sprintf("minp alone                                               : %d\n", sum(flag$minp)))
cat(sprintf("sequential (empty -> minp -> rmin -> size), as the wrapper: %d / %d / %d / %d; kept %d\n",
            sum(flag$empty), sum(!flag$empty & flag$minp),
            sum(!flag$empty & !flag$minp & flag$rmin_prop),
            sum(!flag$empty & !flag$minp & !flag$rmin_prop & flag$size),
            sum(!flag$empty & !flag$minp & !flag$rmin_prop & !flag$size)))
stopifnot(sum(!flag$empty & !flag$minp & !flag$rmin_prop & !flag$size) == fam$counts[["kept"]])

# ---- 4. the prediction ------------------------------------------------------
SEED <- 20260825L; DRAWS <- 2e5
t0 <- proc.time()[3]
pred <- fs_oc_predict(family = fam, n = n_trial,
                      c1 = fs_args$effect.threshold,
                      c2 = fs_args$consistency.threshold,
                      consistency_method = "split", draws = DRAWS, seed = SEED)
cat(sprintf("\nfs_oc_predict on M = %d: %.0f s\n", fam$M, proc.time()[3] - t0))
print(pred)

# ---- 5. the document's own figures, from its own chunks ---------------------
.chunk_code <- function(lines, label) {
  open  <- grep(sprintf("^```\\{r %s[ ,}]", label), lines)
  if (length(open) != 1L) stop("chunk '", label, "' not found exactly once")
  close <- grep("^```\\s*$", lines); close <- close[close > open][1L]
  lines[(open + 1L):(close - 1L)]
}
doc_lines <- readLines(DOC, warn = FALSE)
old_wd <- setwd(dirname(DOC))
env <- new.env(parent = globalenv())
invisible(capture.output({
  for (ch in c("anchor", "worked-scenario", "worked-predictions", "appendix-n500"))
    eval(parse(text = .chunk_code(doc_lines, ch)), envir = env)
}))
setwd(old_wd)
stopifnot(env$M == 16L, env$Rdraw == 2e5)

# The document's own row mapping for the measured column (appendix-n500):
# detection, nH, sens, spec, ppv, npv, EbetaH = oracle - oracle_bias
# (bH500["oracle"]), naive_bias.
rows <- c("det_rate", "EnH", "Esens", "Espec", "Eppv", "Enpv", "EbetaH", "Enaive_bias")
doc_pred <- unlist(mget(rows, envir = env))
measured <- c(env$m500$detection, env$m500$nH, env$m500$sens, env$m500$spec,
              env$m500$ppv, env$m500$npv, unname(env$bH500["oracle"]),
              env$m500$naive_bias)
wrap <- unlist(pred[rows])

tab <- data.frame(
  quantity = c(rows, "M"),
  wrapper_M1601_split = c(round(wrap, 3), fam$M),
  document_M16_split  = c(round(doc_pred, 3), env$M),
  measured_n500_resample = c(measured, NA),
  check.names = FALSE)
cat("\n=== fs_oc_predict on the MD40 fixture vs the document vs the n = 500 record ===\n")
cat(sprintf("wrapper: n = %d, c1 = %g, c2 = %g, gate = split, draws = %g, seed = %d, MC SE(det_rate) = %.4f\n",
            n_trial, fs_args$effect.threshold, fs_args$consistency.threshold,
            DRAWS, SEED, pred$det_rate_se))
cat("document: worked-predictions chunk (M = 16, split, Rdraw = 2e5, seed 20260825)\n")
cat("measured: appendix-n500 m500 block (driver: consistency_method = 'resample', pcons 0.90, fs.splits 400)\n\n")
print(tab, row.names = FALSE)

cat("\nwrapper top-8 selections (label, P_selected, sel | det, Pg, PQg, beta_g):\n")
top <- order(-pred$p_sel)[1:8]
print(data.frame(lab = pred$lab[top], p_sel = round(pred$p_sel[top], 4),
                 sel_c = round(pred$sel_c[top], 4), Pg = round(fam$Pg[top], 3),
                 PQg = round(fam$PQg[top], 3), beta = round(fam$beta_g[top], 1)),
      row.names = FALSE)
cat(sprintf("\nselection mass on rules with true mean below c1: %.3f;  max P1 = %.3f;  candidates with P1 > 0.5: %d\n",
            pred$mass_below, max(pred$P1), sum(pred$P1 > 0.5)))
