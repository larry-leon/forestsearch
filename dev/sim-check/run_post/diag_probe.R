# Diagnostic probe for CC_BRIEF_sim_integrity.md sections 2, 3, 5.
# Reproduces the post-edit qmd's setup + one replicate (sim_id = 1), then:
#   §2 resolved parameters in fs.est$args_call_all
#   §5 field paths for fb_* / mr_*
#   §3 FB bootstrap with mr_in_replicates = FALSE vs TRUE
# Run from dev/sim-check/run_post (needs betaHhat_truth.R in wd).

library(forestsearch); library(survival); library(data.table)
library(foreach); library(doFuture); library(future)
source("betaHhat_truth.R")

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

## ---- knobs copied verbatim from the post-edit setup chunk -------------------
parallel_mode <- "boots"
nb_boots <- 20L; k_random_noise <- 0
subgroup_method <- "consistency"
sg_focus <- "maxcons"
target_hr_harm <- 1.0; n_sample <- 500L
method_tag <- if (identical(subgroup_method, "consistency")) "fs" else subgroup_method
# One source of truth: fs_focus_tag() holds the (subgroup_method, sg_focus)
# -> stem tag map for every driver.  It was pasted into each one and drifted;
# the method-blind copies tagged DINA/GRF runs with the consistency name for a
# rule those engines never ran.
focus_tag <- forestsearch::fs_focus_tag(subgroup_method, sg_focus)
if (!identical(focus_tag, sg_focus))
  cat(sprintf("NOTE: sg_focus '%s' is an alias on the %s path; output stem uses focus tag '%s'.
",
              sg_focus, subgroup_method, focus_tag))
rds_stem   <- sprintf("%s_%s_fb_mr_m1_h%02d_knoise%d_n%d",
                      method_tag, focus_tag, round(10 * target_hr_harm),
                      k_random_noise, n_sample)
sim_id_start <- 1L; n_sims <- 5L
sim_id_end <- sim_id_start + n_sims - 1L
rds_path   <- sprintf("%s_res_%d_%d.rds", rds_stem, sim_id_start, sim_id_end)
combine_glob <- sprintf("%s_res_*.rds", rds_stem)
mr_draws <- 500L
dgm_model <- "alt"; analysis_time <- 84; cens_adjust <- log(1.5); n_super <- 100000L
consistency_method <- "resample"
selection_rule <- "neighborhood"; effect_neighborhood <- 0.10
hr_threshold <- 0.90; hr_consistency <- 0.80; pconsistency <- 0.90
fs_splits <- 400L; maxk <- 2L; n_min <- NULL; d0_min <- 10L; d1_min <- 10L
use_lasso <- FALSE; use_dina <- FALSE; use_grf <- FALSE; use_twostage <- TRUE
fs_conf_force <- c("meno == 0", "er <= 0", "pgr <= 0")
fs_conf.cont_jcuts <- list(er = 10)
outcome_name <- "y_sim"; event_name <- "event_sim"; treat_name <- "treat_sim"
id_name <- "id"; harm_col <- "flag_harm"
confounders_base <- c("er", "age", "meno", "pgr", "nodes", "size", "grade")
seed_base <- 8316951L
n_workers <- ceiling(0.90 * max(1L, min(parallel::detectCores(logical = FALSE) - 1L)))
inner_parallel <- if (identical(parallel_mode, "boots"))
  list(plan = "multisession", workers = n_workers) else list(plan = "sequential")

cat("### STRUCTURAL CHECK 1/2 -- setup chunk symbols resolve\n")
cat(sprintf("  rds_stem     = %s\n", rds_stem))
cat(sprintf("  rds_path     = %s\n", rds_path))
cat(sprintf("  combine_glob = %s\n", combine_glob))
cat(sprintf("  n_workers    = %d\n\n", n_workers))

## ---- DGM + one replicate ----------------------------------------------------
k_inter <- calibrate_k_inter(target_hr_harm = target_hr_harm, model = dgm_model,
                             use_ahr = FALSE)
dgm <- setup_gbsg_dgm(model = dgm_model, k_inter = k_inter, n_super = n_super,
                      seed = seed_base)
dgm <- compute_dgm_cde(dgm)

sim_id <- 1L
df <- simulate_from_dgm(dgm, n = n_sample, analysis_time = analysis_time,
                        cens_adjust = cens_adjust, seed = seed_base + sim_id)
df[[id_name]] <- seq_len(nrow(df))
confs <- intersect(confounders_base, names(df))

base_args <- list(
  df.analysis = df, outcome.name = outcome_name, event.name = event_name,
  treat.name = treat_name, id.name = id_name,
  flag_harm.name = harm_col, confounders.name = confs,
  is.RCT = TRUE, seedit = seed_base + sim_id, quiet = TRUE,
  sg_focus = sg_focus, subgroup_method = subgroup_method,
  hr.threshold = hr_threshold, hr.consistency = hr_consistency,
  pconsistency.threshold = pconsistency, n.min = n_min,
  selection_rule = selection_rule,
  effect_neighborhood = effect_neighborhood,
  stop_threshold = NULL,
  parallel_args = inner_parallel,
  mr_inference = TRUE,
  mr_inference_args = list(ci_method = "ij", draws = mr_draws,
                           include_complement = TRUE))
method_args <- list(
  consistency_method = consistency_method,
  use_lasso = use_lasso, use_grf = use_grf, use_twostage = use_twostage,
  use_dina = use_dina, conf_force = fs_conf_force,
  conf.cont_jcuts = fs_conf.cont_jcuts,
  fs.splits = fs_splits, maxk = maxk, d0.min = d0_min, d1.min = d1_min)

cat("### SECTION 2 -- resolved parameters (post-edit call, stop_threshold = NULL pinned)\n")
warns_post <- character(0)
fs.est <- withCallingHandlers(
  do.call(forestsearch, c(base_args, method_args)),
  warning = function(w) { warns_post <<- c(warns_post, conditionMessage(w))
                          invokeRestart("muffleWarning") })

a <- fs.est$args_call_all
show1 <- function(nm) {
  v <- a[[nm]]
  cat(sprintf("  %-22s = %s\n", nm,
              if (is.null(v)) "NULL" else paste(deparse(v), collapse = " ")))
}
for (nm in c("sg_focus", "subgroup_method", "consistency_method",
             "stop_threshold", "selection_rule", "effect_neighborhood",
             "mr_inference")) show1(nm)
cat(sprintf("  %-22s = %s\n", "stop_threshold in list?",
            "stop_threshold" %in% names(a)))
cat(sprintf("  warnings emitted (post-edit route): %d\n", length(warns_post)))
if (length(warns_post)) cat(paste0("    - ", warns_post, collapse = "\n"), "\n")
cat(sprintf("  any mentioning stop_threshold: %s\n\n",
            any(grepl("stop_threshold", warns_post))))

## same fit, pre-edit route: stop_threshold omitted entirely
cat("### SECTION 2b -- same call with stop_threshold OMITTED (pre-edit route)\n")
base_args_pre <- base_args
base_args_pre$stop_threshold <- NULL          # remove the element entirely
base_args_pre$effect_neighborhood <- NULL
base_args_pre <- base_args_pre[!vapply(names(base_args_pre),
  function(n) n %in% c("stop_threshold", "effect_neighborhood"), logical(1))]
base_args_pre$sg_focus <- "eff"               # pre-edit spelling
warns_pre <- character(0)
fs.pre <- withCallingHandlers(
  do.call(forestsearch, c(base_args_pre, method_args)),
  warning = function(w) { warns_pre <<- c(warns_pre, conditionMessage(w))
                          invokeRestart("muffleWarning") })
ap <- fs.pre$args_call_all
for (nm in c("sg_focus", "subgroup_method", "consistency_method",
             "stop_threshold", "selection_rule", "effect_neighborhood",
             "mr_inference")) {
  v <- ap[[nm]]
  cat(sprintf("  %-22s = %s\n", nm,
              if (is.null(v)) "NULL" else paste(deparse(v), collapse = " ")))
}
cat(sprintf("  warnings emitted (pre-edit route): %d\n", length(warns_pre)))
if (length(warns_pre)) cat(paste0("    - ", warns_pre, collapse = "\n"), "\n")
cat(sprintf("  any mentioning stop_threshold: %s\n", any(grepl("stop_threshold", warns_pre))))
cat(sprintf("  sg.harm identical pre vs post: %s\n",
            identical(fs.pre$sg.harm, fs.est$sg.harm)))
cat(sprintf("  mr_inference$debiased identical: %s\n\n",
            isTRUE(all.equal(fs.pre$mr_inference$debiased,
                             fs.est$mr_inference$debiased))))

## ---- §5 field paths ---------------------------------------------------------
cat("### SECTION 5 -- field paths on the real fit\n")
g <- fs.est$mr_inference
cat("  names(fs.est$mr_inference): ", paste(names(g), collapse = ", "), "\n")
cat("  names(...$debiased):        ", paste(names(g$debiased), collapse = ", "), "\n")
for (f in c("est", "lower", "upper", "se_ij")) {
  v <- g$debiased[[f]]
  cat(sprintf("    mr_inference$debiased$%-6s = %s\n", f,
              if (is.null(v)) "<<ABSENT>>" else format(v, digits = 8)))
}
cat("  complement$debiased present: ",
    !is.null(g$complement$debiased), "\n\n")

## ---- §3 FB bootstrap, mr_in_replicates FALSE vs TRUE ------------------------
cat("### SECTION 3 -- mr_in_replicates FALSE vs TRUE, same seed\n")
t0 <- proc.time()[3]
boot_F <- forestsearch_bootstrap_dofuture(
  fs.est, nb_boots = nb_boots, seed = seed_base + sim_id,
  parallel_args = inner_parallel, mr_in_replicates = FALSE)
t_F <- proc.time()[3] - t0

t0 <- proc.time()[3]
boot_T <- forestsearch_bootstrap_dofuture(
  fs.est, nb_boots = nb_boots, seed = seed_base + sim_id,
  parallel_args = inner_parallel, mr_in_replicates = TRUE)
t_T <- proc.time()[3] - t0

cat("  names(boot): ", paste(names(boot_F), collapse = ", "), "\n")
cat("  names(H_estimates): ", paste(names(boot_F$H_estimates), collapse = ", "), "\n")
for (f in c("H2", "H2_lower", "H2_upper", "sdH2")) {
  v <- boot_F$H_estimates[[f]]
  cat(sprintf("    H_estimates$%-9s = %s\n", f,
              if (is.null(v)) "<<ABSENT>>" else format(v, digits = 8)))
}
for (f in c("H2", "H2_lower", "H2_upper", "sdH2")) {
  v <- boot_F$Hc_estimates[[f]]
  cat(sprintf("    Hc_estimates$%-8s = %s\n", f,
              if (is.null(v)) "<<ABSENT>>" else format(v, digits = 8)))
}
cat(sprintf("\n  H_estimates  identical: %s\n",
            identical(boot_F$H_estimates,  boot_T$H_estimates)))
cat(sprintf("  Hc_estimates identical: %s\n",
            identical(boot_F$Hc_estimates, boot_T$Hc_estimates)))
per_rep_names <- intersect(names(boot_F),
  c("df.est", "dfRes", "boot_results", "result_boot", "df_boot", "out.boots"))
for (nm in per_rep_names)
  cat(sprintf("  per-replicate table '%s' identical: %s\n", nm,
              identical(boot_F[[nm]], boot_T[[nm]])))
cat(sprintf("  mr_replicates: FALSE-> %s ; TRUE-> %s\n",
            if (is.null(boot_F$mr_replicates)) "NULL" else
              sprintf("list of %d", length(boot_F$mr_replicates)),
            if (is.null(boot_T$mr_replicates)) "NULL" else
              sprintf("list of %d", length(boot_T$mr_replicates))))
setdiff_names <- setdiff(names(boot_T), names(boot_F))
cat("  names present only with TRUE: ",
    if (length(setdiff_names)) paste(setdiff_names, collapse = ", ") else "(none)", "\n")
# full-object comparison ignoring the MR payload and timings
drop <- c("mr_replicates")
bF <- boot_F[setdiff(names(boot_F), drop)]
bT <- boot_T[setdiff(names(boot_T), drop)]
cmp <- all.equal(bF, bT)
cat("  all.equal(boot_F, boot_T) ex mr_replicates: ",
    if (isTRUE(cmp)) "TRUE" else paste(cmp, collapse = " | "), "\n")
cat(sprintf("\n  wall-clock: FALSE %.1f s ; TRUE %.1f s ; ratio %.2fx\n",
            t_F, t_T, t_T / t_F))

plan("sequential")
cat("\n### DONE\n")
