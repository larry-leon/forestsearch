# =============================================================================
# dina_grf_routing_diagnostic.R
# -----------------------------------------------------------------------------
# Purpose: isolate WHY the dina-* and grf-* arms detect nothing in the methods
# template, by running -- on ONE identical simulated trial -- the legacy GRF
# path (.run_grf_analysis_gen -> grf.subg.harm.glm, the prior vignette's 98%
# path) against the forestsearch-routed arms (forestsearch(subgroup_method=...))
# and printing details = TRUE traces for GRF and DINA.
#
# Run:  Rscript dina_grf_routing_diagnostic.R     (or source() in RStudio)
# Requires: forestsearch installed (devtools::install(), NOT load_all).
# =============================================================================

suppressPackageStartupMessages({
  library(forestsearch)
  library(speff2trial)
  library(data.table)
})

cal_target_or <- 2.0
sim_n         <- 700L
SIM_ID        <- 42L          # fixed replicate index
seed_base     <- 8316951L     # run_simulation_analysis() default; the trial it
                              # analyzes for sim_id is drawn at seed_base + sim_id

# ---- 1. ACTG-175 data prep (verbatim from the template) ---------------------
actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
actg_df$id    <- seq_len(nrow(actg_df))
actg_df$treat <- ifelse(actg_df$arms == 1L, 1L, 0L)
actg_df$y_binary <- as.integer(actg_df$cd420 <= actg_df$cd40)
actg_df <- actg_df[!is.na(actg_df$cd420), ]
actg_df$z1 <- as.factor(ifelse(actg_df$cd40 > 390, 1L, 0L))
actg_df$z2 <- as.factor(ifelse(actg_df$wtkg > 80,  1L, 0L))
actg_df$z7  <- as.factor(actg_df$hemo);    actg_df$z8  <- as.factor(actg_df$homo)
actg_df$z9  <- as.factor(actg_df$drugs);   actg_df$z10 <- as.factor(actg_df$race)
actg_df$z11 <- as.factor(actg_df$gender);  actg_df$z12 <- as.factor(actg_df$symptom)
actg_df$ar_naive  <- as.integer(ifelse(actg_df$str2 == 0, 1L, 0L))
actg_df$prior_6mo <- as.integer(ifelse(actg_df$preanti <= 6 * 30.4375, 1L, 0L))

cont_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
bin_vars  <- c("hemo", "homo", "drugs", "race", "gender", "symptom",
               "ar_naive", "prior_6mo")
# Binary confounders kept NUMERIC 0/1: DINA requires numeric covariates and
# the GRF policy-tree emits numeric cut-strings (e.g. "symptom <= 0") that are
# invalid for factors.  z7-z12 (created above) stay factors for the DGM.
for (v in bin_vars) actg_df[[v]] <- as.integer(actg_df[[v]])
confounders <- c(cont_vars, bin_vars)

dgm_factor_vars   <- c("z1", "z2", paste0("z", 7:12))
dgm_subgroup_vars <- c("z1", "z2")
dgm_subgroup_cuts <- list(z1 = 1L, z2 = 1L)

# ---- 2. Calibrated DGM (OR(Q) = 2.0 harm signal) ----------------------------
dgm_calibrated <- calibrate_glm_interaction(
  data = actg_df, factor_vars = dgm_factor_vars, continuous_vars = cont_vars,
  outcome_var = "y_binary", treatment_var = "treat",
  target_effect = cal_target_or, outcome_type = "binary", effect_measure = "OR",
  subgroup_vars = dgm_subgroup_vars, subgroup_cuts = dgm_subgroup_cuts,
  k_inter_range = c(0, 3), grid_step = 0.05, n_super = 5000L,
  seed = 8316951L, verbose = FALSE
)
cat(sprintf("DGM: theta-dagger(Q) = %.3f (target OR = %.1f)\n\n",
            dgm_calibrated$hazard_ratios$harm_subgroup, cal_target_or))

# ---- 3. Shared base params (verbatim from the template) ---------------------
fs_base <- list(
  outcome.name = "y_sim", event.name = "y_sim", treat.name = "treat_sim",
  id.name = "id", outcome_type = "binary", effect_measure = "OR",
  hr.threshold = 1.25, hr.consistency = 1.0, pconsistency.threshold = 0.90,
  stop_threshold = NULL,
  conf.cont_jcuts = list(cd40 = 10, wtkg = 10),
  conf_force = c("ar_naive == 1", "prior_6mo == 1"),
  maxk = 2L, fs.splits = 1000L,
  n.min = NULL, d0.min = 10L, d1.min = 10L,
  use_grf = TRUE, use_dina = TRUE, use_lasso = FALSE, vi.grf.min = -0.2,
  is.RCT = TRUE, seedit = 8316951L,
  parallel_args = list(plan = "sequential", workers = 1, show_message = FALSE)
)

# GRF SELECTION floors set to the permissive baseline for this test.
grf_select_dmin     <- 0.0
# Fairness: every engine searches the same two-factor space (matches the
# signature-recovery vignette and make_methods() in the methods vignette).
dina_max_depth      <- 2L
dina_grid_probs     <- seq(0.1, 0.9, by = 0.1)
effect_neighborhood <- 0.10

# ---- 4. forestsearch-routed methods spec (matches make_methods()) -----------
METHODS_SPEC <- list(
  "fs-effMaxSG"   = list(subgroup_method = "consistency", sg_focus = "effMaxSG",
                         selection_rule = "neighborhood",
                         dmin.grf = grf_select_dmin),
  "dina-effMaxSG" = list(subgroup_method = "dina", sg_focus = "effMaxSG",
                         selection_rule = "neighborhood",
                         effect_neighborhood = effect_neighborhood,
                         dina_args = list(max_depth = dina_max_depth,
                                          grid_probs = dina_grid_probs)),
  "grf-effMaxSG"  = list(subgroup_method = "grf", grf_selection = "frontier",
                         sg_focus = "effMaxSG", grf_depth = 2L,
                         dmin.grf = grf_select_dmin,
                         effect_neighborhood = effect_neighborhood,
                         use_grf = FALSE, use_lasso = FALSE),
  "grf-tree"      = list(subgroup_method = "grf", grf_selection = "tree",
                         grf_depth = 2L, dmin.grf = grf_select_dmin,
                         use_grf = FALSE, use_lasso = FALSE)
)

`%||%` <- function(a, b) if (is.null(a)) b else a
tidy_sg_def <- function(s, digits = 4L) {
  if (length(s) != 1L || is.na(s) || !nzchar(s)) return(s)
  s <- gsub("[[:space:]]+", " ", s)
  m <- gregexpr("-?[0-9]+\\.[0-9]+", s)[[1L]]
  if (m[1L] == -1L) return(s)
  L <- attr(m, "match.length")
  for (i in rev(seq_along(m))) {
    num <- as.numeric(substr(s, m[i], m[i] + L[i] - 1L))
    rep <- sub("\\.$", "", formatC(signif(num, digits), format = "fg", flag = "#"))
    s <- paste0(substr(s, 1L, m[i] - 1L), rep, substr(s, m[i] + L[i], nchar(s)))
  }
  s
}
show_det <- function(res, tag) {
  if (is.null(res) || !nrow(res)) { cat(sprintf("  [%s] no rows\n", tag)); return(invisible()) }
  for (i in seq_len(nrow(res))) {
    sgd <- tidy_sg_def(res$sg.def[i] %||% "")
    cat(sprintf("  %-14s found=%s  sg=%s  |H|=%s\n",
                res$analysis[i],
                res$any.H[i] %||% NA,
                if (!nzchar(sgd)) "<none>" else sgd,
                res$size.H[i] %||% NA))
  }
}

# ---- 5a. LEGACY path: FS (consistency) + standalone GRF ---------------------
cat("================ LEGACY path (methods = NULL) ================\n")
cat("FS via consistency; GRF via .run_grf_analysis_gen -> grf.subg.harm.glm\n")
res_legacy <- tryCatch(
  run_simulation_analysis(
    sim_id = SIM_ID, dgm = dgm_calibrated, n_sample = sim_n,
    confounders_base = confounders,
    methods = NULL, run_fs = TRUE, run_fs_grf = FALSE, run_grf = TRUE,
    fs_params = fs_base, verbose = FALSE
  ),
  error = function(e) { cat("  LEGACY run errored:", conditionMessage(e), "\n"); NULL }
)
show_det(res_legacy, "legacy")

# ---- 5b. METHODS path: forestsearch-routed fs / dina / grf ------------------
cat("\n============ METHODS path (forestsearch-routed) ============\n")
res_methods <- tryCatch(
  run_simulation_analysis(
    sim_id = SIM_ID, dgm = dgm_calibrated, n_sample = sim_n,
    confounders_base = confounders,
    methods = METHODS_SPEC, fs_params = fs_base, verbose = FALSE
  ),
  error = function(e) { cat("  METHODS run errored:", conditionMessage(e), "\n"); NULL }
)
show_det(res_methods, "methods")

# ---- 6. Verbose traces for the two failing engines on one trial -------------
# Reproduce the SAME trial the methods/legacy paths analyzed: run_simulation_analysis()
# simulates with seed = seed_base + sim_id, so the trace must use that seed too
# (using bare SIM_ID would draw a DIFFERENT trial and the subgroups would not match).
df_sim     <- simulate_from_glm_dgm(dgm_calibrated, n = sim_n, seed = seed_base + SIM_ID)
fs_formals <- names(formals(forestsearch))
.call_fs <- function(overrides) {
  args <- c(list(df.analysis = df_sim, confounders.name = confounders,
                 details = TRUE, plot.sg = FALSE, quiet = FALSE),
            utils::modifyList(fs_base, overrides)[
              intersect(names(utils::modifyList(fs_base, overrides)), fs_formals)])
  args$outcome.name <- "y_sim"; args$event.name <- "y_sim"
  do.call(forestsearch, args)
}

cat("\n================ TRACE: forestsearch GRF (frontier) ================\n")
g <- tryCatch(.call_fs(list(subgroup_method = "grf", grf_selection = "frontier",
                            sg_focus = "effMaxSG", grf_depth = 2L,
                            dmin.grf = grf_select_dmin,
                            effect_neighborhood = effect_neighborhood,
                            use_grf = FALSE, use_lasso = FALSE)),
              error = function(e) { cat("  errored:", conditionMessage(e), "\n"); NULL })
cat("  -> sg.harm:", if (is.null(g$sg.harm)) "<none>" else paste(g$sg.harm, collapse = " & "), "\n")

cat("\n================ TRACE: forestsearch DINA ================\n")
d <- tryCatch(.call_fs(list(subgroup_method = "dina", sg_focus = "effMaxSG",
                            selection_rule = "neighborhood",
                            effect_neighborhood = effect_neighborhood,
                            dina_args = list(max_depth = dina_max_depth,
                                             grid_probs = dina_grid_probs))),
              error = function(e) { cat("  errored:", conditionMessage(e), "\n"); NULL })
cat("  -> sg.harm:", if (is.null(d$sg.harm)) "<none>" else paste(d$sg.harm, collapse = " & "), "\n")

cat("\n================ SUMMARY ================\n")
cat("Compare legacy GRF detection vs methods grf-* detection on the SAME trial.\n")
cat("If legacy GRF finds a subgroup but grf-effMaxSG/grf-tree do not, the gap\n")
cat("is in forestsearch's GRF selection layer, not the data. The DINA trace\n")
cat("above shows why DINA selects nothing (it ignores dmin.grf).\n")
