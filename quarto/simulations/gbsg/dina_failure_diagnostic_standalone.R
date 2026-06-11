# ============================================================================
# Self-contained DINA single-trial diagnostic
# ----------------------------------------------------------------------------
# Reproduces ONE replicate of dina-tier2_gate_standalone_simulation.qmd and
# runs the DINA forestsearch() call WITHOUT the harness's error-swallowing
# tryCatch (and with details = TRUE / quiet = FALSE), so the real error
# message -- or the DINA selection trace -- prints to the console instead of
# collapsing to an all-NA record.
#
# HOW TO RUN
#   - RStudio: open this file, Source it (Cmd/Ctrl+Shift+S), read the console;
#     or  source("dina_failure_diagnostic_standalone.R").
#   - Terminal: Rscript dina_failure_diagnostic_standalone.R
#
# REQUIRES: devtools::install() of the current forestsearch build first
#   (NOT load_all(); the gate-wrap fix must be in the installed package).
#
# TO LOCALISE: run once as-is (dina_select_statistic = "effect"), then flip
#   it to "dina" and run again.  If "dina" works, the effect re-rank is the
#   culprit; if it still errors, the effect path is innocent.
# ============================================================================

suppressMessages({
  library(forestsearch)
  library(survival)
  library(future)
})

# ── Knobs (mirrors the qmd setup chunk) ─────────────────────────────────────
# DGM
dgm_model      <- "alt"
target_hr_harm <- 2.0
n_sample       <- 700L
analysis_time  <- 84
cens_adjust    <- log(1.5)
n_super        <- 5000L
seed_base      <- 8316951L

# Selection engine + selector under test
subgroup_method       <- "dina"
dina_args             <- list()
dina_select_statistic <- "effect"      # <- flip to "dina" to localise

# Search / gate configuration
sg_focus       <- "effMinSG"
selection_rule <- "neighborhood"
hr_threshold   <- 1.0
hr_consistency <- 1.0
pconsistency   <- 0.90
n_min          <- 60L
gate_draws     <- 2000L

# Column names emitted by simulate_from_dgm()
outcome_name <- "y_sim"
event_name   <- "event_sim"
treat_name   <- "treat_sim"
id_name      <- "id"
harm_col     <- "flag_harm"
confounders_base <- c("er", "age", "meno", "pgr", "nodes", "size", "grade")

# ── Build the DGM (same calls as the build-dgm chunk) ───────────────────────
k_inter <- calibrate_k_inter(target_hr_harm = target_hr_harm,
                             model = dgm_model, use_ahr = FALSE)
dgm <- setup_gbsg_dgm(model = dgm_model, k_inter = k_inter,
                      n_super = n_super, seed = seed_base)
dgm <- compute_dgm_cde(dgm)

# ── Simulate ONE trial ──────────────────────────────────────────────────────
df    <- simulate_from_dgm(dgm, n = n_sample, analysis_time = analysis_time,
                           cens_adjust = cens_adjust, seed = seed_base + 1L)
confs <- intersect(confounders_base, names(df))
cat(sprintf("n = %d | events = %d | arms = %d | confs = %s\n",
            nrow(df), sum(df[[event_name]]),
            length(unique(df[[treat_name]])), paste(confs, collapse = ", ")))
cat(sprintf("subgroup_method = %s | dina_select_statistic = %s | sg_focus = %s\n\n",
            subgroup_method, dina_select_statistic, sg_focus))

# ── ONE DINA forestsearch() call, error EXPOSED ─────────────────────────────
fs <- tryCatch(
  forestsearch(
    df.analysis            = df,
    outcome.name           = outcome_name, event.name = event_name,
    treat.name             = treat_name,   id.name    = id_name,
    flag_harm.name         = harm_col,     confounders.name = confs,
    is.RCT                 = TRUE,         seedit = seed_base + 1L,
    sg_focus               = sg_focus,    subgroup_method = subgroup_method,
    hr.threshold           = hr_threshold, hr.consistency = hr_consistency,
    pconsistency.threshold = pconsistency, n.min = n_min,
    selection_rule         = selection_rule,
    parallel_args          = list(plan = "sequential"),
    debias_gate            = TRUE,
    debias_gate_args       = list(ci_method = "ij", draws = gate_draws,
                                  include_complement = TRUE),
    dina_args = modifyList(dina_args, list(select_statistic = dina_select_statistic)),
    details = TRUE, quiet = FALSE),                # <- expose DINA's decisions
  error = function(e) {
    message("\n>>> forestsearch() ERROR: ", conditionMessage(e))
    cc <- conditionCall(e)
    if (!is.null(cc)) message(">>> call: ", deparse(cc)[1])
    NULL
  })

# ── Report ──────────────────────────────────────────────────────────────────
if (!is.null(fs)) {
  cat("\nsg.harm :", if (is.null(fs$sg.harm)) "NULL (no subgroup)"
                     else paste(fs$sg.harm, collapse = " & "), "\n")
  cat("gate    :", if (is.null(fs$debias_gate)) "NULL" else "present", "\n")
  if (!is.null(fs$debias_gate)) {
    g <- fs$debias_gate
    cat(sprintf("  naive HR = %.3f | debiased HR = %.3f [%.3f, %.3f] | n_family = %d\n",
                g$naive$est, g$debiased$est, g$debiased$lower, g$debiased$upper,
                g$n_family))
  }
}
