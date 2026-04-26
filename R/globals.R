#' Global Variables Declaration for ForestSearch Package
#'
#' Consolidated declarations to suppress R CMD check NOTEs for
#' non-standard evaluation (NSE) in data.table, ggplot2, and dplyr operations.
#'
#' These variables are used in contexts where column names are referenced
#' directly without quotes. R CMD check cannot detect these are valid
#' column references, so we declare them here to suppress false positive warnings.
#'
#' @keywords internal
#' @noRd

utils::globalVariables(c(


# ============================================================================
# Data.table special symbols
# ============================================================================
  ".",
  ".N",
  ".SD",
  ".I",
  ".GRP",
  ".BY",
  "..col_names",

# ============================================================================
# Future package
# ============================================================================
  ".Last.future.plan",

# ============================================================================
# Bootstrap & bias correction
# ============================================================================
  "boot",
  "boot_id",
  "iteration",
  "converged",
  "error_message",
  "time_elapsed",
  "H",
  "H_star",
  "H_bootstrap",
  "H_biasadj_1",
  "H_biasadj_2",
  "Hc_biasadj_1",
  "Hc_biasadj_2",
  "grf_cuts_b",
  "grf_cut",

# ============================================================================
# Subgroup identification
# ============================================================================
  "Subgroup",
  "match_string",
  "treat.recommend",
  "flag.harm",
  "flag_harm",
  "M.1",
  "M.2",
  "M.3",
  "M.4",
  "M.5",
  "M.6",
  "M.7",
  "grp",
  "g_sg",
  "m_sg",

# ============================================================================
# Metrics: HR, consistency, sample size
# ============================================================================
  "HR",
  "hr",
  "hr_sg",
  "hr_individual",
  "Pcons",
  "Pcons_bin",
  "N",
  "N_sg",
  "N_sg_bin",
  "K",
  "K_sg",
  "E_sg",
  "Percent",
  "Percent_of_successful",
  "Rank",

# ============================================================================
# Summary table columns
# ============================================================================
  "Category",
  "Metric",
  "Value",
  "Count",
  "Base_Factor",
  "Factor_Definition",
  "Factor",
  "Consistency Range",
  "Size Range",
  "N_positions",
  "N_total",
  "Positions",
  "Position",

# ============================================================================
# Cross-validation (includes foreach loop variable from Kfold/tenfold)
# ============================================================================
  "cv_index",
  "cvindex",
  "ksim",

# ============================================================================
# ForestSearch parameters (used in NSE contexts)
# ============================================================================
  "est.scale",
  "insplit1",
  "confounders.name",
  "outcome.name",
  "event.name",
  "treat.name",

# ============================================================================
# Duplicate detection
# ============================================================================
  "dup_key",
  "dup_group",
  "dup_count",
  "dup_rank",

# ============================================================================
# Survival analysis & Cox models
# ============================================================================
  "treat",
  "event",
  "Y",
  "E",
  "Treat",
  "stratum",
  "z",
  "loghr_po",
  "theta_1",
  "theta_0",
  "treat_harm",

# ============================================================================
# AFT DGM / Simulation variables
# ============================================================================
  "y",
  "y_sim",
  "event_sim",
  "t_true",
  "c_time",
  "lin_pred_1",
  "lin_pred_0",
  "lin_pred_obs",
  "lin_pred_cens_1",
  "lin_pred_cens_0",
  "z_age",
  "z_size",
  "z_nodes",
  "z_pgr",
  "z_er",
  "z_meno",
  "z_grade_1",
  "z_grade_2",
  "z_hormon",

# ============================================================================
# GBSG dataset column names (used in examples/tests)
# ============================================================================
  "id",
  "pid",
  "status",
  "rfstime",
  "age",
  "size",
  "nodes",
  "pgr",
  "er",
  "meno",
  "grade",
  "hormon",
  "desc",

# ============================================================================
# Model comparison
# ============================================================================
  "By_AIC",
  "By_BIC",
  "By_LogLik",
  "AIC_value",
  "BIC_value",
  "LogLik_value",

# ============================================================================
# GBSG DGM internal variables (from sim_aft_gbsg_refactored.R)
# ============================================================================
  "z1", "z2", "z3", "z4", "z5",
  "zh",
  "v1", "v2", "v3", "v4", "v5", "v6", "v7",
  "grade3",
  "lin.conf.true", "lin1.conf", "lin0.conf",
  "linC1.conf", "linC0.conf",
  "hlin.conf.1", "hlin.conf.0", "hlin.ratio",
  "h1.potential", "h0.potential",
  "Ts", "es",
  "t.sim",

# ============================================================================
# Operating characteristics / simulation analysis
# (from oc_analyses.R, simulation_tables.R)
# ============================================================================
  "any.H",
  "ppv", "npv", "sens", "spec",
  "size.H", "size.Hc",
  "hr.H.true", "hr.H.hat",
  "hr.Hc.true", "hr.Hc.hat",
  "hr.itt", "hr.adj.itt",
  "p.cens", "taumax",
  "sim_id",
  "analysis",
  "trimmed",
  ".trimmed",
  "aa",
  "sg_hat",
  "ahr.H.true", "ahr.Hc.true",
  "ahr.H.hat", "ahr.Hc.hat",
  "cde.H.true", "cde.Hc.true",
  "cde.H.hat", "cde.Hc.hat",
  "CDE", "CDE_H_true", "CDE_Hc_true",
  "CDE_harm", "CDE_no_harm",
  "Group",
  "Scenario",
  "Estimator",

# ============================================================================
# MRCT simulation variables (from mrct_simulation.R, plot_sg_distribution.R)
# ============================================================================
  "sim",
  "hr_test",
  "any_found", "sg_found", "hr_sg_null",
  "hr_sg_train", "POhr_sg_train", "n_sg_train",
  "regAflag", "sg_le85", "regAflag2", "regAflag3",
  "found",
  "sg_biomarker", "sg_age", "sg_male", "sg_ecog",
  "sg_histology", "sg_CTregimen", "sg_region",
  "sg_surgery", "sg_prior_treat",
  "est",
  "region_var", "z_regA",
  "sg_label", "count", "bar_label",

# ============================================================================
# gg_forest() -- ggplot2 aes() column names
# ============================================================================
  "subgroup",   # y aesthetic in label and CI panels
  "row_col",    # colour/fill aesthetic mapped from data frame column
  "value",      # label aesthetic in annotation panels
  "col",        # colour aesthetic in label panel
  "lo",         # xend of CI segment
  "hi"          # x of CI upper bound (also used in geom_point)
))


# ============================================================================
# Additional package-level imports
# ============================================================================
# These imports supplement those in forestsearch-package.R and ensure
# all commonly used functions are available without explicit namespacing.

#' @import data.table
#' @import survival
#' @importFrom stats AIC BIC model.frame vcov sd median quantile
#' @importFrom stats cor coef predict residuals fitted formula terms update
#' @importFrom stats pchisq qnorm qt pnorm pt
#' @importFrom stats rexp runif rnorm rbinom rgamma
#' @importFrom stats na.omit complete.cases aggregate
#' @importFrom graphics hist barplot mtext axis title rug text rect plot.window
#' @importFrom utils data object.size head tail str
#' @importFrom utils capture.output write.table read.table
#' @importFrom utils packageVersion sessionInfo
#' @importFrom splines ns
NULL

# ============================================================================
# Notes on base R functions (no import needed)
# ============================================================================
# The following functions are from base R and don't need importing:
# - Math: mean, sum, min, max, exp, log, sqrt, abs, round, floor, ceiling
# - Structure: length, nrow, ncol, dim, names, colnames, rownames
# - Creation: c, list, data.frame, matrix, as.matrix, as.data.frame
# - Manipulation: subset, which, order, sort, unique, duplicated
# - Text: paste, paste0, sprintf, cat, print
# - Control: if, else, for, while, return, stop, warning, message
# - Constants: TRUE, FALSE, NULL, NA, NaN, Inf
