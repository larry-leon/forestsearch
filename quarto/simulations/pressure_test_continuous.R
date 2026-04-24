# =============================================================================
# Pressure test: continuous vignette with GRF re-enabled
# =============================================================================
# Exercises:
#   - setup chunk variables (new grf_* block, grf_params list)
#   - the updated chunks (oc-alt with dual build_estimation_table calls,
#     classification, oc-sgdist-build, effect-distribution violins)
# using realistic fake run_simulation_analysis() output that includes
# BOTH FS and GRF rows.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

# source("plot_sg_distribution.R")
# source("plot_effect_distribution.R")

PASS <- 0L; FAIL <- 0L
ok <- function(label, expr) {
  res <- tryCatch(expr, error = function(e) e)
  if (inherits(res, "error")) {
    FAIL <<- FAIL + 1L
    cat(sprintf("  [FAIL] %s\n         %s\n", label, conditionMessage(res)))
    return(invisible(NULL))
  }
  if (isTRUE(res)) {
    PASS <<- PASS + 1L
    cat(sprintf("  [ OK ] %s\n", label))
  } else {
    FAIL <<- FAIL + 1L
    cat(sprintf("  [FAIL] %s (got: %s)\n", label,
                deparse(res, width.cutoff = 60)[1]))
  }
}

set.seed(8316951L)

# =============================================================================
# 1. Emulate the continuous vignette setup chunk
# =============================================================================
cat("\n=== Continuous vignette setup chunk ===\n\n")

# Simulation control
nsims_alt   <- 25L
nsims_null  <- 50L
sim_n_sample <- 1000L

# Analysis flags (the edit)
run_fs     <- TRUE
run_fs_grf <- FALSE
run_grf    <- TRUE        # <-- flipped

# GRF parameters (the new block)
grf_dmin         <- 10
grf_maxdepth     <- 2L
grf_sg_criterion <- "mDiff"
grf_tune         <- FALSE

grf_params <- list(
  dmin.grf     = grf_dmin,
  maxdepth     = grf_maxdepth,
  sg.criterion = grf_sg_criterion,
  tune_grf     = grf_tune
)

# Reporting + visualization
sg_notation            <- "benefit"
sgdist_min_pct_alt     <- 5
sgdist_min_pct_null    <- 5
sgdist_top_k           <- 15L
sgdist_show_other_alt  <- FALSE
sgdist_show_other_null <- FALSE

ok("1.1 run_grf flipped to TRUE",               isTRUE(run_grf))
ok("1.2 grf_params has 4 elements",             length(grf_params) == 4L)
ok("1.3 grf_params has required names",
   identical(sort(names(grf_params)),
             sort(c("dmin.grf", "maxdepth", "sg.criterion", "tune_grf"))))
ok("1.4 grf_sg_criterion = 'mDiff' (identity)",
   identical(grf_sg_criterion, "mDiff"))

# =============================================================================
# 2. Fake run_simulation_analysis() output with BOTH FS and GRF rows
# =============================================================================
cat("\n=== Fake results with FS + GRF rows ===\n\n")

# Build one row-per-(sim, method) block mirroring the real shape.
# For MD/continuous: scale values are mean differences (identity scale).
# Detection rates and typical magnitudes emulated from prior runs:
#   FS:  det ~ 60%, H-hat switched MD ~ +40
#   GRF: det ~ 80%, H-hat switched MD ~ +35 (lower power, similar mean)
build_continuous_results <- function(nsim, method, p_detect,
                                     h_mean, hc_mean = 0,
                                     itt_mean = -5) {
  det <- rbinom(nsim, 1L, p_detect)
  labels_pool <- c("{z1 = 1}", "{z2 = 1}", "{z1 = 1 & z2 = 1}",
                   paste0("{z", 3:10, "_high = 1}"))
  data.frame(
    analysis  = method,
    any.H     = det,
    hr.itt    = rnorm(nsim, mean = itt_mean, sd = 2),
    hr.H.hat  = ifelse(det == 1L, rnorm(nsim, h_mean, 5), NA_real_),
    hr.Hc.hat = ifelse(det == 1L, rnorm(nsim, hc_mean, 3), NA_real_),
    sg.def    = ifelse(det == 1L,
      sample(labels_pool, nsim, replace = TRUE,
             prob = c(0.25, 0.20, 0.15, rep(0.05, 8))), ""),
    stringsAsFactors = FALSE
  )
}

results_alt <- rbind(
  build_continuous_results(nsims_alt, "FS",  p_detect = 0.60,
                           h_mean = 40, hc_mean = 0,  itt_mean = -5),
  build_continuous_results(nsims_alt, "GRF", p_detect = 0.80,
                           h_mean = 35, hc_mean = 0,  itt_mean = -5)
)
results_null <- rbind(
  build_continuous_results(nsims_null, "FS",  p_detect = 0.04,
                           h_mean = 10, hc_mean = 0,  itt_mean = 0),
  build_continuous_results(nsims_null, "GRF", p_detect = 0.10,
                           h_mean = 10, hc_mean = 0,  itt_mean = 0)
)

ok("2.1 results_alt has FS and GRF rows",
   all(c("FS", "GRF") %in% unique(results_alt$analysis)))
ok("2.2 results_null has FS and GRF rows",
   all(c("FS", "GRF") %in% unique(results_null$analysis)))
ok("2.3 alt schema has required columns",
   all(c("analysis", "any.H", "hr.itt", "hr.H.hat", "hr.Hc.hat",
         "sg.def") %in% names(results_alt)))

# =============================================================================
# 3. Chunk oc-sgdist-build — must work on FS-only (the vignette filters)
# =============================================================================
cat("\n=== Chunk: oc-sgdist-build ===\n\n")

plot_sgdist_alt <- plot_sg_distribution(
  results              = results_alt[results_alt$analysis == "FS", ],
  any_col              = "any.H",
  label_col            = "sg.def",
  min_pct              = sgdist_min_pct_alt,
  top_k                = sgdist_top_k,
  show_other           = sgdist_show_other_alt,
  max_halvings         = 2L,
  bar_label_inside     = TRUE,
  placeholder_on_empty = TRUE,
  title                = "HTE: Subgroups Identified (H-hat)",
  wrap_width           = 30
)

plot_sgdist_null <- plot_sg_distribution(
  results              = results_null[results_null$analysis == "FS", ],
  any_col              = "any.H",
  label_col            = "sg.def",
  min_pct              = sgdist_min_pct_null,
  top_k                = sgdist_top_k,
  show_other           = sgdist_show_other_null,
  max_halvings         = 2L,
  bar_label_inside     = TRUE,
  placeholder_on_empty = TRUE,
  title                = "Null: Subgroups Identified (H-hat)",
  wrap_width           = 30
)

h_alt  <- sgdist_fig_height(attr(plot_sgdist_alt,  "n_bars"))
h_null <- sgdist_fig_height(attr(plot_sgdist_null, "n_bars"))

ok("3.1 sgdist H1 returns ggplot",  inherits(plot_sgdist_alt, "ggplot"))
ok("3.2 sgdist H0 returns ggplot",  inherits(plot_sgdist_null, "ggplot"))
ok("3.3 h_alt in [4, 10]",          h_alt  >= 4 && h_alt  <= 10)
ok("3.4 h_null in [4, 10]",         h_null >= 4 && h_null <= 10)

# =============================================================================
# 4. Chunk oc-alt-mddist and oc-null-mddist (FS only, MD scale, benefit)
# =============================================================================
cat("\n=== Chunks: oc-alt-mddist, oc-null-mddist ===\n\n")

# Fake DGM for reference effect
dgm_calibrated <- list(hazard_ratios = list(
  harm_subgroup = 40, no_harm_subgroup = 0, overall = -5))
dgm_null <- list(hazard_ratios = list(overall = 0))

md_Q_true <- -dgm_calibrated$hazard_ratios$harm_subgroup

p_alt_md <- plot_effect_distribution(
  results           = results_alt,
  analysis_method   = "FS",
  effect_measure    = "MD",
  subgroup_notation = sg_notation,
  reference_effect  = md_Q_true,
  title             = "ACTG175 Continuous (HTE): MD Estimates",
  subtitle          = sprintf("n=%d, %d sims | MD(Q) = %.1f",
                              sim_n_sample, nsims_alt, md_Q_true)
)

p_null_md <- plot_effect_distribution(
  results           = results_null,
  analysis_method   = "FS",
  effect_measure    = "MD",
  subgroup_notation = sg_notation,
  title             = "ACTG175 Continuous (Null): MD Estimates",
  subtitle          = sprintf("n=%d, %d sims", sim_n_sample, nsims_null)
)

ok("4.1 H1 mddist returns ggplot",  inherits(p_alt_md,  "ggplot"))
ok("4.2 H0 mddist returns ggplot",  inherits(p_null_md, "ggplot"))
# Benefit-flip check: H-hat mean on benefit scale should be near -40
pd_alt <- attr(p_alt_md, "panel_data")
H_mean_benefit <- mean(
  pd_alt$est[grepl("^Identified H-hat", as.character(pd_alt$group))],
  na.rm = TRUE)
ok("4.3 FS H-hat benefit-scale mean near -40",
   abs(H_mean_benefit - (-40)) < 5)

# =============================================================================
# 5. Verify plot_effect_distribution works on GRF rows too
#    (the vignette doesn't call this yet, but we pressure-test it in case
#     a future edit adds a GRF violin)
# =============================================================================
cat("\n=== GRF slicing robustness ===\n\n")

p_grf <- plot_effect_distribution(
  results           = results_alt,
  analysis_method   = "GRF",
  effect_measure    = "MD",
  subgroup_notation = sg_notation,
  reference_effect  = md_Q_true
)
ok("5.1 GRF slice returns ggplot", inherits(p_grf, "ggplot"))

# The H-hat mean for GRF on benefit scale should be near -35
pd_grf <- attr(p_grf, "panel_data")
H_mean_grf <- mean(
  pd_grf$est[grepl("^Identified H-hat", as.character(pd_grf$group))],
  na.rm = TRUE)
ok("5.2 GRF H-hat benefit mean near -35",
   abs(H_mean_grf - (-35)) < 5)

# =============================================================================
# 6. Grob materialization — catches ggplot-layer errors
# =============================================================================
cat("\n=== Grob materialization ===\n\n")

all_plots <- list(
  sgdist_alt  = plot_sgdist_alt,
  sgdist_null = plot_sgdist_null,
  md_alt      = p_alt_md,
  md_null     = p_null_md,
  md_grf      = p_grf
)
for (nm in names(all_plots)) {
  ok(sprintf("6.x grob: %s", nm),
     inherits(ggplot_gtable(ggplot_build(all_plots[[nm]])), "gtable"))
}

# =============================================================================
# 7. Simulate classification/estimation-table consumption of FS+GRF
# =============================================================================
cat("\n=== Classification / estimation table inputs ===\n\n")

# build_classification_table uses sort(unique(results$analysis))
analyses <- sort(unique(c(unique(results_null$analysis),
                          unique(results_alt$analysis))))
ok("7.1 classification analyses contains FS+GRF",
   identical(analyses, c("FS", "GRF")))

# build_estimation_table(analysis_method = "FS") should find FS rows
fs_rows <- results_alt[results_alt$analysis == "FS", ]
grf_rows <- results_alt[results_alt$analysis == "GRF", ]
ok("7.2 FS slice has >=1 row",  nrow(fs_rows)  >= 1L)
ok("7.3 GRF slice has >=1 row", nrow(grf_rows) >= 1L)

# Sanity: detection rates differ between FS and GRF (they should)
ok("7.4 GRF detection rate > FS under H1",
   mean(grf_rows$any.H) > mean(fs_rows$any.H))

# =============================================================================
cat(sprintf("\n=== TOTAL: %d passed, %d failed ===\n", PASS, FAIL))
if (FAIL > 0L) quit(save = "no", status = 1L)
