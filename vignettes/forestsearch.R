## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(

  echo = TRUE,
  warning = FALSE,
  message = FALSE,

  fig.align = "center",
  fig.retina = 2,
  fig.width = 8,
  fig.height = 6
)

# Bootstrap 5 code-folding hook for pkgdown compatibility
# Usage: add `echo=FALSE, code_fold=TRUE` to any chunk header
knitr::knit_hooks$set(code_fold = function(before, options, envir) {
  if (!before && isTRUE(options$code_fold)) {
    id <- gsub("[^a-zA-Z0-9]", "", options$label)
    code_lines <- knitr::knit_code$get(options$label)
    code_text  <- htmltools::htmlEscape(paste(code_lines, collapse = "\n"))
    sprintf(
      '<p><a class="btn btn-outline-secondary btn-sm" data-bs-toggle="collapse" href="#collapse-%s" role="button" aria-expanded="false" aria-controls="collapse-%s"><i class="bi bi-code-slash"></i> Code</a></p>\n<div class="collapse" id="collapse-%s"><div class="card card-body p-0"><pre class="r"><code class="hljs">%s</code></pre></div></div>',
      id, id, id, code_text
    )
  }
})

# Initialize timing
timings <- list()
t_vignette_start <- proc.time()


## ----load-packages------------------------------------------------------------
library(forestsearch)
library(survival)
library(data.table)
library(ggplot2)
library(gt)
library(grf)
library(policytree)
library(doFuture)

# Optional packages for enhanced output
library(patchwork)
library(weightedsurv)

# Set ggplot theme
theme_set(theme_minimal(base_size = 12))


## ----data-setup---------------------------------------------------------------
# Load GBSG data (included in forestsearch package)
df.analysis <- gbsg

# Prepare analysis variables
df.analysis <- within(df.analysis, {
  id <- seq_len(nrow(df.analysis))
  time_months <- rfstime / 30.4375
  grade3 <- ifelse(grade == "3", 1, 0)
  treat <- hormon
})

# Define variable roles
confounders.name <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")
outcome.name <- "time_months"
event.name <- "status"
id.name <- "id"
treat.name <- "hormon"

# Display data structure
cat("Sample size:", nrow(df.analysis), "\n")
cat("Events:", sum(df.analysis[[event.name]]), 
    sprintf("(%.1f%%)\n", 100 * mean(df.analysis[[event.name]])))
cat("Baseline factors:", paste(confounders.name, collapse = ", "), "\n")


## ----baseline-table-----------------------------------------------------------
create_summary_table(
  data = df.analysis,
  treat_var = treat.name,
  table_title = "GBSG Baseline Characteristics by Treatment Arm",
  vars_continuous = c("age", "nodes", "size", "er", "pgr"),
  vars_categorical = c("grade", "meno"),
  font_size = 12
)


## ----km-itt, fig.width=8, fig.height=5----------------------------------------
# Prepare counting process data for KM plot
dfcount <- df_counting(
  df = df.analysis,
  by.risk = 6,
  tte.name = outcome.name,
  event.name = event.name,
  treat.name = treat.name
)

# Plot with confidence intervals and log-rank test
plot_weighted_km(
  dfcount,
  conf.int = TRUE,
  show.logrank = TRUE,
  ymax = 1.05,
  xmed.fraction = 0.775,
  ymed.offset = 0.125
)


## ----grf-analysis-------------------------------------------------------------
t0 <- proc.time()

grf_est <- grf.subg.harm.survival(
  data = df.analysis,
  confounders.name = confounders.name,
  outcome.name = outcome.name,
  event.name = event.name,
  id.name = id.name,
  treat.name = treat.name,
  maxdepth = 2,
  n.min = 60,
  dmin.grf = 12,
  frac.tau = 0.6,
  details = TRUE,
  return_selected_cuts_only = FALSE
)

timings$grf <- (proc.time() - t0)["elapsed"]


## ----grf-trees, fig.width=10, fig.height=4------------------------------------
# Display policy trees
# leaf1 = recommend control, leaf2 = recommend treatment
oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar), add = TRUE)
par(mfrow = c(1, 2))
plot(grf_est$tree1, leaf.labels = c("Control", "Treat"), main = "Depth 1")
plot(grf_est$tree2, leaf.labels = c("Control", "Treat"), main = "Depth 2")


## ----parallel-setup-----------------------------------------------------------
# Detect available cores (limited to 2 cores for CRAN checks)
n_cores <- 2
n_cores_total <- parallel::detectCores()
cat("Using", n_cores, "of", n_cores_total, "total cores for parallel processing")


## ----forestsearch-main, fig.width=10, fig.height=7----------------------------
t0 <- proc.time()

fs <- forestsearch(
  df.analysis,
  confounders.name = confounders.name,
  outcome.name = outcome.name,
  treat.name = treat.name,
  event.name = event.name,
  id.name = id.name,
  # Threshold parameters (per León et al. 2024)
  hr.threshold = 1.25,
  hr.consistency = 1.0,
  pconsistency.threshold = 0.80,
  stop_threshold = 0.80,
  # Search configuration
  sg_focus = "hr",
  max_subgroups_search = 3,
  use_twostage = TRUE,
  # Factor selection
  use_grf = TRUE, 
  return_selected_cuts_only = TRUE,
  use_lasso = TRUE,
  cut_type = "default",
  # Subgroup constraints
  maxk = 2,
  n.min = 60,
  d0.min = 12,
  d1.min = 12,
  # Consistency evaluation
  fs.splits = 100,
  # Parallel processing
  parallel_args = list(
    plan = "multisession",
    workers = n_cores,
    show_message = TRUE
  ),
  # Output options
  showten_subgroups = TRUE,
  details = TRUE,
  plot.sg = TRUE
)

plan("sequential")
timings$forestsearch <- (proc.time() - t0)["elapsed"]

cat("\nForestSearch completed in", 
    round(timings$forestsearch, 1), "seconds\n")


## ----fs-results---------------------------------------------------------------
# Generate results tables
res_tabs <- sg_tables(fs, ndecimals = 3, which_df = "est")

# Display top subgroups meeting criteria
res_tabs$sg10_out


## ----fs-estimates-------------------------------------------------------------
# ITT and subgroup estimates
res_tabs$tab_estimates


## ----fs-subgroup--------------------------------------------------------------
cat("Identified subgroup (H):", paste(fs$sg.harm, collapse = " & "), "\n")
cat("Subgroup size:", sum(fs$df.est$treat.recommend == 0), 
    sprintf("(%.1f%% of ITT)\n", 
            100 * mean(fs$df.est$treat.recommend == 0)))


## ----bootstrap, eval=TRUE-----------------------------------------------------
# Number of bootstrap iterations
# Use 500-2000 for production; reduced here for vignette
NB <- 2

t0 <- proc.time()

fs_bc <- forestsearch_bootstrap_dofuture(
  fs.est = fs,
  nb_boots = NB,
  show_three = FALSE,
  details = FALSE,
  parallel_args = list(
    plan = "multisession",
    workers = n_cores,
    show_message = TRUE
  )
)

plan("sequential")
timings$bootstrap <- (proc.time() - t0)["elapsed"]

cat("\nBootstrap completed in", 
    round(timings$bootstrap / 60, 1), "minutes\n")


## ----bootstrap-summary--------------------------------------------------------
# Comprehensive summary with diagnostics
summaries <- summarize_bootstrap_results(
  sgharm = fs$sg.harm,
  boot_results = fs_bc,
  create_plots = TRUE,
  est.scale = "hr"
)

# Display bias-corrected estimates table
summaries$table


## ----fs-results_figure, fig.width=10, fig.height=7, fig.cap="Kaplan-Meier survival curves by identified subgroup"----

 km_result <- plot_sg_weighted_km(
   fs.est = fs,
   outcome.name = "time_months",
   event.name = "status",
   treat.name = "hormon",
   show.logrank = FALSE,
   conf.int = TRUE,
   by.risk = 12,
   show.cox = FALSE, show.cox.bc = TRUE,
   fs_bc = fs_bc,
   hr_bc_position = "topright"
 )


## ----event-summary------------------------------------------------------------
# note that default required minimum events is 12 for subgroup candidate
# Here we evaluate frequency of subgroup candidates in bootstrap samples less than 15
event_summary <- summarize_bootstrap_events(fs_bc, threshold = 15)


## ----bootstrap-diagnostics----------------------------------------------------
# Quality metrics
summaries$diagnostics_table_gt


## ----subgroup-agreement-------------------------------------------------------
# Agreement with original analysis
if (!is.null(summaries$subgroup_summary$original_agreement)) {
  summaries$subgroup_summary$original_agreement
}

# Factor presence across bootstrap iterations
if (!is.null(summaries$subgroup_summary$factor_presence)) {
  summaries$subgroup_summary$factor_presence
}


## ----bootstrap-plots, fig.width=10, fig.height=4, eval=TRUE-------------------
if (!is.null(summaries$plots)) {
  summaries$plots$H_distribution + summaries$plots$Hc_distribution
}


## ----kfold-cv, eval = TRUE----------------------------------------------------
# 10-fold CV with multiple iterations
# Use Ksims >= 50 for production
Ksims <- 1

t0 <- proc.time()

fs_kfold <- forestsearch_tenfold(
  fs.est = fs,
  sims = Ksims,
  Kfolds = 2,
  details = FALSE,
  parallel_args = list(
    plan = "multisession",
    workers = n_cores,
    show_message = FALSE
  )
)

plan("sequential")
timings$kfold <- (proc.time() - t0)["elapsed"]
metrics_tables <- cv_metrics_tables(fs_kfold)
metrics_tables



## ----oob-cv, eval = FALSE-----------------------------------------------------
# t0 <- proc.time()
# 
# fs_OOB <- forestsearch_Kfold(
#   fs.est = fs,
#   details = FALSE,
#   Kfolds = round(nrow(df.analysis)/100,0),  # N-fold = leave-one-out
#   parallel_args = list(
#     plan = "multisession",
#     workers = n_cores,
#     show_message = TRUE
#   )
# )
# 
# plan("sequential")
# timings$oob <- (proc.time() - t0)["elapsed"]
# 
# # Summarize OOB results
# cv_out <- forestsearch_KfoldOut(
#   res = fs_OOB,
#   details = FALSE,
#   outall = TRUE
# )
# 
# tables <- cv_summary_tables(cv_out)
# 
# tables$combined_table
# 
# tables$metrics_table
# 
# 
# 
# 


## ----forest-plot, fig.width=18, fig.height=12, fig.cap="Subgroup forest plot including identified subgroups"----


# Define reference subgroups for comparison
subgroups <- list(
  age_gt65 = list(
    subset_expr = "age > 65",
    name = "Age > 65",
    type = "reference"
  ),
  age_le65 = list(
    subset_expr = "age <= 65",
    name = "Age ≤ 65",
    type = "reference"
  ),
  pgr_positive = list(
    subset_expr = "pgr > 0",
    name = "PgR > 0",
    type = "reference"
  ),
  pgr_negative = list(
    subset_expr = "pgr <= 0",
    name = "PgR ≤ 0",
    type = "reference"
  )
)


my_theme <- create_forest_theme(base_size = 24, 
footnote_fontsize = 17, cv_fontsize = 22)


# Create forest plot
# Include fs_kfold and fs_OOB if available for CV metrics
result <- plot_subgroup_results_forestplot(
  fs_results = list(
    fs.est = fs,
    fs_bc = fs_bc,
    fs_OOB = NULL,
    fs_kfold = fs_kfold
  ),
  df_analysis = df.analysis,
  subgroup_list = subgroups,
  outcome.name = outcome.name,
  event.name = event.name,
  treat.name = treat.name,
  E.name = "Hormonal",
  C.name = "Chemo",
  ci_column_spaces = 25,
  xlog = TRUE,
  theme = my_theme
)


# Option 2: Custom sizing
render_forestplot(result)   



## ----KMdiffs, fig.width = 8, fig.height = 6, fig.align="center"---------------

# Add additional subgroups along with ITT and identified subgroups
ref_sgs <- list(
age_young = list(subset_expr = "age < 65", color = "brown"),
age_old = list(subset_expr = "age >= 65", color = "orange")
)

plot_km_band_forestsearch(
 df = df.analysis,
   fs.est = fs,
 ref_subgroups = ref_sgs,
 outcome.name = outcome.name,
   event.name = event.name,
   treat.name = treat.name,
 draws_band = 20
)
 
# # Example with more subgroups
# ref_sgs <- list(
# pgr_positive = list(subset_expr = "pgr > 0", color ="green"),
# pgr_negative = list(subset_expr = "pgr <= 0", color = "purple"),
# age_young = list(subset_expr = "age < 65", color = "brown"),
# age_old = list(subset_expr = "age >= 65", color = "orange")
# )



## ----summary-findings---------------------------------------------------------
# Extract key results
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("FORESTSEARCH ANALYSIS SUMMARY\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n\n")

cat("Dataset: GBSG (N =", nrow(df.analysis), ")\n")
cat("Outcome: Recurrence-free survival\n\n")

cat("ITT Analysis:\n")
cat("  HR (95% CI): 0.69 (0.54, 0.89)\n\n")

cat("Identified Subgroup (H):\n")
cat("  Definition:", paste(fs$sg.harm, collapse = " & "), "\n")
cat("  Size:", sum(fs$df.est$treat.recommend == 0), 
    sprintf("(%.1f%%)\n", 100 * mean(fs$df.est$treat.recommend == 0)))
cat("  Unadjusted HR:", sprintf("%.2f", exp(fs$grp.consistency$out_sg$result$hr[1])), "\n")

cat("\nComplement Subgroup (Hc):\n")
cat("  Size:", sum(fs$df.est$treat.recommend == 1),
    sprintf("(%.1f%%)\n", 100 * mean(fs$df.est$treat.recommend == 1)))


## ----timing-summary, echo=FALSE, code_fold=TRUE-------------------------------
timings$total <- (proc.time() - t_vignette_start)["elapsed"]

timing_df <- data.frame(
  Analysis = c("GRF", "ForestSearch", "Bootstrap", "Total"),
  Seconds = c(
    timings$grf,
    timings$forestsearch,
    timings$bootstrap,
    timings$total
  )
)
timing_df$Minutes <- timing_df$Seconds / 60

gt(timing_df) |>
  tab_header(title = "Computational Timing") |>
  fmt_number(columns = c(Seconds, Minutes), decimals = 1) |>
  cols_label(
    Analysis = "Component",
    Seconds = "Time (sec)",
    Minutes = "Time (min)"
  )


## ----session-info-------------------------------------------------------------
sessionInfo()

