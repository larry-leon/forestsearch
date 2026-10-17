# Comprehensive Wrapper for Cox Spline Analysis with AHR and CDE Plotting

This wrapper function combines Cox spline fitting with comprehensive
visualization of Average Hazard Ratios (AHRs) and Controlled Direct
Effects (CDEs) as described in the MRCT subgroups analysis
documentation.

## Usage

``` r
cox_ahr_cde_analysis(
  df,
  tte_name = "os_time",
  event_name = "os_event",
  treat_name = "treat",
  z_name = "biomarker",
  loghr_po_name = "loghr_po",
  theta1_name = "theta_1",
  theta0_name = "theta_0",
  spline_df = 3,
  alpha = 0.2,
  hr_threshold = 0.7,
  plot_style = c("combined", "separate", "grid"),
  plot_select = c("all", "profile_ahr", "ahr_only"),
  save_plots = FALSE,
  output_dir = "plots",
  verbose = TRUE
)
```

## Arguments

- df:

  Data frame containing survival data with potential outcomes.

- tte_name:

  Character string specifying time-to-event variable name. Default:
  `"os_time"`.

- event_name:

  Character string specifying event indicator variable name. Default:
  `"os_event"`.

- treat_name:

  Character string specifying treatment variable name. Default:
  `"treat"`.

- z_name:

  Character string specifying continuous covariate/biomarker name.
  Default: `"biomarker"`.

- loghr_po_name:

  Character string specifying potential outcome log HR variable.
  Default: `"loghr_po"`.

- theta1_name:

  Optional: variable name for theta_1 (treated potential outcome).
  Default: `"theta_1"`.

- theta0_name:

  Optional: variable name for theta_0 (control potential outcome).
  Default: `"theta_0"`.

- spline_df:

  Integer degrees of freedom for spline fitting. Default: 3.

- alpha:

  Numeric significance level for confidence intervals. Default: 0.20.

- hr_threshold:

  Numeric hazard ratio threshold for subgroup identification, or `NULL`
  to suppress threshold-based elements. When `NULL`, plots omit the HR
  threshold line, optimal cutpoint line, and subgroup-dependent panels
  (Plots 5–6); subgroup statistics report only the overall population.
  Default: 0.7.

- plot_style:

  Character: `"combined"`, `"separate"`, or `"grid"` for plot layout.
  Default: `"combined"`.

- plot_select:

  Character controlling which panels to display: `"all"` (default) shows
  the full grid layout; `"profile_ahr"` shows only the treatment effect
  profile (top-left) and the AHR curve for \\z \geq\\ threshold
  (top-middle) side by side, with the AHR panel's y-axis scaled to the
  data range; `"ahr_only"` shows only the AHR curve for \\z \geq\\
  threshold as a single panel with self-contained y-axis.

- save_plots:

  Logical whether to save plots to file. Default: FALSE.

- output_dir:

  Character directory for saving plots. Default: `"plots"`.

- verbose:

  Logical for diagnostic output. Default: TRUE.

## Value

List of class `"cox_ahr_cde"` containing:

- cox_fit:

  Results from
  [`cox_cs_fit`](https://larry-leon.github.io/forestsearch/reference/cox_cs_fit.md)
  function.

- ahr_results:

  AHR calculations for different subgroup definitions.

- cde_results:

  CDE calculations if theta variables available.

- optimal_cutpoint:

  Optimal biomarker cutpoint, or `NULL` when `hr_threshold` is `NULL`.

- subgroup_stats:

  Statistics for recommended and questionable subgroups, or overall-only
  when `hr_threshold` is `NULL`.

- data:

  List with z_values, loghr_po, and subgroup assignments.

## Examples

``` r
if (FALSE) { # \dontrun{
# With threshold - full subgroup analysis
results <- cox_ahr_cde_analysis(
  df = df_large, z_name = "z_bm",
  hr_threshold = 1.25, plot_style = "grid"
)

# Without threshold - pure AHR/CDE curves
results <- cox_ahr_cde_analysis(
  df = df_large, z_name = "z_bm",
  hr_threshold = NULL, plot_style = "grid"
)

# Compact two-panel without threshold
results <- cox_ahr_cde_analysis(
  df = df_large, z_name = "z_bm",
  hr_threshold = NULL,
  plot_select = "profile_ahr"
)

# Single AHR panel only
results <- cox_ahr_cde_analysis(
  df = df_large, z_name = "z_bm",
  hr_threshold = 1.25,
  plot_select = "ahr_only"
)
} # }
```
