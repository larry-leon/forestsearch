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
  save_plots = FALSE,
  output_dir = "plots",
  verbose = TRUE
)
```

## Arguments

- df:

  Data frame containing survival data with potential outcomes

- tte_name:

  Character string specifying time-to-event variable name. Default:
  "os_time"

- event_name:

  Character string specifying event indicator variable name. Default:
  "os_event"

- treat_name:

  Character string specifying treatment variable name. Default: "treat"

- z_name:

  Character string specifying continuous covariate/biomarker name.
  Default: "biomarker"

- loghr_po_name:

  Character string specifying potential outcome log HR variable.
  Default: "loghr_po"

- theta1_name:

  Optional: variable name for theta_1 (treated potential outcome).
  Default: "theta_1"

- theta0_name:

  Optional: variable name for theta_0 (control potential outcome).
  Default: "theta_0"

- spline_df:

  Integer degrees of freedom for spline fitting. Default: 3

- alpha:

  Numeric significance level for confidence intervals. Default: 0.20

- hr_threshold:

  Numeric hazard ratio threshold for subgroup identification. Default:
  0.7

- plot_style:

  Character: "combined", "separate", or "grid" for plot layout. Default:
  "combined"

- save_plots:

  Logical whether to save plots to file. Default: FALSE

- output_dir:

  Character directory for saving plots. Default: "plots"

- verbose:

  Logical for diagnostic output. Default: TRUE

## Value

List containing:

- cox_fit:

  Results from cox_cs_fit function

- ahr_results:

  AHR calculations for different subgroup definitions

- cde_results:

  CDE calculations if theta variables available

- optimal_cutpoint:

  Optimal biomarker cutpoint for subgroup

- subgroup_stats:

  Statistics for recommended and questionable subgroups

- plots:

  List of generated plot objects (if using ggplot2)

## Examples

``` r
if (FALSE) { # \dontrun{
# Load your ForestSearch output data
results <- cox_ahr_cde_analysis(
  df = df_nonAP,
  z_name = "biomarker",
  hr_threshold = 0.7,
  plot_style = "grid"
)

# Access specific results
print(results$subgroup_stats)
print(results$optimal_cutpoint)
} # }
```
