# Plot Subgroup Results Forest Plot

Generates a comprehensive forest plot showing:

- ITT (Intent-to-Treat) population estimate

- Reference subgroups (e.g., by biomarker levels)

- Post-hoc identified subgroups with bias-corrected estimates

- Cross-validation agreement metrics as annotations

## Usage

``` r
plot_subgroup_results_forestplot(
  fs_results,
  df_analysis,
  subgroup_list = NULL,
  outcome.name,
  event.name,
  treat.name,
  E.name = "Experimental",
  C.name = "Control",
  est.scale = "hr",
  xlog = TRUE,
  title_text = NULL,
  arrow_text = c("Favors Experimental", "Favors Control"),
  footnote_text = c("Eg 80% of training found SG: 70% of B (+) also B in CV testing"),
  xlim = c(0.25, 1.5),
  ticks_at = c(0.25, 0.7, 1, 1.5),
  show_cv_metrics = TRUE,
  cv_source = c("auto", "kfold", "oob", "both"),
  posthoc_colors = c("powderblue", "beige"),
  reference_colors = c("yellow", "powderblue"),
  ci_column_spaces = 20,
  conf.level = 0.95,
  theme = NULL
)
```

## Arguments

- fs_results:

  List. A list containing ForestSearch analysis results with elements:

  - `fs.est`: ForestSearch estimation object from
    [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)

  - `fs_bc`: Bootstrap bias-corrected results from
    [`forestsearch_bootstrap_dofuture`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)

  - `fs_kfold`: K-fold cross-validation results from
    [`forestsearch_Kfold`](https://larry-leon.github.io/forestsearch/reference/forestsearch_Kfold.md)
    or
    [`forestsearch_tenfold`](https://larry-leon.github.io/forestsearch/reference/forestsearch_tenfold.md)
    (optional)

  - `fs_OOB`: Out-of-bag cross-validation results (optional, alternative
    to fs_kfold)

- df_analysis:

  Data frame. The analysis dataset with outcome, event, and treatment
  variables.

- subgroup_list:

  List. Named list of subgroup definitions to include in the plot. Each
  element should be a list with:

  - `subset_expr`: Character string for subsetting (e.g., "BM\> 1")

  - `name`: Display name for the subgroup

  - `type`: Either "reference" for pre-specified or "posthoc" for
    identified

- outcome.name:

  Character. Name of the survival time variable.

- event.name:

  Character. Name of the event indicator variable.

- treat.name:

  Character. Name of the treatment variable.

- E.name:

  Character. Label for experimental arm (default: "Experimental").

- C.name:

  Character. Label for control arm (default: "Control").

- est.scale:

  Character. Estimate scale: "hr" or "1/hr" (default: "hr").

- xlog:

  Logical. If TRUE (default), the x-axis is plotted on a logarithmic
  scale. This is standard for hazard ratio forest plots where equal
  distances represent equal relative effects.

- title_text:

  Character. Plot title (default: NULL).

- arrow_text:

  Character vector of length 2. Arrow labels for forest plot (default:
  c("Favors Experimental", "Favors Control")).

- footnote_text:

  Character vector. Footnote text for the plot explaining CV metrics
  (default provides CV interpretation guidance; set to NULL to omit).

- xlim:

  Numeric vector of length 2. X-axis limits (default: c(0.25, 1.5)).

- ticks_at:

  Numeric vector. X-axis tick positions (default: c(0.25, 0.70, 1.0,
  1.5)).

- show_cv_metrics:

  Logical. Whether to show cross-validation metrics (default: TRUE if
  fs_kfold or fs_OOB available).

- cv_source:

  Character. Source for CV metrics: "auto" (default, uses both if
  available, otherwise whichever is present), "kfold" (use fs_kfold
  only), "oob" (use fs_OOB only), or "both" (explicitly use both
  fs_kfold and fs_OOB, with K-fold first then OOB).

- posthoc_colors:

  Character vector. Colors for post-hoc subgroup rows (default:
  c("powderblue", "beige")).

- reference_colors:

  Character vector. Colors for reference subgroup rows (default:
  c("yellow", "powderblue")).

- ci_column_spaces:

  Integer. Number of spaces for the CI plot column width. More spaces =
  wider CI column (default: 20).

- conf.level:

  Numeric. Confidence level for intervals (default: 0.95 for 95% CI).
  Used to calculate the z-multiplier as qnorm(1 - (1 - conf.level)/2).

- theme:

  An fs_forest_theme object from
  [`create_forest_theme()`](https://larry-leon.github.io/forestsearch/reference/create_forest_theme.md).
  Use this to control plot sizing (fonts, row height, CI appearance, CV
  annotation font size). Default: NULL (uses default theme).

## Value

A list containing:

- plot:

  The forestploter grob object (can be rendered with plot())

- data:

  The data frame used for the forest plot

- row_types:

  Character vector of row types for styling reference

- cv_metrics:

  Cross-validation metrics text (if available)

## Details

Creates a publication-ready forest plot displaying identified subgroups
with hazard ratios, bias-corrected estimates, and cross-validation
metrics. This wrapper integrates ForestSearch results with the
forestploter package.

### ForestSearch Labeling Convention

ForestSearch identifies subgroups based on hazard ratio thresholds:

- `sg.harm`: Contains the definition of the "harm" or "questionable"
  subgroup (H)

- `treat.recommend == 0`: Patient is IN the harm subgroup (H)

- `treat.recommend == 1`: Patient is in the COMPLEMENT subgroup (Hc,
  typically benefit)

For `est.scale = "hr"` (searching for harm):

- H (treat.recommend=0): Subgroup defined by sg.harm with elevated HR
  (harm/questionable)

- Hc (treat.recommend=1): Complement of sg.harm (potential benefit)

For `est.scale = "1/hr"` (searching for benefit):

- Roles are reversed: H becomes the benefit group

## See also

[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
for running the subgroup analysis
[`forestsearch_bootstrap_dofuture`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)
for bootstrap bias correction
[`forestsearch_Kfold`](https://larry-leon.github.io/forestsearch/reference/forestsearch_Kfold.md)
for cross-validation
[`create_forest_theme`](https://larry-leon.github.io/forestsearch/reference/create_forest_theme.md)
for customizing plot appearance
[`render_forestplot`](https://larry-leon.github.io/forestsearch/reference/render_forestplot.md)
for rendering the plot

## Examples

``` r
if (FALSE) { # \dontrun{
# Load ForestSearch results
load("fs_analysis_results.Rdata")  # Contains fs.est, fs_bc, fs_kfold

# Define subgroups to display
subgroups <- list(
  bm_gt1 = list(
    subset_expr = "BM > 1",
    name = "BM > 1",
    type = "reference"
  ),
  bm_gt1_size_gt19 = list(
    subset_expr = "BM > 1 & tmrsize > 19",
    name = "BM > 1 & Tumor Size > 19",
    type = "posthoc"
  )
)

# Create the forest plot with default theme
result <- plot_subgroup_results_forestplot(
  fs_results = list(fs.est = fs.est, fs_bc = fs_bc, fs_kfold = fs_kfold),
  df_analysis = df_itt,
  subgroup_list = subgroups,
  outcome.name = "os_time",
  event.name = "os_event",
  treat.name = "combo",
  E.name = "Experimental+CT",
  C.name = "CT"
)

# Create with custom theme for larger plot
large_theme <- create_forest_theme(
  base_size = 14,
  row_padding = c(6, 4),
  cv_fontsize = 11
)

result <- plot_subgroup_results_forestplot(
  fs_results = list(fs.est = fs.est, fs_bc = fs_bc, fs_kfold = fs_kfold),
  df_analysis = df_itt,
  subgroup_list = subgroups,
  outcome.name = "os_time",
  event.name = "os_event",
  treat.name = "combo",
  theme = large_theme
)

# Display the plot
render_forestplot(result)
} # }
```
