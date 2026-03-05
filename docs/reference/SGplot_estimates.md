# Violin/Boxplot Visualization of HR Estimates

Creates violin plots with embedded boxplots showing the distribution of
hazard ratio estimates across simulations for different analysis
populations. Supports symmetric trimming to handle extreme values that
can distort the display when small subgroups produce very large HR
estimates.

## Usage

``` r
SGplot_estimates(
  df,
  label_training = "Training",
  label_testing = "Testing",
  label_itt = "ITT (stratified)",
  label_sg = "Testing (subgroup)",
  trim_fraction = NULL,
  ylim = NULL,
  show_summary = NULL,
  title = "Distribution of HR Estimates Across Simulations",
  subtitle = NULL
)
```

## Arguments

- df:

  data.frame or data.table. Simulation results from
  [`mrct_region_sims`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)

- label_training:

  Character. Label for training data estimates. Default: "Training"

- label_testing:

  Character. Label for testing data estimates. Default: "Testing"

- label_itt:

  Character. Label for ITT estimates. Default: "ITT (stratified)"

- label_sg:

  Character. Label for subgroup estimates. Default: "Testing (subgroup)"

- trim_fraction:

  Numeric or NULL. Fraction of observations to trim from each tail
  (e.g., 0.01 trims the lowest 1\\ When non-NULL, trimmed means and SDs
  are computed for each group, extreme observations are flagged, and the
  y-axis is clipped to the trimmed data range. Set to NULL (default) for
  no trimming (backward compatible).

- ylim:

  Numeric vector of length 2 or NULL. Explicit y-axis limits as
  `c(lower, upper)`. Overrides automatic limits from trimming. Default:
  NULL (auto).

- show_summary:

  Logical. Annotate each violin with mean (SD) below the x-axis labels.
  When trimming is active, displays trimmed statistics. Default: TRUE
  when `trim_fraction` is non-NULL, FALSE otherwise.

- title:

  Character. Plot title. Default: "Distribution of HR Estimates Across
  Simulations".

- subtitle:

  Character or NULL. Plot subtitle. When trimming is active and subtitle
  is NULL, an auto-generated note indicating the trim fraction and
  number of flagged observations is shown. Default: NULL.

## Value

List with components:

- dfPlot_estimates:

  data.table formatted for plotting, with a `trimmed` logical column
  when trimming is active

- plot_estimates:

  ggplot2 object

- trim_info:

  List of per-group trimming diagnostics (NULL when no trimming). Each
  element contains: `n_total`, `n_trimmed`, `n_flagged`, `raw_mean`,
  `raw_sd`, `trimmed_mean`, `trimmed_sd`, `lower_bound`, `upper_bound`.

## See also

[`mrct_region_sims`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)
for generating simulation results,
[`summaryout_mrct`](https://larry-leon.github.io/forestsearch/reference/summaryout_mrct.md)
for tabular summaries with trimming

## Examples

``` r
if (FALSE) { # \dontrun{
# Default (no trimming) — backward compatible
plot_results <- SGplot_estimates(
  results_alt,
  label_training = "Non-Region A, ITT",
  label_itt = "Overall, ITT",
  label_testing = "Region A, ITT",
  label_sg = "Region A, identified subgroup"
)
print(plot_results$plot_estimates)

# With 1% symmetric trimming
plot_results <- SGplot_estimates(
  results_alt,
  label_training = "Non-AP, ITT",
  label_itt = "Overall, ITT",
  label_testing = "AP, ITT",
  label_sg = "AP, identified subgroup",
  trim_fraction = 0.01
)
print(plot_results$plot_estimates)
# Inspect trimming diagnostics
print(plot_results$trim_info)

# With explicit y-axis limits
plot_results <- SGplot_estimates(
  results_alt,
  trim_fraction = 0.01,
  ylim = c(0.2, 3.0)
)
} # }
```
