# Violin/Boxplot Visualization of HR Estimates

Creates violin plots with embedded boxplots showing the distribution of
hazard ratio estimates across simulations for different analysis
populations.

## Usage

``` r
SGplot_estimates(
  df,
  label_training = "Training",
  label_testing = "Testing",
  label_itt = "ITT (stratified)",
  label_sg = "Testing (subgroup)"
)
```

## Arguments

- df:

  data.frame or data.table. Simulation results from
  [`mrct_region_sims`](https://github.com/larry-leon/forestsearch/reference/mrct_region_sims.md)

- label_training:

  Character. Label for training data estimates. Default: "Training"

- label_testing:

  Character. Label for testing data estimates. Default: "Testing"

- label_itt:

  Character. Label for ITT estimates. Default: "ITT (stratified)"

- label_sg:

  Character. Label for subgroup estimates. Default: "Testing (subgroup)"

## Value

List with components:

- dfPlot_estimates:

  data.table formatted for plotting

- plot_estimates:

  ggplot2 object

## See also

[`mrct_region_sims`](https://github.com/larry-leon/forestsearch/reference/mrct_region_sims.md)
for generating simulation results

## Examples

``` r
if (FALSE) { # \dontrun{
# After running simulations
plot_results <- SGplot_estimates(
  results_alt,
  label_training = "Non-Region A, ITT",
  label_itt = "Overall, ITT",
  label_testing = "Region A, ITT",
  label_sg = "Region A, identified subgroup"
)
print(plot_results$plot_estimates)
} # }
```
