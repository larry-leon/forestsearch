# ForestSearch Data Preparation and Feature Selection

Prepares a dataset for ForestSearch, including options for LASSO-based
dimension reduction, GRF cuts, forced cuts, and flexible cut strategies.
Returns a list with the processed data, subgroup factor names, cut
expressions, and LASSO selection results.

## Usage

``` r
get_FSdata(
  df.analysis,
  use_lasso = FALSE,
  use_grf = FALSE,
  grf_cuts = NULL,
  confounders.name,
  cont.cutoff = 4,
  conf_force = NULL,
  conf.cont_medians = NULL,
  conf.cont_medians_force = NULL,
  replace_med_grf = TRUE,
  defaultcut_names = NULL,
  cut_type = "default",
  exclude_cuts = NULL,
  outcome.name = "tte",
  event.name = "event",
  details = TRUE
)
```

## Arguments

- df.analysis:

  Data frame containing the data.

- use_lasso:

  Logical. Whether to use LASSO for dimension reduction.

- use_grf:

  Logical. Whether to use GRF cuts.

- grf_cuts:

  Character vector of GRF cut expressions.

- confounders.name:

  Character vector of confounder variable names.

- cont.cutoff:

  Integer. Cutoff for continuous variable determination.

- conf_force:

  Character vector of forced cut expressions.

- conf.cont_medians:

  Character vector of continuous confounders to cut at median.

- conf.cont_medians_force:

  Character vector of additional continuous confounders to force median
  cut.

- replace_med_grf:

  Logical. If TRUE, removes median cuts that overlap with GRF cuts.

- defaultcut_names:

  Character vector of confounders to force default cuts.

- cut_type:

  Character. "default" or "median" for cut strategy.

- exclude_cuts:

  Character vector of cut expressions to exclude.

- outcome.name:

  Character. Name of outcome variable.

- event.name:

  Character. Name of event indicator variable.

- details:

  Logical. If TRUE, prints details during execution.
