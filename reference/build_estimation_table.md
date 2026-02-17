# Build Estimation Properties Table from Simulation Results

Constructs a publication-quality `gt` table summarizing estimation
properties for hazard ratios in the identified subgroup and its
complement. The layout mirrors the estimation table in Leon et al.
(2024), showing average estimate, empirical SD, min, max, and relative
bias for each estimator.

## Usage

``` r
build_estimation_table(
  results,
  dgm,
  analysis_method = "FSlg",
  n_boots = NULL,
  digits = 2,
  title = "Estimation Properties",
  font_size = 12
)
```

## Arguments

- results:

  `data.table` or `data.frame`. Simulation results from
  [`run_simulation_analysis`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md),
  optionally enriched with bootstrap bias-corrected columns (`hr.H.bc`,
  `hr.Hc.bc`).

- dgm:

  DGM object. Used for true parameter values (`hr_H_true`, `hr_Hc_true`,
  and AHR truth via
  [`get_dgm_hr`](https://larry-leon.github.io/forestsearch/reference/get_dgm_hr.md)).

- analysis_method:

  Character. Which analysis method to tabulate (e.g., `"FSlg"`).
  Default: `"FSlg"`.

- n_boots:

  Integer or `NULL`. Number of bootstraps. When non-NULL, appended to
  the subtitle as "(B = n_boots bootstraps)". Default: `NULL`.

- digits:

  Integer. Decimal places. Default: 2.

- title:

  Character. Table title.

- font_size:

  Numeric. Font size in pixels for table text. Default: 12. Increase to
  14 or 16 for larger display.

## Value

A `gt` table object, or `NULL` if no estimable realizations exist.

## Details

Includes both Cox-based HR and AHR (Average Hazard Ratio from loghr_po)
estimators when AHR columns are present in the results.

For each subgroup (H and Hc) the function reports:

- **Avg**: Mean of the estimates across estimable simulations.

- **SD**: Empirical standard deviation.

- **Min / Max**: Range.

- **Rel Bias**: Relative bias in percent, `100 * (Avg - true) / true`.

When bootstrap-corrected columns (`hr.H.bc`, `hr.Hc.bc`) are present in
`results`, an additional bias-corrected row (`"theta-hat*(H)"` or
`"theta-hat*(Hc)"`) is added per subgroup.

When AHR columns (`ahr.H.hat`, `ahr.Hc.hat`) are present, AHR estimation
rows are appended using the DGM's true AHR values for relative bias
calculation.

## See also

[`build_classification_table`](https://larry-leon.github.io/forestsearch/reference/build_classification_table.md),
[`format_oc_results`](https://larry-leon.github.io/forestsearch/reference/format_oc_results.md),
[`get_dgm_hr`](https://larry-leon.github.io/forestsearch/reference/get_dgm_hr.md)

## Examples

``` r
if (FALSE) { # \dontrun{
build_estimation_table(
  results = results_alt,
  dgm = dgm_calibrated,
  analysis_method = "FSlg"
)
} # }
```
