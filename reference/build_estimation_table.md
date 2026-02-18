# Build Estimation Properties Table from Simulation Results

Constructs a publication-quality `gt` table summarizing estimation
properties for hazard ratios in the identified subgroup and its
complement. The layout mirrors Table 5 of Leon et al. (2024), showing
average estimate, empirical SD, min, max, and relative bias for each
estimator.

## Usage

``` r
build_estimation_table(
  results,
  dgm,
  analysis_method = "FSlg",
  n_boots = NULL,
  digits = 2,
  title = "Estimation Properties",
  subtitle = NULL,
  font_size = 12,
  cde_H = NULL,
  cde_Hc = NULL
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

- subtitle:

  Character or `NULL`. Optional user-supplied subtitle text.
  Auto-populated statistics (method, estimable count, proportion) are
  always shown. When non-NULL, the custom text is displayed first with
  stats appended in parentheses, e.g. "My title (FSlg: 18/20 (90%)
  estimable)".

- font_size:

  Numeric. Font size in pixels for table text. Default: 12. Increase to
  14 or 16 for larger display.

- cde_H:

  Numeric or `NULL`. Controlled direct effect (theta-ddagger(H)) for the
  true harm subgroup. When non-NULL, an additional b-ddagger bias column
  is shown. When `NULL` (default), auto-detected from `dgm$cde_H` or
  `dgm$hazard_ratios$CDE_harm`.

- cde_Hc:

  Numeric or `NULL`. Controlled direct effect for the complement.
  Auto-detected analogously.

## Value

A `gt` table object, or `NULL` if no estimable realizations exist.

## Details

Uses the paper's notation conventions:

- theta-dagger: Marginal (causal) HR truth

- theta-ddagger: Controlled direct effect (CDE) truth

- theta-hat(H-hat): Plugin Cox estimate in identified subgroup

- theta-hat\*(H-hat): Bootstrap bias-corrected estimate

Includes both Cox-based HR and AHR (Average Hazard Ratio from loghr_po)
estimators when AHR columns are present in the results.

For each subgroup (H and Hc) the function reports:

- **Avg**: Mean of the estimates across estimable simulations.

- **SD**: Empirical standard deviation.

- **Min / Max**: Range.

- **b-dagger**: Relative bias (percent) vs marginal truth,
  `100 * (Avg - theta_dagger) / theta_dagger`.

- **b-ddagger** (conditional): Relative bias (percent) vs CDE truth,
  shown when CDE values are available.

When bootstrap-corrected columns (`hr.H.bc`, `hr.Hc.bc`) are present in
`results`, an additional bias-corrected row (theta-hat\*(H-hat)) is
added per subgroup.

When AHR columns (`ahr.H.hat`, `ahr.Hc.hat`) are present, AHR estimation
rows are appended using the DGM's true AHR values for relative bias
calculation.

When CDE columns (`cde.H.hat`, `cde.Hc.hat`) are present and CDE truth
values are available, CDE estimation rows (theta-ddagger(H-hat)) are
appended. The b-dagger column for CDE rows reports bias relative to the
CDE truth rather than the marginal HR.

## See also

[`build_classification_table`](https://larry-leon.github.io/forestsearch/reference/build_classification_table.md),
[`format_oc_results`](https://larry-leon.github.io/forestsearch/reference/format_oc_results.md),
[`get_dgm_hr`](https://larry-leon.github.io/forestsearch/reference/get_dgm_hr.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage (marginal truth only)
build_estimation_table(
  results = results_alt,
  dgm = dgm_calibrated,
  analysis_method = "FSlg"
)

# With CDE truth (full Table 5 alignment)
build_estimation_table(
  results = results_alt,
  dgm = dgm_calibrated,
  cde_H = 2.25,
  cde_Hc = 0.60
)
} # }
```
