# Build Classification Rate Table from Simulation Results

Constructs a publication-quality `gt` table summarizing subgroup
identification and classification rates across one or more data
generation scenarios and analysis methods. The layout mirrors Table 4 of
Leon et al. (2024) with metrics grouped by model scenario (null / alt)
and columns for each analysis method.

## Usage

``` r
build_classification_table(
  scenario_results,
  analyses = NULL,
  digits = 2,
  title = "Subgroup Identification and Classification Rates",
  n_sims = NULL,
  bold_threshold = 0.05,
  font_size = 12
)
```

## Arguments

- scenario_results:

  Named list. Each element is itself a list with:

  results

  :   `data.table` from
      [`run_simulation_analysis`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md).

  label

  :   Character scenario label, e.g., `"M1"`.

  n_sample

  :   Integer sample size.

  dgm

  :   DGM object (for true HRs and subgroup prevalence).

  hypothesis

  :   Character: `"null"` or `"alt"`.

- analyses:

  Character vector of analysis labels to include (e.g.,
  `c("FS", "FSlg", "GRF")`). When `NULL`, all unique values of
  `results$analysis` across scenarios are used.

- digits:

  Integer. Decimal places for proportions. Default: 2.

- title:

  Character. Table title. Default:
  `"Subgroup Identification and Classification Rates"`.

- n_sims:

  Integer. Number of simulations (for subtitle). Default: `NULL`.

- bold_threshold:

  Numeric. Type I error threshold above which the `any(H)` value is
  shown in bold. Set `NULL` to disable. Default: 0.05.

- font_size:

  Numeric. Font size in pixels for table text. Default: 12. Increase to
  14 or 16 for larger display.

## Value

A `gt` table object.

## Details

For each scenario the function computes:

- `any(H)`: Proportion of simulations identifying any subgroup.

- `sens(H)`: Mean sensitivity (only under alternative).

- `sens(Hc)`: Mean specificity.

- `ppv(H)`: Mean positive predictive value (only under alternative).

- `ppv(Hc)`: Mean negative predictive value.

- `avg|H|`: Mean size of identified subgroup (when found).

Under the null hypothesis the rows are reduced to `any(H)`, `sens(Hc)`,
`ppv(Hc)`, and `avg|H|`.

## See also

[`format_oc_results`](https://larry-leon.github.io/forestsearch/reference/format_oc_results.md),
[`summarize_simulation_results`](https://larry-leon.github.io/forestsearch/reference/summarize_simulation_results.md)
