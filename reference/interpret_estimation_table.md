# Generate Narrative Interpretation of Estimation Properties

Produces a templated text summary of the estimation properties table,
automatically populating numerical results from the simulation output.
Useful for reproducible vignettes where interpretation paragraphs should
update when simulations are re-run.

## Usage

``` r
interpret_estimation_table(
  results,
  dgm,
  analysis_method = "FSlg",
  n_sims = NULL,
  n_boots = 300,
  digits = 2,
  scenario = NULL,
  cat = TRUE
)
```

## Arguments

- results:

  Data frame of simulation results (same as for
  [`build_estimation_table`](https://larry-leon.github.io/forestsearch/reference/build_estimation_table.md)).

- dgm:

  DGM object with true parameter values.

- analysis_method:

  Character. Which analysis method to summarise. Default: `"FSlg"`.

- n_sims:

  Integer. Total number of simulations (for detection rate). If `NULL`
  (default), derived from `nrow(results)` after filtering to the
  analysis method.

- n_boots:

  Integer. Number of bootstraps (for narrative). Default: 300.

- digits:

  Integer. Decimal places for reported values. Default: 2.

- scenario:

  Character. One of `"null"` or `"alt"` (default). Controls the
  interpretive framing:

  - `"null"`: emphasises false-positive rate and selection bias

  - `"alt"`: emphasises power, bias relative to true effect

  If `NULL`, inferred from the DGM (`hr_H_true == hr_Hc_true` implies
  null).

- cat:

  Logical. If `TRUE` (default), prints the paragraph via
  [`cat()`](https://rdrr.io/r/base/cat.html). If `FALSE`, returns it
  invisibly as a character string (useful for programmatic insertion
  into Rmd via `results = "asis"`).

## Value

Invisibly returns the interpretation as a character string.

## See also

[`build_estimation_table`](https://larry-leon.github.io/forestsearch/reference/build_estimation_table.md),
[`format_oc_results`](https://larry-leon.github.io/forestsearch/reference/format_oc_results.md),
[`get_dgm_hr`](https://larry-leon.github.io/forestsearch/reference/get_dgm_hr.md)
