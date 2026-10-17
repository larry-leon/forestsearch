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
  cat = TRUE,
  subgroup_notation = c("harm", "benefit"),
  trim_threshold = 1000,
  trim_fraction = 0.01
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

- subgroup_notation:

  Character. `"harm"` (default) labels subgroups as H/Hc; `"benefit"`
  labels as G/Gc for benefit-search analyses (treatment switching).

- trim_threshold:

  Numeric or `NULL`. When the raw mean of an estimate vector exceeds
  this absolute value, the narrative reports averages and SDs computed
  on the central `1 - 2 * trim_fraction` of estimates and adds a closing
  note disclosing the trimming. Set `NULL` to disable. Default: `1000`
  (matches `build_estimation_table`).

- trim_fraction:

  Numeric in `(0, 0.5)`. Fraction of estimates to trim from each tail
  when trimming triggers. Default: `0.01`. Has no effect when
  `trim_threshold = NULL`.

## Value

Invisibly returns the interpretation as a character string.

## See also

[`build_estimation_table`](https://larry-leon.github.io/forestsearch/reference/build_estimation_table.md),
[`format_oc_results`](https://larry-leon.github.io/forestsearch/reference/format_oc_results.md),
[`get_dgm_hr`](https://larry-leon.github.io/forestsearch/reference/get_dgm_hr.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# In a vignette chunk with results = "asis":
interpret_estimation_table(results_null, dgm_null, scenario = "null")

# Capture for further processing:
txt <- interpret_estimation_table(results_alt, dgm_alt, cat = FALSE)
} # }
```
