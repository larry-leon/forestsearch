# Build the Cox/GLM `spec` for MR from forestsearch()'s resolved columns.

DINA/GRF operate on the user-facing analysis frame `df` (not the
standardized `df.fs`), so survival uses the user column names here –
unlike the consistency hook, which uses the internal Y/Event/Treat names
on `df.fs`.

## Usage

``` r
.fs_mr_spec(
  outcome_type,
  effect_measure,
  treat.name,
  outcome.name,
  event.name,
  offset.name,
  adjust_covariates,
  adverse_outcome,
  df
)
```
