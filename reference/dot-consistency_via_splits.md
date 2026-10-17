# Consistency proportion via literal splitting (single-stage)

Runs `n.splits` Bernoulli splits through
[`run_single_consistency_split()`](https://larry-leon.github.io/forestsearch/reference/run_single_consistency_split.md)
and returns the rounded proportion of consistent splits, or `NA_real_`
when too few valid splits are obtained. Survival and GLM
(`estimator_fn`) paths.

## Usage

``` r
.consistency_via_splits(
  df.x,
  N.x,
  n.splits,
  hr.consistency,
  cox_init,
  estimator_fn,
  consistency_threshold,
  adjust_covariates,
  pconsistency.digits,
  m = NA,
  details = FALSE
)
```
