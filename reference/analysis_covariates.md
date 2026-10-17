# Covariate set handed to forestsearch()/mrct_region_sims() for the paper's observed vs unobserved analysis cases.

Covariate set handed to forestsearch()/mrct_region_sims() for the
paper's observed vs unobserved analysis cases.

## Usage

``` r
analysis_covariates(
  dgm,
  case = c("observed", "unobserved"),
  drop = character(0)
)
```

## Arguments

- dgm:

  An MRCT data-generating mechanism returned by
  [`inject_mrct_structure()`](https://larry-leon.github.io/forestsearch/reference/inject_mrct_structure.md).
  Its `df_super` supplies the `z_` covariate columns and its `mrct`
  component supplies the region indicator name and the effect-modifier
  (`X2`) variables.

- case:

  Character. `"observed"` (default) keeps every `z_` covariate except
  the region indicator and spline helper columns; `"unobserved"`
  additionally drops the effect-modifier (`X2`) covariates.

- drop:

  Character vector of covariate names (without the `z_` prefix) to
  remove from the returned set. Default `character(0)` removes nothing.

## Value

Character vector of `z_`-prefixed covariate column names.
