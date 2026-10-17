# Per-candidate influence pieces (Cox / GLM dispatch)

Thin dispatch onto the existing forestsearch internals so MR uses the
identical treatment `dfbeta`, `beta_hat`, and robust `sigma_D` the
resample consistency engine uses.

## Usage

``` r
.fs_mr_pieces(df_sub, spec)
```

## Arguments

- df_sub:

  Data frame for one candidate subgroup's members.

- spec:

  List with `outcome_type`, `effect_measure`, `treat.name`,
  `outcome.name`, `event.name`, `offset.name`, `adjust_covariates`,
  `adverse_outcome`.

## Value

Output of
[`.consistency_glm_pieces()`](https://larry-leon.github.io/forestsearch/reference/dot-consistency_glm_pieces.md)
/
[`.consistency_cox_pieces()`](https://larry-leon.github.io/forestsearch/reference/dot-consistency_cox_pieces.md),
or `NULL`.
