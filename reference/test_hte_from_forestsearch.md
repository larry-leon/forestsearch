# HTE Test Using ForestSearch-Selected Covariates

HTE Test Using ForestSearch-Selected Covariates

## Usage

``` r
test_hte_from_forestsearch(
  fs_result,
  df,
  outcome.name,
  treat.name,
  regression = c("ols", "logistic", "poisson"),
  offset = NULL,
  poly_order = 1L
)
```

## Arguments

- fs_result:

  A `forestsearch` result object.

- df:

  Data frame used in the
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  call.

- outcome.name:

  Character. Outcome variable name.

- treat.name:

  Character. Treatment variable name.

- regression:

  Character: "ols", "logistic", or "poisson".

- offset:

  Numeric offset vector (Poisson only).

- poly_order:

  Integer. Default: 1.

## Value

List with `hte_test` (the
[`test_hte`](https://larry-leon.github.io/forestsearch/reference/test_hte.md)
result), `covariates_used` (the `df` columns entered in the test),
`n_covariates` (their number), `n_candidates`
(`length(fs_result$find.grps$max_sg_est)`, `NA` when that is absent) and
`source`; or `NULL` with a warning when no covariate is found.

The covariates come from `fs_result$confounders.evaluated` when it is
non-empty, else from `fs_result$confounders.candidate`, taken as column
names. Each `confounders.evaluated` label is mapped to columns of `df`:
a label that is exactly a column name (a categorical confounder such as
`"meno"`, including a non-syntactic name such as `"ER-pos"`) is that
column; any other label (a cut such as `"age <= 53"`) is parsed, and its
variables that are columns of `df` are used. Each column is used once,
in the order first met. Two warnings: a label that resolves to no column
of `df` warns that it "cannot be resolved to a column of df" and is
skipped; and when no covariate column is left, the function warns "No
covariate found in df for the ForestSearch covariates" and returns
`NULL`. (When neither field is present it warns "No covariates found."
and returns `NULL`.) Factor and character columns enter the test as
their level codes minus 1, other columns as numbers.

`source` names the field used:
`"confounders.evaluated (GRF/LASSO selected)"` or
`"confounders.candidate (all candidates)"`. The first string is fixed.
On the default (consistency) route `confounders.evaluated` holds the
fit's candidate factor labels, which are LASSO-screened or include GRF
cuts only when the fit used `use_lasso` or `use_grf` (both `FALSE` by
default); on the `subgroup_method = "dina"` and `"grf"` routes it is the
fit's `confounders.name`.
