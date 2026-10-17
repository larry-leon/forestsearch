# Get J-quantile cut expressions for variables

For each named variable in `conf_jcuts`, returns `J` deferred cut
expressions of the form `"X <= qj(X, k, J + 1)"` via
[`cut_var_jq()`](https://larry-leon.github.io/forestsearch/reference/cut_var_jq.md).
Variables that are not continuous (per `cont.cutoff`) are skipped with a
warning. This is the `cut_var_jq`/J-quantile analog of
[`get_conf_force()`](https://larry-leon.github.io/forestsearch/reference/get_conf_force.md).

## Usage

``` r
get_conf_force_jq(df, conf_jcuts, cont.cutoff = 4)
```

## Arguments

- df:

  Data frame.

- conf_jcuts:

  Named list of positive integers. Each name is a continuous-variable
  column in `df`; each value `J` is the number of binary cut expressions
  to emit for that variable. The cuts are placed at the (k/(J+1))-th
  empirical quantiles for `k = 1, ..., J`, partitioning the variable's
  range into `J + 1` non-overlapping intervals.

- cont.cutoff:

  Integer. Cutoff for continuous determination (passed through to
  [`is.continuous()`](https://larry-leon.github.io/forestsearch/reference/is.continuous.md)).

## Value

Character vector of cut expressions (length `sum(J_v)` over continuous
`v`), or `character(0)` if none.

## Details

Conflicts with other forced-cut mechanisms (`defaultcut_names`,
`conf.cont_medians_force`) are validated upstream by
[`get_FSdata()`](https://larry-leon.github.io/forestsearch/reference/get_FSdata.md);
this helper performs only structural checks.

## Examples

``` r
df <- data.frame(age = c(45, 60, 35, 50, 70, 25, 80, 55))
# 5 binary cuts at the 1/6, 2/6, ..., 5/6 quantiles of age:
get_conf_force_jq(df, conf_jcuts = list(age = 5))
#> [1] "age <= qj(age, 1, 6)" "age <= qj(age, 2, 6)" "age <= qj(age, 3, 6)"
#> [4] "age <= qj(age, 4, 6)" "age <= qj(age, 5, 6)"
```
