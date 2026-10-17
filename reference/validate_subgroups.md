# Validate subgroup definitions against a data frame

Evaluates every subgroup `id` expression against `data` (augmented with
`flag_itt = 1L`, the scalar `cutpoints` as columns, and 0-valued
benchmark membership columns) and fails with a named, collected error if
any expression errors, is non-logical, or has the wrong length. Run
automatically at the top of
[`run_subgroup_sims()`](https://larry-leon.github.io/forestsearch/reference/run_subgroup_sims.md)
so that a typo in a subgroup definition stops the study in milliseconds
instead of surfacing as a silent all-`NA` column after the simulation
loop.

## Usage

``` r
validate_subgroups(
  subgroups,
  data,
  cutpoints = list(),
  benchmarks = benchmark_spec()
)
```

## Arguments

- subgroups:

  List of subgroup definitions; each element a list with character
  scalars `id` (an R expression over columns), `name`, and `grp`.

- data:

  A data frame carrying the covariates the `id` expressions reference
  (e.g. `dgm$df_source` or `dgm$df_super`).

- cutpoints:

  Named list of scalars exposed as columns (e.g. `list(age_med = 53)`),
  matching the runner's `cutpoints` argument.

- benchmarks:

  A
  [`benchmark_spec()`](https://larry-leon.github.io/forestsearch/reference/benchmark_spec.md)
  (or `NULL`) whose membership columns are stubbed as `0L` during
  validation.

## Value

Invisibly `TRUE`; errors otherwise.
