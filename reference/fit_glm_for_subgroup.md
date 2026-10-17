# Fit GLM for Subgroup via Estimator Closure

Subsets the analysis data frame to the candidate subgroup and calls the
pre-built estimator closure. Returns a list compatible with
[`create_result_row`](https://larry-leon.github.io/forestsearch/reference/create_result_row.md)
so the downstream pipeline is unchanged.

## Usage

``` r
fit_glm_for_subgroup(df_clean, id.x, estimator_fn)
```

## Arguments

- df_clean:

  Data frame. The analysis data, row-aligned with the cleaned vectors.

- id.x:

  Integer vector. 1 = subject in this subgroup, 0 = not.

- estimator_fn:

  Closure from
  [`make_effect_estimator`](https://larry-leon.github.io/forestsearch/reference/make_effect_estimator.md).

## Value

A list with components `hr`, `lower`, `upper`, `med0`, `med1`, or `NULL`
on failure.

## Details

For binary outcomes with `effect_measure = "RD"`, the `hr` slot contains
the risk difference (not a hazard ratio) – the name is retained for
pipeline compatibility. Similarly, `lower`/`upper` are Wald-type CI
bounds and `med0`/`med1` are `NA`.
