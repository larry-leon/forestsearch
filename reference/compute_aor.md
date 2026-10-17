# Compute Marginal Causal Effect from GLM Potential Outcomes

Computes the marginal (average) treatment effect from individual-level
potential outcomes stored in GLM simulation data. This is the GLM
analogue of
[`compute_ahr`](https://larry-leon.github.io/forestsearch/reference/compute_ahr.md)
for survival data.

## Usage

``` r
compute_aor(df, subset_indicator = NULL, effect_measure = "OR")
```

## Arguments

- df:

  Data frame with potential outcome columns (`p0`/`p1` or `mu0`/`mu1`).

- subset_indicator:

  Optional logical/integer vector for subsetting. If provided, only rows
  where `subset_indicator == 1` are used.

- effect_measure:

  Character. One of `"OR"`, `"RD"`, `"RR"`, `"MD"`, `"IRR"`, `"IRD"`.
  Default `"OR"`.

## Value

Numeric effect value, or `NA_real_` if required columns are missing.

## Details

For binary outcomes (`p0`, `p1`):

- `OR`: `odds(mean(p1)) / odds(mean(p0))`

- `RD`: `mean(p1) - mean(p0)`

- `RR`: `mean(p1) / mean(p0)`

For continuous or count outcomes (`mu0`, `mu1`):

- `MD`: `mean(mu1) - mean(mu0)`

- `IRR`: `mean(mu1) / mean(mu0)`

- `IRD`: `mean(mu1) - mean(mu0)`
