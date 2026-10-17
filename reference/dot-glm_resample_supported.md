# Can the GLM resampling approximation represent this effect measure?

The resampling approximation requires the treatment effect to be a
single model coefficient with a well-defined influence function. That
holds for OR/RR (logistic / log-binomial), RD (identity-link binomial),
MD (OLS), and IRR (Poisson with offset). It does not hold for IRD (a
delta-method rate difference) or propensity-adjusted (IPTW /
G-computation) effects, for which the caller falls back to literal
splitting.

## Usage

``` r
.glm_resample_supported(spec)
```

## Arguments

- spec:

  List; the `glm_resample_spec` threaded from
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md).

## Value

Logical scalar.
