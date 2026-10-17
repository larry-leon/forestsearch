# Compute CDE Analogue from GLM Potential Outcomes

Computes the Controlled Direct Effect analogue for GLM outcomes. For
binary (logit link), this averages the individual-level odds first, then
takes the ratio — analogous to `mean(exp(theta_1)) / mean(exp(theta_0))`
in survival:

## Usage

``` r
compute_cde_glm(df, subset_indicator = NULL, effect_measure = "OR")
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

Numeric CDE value (binary OR only), or `NA_real_`.

## Details

`CDE(S) = mean(p1[S]/(1-p1[S])) / mean(p0[S]/(1-p0[S]))`

This differs from the AOR (`odds(mean(p1)) / odds(mean(p0))`) due to
Jensen's inequality whenever individual probabilities are heterogeneous.
For identity-link (MD) and log-link (IRR) outcomes, CDE equals AOR and
this function returns `NA_real_`.
