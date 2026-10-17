# Effective Information for Survival (Cox PH)

Computes \\d\_{\text{eff}} = n_H (1 - p_c)\\ where \\p_c\\ is the
censoring proportion.

## Usage

``` r
d_eff_survival(n_sg, prop_cens)
```

## Arguments

- n_sg:

  Integer. Subgroup sample size.

- prop_cens:

  Numeric. Proportion censored (0–1).

## Value

Numeric. Effective information (expected events).

## Examples

``` r
d_eff_survival(n_sg = 60, prop_cens = 0.45)
#> [1] 33
# 33
```
