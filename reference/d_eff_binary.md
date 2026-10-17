# Effective Information for Binary Outcomes (Logistic)

Computes \\d\_{\text{eff}} = n_H \bar{p} (1 - \bar{p})\\ from the Fisher
information of the logistic regression treatment coefficient.

## Usage

``` r
d_eff_binary(n_sg, p_event)
```

## Arguments

- n_sg:

  Integer. Subgroup sample size.

- p_event:

  Numeric. Event probability (0–1), averaged across arms.

## Value

Numeric. Effective information.

## Details

For the model \\\text{logit}(p_i) = \beta_0 + \beta V_i\\ with balanced
arms (\\n_1 = n_0 = n_H/2\\): \$\$\text{Var}(\hat\beta) = \frac{4}{n_H
\bar{p}(1-\bar{p})}\$\$

## Examples

``` r
d_eff_binary(n_sg = 100, p_event = 0.30)
#> [1] 21
# 21
```
