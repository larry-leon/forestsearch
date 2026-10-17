# Effective Information for Continuous Outcomes (Gaussian)

Computes \\d\_{\text{eff}} = n_H / \sigma_Y^2\\ from the Fisher
information of the Gaussian linear model treatment coefficient.

## Usage

``` r
d_eff_continuous(n_sg, sigma_y)
```

## Arguments

- n_sg:

  Integer. Subgroup sample size.

- sigma_y:

  Numeric. Residual standard deviation of the outcome.

## Value

Numeric. Effective information.

## Details

For the model \\Y_i = \beta_0 + \beta V_i + \epsilon_i\\, \\\epsilon_i
\sim N(0, \sigma^2)\\, with balanced arms: \$\$\text{Var}(\hat\beta) =
\frac{4\sigma^2}{n_H}\$\$

## Examples

``` r
d_eff_continuous(n_sg = 100, sigma_y = 1.5)
#> [1] 44.44444
# 44.4
```
