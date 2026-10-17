# Effective Information for Count Outcomes (Poisson)

Computes \\d\_{\text{eff}} = D / \phi\\ where \\D\\ is the total
expected events and \\\phi\\ is the dispersion parameter.

## Usage

``` r
d_eff_count(
  n_sg = NULL,
  event_rate = NULL,
  total_events = NULL,
  dispersion = 1
)
```

## Arguments

- n_sg:

  Integer. Subgroup sample size. Required if `total_events` is not
  provided.

- event_rate:

  Numeric. Mean events per patient. Required if `total_events` is not
  provided.

- total_events:

  Numeric. Total observed events \\D\\. Overrides `n_sg * event_rate` if
  provided.

- dispersion:

  Numeric. Overdispersion parameter (1.0 for Poisson, greater than 1 for
  quasi-Poisson). Default: 1.0.

## Value

Numeric. Effective information.

## Details

For Poisson regression with log-offset: \$\$\text{Var}(\hat\beta) = \phi
\left(\frac{1}{D_0} + \frac{1}{D_1}\right) \approx \frac{4\phi}{D}\$\$
Under proportional hazards, \\D = d\\ (Cox events), so
\\d\_{\text{eff}}^{\text{Poisson}} = d\_{\text{eff}}^{\text{Cox}}\\
(Laird and Olivier, 1981).

## Examples

``` r
# From total events
d_eff_count(total_events = 55)
#> [1] 55

# From sample size and rate
d_eff_count(n_sg = 100, event_rate = 0.8)
#> [1] 80

# Quasi-Poisson with overdispersion
d_eff_count(total_events = 55, dispersion = 1.5)
#> [1] 36.66667
```
