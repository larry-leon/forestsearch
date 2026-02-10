# Calculate Hazard Ratios from Potential Outcomes

Calculate Hazard Ratios from Potential Outcomes

## Usage

``` r
calculate_hazard_ratios(df_super, n_super, mu, tau, model, verbose)
```

## Arguments

- df_super:

  Data frame with super population

- n_super:

  Size of super population

- mu:

  Intercept parameter

- tau:

  Scale parameter

- model:

  Model type ("alt" or "null")

- verbose:

  Logical for verbose output

## Value

List of hazard ratios
