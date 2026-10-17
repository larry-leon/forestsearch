# Get L_eff at a Given Sample Size

Convenience function: loads calibration and computes
\\L\_{\text{eff}}(N) = C \cdot (N / n\_{\min})^\alpha\\.

## Usage

``` r
get_leff(N, outcome_type, output_dir = NULL)
```

## Arguments

- N:

  Integer. Sample size.

- outcome_type:

  Character. Outcome type to load.

- output_dir:

  Character. Default: NULL (auto).

## Value

Numeric L_eff value.
