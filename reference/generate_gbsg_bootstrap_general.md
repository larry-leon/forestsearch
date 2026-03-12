# Generate Synthetic GBSG Data using Generalized Bootstrap

Generate Synthetic GBSG Data using Generalized Bootstrap

## Usage

``` r
generate_gbsg_bootstrap_general(n = 686, seed = 123, noise_level = 0.1)
```

## Arguments

- n:

  Number of observations

- seed:

  Random seed

- noise_level:

  Noise level for perturbation

## Value

Synthetic GBSG dataset

## Examples

``` r
if (FALSE) { # \dontrun{
df_synth <- generate_gbsg_bootstrap_general(n = 200, seed = 42)
dim(df_synth)
names(df_synth)
} # }
```
