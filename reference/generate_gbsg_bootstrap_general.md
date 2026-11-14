# Generate Synthetic GBSG Data using Generalized Bootstrap

Calls
[`generate_bootstrap_synthetic()`](https://larry-leon.github.io/forestsearch/reference/generate_bootstrap_synthetic.md)
on [`survival::gbsg`](https://rdrr.io/pkg/survival/man/gbsg.html) with
fixed variable lists. `grade` is passed in `ordinal_vars` only, so it is
resampled but not perturbed.

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

## Superseded

This GBSG-specific wrapper is kept for existing code. New work should
call
[`generate_bootstrap_synthetic()`](https://larry-leon.github.io/forestsearch/reference/generate_bootstrap_synthetic.md)
directly: the call in its first example, with the same `n`, `seed` and
`noise_level`, returns the same data. This wrapper will be deprecated in
a future release.

## Examples

``` r
if (FALSE) { # \dontrun{
df_synth <- generate_gbsg_bootstrap_general(n = 200, seed = 42)
dim(df_synth)
names(df_synth)
} # }
```
