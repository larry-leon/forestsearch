# Validate input data for GRF analysis

Checks that input data meets requirements for GRF analysis

## Usage

``` r
validate_grf_data(W, D, n.min)
```

## Arguments

- W:

  Numeric vector. Treatment indicator

- D:

  Numeric vector. Event indicator

- n.min:

  Integer. Minimum subgroup size

## Value

Logical. TRUE if data is valid, FALSE with warning otherwise

## Examples

``` r
W <- rep(0:1, each = 50)
D <- rbinom(100, 1, 0.6)
validate_grf_data(W, D, n.min = 60)
#> Warning: Insufficient sample size: treatment=50, control=50, required=60 per arm.
#> [1] FALSE
```
