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
