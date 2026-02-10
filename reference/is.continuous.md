# Check if a variable is continuous

Determines if a variable is continuous based on the number of unique
values.

## Usage

``` r
is.continuous(x, cutoff = 4)
```

## Arguments

- x:

  A vector.

- cutoff:

  Integer. Minimum number of unique values to be considered continuous.

## Value

1 if continuous, 2 if not.
