# Find Covariate Any Match

Helper function to determine if any CV fold found a subgroup involving
the same covariate (not necessarily same cut).

## Usage

``` r
find_covariate_any_match(sg_target, sg1, sg2, confs)
```

## Arguments

- sg_target:

  Character. Target subgroup definition to match.

- sg1:

  Character vector. Subgroup 1 labels for each fold.

- sg2:

  Character vector. Subgroup 2 labels for each fold.

- confs:

  Character vector. Confounder names.

## Value

Numeric vector (0/1) indicating match for each fold.
