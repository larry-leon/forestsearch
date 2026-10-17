# Cache and validate cut expressions efficiently

Evaluates all cut expressions once and caches results to avoid redundant
evaluation. Much faster than evaluating repeatedly.

## Usage

``` r
evaluate_cuts_once(confs, df, details = FALSE)
```

## Arguments

- confs:

  Character vector of cut expressions.

- df:

  Data frame to evaluate expressions against.

- details:

  Logical. Print details during execution.

## Value

List with:

- evaluations: List of evaluated vectors (logical TRUE/FALSE) for each
  cut

- is_valid: Logical vector indicating which cuts produced \>1 unique
  value

- has_error: Logical vector indicating which cuts failed to evaluate

## Details

This replaces multiple eval(parse()) calls scattered throughout
get_FSdata. By caching results, we avoid:

1.  Repeated parsing of expressions

2.  Repeated evaluation on dataframe

3.  Redundant uniqueness checks

## Examples

``` r
df <- data.frame(age = c(40, 55, 70), size = c(20, 30, 25))
confs <- c("age <= 50", "size > 25")
evaluate_cuts_once(confs, df)
#> $evaluations
#> $evaluations[[1]]
#> [1]  TRUE FALSE FALSE
#> 
#> $evaluations[[2]]
#> [1] FALSE  TRUE FALSE
#> 
#> 
#> $is_valid
#> [1] TRUE TRUE
#> 
#> $has_error
#> [1] FALSE FALSE
#> 
```
