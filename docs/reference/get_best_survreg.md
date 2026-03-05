# Get Best Model from Comparison

Extracts the best fitting model object from a comparison result. If no
single best model can be determined, returns the Weibull model if
selected by either AIC or BIC. Defaults to Weibull0 model if no model
can be determined.

## Usage

``` r
get_best_survreg(comparison_result)
```

## Arguments

- comparison_result:

  Output from compare_survreg_models or compare_multiple_survreg

## Value

A survreg model object (defaults to Weibull0 model)
