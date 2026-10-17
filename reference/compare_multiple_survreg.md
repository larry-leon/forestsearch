# Compare Multiple Survival Regression Models

Performs comprehensive comparison of multiple survreg models including
convergence checking, information criteria comparison, and model
selection.

## Usage

``` r
compare_multiple_survreg(
  ...,
  model_names = NULL,
  verbose = TRUE,
  criteria = c("AIC", "BIC")
)
```

## Arguments

- ...:

  survreg model objects to compare

- model_names:

  Optional character vector of model names

- verbose:

  Logical, whether to print detailed output (default: TRUE)

- criteria:

  Character vector of criteria to use ("AIC", "BIC", or both)

## Value

A list of class "multi_survreg_comparison" containing:

- models:

  Named list of input models

- convergence:

  Convergence status for each model

- comparison:

  Model comparison statistics

- rankings:

  Model rankings by different criteria

- best_model:

  Name of the best model

- recommendation:

  Text recommendation

## Examples

``` r
if (FALSE) { # \dontrun{
fit1 <- survreg(Surv(time, status) ~ x, dist = "weibull")
fit2 <- survreg(Surv(time, status) ~ x, dist = "lognormal")
comparison <- compare_multiple_survreg(fit1, fit2)
} # }
```
