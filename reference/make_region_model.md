# Build the logistic region model and draw Region on a population.

Build the logistic region model and draw Region on a population.

## Usage

``` r
make_region_model(
  pop,
  x_pred = NULL,
  or_pred = 1,
  x1_vars = character(0),
  or_x1 = 1,
  prevalence = 0.2,
  scale_ref = NULL,
  scale = c("rank", "minmax"),
  band = NULL
)
```

## Arguments

- pop:

  data.frame holding x_pred and x1_vars (natural scale)

- x_pred:

  name of the region-associated effect modifier (X2); NULL for no X2
  term (regional imbalance only through x1_vars)

- or_pred:

  odds ratio for s(x_pred) -\> Region (OR = 1: none)

- x1_vars:

  names of imbalanced-only covariates (X1); character(0) for none

- or_x1:

  odds ratios for the X1 terms (recycled)

- prevalence:

  target marginal P(Region = 1)

- scale_ref:

  data.frame used to define the scaling (defaults to pop)

- scale:

  "rank" (ECDF, default) or "minmax" (the paper's scaling)

- band:

  NULL (default, no band term), or list(var, cut, or): adds
  `log(or) * 1(df[[var]] <= cut)` to the logit (var on its natural
  scale, not passed through s())

## Value

list(alpha0, coefs, band, prob = function(df), draw = function(df,
seed))
