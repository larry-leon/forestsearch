# Simulate Trial Data from GBSG DGM

Generates simulated clinical trial data from a GBSG-based data
generating mechanism.

## Usage

``` r
simulate_from_gbsg_dgm(
  dgm,
  n = NULL,
  rand_ratio = 1,
  sim_id = 1,
  max_follow = Inf,
  muC_adj = 0,
  min_cens = NULL,
  max_cens = NULL,
  draw_treatment = TRUE
)
```

## Arguments

- dgm:

  A "gbsg_dgm" object from
  [`create_gbsg_dgm`](https://larry-leon.github.io/forestsearch/reference/create_gbsg_dgm.md)

- n:

  Integer. Sample size. If NULL, uses full super-population. Default:
  NULL

- rand_ratio:

  Numeric. Randomization ratio (treatment:control). Default: 1 (1:1
  randomization)

- sim_id:

  Integer. Simulation ID used for seed offset. Default: 1

- max_follow:

  Numeric. Administrative censoring time (months). Default: Inf (no
  administrative censoring)

- muC_adj:

  Numeric. Adjustment to censoring distribution location parameter.
  Positive values increase censoring. Default: 0

- min_cens:

  Numeric. Minimum censoring time for uniform censoring. Required if
  cens_type = "uniform"

- max_cens:

  Numeric. Maximum censoring time for uniform censoring. Required if
  cens_type = "uniform"

- draw_treatment:

  Logical. If TRUE, randomly assigns treatment. If FALSE, samples from
  existing treatment arms. Default: TRUE

## Value

Data frame with simulated trial data including:

- id:

  Subject identifier

- y.sim:

  Observed follow-up time

- event.sim:

  Event indicator (1 = event, 0 = censored)

- t.sim:

  True event time (before censoring)

- treat:

  Treatment indicator

- flag.harm:

  Harm subgroup indicator

- loghr_po:

  Individual log hazard ratio (potential outcome)

- v1-v7:

  Analysis factors
