# Internal Workhorse for Single-Theta Detection Probability

Internal Workhorse for Single-Theta Detection Probability

## Usage

``` r
.detect_prob_glm_single(beta, sigma2_s, k1, k2, method, n_mc, seed)
```

## Arguments

- beta:

  Numeric. True effect on natural parameter scale.

- sigma2_s:

  Numeric. Split-half variance (8 / d_eff).

- k1:

  Numeric. Screening threshold (natural scale).

- k2:

  Numeric. Consistency threshold (natural scale).

- method:

  Character. Integration method.

- n_mc:

  Integer. Monte Carlo samples.

- seed:

  Integer. RNG seed.

## Value

Numeric probability.
