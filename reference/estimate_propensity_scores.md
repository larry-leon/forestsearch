# Estimate Propensity Scores

Computes \\P(W = 1 \mid X)\\ using the requested method, then constructs
both stabilized IPTW weights and the Bang & Robins (2005) inverse
propensity score (IPS) covariate. The data frame is returned with three
new columns: `ps_hat` (raw PS), `sw` (stabilized IPTW weights), and
`ips_covar` (IPS covariate).

## Usage

``` r
estimate_propensity_scores(
  data,
  treat.name,
  confounders.name,
  method = c("grf", "lasso", "logistic", "none"),
  grf_forest = NULL,
  seed = 8316951L,
  trim = c(0.025, 0.975)
)
```

## Arguments

- data:

  Data frame.

- treat.name:

  Character. Name of binary treatment column (0/1).

- confounders.name:

  Character vector. Covariate names for the propensity model.

- method:

  Character. One of `"grf"` (default), `"lasso"`, `"logistic"`, or
  `"none"`.

- grf_forest:

  Optional. A fitted `causal_forest` object from which to extract
  `W.hat`.

- seed:

  Integer. Random seed for GRF.

- trim:

  Numeric vector of length 2. Quantile bounds for PS trimming. Default:
  `c(0.025, 0.975)`.

## Value

A list with:

- data:

  Data frame with `ps_hat`, `sw`, and `ips_covar` columns appended.

- ps_hat:

  Numeric vector of estimated propensity scores.

- sw:

  Numeric vector of stabilized IPTW weights.

- ips_covar:

  Numeric vector of inverse PS covariates (Bang & Robins 2005).

- method:

  Character. Method actually used.

- trimmed:

  Integer. Number of observations trimmed.

## References

Bang H, Robins JM (2005). "Doubly Robust Estimation in Missing Data and
Causal Inference Models." *Biometrics* 61(4): 962–973.

Hernan MA, Robins JM, Brumback B (2000). "Marginal Structural Models and
Causal Inference in Epidemiology." *Epidemiology* 11(5): 550–560.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- survival::gbsg
df$grade3 <- as.integer(df$grade == "3")
ps_result <- estimate_propensity_scores(
  data = df,
  treat.name = "hormon",
  confounders.name = c("age", "meno", "size", "grade3", "nodes",
                       "pgr", "er"),
  method = "grf"
)
summary(ps_result$ps_hat)
summary(ps_result$sw)
} # }
```
