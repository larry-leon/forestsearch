# Make an Effect-Estimator Closure

Returns a closure `function(data_slice)` that estimates the
within-subgroup treatment effect for the requested outcome type and
effect measure. The closure is called repeatedly inside the
splitting-consistency loop and the bootstrap worker, so it must be
lightweight and handle small-sample edge cases gracefully.

## Usage

``` r
make_effect_estimator(
  outcome_type,
  treat.name,
  outcome.name,
  event.name = NULL,
  offset.name = NULL,
  effect_measure = NULL,
  adverse_outcome = TRUE,
  robust_se = TRUE,
  adjust_covariates = NULL,
  ps_adjust_method = c("none", "iptw", "dr_gcomp"),
  ...
)
```

## Arguments

- outcome_type:

  Character. One of `"survival"`, `"binary"`, `"continuous"`, or
  `"count"`.

- treat.name:

  Character. Name of the binary treatment column (0/1).

- outcome.name:

  Character. Name of the outcome column.

- event.name:

  Character or `NULL`. Name of the event indicator (survival only).

- offset.name:

  Character or `NULL`. Name of the follow-up time column for rate-based
  measures (IRR, IRD).

- effect_measure:

  Character or `NULL`. For binary: `"OR"` (default), `"RR"`, `"RD"`,
  `"IRR"`, `"IRD"`. For continuous: `"MD"` (default). For survival:
  `"HR"` (default). `NULL` uses the default for the outcome type, which
  for binary is `"OR"` – the same default
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  uses.

- adverse_outcome:

  Logical. If `TRUE` (default), higher outcome values indicate harm and
  the effect is computed on the raw outcome. If `FALSE`, the outcome is
  beneficial (e.g., improvement) and the effect is computed on `1 - Y`
  internally so that effect \> 0 (or OR \> 1) consistently indicates
  treatment harm. Only applies to binary outcomes; ignored for survival,
  continuous, and count.

- robust_se:

  Logical. Use sandwich robust SE when available. Default `TRUE`.

- adjust_covariates:

  This option is implemented but has not been tested against the
  manuscript analyses or the simulation studies of León and Anderson
  (2026a, 2026b), which use the unadjusted treatment-only fits.
  Character vector or `NULL`. Additional covariate names to include in
  the outcome model formula.

- ps_adjust_method:

  This option is implemented but has not been tested against the
  manuscript analyses or the simulation studies of León and Anderson
  (2026a, 2026b), which use the unadjusted treatment-only fits.
  Character. Propensity score adjustment method. One of `"none"`
  (default), `"iptw"` (stabilized inverse probability of treatment
  weights), or `"dr_gcomp"` (Bang & Robins 2005 IPS covariate with
  G-computation).

- ...:

  Additional arguments (reserved for future extensions).

## Value

A closure `function(data_slice)` returning a list with components:
`estimate`, `se`, `converged`, `n0`, `n1`, `measure`, `method_used`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Binary outcome with odds ratio
fn <- make_effect_estimator("binary", "treat", "event",
  effect_measure = "OR")
result <- fn(my_data)
exp(result$estimate)  # odds ratio

# Poisson rate with IPTW adjustment
fn_iptw <- make_effect_estimator("binary", "treat", "event",
  offset.name = "time", effect_measure = "IRR",
  ps_adjust_method = "iptw")
} # }
```
