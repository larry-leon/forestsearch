# GLM per-subgroup fit function

Returns a function `data -> c(estimate, UB)` for non-survival outcomes
by delegating estimation to the package's native effect-estimator layer:
the closure from
[`make_effect_estimator()`](https://larry-leon.github.io/forestsearch/reference/make_effect_estimator.md)
(the same code the
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
search, bootstrap, and consistency loops use), including the
bit-identical unadjusted fast path (`.make_lm_estimator_fast()`) under
the same routing predicate
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
applies. The estimator's internal-scale `estimate`/`se` are mapped to
the wrapper contract as
`upper = estimate + qnorm(1 - (1 - level)/2) * se`, back-transformed
([`exp()`](https://rdrr.io/r/base/Log.html)) for ratio measures – the
package's own display convention – so ratio-measure matrices land on the
natural scale exactly as
[`subgroup_cox()`](https://larry-leon.github.io/forestsearch/reference/subgroup_cox.md)'s
`exp(coef)` / `upper .95` do. Estimation failure, non-convergence, `NA`
estimate/SE (e.g. a single-arm subgroup, where the treatment coefficient
is aliased), or a non-finite result (e.g. a zero-residual-df fit)
returns `c(NA, NA)`.

## Usage

``` r
subgroup_glm(
  outcome.name = "y_sim",
  treat.name = "treat_sim",
  outcome_type = "continuous",
  effect_measure = NULL,
  adverse_outcome = TRUE,
  robust_se = TRUE,
  offset.name = NULL,
  adjust_covariates = NULL,
  level = 0.95,
  lean = TRUE
)
```

## Arguments

- outcome.name, treat.name:

  Column names of the outcome and the 0/1 treatment indicator; defaults
  match
  [`simulate_from_glm_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_glm_dgm.md)'s
  `y_sim` / `treat_sim`.

- outcome_type:

  `"continuous"`, `"binary"`, or `"count"` (all validated; count
  requires `offset.name`).

- effect_measure:

  `NULL` resolves to the wrapper default for the outcome type: `"MD"`
  (continuous), `"OR"` (binary – the
  [`generate_glm_dgm()`](https://larry-leon.github.io/forestsearch/reference/generate_glm_dgm.md)
  default and the calibration target; note that
  [`make_effect_estimator()`](https://larry-leon.github.io/forestsearch/reference/make_effect_estimator.md)'s
  own binary default is `"RD"`), or `"IRR"` (count – all layers agree).
  Binary also accepts `"RD"`; count also accepts `"IRD"`; see the
  supported-settings section for the gated measures.

- adverse_outcome:

  Higher outcome = harm (`TRUE`, default). See the direction convention
  above.

- robust_se, offset.name, adjust_covariates:

  Passed to
  [`make_effect_estimator()`](https://larry-leon.github.io/forestsearch/reference/make_effect_estimator.md).
  `robust_se` is consulted only where the estimator layer uses
  sandwich-style SEs (the count IRR path, the RR modified-Poisson
  fallback, and other rate estimators); the unadjusted OR path, the RD
  tiers, the unadjusted IRD delta chain, and MD use model-based/analytic
  SEs, so the flag is a no-op for those measures (asserted in the phase
  gates). `offset.name` is REQUIRED for `outcome_type = "count"`: it
  names a strictly positive exposure / follow-up time column in the
  analysis data (the estimator adds an internal `.log_offset` column,
  clamping non-positive values at `.Machine$double.eps`); for plain
  counts supply a unit column of 1s – `log(1) = 0` reproduces the
  no-offset Poisson model exactly. Supplying `adjust_covariates`
  disables the unadjusted fast path (matching
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s
  routing).

- level:

  Confidence level for the upper bound (default `0.95`).

- lean:

  Re-home the provenance formula's environment (default `TRUE`).

## Value

A function of one argument (`data`) returning `c(estimate, ub)`; carries
the provenance formula in `attr(, "formula")` and scale metadata in
`attr(, "effect")` – a list with `measure`, `outcome_type`, `log_scale`,
`null_value`, `adverse_outcome`, `est_label`, `ub_label`,
`est_thresholds`, `ub_thresholds`, `level` – which
[`run_subgroup_sims()`](https://larry-leon.github.io/forestsearch/reference/run_subgroup_sims.md)
stamps onto GLM results for scale-aware summaries and plots. A custom
fitter can opt into the same behavior by attaching an `"effect"`
attribute of this shape.

## Direction convention

Estimates follow the estimator layer's harm-positive normalization: with
`adverse_outcome = FALSE` (higher outcome = better) the outcome is
negated internally so that `estimate > 0` consistently indicates
treatment harm, aligning with the survival convention where HR \> 1 =
harm. Consequently the null value is always 0 (identity measures) or 1
(ratio measures), whatever the raw outcome's direction.

## Supported settings

Phase 4.5 validates `outcome_type = "continuous"` (MD),
`outcome_type = "binary"` (`"OR"` – the binary default – and `"RD"`),
and `outcome_type = "count"` (`"IRR"` – the count default – and `"IRD"`;
`offset.name` is required, see below). Binary `"RR"`
(estimator-supported; wrapper parity pending) still stops with an
informative error; the estimator layer
([`make_effect_estimator()`](https://larry-leon.github.io/forestsearch/reference/make_effect_estimator.md))
already supports it and the gate lifts with its parity phase – no
interface change.

## Serialization (`lean`)

The returned closure is stripped of source references and given a
minimal environment (the estimator closure and two scalars), and the
estimator closure itself is rebuilt onto a fresh environment holding
only its forced argument values (parented on the package namespace,
which serializes by reference) with its own source references stripped.
This severs shared srcfile environments (keep-source installs) and the
constructor's frame chain from the serialized graph, so the fitter ships
to parallel workers in the low single-digit KB regardless of session
state. The provenance formula attached in `attr(, "formula")` (built
from the column names, e.g. `y_sim ~ treat_sim`) is re-homed into an
empty stats-parented environment when `lean = TRUE`, keeping the
parallel-shipment guard in
[`run_subgroup_sims()`](https://larry-leon.github.io/forestsearch/reference/run_subgroup_sims.md)
meaningful.

## See also

[`subgroup_cox()`](https://larry-leon.github.io/forestsearch/reference/subgroup_cox.md),
[`run_subgroup_sims()`](https://larry-leon.github.io/forestsearch/reference/run_subgroup_sims.md),
[`make_effect_estimator()`](https://larry-leon.github.io/forestsearch/reference/make_effect_estimator.md),
[`generate_glm_dgm()`](https://larry-leon.github.io/forestsearch/reference/generate_glm_dgm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- subgroup_glm()                      # y_sim ~ treat_sim, MD
fit(simulate_from_glm_dgm(dgm, n = 500, seed = 1))

# Binary (Phase 4.4): OR by default, RD available
fit_or <- subgroup_glm(outcome_type = "binary")
fit_rd <- subgroup_glm(outcome_type = "binary", effect_measure = "RD")

# Count (Phase 4.5): IRR by default, IRD available; offset required
fit_irr <- subgroup_glm(outcome_type = "count", offset.name = "t_exp")
fit_ird <- subgroup_glm(outcome_type = "count", offset.name = "t_exp",
                        effect_measure = "IRD")
} # }
```
