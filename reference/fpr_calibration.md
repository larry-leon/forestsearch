# Direct Simulation-Based FPR Calibration at a Specific Threshold Pair

Estimates the procedure-level false-positive rate (FPR) for a
[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
analysis at a **specific** \\(c_1, c_2)\\ threshold pair by running
`nsims` null simulations. No power-law model or threshold-independence
assumption is required.

## Usage

``` r
fpr_calibration(
  c1,
  c2,
  df.analysis = NULL,
  fs_params = list(),
  nsims = 200L,
  nsims_verify = 0L,
  null_fn = NULL,
  n_sample = NULL,
  n_min = 60L,
  outcome_type_approx = NULL,
  d_eff_approx = NULL,
  prop_cens = 0,
  p_event = 0.5,
  sigma_y = 1,
  sim_seed = 42L,
  n_workers = 1L,
  quiet = FALSE
)
```

## Arguments

- c1:

  Numeric scalar. The **screening threshold** on a candidate subgroup's
  own effect. `c1` is this function's name for
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s
  `effect.threshold` (legacy `hr.threshold`), and is passed through to
  it under that name, on the natural scale –
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  resolves the comparison scale per estimand.

- c2:

  Numeric scalar. The **per-split effect threshold**. `c2` is this
  function's name for
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s
  `consistency.threshold` (legacy `hr.consistency`), passed through
  under that name on the natural scale. It is an effect threshold
  applied within each split half, **not** the consistency rate: the rate
  is `pconsistency.threshold` and travels in `fs_params`. Must satisfy
  `c2 <= c1`.

- df.analysis:

  Data frame. The analyst's dataset. Used directly in permutation mode;
  also used to extract N for the per-subgroup approximation. Can be
  `NULL` if `null_fn` is supplied and `n_sample` is given explicitly.

- fs_params:

  Named list. All forestsearch() arguments except `hr.threshold` and
  `hr.consistency` (which are set to `c1`/`c2`) and `df.analysis`
  (supplied separately). Must include at minimum `outcome.name`,
  `treat.name`, `id.name`, and — for binary/continuous/count —
  `outcome_type` and `effect_measure`.

- nsims:

  Integer. Number of null replicates. Default 200. At FPR ~ 10%,
  `nsims = 200` gives a 95% CI half-width of ~±2.1 pp; `nsims = 500`
  gives ~±1.3 pp.

- nsims_verify:

  Integer. Number of independent verification replicates. Uses a
  separate seed range from the calibration batch to provide an honest
  check of the corrected FPR. Default 0 (skip verification). Set to
  `nsims` or larger for production use.

- null_fn:

  Function or `NULL`. If provided, called as `null_fn(seed)` and must
  return a data frame with the same column structure as `df.analysis`.
  Overrides permutation mode. Default `NULL` (permutation mode).

- n_sample:

  Integer or `NULL`. Sample size used when `null_fn` is provided.
  Ignored in permutation mode (N is taken from `nrow(df.analysis)`).

- n_min:

  Integer. Minimum subgroup size for the per-subgroup approximation
  \\P_1\\. Should match `fs_params$n.min`. Default 60.

- outcome_type_approx:

  Character. One of `"binary"`, `"survival"`, `"count"`, `"continuous"`.
  Used only for computing the asymptotic \\P_1\\ approximation. If
  `NULL` (default), inferred from `fs_params$outcome_type`.

- d_eff_approx:

  Numeric or `NULL`. Effective information for \\P_1\\. If `NULL`
  (default), estimated from `df.analysis` (or `n_sample`) using
  `d_eff_*()`. Supply explicitly if you want to override (e.g. use an
  expected event count rather than observed).

- prop_cens:

  Numeric. Censoring proportion for survival outcomes. Used only when
  `d_eff_approx` is not supplied and `outcome_type_approx = "survival"`.
  Default 0.

- p_event:

  Numeric. Event probability for binary outcomes. Used only when
  `d_eff_approx` is not supplied and `outcome_type_approx = "binary"`.
  Default 0.5.

- sigma_y:

  Numeric. Residual SD for continuous outcomes. Default 1.

- sim_seed:

  Integer. Base random seed. Each replicate uses `sim_seed + i`. Default
  42.

- n_workers:

  Integer. Parallel workers for the foreach loop. Default 1
  (sequential). Set to `parallel::detectCores() - 1` for production
  runs.

- quiet:

  Logical. Suppress progress messages. Default `FALSE`.

## Value

An object of class `"fpr_calibration"` with components:

- `fpr`:

  Numeric. Estimated procedure-level FPR.

- `fpr_ci`:

  Numeric vector of length 2. 95% Wilson CI.

- `n_found`:

  Integer. Number of replicates where a subgroup was found.

- `nsims`:

  Integer. Total replicates attempted.

- `n_failed`:

  Integer. Replicates that errored.

- `P1`:

  Numeric. Per-subgroup approximation \\P_1\\ at `(c1, c2)`.

- `d_eff`:

  Numeric. Effective information used for \\P_1\\.

- `L_eff_implied`:

  Numeric. Implied \\L\_{\text{eff}}\\ from the simulation result:
  \\\log(1 - \text{FPR}) / \log(1 - P_1)\\.

- `c1, c2`:

  Numeric. The threshold pair.

- `N`:

  Integer. Sample size.

- `mode`:

  Character. `"permutation"` or `"dgm"`.

- `call_args`:

  List. Key call arguments for reproducibility.

## Details

Two null-generation modes are supported:

- `null_fn` provided:

  Each replicate draws a fresh null dataset from `null_fn(seed)`. Use
  this when you have a parametric data generating mechanism (DGM) or
  want to mimic a specific study design.

- `df.analysis` only (permutation mode):

  Each replicate permutes the treatment column within the analyst's
  actual dataset. This preserves the observed covariate and outcome
  distributions and is the natural H0 reference for a real analysis.
  Requires the `treat.name` column to be present in `df.analysis`.

### Permutation mode

Under the null hypothesis of no subgroup heterogeneity, the treatment
assignment is exchangeable with the observed outcomes. Permuting
`treat.name` within each replicate breaks any true treatment-by-
subgroup interaction while preserving the marginal distributions of all
other variables. This is a model-free H0 reference that adapts
automatically to the actual dataset structure.

### P1 approximation

The asymptotic per-subgroup \\P_1\\ at `(c1, c2)` is computed via
[`compute_detection_probability_glm()`](https://larry-leon.github.io/forestsearch/reference/compute_detection_probability_glm.md)
using the effective information \\d\_{\text{eff}}\\ at `n.min`. The
ratio \\L\_{\text{eff}} = \log(1 - \hat{\text{FPR}}) / \log(1 - P_1)\\
is the implied multiplicity factor specific to this threshold pair and
dataset — it replaces the universal power-law model and does not assume
threshold-independence.

## See also

[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md),
[`compute_detection_probability_glm`](https://larry-leon.github.io/forestsearch/reference/compute_detection_probability_glm.md),
[`d_eff_binary`](https://larry-leon.github.io/forestsearch/reference/d_eff_binary.md),
[`d_eff_survival`](https://larry-leon.github.io/forestsearch/reference/d_eff_survival.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# ---- Binary example (permutation mode) ----
set.seed(1)
df <- data.frame(
  id = 1:500, treat = rbinom(500, 1, 0.5),
  bm1 = as.factor(rbinom(500, 1, 0.7)),
  bm2 = as.factor(rbinom(500, 1, 0.5)),
  age = round(rnorm(500, 55, 10)),
  y   = rbinom(500, 1, 0.30)   # null: no treatment effect
)

fs_params <- list(
  confounders.name = c("bm1", "bm2", "age"),
  outcome.name = "y", event.name = "y",
  treat.name = "treat", id.name = "id",
  outcome_type = "binary", effect_measure = "OR",
  adverse_outcome = TRUE,
  pconsistency.threshold = 0.90, fs.splits = 400,
  n.min = 60, maxk = 2,
  use_lasso = TRUE, use_grf = TRUE,
  is.RCT = TRUE, details = FALSE, quiet = TRUE,
  parallel_args = list(plan = "sequential", workers = 1)
)

res <- fpr_calibration(
  c1 = 1.70, c2 = 1.70,
  df.analysis = df,
  fs_params = fs_params,
  nsims = 200,
  n_min = 60, outcome_type_approx = "binary", p_event = 0.30,
  sim_seed = 42, n_workers = 4
)
print(res)
} # }
```
