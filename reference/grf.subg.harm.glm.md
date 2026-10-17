# GRF Subgroup Identification for GLM Outcomes

Identifies treatment effect subgroups for binary, continuous, and count
outcomes using Generalized Random Forests (GRF) for candidate factor
screening, followed by exhaustive subgroup enumeration with GLM-based
effect estimation.

## Usage

``` r
grf.subg.harm.glm(
  data,
  confounders.name,
  outcome.name,
  treat.name = "treat",
  id.name = "id",
  outcome_type = c("binary", "continuous", "count"),
  effect_measure = NULL,
  offset.name = NULL,
  overdispersion = c("none", "quasi", "negbin"),
  grf_count_transform = c("log", "identity"),
  n.min = 60L,
  dmin.grf = 0,
  frac.tau = 0.5,
  maxdepth = 2L,
  RCT = TRUE,
  sg.criterion = c("mDiff", "Nsg"),
  conf.level = 0.95,
  seedit = 8316951L,
  return_selected_cuts_only = FALSE,
  adverse_outcome = FALSE,
  tune_grf = FALSE,
  grf_selection = c("tree", "frontier"),
  frontier_rule = c("effMaxSG", "eff", "maxSG", "minSG", "effMinSG"),
  effect_neighborhood = 0.1,
  selection_rule = c("neighborhood", "pareto", "both"),
  details = FALSE,
  verbose = FALSE
)
```

## Arguments

- data:

  Data frame containing the analysis data.

- confounders.name:

  Character vector of potential effect modifier (confounder) names.
  These are the candidate split variables for GRF and the subgroup
  enumeration.

- outcome.name:

  Character. Name of the outcome column.

- treat.name:

  Character. Name of the binary treatment indicator (0/1). Default:
  `"treat"`.

- id.name:

  Character. Name of the subject ID column. Default: `"id"`.

- outcome_type:

  Character. Type of outcome. One of `"binary"` (outcome is 0/1; use
  `effect_measure` to select log-OR, log-RR, or RD), `"continuous"`
  (estimates mean difference), or `"count"` (count or rate outcome; use
  `effect_measure = "log_IRR"` with `offset.name` for incidence rate
  ratios). Default: `"binary"`.

- effect_measure:

  Character. Effect measure to estimate. One of `"log_OR"` (log odds
  ratio, binary/logistic), `"log_RR"` (log relative risk,
  binary/modified Poisson), `"RD"` (risk difference, binary/identity
  link), `"MD"` (mean difference, continuous/Gaussian), or `"log_IRR"`
  (log incidence rate ratio, count/Poisson with offset). Default is
  `"log_OR"` for binary, `"MD"` for continuous, and `"log_IRR"` for
  count.

- offset.name:

  Character or `NULL`. Name of the person-time (exposure) column. When
  supplied, `log(offset)` enters the linear predictor of the GLM as an
  offset term. Required when `outcome_type = "count"` and
  `effect_measure = "log_IRR"`. Ignored for binary and continuous
  outcomes. Default: `NULL`.

- overdispersion:

  Character. Overdispersion correction for count models. One of `"none"`
  (standard Poisson, model-based SEs), `"quasi"` (quasi-Poisson,
  dispersion-corrected SEs), or `"negbin"` (negative binomial via
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html)). Ignored
  for binary and continuous outcomes. Default: `"none"`.

- grf_count_transform:

  Character. Transformation applied to the count outcome *before*
  passing it to
  [`grf::causal_forest()`](https://rdrr.io/pkg/grf/man/causal_forest.html)
  for factor screening. Has no effect on binary or continuous outcomes.
  `"log"` (default) applies \\Y^\* = \log(Y + 0.5)\\, which prevents
  \\\log(0)\\ when zero counts are present and approximates variance
  stabilisation for Poisson data. `"identity"` passes raw counts \\Y^\*
  = Y\\ directly, appropriate when counts are large and rarely zero.
  Both options affect only GRF factor screening; all subgroup effect
  estimates use the correctly specified GLM regardless. Default:
  `"log"`.

- n.min:

  Integer. Minimum subgroup sample size for a valid split. Default:
  `60`.

- dmin.grf:

  Numeric. Minimum absolute CATE magnitude (on the doubly robust score
  scale) required for the identified subgroup to be declared valid. The
  score equals `-CATE(H)`, so a positive threshold requires the subgroup
  CATE to be sufficiently negative (treatment hurts). Default: `0`.

- frac.tau:

  Numeric. Fraction of the sample used for the GRF horizon (time-horizon
  analogue for non-survival outcomes). Default: `0.5`.

- maxdepth:

  Integer. Maximum depth for GRF policy trees. Default: `2`.

- RCT:

  Logical. Whether data come from a randomised trial. When `TRUE`,
  propensity scores are fixed at 0.5. Default: `TRUE`.

- sg.criterion:

  Character. Subgroup selection criterion. One of `"mDiff"` (maximum
  score difference) or `"Nsg"` (largest valid subgroup). Default:
  `"mDiff"`.

- conf.level:

  Numeric. Confidence level for subgroup effect estimates. Default:
  `0.95`.

- seedit:

  Integer. Random seed for GRF. Default: `8316951`.

- return_selected_cuts_only:

  Logical. If `TRUE`, returns only cuts from the tree depth that
  identified the selected subgroup. If `FALSE` (default), returns all
  cuts from all fitted trees. Matches the interface of
  [`grf.subg.harm.survival`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.survival.md).
  Default: `FALSE`.

- adverse_outcome:

  Logical. If `TRUE`, a higher outcome value indicates a worse result
  (e.g., adverse event count), so the "harm" subgroup is where treatment
  increases the outcome. If `FALSE` (default), higher values are better.
  Default: `FALSE`.

- tune_grf:

  Logical. If `TRUE`, enables cross-validated hyperparameter tuning for
  the causal forest via `tune.parameters = "all"`, which tunes
  `min.node.size`, `mtry`, `sample.fraction`, `alpha`, and
  `imbalance.penalty`. If `FALSE` (default), all GRF parameters remain
  at their defaults, preserving existing behavior. Most impactful for
  observational data where nuisance models benefit from tuning; defaults
  are near-optimal for RCTs (Dandl et al., 2024).

- grf_selection:

  Character, one of `"tree"` (default) or `"frontier"`. `"tree"` uses
  the GRF policy tree (standard behavior). `"frontier"` is an
  **experimental** alternative that, on the same doubly-robust scores,
  enumerates single-covariate thresholds (and depth-2 covariate-pair
  conjunctions when `maxdepth >= 2`), takes the Pareto frontier of
  (harm-effect, size), and selects one subgroup under `frontier_rule`.
  The selected subgroup is always a single conjunction (never the
  disjunctive union the tree path can return). Provided for comparison
  and exploration, not as a recommended default.

- frontier_rule:

  Character, one of `"effMaxSG"` (default), `"eff"`, `"maxSG"`,
  `"minSG"`, or `"effMinSG"`; the rule applied to the frontier when
  `grf_selection = "frontier"`. See
  [`grf.subg.harm.survival`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.survival.md).

- effect_neighborhood:

  Numeric in (0, 1); relative neighborhood for the `"effMaxSG"` /
  `"effMinSG"` rules. Default 0.10. Used only when
  `grf_selection = "frontier"`.

- selection_rule:

  Character, one of `"neighborhood"` (default), `"pareto"`, `"both"`;
  the inclusion band for the `eff*SG` frontier rules. Inert for the
  other rules.

- details:

  Logical. Print GRF diagnostic information. Default: `FALSE`.

- verbose:

  Logical. Print progress messages. Default: `FALSE`.

## Value

A list of class `"grf_glm_result"` containing:

- sg.harm.id:

  **Character vector of cut expressions** defining the identified
  subgroup (e.g., `c("v1=1", "v3=1")`), length equal to the depth of the
  selected policy tree. `NULL` if no subgroup was found. **Not** a
  per-subject membership indicator; see the *Field naming collision*
  section below.

- data:

  The input data frame with a `treat.recommend` column added (`0` = in
  harm/questionable subgroup, `1` = complement).

- tree.cuts:

  Named list of GRF policy-tree split points.

- grf_varimp:

  Named numeric vector of GRF variable importances.

- outcome.name:

  Character: name of the outcome column.

- treat.name:

  Character: name of the treatment column.

- id.name:

  Character: name of the subject-identifier column.

- confounders.name:

  Character vector of covariate column names.

- effect_measure:

  Character: the effect measure used.

- outcome_type:

  Character: the outcome type.

- overdispersion:

  Character: the overdispersion correction applied.

- offset.name:

  Character or `NULL`: the offset variable name.

- grf_count_transform:

  Character: the transform argument as matched (`"log"` or
  `"identity"`).

- grf_y_transform:

  Character: transformation actually applied to Y for GRF screening
  (`"log(Y + 0.5)"`, `"identity"`, or `"none"`).

## Details

For count / rate outcomes, the GRF screening step applies a
variance-stabilising \\\log(Y + 0.5)\\ transformation so that
[`grf::causal_forest()`](https://rdrr.io/pkg/grf/man/causal_forest.html)
can be applied without dedicated count-specific infrastructure. All
subgroup effect estimates are computed from correctly specified Poisson,
quasi-Poisson, or negative-binomial GLMs via
[`create_glm_row`](https://larry-leon.github.io/forestsearch/reference/create_glm_row.md).

### GRF Screening for Count Outcomes

[`grf::causal_forest()`](https://rdrr.io/pkg/grf/man/causal_forest.html)
is designed for continuous outcomes. For count data, this function
transforms \\Y\\ before passing it to GRF according to the
`grf_count_transform` argument:

- `"log"` (default):

  Applies `Y_star = log(Y + 0.5)`. The +0.5 shift avoids `log(0)` for
  zero-count observations; the log scale approximates variance
  stabilisation for Poisson data and is recommended when any counts are
  zero or small.

- `"identity"`:

  Passes raw counts `Y_star = Y` directly. Appropriate when counts are
  large, rarely zero, and the raw scale is near-continuous (e.g., daily
  hospitalisation volumes in a large centre).

In both cases the transformation applies *only* to GRF factor screening
and policy-tree cut generation – it determines *which variables* and
*which cut-points* are candidate subgroup definitions. The actual effect
estimates in each subgroup are always computed from a correctly
specified Poisson (or quasi-Poisson / negative-binomial) GLM via
[`create_glm_row`](https://larry-leon.github.io/forestsearch/reference/create_glm_row.md),
regardless of `grf_count_transform`.

Alternatively, setting `use_grf = FALSE` in the parent
[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
call routes to LASSO-based factor selection, which supports
`family = "poisson"` natively via glmnet.

### Offset and Person-Time

When `offset.name` is supplied, \\\log(\text{exposure}\_i)\\ enters
every GLM fit in the subgroup enumeration loop, so that subgroup effects
are on the incidence-rate-ratio scale adjusted for differential
follow-up.

## Field naming collision with forestsearch / subgroup.consistency

The field name `sg.harm.id` has **different semantics** on this object
versus on the `grp.consistency` list returned by
[`subgroup.consistency`](https://larry-leon.github.io/forestsearch/reference/subgroup.consistency.md)
(and nested inside the
[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
result):

|  |  |  |
|----|----|----|
| **Object** | **`sg.harm.id` contains** | **Length / type** |
| `grf.subg.harm.glm()` result (this function) | character vector of cut expressions | character, length = depth of selected tree |
| [`subgroup.consistency()`](https://larry-leon.github.io/forestsearch/reference/subgroup.consistency.md) result | per-subject 0/1 membership indicator | integer, length `nrow(df)` |

**Practical consequence.** On this object, pasting
`paste(obj$sg.harm.id, collapse = " & ")` correctly renders the
identified subgroup. On a `grp.consistency` list the same expression
silently concatenates a long 0/1 indicator and produces output like
`"0 & 0 & 1 & 0 & ..."`. When writing code that must handle both object
types, dispatch on class first, or use the top-level `sg.harm` field
that the
[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
result exposes.

This naming collision is a documented CRAN-stable API for v0.1.x and
v0.2.x; it is expected to be resolved via a deprecation cycle in a
future minor release.

## See also

[`grf.subg.harm.survival`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.survival.md)
for time-to-event outcomes,
[`create_glm_row`](https://larry-leon.github.io/forestsearch/reference/create_glm_row.md)
for per-subgroup GLM estimation,
[`glm_effect_profile`](https://larry-leon.github.io/forestsearch/reference/glm_effect_profile.md)
for continuous biomarker profiles.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate count data with heterogeneous treatment effect
set.seed(123)
n <- 600
z1 <- rbinom(n, 1, 0.4)   # binary subgroup factor
df_count <- data.frame(
  id          = seq_len(n),
  treat       = rbinom(n, 1, 0.5),
  v1          = as.factor(z1),
  v2          = as.factor(rbinom(n, 1, 0.5)),
  person_time = runif(n, 0.5, 2.0),
  age         = rnorm(n, 60, 10)
)
# Generate outcome: higher event rate for treated in z1 = 1 subgroup
log_rate <- -1.5 +
  0.8 * df_count$treat * (df_count$v1 == 1) -
  0.3 * df_count$treat * (df_count$v1 == 0) +
  log(df_count$person_time)
df_count$events <- rpois(n, lambda = exp(log_rate))

# Run GRF subgroup identification for count outcome (log transform, default)
grf_glm <- grf.subg.harm.glm(
  data                = df_count,
  confounders.name    = c("v1", "v2", "age"),
  outcome.name        = "events",
  treat.name          = "treat",
  id.name             = "id",
  outcome_type        = "count",
  effect_measure      = "log_IRR",
  offset.name         = "person_time",
  overdispersion      = "quasi",
  grf_count_transform = "log",      # default: log(Y + 0.5)
  n.min               = 60,
  verbose             = TRUE
)
print(grf_glm$sg.harm.id)

# Alternative: identity transform (large counts, no zeros)
grf_glm_id <- grf.subg.harm.glm(
  data                = df_count,
  confounders.name    = c("v1", "v2", "age"),
  outcome.name        = "events",
  treat.name          = "treat",
  id.name             = "id",
  outcome_type        = "count",
  effect_measure      = "log_IRR",
  offset.name         = "person_time",
  overdispersion      = "none",
  grf_count_transform = "identity",  # raw counts passed to GRF
  n.min               = 60
)
} # }
```
