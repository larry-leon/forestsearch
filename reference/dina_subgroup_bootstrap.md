# Selection-adjusted bootstrap CI for dina_subgroup()

Performs selection-adjusted inference for the harm subgroup discovered
by
[`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md).
For each of `n_boot` bootstrap resamples of the input data frame, the
function refits DINA via
[`dina()`](https://larry-leon.github.io/forestsearch/reference/dina.md)
and re-runs the subgroup search via
[`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md),
then reports three classes of output:

## Usage

``` r
dina_subgroup_bootstrap(
  df,
  outcome,
  treatment,
  covariates,
  family = c("gaussian", "binomial", "poisson", "cox"),
  status = NULL,
  m_diff,
  n_min = 60L,
  n_min.frac = 0.1,
  direction = c("both", "left", "right"),
  sg_focus = "maxSG",
  effect_neighborhood = 0.1,
  alpha = 0.05,
  n_boot = 200L,
  parallel = c("none", "boots"),
  seed = NULL,
  verbose = FALSE,
  fit = NULL,
  sg = NULL,
  refit = TRUE,
  refit_strata = NULL,
  refit_confounders = "none",
  ...
)
```

## Arguments

- df:

  data frame containing the outcome, treatment, and covariate columns.

- outcome:

  character(1); name of the response column in `df`. For
  `family = "cox"` this is the event/censoring time column.

- treatment:

  character(1); name of the binary treatment indicator column.

- covariates:

  character vector; names of the covariate columns to search over. All
  must be numeric.

- family:

  one of `"gaussian"`, `"binomial"`, `"poisson"`, `"cox"`.

- status:

  character(1) or `NULL`; for `family = "cox"` only, the event indicator
  column name.

- m_diff:

  scalar harm threshold on the natural-parameter scale.

- n_min:

  positive integer, or `NULL`; minimum subgroup size. Default `60L`.
  Passing `n_min = NULL` opts into the sample-size-adaptive floor
  `max(60, ceiling(n_min.frac * nrow(df)))`, matching
  [`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md).
  When supplying a pre-computed `sg`, use the same `n_min` (or
  `n_min.frac`) it was found with, since the resolved value must match
  `sg$n_min`.

- n_min.frac:

  numeric in (0, 1); fraction of `nrow(df)` used for the adaptive
  `n_min` floor when `n_min = NULL`. Default `0.10`. Ignored when
  `n_min` is supplied.

- direction:

  one of `"both"` (default), `"left"`, `"right"`.

- sg_focus:

  character; subgroup selection criterion, forwarded to
  [`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md)
  for both the original-data point estimate and every bootstrap
  iteration. Accepts the same set as
  [`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md):
  one of `"maxSG"` (default), `"minSG"`, `"eff"`, `"effMaxSG"`,
  `"effMinSG"`, `"maxeff"`, `"maxeffCons"` (canonical `"hr"`,
  `"hrMaxSG"`, `"hrMinSG"` and the alias `"maxcons"` also accepted).
  `"maxeff"` and `"maxeffCons"` were previously rejected here even
  though
  [`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md)
  and
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  both accept them. On this path all five of `"eff"`, `"hr"`,
  `"maxcons"`, `"maxeff"` and `"maxeffCons"` rank identically – see
  [`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md)
  and
  [`fs_focus_tag()`](https://larry-leon.github.io/forestsearch/reference/fs_focus_tag.md).

- effect_neighborhood:

  numeric in `[0, 1)`; relative tolerance for the `"effMaxSG"` /
  `"effMinSG"` effect band, forwarded to
  [`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md).
  Default `0.10`. Ignored for the non-band foci.

- alpha:

  confidence level for the percentile CI. Default `0.05` (percentile
  endpoints at the `alpha/2` and `1 - alpha/2` quantiles).

- n_boot:

  positive integer; number of bootstrap iterations. Default `200L`.
  Warns if `n_boot < 100`.

- parallel:

  one of `"none"` (default) or `"boots"`. When `"boots"`, bootstrap
  iterations dispatch via `foreach::foreach() %dofuture%` against the
  active
  [`future::plan()`](https://future.futureverse.org/reference/plan.html).

- seed:

  optional integer seed for reproducibility of the bootstrap indices.
  Default `NULL`.

- verbose:

  logical; if `TRUE`, print a progress message every 10 iterations
  (serial mode) or a start summary (parallel mode). Default `FALSE`.

- fit:

  optional precomputed DINA fit (`"dina"` or `"dina_bagged"`) on the
  original `df`. When supplied (and `sg` is `NULL`), the function reuses
  it for the original-data point estimate instead of refitting, then
  runs only the subgroup search. Its `family` must match `family`.
  Default `NULL`.

- sg:

  optional precomputed `"dina_subgroup"` result on the original `df`.
  When supplied, it is used directly as the point estimate with no refit
  and no re-search; `fit` is then ignored. Its stored `m_diff`, `n_min`,
  `alpha`, `family`, and row count must match the corresponding
  arguments / `df`, or an error is raised. Supplying the upstream `sg`
  is the recommended way to guarantee the bootstrap's point estimate is
  identical to a standalone
  [`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md)
  result. Default `NULL`.

- refit:

  logical; if `TRUE` (default), also fit the *standard* within-subgroup
  treatment-effect model (Cox for survival, GLM otherwise) via
  [`dina_subgroup_refit`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup_refit.md)
  on the original data, and bootstrap it within the fixed discovered
  signature. Reported alongside the DINA effect. Skipped (with a
  warning) if the original-data subgroup is not found or the
  within-subgroup fit errors.

- refit_strata:

  `NULL` (default) or a character vector of column names entered as Cox
  [`strata()`](https://rdrr.io/pkg/survival/man/strata.html) terms in
  the within-subgroup model. Cox family only.

- refit_confounders:

  adjustment set for the within-subgroup model: `"none"` (default,
  unadjusted), `NULL` (automatic: the DINA `covariates` with the
  subgroup-defining covariate omitted), or a character vector of column
  names. See
  [`dina_subgroup_refit`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup_refit.md).

- ...:

  further arguments forwarded to
  [`dina()`](https://larry-leon.github.io/forestsearch/reference/dina.md)
  for each per-iteration DINA fit (e.g., `n_folds`, `cens_type`,
  `cens_params`, `propensity_method`, `baseline_method`).

## Value

An object of class `"dina_subgroup_bootstrap"`, a list with components:

- point:

  the `"dina_subgroup"` object computed on the original (unbootstrapped)
  data.

- effect_ci:

  named length-2 numeric vector (`lower`, `upper`) giving the bootstrap
  percentile CI on the fixed-subgroup linear functional `a*^T beta_b`,
  where `a* = (1, x_bar_S)` is the covariate-mean vector of the
  original-data subgroup and `beta_b` is each bootstrap iteration's DINA
  coefficient vector. This is the meaningful effect-uncertainty CI;
  width reflects sampling variance of `beta-hat`.

- effect_dist:

  numeric vector of length `n_boot` containing `a*^T beta_b` for each
  iteration (`NA` for iterations whose DINA fit failed).

- refit:

  the `"dina_subgroup_refit"` object (standard within-subgroup model on
  the original data), or `NULL` if `refit = FALSE`, no subgroup was
  found, or the fit failed.

- refit_effect_ci:

  named length-2 numeric vector (`lower`, `upper`); bootstrap percentile
  CI on the within-subgroup standard-model treatment effect, refit
  within the FIXED discovered signature on each resample. Conditional on
  the signature; comparable to `effect_ci`. `(NA, NA)` if `refit` was
  not performed.

- refit_effect_dist:

  numeric vector of length `n_boot` of per-iteration within-signature
  standard-model treatment effects (`NA` where not computed or the fit
  failed).

- n_subgroup_ci, threshold_ci:

  named length-2 numeric vectors giving percentile CIs for the
  structural quantities of bootstrap-selected subgroups, restricted to
  iterations that selected the same `(covariate, direction)` as the
  original-data subgroup. Subgroup-size and boundary-location stability
  diagnostics.

- n_modal_match:

  integer; number of bootstrap iterations whose selection matched the
  original-data `(covariate, direction)`. Denominator for the stability
  CIs.

- selection_frequency:

  a `table` of `(covariate, direction)` selection counts across found
  iterations. Covariate-selection stability diagnostic.

- boot_results:

  data frame with one row per bootstrap iteration; columns `found`,
  `covariate`, `direction`, `threshold`, `n_subgroup`, `mean_tau_hat`,
  `failed`. The `mean_tau_hat` column records the bootstrap-selected
  subgroup's boundary value (pinned near `m_diff` by the search rule)
  and is retained for diagnostic transparency rather than as a
  meaningful effect estimate.

- n_boot, n_boot_found, n_boot_failed:

  convergence diagnostics.

- alpha, m_diff, n_min, sg_focus, effect_neighborhood, family:

  echoed inputs.

- call:

  the matched call.

## Details

- An **effect CI** on `a*^T beta_b`, where `a* = (1, x_bar_S)` is held
  FIXED from the original-data subgroup `S*` and `beta_b` is each
  iteration's DINA coefficient vector. This is the bootstrap analog of
  the Wald CI returned by
  [`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md)
  – conditional on the discovered subgroup definition, with width
  reflecting sampling variance of `beta-hat`.

- **Stability CIs** on the bootstrap-selected `n_subgroup` and
  `threshold`, restricted to iterations that selected the same
  `(covariate, direction)` as the original-data subgroup.

- A **selection-frequency table** counting how often each
  `(covariate, direction)` was chosen. The covariate-selection stability
  diagnostic.

Each bootstrap iteration calls
[`dina()`](https://larry-leon.github.io/forestsearch/reference/dina.md)
(single-pass, sandwich variance) under the hood, not
[`dina_bagged()`](https://larry-leon.github.io/forestsearch/reference/dina_bagged.md).

If you have already run
[`dina()`](https://larry-leon.github.io/forestsearch/reference/dina.md)
and
[`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md)
on the original data, pass those objects via `fit` and/or `sg` so the
original-data point estimate reuses them rather than being recomputed
here. This guarantees the reported point estimate matches your
standalone result exactly; otherwise the internal refit uses its own
cross-fitting fold assignment, which at a boundary-case `m_diff` can
disagree with the upstream fit on whether a subgroup qualifies.

## Why the bootstrap target is the *fixed-subgroup* effect

A naive bootstrap of
[`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md)'s
`mean_tau_hat` is uninformative. Under the default `sg_focus = "maxSG"`
the search selects the *largest* subgroup whose mean tau-hat exceeds
`m_diff`, so the chosen subgroup lies at the boundary by construction:
adding one more patient would push the mean below `m_diff`. Each
bootstrap iteration's selected subgroup is similarly pinned, so the
bootstrap distribution of `mean_tau_hat` clusters tightly above `m_diff`
and the resulting CI quantifies threshold-grid granularity rather than
effect uncertainty. (Under the band foci `"effMaxSG"` / `"effMinSG"` the
per-iteration subgroup is instead anchored to the maximum effect rather
than the `m_diff` boundary, but the same caveat applies: `mean_tau_hat`
is a selection-pinned quantity, not an effect estimate.) The
fixed-subgroup linear functional `a*^T beta_b` instead targets the
clinically meaningful estimand "treatment effect for patients whose
covariates lie in the discovered subgroup". Because DINA's `beta-hat` is
cross-fit (and hence asymptotically unbiased even for data-driven
contrasts), the percentile distribution of `a*^T beta_b` gives a
non-parametric CI that complements the sandwich Wald CI of
[`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md).

## Two complementary effect estimates

When `refit = TRUE` (default), the function reports two treatment
effects for the discovered subgroup, on the same natural-parameter
scale:

- the **DINA effect** (`effect_ci`): the BLP-analog `a*^T beta`, a
  projection of the cross-fitted CATE surface averaged over the
  subgroup; and

- the **within-subgroup standard-model effect** (`refit_effect_ci`): the
  treatment coefficient from a plain Cox model (or GLM) fit on the
  subgroup, refit within the FIXED discovered signature on each
  resample. This is the contrast a clinical-trial reader usually
  expects.

The two coincide only under correct linearity; reporting both is
informative. Crucially, *both* CIs are conditional on the discovered
signature treated as pre-specified – neither corrects for the
data-driven selection of the signature. A selection-adjusted (de-biased)
interval requires the more involved bias-corrected bootstrap of Leon et
al. (2024, Section 3.2) and is not provided here; selection *stability*
is instead summarized by `selection_frequency`.

## Performance

With `parallel = "boots"`, each bootstrap iteration runs in a worker via
`foreach::foreach() %dofuture%`. Set
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
to a parallel backend (e.g.,
[`future::multisession`](https://future.futureverse.org/reference/multisession.html))
before calling, and reset to `"sequential"` afterwards. Workers run from
the installed package; `devtools::install()` is required after code
edits before the parallel branch will see them.

## See also

[`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md)
for the underlying threshold search;
[`dina()`](https://larry-leon.github.io/forestsearch/reference/dina.md)
for the per-iteration DINA fit;
[`dina_subgroup_refit`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup_refit.md)
for the standard within-subgroup model reported when `refit = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 400
df_demo <- data.frame(
  w  = stats::rbinom(n, 1, 0.5),
  x1 = stats::runif(n, -1, 1),
  x2 = stats::runif(n, -1, 1)
)
tau_x      <- 0.4 + 1.2 * df_demo$x1
df_demo$y  <- 0.5 * df_demo$x1 + df_demo$w * tau_x + stats::rnorm(n)

bs <- dina_subgroup_bootstrap(
  df         = df_demo,
  outcome    = "y",
  treatment  = "w",
  covariates = c("x1", "x2"),
  family     = "gaussian",
  m_diff     = 0.5,
  n_boot     = 50L,
  seed       = 1L
)
print(bs)
} # }
```
