# Identify a harm subgroup from a DINA fit via univariate threshold search

Given a fitted DINA object (`"dina"` or `"dina_bagged"`) and a data
frame containing the covariates the model was fit on, search over
(covariate, threshold, direction) candidates for subgroups whose mean
per-patient HTE estimate meets the harm threshold `m_diff`, then select
one according to `sg_focus` (default `"maxSG"`: the largest qualifying
subgroup). This is the simplest bridge between DINA's continuous
estimator \\\hat\tau(x) = \hat\beta_0 + x'\hat\beta\\ and forestsearch's
interpretable signature framework.

## Usage

``` r
dina_subgroup(
  fit,
  df,
  covariates,
  m_diff,
  n_min = 60L,
  n_min.frac = 0.1,
  direction = c("both", "left", "right"),
  max_depth = 2L,
  grid_probs = seq(0.1, 0.9, by = 0.1),
  sg_focus = "maxSG",
  selection_rule = "neighborhood",
  effect_neighborhood = 0.1,
  alpha = 0.05,
  tau_sign = 1
)
```

## Arguments

- fit:

  a fitted DINA object of class `"dina"` or `"dina_bagged"`, as returned
  by
  [`dina_fit()`](https://larry-leon.github.io/forestsearch/reference/dina_fit.md),
  [`dina_fit_bagged()`](https://larry-leon.github.io/forestsearch/reference/dina_fit_bagged.md),
  [`dina()`](https://larry-leon.github.io/forestsearch/reference/dina.md),
  or
  [`dina_bagged()`](https://larry-leon.github.io/forestsearch/reference/dina_bagged.md).

- df:

  data frame containing the covariate columns referenced by
  `covariates`. Typically the same data frame used to fit `fit`.

- covariates:

  character vector of column names in `df` to search over. All must be
  numeric.

- m_diff:

  scalar harm threshold on the natural-parameter scale of `fit$family`.
  Subgroups are kept only if their mean tau-hat meets or exceeds this
  value.

  **It is a candidacy floor, not a ranking key.** It decides which
  candidates enter the qualifying family; `sg_focus` then orders that
  family. Every focus ranks only over survivors of this floor, so
  widening or tightening `m_diff` changes the selection by changing the
  family, not by changing the sort.

  **Scale.** The link (tau) scale of `fit$family`, which is not the
  scale the effect is usually reported on: `log(HR)` for `cox`,
  `log(OR)` for `binomial`, `log(IRR)` for `poisson`, and the raw mean
  difference for `gaussian`. A harm threshold of HR 1.25 is therefore
  `m_diff = log(1.25)`, not `1.25`.

  **Under
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  with `subgroup_method = "dina"` this argument is derived, not passed
  through.** It is computed as `log(hr.threshold)` (identity for
  `gaussian`) from the resolved `effect.threshold`, and an `m_diff`
  supplied in `dina_args` – a
  [`dina_frontier()`](https://larry-leon.github.io/forestsearch/reference/dina_frontier.md)
  key – is **ignored** on that path. Set the floor with
  `effect.threshold` / `hr.threshold` instead.

  **Under `dina_select_statistic = "effect"`** the floor still defines
  the qualifying family on the tau-hat scale, while selection re-ranks
  that family on the refit inferential effect. The two scales therefore
  do different jobs in one selection: `m_diff` shapes the family, the
  refit effect orders it. Since multiplier resampling de-biases over the
  realized selection event, `m_diff` is what fixes the family MR
  de-biases over – changing it changes the inferential target, not just
  the winner.

- n_min:

  positive integer, or `NULL`; minimum subgroup size. Default `60L`.
  Supplying a value (or omitting it) uses that fixed floor. Passing
  `n_min = NULL` opts into a sample-size-adaptive floor
  `max(60, ceiling(n_min.frac * n))`, where `n = nrow(df)` – at least
  60, and at least `n_min.frac` of `n`. The resolved integer is stored
  on the returned object (`$n_min`), so
  [`dina_subgroup_bootstrap()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup_bootstrap.md)
  and
  [`dina_subgroup_refit()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup_refit.md)
  inherit it.

- n_min.frac:

  numeric in (0, 1); fraction of `nrow(df)` used for the adaptive
  `n_min` floor when `n_min = NULL`. Default `0.10`. Ignored when
  `n_min` is supplied. Default `60L`, matching the convention in GRF and
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md).

- direction:

  one of `"both"` (default), `"left"`, or `"right"`. For each covariate
  `x_j` controls whether to search subgroups of the form `x_j <= q`
  (left), `x_j >= q` (right), or both.

- max_depth:

  integer, `1L` or `2L` (default `2L`, aligning with the consistency
  search `maxk = 2` and GRF `grf_depth = 2`). `1L` searches
  single-covariate cuts only (the legacy behaviour). `2L` additionally
  searches AND-conjunctions of two distinct covariates and ranks them
  jointly with the single-covariate cuts; see the *Depth-2 conjunctions*
  section.

- grid_probs:

  numeric vector of probabilities in `[0, 1]` giving the per-covariate
  quantile grid used to generate depth-2 pair thresholds. Only consulted
  when `max_depth = 2L`. Default interior deciles `seq(0.1, 0.9, 0.1)`.
  A finer grid recovers more of the exhaustive pair space at higher
  cost; covariates with fewer unique values than `length(grid_probs)`
  are searched at their unique values instead.

- sg_focus:

  character; subgroup selection criterion applied to the qualifying
  candidates. One of `"maxSG"` (default), `"minSG"`, `"eff"`,
  `"effMaxSG"`, `"effMinSG"`, `"maxeff"`, `"maxeffCons"`. The canonical
  forms `"hr"`, `"hrMaxSG"`, `"hrMinSG"` (used internally by
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md))
  and the alias `"maxcons"` for `"eff"` are also accepted. See
  Description for the semantics of each.

  **Five spellings, one rule.** On this path `"eff"`, `"hr"`,
  `"maxcons"`, `"maxeff"` and `"maxeffCons"` all resolve to the same
  ordering, `order(-eff)`. Two separate collapses produce that: `"eff"`
  and `"maxcons"` are aliases of `"hr"` everywhere (see
  `.normalize_sg_focus()`), and `"maxeff"` / `"maxeffCons"` are
  behavioural synonyms of `"hr"` here specifically. Use
  [`fs_focus_tag()`](https://larry-leon.github.io/forestsearch/reference/fs_focus_tag.md)
  to label output by the rule that ran; a hand-written tag reading
  `"maxeffCons"` on this path names a rule that did not execute.
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  announces the collapse at run time, gated on its `quiet` argument;
  `dina_subgroup()` called directly does not, having no `quiet` of its
  own.

  `"maxeff"` and `"maxeffCons"` are **synonyms here**, both ranking the
  qualifying candidates by effect (ties broken by insertion order)
  exactly as `"eff"` / `"hr"` does. The `Cons` suffix names the
  consistency floor, and DINA has none – it computes no `Pcons`, and no
  ordering key below contains one – so there is nothing for the
  qualifier to bind to. This is engine-specific: under
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  with `subgroup_method = "consistency"` the two are genuinely distinct,
  since `"maxeffCons"` retains both the consistency and effect floors
  while `"maxeff"` disables both. The synonymy disables no floor –
  candidates are qualified before this ordering is applied – and raises
  no condition, by design.

- selection_rule:

  character; rule defining the candidate inclusion band for `"effMaxSG"`
  / `"effMinSG"`. One of `"neighborhood"` (default; 1-D effect band),
  `"pareto"` (2-D non-dominated frontier in (effect, N)), or `"both"`
  (intersection). Must be `"neighborhood"` for the other foci. Forwarded
  unchanged to the shared band computation, so semantics match
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md).

- effect_neighborhood:

  numeric in `[0, 1)`; relative tolerance for the effect band. A
  candidate is in the neighborhood iff its (natural-scale) effect is at
  least `(1 - effect_neighborhood) * max(effect)` over the qualifying
  set. Default `0.10` (within 10% of the strongest effect). Consulted by
  `selection_rule = "neighborhood"` and `"both"`; ignored by `"pareto"`
  (and must be left at its default there) and for `"maxSG"`, `"minSG"`,
  `"eff"`.

- alpha:

  confidence level for the Wald interval on the subgroup-mean tau-hat.
  Default `0.05`.

- tau_sign:

  `1` (default) or `-1`; the orientation applied to the per-patient
  tau-hat before the `m_diff` floor, the `sg_focus` ranking and the
  reported mean and interval. `1` uses the fit's own scale. `-1` negates
  it, for a fit whose harm direction is the negative effect: the floor
  then reads "harm of at least `m_diff`".
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  passes `-1` for continuous and binary outcomes with
  `adverse_outcome = FALSE`, the same outcomes its effect estimator
  flips (`-Y`, `1 - Y`), so DINA's proposal floor and the admission
  floor point the same way. With `-1`, `mean_tau_hat`, `ci` and
  `candidates$tau_hat` are on the negated scale. The recorded `call`
  carries `tau_sign` only when it is not `1`.

## Value

An object of class `"dina_subgroup"`, a list with components:

- found:

  logical; `TRUE` if any candidate satisfied `m_diff` and `n_min`,
  otherwise `FALSE`. When `FALSE`, the subgroup-specific components
  below are `NULL`.

- covariate, direction, threshold:

  the chosen cut(s). Length 1 for a single-covariate subgroup; length 2
  (one entry per covariate) when `max_depth = 2L` selects a conjunction.

- labels:

  character vector of the chosen `{var op q}` factor label(s), in the
  form consumed by
  [`get_dfpred()`](https://larry-leon.github.io/forestsearch/reference/get_dfpred.md)
  (AND-composed). Length matches `covariate`.

- label:

  a single combined display string, e.g. `"{nodes >= 10 & age >= 60}"`.

- depth:

  integer; number of covariates in the selected subgroup (`1` or `2`).

- n_subgroup:

  integer size of the chosen subgroup.

- mean_tau_hat:

  scalar subgroup-mean tau-hat (in the `tau_sign` orientation).

- se_mean_tau_hat:

  Wald standard error, computed as `sqrt(a_S^T vcov(fit) a_S)`
  (CONDITIONAL on the chosen subgroup – not selection-adjusted).

- ci:

  named length-2 vector `(lower, upper)` giving the Wald 1 - alpha CI.

- mask:

  logical vector of length `nrow(df)` marking which rows are in the
  chosen subgroup.

- m_diff, n_min, sg_focus, selection_rule, effect_neighborhood, alpha,
  family, max_depth:

  the inputs / fit family, echoed for reproducibility.

- n_total:

  `nrow(df)`.

- n_candidates_searched:

  total number of `(j, dir, q)` triples evaluated.

- n_candidates_qualifying:

  number of those triples that met both the size and `m_diff`
  constraints.

- call:

  the matched call.

## Details

For each covariate `x_j` in `covariates`, the search visits the sorted
unique values of `df[[x_j]]` as candidate thresholds `q` and for each
threshold considers the left-tail subgroup of patients with \\x\_{j,i}
\le q\\ (and optionally also the right-tail subgroup with \\x\_{j,i} \ge
q\\). A candidate subgroup must satisfy:

- subgroup size `|S| >= n_min`, AND

- mean per-patient tau-hat over `S` `>= m_diff`.

Among the qualifying candidates, `sg_focus` selects the returned
subgroup, mirroring
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s
selection vocabulary:

- `"maxSG"` (default) – largest qualifying subgroup (legacy behaviour;
  size ties broken by the most extreme mean tau-hat).

- `"minSG"` – smallest qualifying subgroup.

- `"eff"` – most extreme effect, ignoring size.

- `"effMaxSG"` – largest subgroup in the candidate inclusion band (see
  `selection_rule`).

- `"effMinSG"` – smallest subgroup in the inclusion band. The `"eff*"`
  band foci reuse the inclusion-band logic shared with
  [`sort_subgroups()`](https://larry-leon.github.io/forestsearch/reference/sort_subgroups.md);
  for ratio families (binomial, poisson, cox) the band is applied on the
  natural (exponentiated) effect scale. `selection_rule` defines the
  band, exactly as in
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md):
  `"neighborhood"` (1-D effect band within `effect_neighborhood` of the
  maximum qualifying effect), `"pareto"` (the 2-D non-dominated frontier
  in (effect, N) space), or `"both"` (their intersection). DINA's
  univariate (covariate, threshold, direction) candidates populate the
  (effect, N) plane just as forestsearch's signatures do, so the
  frontier is defined analogously.

## Depth-2 conjunctions

With `max_depth = 2` the candidate set is extended to AND-conjunctions
of two thresholds on *distinct* covariates, e.g.
`nodes >= q1 & age >= q2`. Because tau-hat is linear in X, the mean over
the intersection mask is scored exactly like a depth-1 cut, so the
`m_diff`, `n_min`, and `sg_focus` machinery is unchanged. To stay
tractable, depth-2 pairs are generated over a per-covariate quantile
grid (`grid_probs`) rather than every unique value, while depth-1
singletons retain full resolution. Singletons and pairs are ranked
*jointly* under `sg_focus`, so a depth-2 conjunction is returned only
when it outranks every depth-1 candidate; otherwise the depth-1
selection is recovered exactly. When the selected subgroup is a
conjunction, `covariate`, `direction`, and `threshold` are length-2
vectors and `labels` holds the two `{var op q}` factors in the form
consumed by
[`get_dfpred()`](https://larry-leon.github.io/forestsearch/reference/get_dfpred.md).

The variance of the subgroup-mean tau-hat is computed as
\$\$\widehat{\mathrm{Var}}(\bar{\hat\tau}\_S) = a_S^\top
\widehat{\mathrm{Var}}(\hat\beta) a_S,\$\$ with \\a_S = (1,
\bar{x}\_S^\top)\\ and \\\widehat{\mathrm{Var}}(\hat\beta) =
\mathrm{vcov}(\mathtt{fit})\\. The variance source therefore tracks the
fit class: sandwich for `"dina"` and infinitesimal jackknife for
`"dina_bagged"`.

## Inference caveat

The returned Wald confidence interval is *conditional on the selected
(covariate, threshold) being treated as pre-specified*. It does not
adjust for the selection across the candidate set and will undercover
when the search space is large. Use this CI for descriptive reporting;
for hypothesis testing against `m_diff`, a selection-adjusted interval
via bootstrap of the full search procedure is required. The
undercoverage is more pronounced with `max_depth = 2L`, whose candidate
space is substantially larger.

## See also

[`dina_fit()`](https://larry-leon.github.io/forestsearch/reference/dina_fit.md)
/
[`dina()`](https://larry-leon.github.io/forestsearch/reference/dina.md)
for the underlying estimator;
[`dina_fit_bagged()`](https://larry-leon.github.io/forestsearch/reference/dina_fit_bagged.md)
/
[`dina_bagged()`](https://larry-leon.github.io/forestsearch/reference/dina_bagged.md)
for the bagged variant whose IJ variance propagates through the
subgroup-mean SE here;
[`sort_subgroups()`](https://larry-leon.github.io/forestsearch/reference/sort_subgroups.md)
and
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
for the selection vocabulary mirrored by `sg_focus`.

## Examples

``` r
set.seed(1)
n <- 400
df_demo <- data.frame(
  w  = stats::rbinom(n, 1, 0.5),
  x1 = stats::runif(n, -1, 1),
  x2 = stats::runif(n, -1, 1),
  x3 = stats::runif(n, -1, 1)
)
tau_x      <- 0.3 + 1.2 * df_demo$x1 - 0.4 * df_demo$x2
df_demo$y  <- 0.5 * df_demo$x1 + df_demo$w * tau_x + stats::rnorm(n)

fit <- dina(df_demo, outcome = "y", treatment = "w",
            covariates = c("x1", "x2", "x3"),
            family = "gaussian", seed = 1L)

sg <- dina_subgroup(fit, df_demo,
                    covariates = c("x1", "x2", "x3"),
                    m_diff = 0.5, n_min = 60L)
sg
#> DINA subgroup search
#>   Family:      gaussian
#>   Focus:       maxSG (largest qualifying subgroup)
#>   m_diff:      0.5
#>   n_min:       60
#>   Candidates:  3372 searched, 602 satisfied m_diff and size
#> 
#>   Signature:   x1 >= -0.6892
#>   Subgroup n:  341 / 400 (85.2%)
#>   mean tau-hat:0.503 (SE 0.1066)
#>   Wald 95% CI:  [0.2941, 0.712]
#> 
#>   Note: CI is conditional on the selected (covariate, threshold)
#>   and does NOT adjust for the selection across candidates.

# Tighter subgroup via the effect-neighborhood band:
dina_subgroup(fit, df_demo, covariates = c("x1", "x2", "x3"),
              m_diff = 0.5, n_min = 60L, sg_focus = "effMaxSG")
#> DINA subgroup search
#>   Family:      gaussian
#>   Focus:       effMaxSG (largest within effect band)
#>   Inclusion:   within 10% of max effect
#>   m_diff:      0.5
#>   n_min:       60
#>   Candidates:  3372 searched, 602 satisfied m_diff and size
#> 
#>   Signature:   x1 >= 0.4627
#>   Subgroup n:  109 / 400 (27.2%)
#>   mean tau-hat:1.147 (SE 0.1612)
#>   Wald 95% CI:  [0.8308, 1.463]
#> 
#>   Note: CI is conditional on the selected (covariate, threshold)
#>   and does NOT adjust for the selection across candidates.

# Frontier-aware selection (2-D non-dominated set in (effect, N)):
dina_subgroup(fit, df_demo, covariates = c("x1", "x2", "x3"),
              m_diff = 0.5, n_min = 60L, sg_focus = "effMaxSG",
              selection_rule = "both")
#> DINA subgroup search
#>   Family:      gaussian
#>   Focus:       effMaxSG (largest within effect band)
#>   Inclusion:   within 10% of max effect AND Pareto frontier
#>   m_diff:      0.5
#>   n_min:       60
#>   Candidates:  3372 searched, 602 satisfied m_diff and size
#> 
#>   Signature:   x1 >= 0.4627
#>   Subgroup n:  109 / 400 (27.2%)
#>   mean tau-hat:1.147 (SE 0.1612)
#>   Wald 95% CI:  [0.8308, 1.463]
#> 
#>   Note: CI is conditional on the selected (covariate, threshold)
#>   and does NOT adjust for the selection across candidates.

# Depth-2: allow AND-conjunctions of two covariates.
dina_subgroup(fit, df_demo, covariates = c("x1", "x2", "x3"),
              m_diff = 0.5, n_min = 60L, max_depth = 2L)
#> DINA subgroup search
#>   Family:      gaussian
#>   Focus:       maxSG (largest qualifying subgroup)
#>   m_diff:      0.5
#>   n_min:       60
#>   Candidates:  3372 searched, 602 satisfied m_diff and size
#> 
#>   Signature:   x1 >= -0.6892
#>   Subgroup n:  341 / 400 (85.2%)
#>   mean tau-hat:0.503 (SE 0.1066)
#>   Wald 95% CI:  [0.2941, 0.712]
#> 
#>   Note: CI is conditional on the selected (covariate, threshold)
#>   and does NOT adjust for the selection across candidates.
```
