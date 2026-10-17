# Evaluate Subgroup Consistency

Evaluates candidate subgroups using split-sample consistency validation.
For each candidate, repeatedly splits the data and checks whether the
treatment effect direction is consistent across splits.

## Usage

``` r
subgroup.consistency(
  df,
  hr.subgroups,
  hr.threshold = 1,
  hr.consistency = 1,
  pconsistency.threshold = 0.9,
  m1.threshold = Inf,
  n.splits = 100,
  details = FALSE,
  by.risk = 12,
  plot.sg = FALSE,
  maxk = 7,
  Lsg,
  confs_labels,
  sg_focus = "hr",
  selection_rule = "neighborhood",
  effect_neighborhood = 0.1,
  stop_Kgroups = 200,
  stop_threshold = NULL,
  show_candidate_summary = FALSE,
  max_print = 10,
  pconsistency.digits = 2,
  seed = 8316951,
  checking = FALSE,
  use_twostage = FALSE,
  twostage_args = list(),
  parallel_args = list(),
  estimator_fn = NULL,
  consistency_threshold = NULL,
  effect_label = "HR",
  effect_log_scale = FALSE,
  adjust_covariates = NULL,
  consistency_method = c("resample", "split"),
  glm_resample_spec = NULL
)
```

## Arguments

- df:

  Data frame containing the analysis dataset. Must include columns for
  outcome (Y), event indicator (Event), and treatment (Treat).

- hr.subgroups:

  Data.table of candidate subgroups from subgroup search, containing
  columns: HR, n, E, K, d0, d1, m0, m1, grp, and factor indicators.

- hr.threshold:

  Numeric. `c1`, the **screening threshold** on a candidate's own
  effect, re-applied here to the supplied `hr.subgroups` table. On the
  comparison scale of the estimand
  ([`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  passes the already-resolved value for GLM outcomes). Default: 1.0

- hr.consistency:

  Numeric. `c2`, **an effect threshold** – the minimum effect required
  in *each* half of a split for that split to count as consistent.
  **Not** a proportion and **not** the consistency rate; for the
  proportion see `pconsistency.threshold`. Default: 1.0

- pconsistency.threshold:

  Numeric. `p*`, **a proportion** in `[0, 1]` – the minimum fraction of
  splits that must be consistent, where consistency of one split is
  decided by `hr.consistency` (`c2`). `c2` sets the bar; `p*` counts how
  often it is met. Both names contain "consistency" and they are
  different quantities. Default: 0.9

- m1.threshold:

  Numeric. Maximum m1 threshold for filtering. Default: Inf

- n.splits:

  Integer. Number of splits for consistency evaluation. Default: 100

- details:

  Logical. Print progress details. Default: FALSE

- by.risk:

  Numeric. Risk interval for KM plots. Default: 12

- plot.sg:

  Logical. Generate subgroup plots. Default: FALSE

- maxk:

  Integer. Maximum number of factors in subgroup. Default: 7

- Lsg:

  List of subgroup parameters.

- confs_labels:

  Character vector mapping factor names to labels.

- sg_focus:

  Character. Subgroup selection criterion. One of:

  `"hr"` (or `"eff"`)

  :   Sort by \\(-Pcons, -hr, K)\\; pick top.

  `"maxeff"`

  :   Sort by \\(-hr, K)\\; pick top. No consistency filter. Consistency
      (`Pcons`) does not enter the sort key, so it is evaluated for the
      SELECTED winner only (reported in the footer) and is `NA` for
      every other candidate – the effect (from the search) is carried
      for all, but the per-candidate consistency refits are skipped.
      Selection is identical to evaluating every candidate.

  `"maxSG"`

  :   Sort by \\(-N, -Pcons, K)\\; pick top.

  `"minSG"`

  :   Sort by \\(N, -Pcons, K)\\; pick top.

  `"hrMaxSG"` (or `"effMaxSG"`)

  :   Among candidates with effect size within `effect_neighborhood` of
      the maximum, pick the one with the largest \\N\\.

  `"hrMinSG"` (or `"effMinSG"`)

  :   Among candidates with effect size within `effect_neighborhood` of
      the maximum, pick the one with the smallest \\N\\.

  Default: `"hr"`. The `"eff*"` forms are aliases for the `"hr*"` forms
  and read more naturally in GLM contexts (continuous MD, binary
  OR/RR/RD, count IRR). Both vocabularies produce identical results.

- selection_rule:

  Character. Rule defining the candidate inclusion set for `"hrMaxSG"` /
  `"hrMinSG"`. One of `"neighborhood"` (default; current behaviour),
  `"pareto"` (2D Pareto-non-dominated set in (effect, N)), or `"both"`
  (intersection of the two). Forwarded to
  [`sort_subgroups_preview`](https://larry-leon.github.io/forestsearch/reference/sort_subgroups_preview.md)
  (pre-consistency triage) and
  [`sort_subgroups`](https://larry-leon.github.io/forestsearch/reference/sort_subgroups.md)
  (post-consistency selection). Must be `"neighborhood"` for other
  `sg_focus` values.

- effect_neighborhood:

  Numeric in `[0, 1)`. Relative tolerance for `"hrMaxSG"` and
  `"hrMinSG"`. A candidate is in the *effect-size neighborhood* iff its
  (natural-scale) effect is at least
  `(1 - effect_neighborhood) * max(effect)`. Default `0.10` (i.e.,
  within 10\\ For GLM ratio measures (OR, IRR), the neighborhood test is
  applied on the natural scale (after exponentiation of `hr`). Ignored
  when `sg_focus` is `"hr"`, `"maxSG"`, or `"minSG"`. When
  `selection_rule = "pareto"`, must be left at its default (the rule
  does not consult it).

- stop_Kgroups:

  Integer. Maximum number of candidates to evaluate. When the ranked
  pool exceeds this cap it is truncated and a warning naming both counts
  and the `sg_focus` ordering is emitted. Default: 200

- stop_threshold:

  Numeric in `[0, 1]` or `NULL`. When a candidate subgroup's consistency
  probability (`Pcons`) meets or exceeds this threshold, evaluation
  stops early – remaining candidates are skipped. Set to `NULL` to
  disable early stopping and evaluate all candidates up to
  `stop_Kgroups`. Default: `NULL`.

  **Note:** Values \> 1.0 are not permitted. To disable early stopping,
  use `stop_threshold = NULL`, not a value above 1.

  **Interaction with `sg_focus`:**

  `"hr"`, `"maxSG"`, `"minSG"`

  :   Early stopping is valid because candidates are sorted by a single
      criterion. The first candidate passing the threshold is optimal
      under that criterion.

  `"hrMaxSG"`, `"hrMinSG"`

  :   Should generally be `NULL`, because these compound criteria
      require comparing HR *and* size across all candidates.
      [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
      automatically resets to `NULL` with a warning for these.

  For parallel execution, early stopping is checked after each batch
  completes, so some additional candidates beyond the first meeting the
  threshold may be evaluated. Use a smaller `batch_size` in
  `parallel_args` for finer-grained early stopping.

- show_candidate_summary:

  Logical. If `TRUE`, prints two diagnostic tables: a
  **pre-consistency** candidate preview (all top-`stop_Kgroups`
  candidates entering consistency evaluation, with Frontier and InBand
  flags) and a **post-consistency** summary (passing candidates only,
  with Pcons, Frontier, InBand, and Selected flags). The two views
  together make the rule's filter visible end-to-end. Pcons is shown
  only for passing candidates (the consistency evaluator returns NULL
  for non-passing ones). Default: FALSE.

- max_print:

  Integer. Maximum number of candidate rows printed in each
  `show_candidate_summary` table (forwarded to the internal
  `print_candidate_summary()`). Candidates are sorted selection-first,
  so the selected subgroup is always shown; any remaining candidates
  beyond `max_print` are collapsed into a trailing "... N more
  candidate(s) not shown" line. Affects console output only, not the
  returned object. Default: 10.

- pconsistency.digits:

  Integer. Number of decimal places to which the consistency proportion
  is rounded **before** it is compared with `pconsistency.threshold`.
  The rounded value decides admission and is the value the selection
  sort orders on, so this setting affects which subgroups are admitted,
  not only what is displayed. At the default of 2, a threshold of `p*`
  admits any candidate whose unrounded proportion is at or above the
  largest 0.01 grid point strictly below `p*` plus half a step (0.005) –
  so `p* = 0.99` admits from about 0.985. Raising the setting makes
  admission finer-grained and produces fewer exact ties for the
  selection rule to break. Under `consistency_method = "resample"` the
  proportion is a continuous closed-form quantity, so the rounding is a
  deliberate coarsening rather than a formatting step. Default: 2

- seed:

  Integer. Random seed for reproducible consistency splits.
  Default: 8316951. Set to NULL for non-reproducible random splits. The
  seed is used both for sequential execution (via set.seed()) and
  parallel execution (via future.seed).

- checking:

  Logical. Enable additional validation checks. Default: FALSE

- use_twostage:

  Logical. Use two-stage adaptive algorithm. Default: FALSE

- twostage_args:

  List. Parameters for two-stage algorithm:

  n.splits.screen

  :   Splits for Stage 1 screening. Default: 30

  screen.threshold

  :   Consistency threshold for Stage 1. Default: auto

  batch.size

  :   Splits per batch in Stage 2. Default: 20

  conf.level

  :   Confidence level for early stopping. Default: 0.95

  min.valid.screen

  :   Minimum valid Stage 1 splits. Default: 10

- parallel_args:

  List. Parallel processing configuration:

  plan

  :   Future plan: "multisession", "multicore", "callr", or
      "sequential". A sequential plan runs the consistency screen as a
      plain in-process loop regardless of `workers`: no future backend
      is set up and no per-batch globals resolution is paid. A parallel
      plan with `workers = 1` is honoured as parallel (one-worker
      batched backend); pass `plan = "sequential"` to run the plain
      loop. Results under `consistency_method = "resample"` are
      identical on every path (the closed form consumes no RNG). Under
      `"split"` the plain loop draws from a single RNG stream seeded
      once, whereas the batched backend gives each candidate its own
      L'Ecuyer-CMRG stream indexed by global position, so a sequential
      plan's `p.consistency` values differ from a parallel plan's within
      Monte-Carlo noise at the stated `fs.splits`; invariance to worker
      count holds within any plan.

  workers

  :   Number of parallel workers on a parallel plan.

  batch_size

  :   Number of candidates to evaluate per batch on a genuinely parallel
      plan – the knob for the per-batch overhead described above.
      Smaller values provide finer-grained early stopping but may
      increase overhead. Default: When stop_threshold is set and
      sg_focus is "hr" or "minSG", defaults to 1 (stop immediately when
      first candidate passes). For other sg_focus values with
      stop_threshold, defaults to min(workers, n_candidates/4). When
      stop_threshold is NULL, defaults to workers\*2 for efficiency.

  show_message

  :   Print parallel config messages

- estimator_fn:

  Closure or `NULL`. Effect-estimator closure from
  [`make_effect_estimator`](https://larry-leon.github.io/forestsearch/reference/make_effect_estimator.md)
  for GLM outcomes. Default `NULL` (uses Cox model).

- consistency_threshold:

  Numeric or `NULL`. Threshold for GLM consistency evaluation. Default
  `NULL` (uses `hr.consistency`).

- effect_label:

  Character. Column label for the effect measure in diagnostic output
  (e.g., candidate subgroup table). Default `"HR"`. Set automatically by
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  to the resolved `effect_measure` (e.g., `"IRR"`, `"OR"`) for GLM
  outcomes.

- effect_log_scale:

  Logical. If `TRUE`, the `HR` column in `hr.subgroups` stores log-scale
  values (log-OR, log-IRR) and will be exponentiated for display in
  diagnostic output. Default `FALSE` (Cox stores natural-scale HR).

- adjust_covariates:

  Character vector or `NULL`. Additional Cox model terms used when
  scoring survival subgroups during consistency evaluation, e.g.
  `"strata(x1)"` for a stratified model or `"x2"` for linear adjustment.
  Referenced columns must be present in `df`. Ignored on the GLM path.
  Default `NULL` (treatment-only, unadjusted).

- consistency_method:

  Character. `"resample"` (default) for the multiplier
  (influence-function) approximation via
  [`consistency_resample`](https://larry-leon.github.io/forestsearch/reference/consistency_resample.md),
  or `"split"` – the retained fallback – for the literal repeated 50/50
  split-and-refit calculation, which costs `fs.splits * 2` refits per
  candidate evaluation. On the GLM path `"resample"` requires a
  `glm_resample_spec`; without one (or when `consistency_resample` is
  unavailable) it falls back to `"split"` with a warning.
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  builds that spec automatically; a direct caller on a GLM outcome must
  supply it.

- glm_resample_spec:

  List or `NULL`. GLM resampling specification threaded from
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  when `consistency_method = "resample"` on a GLM outcome: components
  `outcome_type`, `effect_measure`, `treat.name`, `outcome.name`,
  `offset.name`, `adjust_covariates`, `adverse_outcome`, and
  `comparison_threshold` (the consistency threshold already on the
  comparison scale, i.e. `log(threshold)` for ratio measures). `NULL` on
  the survival path.

## Value

A list containing:

- out_sg:

  Selected subgroup results. When non-`NULL`, contains `result` (sorted
  candidate table; top row is the selected subgroup) and
  `pareto_frontier` (data.table of non-dominated candidates on (effect,
  N), both maximized – a post-hoc diagnostic, not used for selection).

- sg_focus:

  Selection criterion used

- df_flag:

  Data frame with treatment recommendations

- sg.harm:

  **Subgroup definition labels** – character vector of factor-level cut
  names (e.g., `c("z1.1", "z2.1")`), of length equal to the number of
  cuts defining the subgroup. `NULL` if no subgroup was identified.

- sg.harm.id:

  **Per-subject subgroup-membership indicator** – integer vector of
  length `nrow(df)` with `1` if subject \\i\\ is in the identified
  subgroup and `0` otherwise. `NULL` if no subgroup was identified.
  **Not** a character vector of cut expressions; see the *Field naming
  collision* section below.

- algorithm:

  "twostage" or "fixed"

- n_candidates_evaluated:

  Number of candidates actually evaluated

- n_candidates_total:

  Total candidates available

- n_passed:

  Number meeting consistency threshold

- early_stop_triggered:

  Logical indicating if early stop occurred

- early_stop_candidate:

  Index of candidate triggering early stop

- stop_threshold:

  Threshold used for early stopping

- seed:

  Random seed used for reproducibility (NULL if not set)

## Field naming collision with GRF results

The field name `sg.harm.id` has **different semantics** on this object
versus on the result objects returned by
[`grf.subg.harm.glm`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.glm.md)
and
[`grf.subg.harm.survival`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.survival.md):

|  |  |  |
|----|----|----|
| **Object** | **`sg.harm.id` contains** | **Length / type** |
| `subgroup.consistency()` result (this function) | per-subject 0/1 membership indicator | integer, length `nrow(df)` |
| [`grf.subg.harm.glm()`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.glm.md) result | character vector of cut expressions | character, length = depth of selected tree |
| [`grf.subg.harm.survival()`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.survival.md) result | character vector of cut expressions | character, length = depth of selected tree |

**Practical consequence.** Do not paste `sg.harm.id` with
`paste(..., collapse = " & ")` to print "the identified subgroup"
without first checking object class. For subgroup labels from a
[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
or `subgroup.consistency()` result, use `sg.harm` (character vector of
cut names) – the FS main result exposes `sg.harm` at the top level.

This naming collision is a documented CRAN-stable API for v0.1.x and
v0.2.x; it is expected to be resolved via a deprecation cycle in a
future minor release.

## Examples

``` r
if (FALSE) { # \dontrun{
# Standard evaluation
result <- subgroup.consistency(
  df = trial_data,
  hr.subgroups = candidates,
  sg_focus = "hr",
  n.splits = 400,
  parallel_args = list(plan = "multisession", workers = 6)
)

# Show top 10 candidates before evaluation
result <- subgroup.consistency(
  df = trial_data,
  hr.subgroups = candidates,
  sg_focus = "hr",
  show_candidate_summary = TRUE,  # Post-consistency summary
  n.splits = 400
)

# With early stopping and custom batch size
result <- subgroup.consistency(
  df = trial_data,
  hr.subgroups = candidates,
  sg_focus = "hr",
  stop_threshold = 0.95,
  show_candidate_summary = TRUE,
  parallel_args = list(
    plan = "multisession",
    workers = 6,
    batch_size = 2  # Check early stopping after every 2 candidates
  )
)
} # }
```
