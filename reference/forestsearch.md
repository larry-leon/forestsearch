# ForestSearch: Exploratory Subgroup Identification

Implements advanced statistical methods for exploratory subgroup
identification in clinical trials. Provides tools for identifying
patient subgroups with differential treatment effects using machine
learning approaches including Generalized Random Forests ('GRF'), LASSO
regularization, and exhaustive combinatorial search algorithms. Supports
survival endpoints (Cox proportional hazards), binary outcomes (log odds
ratio, log relative risk, risk difference), continuous outcomes (mean
difference), and count / rate outcomes (log incidence rate ratio via
Poisson, quasi-Poisson, or negative-binomial GLMs with optional
person-time offset). Features bootstrap bias correction using
infinitesimal jackknife methods to address selection bias in post-hoc
analyses. Designed for clinical researchers conducting exploratory
subgroup analyses in randomized controlled trials, particularly for
multi-regional clinical trials ('MRCT') requiring regional consistency
evaluation. Methods are described in Leon et al. (2024)
[doi:10.1002/sim.10163](https://doi.org/10.1002/sim.10163) .

## Usage

``` r
forestsearch(
  df.analysis,
  outcome.name = "tte",
  event.name = "event",
  treat.name = "treat",
  id.name = "id",
  potentialOutcome.name = NULL,
  flag_harm.name = NULL,
  confounders.name = NULL,
  parallel_args = list(plan = "multisession", workers = .default_parallel_workers(),
    show_message = TRUE),
  df.predict = NULL,
  df.test = NULL,
  is.RCT = TRUE,
  seedit = 8316951,
  est.scale = "hr",
  use_lasso = FALSE,
  use_grf = FALSE,
  grf_res = NULL,
  grf_cuts = NULL,
  use_dina = FALSE,
  dina_res = NULL,
  dina_cuts = NULL,
  dina_args = list(),
  dina_select_statistic = c("effect", "dina"),
  subgroup_method = c("consistency", "dina", "grf"),
  max_n_confounders = 1000,
  grf_depth = 2,
  grf_selection = c("frontier", "tree"),
  grf_select_statistic = c("effect", "dr"),
  dmin.grf = 0,
  frac.tau = 0.8,
  return_selected_cuts_only = TRUE,
  conf_force = NULL,
  defaultcut_names = NULL,
  cut_type = "default",
  exclude_cuts = NULL,
  replace_med_grf = FALSE,
  cont.cutoff = 4,
  conf.cont_medians = NULL,
  conf.cont_medians_force = NULL,
  conf.cont_jcuts = NULL,
  collapse_cuts = TRUE,
  collapse_cuts_args = list(),
  n.min = 60,
  n.min.frac = 0.1,
  effect.threshold = NULL,
  consistency.threshold = NULL,
  hr.threshold = 1.25,
  hr.consistency = 1,
  sg_focus = "hr",
  selection_rule = "neighborhood",
  effect_neighborhood = 0.1,
  fs.splits = 1000,
  m1.threshold = Inf,
  pconsistency.threshold = 0.9,
  stop_threshold = pconsistency.threshold,
  pconsistency.digits = 2,
  show_candidate_summary = FALSE,
  max_print = 10,
  d0.min = 10,
  d1.min = 10,
  max.minutes = 3,
  minp = 0.025,
  details = FALSE,
  quiet = FALSE,
  maxk = 2,
  by.risk = 12,
  plot.sg = FALSE,
  plot.grf = FALSE,
  max_subgroups_search = Inf,
  vi.grf.min = NULL,
  use_twostage = TRUE,
  twostage_args = list(),
  consistency_method = c("resample", "split"),
  outcome_type = c("survival", "binary", "continuous", "count"),
  effect_measure = NULL,
  offset.name = NULL,
  adverse_outcome = NULL,
  overdispersion = c("none", "quasi", "negbin"),
  grf_count_transform = c("log", "identity"),
  tune_grf = FALSE,
  adjust_covariates = NULL,
  ps_method = NULL,
  ps_adjust_method = c("none", "iptw", "dr_gcomp"),
  ps_hat = NULL,
  mr_inference = FALSE,
  mr_inference_args = list()
)
```

## Arguments

- df.analysis:

  Data frame. Analysis dataset with required columns.

- outcome.name:

  Character. Name of time-to-event outcome variable. Default "tte".

- event.name:

  Character. Name of event indicator (1=event, 0=censored). Default
  "event".

- treat.name:

  Character. Name of treatment variable (1=treatment, 0=control).
  Default "treat".

- id.name:

  Character. Name of subject ID variable. Default "id".

- potentialOutcome.name:

  Character. Name of potential outcome variable (optional).

- flag_harm.name:

  Character. Name of true harm flag for simulations (optional).

- confounders.name:

  Character vector. Names of candidate subgroup-defining variables.

- parallel_args:

  List. Parallel processing configuration:

  plan

  :   Character. One of `"multisession"` (default; cross-platform,
      reuses workers across calls), `"multicore"` (Linux/macOS only,
      fastest startup but fails in RStudio interactive sessions),
      `"callr"` (slower per-call, fully isolated R processes), or
      `"sequential"` (no parallelism).

  workers

  :   Integer. Number of parallel workers. Default: `floor(0.80 * N)`
      where `N` is the number of physical (non-hyperthreaded) cores,
      capped to 2 when `_R_CHECK_LIMIT_CORES_` is set (e.g., during
      `R CMD check`). On a 4-core laptop the default is 3 workers; on a
      16-core workstation it is 12.

  show_message

  :   Logical. Show parallel setup messages.

  batch_size

  :   Integer. Number of *candidate subgroups* dispatched to the workers
      per batch during consistency evaluation. Default: chosen from the
      worker count and pool size. **This is a performance knob only – it
      does not affect the selected subgroup.** Smaller batches let an
      early stop take effect sooner (relevant only for
      `sg_focus = "maxeffCons"`, the one focus that early-stops); larger
      batches amortise dispatch overhead. Not to be confused with
      `twostage_args$batch.size`, which batches *splits* within Stage 2
      of the two-stage consistency algorithm – a different quantity in
      different units.

  To disable parallelism entirely, pass
  `parallel_args = list(plan = "sequential")`.

- df.predict:

  Data frame. Prediction dataset (optional).

- df.test:

  Data frame. Test dataset (optional).

- is.RCT:

  Logical. Is this a randomized controlled trial? Default TRUE.

- seedit:

  Integer. Random seed. Default 8316951.

- est.scale:

  Character. Estimation scale ("hr" or "rmst"). Default "hr".

- use_lasso:

  Logical. Use LASSO for variable selection. Default FALSE. The default
  changed from `TRUE` in an earlier release: with
  `subgroup_method = "consistency"` a model-based front end makes the
  candidate family data-dependent, which multiplier resampling
  (`mr_inference = TRUE`) requires be fixed. Set `TRUE` to restore the
  prognostic-lasso prefilter used in Leon et al. (2024).

- use_grf:

  Logical. Generate additional candidate cuts from a GRF fit:
  tree-derived cutpoints are extracted via
  [`grf.subg.harm.survival()`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.survival.md)
  /
  [`grf.subg.harm.glm()`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.glm.md)
  (or taken from `grf_res` / `grf_cuts` when supplied) and enter the
  candidate family through
  [`get_FSdata`](https://larry-leon.github.io/forestsearch/reference/get_FSdata.md)
  alongside the quartile (and any DINA) cuts. This flag governs
  candidate-cut generation only: GRF variable-importance ordering of the
  cut columns is a separate step controlled solely by `vi.grf.min` and
  runs, or is skipped, regardless of `use_grf`. Default FALSE. The
  default changed from `TRUE` in an earlier release, for the same
  fixed-candidate-family reason given under `use_lasso`.

- grf_res:

  GRF results object (optional, for reuse).

- grf_cuts:

  List. Custom GRF cut points (optional).

- use_dina:

  Logical. Generate additional screening-stage candidate cuts from a
  DINA per-covariate Pareto frontier, fed into
  [`get_FSdata`](https://larry-leon.github.io/forestsearch/reference/get_FSdata.md)
  alongside the GRF cuts. Default FALSE (opt-in); existing behaviour is
  unchanged when FALSE.

- dina_res:

  A fitted DINA object (optional, for reuse). When NULL and
  `use_dina = TRUE`, a DINA model is fit at the screening stage.

- dina_cuts:

  Character vector of pre-supplied DINA cut expressions (optional). When
  supplied, the DINA fit / frontier step is skipped and these are used
  directly.

- dina_args:

  Named list of DINA tuning options, all optional. Recognised keys:
  `family` (defaults to the mapping from `outcome_type`: survival -\>
  cox, binary -\> binomial, continuous -\> gaussian, count -\> poisson),
  `seed` (defaults to `seedit`), and the fit-only keys `n_folds`,
  `cens_type`, `cens_params` (passed to
  [`dina`](https://larry-leon.github.io/forestsearch/reference/dina.md)
  only when set); plus the frontier keys forwarded to
  [`dina_frontier`](https://larry-leon.github.io/forestsearch/reference/dina_frontier.md):
  `scope` (default `"wide"`), `m_diff`, `n_min` (defaults to the `n.min`
  of this call), `direction`, `max_per_covariate`, `max_subgroups`, and
  `digits`; plus the screening-behavior keys `selected_only` (default
  `TRUE`, matching GRF's `return_selected_cuts_only`), `max_depth` (`2L`
  default – aligning DINA with `maxk = 2` and `grf_depth = 2` so the
  selected cut may be an AND-conjunction of two covariates; set `1L` for
  single-covariate cuts only), and `grid_probs` (the per-covariate
  quantile grid for depth-2 pair thresholds; default interior deciles).
  `max_depth` / `grid_probs` are forwarded to
  [`dina_subgroup`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md)
  and apply to both the selected-cut screening and
  `subgroup_method = "dina"`. (The selector statistic formerly accepted
  here as `select_statistic` is now the top-level
  `dina_select_statistic` argument; a value supplied inside `dina_args`
  is deprecated and overridden with a warning.) When
  `selected_only = TRUE` (the default), `use_dina` screening contributes
  the cut(s) chosen by
  [`dina_subgroup`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md)
  – using this call's `sg_focus` / `selection_rule` /
  `effect_neighborhood` / `n.min` / `hr.threshold`, i.e. the same cut
  `subgroup_method = "dina"` would select. Set `FALSE` to contribute the
  full frontier candidate set instead.

  **`m_diff` is the candidacy floor, on the link scale, and is derived
  rather than taken from this list.** It is the harm floor a candidate
  must clear to enter the qualifying family that `sg_focus` then ranks –
  `log(HR)` for `cox`, `log(OR)` for `binomial`, `log(IRR)` for
  `poisson`, and the raw mean difference for `gaussian`. Under
  `subgroup_method = "dina"` (and under `selected_only = TRUE`
  screening) it is computed as `log(hr.threshold)` – identity for
  `gaussian` – from the resolved `effect.threshold`, and an `m_diff`
  supplied here is **ignored**. Set the floor with `effect.threshold` /
  `hr.threshold`. The key remains meaningful only where the raw
  [`dina_frontier`](https://larry-leon.github.io/forestsearch/reference/dina_frontier.md)
  candidate set is used directly.

  **The same holds for every frontier key, not only `m_diff`.** All
  seven – `scope`, `m_diff`, `n_min`, `direction`, `max_per_covariate`,
  `max_subgroups` and `digits` – are **inert under**
  `subgroup_method = "dina"`, which selects with
  [`dina_subgroup`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md)
  and never consults a frontier table. Under `use_dina = TRUE` screening
  with `selected_only = TRUE` (the default) a frontier is computed and
  then discarded in favour of the selected cut, so six of the seven
  change nothing; only `digits` still acts, rounding the threshold in
  the emitted cut expression. The seven keys govern the candidate pool
  only under `selected_only = FALSE`, where
  [`dina_frontier`](https://larry-leon.github.io/forestsearch/reference/dina_frontier.md)'s
  `cut_expr` column *is* the pool. Supplying any of the seven under
  `subgroup_method = "dina"` now raises a single
  [`warning()`](https://rdrr.io/r/base/warning.html) naming every
  offending key.

  Under `dina_select_statistic = "effect"` the floor still defines the
  qualifying family on the tau-hat scale while selection re-ranks that
  family on the refit inferential effect; because multiplier resampling
  de-biases over the realized selection event, `m_diff` fixes the family
  MR de-biases over. See
  [`dina_subgroup`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md).
  Unknown keys raise an error.

- dina_select_statistic:

  Character, one of `"effect"` (default) or `"dina"`, used only when
  `subgroup_method = "dina"`. The DINA analogue of
  `grf_select_statistic`. `"dina"` ranks DINA's qualifying family on its
  native subgroup-mean tau-hat (unchanged legacy behaviour). `"effect"`
  re-ranks the same family on the inferential effect measure (Cox HR for
  survival; the resolved GLM effect otherwise), scored with the same
  per-candidate estimator multiplier resampling (MR) uses, so the
  realized selection is the exact event MR de-biases (and MR's
  native-family restriction follows the effect rather than tau-hat).
  `"effect"` is the default because only it ranks on the inferential
  coefficient; `"dina"` is rejected under `mr_inference = TRUE`.

- subgroup_method:

  Character, one of `"consistency"` (default), `"dina"`, or `"grf"`.
  `"consistency"` is the standard ForestSearch pipeline (GRF/LASSO
  screening then the consistency search). `"dina"` delegates subgroup
  *selection* to
  [`dina_subgroup`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md)
  and bypasses GRF, LASSO and the consistency search: a DINA model is
  fit, a single subgroup is selected using this call's `sg_focus` /
  `selection_rule` / `effect_neighborhood` / `n.min` / `hr.threshold`,
  and that subgroup is returned as `sg.harm`. `dina_args` supplies only
  the DINA *fit* tuning in this mode (frontier keys are ignored).
  `"grf"` is the GRF analogue: it delegates selection to the GRF
  subgroup-identification routine
  ([`grf.subg.harm.survival`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.survival.md)
  for survival,
  [`grf.subg.harm.glm`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.glm.md)
  for GLM outcomes), again bypassing LASSO, DINA and the consistency
  search, and returns the policy-tree harm leaf as `sg.harm`.
  GRF-selection mode is configured from this call's existing GRF
  arguments (`frac.tau`, `dmin.grf`, `grf_depth`, `is.RCT`, `seedit`,
  and the GLM extras) – there is no separate argument list – so it uses
  the same settings GRF screening would. For both `"dina"` and `"grf"`,
  the existing bootstrap de-biasing applies unchanged, re-fitting and
  re-selecting on each replicate. The selected subgroup may be a single
  covariate or, where the method's depth allows (DINA `max_depth = 2`;
  GRF `grf_depth >= 2`), an AND-conjunction.

- max_n_confounders:

  Integer. Maximum number of candidate cut columns retained after GRF
  variable-importance screening, keeping the highest-importance columns.
  Default 1000. **Applied only when `vi.grf.min` is numeric**: the
  truncation lives inside the variable-importance block, so with
  `vi.grf.min = NULL` (the default) this argument has no effect, and a
  finite value below the number of cut columns triggers a warning saying
  so.

- grf_depth:

  Integer. GRF tree depth. Default 2.

- grf_selection:

  Character, one of `"frontier"` (default) or `"tree"`, used only when
  `subgroup_method = "grf"`. `"tree"` selects via the GRF policy tree.
  `"frontier"` selects from the Pareto frontier of the doubly-robust
  scores (see
  [`grf.subg.harm.survival`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.survival.md));
  the frontier rule is taken from `sg_focus` – `"eff"`, `"effMaxSG"`,
  `"effMinSG"`, `"maxSG"` and `"minSG"` each map to the aligned frontier
  rule, as in
  [`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md)
  – and the relative band from `effect_neighborhood`. `"maxeff"` and
  `"maxeffCons"` map **explicitly** to the plain effect-argmax rule
  (`"eff"`), the same rule `"hr"` takes: GRF computes no `Pcons`, so the
  `Cons` qualifier has nothing to bind to. The mapping is exhaustive
  over the canonical foci and raises an error rather than falling
  through, so no focus can silently select under a rule the caller did
  not request. (Earlier versions had no branch for these two and
  defaulted them to `"effMaxSG"` silently; pairing them with
  `grf_selection = "frontier"` is now well defined.) `"frontier"` is the
  default because it enumerates a family that can be ranked on the
  inferential effect; the policy tree exposes no such family, so
  `"tree"` is rejected under `mr_inference = TRUE`.

- grf_select_statistic:

  Character, one of `"effect"` (default) or `"dr"`, used only when
  `subgroup_method = "grf"` and `grf_selection = "frontier"`. `"dr"`
  selects the frontier winner on GRF's native doubly-robust score
  (legacy behaviour). `"effect"` re-ranks the same DR-candidate frontier
  on the inferential effect measure (Cox HR for survival; the resolved
  GLM effect otherwise), scored with the same per-candidate estimator
  multiplier resampling (MR) uses, so the realized selection is the
  exact event MR de-biases (and MR's native-family restriction follows
  the effect rather than the DR score). `"effect"` is the default
  because only it ranks on the inferential coefficient; `"dr"` is
  rejected under `mr_inference = TRUE`. Ignored for
  `grf_selection = "tree"` (the policy-tree leaf is the selection, with
  no enumerated family to rank).

- dmin.grf:

  Numeric. Minimum events for GRF. Default 0.0.

- frac.tau:

  Numeric in (0, 1\]. Multiplier on the GRF time horizon passed to
  [`grf::causal_survival_forest()`](https://rdrr.io/pkg/grf/man/causal_survival_forest.html).
  The effective horizon is
  `frac.tau * min(max(Y[treated, event]), max(Y[control, event]))`, i.e.
  `frac.tau` times the smaller of the two arm-specific maximum event
  times. Values \< 1 truncate the horizon inward of that maximum (the
  conservative choice when events are sparse near the tail). Applies to
  both the GRF variable-importance pre-screen inside `forestsearch()`
  and the downstream call to
  [`grf.subg.harm.survival`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.survival.md).
  Used only when `outcome_type = "survival"`; ignored for GLM outcomes
  (binary, continuous, count). Default 0.8.

- return_selected_cuts_only:

  Logical. If TRUE (default), GRF returns only cuts from the tree depth
  that identified the selected subgroup meeting `dmin.grf`. If FALSE
  returns all cuts from all fitted trees (depths 1 through `grf_depth`).
  See
  [`grf.subg.harm.survival`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.survival.md)
  for details.

- conf_force:

  Character vector of forced cut expressions (optional), e.g.
  `c("er <= 0", "pgr <= 10")`. Literal cut locations: unlike the
  `conf.cont_*` arguments these do not move with a resample.

- defaultcut_names:

  Character vector. Default cut variable names (optional).

- cut_type:

  Character. Cut strategy: `"default"` or `"median"`. Default
  `"default"`. These are the only permitted values –
  [`get_FSdata()`](https://larry-leon.github.io/forestsearch/reference/get_FSdata.md)
  validates against them and errors otherwise.

- exclude_cuts:

  Character vector. Variables to exclude from cutting (optional).

- replace_med_grf:

  Logical. Replace median with GRF cuts. Default FALSE.

- cont.cutoff:

  Integer. Cutoff for continuous vs categorical. Default 4.

- conf.cont_medians:

  Character vector of continuous confounders to cut at their median
  (optional). **Names, not values** – the cut location is the median
  computed from the supplied data, so it is not fixed by this argument
  and moves with a resample. To pin a cut location, pass a literal cut
  expression via `conf_force`.

- conf.cont_medians_force:

  Character vector of additional continuous confounders to force a
  median cut on (optional). Names, not values, with the same consequence
  as `conf.cont_medians` above.

- conf.cont_jcuts:

  Named integer list (each value \>= 1). Per-variable J-quantile cut
  override for continuous covariates: for an entry `X = J`, the default
  cut set for X is replaced by `J` binary cut points at the (k/(J+1))-th
  quantiles of X (k = 1, ..., J), defining J+1 non-overlapping
  intervals. Names must be in `confounders.name` and must not overlap
  with `defaultcut_names` or `conf.cont_medians_force`. J-cut variables
  bypass LASSO filtering, matching `defaultcut_names` semantics.
  Variables not listed retain default behaviour. Forwarded to
  [`get_FSdata`](https://larry-leon.github.io/forestsearch/reference/get_FSdata.md)
  via
  [`filter_call_args()`](https://larry-leon.github.io/forestsearch/reference/filter_call_args.md);
  see there for full semantics.

- collapse_cuts:

  Logical. If TRUE, collapse near-redundant continuous candidate cuts
  before the search. Continuous covariates can generate many thresholds
  that are practically redundant (e.g. `"age <= 35.7"` and
  `"age <= 35"`); when enabled, cuts on the same variable and operator
  within a per-variable standard-error band are merged to a single
  rounded-centroid threshold, subject to a membership safety check.
  Forwarded to
  [`get_FSdata`](https://larry-leon.github.io/forestsearch/reference/get_FSdata.md);
  see
  [`collapse_redundant_cuts`](https://larry-leon.github.io/forestsearch/reference/collapse_redundant_cuts.md)
  for the algorithm. Default TRUE. Pass `collapse_cuts = FALSE` to
  recover the un-coarsened candidate pool used by earlier package
  versions.

- collapse_cuts_args:

  List of overrides merged onto the defaults
  `list(c = 1.0, tol = 0.05, digits = 0L)`: band multiplier `c`
  (`band = c * sd(x)/sqrt(n)`), membership safety tolerance `tol`
  (fraction of n when \< 1, absolute count when \>= 1), and
  representative rounding `digits`. Ignored when
  `collapse_cuts = FALSE`.

- n.min:

  Integer or `NULL`. Minimum subgroup size. Default 60. Supplying a
  value (or omitting the argument) uses that fixed floor, as in prior
  versions. Passing `n.min = NULL` opts into a sample-size-adaptive
  floor `max(60, ceiling(n.min.frac * N))`, where `N` is the analysis
  sample size – at least 60, and at least `n.min.frac` of `N`. Growing
  the minimum subgroup size with `N` curbs false discovery from very
  small subgroups in large samples. The resolved integer is used
  everywhere downstream (screening, DINA, GRF, the consistency search,
  and the bootstrap).

- n.min.frac:

  Numeric in (0, 1). Fraction of the analysis sample size used for the
  adaptive `n.min` floor when `n.min = NULL`. Default 0.10. Ignored when
  `n.min` is supplied.

- effect.threshold:

  Numeric or NULL. `c1`, the **screening threshold on a candidate
  subgroup's own effect**, estimated on the full sample. A candidate
  must clear it to enter the consistency stage. For ratio-scale measures
  (OR, RR, IRR, HR): on the ratio scale (e.g., 1.5 means OR \>= 1.5).
  For identity-scale measures (RD, IRD, MD): on the identity scale
  (e.g., 0.07 means RD \>= 7 pct points). If NULL (default), falls back
  to `hr.threshold`. Distinct from `consistency.threshold` (`c2`), which
  is applied within each split half, and from `pconsistency.threshold`
  (`p*`), which is a proportion rather than an effect. On the
  consistency path with a ratio estimand (survival HR, binary `"OR"`),
  supplying `c1` and **not** `c2` derives `c2 = 0.80 * c1` on the ratio
  scale, with a message; see the *Threshold vocabulary and resolved
  defaults* section for that rule and for the value this resolves to
  under each estimand.

- consistency.threshold:

  Numeric or NULL. `c2`, **an effect threshold** – **not** a proportion,
  and **not** the consistency rate. It is applied to the subgroup's
  effect estimated *within each split half*; a split counts as
  consistent when both of its halves clear `c2`. Same scale conventions
  as `effect.threshold`. If NULL (default), falls back to
  `hr.consistency`. Its companion is `pconsistency.threshold` (`p*`),
  the *fraction* of splits that must clear `c2`: this argument sets the
  bar, that one counts how often the bar is met. Both names contain the
  word "consistency" and they are different quantities. On the
  consistency path with a ratio estimand (survival HR, binary `"OR"`)
  `c2` must be **less than or equal to** `c1`; `c2 > c1` is an error.
  Left unsupplied there, it is derived as `0.80 * c1` rather than left
  at its default. See the *Threshold vocabulary and resolved defaults*
  section for both rules and for the value this resolves to under each
  estimand.

- hr.threshold:

  Numeric. Legacy name for `effect.threshold` – `c1`, the screening
  threshold on a candidate subgroup's own effect. Retained for backward
  compatibility. Default 1.25. When `effect.threshold` is provided,
  `hr.threshold` is ignored; supplying both at disagreeing values is an
  error. Carries the same `c2` derivation as `effect.threshold`.

- hr.consistency:

  Numeric. Legacy name for `consistency.threshold` – `c2`, the per-split
  **effect** threshold, not a proportion and not the consistency rate.
  Retained for backward compatibility. Default 1.0. When
  `consistency.threshold` is provided, `hr.consistency` is ignored;
  supplying both at disagreeing values is an error. Carries the same
  `c2 <= c1` requirement as `consistency.threshold`.

- sg_focus:

  Character. Subgroup selection focus – *what* is selected, not merely
  how candidates are sorted. Except for `"maxeff"`, every focus selects
  among **qualifiers**: candidates whose consistency rate reaches
  `pconsistency.threshold`. One of:

  `"hr"` (or `"eff"`, `"maxcons"`)

  :   Pick the qualifier with the highest consistency (`Pcons`); ties
      broken by largest effect. **This is consistency-primary with
      effect only as a tiebreak** – despite the name, it does not
      maximise the effect. Use `"maxeffCons"` for that. `"maxcons"` is
      the recommended spelling, since it says what the rule does.

  `"maxeffCons"`

  :   Pick the qualifier with the largest effect; ties broken by fewest
      cuts (`K`). The effect maximiser *subject to* the consistency
      floor – as against `"maxeff"`, which applies no floor. This is the
      only focus for which early stopping is sound; see
      `stop_threshold`.

  `"maxSG"`

  :   Pick the qualifier with the largest sample size; ties broken by
      consistency. Not a spelling of `"hrMaxSG"`: it applies no effect
      band.

  `"minSG"`

  :   Pick the qualifier with the smallest sample size; ties broken by
      consistency. Not a spelling of `"hrMinSG"`.

  `"hrMaxSG"` (or `"effMaxSG"`)

  :   Among candidates with effect size within `effect_neighborhood` of
      the maximum, pick the one with the *largest* sample size.

  `"hrMinSG"` (or `"effMinSG"`)

  :   Among candidates with effect size within `effect_neighborhood` of
      the maximum, pick the one with the *smallest* sample size.

  `"maxeff"`

  :   The effect maximiser over the enumerated family, with no
      consistency threshold, effect threshold, early stopping, two-stage
      screening or candidate-pool truncation. Implements the argmax
      primitive of Guo and He (2021). Size/event minima (`n.min`,
      `d0.min`, `d1.min`) still apply as estimability constraints;
      prevalence and redundancy pruning are relaxed to their most
      permissive setting (`minp` -\> 0, `rmin` -\> 0) and duplicate
      removal collapses identical-membership candidates to their
      fewest-cut representative, so only degenerate or
      exactly-duplicated candidates are dropped and the argmax is
      invariant. Because consistency (`Pcons`) never enters the argmax,
      it is computed for the SELECTED subgroup only (reported alongside
      it) rather than for every candidate, avoiding the per-candidate
      consistency refits. Contrast `"maxeffCons"`, which applies the
      same argmax but only among candidates clearing
      `pconsistency.threshold`, and `"hr"`, which ranks by consistency
      first and uses effect only as a tiebreaker.

  **The descriptions above are for `subgroup_method = "consistency"`.**
  There, `"maxeffCons"` retains *both* the consistency floor
  (`pconsistency.threshold`) and the effect floor, while `"maxeff"`
  disables both – the distinction the `Cons` suffix names.

  For `subgroup_method = "dina"` and `"grf"`, `"maxeff"` and
  `"maxeffCons"` are **synonyms**, both selecting the effect maximiser
  over the identifier's qualifiers with the harm floor active. Neither
  engine has a consistency floor – neither computes `Pcons`, and no DINA
  or GRF sort key contains one – so the `Cons` qualifier has nothing to
  bind to. The two could in principle have been separated by the effect
  floor alone; they deliberately are not, because the manuscript keeps
  the harm threshold in a common form across all three identifiers, and
  `"maxeff"` exists for `"consistency"` as a comparison against the
  argmax primitive of Guo and He (2021), which has no DINA or GRF
  counterpart. No floor is disabled by this on any path, and the
  collapse raises no condition: it is intended behaviour, and a warning
  would fire on every such run.

  **Per-engine resolution.** Every accepted spelling, and the rule that
  actually runs for it on each engine:

  |  |  |  |
  |----|----|----|
  | **sg_focus** | **consistency** | **dina / grf** |
  | `eff`, `hr`, `maxcons` | consistency argmax (`Pcons`-primary, effect breaks ties) | effect argmax, `order(-eff)` |
  | `maxeff` | effect argmax, no floors | effect argmax, `order(-eff)` |
  | `maxeffCons` | effect argmax with the consistency floor | effect argmax, `order(-eff)` |
  | `effMaxSG`, `hrMaxSG` | largest N within the effect band, per `selection_rule` | as documented in [`dina_subgroup`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md) |
  | `effMinSG`, `hrMinSG` | smallest N within the effect band, per `selection_rule` | as documented in [`dina_subgroup`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md) |
  | `maxSG` | largest N among qualifiers | largest N |
  | `minSG` | smallest N among qualifiers | smallest N |

  The collapse sets, stated explicitly. On `"consistency"`, {`eff`,
  `hr`, `maxcons`} is one rule, and `maxeff` and `maxeffCons` are two
  further *distinct* rules. On `"dina"` and `"grf"`, {`eff`, `hr`,
  `maxcons`, `maxeff`, `maxeffCons`} is **one** rule – five spellings,
  one `order(-eff)`.

  `forestsearch()` announces both collapses at run time: when the
  canonical rule differs from the spelling passed, it emits one
  [`message()`](https://rdrr.io/r/base/message.html) (suppressed by
  `quiet`) naming the rule that will actually run.
  [`fs_focus_tag`](https://larry-leon.github.io/forestsearch/reference/fs_focus_tag.md)
  returns the stem tag for each cell of the table above, and is the
  supported way to label output files by the rule that ran rather than
  by the spelling requested.

  Default `"hr"`. The `"eff*"` forms are aliases for the `"hr*"` forms
  and read more naturally in GLM contexts (continuous MD, binary
  OR/RR/RD, count IRR) where there is no hazard ratio. Both vocabularies
  produce identical results; pick whichever fits the outcome type. See
  [`sort_subgroups`](https://larry-leon.github.io/forestsearch/reference/sort_subgroups.md)
  for full sort keys and tiebreakers.

  **Scope of the band arguments.** `selection_rule` and
  `effect_neighborhood` configure the effect band, which only
  `"hrMaxSG"` / `"hrMinSG"` (equivalently `"effMaxSG"` / `"effMinSG"`)
  consult. They are inert for every other focus, and supplying
  `selection_rule` other than `"neighborhood"` with a non-band focus is
  an error rather than a silent no-op.

  **Choosing between the three effect-oriented rules.** `"maxeff"`
  answers "what is the most extreme subgroup in this family?";
  `"maxeffCons"` answers "what is the most extreme subgroup that also
  replicates across splits?"; `"hr"`/`"maxcons"` answers "which subgroup
  replicates most reliably?". They can and do select different subgroups
  on the same data.

- selection_rule:

  Character. Rule defining the candidate inclusion set for `"hrMaxSG"` /
  `"hrMinSG"`. One of `"neighborhood"` (default; current behaviour),
  `"pareto"` (2D Pareto-non-dominated set in (effect, N)), or `"both"`
  (intersection of the two). Forwarded to
  [`sort_subgroups`](https://larry-leon.github.io/forestsearch/reference/sort_subgroups.md)
  via
  [`filter_call_args()`](https://larry-leon.github.io/forestsearch/reference/filter_call_args.md);
  see there for full semantics. Must be `"neighborhood"` for other
  `sg_focus` values.

- effect_neighborhood:

  Numeric in `[0, 1)`. Relative tolerance for `"hrMaxSG"` and
  `"hrMinSG"` selection. A candidate is in the *effect-size
  neighborhood* iff its (natural-scale) effect is at least
  `(1 - effect_neighborhood) * max(effect)`. Default `0.10` (within 10\\
  `effect_neighborhood = 0` reduces `"hrMaxSG"`/ `"hrMinSG"` to a strict
  max-effect filter (only candidates tied at the maximum qualify). For
  GLM ratio measures (OR, IRR), the neighborhood test is applied on the
  natural scale. Ignored when `sg_focus` is `"hr"`, `"maxSG"`, or
  `"minSG"`.

- fs.splits:

  Integer. Number of splits for consistency evaluation (or maximum
  splits when `use_twostage = TRUE`). Default 1000.

- m1.threshold:

  Numeric. Maximum median survival threshold. Default Inf.

- pconsistency.threshold:

  Numeric. `p*`, a **proportion** in `[0, 1]` – the fraction of random
  splits that must be consistent, where a split is consistent when both
  of its halves clear `consistency.threshold` (`c2`). It is a rate, not
  an effect threshold: it never enters an effect comparison, and unlike
  `c1` and `c2` it is **not** remapped per estimand – the same number
  means the same thing for every `outcome_type`. Default 0.90.
  Overridden to `0` under `sg_focus = "maxeff"` (no consistency filter).

- stop_threshold:

  Numeric in `[0, 1]` or `NULL`. Early stopping threshold for
  consistency evaluation. When a candidate subgroup's consistency
  probability (`Pcons`) meets or exceeds this threshold, evaluation
  stops early – remaining candidates are skipped. Set to `NULL` to
  disable early stopping and force evaluation of all candidates up to
  `max_subgroups_search`. Defaults to `pconsistency.threshold`, so it
  tracks the consistency floor automatically if that is changed.

  **Applies to `sg_focus = "maxeffCons"` only.** Early stopping halts
  the scan and then applies the focus's selection key to the resulting
  prefix, so it is sound only when that key is fully determined by the
  enumeration order. `Pcons` is not – it is unknown until a candidate is
  evaluated – so any focus whose key contains a `Pcons` term can be
  beaten by a candidate the scan never reached. `"maxeffCons"` is the
  only focus with no such term: its key `(-hr, K)` is exactly the
  enumeration order `(-HR, K)`, ties included, so the first qualifying
  candidate is provably the winner.

  For every other focus `stop_threshold` is reset to `NULL`, with a
  warning if it was supplied explicitly:

  - `"hr"` (`"eff"`, `"maxcons"`) selects on `Pcons` directly;

  - `"maxSG"` / `"minSG"` use `Pcons` to break ties in subgroup size,
    which a truncated scan can miss;

  - `"hrMaxSG"` / `"hrMinSG"` additionally need the whole pool to
    compute the effect neighborhood;

  - `"maxeff"` evaluates only its winner, so early stopping has nothing
    to skip.

  **Note:** Values \> 1.0 are not permitted. To disable early stopping,
  use `stop_threshold = NULL`, not a value above 1.

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
  deliberate coarsening rather than a formatting step. Passed to
  [`subgroup.consistency`](https://larry-leon.github.io/forestsearch/reference/subgroup.consistency.md).
  Default 2.

- show_candidate_summary:

  Logical. If `TRUE`, prints a post-consistency summary table of all
  passing candidates with Frontier/InBand/Selected flags. See
  [`subgroup.consistency`](https://larry-leon.github.io/forestsearch/reference/subgroup.consistency.md)
  for details. Default `FALSE`.

- max_print:

  Integer. Maximum number of candidate subgroups to print in the pre-
  and post-consistency candidate summaries when
  `show_candidate_summary = TRUE`. Candidates are printed in sorted
  order with the selected subgroup first, so `max_print = 10` shows the
  selected subgroup and the next nine, followed by a "... N more not
  shown" line. Only the per-candidate rows are capped; the banner,
  header, legend, and summary footer are always printed. Default `10`
  (the selected subgroup and the next nine); use `Inf` to print all
  candidates.

- d0.min:

  Integer. Minimum per-arm filter for candidate subgroups. For
  `"survival"` and `"binary"` outcomes, this is the minimum number of
  events in the control arm (for binary, events are Y = 1). Ignored for
  `"continuous"` outcomes (only `n.min` applies). Default 10.

- d1.min:

  Integer. Same as `d0.min` for the treatment arm. Default 10.

- max.minutes:

  Numeric. **Currently inert; scheduled for deprecation in v0.3.0.**
  Previously intended as a wall-clock time budget for the combination
  search, this argument is no longer enforced in the parallelized search
  path and has no effect on behavior. Search scope is governed by
  `maxk`, the candidate-factor front ends (`use_grf`, `use_lasso`,
  `conf_force`), and the number of parallel workers. The default of 3 is
  retained only for signature compatibility.

- minp:

  Numeric. Minimum prevalence threshold. Default 0.025.

- details:

  Logical. Print progress details. Default FALSE.

- quiet:

  Logical. If TRUE, suppress the configuration summary message printed
  at startup. Useful when `forestsearch()` is called repeatedly inside
  bootstrap or cross-validation loops. Default FALSE.

- maxk:

  Integer. Maximum number of factors per subgroup. Default 2.

- by.risk:

  Integer. Risk table interval. Default 12.

- plot.sg:

  Logical. Plot subgroup survival curves. Default FALSE.

- plot.grf:

  Logical. Plot GRF results. Default FALSE.

- max_subgroups_search:

  Numeric. Maximum subgroups to evaluate. **Default `Inf`** – no
  truncation. The cap has no counterpart in the method: the candidate
  family is every conjunction of at most two conditions, less empty
  intersections and size/event failures, and selection ranges over that
  family. A finite value truncates the pool under a *preview* ordering,
  so the subgroup that is optimal under the actual criterion can be
  discarded before it is ever evaluated – most acutely for the size foci
  `"maxSG"` / `"minSG"`, where small high-effect subgroups sort last,
  and for the band foci `"effMaxSG"` / `"effMinSG"`. When a finite value
  does truncate the pool a warning naming both counts is emitted.

  Set a finite value only to bound runtime, knowing it can change which
  subgroup is selected. The default was `200` in earlier versions.

- vi.grf.min:

  Numeric or `NULL`. GRF variable-importance screening of the candidate
  cut columns. `NULL` (the default since 0.3.0) skips the screening
  entirely: no causal forest is fitted and the cut columns are searched
  in the order
  [`get_FSdata()`](https://larry-leon.github.io/forestsearch/reference/get_FSdata.md)
  builds them. A numeric value fits a causal forest
  ([`grf::causal_forest()`](https://rdrr.io/pkg/grf/man/causal_forest.html)
  for GLM outcomes,
  [`grf::causal_survival_forest()`](https://rdrr.io/pkg/grf/man/causal_survival_forest.html)
  for survival), **orders** the cut columns by variable importance, and
  retains those with `vi_ratio > vi.grf.min`, where `vi_ratio` is each
  column's importance divided by the maximum. At values `<= 0` nothing
  is filtered out –
  [`grf::variable_importance()`](https://rdrr.io/pkg/grf/man/variable_importance.html)
  is non-negative and the block is guarded by `vi_max > 0` – so the
  effect is ordering only. That ordering sets the candidate matrix's
  column order; it does **not** reach candidate selection, because the
  candidate table is re-sorted by `(-HR, K)` before the consistency
  loop. Measured effect on `n_candidates_evaluated` and
  `early_stop_candidate`: zero; the selected subgroup's subject
  membership is unchanged and only the clause order of `sg.harm` may
  differ. Versions before 0.3.0 defaulted to `-0.2`; pass
  `vi.grf.min = -0.2` explicitly to reproduce their output.
  `max_n_confounders` is applied **only when this screening runs** – see
  that argument.

- use_twostage:

  Logical. Use two-stage sequential consistency algorithm for improved
  performance. Default FALSE for backward compatibility. When TRUE,
  `fs.splits` becomes the maximum number of splits, and early stopping
  is enabled. See Details. Overridden to `FALSE` under
  `sg_focus = "maxeff"` (no two-stage screening).

- twostage_args:

  List. Parameters for two-stage algorithm (only used when
  `use_twostage = TRUE`):

  n.splits.screen

  :   Integer. Splits for Stage 1 screening. Default 30.

  screen.threshold

  :   Numeric. Consistency threshold for Stage 1. Default is
      automatically calculated to provide ~2.5 SE margin.

  batch.size

  :   Integer. Splits per batch in Stage 2. Default 20.

  conf.level

  :   Numeric. Confidence level for early stopping. Default 0.95.

  min.valid.screen

  :   Integer. Minimum valid Stage 1 splits. Default 10.

- consistency_method:

  Character. `"resample"` (default) uses the multiplier
  (influence-function / `dfbeta`) approximation
  ([`consistency_resample`](https://larry-leon.github.io/forestsearch/reference/consistency_resample.md)),
  returning the rate from a single subgroup fit; `"split"` is retained
  as the fallback and runs the literal repeated 50/50 split-and-refit
  calculation, which costs `fs.splits * 2` refits per candidate
  evaluation. Despite the name, `"resample"` performs no numerical
  resampling at evaluation time – the rate is a closed-form expression.
  The name reflects the method's derivation: the repeated-split
  consistency rate admits an equivalent multiplier-resampling (MR)
  representation, and taking that representation to its analytic limit
  yields the closed form implemented here. The `"resample"` path applies
  to the Cox (survival) outcome and to GLM outcomes whose effect is a
  single model coefficient (OR, RR, RD, MD, IRR); configurations it
  cannot represent that way (IRD, propensity-adjusted effects, or a
  non-convergent fit on a given candidate) fall back to `"split"`
  automatically. When `use_twostage = TRUE`, `"resample"` bypasses the
  two-stage split screening entirely (the rate is computed directly).

- outcome_type:

  Character. One of `"survival"` (default), `"binary"`, `"continuous"`,
  or `"count"`.

- effect_measure:

  Character or `NULL`. Effect measure for GLM outcomes. For binary:
  `"OR"` (default), `"RD"`, `"RR"`, `"IRR"`, `"IRD"`. For continuous:
  `"MD"` (default). For count: `"IRR"` (default), `"IRD"`. `NULL`
  resolves to the per-outcome default shown here; the resolved measure
  then fixes the scale on which `effect.threshold` and
  `consistency.threshold` are compared – log for the ratio measures
  (`"OR"`, `"RR"`, `"IRR"`), identity for `"RD"`, `"IRD"` and `"MD"` –
  so an unset binary call screens on log(1.25) rather than 0.05. See the
  *Threshold vocabulary and resolved defaults* section. Ignored for
  `outcome_type = "survival"`, which always uses the hazard ratio.

- offset.name:

  Character or `NULL`. Name of the follow-up time column for rate-based
  measures (IRR, IRD).

- adverse_outcome:

  Logical or `NULL`. If `TRUE`, higher outcome values indicate harm
  (default for binary and count). If `FALSE`, higher values indicate
  benefit (default for continuous).

- overdispersion:

  Character. Overdispersion correction for `outcome_type = "count"`. One
  of `"none"` (default, standard Poisson), `"quasi"` (quasi-Poisson), or
  `"negbin"` (negative binomial via
  [`MASS::glm.nb`](https://rdrr.io/pkg/MASS/man/glm.nb.html)). Ignored
  for other outcome types.

- grf_count_transform:

  Character. Transformation applied to the count outcome before passing
  to
  [`grf::causal_forest()`](https://rdrr.io/pkg/grf/man/causal_forest.html)
  for factor screening. `"log"` (default) applies \\\log(Y + 0.5)\\;
  `"identity"` passes raw counts. Ignored for non-count outcomes.

- tune_grf:

  Logical. If `TRUE`, enables cross-validated hyperparameter tuning for
  GRF's
  [`causal_forest()`](https://rdrr.io/pkg/grf/man/causal_forest.html)
  via `tune.parameters = "all"`. If `FALSE` (default), all GRF
  parameters remain at their defaults, preserving existing behavior.
  Applies to both internal GRF pre-screening and standalone GRF analysis
  when passed through
  [`run_simulation_analysis()`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md).

- adjust_covariates:

  Character vector or `NULL`. Additional terms used to adjust the Cox
  model that scores survival subgroups during consistency evaluation.
  Terms are appended verbatim to the model right-hand side after the
  treatment term, so survival modelling helpers may be used directly:
  `adjust_covariates = "strata(x1)"` fits a model stratified by `x1`,
  whereas `adjust_covariates = "x1"` includes `x1` as a linear
  covariate. The referenced raw columns must be present in
  `df.analysis`; they are carried into the internal scoring frame
  automatically. Applies to the survival (Cox) consistency path and to
  GLM outcomes (the effect-estimator closure adds the terms to the
  outcome model, e.g. logistic for binary/OR). The candidate-search
  ranking is adjusted on both paths. Mutually exclusive with
  `ps_method != "none"`. Default `NULL` (treatment-only, unadjusted –
  the previous behaviour).

- ps_method:

  Character or `NULL`. Propensity score estimation method: `"grf"`,
  `"lasso"`, `"logistic"`, or `"none"`. Default: `"none"` for RCT,
  `"grf"` for observational.

  **GLM outcome types only** (`"binary"`, `"continuous"`, `"count"`).
  Propensity adjustment is **not implemented for
  `outcome_type = "survival"`**: no Cox estimation path reads the
  resulting weights, so any value other than `"none"` is an error there
  rather than a silent no-op. For confounding control on a survival
  outcome use `adjust_covariates`. Note the default resolves to `"grf"`
  when `is.RCT = FALSE`, so an observational survival analysis must pass
  `ps_method = "none"` explicitly.

- ps_adjust_method:

  Character. PS adjustment method: `"none"` (default), `"iptw"`
  (stabilized inverse probability of treatment weights), or `"dr_gcomp"`
  (IPS covariate).

  **GLM outcome types only**, and consumed only when `ps_method` is not
  `"none"` – it selects how the estimator closure is rebuilt, and that
  rebuild is GLM-gated. It has no effect on a survival outcome, which is
  refused at `ps_method`.

- ps_hat:

  Numeric vector or `NULL`. User-supplied propensity scores. Must have
  length equal to `nrow(df.analysis)`, and is applied **positionally**
  to `df.analysis`.

  Not carried into bootstrap replicates or cross-validation folds: a
  resample has the same row count but different subjects, so a supplied
  score would be attached to the wrong people. Both re-estimate their
  own score per replicate/fold under `ps_method`.

- mr_inference:

  Logical. Whether to follow the search with multiplier resampling (MR):
  a refit-free de-biased estimate of the selected subgroup's treatment
  effect, an interval for it, and a harm-confirmation flag. **Default
  `FALSE`**, so by default this function performs *identification only*
  – no de-biasing, no harm confirmation, and `mr_inference` /
  `mr_harm_confirmed` come back `NULL` / `NA`. Setting it `TRUE` adds a
  step that runs *after* the subgroup has been chosen and cannot change
  that choice. Supported on the GLM resample-consistency path
  (`consistency_method = "resample"` with a GLM effect) and on the
  survival/Cox path (`outcome_type = "survival"`), where it runs under
  the default split-consistency search. MR approximates the full
  bootstrap
  ([`forestsearch_bootstrap_dofuture`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md))
  to leading order and does not replace it. See the vocabulary section.

- mr_inference_args:

  List of optional MR controls; used only when `mr_inference = TRUE`.

  `confirm_rule`

  :   Which harm-confirmation rule to apply to the de-biased estimate:
      `"point"` (default) compares the de-biased point estimate to
      `t_confirm`; `"ci"` compares the one-sided 95% selection-adjusted
      lower bound, and is therefore the stricter requirement.

  `t_confirm`

  :   The threshold `confirm_rule` compares against, on the **effect**
      scale (not the working log scale). `NULL` (default) resolves to
      the null effect: `1` for ratio measures, `0` for differences. Set
      near the null rather than at the screening threshold, because the
      correction over-shrinks true effects.

  `reselection`

  :   Which selection rule to replay on each perturbed draw when forming
      the bias term. Defaults to the rule `sg_focus` implies, so it
      normally needs no setting.

  `draws`

  :   Number of multiplier draws. Default `2000`.

  `multiplier`

  :   Multiplier distribution, default `"poisson"` (centred Poisson(1)
      weights, the multiplier-bootstrap analogue of
      nonparametric-bootstrap multiplicities); also `"gaussian"`,
      `"rademacher"`.

  `ci_method`

  :   `"field"` (default) additionally runs the field-calibrated block,
      whose one-sided bounds are the primary reported products; the
      de-biased element, its infinitesimal-jackknife SE and its IJ
      intervals are returned in every call regardless, so nothing is
      removed – only a per-fit Monte Carlo cost is added. `"ij"`
      restores the previous default: the de-biased interval from the
      infinitesimal-jackknife variance, the leading-order analogue of
      the full-bootstrap interval, with the field block omitted.
      `"wald"` uses the subgroup robust SE and also omits the field
      block. The naive interval always uses the robust SE.

  `include_complement`

  :   Also de-bias the complement subgroup. Passed as `TRUE` by every
      internal caller.

  `seed`

  :   Seed for the multiplier draws; defaults to `seedit`.

  `return_reselection`

  :   Logical, default `FALSE`. When `TRUE` the gate's return gains the
      `reselection` element (per-draw winners and the re-selection
      frequencies `p_hat`); add-only pass-through, the default
      reproduces prior output exactly.

  `field_M_cap`

  :   Optional override of the uniform (kappa) sweep's mass-carrying
      cap; `NULL` (default) keeps
      [`fs_mr_field_uniform()`](https://larry-leon.github.io/forestsearch/reference/fs_mr_field_uniform.md)'s
      default. Consulted only when `field_uniform = TRUE`.

  `field_complement`

  :   Logical, default `FALSE`. Under `ci_method = "field"` with
      `include_complement = TRUE`, also evaluate the field block for the
      complement subgroup, attached as `field$complement` with the
      one-sided 95\\ bound as its primary product (the benefit-claim
      orientation; see
      [`fs_mr_inference`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)).
      Add-only pass-through: the default reproduces prior output
      exactly.

  `ij_residual`

  :   `"two_term"` (default), `"winner"` or `"winner_floor"`: which
      infinitesimal-jackknife residual populates the reported de-biased
      SE and interval (both blocks); all three SEs are returned side by
      side regardless (see
      [`fs_mr_inference`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)).
      Add-only pass-through: the default reproduces prior output
      exactly.

  `keep_declaration_field`, `keep_field_matrix`

  :   Logical, both default `FALSE`. Retain the standardized
      perturbation field that
      [`fs_declaration_calibration`](https://larry-leon.github.io/forestsearch/reference/fs_declaration_calibration.md)
      reads, as `mr_inference$declaration_field`; `keep_field_matrix`
      additionally keeps the full draws-by-candidates matrix (needed for
      `family = "reduced"`). Add-only pass-through: unset, nothing is
      computed or stored and the default reproduces prior output
      exactly.

  `declaration_c0`

  :   `NULL` (default) or a numeric vector of pre-specified protected
      null levels `c0` on the natural scale of `hr.consistency`
      (`c0 <= c2`); consulted only with `keep_declaration_field = TRUE`.
      The capture then also stores the shifted family maxima `Mstar_c0`
      that `fs_declaration_calibration(fit, c0 = )` reads (see
      [`fs_mr_inference`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)).
      Add-only pass-through: unset, nothing is computed or stored.

## Value

A list of class "forestsearch" containing:

- grp.consistency:

  Consistency evaluation results including:

  - out_sg: Selected subgroup based on sg_focus

  - sg_focus: Focus criterion used

  - df_flag: Treatment recommendations

  - sg.harm: Subgroup definition labels (character vector of cut names)
    – same as the top-level `sg.harm`

  - sg.harm.id: per-subject 0/1 membership indicator (*not* a character
    vector of cut expressions; see the *Field naming collision* section
    below)

  - algorithm: "twostage" or "fixed"

  - n_candidates_evaluated: Number evaluated

  - n_passed: Number passing threshold

- find.grps:

  Subgroup search results

- confounders.candidate:

  Candidate confounders considered

- confounders.evaluated:

  Confounders after variable selection

- df.est:

  Analysis data with treatment recommendations

- df.predict:

  Prediction data with recommendations (if provided)

- df.test:

  Test data with recommendations (if provided)

- minutes_all:

  Total computation time

- grf_res:

  GRF results object

- sg_focus:

  Subgroup focus criterion used

- sg.harm:

  Selected subgroup definition – character vector of factor-level cut
  names. **Recommended accessor** for displaying the identified subgroup
  (e.g., `paste(obj$sg.harm, collapse = " & ")`). `NULL` if no subgroup
  was identified.

- grf_cuts:

  GRF cut points used

- prop_maxk:

  Proportion of max combinations searched

- max_sg_est:

  Maximum subgroup HR estimate

- grf_plot:

  GRF plot object (if plot.grf = TRUE)

- args_call_all:

  All arguments for reproducibility

- threshold_config:

  The resolved threshold configuration: `outcome_type`,
  `effect_measure`, `screening`, `consistency`, `screening_natural`,
  `consistency_natural`, `scale`, `pconsistency` and two description
  strings. **Read `$scale` before reading `$screening`.** The phrase
  "the screening threshold" names two different numbers depending on
  which object is inspected:

  - On the **survival** path `$screening` is `log(hr.threshold)` – the
    scale the admission set and multiplier resampling work on – while
    the candidate search itself compares the fitted hazard ratio against
    `hr.threshold` on the **natural** scale. The two are the same
    threshold on two scales, and `$screening_natural` holds the search's
    value.

  - On the **GLM** paths `$screening` is the value the search compares
    against directly: the link scale for ratio measures (`log` of the
    supplied ratio) and the identity scale for `RD`, `IRD` and `MD`.

  `$consistency` follows the same convention, and `$pconsistency` is the
  rate `pconsistency.threshold`, which has no scale.

- family_status:

  Character scalar recording whether the candidate family the identifier
  ranked over is fixed – the condition multiplier resampling requires
  (see `mr_inference`). `"no-front-end"`:
  `subgroup_method = "consistency"` with `use_lasso`, `use_grf` and
  `use_dina` all `FALSE`. Named for what it checks: no fitted model
  shapes the family on the observed data. This is weaker than the
  manuscript's Section 2.1 fixed family, which additionally requires
  resample- invariant cut locations – quantile-derived cuts are not.
  `"conditional-removable"`: the same method with at least one front end
  on – data-dependent, but fixed by turning those off.
  `"conditional-inherent"`: `subgroup_method` of `"dina"` or `"grf"`,
  where the family is generated by fitting a model to the same data and
  cannot be made fixed. Descriptive only; it never raises a condition.

- mr_inference:

  Multiplier-resampling (MR) result, or `NULL` when MR did not run. See
  [`fs_mr_inference`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)
  for fields (`naive`, `debiased`, `selection_bias`, `settings`,
  `harm_flag`, `timing_seconds`). MR approximates the full bootstrap
  (FB) to leading order; its interval is the infinitesimal-jackknife
  analogue of the FB interval under the default `ci_method = "field"`
  and under `"ij"`, and the subgroup robust-SE interval under `"wald"`.
  Under the default the return also carries the `field` element.

- mr_harm_confirmed:

  Logical, three-valued. `TRUE` when MR ran and the de-biased estimate
  still indicates harm under `confirm_rule`; `FALSE` when MR ran and it
  does not; `NA` when MR did **not** run or could not be computed. `NA`
  is therefore *absence of evidence, not evidence against harm* – test
  it with [`is.na()`](https://rdrr.io/r/base/NA.html) before
  [`isTRUE()`](https://rdrr.io/r/base/Logic.html), which collapses `NA`
  to `FALSE` and so reports "harm not confirmed" for an analysis that
  never ran. See the vocabulary section.

## Details

**Algorithm Overview:**

1.  **Variable Selection**: GRF identifies variables with treatment
    effect heterogeneity; LASSO selects most predictive

2.  **Subgroup Discovery**: Exhaustive search over factor combinations
    up to `maxk`

3.  **Consistency Validation**: Split-sample validation ensures
    reproducibility

4.  **Selection**: Choose subgroup based on `sg_focus` criterion

**Two-Stage Consistency Algorithm:** When `use_twostage = TRUE`, the
consistency evaluation uses an optimized algorithm that can provide
3-10x speedup:

- **Stage 1**: Quick screening with `n.splits.screen` splits eliminates
  clearly non-viable candidates

- **Stage 2**: Sequential batched evaluation with early stopping for
  candidates passing Stage 1

The two-stage algorithm is recommended for:

- Exploratory analyses with many candidate subgroups

- Large `fs.splits` values (\>200)

- Iterative model development

For final regulatory submissions, `use_twostage = FALSE` may be
preferred for exact reproducibility.

## Threshold vocabulary and resolved defaults

Three arguments carry thresholds, and two of them have the word
"consistency" in the name. They are different quantities:

- `c1` – `effect.threshold` (legacy `hr.threshold`):

  The **screening** threshold on a candidate subgroup's own effect,
  estimated on the full sample. A candidate must clear it to enter the
  consistency stage.

- `c2` – `consistency.threshold` (legacy `hr.consistency`):

  An **effect** threshold, applied to the subgroup's effect estimated
  within each split half. **Not** a proportion, and **not** the
  consistency rate.

- `p*` – `pconsistency.threshold`:

  A **proportion** in `[0, 1]`: the fraction of splits that must clear
  `c2`. **Not** an effect threshold.

In one line: `c2` sets the bar a split half must clear, and `p*` counts
how often that bar is met.

**Resolved defaults, per estimand.** Reading `hr.consistency`'s default
as "1.0" is misleading on the GLM paths: when neither threshold is
supplied, the value a run actually uses is resolved from the estimand.
The resolved pairs are:

|  |  |  |
|----|----|----|
| **resolved estimand** | **c1** | **c2** |
| `HR` (survival) | 1.25, compared by the search on the natural HR scale | 1.0 |
| `OR`, `RR`, `IRR` | log(1.25) | log(1.0), i.e. 0 |
| `RD` | 0.05 | 0.0 |
| `IRD` | 0.01 | 0.0 |
| `MD` | 0.0 | 0.0 |

Supplying `effect.threshold` or `consistency.threshold` (or the legacy
names) suppresses the corresponding remap and the supplied value is used
on the scale shown above.

**The pair rule, on the consistency path with a ratio estimand.** For
`subgroup_method = "consistency"` with survival (HR) or binary
`effect_measure = "OR"` – the two estimands whose `c1` and `c2` sit on
one comparable ratio scale – the two thresholds are treated as a pair:

- `c2 > c1` is an **error**. It is degenerate: the consistency-stage
  entry condition re-screens the admitted family on `c2`, so candidates
  the screen let through are discarded at the stage boundary, and when
  none survives the consistency stage does not run at all and the fit
  reads as a finding of no subgroup. `c2 = c1` is allowed.

- `c1` supplied **without** `c2` derives `c2 = 0.80 * c1` on the
  **ratio** scale (on the log scale an additive shift of `log(0.80)`,
  never `0.80 * log(c1)`): 1.25 -\> 1.00, 1.00 -\> 0.80, 0.90 -\> 0.72.
  Announced with a message, and carried into every bootstrap replicate
  and CV fold. The default pair (1.25, 1.0) is the rule's fixed point,
  so a call that sets neither threshold is unaffected. An explicitly
  supplied `c2`, in either spelling, is never overridden.

- Supplying both spellings of one threshold at disagreeing values is an
  error naming both.

Neither rule applies under `subgroup_method = "dina"` or `"grf"` (which
return before the consistency stage and never consult `c2`), nor to the
identity-scale estimands `RD`, `IRD`, `MD`, nor to `IRR`.

`pconsistency.threshold` is **not** remapped. It is a rate, and it is
identical for every `outcome_type` and every estimand.

## Field naming collision with GRF results

The top-level `sg.harm` on this object and the `grp.consistency$sg.harm`
nested field both hold a character vector of cut expressions. The nested
`grp.consistency$sg.harm.id` field holds a *per-subject 0/1 membership
indicator*, not a cut vector – unlike `sg.harm.id` on GRF result
objects, which *does* hold cut expressions. See
[`subgroup.consistency`](https://larry-leon.github.io/forestsearch/reference/subgroup.consistency.md)
and
[`grf.subg.harm.glm`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.glm.md)
for the full discussion of this naming collision. Code that must handle
both object types should prefer the top-level `sg.harm` accessor on
forestsearch results.

## Vocabulary FB MR and harm confirmation

Four different things in this package used to share the word "gate".
They are distinct, they happen at different stages, and only the last
one is a pass/fail decision:

- Admissibility criteria:

  Which candidate subgroups are even eligible to be considered:
  `hr.threshold`, `pconsistency.threshold`, `n.min`, `d0.min`, `d1.min`,
  `maxk`. Applied during the search.

- Selection rule:

  Which admissible candidate is chosen as the estimated subgroup:
  `sg_focus`. Applied during the search.

- Post-selection inference:

  Correcting the chosen subgroup's treatment effect for the optimism of
  having chosen it, and putting an interval around the corrected value.
  This is FB and MR below. Applied *after* the search.

- Harm confirmation:

  The one pass/fail: does the corrected estimate still indicate harm?
  This is `confirm_rule` against `t_confirm`, reported as
  `mr_harm_confirmed`.

**Full bootstrap (FB)** —
[`forestsearch_bootstrap_dofuture`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md).
Resamples *subjects*, and re-runs the entire subgroup search inside
every replicate, so the replicate may select a different subgroup than
the original analysis did. That is what lets it correct the
selection-induced bias in the reported treatment effect, and supply an
infinitesimal-jackknife interval. It is the reference method, and it is
expensive: cost scales with `nb_boots` times the cost of one search.

**Multiplier resampling (MR)** —
[`fs_mr_inference`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md),
turned on by `mr_inference = TRUE`. Resamples nothing; it perturbs the
per-subject influence contributions (`dfbeta`) of a *single* set of fits
with random mean-zero multipliers, and re-applies the selection rule on
each perturbed draw. It corrects the same selection bias FB corrects and
produces the same kind of interval, to leading order, at a small
fraction of the cost — no refits. MR is **post-selection inference on a
completed analysis**: it runs after the subgroup has been chosen and
**cannot change which subgroup is identified**. Running with
`mr_inference = TRUE` and `FALSE` gives identical `sg.harm`, `df.est`
and `max_sg_est`. MR approximates FB; it does not replace it.

**Harm confirmation.** Once MR has produced a de-biased estimate,
`mr_inference_args$confirm_rule` decides whether that estimate still
indicates harm:

- `"point"` (default):

  Compares the de-biased *point estimate* to `t_confirm`.

- `"ci"`:

  Compares the one-sided 95% selection-adjusted *lower bound* to
  `t_confirm`. Strictly the stronger requirement, so anything it
  confirms `"point"` also confirms.

**`t_confirm`** is the threshold those rules compare against, on the
**effect** scale (HR, OR, RR, IRR; or RD, MD for differences) — not the
working log scale used internally. Left `NULL` it resolves to the null
effect: `1` for ratio measures, `0` for differences. It is deliberately
set near the null rather than at `hr.threshold`, because the de-biasing
correction over-shrinks genuine effects, so re-applying the screening
threshold to a corrected estimate would discard true findings.

**`mr_harm_confirmed`** carries the answer, and is three-valued: `TRUE`
= MR ran and harm is confirmed; `FALSE` = MR ran and harm is not
confirmed; `NA` = MR did not run, or could not be computed. **`NA` is
not evidence against harm.** It arises when `mr_inference = FALSE` (the
default), when no subgroup was identified, or when the MR step was
skipped or failed for the outcome/consistency combination in use.
Because `isTRUE(NA)` is `FALSE`, code that reads this field with
[`isTRUE()`](https://rdrr.io/r/base/Logic.html) alone silently reports
"harm not confirmed" for an analysis that was never performed; branch on
[`is.na()`](https://rdrr.io/r/base/NA.html) first.

## References

- FDA Guidance for Industry: Enrichment Strategies for Clinical Trials

- Athey & Imbens (2016). Recursive partitioning for heterogeneous causal
  effects. PNAS.

- Wager & Athey (2018). Estimation and inference of heterogeneous
  treatment effects using random forests. JASA.

## See also

[`subgroup.consistency`](https://larry-leon.github.io/forestsearch/reference/subgroup.consistency.md)
for consistency evaluation details
[`forestsearch_bootstrap_dofuture`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)
for the full bootstrap (FB)
[`fs_mr_inference`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)
for multiplier resampling (MR), the `mr_inference = TRUE` step
[`mr_estimates_table`](https://larry-leon.github.io/forestsearch/reference/mr_estimates_table.md)
to render the MR estimates
[`fs_fdr_report`](https://larry-leon.github.io/forestsearch/reference/fs_fdr_report.md)
to sweep harm-confirmation thresholds (`c_confirm`) and measure what
confirmation buys
[`forestsearch_Kfold`](https://larry-leon.github.io/forestsearch/reference/forestsearch_Kfold.md)
for cross-validation

Package website: <https://larry-leon.github.io/forestsearch/>

Source code: <https://github.com/larry-leon/forestsearch>

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: Standard analysis (backward compatible)
result <- forestsearch(
  df.analysis = trial_data,
  sg_focus = "hr",
  hr.threshold = 1.25,
  pconsistency.threshold = 0.90,
  fs.splits = 400,
  details = TRUE
)

# Example 2: Fast exploratory analysis with two-stage
result_fast <- forestsearch(
  df.analysis = trial_data,
  sg_focus = "maxSG",
  hr.threshold = 1.15,
  pconsistency.threshold = 0.85,
  fs.splits = 500,
  use_twostage = TRUE,
  details = TRUE
)

# Example 3: Two-stage with custom parameters
result_custom <- forestsearch(
  df.analysis = trial_data,
  sg_focus = "hr",
  hr.threshold = 1.3,
  pconsistency.threshold = 0.95,
  fs.splits = 600,
  use_twostage = TRUE,
  twostage_args = list(
    n.splits.screen = 50,
    batch.size = 25,
    conf.level = 0.99
  ),
  parallel_args = list(plan = "multisession", workers = 4),
  details = TRUE
)

# Example 4: Binary outcome with OR thresholds
result_binary <- forestsearch(
  df.analysis = trial_data,
  outcome_type = "binary",
  effect_measure = "OR",
  effect.threshold = 1.5,
  consistency.threshold = 1.3,
  pconsistency.threshold = 0.90
)

# Example 5: Binary outcome with RD thresholds
result_rd <- forestsearch(
  df.analysis = trial_data,
  outcome_type = "binary",
  effect_measure = "RD",
  effect.threshold = 0.07,       # 7 pct-point harm
  consistency.threshold = 0.03   # 3 pct-point per split
)
} # }
```
