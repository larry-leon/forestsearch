# forestsearch 0.5.0

forestsearch 0.5.0 removes `create_gbsg_dgm()` and `simulate_from_gbsg_dgm()`
(use `generate_aft_dgm_flex()` with `simulate_from_dgm()`; the superseded
`setup_gbsg_dgm()` rebuilds the legacy mechanism), needs survival (>= 3.5-8) and
adds observed-trial and analytic tools for MRCT regional consistency. Results
can change where fixes apply: GLM bootstrap, table and plot effects (on the fit's
`adverse_outcome`), consistency-path bias correction, logistic-tier RD SEs, the
RD, RR and IRD scales and the detection, PPV, NPV and specificity of
`run_simulation_analysis()`, `test_hte_from_forestsearch()`, and some GRF and DINA selections.

### Breaking changes

* `forestsearch_bootstrap_dofuture()` on GLM outcomes now builds the estimator
  behind its bias correction with the parent fit's resolved `adverse_outcome`
  and `ps_adjust_method` (it used `adverse_outcome = TRUE` and no propensity
  weighting). For a continuous fit with `adverse_outcome` unset (it resolves to
  `FALSE`), or any fit with `adverse_outcome = FALSE`, `H_estimates`,
  `Hc_estimates`, `SG_CIs`, the `*_biasadj_*` columns and the effect columns of
  `FSsg_tab` were on the raw outcome scale, the opposite sign to the fit; they
  are now on the fit's scale (an exact sign change for MD, RD and log OR; for
  RR the old values were the ratio of 1 - Y risks). With propensity adjustment
  the bootstrap estimates are now weighted like the fit. Fits with
  `adverse_outcome = TRUE` and no propensity adjustment, which include binary
  and count fits at defaults, are unchanged. `forestsearch()` with propensity
  adjustment on a GLM outcome also keeps `adverse_outcome` when it rebuilds its
  estimator (it used `TRUE`). Scripts that negated the bootstrap output by
  `adverse_outcome` should drop that step. (Code audit A-009.)

* `forestsearch()` and the other exported functions that can reach the fixed
  seed in `subgroup.consistency()` (`set.seed(8316951)`, which `forestsearch()`
  does not expose) no longer leave the caller's random-number state changed.
  These are `subgroup.consistency()` itself and the functions that call
  `forestsearch()`: `forestsearch_bootstrap_dofuture()`, `bootstrap_results()`,
  `forestsearch_Kfold()`, `forestsearch_tenfold()`, `run_simulation_analysis()`,
  `mrct_region_sims()`, `mrct_case_study()`, `mrct_trial_analysis()`,
  `compare_selection_rules()`, `fpr_calibration()`, `fpr_calibration_transfer()`,
  `run_null_calibration()` and `fs_fdr_report()`. The outermost such call saves
  `.Random.seed` on entry and restores it on exit (absent again if it was
  absent), whatever seeds are set inside it. Their own results are unchanged:
  every draw inside a call comes from the same stream as before. A later call
  that draws random numbers without its own seed, straight after
  `forestsearch()` (for example `sample()`, or `simulate_from_dgm()` with
  `seed = NULL`), used to continue the stream `forestsearch()` had seeded and
  now depends on the caller's seed, so set one. `forestsearch_bootstrap_dofuture()`
  seeds itself (`seed = 8316951L`) and is not affected. Functions that only set a
  seed passed to them (for example `simulate_from_dgm(seed = )`) still set it.

* `fpr_calibration()` with neither `outcome_type_approx` nor
  `fs_params$outcome_type` given now approximates `P1` for `"survival"`,
  `forestsearch()`'s default (it used `"binary"`), so the printed `P1`,
  `d_eff` and implied `L_eff` change; the simulated false-positive rate does
  not. An outcome type other than `"survival"`, `"binary"`, `"continuous"` or
  `"count"` now stops with an error naming it. (Code audit A-107.)

* `create_gbsg_dgm()` and `simulate_from_gbsg_dgm()` are removed, with the
  internal `.simulate_from_gbsg_dgm_()` and `get_dgm_with_output()`. Both
  exported functions had warned that they were deprecated since before 0.1.0,
  and nothing in the package, its tests, its vignettes and articles, or the
  quarto analyses called any of the four. Use `generate_aft_dgm_flex()` with
  `simulate_from_dgm()`; to reproduce the legacy GBSG mechanism exactly, use
  `setup_gbsg_dgm()`, which builds it with the same internal builder. No
  remaining function changes behaviour.

* `survival (>= 3.5-8)` is now required: the lean unadjusted search-stage
  Cox fit calls `survival::coxph.fit()` with the calling conventions
  verified on survival 3.5-8 and later.

### Bug fixes

* `mrct_region_sims()` now restores the caller's `progressr.handlers` option
  when it returns or stops. It set the progress handler with
  `progressr::handlers("progress")` and left that setting in the session. The
  progress display during the call and the returned values are unchanged.

* `forestsearch()` with a continuous outcome, DINA (`subgroup_method = "dina"`,
  or `use_dina = TRUE`), and neither `effect.threshold` nor `hr.threshold`
  supplied now uses the resolved screening threshold (0 for a mean difference)
  as the DINA harm floor, the value the returned object records and that
  bootstrap and cross-validation replays use. The fit used the raw survival
  default 1.25. The selected subgroup and the DINA screening cut can change for
  such analyses. A supplied threshold, and every non-gaussian DINA family, give
  the same floor as before. (Code audit A-001.)

* Survival GRF (`grf.subg.harm.survival()`, and `forestsearch()` with
  `use_grf = TRUE` or `subgroup_method = "grf"`) now codes the covariate matrix
  column by column as the GLM path does. When a confounder was a factor or
  character column, every numeric column was rounded to 7 significant digits
  and factor labels such as `"EU"`/`"US"` became `NA` (the GRF step then stopped,
  or its cuts were dropped with a warning). Numeric columns now keep their exact
  values, so the matrix GRF sees and the data a selected cut is evaluated on
  agree. Results with only numeric confounders are unchanged. (Code audit
  A-059.)

* `run_simulation_analysis()` and `fs_mr_oc_summary()`: a replicate whose
  identified subgroup is the whole sample has an empty estimated complement, and
  its NPV (0/0) is now 0, by the same convention that gives PPV = 0 for an empty
  identified subgroup. `run_simulation_analysis()` stored `NA` (dropped by later
  averages) in both the ForestSearch and GRF rows; `fs_mr_oc_summary()`'s
  unconditional NPV filled such a detected replicate with 1. Every other missing
  NPV keeps its previous treatment, and stored payloads are not regenerated.
  This affects simulation summaries in which a method can select everyone.

* `generate_gbsg_bootstrap_general()` now takes its data from `survival::gbsg`.
  `data("gbsg", package = "survival")` loads nothing under current survival, so
  every call stopped with an error. (Code audit A-004.)

* Functions that set a `future` plan for their own work now restore the
  caller's plan on exit, including after an error: `subgroup.search()` (and so
  `forestsearch()` and the functions that call it in the session),
  `forestsearch_bootstrap_dofuture()`, `mrct_region_sims()`,
  `fpr_calibration()` (with `n_workers > 1`) and `fs_fdr_report()`. They left
  the session on a sequential plan, replacing a plan the user had set. Results
  are unchanged. `setup_parallel_SGcons()` and `reset_workers()` still set the
  plan, which is their purpose. (Code audit INV-2, plan half.)

* `run_simulation_analysis()`: a replicate whose identified set is empty is now
  scored as a no-find on every path (`any.H` = 0; `ppv` = 0 under an
  alternative). The detection rate and the all-replicate PPV average can change
  where such replicates occur.

* `forestsearch_KfoldOut()`: when the folds identify no subgroup and the
  full-data analysis has one, the agreement metrics (`sens_H`, `ppv_H`,
  `sens_Hc`, `ppv_Hc`) are now 0, 0, 1 and the full-data complement share
  (the sensitivity for H was reported as 1). When the folds put everyone in
  the subgroup, `sens_Hc` and `ppv_Hc` are now 0 (they were NA).

* `fs_mr_oc_summary()` and `build_classification_table()`: a replicate with no
  identified subgroup now contributes NPV = |Hc|/N, the share of that
  replicate's subjects truly in the complement (it was 1). Under the null this
  is still 1. Specificity, sensitivity, PPV and detected-only summaries are
  unchanged. Stored payloads are not regenerated.
  - `fs_mr_oc_summary()` (unconditional convention) takes |Hc|/N as
    `1 - n_true / meta$n_sample`. Its unconditional NPV is `NA` when that is
    not available for every non-detected replicate (no numeric `n_true` or
    usable `meta$n_sample`, or `n_true` missing on a non-detected row).
  - `build_classification_table()` changes only rows stored with `npv = NA`
    and `any.H == 0` under an alternative (output of an older package version
    or another producer). It takes |Hc|/N from the row's `propHc_true`. A row
    without it (no truth column) is still filled with 1.

* `run_simulation_analysis()`: when no subgroup is identified under an
  alternative, `spec` is now 1 and `npv` is |Hc|/N (the two were stored the
  other way round). All-replicate averages of specificity and NPV from
  `build_classification_table()`, `format_oc_results()` and
  `summarize_single_analysis()` change accordingly. Replicates with a subgroup
  found, null scenarios, `sens` and `ppv` are unchanged.

* Bootstrap bias correction on the consistency path now evaluates the
  bootstrap-selected subgroup on the original data and on the bootstrap sample
  as Eq. 7 defines. Previously both terms were computed on the resample, through
  a merge on repeated ids. Bias-corrected H/Hc estimates from
  `forestsearch_bootstrap_dofuture()` after a consistency-path fit change.
  Naive estimates, selection, DINA/GRF and multiplier resampling are unchanged
  (A-002).

* `expand_covariates()` (and so `inject_mrct_structure()`) now warns when
  the latent correlation is not fitted by latentcor (package not installed,
  or its call failed) and the normal-score correlation is used instead. The
  computation is unchanged.

* `summarize_bootstrap_subgroups()` (and so `summarize_bootstrap_results()`)
  now counts a bootstrap consistency or subgroup size that lies exactly on a
  bin edge in the bin whose label starts there: a consistency of 0.95 is in
  ">=0.95" and a size of 50 in "50-99". Both were counted one bin lower.
  Values off an edge are binned as before. (Code audit A-147.)

* `forestsearch_KfoldOut()` no longer stops with "subscript out of bounds"
  when the held-out assignments are one level plus `NA` (`{1, NA}` or
  `{0, NA}`); the agreement values follow the existing single-level rule.
  `forestsearch_tenfold()`, which stopped too, now completes, and
  `forestsearch_Kfold()`, which caught the error and returned all-`NA`
  agreement summaries, now returns the values. Assignments `{0, 1}` and
  `{0, 1, NA}` give the same values as before.

* Missing values in a covariate that defines the subgroup now say so; no
  returned value changes. `run_simulation_analysis()` stops in the same
  cases as before, but with "N of M subjects have no subgroup assignment: a
  covariate that defines the identified subgroup ... has missing values"
  instead of "missing value where TRUE/FALSE needed". The GRF route for
  binary, continuous and count outcomes warns once with the count and the
  defining covariates when the selected cut leaves subjects unassigned, and
  `forestsearch_Kfold()` and `forestsearch_tenfold()` warn once per call with
  the count of held-out assignments that are `NA`.

* Cross-validation "Exact" (`find_metrics["Exact"]`, and `exact_match` of
  `CV_sgs()`) for a one-factor full-data subgroup now counts a fold only when
  its subgroup is that factor alone. A fold that adds a second factor still
  counts for "At least 1" and "Cov1 exact", which keep their values. So the
  Exact rate can fall for such analyses; two-factor subgroups are unchanged.
  (Code audit A-103.)

* `test_hte_from_forestsearch()` now uses the continuous covariates of the
  fit: an evaluated cut label such as `"age <= 55"` is mapped to its column
  `age`. Only labels that were exact column names were used before, so
  continuous covariates were left out, and a fit with only continuous
  confounders returned `NULL`. A label that cannot be resolved is named in a
  warning, and when no covariate remains the function warns before returning
  `NULL`. (Code audit A-151.)

* `lasso_selection()` for binary, continuous and count outcomes no longer
  lists `"(Intercept)"` in `selected` or `omitted` (it stays in
  `coefficients`). When LASSO keeps no covariate, `forestsearch()` therefore
  keeps the full candidate pool, as it already did for survival, instead of
  failing at `get_FSdata` with "no candidate cut expressions could be
  constructed". With `details = TRUE` the printed line names the model by
  outcome type ("Logistic-LASSO", "Linear-LASSO", "Poisson-LASSO") and its
  counts no longer include the intercept. (Code audit A-042.)

* `get_FSdata()`, and so `forestsearch()` with GRF or DINA cuts, drops a cut
  as a duplicate only when its variable is exactly a categorical confounder or
  the whole cut is already in the pool. A cut on `bmi` is no longer dropped
  because a categorical confounder `bm` exists, and `x <= 25.3` is kept when
  the pool holds `x <= 2`. A median cut is removed only for a variable GRF cut
  exactly. (Code audit A-047.)

* `mrct_region_sims()` and `mrct_case_study()` without a `confounders.name`
  and with a mechanism built with a spline no longer include the three internal
  spline columns in the default confounders (`z_bm_treat`, `z_bm_k`,
  `z_bm_k_treat` for `create_dgm_for_mrct()`). Covariates such as
  `z_prior_treat`, and the region indicator, stay. (Code audit A-125.)

* `grf.subg.harm.glm()` with a count outcome and `adverse_outcome = TRUE` (the
  value `forestsearch()` and `run_simulation_analysis()` resolve for counts),
  and so `forestsearch()` with `use_grf = TRUE` or
  `subgroup_method = "grf"` and the GRF arm of `run_simulation_analysis()`,
  now looks for the region where treatment raises the count. It looked for
  the region where treatment lowers it. (Code audit A-064.)

* `forestsearch()` with `use_dina = TRUE` now orients the DINA screening cut
  as `subgroup_method = "dina"` does. For a continuous outcome with the
  default orientation, or a binary outcome with `adverse_outcome = FALSE`, the
  proposed cut was on the benefit side; it now points at the harm region, so
  the candidate pool and the selected subgroup can change. Survival, count and
  `adverse_outcome = TRUE` fits are unchanged.

* `plot_sg_glm_outcomes()`, `plot_subgroup_results_forestplot()` and
  `mr_estimates_table()` without `boot_results` now show a GLM fit's effect
  on the fit's own orientation (its `adverse_outcome`). For a continuous fit
  with the default orientation the displayed effect had the opposite sign of
  the fit's estimate. The raw arm means and their difference in the summary
  table are unchanged.

* Change of method: the risk-difference standard error of
  `make_effect_estimator(outcome_type = "binary", effect_measure = "RD")`
  without adjustment covariates is now the two-sample binomial form
  `sqrt(p1 (1 - p1) / n1 + p0 (1 - p0) / n0)` when the logistic tier is used.
  With 40 patients per arm, control 30% and treated 100%, it is 0.0725 (it was
  about 5e-06, giving a zero-width interval); with control 0% and treated
  80%, 0.0632 (it was 448.5). The estimate is unchanged, and so is the
  covariate-adjusted case, except a direct call with `adjust_covariates` and
  `ps_adjust_method = "dr_gcomp"` (or `"iptw"` without an `sw` column), whose
  covariates are dropped before the fit and which now gets the new form.
  (Code audit A-055.)

* `forestsearch()` with a binary risk difference or a count rate difference,
  `dina_args = list(family = "gaussian")` and a supplied ratio-scale
  threshold that is remapped with a warning (for example 1.5 to 0.05) now
  uses the remapped value as the DINA floor in the fit, the value the object
  records and the bootstrap and cross-validation replay. The fit used the
  supplied 1.5 and so usually found no subgroup. This narrows the A-001 entry
  above: every other supplied threshold gives the same floor as before.

* `sg_tables()` now builds the effect column of its GLM tables with the fit's
  own orientation: the `adverse_outcome` that `forestsearch()` stored (for a
  `grf.subg.harm.glm()` result, the value stored on it, `FALSE` by default),
  or `forestsearch()`'s rule when none is stored; it used `TRUE`. For a
  continuous fit with `adverse_outcome = FALSE` (the default) the MD column
  changes sign: in one example the selected subgroup read MD -1.094 and now
  reads 1.094, the fit's own estimate. Binary OR, RR and RD tables with
  `adverse_outcome = FALSE` (including a `grf.subg.harm.glm()` result at its
  default) now show the effect for the non-event, as the fit does. The arm
  columns and `Diff` stay raw. Survival tables, count and rate-measure (IRR,
  IRD) tables, and fits with `adverse_outcome = TRUE` are unchanged.

* `forestsearch_KfoldOut(outall = TRUE)` on a `forestsearch_Kfold()` result
  now builds the effect column of its GLM tables (ITT, full-data and fold
  subgroups) with the fit's own `adverse_outcome`, or `forestsearch()`'s rule
  when none is stored; it used `TRUE`. For a continuous fit with
  `adverse_outcome = FALSE` (the default) the MD column changes sign: in one
  example the full-data subgroup read MD -1.0941 and now reads 1.0941, the
  sign of the fit's own estimate. Binary OR, RR and RD tables with
  `adverse_outcome = FALSE` now show the effect for the non-event. The arm
  columns and `Diff` stay raw. Survival tables, count and rate-measure (IRR,
  IRD) tables, and fits with `adverse_outcome = TRUE` are unchanged.

* `compute_frontier_cis()` on a fit object that stores no `adverse_outcome`
  now orients the effect by `forestsearch()`'s rule (`TRUE` for binary and
  count outcomes, `FALSE` for continuous); it used `TRUE`. For a continuous
  outcome every estimate and limit changes sign (the limits swap), so the
  table has the orientation of a fit made at the default; in one example the
  selected subgroup's row went from -1.0941 to 1.0941, the fit's own
  estimate. A stored value is used as before; objects returned by
  `forestsearch()` store it and are unchanged.

* Change of method: the risk-difference standard error of
  `make_effect_estimator(outcome_type = "binary", effect_measure = "RD")` with
  `adjust_covariates`, when the logistic g-computation tier is used
  (`method_used = "logistic_gcomp_adj"`), is now the full delta-method SE
  `sqrt(g' V g)` of the averaged predictions over every estimable coefficient
  of the logistic fit. It was the treatment coefficient's SE scaled by the
  average derivative, which left out the other coefficients and degenerated
  when an arm had 0% or 100% events: near 0 with the treated arm there
  (control 14/40, treated 40/40, one covariate: 6.8e-06, an interval printed
  as (0.6125, 0.6125); now 0.0775) and in the hundreds with only the control
  arm there (control 0%, treated 80%, 40 per arm, adjusted for age: the
  interval around 0.8044 was (-830.78, 832.39); now (0.6838, 0.9251)); with
  both arms there it stays near 0, as before. If the new value is not a
  finite, non-negative number, or a column that treatment changes is aliased,
  the previous value is kept. The estimate, the unadjusted risk difference,
  the identity-link tier and every other measure are unchanged. The SE reaches
  the confidence limits of `forestsearch()`'s candidate table (which also
  serve as near-duplicate keys), the bootstrap's observed-estimate interval,
  the GLM tables and plots of a fit with `adjust_covariates`, and the upper
  bound `subgroup_glm()` returns; no default call reaches this tier.

* The internal `format_subgroup_summary_tables()` now highlights the
  `">=0.95"` row of the bootstrap consistency distribution, as its footnote
  says. The label was written `">= 0.95"`, which matched no row.

* `forestsearch()` with `effect_measure = "IRD"` and an `effect.threshold`
  above 1: the validator warning now speaks of a rate-difference threshold
  and suggests `effect_measure = 'IRR'` (it said "A risk difference cannot
  exceed 1.0" and suggested 'OR'), and the remap warning suggests the IRR
  scale, prints the 0.01 it assigns and gives a rate-difference example (it
  suggested the OR scale, printed "Remapping to default: 0.05" and gave a
  7 pct-point example). The risk-difference messages and every remapped
  value are unchanged.

* `sg_tables()` on a `grf.subg.harm.glm()` result with a count outcome no
  longer labels the GRF CATE note "(risk-difference scale)". It names the
  scale of the count outcome given to the forest (before the negation for
  `adverse_outcome = TRUE`), "(log(Y + 0.5) scale, Y = count)" or "(count
  scale)", adding ", not offset-adjusted" when the fit has an offset; the
  sign convention of the printed CATEs is unchanged. Binary and continuous
  results print the same note as before.

* `run_simulation_analysis()` with a GLM outcome now reports the risk
  difference (`fs_params$effect_measure = "RD"`), the risk ratio (`"RR"`) and
  the rate difference (`"IRD"`) in `hr.itt`, `hr.H.hat`, `hr.Hc.hat`,
  `hr.H.true` and `hr.Hc.true` on their own scales: the difference and the
  ratio of the two fitted arm proportions, and the difference of the two
  fitted arm rates, treatment against control. They were the log odds ratio,
  the odds ratio and the log rate ratio. Tables and plots built from these
  columns change accordingly. With an offset (the column named by
  `fs_params$offset.name`, when it is in the data) the rate difference is per
  unit of exposure, the convention of the rate ratio, while the IRD truth in
  the `ahr.*` columns and the DGM's IRD are per patient, so the two are not on
  the same scale. Odds-ratio, rate-ratio, mean-difference and survival results
  are unchanged, and stored payloads are not regenerated. (Code audit A-138,
  in part.)

* `xlog` in `plot_subgroup_results_forestplot()` never reached forestploter
  (it was passed to `forestploter::forest()`, which has no such argument),
  so the HR/OR axis was linear despite the documentation. It is now a real
  log axis via `set_xaxis(x_trans = "log")`. The call uses the forestploter
  1.2.0 API (`set_xaxis()`, `set_labs()`, `style`, `refline_gp`), which
  removes the "will be deprecated" and "Unknown arguments ignored" messages;
  `Suggests: forestploter (>= 1.2.0)`. Reported by adayim (#7).

* `ps_adjust_method = "iptw"` with a continuous outcome: the weighted
  least-squares fit failed for every candidate, because the `weights`
  argument did not resolve inside the estimator closure, so the search
  returned no candidates without an error. The weights now resolve from the
  data slice. The unadjusted and covariate-adjusted continuous fits, and the
  binary, count and survival paths, are unchanged; no committed analysis
  used this path (all are randomized-trial analyses with `is.RCT = TRUE`).
  Regression test `tests/testthat/test-iptw-continuous.R`.

### New features

* Wrappers for the analyses the articles write by hand. Each returns a plain data
  frame and no existing function changes:
  - `mrct_core_summary()`: the core operating characteristics of an
    `mrct_region_sims(record = "full")` campaign.
  - `ahr_truths_table()`, `ahr_design_table()`: the analytic layer's truths and design
    probabilities over cells and sets, with the seven consistency rows.
  - `mr_method_comparison()`, `mr_cost_table()`: naive, full-bootstrap and multiplier
    estimates side by side, method by method, and what the two corrections cost.
  - `fs_declaration_table()`: the declaration calibration at protected levels and the
    claim threshold, from one fit.
  - `mrct_record()`, `fs_in_subgroup()`, `mrct_build_dgm()`: small helpers for the
    shipped records, the declared subgroup's membership, and a DGM build that keeps its
    warnings.

* `mrct_trial_analysis()`: the regional-consistency pipeline on an observed trial,
  the counterpart of `mrct_case_study()` for real data: the all-patient and region
  Cox fits and the trigger; `forestsearch()` on the non-region patients with the
  region as the test set; the declared subgroup applied; the four-set Cox summary;
  and multiplier-resampling inference. `$stage` reports how far it went
  (`not_positive`, `no_trigger`, `no_declaration`, `complete`), and none of these
  is an error. A column of the true subgroup is optional, as in simulation.

* `mrct_four_sets()`: the four-set Cox summary from columns you name, with
  `print()`, `plot(which = c("forest", "km"))` and `as.data.frame()` methods, so
  these no longer need an `mrct_case_study` object.

* `mrct_case_study()` now runs on the same analysis code as `mrct_trial_analysis()`.
  Its arguments, result, printed output and plots are unchanged.

* `mr_values()`: a new tidy accessor for the multiplier-resampling (MR) results of
  a `forestsearch()` fit or a bare `fs_mr_inference()` result. It returns a data
  frame with one row per region (`"H"`, and `"Hc"` when the complement was
  computed) and a fixed set of columns: the naive estimate, the two-term de-biased
  estimate with its IJ interval and SE, the field and field-s one-sided bounds, the
  Bonferroni pairs, and the fit-level diagnostics (selection bias, `p_hat`, family
  size, timing). A quantity the fit did not compute is `NA`. Pure extraction:
  nothing is recomputed, apart from the optional `log()` of `scale = "log"`, and
  no existing function changes.

* The analytic layer of the MRCT regional-consistency work, moved from the
  mrct-consistency repository (`docs/theory/scripts/sk_limit.R`,
  `ahr_setup.R`, `sims/regional_entry_design.R`) in its per-patient-entry
  form. For any patient set of a Weibull AFT DGM built by
  `inject_mrct_structure()`, analysed at a calendar time with uniform
  staggered entry and administrative censoring only, it computes by
  quadrature (nothing is drawn) the exact limit of the two-sample Cox
  estimate, its information and robust variance, the expected events, the
  exact correlation between nested sets' estimates, the set's AHR and CDE,
  and the normal-approximation trigger probabilities at a trial size:
  * `ahr_population()`: the population of a DGM, with an optional
    per-patient or per-group entry delay (`entry_delay`);
  * `ahr_limit()`: the limit and its companions for a named or logical set;
  * `ahr_design()`: SEs, events, correlations and `p_itt_ub_lt1`,
    `p_region_itt_gt1`, `p_region_Gstar_lt1`, `p_trigger`,
    `p_trigger_and_Gstar_lt1` at a trial size `n`;
  * `ahr_solve()`: the scan-and-`uniroot()` solve of P(trigger) = target in
    one free parameter, returning every bracketed root.
  With no delay the layer reproduces the committed mrct-consistency records
  to below 1e-14.

* `scale_dgm_hazard(dgm, m)`: multiplies every hazard of a DGM by `m`
  (`mu <- mu_fitted - tau * log(m)`, hazard ratios unchanged), recording `m`
  in `dgm$mrct$hazard_multiplier` and `attr(dgm, "m")`.

* `simulate_from_dgm()` gains `entry_delay` (default `NULL`, the current
  behaviour): `list(var = <column>, delay = a)` makes patients with
  `var == 1` enter at `a + U(0, max_entry)`, redrawn every call; a numeric
  vector gives each sampled patient's delay. `entry_var` keeps precedence.
  `mrct_region_sims()` passes it through `sim_args`; its documentation now
  states the `entrytime` column convention.

* `expand_covariates()` and `inject_mrct_structure()` gain `latent_R`
  (default `NULL`, the current behaviour): a supplied latent correlation
  matrix is validated and used as given, with no latentcor fit, so a DGM can
  be rebuilt identically on any platform from a stored `dgm$mrct$latent_R`.

* `mrct_region_sims()` gains `record = c("point", "full")`. `"point"` (the
  default) returns exactly the current table. `"full"` keeps the
  per-replicate record of the same trials (nothing is redrawn): the table
  gains the 2 x 2 counts of Ĥ membership against `flag_harm` over the
  non-region and the region patients (`TP_nonregion`, ..., `TN_region`), and
  carries `attr(res, "estimators")`, one row per trial and patient set
  (`itt_all`, `itt_nonregion`, `itt_region`, `region_Ghat`, `region_Gstar`,
  `Hstar_nonregion`) with the unadjusted Cox fit (`hr`, Wald `lo`/`hi`,
  `se_loghr`), patients and events by arm, the set's `ahr` and `cde`, and
  `declared`. It reproduces the redraw-based records of mrct-consistency
  (`sims/123_ident_rlog4_full_records.R`) exactly.

* `mrct_region_sims(record = "full")` adds, after the counts, the trial's
  multiplier-resampling bounds from `fs$mr_inference` when `fs_args`
  carries `mr_inference = TRUE`: `mr_ci_method`, `mr_debiased_H`,
  `mr_ij_lower_H`, `mr_ij_upper_H`, `mr_field_lower_1s_H`,
  `mr_debiased_Hc`, `mr_ij_lower_Hc`, `mr_ij_upper_Hc`,
  `mr_field_upper_1s_Hc`, `mr_field_upper_1s_s_Hc`, `mr_bonf_lower_H`,
  `mr_bonf_upper_Hc`, `mr_p_hat_H`, `mr_selection_rate`; `NA` where the
  element was not computed (no MR, no declaration, complement or field
  options off, `return_reselection = FALSE` for `mr_p_hat_H`).
  `record = "point"` is unchanged.

* `fs_attach_betaHhat()`, `fs_betaHhat_table()` and `fs_betaHhat_one()`
  gain `method = c("frame", "exact")` (with `dgm`, `analysis_time`,
  `max_entry`, `cens_adjust`). `"frame"` (the default) is the current route.
  `"exact"` (survival only) computes beta(Ĥ) and beta(Ĥᶜ) without a draw,
  from `ahr_population()` / `ahr_limit()` on the realised membership in
  `dgm$df_super`: the same population quantity as the frame route, without
  the Monte Carlo error of its single evaluation draw. It covers
  administrative censoring only and stops when `cens_adjust` leaves the
  DGM's random censoring in place (as `0` and the frame default `log(1.5)`
  do; `50` or `Inf` remove it).

* Two helpers find and examine a case-study trial of an MRCT simulation:
  * `mrct_select_case(res)`: on a `mrct_region_sims(record = "full")`
    result, flags each trial for the trigger (ITT upper 95% limit <
    `itt_ub` and region HR > `region_itt`), a declared subgroup, and
    region(Ĝ) < 1.0, and returns the first trial meeting the chosen
    conditions (`which = "first"`, the rule that chose replicate 5 of S2 in
    the mrct-consistency report) or all of them, with the trial's hazard
    ratios. A missing HR never satisfies a condition.
  * `mrct_case_study(dgm, seed, k, ...)`: redraws trial `k` exactly as
    `mrct_region_sims()` drew it (`simulate_from_dgm(seed = seed + k)` under
    L'Ecuyer-CMRG, the caller's generator restored), runs `forestsearch()`
    on the non-region patients with the run's arguments and multiplier
    resampling (seeded `seed + k`), and returns the four sets (all, region,
    region(G*), region(Ĝ)) with their Cox fits, Ĥ membership, the 2 × 2
    tables of Ĥ against H* and the MR result. `print()` and `plot()`
    (forest plot; Kaplan–Meier panels with weightedsurv) methods;
    `mr_estimates_table(x$fs)` gives the inference table.

* `inject_mrct_structure()` has an observed-region mode:
  `region = list(observed = "<col>")` takes the region from a 0/1 column of
  `seed_data`, expands it with the covariates (no region model is fitted or
  drawn), and estimates its prognostic log-HR by the AFT fit unless
  `region$loghr` fixes it. The planted-model elements (`prevalence`,
  `or_pred`, `x1_vars`, `or_x1`, `scale`, `band`) are errors in this mode;
  `dgm$mrct$region$observed` records the mode. Defaults are unchanged.

* Region-dependent treatment effect. `inject_mrct_structure()` gains
  `region_loghr_shift` and `region_harm_loghr_shift` (default `0`, the
  current build exactly): the treated log hazard of region patients is
  shifted by γ_R, plus γ_H\*R for region patients in
  H\* (the `x3` subgroup), so with a base effect of `b` outside H\* and `h`
  inside it the region hazard ratios are b_R = b·exp(γ_R) and
  h_R = h·exp(γ_R + γ_H\*R). The shift is applied to
  `theta_1`, `lin_pred_1`, `loghr_po` and treated `lin_pred_obs`; the control
  arm is unchanged; non-zero shifts are recorded in `dgm$mrct$region_effect`.
  `ahr_limit()` and `ahr_design()` gain `b_R` and `h_R` (default `b`, `h`):
  the four-cell (H\* × region) parametric effect, so the analytic layer's
  (`b`, `h`) route matches the extended DGM. `ahr_design()` returns `b_R`
  and `h_R`.

* `mr_estimates_table()` no longer needs a bootstrap template:
  `boot_results` now defaults to `NULL`, and the table is then built from
  the `forestsearch(mr_inference = TRUE)` result alone. The descriptive
  columns come from `fs$df.est` (by `SG_tab_estimates()`, as the bootstrap
  lays them out), the naive effect is `fs$mr_inference$naive` with its
  robust limits (the header says "robust"), followed by the MR de-biased
  estimate with its IJ interval and, when the field quantities exist, the
  field one-sided 95% bound (lower for the selected subgroup, upper for its
  complement) and the Bonferroni pair; complement cells carry MR values only
  under `include_complement = TRUE`. The returned table carries the numbers
  it displays as `attr(, "mr_values")`. A call that supplies
  `boot_results` returns the same numbers; its bias-corrected footnote now
  names the method used (the MR de-biased estimate with the number of
  draws, and its interval from the IJ variance on the multiplier draws, or
  the robust SE under `ci_method = "wald"`) instead of the bootstrap's
  "infinitesimal jackknife method (draws iterations)". `format_bootstrap_table()`
  gains `bc_footnote` (default `NULL`, the bootstrap footnote) for this.

* The lean unadjusted Cox fit lives in the internal helper
  `.fs_cox_two_group()` and is guarded by a once-per-session probe: if
  `coxph.fit()` does not reproduce `coxph()` on a near-tie fixture (for
  example after a change in a future survival), the search-stage fit falls
  back to the formula interface for the session, with one message. Results
  are identical either way (`tests/testthat/test-lean-fits.R`; record
  `dev/reports/REPORT_robustness_041_2026-10-28.md`).

* Two per-candidate search-stage fits now call the fitting engines directly
  instead of going through the formula interface: the unadjusted survival
  fit in `fit_cox_for_subgroup()` calls `survival::coxph.fit()` with exactly
  the arguments `coxph()` passes, and the unadjusted binary odds-ratio fit
  runs `stats::glm.fit()` through a tagged closure and vector path (the twin
  of the continuous fast path). Results are bit-identical (identity tests in
  `tests/testthat/test-lean-fits.R`; verification record
  `dev/reports/REPORT_lean_fits_2026-10-26.md`).

* Measured per-replicate time of a sequential `nb_boots = 20` bootstrap on
  the committed applications: survival (gbsg) 3.60 s -> 2.02 s (-44%),
  binary OR (ACTG175) 5.93 s -> 2.75 s (-54%), continuous (ACTG175)
  unchanged (1.87 s -> 1.84 s).

* Nothing changes for adjusted, weighted (propensity), RD or RR fits: they
  keep their formula paths.

### Superseded and deprecated

* `setup_gbsg_dgm()`, `calibrate_k_inter()` and
  `generate_gbsg_bootstrap_general()` are now documented as superseded (a
  "Superseded" section on each help page; no run-time message). They are kept
  so that published simulations reproduce and will be deprecated in a future
  release. Replacements: `generate_aft_dgm_flex()` with `simulate_from_dgm()`,
  `find_k_inter_for_target_hr()`, and `generate_bootstrap_synthetic()`
  (whose first example returns the same data as
  `generate_gbsg_bootstrap_general()` with the same `n`, `seed` and
  `noise_level`).

* The censoring step of the legacy GBSG wrapper is documented on
  `setup_gbsg_dgm()`, with a note on `fs_fdr_report()` and
  `validate_k_inter_effect()`: the Weibull censoring model is fitted to the
  resampled super-population with the arm assigned for simulation, and its
  covariate term is the fitted term multiplied by the ratio of the event-model
  scale to the censoring-model scale. `generate_aft_dgm_flex()` uses the
  censoring fit's own term. The true hazard ratios and the `k_inter`
  calibration do not depend on it; simulated censoring times do. The step is
  unchanged so that published simulations reproduce.

* `vi.grf.min` and `max_n_confounders` in `forestsearch()` are deprecated.
  Both formals are kept (same positions, same defaults), so no existing call
  breaks, but a non-default value now emits one deprecation warning and has
  no effect. The GRF variable-importance screening/ordering of the candidate
  cut columns is removed. The candidate family is now the cut grid alone,
  and no model is fitted in the forest-search path unless `use_lasso`,
  `use_grf` or `use_dina` is `TRUE`.

* `max.minutes` in `forestsearch()` and `subgroup.search()` is deprecated.
  It was documented as inert ("scheduled for deprecation in v0.3.0") and
  that was never actioned. Both formals are kept (same positions, same
  defaults); the argument has no effect. `forestsearch()` emits one
  deprecation warning for a non-default value, and `subgroup.search()` one
  warning whenever it is passed. The formals are removed in a later release.

* What changes for callers who passed a non-NULL `vi.grf.min`: at values
  `<= 0` the screen never filtered anything, so the selected subgroup's
  membership is unchanged. The clause order of the `sg.harm` label (and the
  redundancy column-walk) can still differ. A binding value `> 0` no longer
  screens.

* `run_simulation_analysis()`'s default argument list, the test helper
  `.fs_args_for()` and the `forestsearch` vignette no longer pass either
  argument. `fs_family_report()` reports both on a single
  `"GRF variable-importance screening (deprecated)"` row with status
  `inert`, echoing the values that were passed.

* Both formals will be removed in a later release.

* The internal helper `fit_causal_forest_glm()` is removed (unused since
  the `vi.grf.min` deprecation).

### Documentation

* Help pages corrected to say what the code does (no behaviour change):
  - `generate_bootstrap_with_noise()`: by default every numeric column except
    `id_var` is treated as continuous (by design), so a numeric-coded category
    such as a 0/1 indicator gets noise, rounded back to its value at a small
    `noise_level`, unless `continuous_vars` is given without it; the
    `cat_vars` auto-detection and the `grade` special case are stated
    exactly.
  - `generate_bootstrap_synthetic()`: `ordinal_vars` only marks which of the
    `cat_vars` are ordinal; a variable listed only there is kept as
    resampled; a factor listed in both becomes NA at a changed row unless its
    levels are "1", "2", ...; and at the lowest or highest level of a
    numeric-coded variable the new value is not always the adjacent level.
  - `generate_gbsg_bootstrap_general()`: `grade` is resampled, not perturbed.
  - `forestsearch()`, `dina_subgroup()`, `dina_subgroup_refit()`: DINA's
    Poisson quantity is the log mean-count ratio (DINA takes no exposure
    offset), which equals log(IRR) only when mean follow-up does not depend
    on arm. The printed `"log-IRR"` label is unchanged.
  - `forestsearch()`: `dmin.grf` is a harm threshold on GRF's doubly robust
    scale, not a number of events; under `grf_select_statistic = "effect"` it
    is an existence gate and is then not applied per candidate.
  - `forestsearch()` and `subgroup.consistency()`: a sequential and a
    parallel plan can give different consistency proportions, and possibly a
    different selection, whenever literal splitting runs
    (`consistency_method = "split"`, an IRD effect under `"resample"`, or a
    candidate whose closed form is not available), except under
    `sg_focus = "maxeff"`.
  - `fs_fdr_report()`: `k_inter` is calibrated on the builder's default
    super-population size, not on `n_super`; the calibrated value does not
    depend on that size beyond the root-finding tolerance.
  - `setup_gbsg_dgm()` documents its own arguments (they were inherited from
    the removed `create_gbsg_dgm()`), and the examples that called
    `create_gbsg_dgm()` now call `setup_gbsg_dgm()`.
  - Examples load `survival::gbsg` directly instead of
    `data(gbsg, package = "survival")`, which loads nothing.

* Documentation pass (part 5 of the CRAN preparation; no behaviour change):
  - The display help pages (`sg_tables()`, `forestsearch_KfoldOut()`,
    `cv_summary_tables()`, `plot_sg_glm_outcomes()`,
    `plot_subgroup_results_forestplot()`, `mr_estimates_table()`) say that the
    effect column follows the fit's orientation while the raw arm columns
    (and `Diff`, where shown) do not, so for an unadjusted fit where a lower
    outcome is adverse the effect and the raw difference can have opposite
    signs. `?make_effect_estimator` says what `adverse_outcome` does for each
    outcome type.
  - `?forestsearch`: a supplied risk- or rate-difference threshold above 1
    is remapped, with two warnings for `effect.threshold` (one for
    `consistency.threshold`), and the remapped
    effect threshold is also the DINA floor where DINA selects a cut; under
    the default `dina_select_statistic = "effect"` the `use_dina` screening
    cut can differ from the one `subgroup_method = "dina"` would select; what
    `consistency_method` does under propensity adjustment.
  - `CV_sgs()` defines "Exact" for a one-factor subgroup;
    `summarize_bootstrap_subgroups()` says its bins are closed on the left;
    `test_hte_from_forestsearch()` describes how cut labels map to columns and
    its warnings; `calibrate_k_treat()`, `find_k_inter_for_target_hr()`,
    `compute_dgm_cde()` and `generate_bootstrap_with_noise()` are corrected to
    the code; `setup_gbsg_dgm()` carries the model description of the removed
    `create_gbsg_dgm()` page.
  - Sixteen `\dontrun{}` examples are corrected so that they run as written
    after `library(forestsearch)`; fifteen were run to the end, and
    `fs_fdr_report()`'s, a long simulation, was not. They used, among other
    defects, a bare `gbsg`, the missing `grade3` and id columns, objects or
    list elements that do not exist, a wrong or placeholder call or one
    missing a required argument, a nonexistent `nsim` argument, or
    `n_workers = 100`.
  - The "Biomarker Treatment Effects" article names `setup_gbsg_dgm()`. The
    simulation article no longer calls its `generate_aft_dgm_flex()` block
    equivalent to `setup_gbsg_dgm()`, says how the two differ, and says that
    `setup_gbsg_dgm()` is used there to reproduce the published simulations
    and is superseded for new work.
  - Nine test files of part 4's fixes carry the comments their verifiers
    asked for (how expected values are derived, and which checks pass before
    and after the fix); comments only.

* New test file `test-synthetic-data-documented-behaviour.R` pins what the
  synthetic-data help pages now say: `grade` is resampled but not perturbed by
  `generate_gbsg_bootstrap_general()`, which returns the data of
  `generate_bootstrap_synthetic()`'s first example call; under auto-detection
  `generate_bootstrap_with_noise()` leaves a numeric 0/1 column as resampled and
  perturbs a factor; named in `cat_vars` with `continuous_vars` given, the 0/1
  column is perturbed.

* Test-suite reconciliation of the estimability-boundary invariant with
  `glm_effect_estimators` behaviour; no functional change. The source-level
  test in `test-estimability-boundary.R` counted exactly one
  non-convergence catch. A second, verbatim catch was added with the
  OR-only fast closure `.make_glm_or_estimator_fast()` (the lean binary OR
  fit), and the test now asserts that structure: one OR/RR-scoped catch in
  the slow closure, one in the fast closure, and neither reachable by RD.

* New getting-started vignette, `vignette("forestsearch")`: the German Breast Cancer Study Group analysis of
  León and Anderson (2026b), Table 1, run end to end — the three identifiers (forest search, DINA, GRF) and
  the selection-adjusted estimates and bounds of León and Anderson (2026a) — at the manuscript's settings,
  with the declaration calibration at a protected level.

* The previous getting-started vignette, which documents the workflow of León et al. (2024), is now the
  pkgdown article "Exploratory subgroup identification with forestsearch (the 2024 workflow)" and is no
  longer built on CRAN.

* The post-selection analysis is now an article on the package website rather than a CRAN vignette, because the full analysis takes a few minutes to run: <https://larry-leon.github.io/forestsearch/articles/survival-post-selection.html>. It reproduces the GBSG analysis of León and Anderson (2026a), illustrates the full bootstrap at a small number of replicates, and shows `mr_estimates_table()`.

* The methodology and simulation articles on the package website keep their text as the record of the 2024 workflow, with the deprecated arguments removed from their listings and a dated note at each point where the current method differs. The articles on treatment-effect definitions, biomarker effects and the causal interpretation of hazard ratios now state the estimand the package reports for an identified subgroup, as in León and Anderson (2026a, 2026b). All five articles are now Quarto documents.

* New article on the package website, "Post-selection bounds on a nested continuum: the Guo and He (2021) design": Section S4 of León and Anderson (2026b), with the complete code, a small default run, and the paper's results beside it.

* New article on the package website, "Binary and continuous outcomes: the ACTG 175 analyses": Section S3 of León and Anderson (2026b) — forest search on the binary and continuous endpoints, with the selection-adjusted bounds and the declaration calibration — and Table S3 of León and Anderson (2026a), the three identifiers on the binary endpoint with the full bootstrap, at a small number of replicates, beside multiplier resampling.

* The README describes the current method and shows the analysis of `vignette("forestsearch")`.

* `?forestsearch`: Example 1 passes a bootstrap result to `mr_estimates_table()`, which needs one as its template (the example called it without one); Example 2 names the analysis it runs, León and Anderson (2026a, Table S3); a section reference is corrected.

* The `run_null_calibration()` example runs at an illustrative size (two trial sizes, ten simulations, no LASSO/GRF screening); a calibration for use takes more simulations, as the example now says.

# forestsearch 0.4.0

## Changes since the last CRAN release (0.1.0)

### Breaking changes

* `showten_subgroups` is removed from `forestsearch()` and
  `subgroup.consistency()`. Neither function has `...`, so a call that passes
  it now fails with "unused argument". Use `show_candidate_summary` instead.
* Changed defaults in `forestsearch()` (old -> new):
  * `use_lasso`: `TRUE` -> `FALSE`
  * `use_grf`: `TRUE` -> `FALSE`
  * `dmin.grf`: `12` -> `0`
  * `frac.tau`: `0.6` -> `0.8`
  * `d0.min`: `12` -> `10`
  * `d1.min`: `12` -> `10`
  * `max_subgroups_search`: `10` -> `Inf`
  * `vi.grf.min`: `-0.2` -> `NULL`
  * `stop_threshold`: `0.95` -> `pconsistency.threshold`
  * `parallel_args`: `list(plan = "callr", workers = 6, show_message = TRUE)`
    -> `list(plan = "multisession", workers = .default_parallel_workers(),
    show_message = TRUE)`
* Changed default in `subgroup.consistency()`: `stop_Kgroups`: `10` -> `200`.
* Changed default in `forestsearch_Kfold()`: `Kfolds`: `nrow(fs.est$df.est)`
  -> `.fs_cv_nobs(fs.est)`.
* Changed default in `cox_ahr_cde_analysis()`: `output_dir`: `tempdir()` ->
  `"plots"` (files are written only when `save_plots = TRUE`).
* `save_leff_calibration()` now requires `output_dir`. Called without it, the
  function stops instead of writing to `"_data"` in the working directory;
  pass `output_dir = "_data"` to keep the former behaviour.

### New features

* GLM outcomes: binary (odds ratio), continuous (mean difference) and
  count/rate (incidence rate ratio) endpoints, selected with the
  `outcome_type` and `effect_measure` arguments of `forestsearch()`.
* DINA and GRF as subgroup identifiers alongside the consistency search,
  selected with `forestsearch(subgroup_method = c("consistency", "dina",
  "grf"))`; see `dina_subgroup()` and, for GLM outcomes, `grf.subg.harm.glm()`.
* Post-selection inference: multiplier resampling (`fs_mr_inference()`),
  with infinitesimal-jackknife (IJ) standard errors and field bounds.
* Operating-characteristics functions: `fs_oc_family_enumerate()`,
  `fs_oc_predict()`, `fs_oc_grid()` and `fs_oc_invert()`.
* Threshold rules on the forest-search consistency path (survival HR and
  binary OR only): when only `c1` is given, `c2 = 0.8 * c1` on the ratio
  scale; `c2 > c1` is an error.
* Calibrated declaration threshold: `fs_declaration_calibration()`.
* `forestsearch()` gains `pconsistency.digits` (default 2): the consistency
  proportion is rounded to this many digits before it is compared with
  `pconsistency.threshold`.
* grf forest training is capped at 2 threads when `_R_CHECK_LIMIT_CORES_` is
  set (as under `R CMD check --as-cran`); otherwise grf's own default applies
  unchanged.
* Printed post-selection output and documentation simplified; methodology
  reference León and Anderson (2026) <doi:10.48550/arXiv.2609.38361>.

## Development changes after 0.3.5 (not released on CRAN)

* **MR's admission set now uses the screen's rounded consistency threshold
  (forest search only).** The search admits on
  `round(Pcons, pconsistency.digits) >= p*`, but MR's consistency floor was
  built from the exact cutoff `qnorm((1 + p*) / 2)`, so MR re-selected under a
  different admission rule from the one that produced the selected subgroup,
  in the under-adjusting direction. The floor is now built from the effective
  threshold `.fs_pcons_eff(p*, digits)` on both consistency paths.
  `fs_mr_inference()` gains `pconsistency.digits` (default `NULL`, the
  `subgroup.consistency()` default); `forestsearch()` passes its own value.
  **Corrected estimates, field bounds and the Bonferroni pair can move for
  forest search.** On the GBSG application (seed 8316951) the corrected HR
  moves from 1.2806 to 1.2762, the field lower bound from 0.6365 to 0.6336,
  and the Bonferroni pair from 0.5099 / 0.8397 to 0.5087 / 0.8396. GRF and
  DINA carry no consistency screen and are unaffected (identical output).

* **`fs_declaration_calibration()`: FW-hat and the implied p-star are now
  computed at the screen's rounded admission rule.** The screen admits on
  `round(Pcons, pconsistency.digits) >= p*`, which is `Pcons >= 0.895` for
  `p* = 0.90` at 2 digits, but `fw_size` was taken at the exact
  `qnorm((1 + p*) / 2)`. It is now taken at the effective threshold
  `qnorm((1 + pcons_eff) / 2)` (returned as `pcons_eff` and `z_pstar`), with
  `digits` read from the fit's `args_call_all` (default 2, source recorded), so
  **`fw_size` values rise relative to previous releases** and now depend on
  `pconsistency.digits`; `admitted_pstar` uses the same threshold. On the GBSG
  application (`c0 = 0.75`, `alpha = 0.10`) `fw_size(c0)` moves from 0.6508 to
  0.6672. **`kappa_hat` is unchanged.** The c0 table's `pstar_implied` column
  is replaced by what to set to run `kappa_hat` as a p-star screen:
  `pstar_settable` (smallest settable p* at the fit's digits whose rounded
  screen is at or above `kappa_hat`; `pstar_achievable = FALSE` when none is),
  `pcons_eff_settable`, `z_eff_settable`, `z_gap`, and `digits_fine` /
  `pstar_fine` / `z_gap_fine` (smallest digits with a gap under 0.01 in z).
  Under `consistency_method = "split"` the screen is not a threshold on `T`,
  so `fw_size`, `admitted_pstar` and the settable columns are `NA` / `NULL`.

* The `fs_dgm_feasibility()` example now runs (it called `generate_glm_dgm()`
  without its required data arguments); `fs_plot_bias_coverage()` captions use
  `\u` escapes, so the R source is ASCII -- the rendered text is unchanged.

* **`forestsearch()` gains `pconsistency.digits` (default 2).** It was
  previously reachable only on the internal consistency functions
  (`subgroup.consistency()`, `evaluate_subgroup_consistency()`,
  `evaluate_consistency_twostage()`); every `forestsearch()` run was pinned to
  2 digits. The consistency proportion is rounded to this many digits *before*
  it is compared with `pconsistency.threshold`, so the setting governs
  admission, not only display; the roxygen for all four functions now says so.
  Default behaviour is unchanged: at 2 the GBSG application's consistency table
  is identical to the previous package's. `fs_family_report()` lists the
  argument in its consistency-screen row.

* **New `fs_dgm_feasibility()`: a design-time check that the planted region
  can be declared at all.** Given a `glm_dgm` -- the interface
  `fs_oc_predict()` and `fs_oc_grid()` accept -- it draws `n_rep` replicates
  per sample size through the DGM's *own* generator (`simulate_from_glm_dgm()`,
  the function the campaign templates call, never a re-implementation) and
  reports, per `n`: the share of replicates on which the planted region is
  **undeclarable** (size `<=` `n.min`, the strict test `subgroup.search()`
  itself applies), the share under the per-arm events floor, the share on
  which the estimand does not exist under the estimability boundary's
  per-estimand condition, and size and per-arm cell summaries. `$feasible` is
  `TRUE` when every undeclarable share is at or below `tolerance` (default
  0.05). The print method is loud and unconditional. On the ACTG175 OR 0.75
  design this returns `FALSE` at n = 500 and n = 750: the planted region is at
  or below `n.min = 60` in 90% of replicates at n = 500, so what the search
  declares there is a larger overlapping region, not the planted one.
  **It imposes nothing** -- the floors are read and reported, never applied,
  and `forestsearch()` never consults it. The function does not change the RNG
  kind and restores the caller's RNG stream on exit; the caller's matching
  obligation, stated in the roxygen, is to build and calibrate the DGM
  *before* any replicate switches the kind, since calibrating afterwards
  yields a different super-population. **Survival DGMs are not accepted**, and
  the error and `?fs_dgm_feasibility` say exactly what supporting them would
  need rather than guessing: `setup_gbsg_dgm()` objects carry a different
  generator (`simulate_from_dgm()`) whose `analysis_time` and `cens_adjust`
  arguments have no GLM counterpart and change the per-arm event counts.

* **The estimator boundary now returns `NA` with a reason where the estimand
  does not exist, instead of a finite divergent number.** `glm()` and
  `coxph()` do not fail on a slice with an empty cell or a zero-event arm:
  they return a large finite coefficient with `converged = TRUE`, which was
  then admitted as an estimate (a zero-cell OR returned 19.97 with SE 4809).
  The existence condition is now checked before the fit, per estimand: **OR**
  requires all four cells (control/treated x events/non-events) `>= 1`;
  **RR**, **IRR** and **HR** require at least one event in each arm. A
  candidate failing it returns `estimate = NA`, `se = NA`,
  `converged = FALSE` and a `reason` naming the empty cell, and takes the
  existing fit-failure status (status 5) -- an `NA` effect cannot rank, so no
  selection logic changes. For the ratio estimands only, a fit reporting
  `converged = FALSE` is likewise non-estimable (`"non-convergent fit"`).
  **Unaffected, byte-identical: RD, IRD and MD**, which exist on any slice --
  in particular RD's tier-3 raw-proportions fallback returns
  `converged = FALSE` by design and is never caught. **This is not an
  admission floor:** `n.min`, `d0.min`/`d1.min` and every DINA/GRF floor are
  unchanged, nothing is forwarded to DINA or GRF, and no minimum count is
  imposed anywhere. Non-estimable candidates are counted, not dropped
  silently: `filter_counts$n_nonestimable` and
  `filter_counts$nonestimable_reasons` carry the tally and its reasons,
  `subgroup.search(details = TRUE)` prints them, and `fs_family_report()`
  gains an "estimability boundary" stage row. Committed results are
  unaffected -- no declared subgroup in the committed binary record violates
  the four-cell condition (3,313 replicates, `a3374e6f`), and no committed
  survival bundle stores a non-finite estimate or SE (930 bundles, 568,757
  declared rows).

* **`dina_frontier()`'s display caps now default to `Inf` and say when they
  trim.** `max_per_covariate` and `max_subgroups` defaulted to `3L` and `10L`,
  which silently dropped covariates whole and truncated the round-robin
  rank-2 pass -- observed in the applied vignettes, where the printed frontier
  was not the frontier. Both now default to `Inf`, aligning with
  `max_subgroups_search`'s own `Inf` doctrine, and a user-supplied finite value
  that actually removes rows raises a warning naming the kept and available
  counts (condition class `dina_frontier_cap_trim`). **Affected: callers who
  relied on the finite defaults** -- they now receive the full table; every
  committed in-repo caller passes both caps explicitly, so no committed
  rendered table changes. **Unaffected: the candidate pool anywhere.** The
  frontier is a display object; `dina_subgroup()` never consults it, and
  `forestsearch(use_dina = TRUE)` still supplies the finite screening caps
  itself, so the screening pool is byte-identical. Under the default
  `selected_only = TRUE` that screening frontier is discarded, so the trim
  warning is muffled there rather than naming a table nobody sees.

* **`dina_args` frontier keys now warn when they are ignored.** The seven
  `dina_frontier()` keys -- `scope`, `m_diff`, `n_min`, `direction`,
  `max_per_covariate`, `max_subgroups`, `digits` -- are inert under
  `subgroup_method = "dina"`, which selects through `dina_subgroup()` and
  never consults a frontier table. Supplying any of them there now raises a
  single `warning()` per fit naming every offending key at once, rather than
  accepting them silently. **Affected: `subgroup_method = "dina"` callers who
  pass frontier keys** -- the selection is unchanged, only the messaging.
  **Unaffected:** the `use_dina` screening path, where the keys still act
  (`digits` under `selected_only = TRUE`, all seven under
  `selected_only = FALSE`), and every other `subgroup_method`.

* **DINA now refuses the identity-scale estimands `"RD"` and `"IRD"`.** DINA
  scores candidates on the family link scale, so its admission floor was
  derived as `log(hr.threshold)` for the ratio families and as the identity
  threshold for `gaussian`. For an explicit `effect_measure = "RD"` that
  derivation was meaningless: a 7-percentage-point threshold became a log-odds
  floor of `log(0.07) = -2.66`, which admitted nearly every candidate,
  silently. There is no identity-scale number DINA could correctly apply on a
  link scale, so the derivation now stops instead, naming the ratio and
  non-DINA alternatives. The guard sits at the `m_diff` derivation itself, so
  both routes to it are covered: `subgroup_method = "dina"` (an error, raised
  before the DINA model is fit) and `use_dina` screening with
  `dina_args$selected_only = TRUE` (where it surfaces through that path's
  existing `"DINA analysis failed"` warning and the run continues with no DINA
  cuts). **Affected: only calls that explicitly request `"RD"` or `"IRD"`
  together with DINA** -- nothing committed does, and the defaults per outcome
  type (`"OR"`, `"MD"`, `"IRR"`) are all outside the guard. **Unaffected:**
  `family = "gaussian"` with `"MD"`, which is a correct identity-scale floor
  and cannot reach the refusal, and every ratio estimand on every family.

* **The two effect thresholds are now a pair on the forest-search consistency
  path.** For `subgroup_method = "consistency"` with survival (`HR`) or binary
  `effect_measure = "OR"` -- the estimands whose `c1` and `c2` sit on one
  comparable ratio scale -- `forestsearch()` now enforces `c2 <= c1` and
  derives a `c2` the caller did not supply. **Affected: callers who set `c1`
  without `c2`** (`effect.threshold` or `hr.threshold` alone), which used to
  leave `c2` at the package default `1.0` whatever `c1` was and now resolves
  `c2 = 0.80 * c1` on the ratio scale -- 1.25 -> 1.00, 1.00 -> 0.80,
  0.90 -> 0.72 -- announced with a message and carried into every bootstrap
  replicate and cross-validation fold. The default pair (1.25, 1.0) is the
  rule's fixed point, so a call that set neither threshold resolves exactly
  what it did before, and a call that passed `hr.threshold = 1.25` explicitly
  is likewise unmoved and only gains the message. **Also affected: callers
  with `c2 > c1`**, previously accepted and degenerate -- the consistency-stage
  entry condition re-screens the admitted family on `c2`, so candidates the
  screen let through were discarded at the stage boundary and the run could
  report no subgroup for a reason nothing in the output named. That
  configuration is now an error, raised before any model is fit, naming both
  values and the spellings the caller used. Supplying both spellings of one
  threshold at disagreeing values is an error too. An explicitly supplied `c2`
  is never overridden, and no selection logic changed. **Unaffected:**
  `subgroup_method = "dina"` and `"grf"`, which return before the consistency
  stage and never consult `c2`, and the identity-scale estimands `"RD"`,
  `"IRD"` and `"MD"` together with `"IRR"`, all of which resolve exactly as
  before in every branch.

* **Bootstrap replicates and cross-validation folds now resolve the thresholds
  the original fit resolved.** `forestsearch()` resolves the screening (`c1`)
  and consistency (`c2`) thresholds through branches guarded by whether the
  caller supplied them, and for the legacy spellings (`hr.threshold`,
  `hr.consistency`) that was detected with `missing()`. Because
  `forestsearch_bootstrap_dofuture()`, `forestsearch_Kfold()` and
  `forestsearch_tenfold()` replay the captured argument list with every formal
  supplied, `missing()` was `FALSE` inside a replicate, and a replicate could
  resolve thresholds its parent never used. **Affected: identity-scale
  analyses (`effect_measure` `"RD"`, `"IRD"` or `"MD"`) that left `c1` or `c2`
  at its default, or supplied only one of the two, and then ran the bootstrap
  or cross-validation** -- those replicates used `c2 = 1.0` where the fit used
  `0.0`, and on `"MD"` also `c1 = 1.25` where the fit used `0.0`. A per-split
  effect of 1.0 is a threshold no candidate can meet, so the affected folds
  and replicates identified nothing. Ratio measures (`"OR"`, `"RR"`, `"IRR"`)
  and survival were never affected, and the fix leaves every original fit's
  resolution byte-identical -- only what a replay sees changes. The resolved
  values are now written back into `effect.threshold` /
  `consistency.threshold`, which are detected by value rather than by
  `missing()` and so survive any wrapper. A spurious
  "`effect.threshold` appears to be on a ratio scale" warning that fired once
  per replicate on `"RD"` / `"IRD"` also stops.

* **The estimation-layer entry points default binary to the odds ratio.**
  `make_effect_estimator()` and `consistency_resample()` resolved an unset
  `effect_measure` to `"RD"` for `outcome_type = "binary"` while
  `forestsearch()` resolves it to `"OR"`; both now resolve `"OR"`. Continuous
  (`"MD"`), count (`"IRR"`) and survival (`"HR"`) defaults are unchanged, and
  an explicitly supplied measure is untouched. No in-package caller reached
  the old default -- `forestsearch()` passes `effect_measure` explicitly -- so
  this changes results only for a direct call that omitted it.

* **Binary outcomes now default to the odds ratio.** `forestsearch()` with
  `outcome_type = "binary"` and `effect_measure` left unset resolves to
  `"OR"`; it previously resolved to `"RD"`. **This changes results for binary
  callers who inherit `effect_measure`**: because the screening and
  consistency thresholds are resolved from the estimand, such a run moves
  from `c1 = 0.05` / `c2 = 0.0` on the identity scale to `c1 = log(1.25)` /
  `c2 = log(1.0)` on the log scale, and the search screens candidates on the
  log odds ratio rather than the risk difference. To keep the old behaviour,
  pass `effect_measure = "RD"` explicitly. **Callers who already pass
  `effect_measure` -- any value, including `"RD"` -- are unaffected**, as are
  survival (`HR`), continuous (`MD`) and count (`IRR`) outcomes, whose
  defaults are unchanged. No selection logic changed. The default also aligns
  the three identifiers: `subgroup_method = "dina"` already screened on the
  log odds ratio (its harm floor is `log(hr.threshold)` for a binomial
  family), so this brings the consistency and GRF identifiers onto the same
  estimand; DINA's own behaviour is unchanged.

* **GRF membership on factor covariates.** GRF's subgroup evaluator compared a
  factor covariate's raw column with a numeric cut, so every candidate on a
  factor covariate got `NA` membership and was dropped from the effect
  re-selection and from the multiplier-resampling family; it now codes each
  covariate exactly as the forest's covariate matrix does, through one shared
  internal helper.

* **DINA proposal floor orientation.** With `subgroup_method = "dina"` on a
  continuous or binary outcome and `adverse_outcome = FALSE`, DINA's proposal
  floor was applied to the raw tau-hat, so it admitted benefit while the
  admission floor required harm; the floor, ranking and reported mean now use
  the negated tau-hat for those outcomes (new `dina_subgroup()` argument
  `tau_sign`, default `1`), and survival, count and `adverse_outcome = TRUE`
  runs are unchanged.

* The internal MR wrapper `.fs_apply_mr()` -- which carries the **DINA** and
  **GRF** branches into `fs_mr_inference()` -- now forwards the **full**
  argument set.  It previously forwarded 15 of 25, dropping
  `return_reselection`, `field_R_out`, `field_R_in`, `field_uniform`,
  `field_M_cap`, `field_complement`, `field_decompose`,
  `field_scale_complement`, `ij_residual` and `field_recovery`.  Two field
  products (`field_decompose`, `field_recovery`) were therefore **unreachable**
  on those two engines, and three more were **inert while their defaults
  coincided** with the intended values -- so a run's meta recorded them as set
  while they controlled nothing, and flipping one to test the alternative did
  nothing.  `mr_inference_args` now controls the DINA and GRF branches exactly
  as it controls the consistency branch.  **Add-only: no default changes on any
  branch.**  Each forwarded default is read from `formals(fs_mr_inference)` at
  call time rather than restated in the wrapper, so a call that asks for nothing
  is byte-identical to the same call before the change (verified on both engines
  against captures from the pre-change install).  `ci_method` is unchanged and
  keeps the wrapper's own `"ij"` default.

* `fs_mr_inference()` gains `field_recovery` (default `FALSE`): add-only,
  default-inert **membership-agreement diagnostics** for the field's
  re-selection, attached as `field$recovery`.  Over the outer draws that
  produced a winner, each re-selected candidate's patient set is compared with
  the observed subgroup's and averaged, giving `sens_H` (the mean share of the
  identified patients the re-selections retain -- the primary quantity),
  `ppv_H`, `sens_Hc`, `ppv_Hc`/`npv_Hc`, the containment quantiles
  `q10`/`q50`/`q90`, `share_equal_1`, and the draw accounting.  The four
  classification metrics reuse the cross-validation membership cross-tab's
  definitions (`forestsearch_Kfold()`'s `sens_H`/`sens_Hc`/`ppv_H`/`ppv_Hc`)
  rather than parallel ones.  They exist because the previously reported
  exact-match frequency `p_hat` cannot distinguish re-selecting near-twins of
  the observed subgroup from re-selecting unrelated regions.  Computed from
  draws already made -- no new fits, no new randomness, no RNG consumption --
  and **no construction reads them**, so every bound, bias and SE is identical
  whether or not the block runs.  They answer a narrower question than the
  bootstrap or cross-validation recovery diagnostics: re-selection within the
  fixed kept family under perturbation, not re-discovery from scratch; the two
  are related but not interchangeable.  `summary.forestsearch()` reports
  `sens_H` beside `p_hat` when present (TASK_field_recovery_2026-09-09).

* The `ci_method` default is now `"field"` (was `"ij"`), in
  `fs_mr_inference()` and in the `forestsearch()` pass-through fallback, so a
  default call produces the recommended constructions: the field-calibrated
  one-sided lower bound on the harm subgroup, the studentized field-s
  one-sided upper bound on the complement, and the Bonferroni joint pair.
  **The IJ SE and its intervals are unchanged and still returned in every
  call** -- the IJ variance is computed before and independent of the field
  gate, so `"field"` adds the field block without removing anything, and the
  `debiased` element is identical to what `"ij"` returns.  Prior behaviour is
  `ci_method = "ij"`, which omits the field block.  The new default adds a
  per-fit Monte Carlo cost: `field_R_out` x `field_R_in` (defaults 1000 / 500)
  outer/inner draws, plus the complement's own block under the (also default)
  `field_complement = TRUE`.  Together with the Part D flips this makes the
  certified products reachable from a default call
  (TASK_cimethod_note_2026-09-09, Part D2).

* `fs_mr_inference()` defaults are now the recommended constructions (field
  complement bound, studentized field-s, re-selection diagnostics):
  `field_complement = TRUE`, `field_scale_complement = "selected"`,
  `return_reselection = TRUE`, mirrored in the `forestsearch()` pass-through
  fallbacks.  Prior behaviour is reachable by setting the arguments
  explicitly (`FALSE` / `"none"` / `FALSE`).  `field_decompose` stays `FALSE`.
  Part D left the `ci_method` default at `"ij"`, so the complement field block
  still ran only when `ci_method = "field"` was set explicitly; that default
  was flipped in the bullet above (TASK_cert20_2026-09-08, Part D;
  superseded on this point by TASK_cimethod_note_2026-09-09, Part D2).

* `fs_mr_inference()` gains `field_decompose` (default `FALSE`): add-only
  complement-field scale diagnostics attached as `field$complement$decomp_fields`
  when the complement field block runs; forwarded through
  `forestsearch(mr_inference_args = list(field_decompose = TRUE))`.  Defaults
  are byte-identical (PROPOSAL_complement_field_scale_2026-09-08_v2, Stage 1).
* `fs_mr_inference()` gains `field_scale_complement = c("none", "selected")`:
  under `"selected"` the complement field is studentized to the selected
  complement's influence-norm scale (the percentile-t root, variant R1) and
  `field$complement` gains the `_s` companions of its bounds with
  `field$joint_s` beside `field$joint`; add-beside, no existing field
  changes, the default reproduces prior output exactly
  (PROPOSAL_complement_field_scale_2026-09-08_v2 s5; TASK_field_studentize_e1).

* **Documentation: the threshold arguments.** The roles of
  `effect.threshold` / `hr.threshold` (`c1`, the screening threshold on a
  candidate subgroup's own effect), `consistency.threshold` /
  `hr.consistency` (`c2`) and `pconsistency.threshold` (`p*`) are now stated
  in every exported function that carries them, together with the defaults
  `c1` and `c2` resolve to under each estimand -- `hr.consistency`'s
  documented default of `1.0` is the survival value and is remapped on the
  GLM paths.  `consistency.threshold` is an **effect** threshold applied
  within each split half, **not** a proportion and not the consistency rate;
  the rate is `pconsistency.threshold`, which is never remapped.
  `fpr_calibration()`'s `c1` / `c2` are cross-referenced to the
  `forestsearch()` names in both directions.  Documentation only: no
  behaviour, default or code change.

## forestsearch 0.3.5 (development; not released on CRAN)

## Performance: unadjusted survival fits run on subset vectors

* The unadjusted arm of the per-candidate Cox fit no longer builds the
  per-candidate data frame: `coxph()` runs on three subset vectors
  (`Y`, `E`, `Treat`) resolved from the fitting environment, with the
  identical formula, solver, coefficient name, and 0.3.3 CI arithmetic.
  The frame assembly is relocated verbatim into the adjusted arm -- the
  only arm that still needs covariate columns.  A relocation plus
  substitution verified bit-identical (`identical()` on every retained
  result component across nine fixtures, plus a 20-replicate gbsg
  bootstrap payload, against the committed 0.3.4 baseline).

## forestsearch 0.3.4 (development; not released on CRAN)

## Performance: survival candidate medians computed only for effect-screen survivors

* The per-candidate median survival computation (`survfit()` +
  `summary.survfit()`) is relocated from inside the per-candidate Cox fit to
  the post-screen point of the same iteration, so it runs only for
  candidates that pass the effect screen -- the only candidates whose
  medians reach any output surface (`hr.subgroups` consumers and display
  tables for the selected subgroup).  A relocation of the existing
  computation, not a change to it: verified bit-identical (`identical()` on
  every retained result component across nine fixtures, plus a 20-replicate
  gbsg bootstrap payload).  ~11.7x fewer `survfit` calls on the gbsg
  application (1,410 fitted candidates vs 121 screen survivors).  Adjusted
  and unadjusted survival arms both take the relocated path; continuous and
  binary outcomes never computed medians and are untouched.

## forestsearch 0.3.3 (development; not released on CRAN)

## Performance: unadjusted fast paths for the per-candidate fit

* The per-candidate effect fit takes a fast path when the configuration is
  unadjusted, routed by a construction-time predicate and verified
  bit-identical (`identical()` on every retained result component across
  twelve fixtures, plus a 20-replicate bootstrap payload):
  * Continuous (MD): when there is no `adjust_covariates` and no propensity
    adjustment, the estimator closure is the fast variant
    (`lm.fit()` on the explicit two-column model matrix plus
    `summary.lm()`/`vcov.lm()`'s own extraction operations), and the search
    loop fits candidates from the outcome/treatment vectors directly -- no
    per-candidate data-frame slice.
  * Survival (Cox): the treatment-only branch of the per-candidate fit
    computes the consumed `conf.int` row directly from the `coxph()` fit,
    replicating `summary.coxph()`'s arithmetic (non-robust, 95%), instead of
    building the full summary table.
* Adjusted, weighted and propensity-adjusted paths are unchanged (preserved
  verbatim as the `else` branches); binary/count outcomes are out of scope.

## forestsearch 0.3.2 (development; not released on CRAN)

## Signed orientation in `fs_oc_family_enumerate()`

* Opposite-sign region effects are now supported in
  `fs_oc_family_enumerate()` via the signed orientation
  `s = sign(m_tau[Q])` (the harm direction): the candidate means are the
  signed mixture `beta_g = s * (m_tau[Qc] + (m_tau[Q] - m_tau[Qc]) * PQg)`,
  replacing the oriented-absolute reading, which was undefined when the two
  regions' effects differed in sign and `stop()`ped. Benefit-direction
  candidates carry oriented-negative `beta_g`.
* Same-sign families are bit-identical to the previous reading (gated
  `identical()` on every stored field against the MD40 corrected family).
* The family object gains an `orientation` element (alternative branch
  only): `s`, the signed `m_tau` pair, and the oriented `tauQc`, `bint`.
  **`tauQc` may be negative** for opposite-sign families.
* A planted effect with `m_tau[Q]` exactly zero is rejected with guidance
  (no harm direction to orient by; plant a nonzero Q effect or use the
  null path). The null branch, `se_g`, and all downstream consumers are
  unchanged.

* Documentation: corrected the `use_grf` `@param`, which described the
  variable-importance role belonging to `vi.grf.min`. `use_grf` governs GRF
  candidate-cut generation only; variable-importance ordering/screening of
  the cut columns is controlled solely by `vi.grf.min` (`NULL` default since
  0.3.0). No behavior change.

## forestsearch 0.3.1 (development; not released on CRAN)

## The consistency screen treats `plan = "sequential"` as sequential

`subgroup.consistency()` validated `parallel_args` against `c("plan", "workers")`
and, when both were present, took the future-backed batching path even for
`plan = "sequential"`: the candidates were split into `future_lapply()` calls of
`min(2 * workers, n)` -- two at a time on one worker -- and every call paid a
full round of globals resolution around an evaluation costing a few
milliseconds.  `forestsearch_bootstrap_dofuture()` forces exactly that
configuration (`list(plan = "sequential", workers = 1L)`) on the inner search
of every replicate.  Measured (`REPORT_bootstrap_profile_2026-08-30.md`):
29.4 s per call, 86% of it in `future_lapply` globals resolution, against 4.5 s
for the identical search under the plain loop.

1. **Dispatch.**  A sequential plan now runs the consistency screen as the plain
   in-process loop **regardless of `workers`**.  No future backend is set up
   and no warning is emitted.  A parallel plan with `workers = 1` is still
   parallel (the one-worker batched backend) -- pass `plan = "sequential"` to
   run the plain loop.  The batch-size logic and `parallel_args$batch_size`,
   the knob for per-batch overhead on genuinely parallel plans, are untouched.
2. **Measured effect.**  Continuous MD fixture (n = 500, 749 candidates
   screened): 31.1 s -> 5.5 s per `forestsearch()` call (5.6x).  Survival gbsg
   fixture (120 candidates screened): 8.9 s -> 7.5 s.  A 20-replicate
   `forestsearch_bootstrap_dofuture()` on the continuous fixture: 32.2 s ->
   4.95 s per replicate; B = 1000 projects to 82 min sequential
   against 490-537 min before.  A `multisession` plan with one worker is
   unchanged (30.8 s -> 30.7 s).
3. **Results under `consistency_method = "resample"` are unchanged** -- the
   closed form consumes no RNG, and six reference configurations (sequential
   with and without early stopping, survival, a two-worker multisession
   control and a one-worker multisession control) are `identical()` before and
   after on the selected subgroup, its membership, the candidate table, the
   consistency counters and the fitted effect.
4. **`consistency_method = "split"` under a sequential plan:** the plain loop
   seeds a single stream once (`set.seed(seed)`) where the batched backend
   assigned one L'Ecuyer-CMRG stream per candidate, so callers who passed
   `plan = "sequential"` *together with* `workers` in 0.3.0 will see each
   candidate's consistency rate move within its own Monte-Carlo noise at the
   chosen `fs.splits` (measured at 50 splits: sd of the per-candidate
   difference 0.028 against a binomial sd of 0.042; selected subgroup
   unchanged).  Callers who passed `plan = "sequential"` without `workers`
   were already on the plain loop and see no change.  Invariance to worker
   count is unchanged within any plan.
5. The driver-style `parallel_args = list(plan = "sequential")` (no `workers`)
   now reaches the plain loop through this rule rather than by tripping the
   `"parallel_args missing required elements. Using sequential."` warning,
   which no longer fires for that input.

## forestsearch 0.3.0 (development; not released on CRAN)

## Default change: `vi.grf.min = NULL` -- GRF variable-importance screening is off unless requested

The default of `forestsearch(vi.grf.min = )` changed from `-0.2` to `NULL`.
With `NULL`, Section 5 (the per-call causal forest and the importance
ordering of the candidate cut columns) is skipped entirely; screening runs
only when a numeric value is passed.

1. **Screening is off by default.**  A numeric `vi.grf.min` fits the forest,
   orders the cut columns by importance and retains those with
   `vi_ratio > vi.grf.min`; at values `<= 0` nothing is filtered (importances
   are non-negative and the block is guarded by `vi_max > 0`), so the
   pre-0.3.0 default was an ordering only.
2. **Reproducing results from versions before 0.3.0 requires passing
   `vi.grf.min = -0.2` explicitly.**  Callers that never set the argument now
   take the `NULL` path.
3. **Measured effect where `max_n_confounders` does not bind** (260 paired
   synthetic runs, `REPORT_vi_grf_smoke_2026-08-30.md`): the selected
   subgroup's subject membership is identical, the fitted effect and `Pcons`
   are identical, the candidate families are identical as sets of subject
   memberships, and `n_candidates_evaluated` / `early_stop_candidate` are
   identical -- the ordering never reached the scan, because the candidate
   table is re-sorted by `(-HR, K)` before the consistency loop.  Only the
   clause order of `sg.harm` may differ (e.g. `{age <= 47} {biomarker_hi}`
   instead of `{biomarker_hi} {age <= 47}`), and each call is ~0.13 s faster.
4. **`max_n_confounders` has no effect when `vi.grf.min` is `NULL`.**  The
   truncation lives inside the variable-importance block.  A finite cap below
   the number of cut columns now triggers one warning saying the cap is
   inert and how to restore it (pass a numeric `vi.grf.min`); the cap's
   placement is unchanged.

`R/run_simulation_analysis.R` (`default_fs_params_*()`) still sets
`vi.grf.min = -0.2` explicitly and is unaffected.  No driver, application
document or payload was touched; results recorded under earlier versions
stand as statements about those versions.

## forestsearch 0.2.7 (development; not released on CRAN)

## New: `fs_family_report()` -- what is, and is not, deterministic about the candidate family

`fs_family_report(x, data = NULL, outcome_type = NULL)` takes a named list of
`forestsearch()` arguments or a fitted `forestsearch` object and returns one
row per stage of candidate-family construction -- cut construction, the LASSO
/ DINA / GRF front ends, GRF variable-importance ordering and the
`max_n_confounders` cap, dummy expansion, combination enumeration, the
`minp` / `rmin` / `n.min` floors, the per-arm floors, the effect screen,
near-duplicate removal, the `max_subgroups_search` cap, the `m1.threshold`
filter, the consistency screen, early stopping, winner selection and the
(inert) time cap -- each classified as `deterministic`, `disabled`, `inert`,
`data-dependent` or `data-dependent (not disableable)`, with the governing
arguments, their resolved values and a one-line note.  The `print()` method
leads with a verdict and ends with the stages no argument switches off: cuts
at sample quantiles, the `rmin` and `n.min` floors on sample counts, and the
statistics-keyed near-duplicate removal.  With `data` supplied it also counts
the cut columns and combinations on that data.

It mirrors `forestsearch()`'s own argument resolution (`sg_focus` aliases; the
`maxeff` overrides; `stop_threshold` reset to `NULL` for every focus other than
`maxeffCons`; `max_n_confounders` inert when `vi.grf.min = NULL`; `d0.min` /
`d1.min` skipped for continuous and count outcomes) and a drift-guard test
holds it to the engine.  **It reports and changes nothing**: it fits no model,
runs no search, and is called by nothing else in the package.  No default,
driver or document changed.

## forestsearch 0.2.6 (development; not released on CRAN)

## `fs_oc_family_enumerate()`: a null branch

A DGM whose `df_super$flag_harm` has no member (Q empty;
`generate_glm_dgm(model = "null")`) is detected structurally -- the same
fact the null driver asserts -- and cross-checked against `dgm$model`.
Under the null every candidate has the same true effect, so the family's
fields become `beta_g` = the common effect (oriented as the alternative
orients), `se_g` from the whole-population effective variance
(`fs_dgm_scale(dgm, regions = list(S = ...))`, the S row -- obtained through
that function's public `regions` argument, without editing it) at `(n, Pg)`,
`PQg = 0`, `sens_g = NA` (0/0, undefined), `spec_g = 1 - Pg`, `PQ = 0`.
Enumeration, floors, overlaps and the covariance are unchanged; `NA`
sensitivity propagates only to `E[sens]` downstream.

## `fs_oc_grid()` / `fs_oc_invert()`: the exact order-statistic reduction

At fixed `(n, gate, c2)` neither the eligible set nor the maxeffCons winner
depends on `c1`: per draw, `T` = the winner's `Bhat` (`-Inf` if nothing is
eligible) and declaration at `c1` is exactly `T >= c1`.  The one-block sweep
now takes the argmax once per `c2` and reads every `c1` off the reduction,
and `fs_oc_invert(solve_for = "c1")` returns the `k`-th largest `T`,
`k = ceiling(target * draws)` -- an order statistic, no search; `target` may
be a vector, served by one draw set.  Results are unchanged: a grid point is
still `identical()` to `fs_oc_predict()`, and the twelve stored slow-path
inversions are reproduced to the order-statistic step.

## Alternative-cell results unchanged

`fs_oc_predict()` at a fixed seed is `identical()` to the 0.2.4 reference on
both gates, and the analytic document's `worked-predictions` chunk is still
reproduced bit-identically.

## forestsearch 0.2.5 (development; not released on CRAN)

## New: `fs_oc_grid()` and `fs_oc_invert()` -- threshold sweeps on one draw set, and the declaration-rate inversion

`fs_oc_grid(dgm, forestsearch_args, n, c1, c2, ...)` evaluates every
`fs_oc_predict()` quantity over the full crossing of `n`, `c1` and `c2`,
enumerating the family and drawing the candidates' fitted effects once per
`n` (and per gate) and sweeping the thresholds against those draws: both
gates are thresholds on draws that do not depend on `c1` or `c2`, so the
sweep costs arithmetic, not draws.  A `block` argument generates the draws in
row blocks and accumulates, making memory O(block x M); `block = Inf` is the
exact one-block path, `identical()` to `fs_oc_predict()` at one grid point.

`fs_oc_invert(dgm, forestsearch_args, n, target, solve_for = c("c1", "c2"))`
finds the threshold at which the family declaration rate equals `target`
by bracketing and bisection on fixed draws (a monotone step function),
reporting the achieved rate with its MC SE, and returning `NA` with the
ceiling and the binding threshold named when the target is unattainable.

## `fs_oc_predict()`: results unchanged

The draw, gate and selection/functional blocks of `fs_oc_predict()` were
moved into internal helpers (`.fs_oc_draw`, `.fs_oc_gate`,
`.fs_oc_functionals`, `.fs_oc_result`) so that `fs_oc_grid()` shares the
gate rather than copying it.  This is a move-only refactor: `fs_oc_predict()`'s
output at a fixed seed is `identical()` to 0.2.4's on both gates, and the
analytic document's `worked-predictions` chunk is still reproduced
bit-identically.

## forestsearch 0.2.4 (development; not released on CRAN)

## `fs_oc_predict()`: the `"resample"` gate is implemented

`consistency_method = "resample"` -- the package's production consistency
screen -- now runs in `fs_oc_predict()` instead of stopping.  The gate is the
closed-form inversion of `consistency_resample()`'s rate,
`2*pnorm((beta_hat - c2)/sigma_D) - 1 >= pconsistency`, i.e.
`Bhat >= c2 + qnorm((1 + pconsistency)/2) * sigma_D`, combined with the effect
screen `Bhat >= c1`.  `sigma_D` -- `sqrt(sum(dfbeta[i, treat]^2))`, the
sandwich SE of the subgroup treatment coefficient -- is identified with the
family's `se_g`: on one simulated MD40 trial it equals the direct dfbeta sum of
squares to 1e-10 and sits within 4% of the model SE and within 7% of the
population `se_g` with no prevalence trend
(`dev/glm-continuous-sims/sigma_d_diagnostic_2026-08-29.R`).  No constant
factor is carried.  New argument `pconsistency` (default
`forestsearch_args$pconsistency.threshold`, then `forestsearch()`'s 0.90).
The branch draws only `Bhat ~ N(beta_g, Sg)` -- one draw matrix, not two.
The `"split"` branch is untouched and remains bit-identical to the analytic
document's `worked-predictions` chunk.

## forestsearch 0.2.3 (development; not released on CRAN)

## New: `fs_oc_family_enumerate()` and `fs_oc_predict()` -- predicted operating characteristics over a population-enumerated family

`fs_oc_family_enumerate(dgm, forestsearch_args, n)` builds the candidate
family the search ranks over as a **population enumeration under the
`forestsearch()` cut specification**: the package's own cut machinery
(`get_FSdata()`, `conf.cont_jcuts`, `cut_type`, `collapse_cuts`, ...) is run
on the DGM's super-population frame `df_super`, so cuts land at population
quantiles; both directions of every cut are enumerated to `maxk` with the
search's own combination helpers; and only the structural floors are applied
(`minp` per factor, the `rmin` redundancy rule, and `Pg >= n.min / n`).
Every mean and standard error comes from `fs_dgm_scale(dgm)`.

`fs_oc_predict(dgm, forestsearch_args, n, c1, c2, family, consistency_method,
draws, seed)` is the `worked-predictions` assembly of the analytic prediction
document with its family construction removed: joint-normal half-sample draws
through `fs_sym_root()`, the declaration gate, the maxeffCons argmax, and the
selection-weighted functionals (detection rate, E|Hhat|, sensitivity,
specificity, PPV, NPV, E[beta(Hhat)], naive bias, mass below the floor), with
Monte-Carlo SEs for the proportions.  Given the document's M = 16 family, seed
and draw count it reproduces the document's figures bit-identically
(`dev/glm-continuous-sims/fidelity_fs_oc_predict_2026-08-28.R`).

`consistency_method` mirrors the package's own argument.  `"split"` is
implemented (the document's single-split gate).  `"resample"` stops with an
explanation: the package's closed-form screen uses `sigma_D`, the robust
sandwich SE of the subgroup fit (`R/consistency_resample.R`), which has no
expression in the population quantities the scale table carries; a stand-in
would be a method decision, not a port.

No existing function changes.

## forestsearch 0.2.2 (development; not released on CRAN)

## Builders own population noise

setup_gbsg_dgm() and generate_glm_dgm() gain k_random_noise/noise_seed:
population-level noise drawn once onto df_super at construction (scheme of
2a4787bc, previously driver-local); recorded as noise_scheme='population'.
Default k_random_noise = 0 is output-identical to 0.2.1.

## New: `fs_dgm_scale()` -- computed sampling scale of the difference-in-means estimator

Enumerates the finite-population sampling scale of the within-region
difference-in-means estimator directly from a DGM's potential outcomes:
`Var[betahat(g)] = V_eff(g) / (n P(g))`, with every component column returned
per region. Replaces the practice of anchoring the analytic documents on a
single *measured* standard error -- the measured anchor of 13.6786 at
n = 1000 implied an effective bracket of 16,119, below the fixture's residual
variance of 16,256, a structurally impossible state that a computed scale
cannot express. Exact for identity-scale effect measures (`MD`, `RD`, `IRD`);
ratio measures are rejected rather than silently approximated.

## New: `fs_scale_se()` -- estimator SD at a given sample size and region

Companion to `fs_dgm_scale()`: the standard deviation of the estimator at
sample size `n` on a named region, `jensen = TRUE` applying the
finite-count (random region size and arm split) corrections for the
unconditional standard deviation.

## New: `fs_mr_oc_summary()` -- operating characteristics from a payload

Summarizes estimation and identification operating characteristics directly
from a saved simulation payload -- only the `results`, `truth`, and `meta`
blocks are required, so no DGM argument is needed. The simulation drivers
attach the summary to every saved bundle as `$oc`, so a document's tables and
its payload cannot disagree.

## New: `fs_sym_root()` -- symmetric, continuous matrix square root

Promoted from the analytic prediction document's anchor chunk, where no test
could assert the property that motivated it. Both the symmetric root
`V D^{1/2} V'` and the asymmetric `V D^{1/2}` satisfy the covariance identity,
but the asymmetric root depends on the eigenvector basis, which is not a
continuous function of its input on the rank-deficient covariances that
complement-pair candidate families produce. The test suite asserts the
continuity difference directly: under a 1e-12 perturbation the symmetric root
moves < 1e-6 where the asymmetric root jumps by order one.

## forestsearch 0.2.1 (development; not released on CRAN)

## `sg_focus` alias resolution is now announced at run time

Five spellings resolve to one rule, and on the DINA and GRF paths two more
collapse onto the effect argmax. All of that was silent: a run requested as
`sg_focus = "eff"` ranked by a rule the caller never named, and nothing in the
transcript recorded the substitution.

`forestsearch()` now emits exactly one `message()` -- suppressed by `quiet`,
never `cat()` -- when the canonical rule differs from the spelling passed:

```
sg_focus 'eff' resolves to canonical rule 'hr' (aliases: eff, maxcons).
sg_focus 'maxeffCons' ranks as the effect argmax on the dina path (no Pcons is computed).
```

The alias relation now lives in one table (`.FS_SG_FOCUS_ALIASES`) that both
the normalizer and the announcement read, so the message cannot claim an alias
set the normalizer does not implement.

**No selection changes.** This is an announcement, not a rule change: the
selected subgroup and every estimate are bit-identical across the change, for
each of `subgroup_method` in `c("consistency", "dina")` crossed with `sg_focus`
in `c("eff", "hr", "maxcons", "maxeffCons")`.

## New: `fs_focus_tag()` -- one source of truth for stem tags

The `(subgroup_method, sg_focus)` -> stem-tag map had been pasted into 51
simulation drivers, and the copies were method-blind. On the DINA and GRF
paths they tagged output with the consistency name for a rule those engines
never ran: a `maxeffCons` DINA run wrote files stemmed `maxeffCons` when the
run had ranked by `order(-eff)`. Files named for a rule that did not run are
not recoverable after the fact.

`fs_focus_tag(subgroup_method, sg_focus)` is exported and tested against the
full 11 x 3 matrix, and every driver now calls it. Five drivers change the
stem they write (the DINA and GRF ones); the other 46 are unaffected.

## The DINA bootstrap `sg_focus` whitelist accepted only five of seven foci

`dina_subgroup_bootstrap()` held a re-inlined literal subset, so `maxeff` and
`maxeffCons` -- accepted by both `forestsearch()` and `dina_subgroup()` -- were
rejected here and the bootstrap stopped before its first iteration. It now
assigns from the shared `.FS_SG_FOCUS_CANONICAL`, and
`.assert_sg_focus_dispatch_complete()` covers this site so the literal cannot
come back. The rejection message enumerated five spellings for a whitelist of
seven; `dina_subgroup()`'s enumerated seven but omitted the accepted alias
`maxcons`. Both are corrected.

## The GRF frontier rule map is exhaustive

The `frontier_rule` switch had a default arm that silently ran `"effMaxSG"` for
any unrecognised focus, and a `tryCatch` above it swallowed normalization
errors. Both are gone: the map is exhaustive over `.FS_SG_FOCUS_CANONICAL` and
errors otherwise. Unreachable from `forestsearch()`, which normalizes and
whitelists first -- it exists so a focus added without a branch fails loudly
rather than selecting under a rule nobody asked for.

## Bug fix: MR harm confirmation wrote to stdout

The per-replicate `MR harm confirmation:` line used `cat()`, violating the
message()-only contract the MR messaging section states for itself, and
flooding parallel simulation renders (116 blocks at `n_sims = 20`).
`suppressMessages()` could not silence it because stdout is not the message
stream. It now routes through `.mr_msg(quiet, ...)` like every other MR report.

## Documentation

- `?forestsearch` gains a per-engine `sg_focus` table giving, for every
  spelling, the rule that actually runs -- and states the collapse sets
  explicitly: three spellings are one rule on `"consistency"`, five are one
  rule on `"dina"` and `"grf"`.
- The synonym collapse is mirrored into `dina_subgroup()` and
  `dina_subgroup_bootstrap()`, where it had lived only in a code comment.
- `m_diff` is documented as the candidacy floor on the link scale
  (`log(HR)`/`log(OR)`/`log(IRR)`, raw MD for gaussian), and as **derived**
  from the resolved `effect.threshold` under `subgroup_method = "dina"` -- a
  value supplied in `dina_args` is ignored there. Under
  `dina_select_statistic = "effect"` it still fixes the qualifying family that
  multiplier resampling de-biases over, while selection re-ranks on the refit
  effect.
- The `grf_selection` note warning that `maxeff`/`maxeffCons` "fall through
  without a condition being raised" is removed; they now map explicitly.

## forestsearch 0.2.0 (development; not released on CRAN)

## `consistency_method` now defaults to `"resample"`

The default is `"resample"`; `"split"` is retained as the fallback. Resample is
the intended method: it applies to the Cox path and to GLM outcomes whose
effect is a single model coefficient, and it is the cheaper of the two --
`"split"` costs `fs.splits * 2` refits per candidate evaluation.

**Behavior note:** callers that omit the argument now get resample-consistency
where they previously got split-consistency. Two consequences worth knowing:
results are not bit-identical across the change, and on a GLM outcome
multiplier-resampling inference (`mr_inference = TRUE`) now runs by default,
having previously been skipped because it requires `"resample"`. Pass
`consistency_method = "split"` explicitly to retain the old behavior.

## Breaking: `ps_method != "none"` is now an error for survival outcomes

Propensity adjustment is **not implemented** for `outcome_type = "survival"`.
The estimator-closure rebuild is gated on `if (is_glm)`, and no Cox estimation
path reads `sw`, `ps_hat` or `ips_covar` -- those names appear in no Cox file.
So the score was estimated, the columns were attached, nothing consumed them,
and the returned object still reported the method that was asked for.

Measured on GBSG: `ps_method = "grf"` at `ps_adjust_method` of `"iptw"`,
`"dr_gcomp"` and `"none"` gives results **bit-identical** to
`ps_method = "none"` -- same selected subgroup, same `treat.recommend`, same
`result_1sided` at tolerance 0.

This errors rather than warning, on the same grounds as the silently ignored
`dmin.grf` on the GRF effect path: an adjustment argument that cannot change
any result is a defect, not a documentation matter. There is no configuration
in which it does something, so a warning would invite the reader to think there
is a setting worth tuning.

**Note the default.** `ps_method` resolves to `"grf"` when `is.RCT = FALSE`, so
an observational **survival** analysis that never mentioned `ps_method` now
errors and must pass `ps_method = "none"` explicitly. The message says so. Such
a run was already getting no adjustment while reporting `ps_method = "grf"`, so
the break surfaces a result that was not what it claimed.

For confounding control on a survival outcome use `adjust_covariates`, which
the Cox path does consume.

## Bug fix: the bootstrap carried a user-supplied `ps_hat` into replicates unaligned

`forestsearch()` assigns a user-supplied `ps_hat` positionally
(`df.analysis$ps_hat <- ps_hat`) and derives the IPTW weights from it. A
bootstrap resample has the **same row count** but different subjects, so the
length check passed while every score belonged to a different person. The
replicate then de-biased using propensity scores attached to the wrong
subjects.

`forestsearch_bootstrap_dofuture()` now nulls `ps_hat` in its
`CATEGORY 2: VARIABLE RE-SELECTION` block, alongside `grf_res`, `grf_cuts`,
`dina_res` and `dina_cuts`. Each replicate estimates its own score under the
replicate's `ps_method`. This matches cross-validation, which nulls `ps_hat` at
both entry points on identical reasoning -- the asymmetry was the defect.

**This changes bootstrap results for anyone passing `ps_hat` explicitly**, and
only for them: with `ps_hat = NULL` (the default) the score was already
estimated per replicate and nothing changes. Measured on ACTG175, binary/OR,
`ps_adjust_method = "iptw"`, 40 replicates on a fixed seed: the selected
subgroup and its size are unchanged, the FB estimate moves 2.0919 -> 2.2434
(+7.2%), and the interval widens from [0.812, 5.389] to [0.625, 8.059].

Note the effect is confined to GLM outcomes. On the survival path the score is
attached and never read, so nothing consumed the misaligned values there.

## New: `mean_r` and `mean_r_c` on the MR object

`fs_mr_inference()` now returns two diagnostics beside `selection_rate`:

- `mean_r` -- the mean of the infinitesimal-jackknife residual over the draws
  the IJ actually uses (`ok_H`), not over all `draws`.
- `mean_r_c` -- the same for the complement, `NA_real_` when no complement was
  fit.

**No estimate changes.** This is exposure of a quantity the correction already
formed; no computation was touched.

They exist to make one invariant checkable: the residual mean is **zero by
construction whenever both bias terms share a denominator**, which is the
convention the package implements -- `selection_bias` and `fixed_bias` both
average over the draws that produced a winner, and the IJ runs on that same
set. A non-zero value means the two terms are being normalised differently
somewhere. That was exactly the defect corrected earlier in this cycle, whose
signature was a non-zero residual mean; until now nothing asserted its absence
directly.

## Breaking: `max_subgroups_search` now defaults to `Inf`

`forestsearch()` capped the candidate pool at the top **200** subgroups before
consistency evaluation. The default is now `Inf` -- no truncation.

The cap has no counterpart in the method. The candidate family is every
conjunction of at most two conditions, less empty intersections and
size/event failures, and selection ranges over that family. Truncation applied
a *preview* ordering and discarded everything below the cut, so the subgroup
that is optimal under the actual selection criterion could be dropped before
it was ever evaluated -- most acutely for the size foci `"maxSG"` / `"minSG"`,
where small high-effect subgroups sort last, and for the band foci
`"effMaxSG"` / `"effMinSG"`.

**This changes FS results wherever the cap previously bound**, i.e. wherever
the post-screen pool exceeded 200 and the truncation warning fired. Where the
pool was smaller than the cap, nothing changes. Runs that relied on the cap for
runtime should set a finite value explicitly, knowing it can change which
subgroup is selected:

```r
forestsearch(..., max_subgroups_search = 200)   # previous behaviour
```

`sg_focus = "maxeff"` already forced `Inf` and is unaffected.

## Identifier alignment: new defaults (BREAKING) and an MR configuration guard

Multiplier resampling (MR) linearizes the selection event, which is valid only
when two conditions hold: selection ranks on the same within-subgroup effect
the analysis reports, and the candidate family is fixed rather than estimated
from the outcome data.  Several package defaults did not satisfy them.  Five
defaults changed and a hard configuration check was added.

### Breaking: `use_lasso` and `use_grf` now default to `FALSE`

`forestsearch()` previously ran a prognostic Cox-LASSO prefilter and a GRF
importance screen before enumerating candidates.  Both are now off by default.

**This changes the candidate family on every `forestsearch()` call with
`subgroup_method = "consistency"`, not only under MR.**  A different family
can yield a different selected subgroup, so results from 0.1.0 and earlier
will not necessarily reproduce at the new defaults.

Reproducing Leon et al. (2024, *Statistics in Medicine*, doi:10.1002/sim.10163),
which used the prognostic lasso prefilter, now requires setting the flags
explicitly:

```r
forestsearch(..., use_lasso = TRUE, use_grf = TRUE)
```

That configuration remains fully supported.  It is simply no longer the
default, because a model-based front end makes the family data-dependent.

### Changed defaults: three `match.arg()` re-orderings

The permitted values of all three are unchanged; only which one is the default
moved.  Each previously defaulted to its engine's native statistic, which is
not the effect the analysis reports.

| argument | was | now |
|---|---|---|
| `dina_select_statistic` | `"dina"` (subgroup-mean tau-hat) | `"effect"` |
| `grf_selection` | `"tree"` (policy-tree objective) | `"frontier"` |
| `grf_select_statistic` | `"dr"` (mean doubly-robust score) | `"effect"` |

`use_dina` is unchanged (`FALSE`).

### New: MR configurations are validated, not silently accepted

A misaligned configuration combined with MR previously ran to completion and
returned a correction that did not de-bias the reported effect.  It now raises
an error naming the offending argument and the setting that resolves it, at
three entry points:

* `forestsearch()` under `mr_inference = TRUE`;
* `forestsearch_bootstrap_dofuture()` under `mr_in_replicates = TRUE`;
* `forestsearch_Kfold()` under `mr_in_replicates = TRUE`.

Rejected: `dina_select_statistic = "dina"`, `grf_select_statistic = "dr"`,
`grf_selection = "tree"`, and -- for `subgroup_method = "consistency"` only --
any of `use_lasso` / `use_grf` / `use_dina` left `TRUE`.  That last check is
**not** applied to `subgroup_method` of `"dina"` or `"grf"`: their families
regenerate under resampling by construction, which is a property of the
methods, not a configuration error.

Errors only -- no warnings, no silent coercion.  All six configurations remain
reachable; the three misaligned ones now require explicit opt-in and are
unavailable under MR.

### New: `$family_status` on the returned object

`forestsearch()` now records the candidate-family status as a character
scalar, so a sweep can tabulate it without parsing text: `"no-front-end"`,
`"conditional-removable"` (a front end is on, and turning it off would reach
`"no-front-end"`), or `"conditional-inherent"` (DINA/GRF).  The first level is
deliberately not called `"fixed"`: it records that no fitted model shapes the
family on the observed data, which is weaker than the manuscript's fixed
family, that additionally requires resample-invariant cut locations.  Shown by
`print()` and `summary()`.  Descriptive only; it never raises a condition.

### Unaffected: the full bootstrap

`forestsearch_bootstrap_dofuture()` remains valid for all six configurations.
It replays the whole pipeline per resample rather than linearizing it, so the
alignment conditions do not bear on its validity -- it is the reference
standard those conditions are stated against.  See its documentation.

## Terminology: "gate" and "tier" are gone

The post-selection-inference surface was named after a "gate" while mostly
doing something else, and the two de-biasing routes were called "Tier 1" and
"Tier 2", which said nothing about what either one does.  Both vocabularies
are retired.  **None of the renamed names appeared in CRAN 0.1.0** -- the
whole family was added after that release -- so these are renames of
unreleased names, not breaking changes for any CRAN user, and no deprecation
forwards are provided.  `?forestsearch` now carries a *Vocabulary* section
defining every term.

* **"Tier 1" is now the full bootstrap (FB)** --
  `forestsearch_bootstrap_dofuture()`.  It resamples subjects and re-runs the
  whole search per replicate.
* **"Tier 2" is now multiplier resampling (MR)** -- `fs_mr_inference()`.  It
  resamples nothing; it perturbs the influence contributions of a single set
  of fits.  MR is post-selection inference on a completed analysis and cannot
  change which subgroup is identified.

| removed | replacement |
|---|---|
| `debias_gate` (`forestsearch()` argument) | `mr_inference` |
| `debias_gate_args` | `mr_inference_args` |
| `out$debias_gate` | `out$mr_inference` |
| `out$harm_flag_debiased` | `out$mr_harm_confirmed` |
| `gate_estimates_table()` | `mr_estimates_table()` |
| `fs_debias_gate()` (internal) | `fs_mr_inference()` |
| `mr_inference_args$gate` | `mr_inference_args$confirm_rule` |
| `mr_inference_args$t_gate` | `mr_inference_args$t_confirm` |
| `c_gate` (`fs_fdr_report()` argument and `rep$fdr` column) | `c_confirm` |
| `out$mr_inference$gate` | `out$mr_inference$settings` |
| `out$mr_inference$gate$type` | `out$mr_inference$settings$confirm_rule` |
| `.fs_dg_*`, `.fs_apply_debias_gate()` (internal) | `.fs_mr_*`, `.fs_apply_mr()` |

The permitted values of `confirm_rule` are unchanged (`"point"`, `"ci"`);
only the argument name moved.  Three files were renamed to match their main
function: `R/fs_mr_inference.R`, `R/fs_mr_inference_methods.R`,
`R/mr_estimates_table.R`.

Unrelated uses of the word were deliberately left alone: the calibrators'
hard tolerance gates (a genuine pass/fail that `stop()`s), and the
`tier1`/`tier2` risk-difference fallback ladder in
`glm_effect_estimators.R`.

**`mr_harm_confirmed` is three-valued and `NA` is not a negative result.**
`TRUE` = MR ran and harm is confirmed; `FALSE` = MR ran and it is not;
`NA` = MR did not run or could not be computed -- including the default
`mr_inference = FALSE`, and the case where no subgroup was identified.
Because `isTRUE(NA)` is `FALSE`, code that reads this field with `isTRUE()`
alone reports "harm not confirmed" for an analysis that never ran.  Branch on
`is.na()` first.

## Known limitations of the MR step

Documented here because they are measurable and were not previously written
down.  None is changed by the rename; all are behaviour, and all remain open.

* **MR is silently skipped for GLM outcomes under the default
  `consistency_method = "split"`.**  The consistency branch runs MR only when
  `consistency_method == "resample"` with a GLM effect, or on the survival/Cox
  path.  For a binary/count/continuous outcome at the default
  `consistency_method`, setting `mr_inference = TRUE` does nothing, emits no
  warning, and leaves `mr_harm_confirmed = NA`.  Pass
  `consistency_method = "resample"` to get MR on a GLM outcome.
* **MR's re-selection family does not fully replay the identifier's
  admissibility.**  The per-arm event minima (`d0.min` / `d1.min`) and
  `max_subgroups_search` are not re-applied, so MR re-selects over a superset
  of the candidates the search actually chose among.  This inflates the
  selection-bias term, i.e. over-corrects the de-biased estimate.
* **MR is recomputed inside every bootstrap replicate and CV fold, and the
  result is discarded.**  `mr_inference` propagates through
  `args_call_all`, and neither `forestsearch_bootstrap_dofuture()` nor
  `forestsearch_Kfold()` overrides it.  This is wasted work only -- the
  reported result comes from the original analysis -- but at the default
  `draws = 2000` it is not cheap.

`tests/testthat/test-mr-inference.R` is new and covers this surface, which
previously had no regression tests at all; the first limitation above is
pinned there as present behaviour.

## Behaviour Changes

**These change which subgroup a search returns.**  0.1.0 is on CRAN; an
analysis re-run under 0.2.0 can select a different subgroup than it did
before.  Anyone reproducing published results should pin the version or
re-record the selection.

* **`sg_focus = "hr"` (and its aliases `"eff"`, `"maxcons"`) no longer stops
  early.**  Early stopping halts the candidate scan at the first candidate
  reaching `stop_threshold` and then applies the focus's selection key to that
  prefix.  That is sound only when the key is fully determined by the
  enumeration order, and `Pcons` never is -- it is unknown until a candidate is
  evaluated.  Any key containing a `Pcons` term can therefore be beaten by a
  candidate the scan never reached.  This affects `"hr"` (`Pcons` primary) and
  also `"maxSG"` / `"minSG"`, where `Pcons` breaks ties in subgroup size and a
  truncated scan can halt before a tied candidate with higher `Pcons` is seen.
  `stop_threshold` is now reset to `NULL` for all of these, with a warning when
  it was supplied explicitly.

  A visible consequence: the selected subgroup no longer depends on
  `parallel_args$batch_size`.  Before this change, that performance knob could
  alter a scientific result.

* **`stop_threshold` now defaults to `pconsistency.threshold`** instead of a
  fixed `0.95`.  Previously the two defaults were `0.95` and `0.90`, so a
  candidate with `Pcons` in `[0.90, 0.95)` qualified but did not halt the scan.
  The default is a promise, so it tracks a user-supplied floor rather than
  reverting to the old gap.  This is an efficiency change, not a correctness
  fix: where early stopping is sound at all, the prefix sort recovers the same
  winner for any `stop_threshold >= pconsistency.threshold`.

## New Features

* **New `sg_focus = "maxeffCons"`**: the effect maximiser among candidates
  clearing `pconsistency.threshold`.  It sits between the two existing
  effect-oriented rules -- `"maxeff"` maximises effect with no consistency
  condition at all, while `"hr"`/`"maxcons"` maximises consistency and uses
  effect only as a tiebreak.  `"maxeffCons"` is the only focus for which early
  stopping is sound, since its key `(-hr, K)` is exactly the enumeration order.
  Multiplier resampling replays this rule correctly (it maps to MR's
  `"maxeff"`, which already ranges only over consistency-qualifying
  candidates).

* **New `sg_focus = "maxcons"`**, an alias for `"hr"`.  `"hr"` selects the
  consistency argmax with effect only as a tiebreak, which is the opposite of
  what the name suggests to most readers; `"maxcons"` names the rule it
  actually implements.  `"hr"` and `"eff"` are unchanged and keep their
  meaning.

* `sort_subgroups()` and `sort_subgroups_preview()` now accept the alias
  vocabulary (`"eff"`, `"effMaxSG"`, `"effMinSG"`, `"maxcons"`) when called
  directly, as the other entry points already did.

* `parallel_args$batch_size` is now documented in `forestsearch()`.  It was
  accepted and validated but appeared in no roxygen block.


### GLM Extension

* Subgroup identification for binary, continuous, and count outcomes
  using odds ratios (OR), risk differences (RD), relative risks (RR),
  incidence rate ratios (IRR), incidence rate differences (IRD), and
  mean differences (MD).

* New `outcome_type` and `effect_measure` parameters in `forestsearch()`
  control the GLM pipeline.  Identification uses outcome-type-aware
  LASSO screening and GLM-based subgroup evaluation; estimation uses
  pluggable estimator closures via `make_effect_estimator()`.

* `generate_glm_dgm()`, `calibrate_glm_interaction()`, and
  `simulate_from_glm_dgm()` provide a simulation suite for binary,
  continuous, and count DGMs with calibrated treatment-by-subgroup
  interactions.

* `run_simulation_analysis()` dispatches automatically on DGM class:
  `"glm_dgm"` objects route to `simulate_from_glm_dgm()` and
  `grf.subg.harm.glm()` without user intervention.

### Consistency Resampling Approximation

* New `consistency_method` parameter in `forestsearch()` and
  `subgroup.consistency()`.  `"split"` (default) runs the original
  repeated 50/50 split-and-refit consistency calculation and is
  byte-identical to prior behavior; `"resample"` replaces it with a
  multiplier (influence-function) approximation that returns the rate
  from a single subgroup fit.  Each candidate's half-split treatment
  effects are represented as `beta_hat +/- D`, where `D` is a
  multiplier sum of the per-subject treatment `dfbeta` contributions,
  giving the closed form `2 * Phi((beta_hat - c) / sigma_D) - 1` with
  `sigma_D` the robust (sandwich) SE of the treatment coefficient and
  `c` the consistency threshold on the comparison scale.

* The `"resample"` path applies to the Cox (survival) outcome and to
  GLM outcomes whose effect is a single model coefficient (OR, RR, RD,
  MD, IRR).  Configurations it cannot represent that way --- IRD,
  propensity-adjusted effects, or a non-convergent fit on a given
  candidate --- fall back to literal splitting automatically, so a
  consistency rate is always produced.  When `use_twostage = TRUE`,
  `"resample"` bypasses the two-stage split screening entirely.

* The method propagates through the bootstrap automatically (each
  resample iteration inherits `consistency_method`), where the saving
  compounds because the inner `forestsearch()` runs sequentially.

* New exported `consistency_resample()` (Cox + GLM rate from a single
  fit) and `consistency_resample_compare()` (validates the
  approximation against literal splitting on the survival path).

### GRF Standalone Subgroup Identification

* `grf.subg.harm.glm()` provides standalone GRF-based subgroup
  identification for binary, continuous, and count outcomes using
  `grf::causal_forest()` and `policytree::policy_tree()`.

* `adverse_outcome` parameter controls the Y-flip before
  `causal_forest()`.  When `TRUE` (the default for binary and count
  outcomes), the outcome is flipped (`1 - Y` for binary, `-Y` for
  continuous) so that the most-negative-CATE subgroup corresponds to
  the harm subgroup.  This aligns the GRF sign convention with
  `forestsearch()`'s internal GRF pre-screening and with
  `grf::causal_survival_forest()`'s RMST convention.

* `tune_grf` parameter added to `forestsearch()`,
  `grf.subg.harm.glm()`, and `grf.subg.harm.survival()`.
  When `TRUE`, enables cross-validated hyperparameter tuning
  (`tune.parameters = "all"`) in `causal_forest()`, tuning
  `min.node.size`, `mtry`, `sample.fraction`, `alpha`, and
  `imbalance.penalty`.  Default `FALSE` preserves existing behavior.
  Most impactful for observational data with propensity score
  adjustment; defaults are near-optimal for RCTs (Dandl et al., 2024).

* `dmin.grf` is on the **risk difference** scale for binary outcomes
  (e.g., `dmin.grf = 0.10` requires at least a 10 percentage point
  absolute difference in event rates).  Calibration documents are
  provided for both synthetic 4-confounder settings and the ACTG175
  12-confounder setting.

### Selection Criteria

* `sg_focus = "hrMaxSG"` and `"hrMinSG"` redefined to use
  effect-size neighborhood selection.  Among candidates whose
  effect is within `effect_neighborhood` (default 10%) of the
  maximum effect, `"hrMaxSG"` picks the largest sample size and
  `"hrMinSG"` picks the smallest.  Setting `effect_neighborhood = 0`
  reduces these to a strict max-effect filter.  Applies to all
  outcome types; for GLM ratio measures (OR, IRR), the neighborhood
  test is on the natural scale.

* New `effect_neighborhood` parameter (default `0.10`) in
  `forestsearch()` and `subgroup.consistency()`.

* New `selection_rule` parameter (default `"neighborhood"`) controls
  the candidate-inclusion logic for `"hrMaxSG"` / `"hrMinSG"`.  One of:
    - `"neighborhood"`: within `effect_neighborhood` of the maximum
      effect (the legacy v0.2.0 behavior).
    - `"pareto"`: Pareto-non-dominated set on (effect, N), both
      maximized.
    - `"both"`: intersection of `"neighborhood"` and `"pareto"`.
  Must be `"neighborhood"` for single-criterion focus values
  (`"hr"`, `"maxSG"`, `"minSG"`).

* GLM-natural `sg_focus` vocabulary added as aliases for the existing
  Cox-flavored names:
    - `"eff"`      is an alias for `"hr"`
    - `"effMaxSG"` is an alias for `"hrMaxSG"`
    - `"effMinSG"` is an alias for `"hrMinSG"`
  Both vocabularies produce identical results.  The `"eff*"` forms
  read more naturally in GLM contexts (continuous MD, binary OR/RR/RD,
  count IRR) where there is no hazard ratio.  Old code using the
  `"hr*"` forms continues to work without changes.

* New canonical threshold names `effect.threshold` and
  `consistency.threshold` in `forestsearch()`.  The legacy
  `hr.threshold` and `hr.consistency` continue to work; the new
  names take precedence when both are provided.  Like the
  `sg_focus` aliases, the new names read more naturally in GLM
  contexts.

* Frontier-preserving preview sort.  The candidate triage that
  determines which subgroups enter consistency evaluation now
  guarantees that the full Pareto frontier appears in the top-K
  candidates, regardless of `selection_rule`.  Previously,
  restrictive rules such as `selection_rule = "both"` could crowd
  low-N frontier members out of the top-K by filling slots with
  higher-N dominated candidates, producing different post-consistency
  diagnostics across rules even on the same data.  Affects candidate
  *filtering* only; the post-consistency winner selection in
  `sort_subgroups()` is unchanged, so the selected subgroup is
  unaffected.

### Dual-View Candidate Diagnostic

* `forestsearch()` and `subgroup.consistency()` accept
  `show_candidate_summary = TRUE` to print two diagnostic tables
  during a run:
    - **Pre-consistency preview**: all candidates entering
      consistency evaluation, with Frontier and InBand flags
      computed from the candidate HR/N values.
    - **Post-consistency summary**: passing candidates with Pcons,
      Frontier, InBand, and Selected flags.
  Together the two views make the rule's filter visible end-to-end.
  Column headers use abbreviated forms (`P` for Pcons, `OF` for
  on-frontier, `IB` for in-band, `S` for selected) with a legend
  printed below each table.

### Pareto-Frontier Diagnostic Suite

* `pareto_frontier_table()` renders the frontier as a formatted
  `gt` table or returns it as a `data.table`.  Works uniformly for
  survival (HR) and GLM (OR, RR, IRR, RD, IRD, MD) outcomes;
  effect-column label and scale handling derived from
  `fs$effect_measure`.  Selected subgroup is marked with a ★ and
  optionally highlighted.  The `include_dominated = TRUE` option
  extends the view to all passing candidates (not just the
  frontier), with `on_frontier` and `in_band` flag columns.

* `plot_pareto_frontier()` renders the frontier as a `ggplot`
  scatter with optional ε-band shading and split-derived CI bars.
  Shared style conventions across all frontier plots (theme_bw,
  steelblue step polyline, `#D55E00` highlight for the selected
  point).

* `plot_pareto_combined()` composes multiple `forestsearch` fits
  onto a single Pareto plot when their consistency-passing sets
  are identical.  Each selected subgroup is annotated with one
  or more `S1: <combo_label>` markers naming the configuration(s)
  that picked it; multiple combos picking the same subgroup yield
  stacked multi-line labels.  Returns `NULL` with a clear warning
  when the equality precondition fails (same row count, same
  subgroup-definition set, hr/N/E/K agreement within tolerance);
  the side-by-side panel composition is the natural alternative.

* `compute_frontier_cis()` computes three CIs per frontier member:
  a naive CI (model-based or robust SE), a half-jackknife split CI
  (Shao 1996), and an FSBC-mimic CI.  Pluggable into the table and
  plot functions via `ci_table = ...`.

* `explain_pareto_selection()` produces a verbal account of why the
  selected subgroup wins on the configured lexicographic criterion
  relative to other non-dominated candidates.  Format and verbosity
  configurable; supports markdown rendering via `results = "asis"`
  in Quarto chunks.

* `frontier_member_flags()` returns a per-subject membership matrix
  indicating which frontier subgroups each subject belongs to;
  useful for downstream within-frontier comparisons.

* `compare_selection_rules()` wrapper runs `forestsearch()` across
  multiple `(sg_focus, selection_rule)` tuple combinations with all
  other parameters held fixed.  Captures each run's stdout, builds
  per-combo Pareto plots, composes them via `patchwork`
  side-by-side, and auto-builds a `plot_pareto_combined()` view
  when the passing sets match.  Returns a structured
  `forestsearch_comparison` object with per-combo fs, ci_tab, plot,
  console, and diagnostics slots.

### Other

* Propensity score adjustment for observational studies via stabilized
  IPTW with GRF, LASSO, or logistic regression PS estimation.

* Forest plot (`plot_subgroup_results_forestplot()`) extended to
  support GLM effect measures with auto-detected display defaults,
  clinically meaningful axis limits, and extreme CI capping.

* `plot_sg_glm_outcomes()`: bar charts (binary) and
  point-and-errorbar (continuous) outcome visualization.

* `glm_effect_profile()`: delta-method treatment effect profiles
  across continuous biomarkers with natural cubic spline interactions.

## Breaking Changes (vs. v0.1.0)

* **`sg_focus = "hrMaxSG"` / `"hrMinSG"` selection rule changed.**
  Previously these used lexicographic sorts (size primary, effect
  secondary); now they use effect-size neighborhood selection (see
  *Selection Criteria* above).  Results from these focus values will
  differ from v0.1.0 runs on the same data.  Single-criterion focus
  values (`"hr"`, `"maxSG"`, `"minSG"`) are unchanged.

* **`showten_subgroups` argument removed.**  Renamed to
  `show_candidate_summary`.  The old name no longer works; callers
  passing `showten_subgroups = TRUE` will see an
  `unused argument` error.  The replacement triggers two
  diagnostic views (preview + summary) rather than the legacy fixed
  "top 10" display; see *Dual-View Candidate Diagnostic* above.

* **`pareto_frontier_table()` `digits` default lowered.**  Master
  `digits` argument now defaults to `2L` (was `3L`).  Per-column
  overrides (`digits_effect`, `digits_pcons`, `digits_ci`) inherit
  from the master when not explicitly set.  Rendered tables are more
  compact; pass `digits = 3L` to recover the prior precision.

## Bug Fixes

* **Fixed an incoherence between the two multiplier-resampling bias terms'
  denominators.** `selection_bias` averaged over the draws that produced a
  re-selected winner while `fixed_bias` averaged over all `B` draws. The
  infinitesimal-jackknife residual combines the two, so it mixed differently
  normalised quantities and its mean over the draws the IJ used was not zero --
  precisely the condition the finite-`B` correction needs in order to be the
  centered quantity it is meant to be. Both terms, and the IJ, are now taken
  over the same draws: those on which a selection occurred.

  This conditions the correction on identification, which is deliberate. The
  reported analysis exists only because a subgroup was identified, the selected
  subgroup's effect is already a conditional estimand, and a bootstrap replicate
  that identifies nothing contributes nothing to the full bootstrap that MR
  approximates.

  **Wherever `selection_rate < 1` the de-biased point estimate shifts by exactly
  the old residual mean** -- the incoherence was surfacing directly as a bias of
  that size. Intervals are essentially unchanged, and **`selection_rate = 1`
  analyses are unaffected to every digit**. `$mr_inference$selection_rate` is
  the field that tells you which case you are in. Measured on ACTG175 at
  `selection_rate = 0.865`:

  | quantity | before | after |
  |---|---:|---:|
  | `selection_bias` | 0.6623937 | 0.6623937 (identical) |
  | `fixed_bias` | 0.0053309 | 0.0908747 |
  | de-biased OR | 1.52988 | 1.40445 |
  | `tilde_V` | 0.2965604 | 0.2965604 (identical) |
  | `se_ij` | 0.4190578 | 0.4217818 |

  Note that `fixed_bias` is now a **conditional** expectation and is no longer
  mean-zero: the draws that clear the admission floor are exactly those where
  the perturbation pushed effects up, so the conditioning event correlates with
  the quantity being averaged. It shrinks monotonically to zero as
  `selection_rate` approaches 1 -- the behaviour of a selection effect rather
  than an artefact. Both the shipped DINA configuration and GRF on ACTG175 ran
  at `selection_rate = 1`.

* Fixed cross-validation reusing a cached `grf_res` / `dina_res` inside every
  training fold.  `forestsearch_Kfold()` and `forestsearch_tenfold()` rebuild
  each fold's call from `args_call_all`, which carries all of `forestsearch()`'s
  formals -- including a user-supplied cached fit.  The bootstrap has always
  nulled `grf_res`, `grf_cuts`, `dina_res` and `dina_cuts` for exactly this
  reason, and CV nulled `ps_hat` on the same reasoning while leaving the other
  four.  A fold therefore ranked over a candidate family estimated on **all**
  the data, including the fold being held out, so `sensCV`, `ppvCV` and
  exact-match were optimistic by construction.  All four are now nulled at both
  entry points.

  **This changes CV metrics for anyone passing a cached fit**, and the effect is
  not small.  Measured on ACTG175 (binary, OR, `Kfolds = 3`), cached versus
  uncached:

  | method | metric | cached, before | uncached | cached, after |
  |---|---|---:|---:|---:|
  | GRF | `ppv_H` | 0.6486 | 0.4896 | 0.4896 |
  | GRF | exact match | 0.3333 | 0.0000 | 0.0000 |
  | DINA | `sens_H` | 0.4688 | 0.2188 | 0.2188 |
  | DINA | `ppv_H` | 0.3000 | 0.1538 | 0.1538 |

  After the fix the cached and uncached runs agree exactly on every metric.
  Runs that passed `grf_res = NULL` / `dina_res = NULL` (the default) are
  unaffected -- there was nothing to reuse.

* Fixed the influence path erroring on every `lm()` fit, which disabled
  multiplier resampling (MR) for continuous outcomes entirely.  The internal
  one-step `dfbeta` read `fit$weights` unconditionally; an unweighted `lm()`
  has `NULL` there rather than a vector of ones, so the computation was
  non-conformable and failed.  `effect_measure = "MD"` is the only measure
  fitted with `lm()`, and the failure was silent at every level:
  `consistency_resample()` returned `NA` rates, the consistency engine fell
  back to literal splitting while `consistency_method = "resample"` was in
  force, and `fs_mr_inference()` could fit no candidate at all -- so
  `mr_inference = TRUE` on a continuous outcome produced no correction and
  reported nothing.  **Any MD analysis run with `mr_inference = TRUE` between
  the introduction of the one-step `dfbeta` and this release has no MR result
  and should be re-run.**  Binary, count and survival outcomes are unaffected:
  all fit with `glm()` or `coxph()`, where the weights are always present.

* Fixed `grf.subg.harm.glm()` ignoring `adverse_outcome`: the Y-flip
  before `causal_forest()` was accepted but never applied, causing the
  policy tree to identify the complement of the true subgroup in
  benefit-search scenarios.

* Fixed `run_simulation_analysis()` hardcoding `simulate_from_dgm()`:
  added class-based dispatch so `"glm_dgm"` objects route to
  `simulate_from_glm_dgm()`.

* Fixed `run_simulation_analysis()` GRF dispatch: binary/continuous
  outcomes now route to `grf.subg.harm.glm()` instead of
  `grf.subg.harm.survival()`.  Also removes `event.name` from
  `grf_args` (not accepted by `grf.subg.harm.glm()`) and propagates
  `adverse_outcome` with the same default as `forestsearch()`.

* Fixed `valid_pnames` whitelist in `run_simulation_analysis()`:
  added `outcome_type`, `effect_measure`, `offset.name`, and
  `tune_grf` so GLM parameters pass through to `forestsearch()`.

* Removed stale duplicate `grf.subg.harm.glm()` from `grf_main.R`
  that conflicted with the canonical implementation in
  `grf_subg_harm_glm.R`.

* Fixed `merge(by = "id")` ignoring `id.name` parameter in
  `forestsearch()` output.

* Fixed `geom_errorbarh(height = 0, ...)` deprecation warning under
  `ggplot2 >= 3.5` in `plot_pareto_frontier()` and
  `plot_pareto_combined()`.  Replaced with `geom_linerange()` which
  produces an identical visual without the warning.

* Fixed non-ASCII characters in `plot_pareto_combined.R` for R CMD
  check compliance: em-dashes in user-visible string literals
  replaced with ASCII hyphens; ε and ≥ glyphs replaced with R
  unicode escapes (`\u03b5`, `\u2265`).  Plot output is unchanged.

* Fixed bare `<<-` assignment leakage inside nested
  `capture.output()` / `tryCatch()` in `compare_selection_rules()`.
  Replaced with an explicit holder environment; the prior pattern
  could silently write to the global environment under certain
  R session configurations.

## Internal

* `fit_causal_forest()` and `fit_causal_forest_glm()` accept
  `tune_grf` and pass `tune.parameters` to `grf::causal_forest()`
  and `grf::causal_survival_forest()`.

* Bootstrap bias correction extended for GLM outcomes with
  log-scale / identity-scale handling.

* Cross-validation extended to handle GLM parameters and propensity
  score re-estimation.

* New internal helper `.normalize_sg_focus()` in
  `forestsearch_helpers.R` translates the `"eff*"` vocabulary to
  the canonical `"hr*"` form at entry to `forestsearch()` and
  `subgroup.consistency()`.  Downstream code is keyed on the
  canonical form; the aliases are purely user-facing.

* New internal helper `extract_candidate_diagnostics()` in
  `compare_selection_rules.R` slices the `CANDIDATE EVALUATION
  PREVIEW` and `CANDIDATE EVALUATION SUMMARY` blocks out of
  captured forestsearch stdout for focused presentation.

* Pareto frontier on (effect, N) -- both maximized -- attached to
  consistency results as `out_sg$pareto_frontier`.  Post-hoc
  diagnostic listing non-dominated alternatives to the selected
  subgroup; not used for selection.

* `Pcons` excluded from the value-tolerance check in
  `plot_pareto_combined()`'s equality precondition.  Pcons can
  legitimately drift across rules because the preview-sort change
  alters queue order and therefore the random-split state each
  candidate consumes; equality is keyed on subgroup definitions
  and hr/N/E/K only.

* GLM consistency resampling is wired through an explicit
  `glm_resample_spec` threaded from `forestsearch()` to the
  consistency evaluators (rather than introspecting the estimator
  closure).  `consistency_resample()` gained a `comparison_threshold`
  argument so the engine passes the already-resolved comparison-scale
  threshold without a second log-transform.  New internal helpers
  `.glm_resample_supported()` and `.consistency_via_splits()` gate the
  supported measures and provide the shared split fallback.

## References

* Dandl, S., Haslinger, C., Hothorn, T., Seibold, H., Sverdrup, E.,
  Wager, S., & Zeileis, A. (2024). What makes forest-based
  heterogeneous treatment effect estimators work? *Annals of Applied
  Statistics*, 18(1), 506-528.

* Rehill, P. (2025). How do applied researchers use the causal forest?
  A methodological review. *International Statistical Review*.

* Shao, J. (1996). Bootstrap model selection. *Journal of the
  American Statistical Association*, 91(434), 655-665.

# forestsearch 0.1.0

* Initial CRAN release.
