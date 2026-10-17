# Multiplier resampling (MR) for a selected forestsearch subgroup

Computes a multiplier-bootstrap approximation of the bootstrap
bias-corrected treatment effect for the selected subgroup and flags
whether it is still consistent with harm. Reuses the resample
consistency engine's per-candidate treatment `dfbeta`; see the file
header for the closed form.

## Usage

``` r
fs_mr_inference(
  df,
  candidates,
  spec,
  selected_members,
  admission,
  t_confirm = NULL,
  confirm_rule = c("point", "ci"),
  reselection = c("maxcons", "maxeff", "maxSG", "minSG", "effMaxSG", "effMinSG"),
  effect_neighborhood = 0.1,
  selection_rule = c("neighborhood", "pareto", "both"),
  draws = 2000L,
  multiplier = c("poisson", "gaussian", "rademacher"),
  include_complement = FALSE,
  ci_method = c("field", "ij", "wald"),
  seed = NULL,
  return_reselection = TRUE,
  field_R_out = 1000L,
  field_R_in = 500L,
  field_uniform = FALSE,
  field_M_cap = NULL,
  field_complement = TRUE,
  field_decompose = FALSE,
  field_scale_complement = c("selected", "none"),
  ij_residual = c("two_term", "winner", "winner_floor"),
  field_recovery = FALSE,
  keep_declaration_field = FALSE,
  keep_field_matrix = FALSE,
  declaration_c0 = NULL,
  pconsistency.digits = NULL
)
```

## Arguments

- df:

  Analysis data frame (the standardized `df.fs` inside
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md));
  must contain the columns named in `spec`.

- candidates:

  Named list of integer row-index vectors, one per screened candidate
  subgroup (the family the selection rule chose among).

- spec:

  List with `outcome_type`, `effect_measure`, `treat.name`,
  `outcome.name`, `event.name`, `offset.name`, `adjust_covariates`,
  `adverse_outcome`.

- selected_members:

  Integer row indices of the observed selected subgroup
  (`which(grp.consistency$sg.harm.id == 1)`).

- admission:

  The resolved admission set, as returned by `.fs_resolve_admission()`:
  a list with `effect_floor` (numeric on the **comparison scale** – log
  for ratio measures – or `NULL`) and `consistency` (a list with
  `c_cons` and `p_star`, or `NULL`).

  This replaces the former `c_screen` / `c_consistency` / `p_star`
  arguments, from which MR used to rebuild
  `t_g <- pmax(c_screen, c_consistency + z * sigma_D)` at this site.
  That reconstruction was the defect: the identifier resolved its floors
  in one place and MR re-derived them in another, so the two could
  disagree about which candidates were ever in contention – and under
  `sg_focus = "maxeff"` they did. MR is a linearization of the selection
  *map*, which is a ranking and a domain; passing the resolved domain
  makes the two equal by construction.

  `NULL` means the floor does not apply, and is handled by a separate
  code path rather than by sentinel arithmetic. With neither floor
  present there is no admission filter at all.

- t_confirm:

  Harm-confirmation threshold on the **effect scale** (HR/OR/ RR/IRR, or
  RD/MD for differences – not the working log scale). `NULL` uses the
  near-null default: `1` for ratio measures, `0` for differences. Set it
  near the null rather than at the screening threshold, because the
  de-biasing correction over-shrinks true effects.

- confirm_rule:

  Which harm-confirmation rule to apply to the de-biased estimate:
  `"point"` requires the de-biased point estimate to clear `t_confirm`;
  `"ci"` requires the one-sided 95% selection-adjusted lower bound to
  clear it, and is therefore the stricter of the two.

- reselection:

  Bootstrap re-selection rule for the bias term; default `"maxcons"`
  (validated). `"effMaxSG"` etc. are available but approximate.

- effect_neighborhood:

  Band for the `eff*SG` re-selection rules.

- draws, multiplier, seed:

  Multiplier-bootstrap controls. `multiplier` defaults to `"poisson"`
  (mimics the nonparametric bootstrap).

- include_complement:

  Logical. When `TRUE`, also de-bias the complement subgroup (everyone
  not in the selected subgroup). The complement is not selected
  independently – its bias is induced by selecting the subgroup – so on
  each multiplier draw the complement of the re-selected winner is
  tracked and de-biased with the complement's own influence. Complements
  are fit only for candidates that win across draws (plus the selected
  one), so the extra cost is small. Default `FALSE`.

- ci_method:

  `"field"` (default, the recommendation) computes everything the `"ij"`
  path computes – the `debiased` element, its IJ SE and its IJ intervals
  are identical, and are returned in every call regardless of this
  argument – and **additionally** runs the field-calibrated interval,
  whose one-sided bounds are the primary reported products (the lower
  bound on the harm subgroup, the studentized complement's upper bound,
  and the Bonferroni joint pair). Nothing is removed by choosing it; it
  adds a per-fit Monte Carlo cost (`field_R_out` x `field_R_in`,
  defaults 1000 / 500, plus the complement's own block under
  `field_complement = TRUE`). `"ij"` restores the previous default: the
  **de-biased** CI from the infinitesimal-jackknife variance (Leon et
  al. 2024, Eq. VInfJ_bc), computed from the same multiplier draws – the
  leading-order analogue of the FB interval – with the field block
  omitted. `"wald"` uses the subgroup robust SE (`sigma_D`) for the
  de-biased CI and also omits the field block. The naive CI always uses
  the robust SE. The field block (León and Anderson, 2026) is returned
  as a `field` element: Gaussian-multiplier perturbations of the shrunk
  candidate-effect field (`w = beta_hat` with the winner's entry
  replaced by the de-biased estimate) are pushed through the configured
  re-selection map, with an inner Monte Carlo estimating the selection
  drift, and the resulting `Lambda*` distribution is inverted around the
  de-biased estimate (basic-bootstrap form). Drawn under
  `seed + 900000L` when `seed` is given, after – never inside – the main
  multiplier stream, so `"ij"` and `"wald"` output is unaffected
  byte-for-byte.

- return_reselection:

  Logical. When `TRUE`, the return list gains a `reselection` element
  exposing the per-draw re-selection already computed for the bias term:
  `winner` (integer vector of length `draws`, the index of the
  re-selected candidate in the kept family on each draw, `NA` where no
  draw winner exists) and `p_hat` (named numeric over the kept family,
  `tabulate(winner) / draws` – so `sum(p_hat) == selection_rate`, and
  the frequencies sum to 1 exactly when every draw produced a winner).
  Default `TRUE`: the re-selection diagnostics are part of the
  recommended output (\\\hat p\\ is reported by print() and summary() as
  a descriptive diagnostic). Nothing in the arithmetic depends on this
  switch, so `FALSE` reproduces the previous, smaller return object
  exactly.

- field_R_out, field_R_in:

  Outer and inner Monte Carlo sizes for `ci_method = "field"` (defaults
  1000 / 500); ignored otherwise. The inner draws are shared across
  outer draws.

- field_uniform:

  Logical (default `FALSE`); only consulted under `ci_method = "field"`.
  When `TRUE`, after the field block completes the gate additionally
  runs the uniform (kappa) calibration
  ([`fs_mr_field_uniform()`](https://larry-leon.github.io/forestsearch/reference/fs_mr_field_uniform.md))
  and attaches its result as `field$uniform`: the smallest widening
  factor `kappa` in `[1, 2]` such that the widened two-sided interval
  attains 95% coverage uniformly over the **winner-profile protection
  family** – hypothetical true fields equal to the shrunk field `w` with
  the winner's entry set to the runner-up level plus `delta` winner-SEs,
  `delta` on a 0-4 grid – computed from the trial's own influence
  structure (`Sigma_hat` via `B_eff`), restricted to the **mass-carrying
  candidate set** (smallest `M <= 12` candidates holding \>= 99% of the
  re-selection mass, winner always included). Experimental option.
  Intended properties: the one-sided field bound is uniformly valid; the
  two-sided interval widened by the reported `kappa` is uniformly valid
  over the winner-profile family; the plain quantile interval is
  approximate. The sweep draws under `seed + 910000L` after – never
  inside – the field's own stream, so every other output, the field
  block included, is byte-identical whether or not it runs.

- field_M_cap:

  Optional override of the uniform sweep's mass-carrying cap
  ([`fs_mr_field_uniform()`](https://larry-leon.github.io/forestsearch/reference/fs_mr_field_uniform.md)'s
  `M_cap`); `NULL` (default) uses that function's default. Pass a value
  at least the family size to lift the cap to the full family (the GBSG
  frozen-intervals illustration does). Consulted only when
  `field_uniform = TRUE`; add-only, default-inert.

- field_complement:

  Logical (default `TRUE`, the recommended construction); consulted only
  under `ci_method = "field"` with `include_complement = TRUE`. When
  `TRUE`, after the harm field block completes, the same field procedure
  is evaluated for the complement subgroup and attached as
  `field$complement`. **Orientation:** selection acts on the
  harm-subgroup effects; the complement is its by-product, and its naive
  estimate is pushed toward *benefit*. For a benefit claim ("the effect
  in the complement is at most U") the exposed limit is therefore the
  **upper** bound on the complement's working-scale effect, and the
  primary product here is the one-sided 95% upper bound `upper_1s`. The
  guarantee is the mirror of the harm side's: an upper bound is
  vulnerable to under-correction of the beneficial optimism (U too low
  over-claims) and immune to over-correction. Selection is the harm
  field's own – every outer and inner draw's winner is re-read, never
  re-selected – and the complement's perturbations use the harm field's
  own multipliers projected through the complement's influence, so its
  noise carries its correct correlation with every candidate's noise;
  complements of draw winners the multiplier stage did not fit are fit
  lazily. The block draws nothing, so every other output – the harm
  `field` block and its `uniform` sub-block included – is byte-identical
  whether or not it runs, so `FALSE` reproduces the pre-complement
  output exactly.

- field_decompose:

  Logical (default `FALSE`); consulted only when the complement field
  block runs. When `TRUE`, `field$complement` gains a `decomp_fields`
  list of add-only scale diagnostics (the selected complement's
  influence-norm scale, the draw-weighted mean and CV of the outer
  winners' scales, their ratio `scale_ratio_c`, and the variance /
  covariance of the two pieces of `Lambda*c`); it reads existing objects
  only, so every pre-existing output is unchanged whether or not it
  runs.

- field_scale_complement:

  `"selected"` (default, the recommended construction: the studentized
  complement field, "field-s") or `"none"`; consulted only when the
  complement field block runs. Under `"selected"` the complement's field
  readings are studentized to the selected complement's own scale before
  differencing – each outer reading `zeta^c_{r, G_r}` is multiplied by
  `s_sel / s_{G_r}` and each inner reading by
  `s_sel / s_{G(v_r + zeta'_j)}`, with `s_g = sqrt(sum(Bc[, g]^2))` the
  candidate's complement influence-norm scale (the percentile-t root) –
  and `field$complement` gains the `_s` companions of its bounds
  (`est2_s`, `upper_1s_s`, `lower_1s_s`, `lower_2s_s`, `upper_2s_s`,
  `lower_se_s`, `upper_se_s`, `se_field_s`, `lambda_mean_s`) with
  `field$joint_s` beside `field$joint`. Add-beside: no existing field
  changes under either setting; the default reproduces prior output
  exactly.

- ij_residual:

  Which IJ residual populates the **reported** de-biased SE and interval
  (`debiased$se_ij`/`lower`/`upper`/`lower_1s`, and the complement's).
  `"two_term"` (default) is the paper's residual \\r_b = (bias\_{sel} +
  bias\_{fix}) - D\_{H^\*\_b}(b) - D\_{\hat H}(b)\\, whose limit
  variance is \\(p + e\_{\hat H})' \Sigma (p + e\_{\hat H})\\ –
  \\4\sigma^2\\ when one candidate dominates, a doubled SE exactly where
  there is least to correct (the complement is effectively always in
  that regime). `"winner"` uses the winner-only (total-derivative)
  residual \\r^w_b = bias\_{sel} - D\_{H^\*\_b}(b)\\, exact when one
  candidate dominates but under-covering at ties on the harm side (the
  delta method does not hold in the local regime), likewise retained for
  comparison only; `"winner_floor"` floors it at the naive robust
  variance, \\\max(V^w, \hat\sigma^2_D)\\, retained for comparison only
  and not recommended. All three SEs are computed from the same draws
  and returned side by side in every call (`se_ij_two_term`,
  `se_ij_winner`, `se_ij_winner_floor` with their bounds `*_w` /
  `*_wf`); the argument only selects the reported one, and the point
  estimates never change. Add-only: the default reproduces prior output
  exactly (the new elements aside).

- field_recovery:

  Logical (default `FALSE`); consulted only under `ci_method = "field"`.
  When `TRUE`, `field` gains a `recovery` list of **descriptive
  membership-agreement diagnostics**: over the outer draws that produced
  a winner, each re-selected candidate's patient set `kept[[G_r]]` is
  compared with the observed pick's `kept[[sel]]`, and the four
  classification metrics of the CV membership cross-tab
  ([`forestsearch_Kfold()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_Kfold.md)'s
  `sens_H` / `sens_Hc` / `ppv_H` / `ppv_Hc`) are averaged over draws,
  with `npv_Hc` reported as the standing-convention name of the same
  quantity as `ppv_Hc`. These are computed **from draws already made**:
  no new fits, no new randomness, no RNG consumption, and **no
  construction reads them** – every bound, bias and SE is byte-identical
  whether or not the block runs. They answer a narrower question than
  the full bootstrap (FB) or cross-validation (CV) recovery diagnostics:
  re-selection *within the fixed kept family* under multiplier
  perturbation, not re-discovery of a subgroup from scratch on resampled
  data. The two are related but not interchangeable, and a report that
  quotes one should say which.

- keep_declaration_field:

  Logical (default `FALSE`). When `TRUE`, the return gains a
  `declaration_field` element: the standardized perturbation field of
  the main multiplier stream, `Zstar[b, g] = D_g(b) / sigma_D(g)`,
  summarized by its family maxima `Mstar` (length `draws`), with the
  per-candidate `beta_hat`, `sigma_D`, `family_id` and a `meta` list
  (the multiplier law, `B`, the resolved `c_cons` / `p_star` / effect
  floor, and the field's scaling diagnostics). It is what
  [`fs_declaration_calibration()`](https://larry-leon.github.io/forestsearch/reference/fs_declaration_calibration.md)
  reads. It re-reads the influence matrix and the multipliers already
  drawn – no fit, no draw, no RNG consumption – so every other output is
  byte-identical whether or not it runs, and with `FALSE` nothing is
  computed or stored.

- keep_field_matrix:

  Logical (default `FALSE`); consulted only when
  `keep_declaration_field = TRUE`. When `TRUE`, `declaration_field` also
  keeps the full `draws x G` matrix `Zstar`, which the reduced-family
  diagnostic of
  [`fs_declaration_calibration()`](https://larry-leon.github.io/forestsearch/reference/fs_declaration_calibration.md)
  needs. Off by default because `G` can be large.

- declaration_c0:

  `NULL` (default) or a numeric vector of clinically specified,
  pre-specified protected null levels `c0`; consulted only when
  `keep_declaration_field = TRUE` (ignored otherwise). Each `c0` is on
  the natural scale of the consistency threshold `c2` (the HR itself on
  the survival path; the ratio for OR / RR / IRR; the difference for RD
  / MD) and is mapped to the comparison scale exactly as `c2` is
  ([`log()`](https://rdrr.io/r/base/Log.html) for ratio measures,
  identity otherwise). `c0 <= c2` is required. With
  `delta_g = (c_cons - c0_cmp) / sigma_D(g)`, `declaration_field` then
  also stores `Mstar_c0`, a `draws x K` matrix of the shifted family
  maxima `max_g { Zstar[b, g] - delta_g }`, column names
  `as.character(c0)`, and `meta$c0` / `meta$c0_cmp`. At `c0 = c2` the
  column equals `Mstar` exactly. On an identity-scale path a `c0` given
  on the wrong scale cannot be detected; on a ratio path a non-positive
  `c0` (a log supplied by mistake) errors. With `NULL` nothing is
  computed or stored.

- pconsistency.digits:

  Integer or `NULL` (default). The digits the search rounded `Pcons` to
  before comparing it with `p_star`; consulted only when the admission
  set carries a consistency floor. MR's consistency floor is built from
  the screen's effective threshold (`.fs_pcons_eff(p_star, digits)`,
  mapped to z by `qnorm((1 + .) / 2)`), not from the exact cutoff
  `qnorm((1 + p_star) / 2)`, so that MR re-selects under the admission
  rule that produced the selected subgroup.
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  passes its own `pconsistency.digits`; `NULL` falls back to the
  [`subgroup.consistency()`](https://larry-leon.github.io/forestsearch/reference/subgroup.consistency.md)
  default.

## Value

List with the selected index/label, `naive` and `debiased` estimates
(effect scale, with approximate 95% CIs), `selection_bias`,
`fixed_bias`, `selection_rate`, `mean_r`, `mean_r_c`, the `settings`
actually used (`t_confirm`, `confirm_rule`, `reselection`,
`selection_rule`, `multiplier`, `draws`), `harm_flag`, family/subgroup
sizes, and `timing_seconds`. The `debiased` element carries `se_ij`,
`se_wald`, `var_ij`, and `ij_source`; its CI uses the IJ SE under
`ci_method = "field"` (the default) and `"ij"` (the FB analogue in both
cases) and the robust SE under `"wald"`. When
`include_complement = TRUE`, a `complement` element carries the
complement subgroup's `naive`/`debiased` estimates and bias terms in the
same form, including its own IJ variance. The complement's de-biased CI
follows the same SE convention as the winner's: the IJ SE under
`ci_method = "ij"` and `"field"` (under `"field"` the complement, like
the `debiased` element, is identical to the `"ij"` output), the robust
SE under `"wald"`.

`mean_r` is the mean of the IJ residual \\r_b\\ over `ok_H` – exactly
the draws entering the variance, not all `draws`, since on an excluded
draw \\r_b\\ is not a meaningful quantity. **The invariant is that it is
zero by construction whenever both bias terms share a denominator**,
which is the convention the package implements: `selection_bias` and
`fixed_bias` both average over the draws that produced a winner, and the
IJ runs on that same set. A non-zero value therefore means the two terms
are being normalised differently somewhere – the defect corrected in
`dad0415`, whose signature was precisely a non-zero residual mean.
`mean_r_c` is the same quantity for the complement over `use_c`, and is
`NA_real_` when no complement was fit. Both are diagnostics of the
correction's internal consistency rather than properties of the
estimate, which is why they sit beside `selection_rate` and not inside
`debiased`.

Note the mixed scales in the flat `debiased` list:
`est`/`lower`/`upper`/ `lower_1s` are on the **effect** scale, while
`se`/`se_ij`/`se_wald`/ `var_ij` and the two bias terms are on the
**working** (log, for ratio measures) scale.

When the selected subgroup cannot be fit in the reconstructed family the
return is a short variant carrying only `selected_index` (`NA`),
`selected_label`, `harm_flag` (`NA`), `settings`, `note` and `n_family`;
consumers must tolerate that shape.

Under `return_reselection = TRUE` the full return additionally carries a
`reselection` element (see that argument); the short variant never does.

Under `ci_method = "field"` the full return additionally carries a
`field` element: `lambda_mean`, `lambda_sd`/`se_field` and the `Lambda*`
quantiles `q05/q25/q50/q75/q95/q025/q975` on the **working** scale;
`n_out_used`, `n_in_used_mean`, `R_out`, `R_in`, `seed_offset`,
`timing_seconds`; and, on the **effect** scale (the `debiased`
convention), the second-order point estimate `est2` (de-biased estimate
minus `lambda_mean`), the primary one-sided 95% bound `lower_1s`
(de-biased minus `q95`), the two-sided quantile interval
`lower_2s`/`upper_2s`, and the supplementary SE-type interval
`lower_se`/`upper_se` around `est2`.

Under `field_recovery = TRUE` the `field` element additionally carries
`recovery`: `sens_H` (the primary quantity – the mean share of the
identified patients that the re-selections retain), `ppv_H`, `sens_Hc`,
`ppv_Hc` and its alias `npv_Hc`; the containment quantiles `q10`, `q50`,
`q90` and `share_equal_1` (the share of draws whose re-selection
contains all of the observed subgroup); and the accounting `n_draws`,
`n_used`, `n_skipped`, `n_selected`, `n_all`, `timing_seconds`.
Descriptive only; absent when `field_recovery = FALSE`.

Under `field_complement = TRUE` (with `include_complement = TRUE`) the
`field` element additionally carries `complement`, the complement's own
field block in the same form and scales, inverted around the
complement's two-term de-biased estimate: `lambda_mean`,
`lambda_sd`/`se_field`, the seven `Lambda*c` quantiles; on the effect
scale `est2`, the **primary** one-sided 95% upper bound `upper_1s`
(de-biased minus `q05`), the one-sided lower bound `lower_1s`, the
two-sided `lower_2s`/`upper_2s` and the SE-type `lower_se`/`upper_se`;
the draw accounting `n_out_used`, `n_in_used_mean`,
`n_out_dropped_unfit`; the fit accounting `n_complement_fits` (distinct
complement fits, the multiplier stage's included), `n_new_fits`,
`share_draws_new_fit` (share of outer + inner draw-winner readings whose
candidate the multiplier stage had not fit); `R_out`, `R_in`,
`timing_seconds`. A `note` replaces the numbers when the selected
complement is unfit or fewer than 2 outer draws survive.

With the complement field on, `field` also carries `joint`: the
simultaneous (harm lower, complement upper) pair from the aligned outer
draws \\(\Lambda^\*\_r, \Lambda^{\*c}\_r)\\ – same multipliers, same
winners, no new draws. `gamma` is the equal-tail level (grid from
`alpha` down to `alpha / 2`, step 0.001; the largest value whose joint
probability \\P^\*(\Lambda^\* \le q\_{1-\gamma}, \Lambda^{\*c} \ge
q\_\gamma) \ge 1 - \alpha\\), `joint_prob` the achieved probability,
`lower_H`/`upper_Hc` the calibrated pair (effect scale),
`bonf_lower_H`/`bonf_upper_Hc` the Bonferroni pair at `gamma = alpha/2`
with its `bonf_joint_prob`, `corr` the draws' correlation,
`n_joint_draws`, and the `grid_gamma`/`grid_joint_prob` profile. The
marginal one-sided bounds are unchanged. The top-level `ij_residual`
records the residual that populated the reported IJ SEs.

## Details

This is **post-selection inference on a completed analysis**: it runs
after
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
has chosen the subgroup and cannot change that choice. It is the fast,
refit-free counterpart of the full bootstrap (FB,
[`forestsearch_bootstrap_dofuture()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)),
which re-runs the entire search in every replicate; MR perturbs the
influence contributions of a single fit instead. MR approximates FB to
leading order and does not replace it.

## Alignment is assumed, not checked here

This is an engine-level entry point. Its arguments are a candidate
family, a specification, and a selected membership vector – the
identifier configuration that produced them is not visible, so the
alignment conditions MR requires (selection ranking on the inferential
coefficient \\\hat\beta(g)\\, and a fixed candidate family) cannot be
verified at this level and are assumed to have been established
upstream.

They are enforced by `.validate_mr_configuration()` at the three
configuration-visible entry points –
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
under `mr_inference = TRUE`, and
[`forestsearch_bootstrap_dofuture()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)
and
[`forestsearch_Kfold()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_Kfold.md)
under `mr_in_replicates = TRUE`. Calling `fs_mr_inference()` directly
bypasses those guards: a family ranked on DINA's native tau-hat, on
GRF's doubly-robust score, or on a GRF policy-tree objective will still
produce numbers, but they do not de-bias the reported effect.

## See also

[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
for the `mr_inference` switch and the vocabulary section;
[`mr_estimates_table()`](https://larry-leon.github.io/forestsearch/reference/mr_estimates_table.md)
to render the result;
[`forestsearch_bootstrap_dofuture()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)
for the full bootstrap (FB) this approximates;
[`fs_fdr_report()`](https://larry-leon.github.io/forestsearch/reference/fs_fdr_report.md)
which sweeps `c_confirm` thresholds.
