---
bibliography: []
---

# Note: the MR admission set, and what its correction changed

**Status:** fixed. One item recorded as **OPEN** at the end.
**Change:** `.fs_resolve_admission()` -- admission resolved once and carried to
MR, replacing MR's reconstruction of `t_g` from raw parameters.
**Validated:** 2026-08-05, ACTG175 binary/OR, fixed family
(`use_lasso`/`use_grf`/`use_dina` all `FALSE`).

---

## What was wrong

Multiplier resampling linearizes the identifier's selection **map**, which is a
ranking *and* a domain. MR is valid only if the set it re-selects over is the
set the identifier selected over. MR was rebuilding that domain at its own call
site:

```r
t_g <- pmax(c_screen, c_consistency + z * sigma_D)
```

so the identifier decided its floors in one place and MR re-derived them in
another. They drifted, in two opposite directions:

- **`sg_focus = "maxeff"`, consistency engine.** The identifier applies no
  auxiliary selection condition (Guo and He's argmax primitive), but MR was
  still passed `c_screen = log(hr.threshold)` -- MR admitted *less* than the
  identifier did.
- **DINA and GRF.** Call sites passed `c_consistency = 0` while leaving
  `p_star = pconsistency.threshold`, so `t_g = pmax(c_screen, 1.645 * sigma_D)`.
  Note that `c_consistency = 0` does **not** neutralise the term: the
  `z * sigma_D` part survives. MR applied a consistency floor to engines that
  have no consistency screen -- MR admitted *less* than those identifiers did,
  by a different route.

## What the fix changed, measured

Before/after on ACTG175, fixed family, MR draws 2000. Full bootstrap was run
only where the comparison required it (see the open item below); everything
else is MR-versus-MR on the same fitted object, where a bootstrap adds nothing.

| config | MR before | MR after | change | sel-rate before -> after |
|---|---|---|---|---|
| consistency/maxeff | 1.82500 | 1.82500 | **0.00000** | 1.000 -> 1.000 |
| consistency/maxeffCons | 1.82242 | 1.82242 | 0.00000 | 0.992 -> 0.992 |
| consistency/hr | 1.78362 | 1.78362 | 0.00000 | 0.992 -> 0.992 |
| consistency/maxSG | 1.19842 | 1.19842 | 0.00000 | 0.992 -> 0.992 |
| consistency/minSG | 1.67356 | 1.67356 | 0.00000 | 0.992 -> 0.992 |
| consistency/hrMaxSG | 1.80932 | 1.80932 | 0.00000 | 0.992 -> 0.992 |
| consistency/hrMinSG | 1.79925 | 1.79925 | 0.00000 | 0.992 -> 0.992 |
| **dina/hr** | 0.35652 | **0.57330** | **+0.21679** | **0.051 -> 0.341** |
| dina/maxeff | 0.57330 | 0.57330 | 0.00000 | 0.341 -> 0.341 |
| **grf/hr** | 1.68887 | **1.94281** | **+0.25394** | **0.766 -> 0.981** |
| grf/maxeff | 1.94281 | 1.94281 | 0.00000 | 0.981 -> 0.981 |

Estimated selection bias fell correspondingly: DINA 1.021 -> 0.546, GRF
0.338 -> 0.198.

**The convergence is the check that could have failed.** After the fix
`dina/hr` equals `dina/maxeff` to all printed digits, and `grf/hr` equals
`grf/maxeff`. Before, those pairs differed *only* by accident: the identifier's
`maxeff` override block sets `pconsistency.threshold <- 0` for every method, so
DINA-under-`maxeff` happened to receive the correct admission set while
DINA-under-`hr` received the spurious floor. Once admission is resolved
structurally, the accident and the defect both disappear and the pairs agree.

## Consequence for DINA results predating this fix

**`dina/hr` had a selection rate of 0.051.** The spurious `z * sigma_D` floor
rejected the candidate in roughly 95% of multiplier draws, so MR's selection
bias -- and therefore its de-biased estimate and its IJ interval -- rested on
about 5% of the resampling distribution. That is not a small perturbation of a
well-estimated quantity; it is an estimate computed from a twentieth of the
draws it was meant to use.

This bears directly on:

- **The ACTG175 DINA rows** in `sim_analyses/actg175_table_payload.rds` and in
  the re-run under `rerun/`. Those runs used `sg_focus = "effMaxSG"`, not
  `hr`, so the exact 0.051 figure is specific to the configuration measured
  here -- but every DINA (and GRF) MR result produced before this fix carried
  the same spurious consistency term, and its severity depends on
  `sigma_D` relative to `c_screen` in that run. Their `mr_*` columns should be
  regenerated rather than trusted.
- **Any DINA or GRF MR result predating this fix**, in the simulations as well
  as the applications. Note this is independent of, and compounds with, the
  R 4.6.0 `stats::dfbeta()` issue recorded in
  `NOTE_actg175_naive_interval_rootcause.md`: that one inflated `sigma_D`
  itself, and an inflated `sigma_D` makes this floor bite *harder*, since the
  offending term is `z * sigma_D`.

The consistency engine is unaffected -- all seven foci are numerically
unchanged.

## OPEN: the maxeff FB/MR gap does not close

Theorem 2 gives first-order agreement between MR and the full bootstrap under a
fixed family. Under `sg_focus = "maxeff"` it is not observed, and the admission
fix did not change it:

| | FB | MR | gap (log scale) |
|---|---|---|---|
| before | 2.6309 | 1.8250 | 0.36574 |
| after | 2.6309 | 1.8250 | **0.36574** |

**Why the fix did not move it** (established, not assumed). Under `maxeff` the
identifier had *already* set `pconsistency.threshold <- 0`, so
`z = qnorm(0.5) = 0` and `t_g = pmax(c_screen, c_cons + 0 * sigma_D)` was
**constant** across candidates at `log(1.25)`. A constant floor cannot move an
argmax-by-effect while the passing set is non-empty, and the selection rate was
already 1.000, so no draw was dropped and the winner was the global argmax on
every draw. The old admission set was wrong in principle and inert in
arithmetic *for this focus, rule and dataset*. It would not be inert with a
selection rate below 1 (draws get dropped -- exactly DINA's case), with a
size-based rule (`maxSG`/`effMaxSG`, where a floor changes the pool rather than
the maximum), or with a non-zero `p_star` making `t_g` candidate-varying.

**So the admission set was not the cause of this gap.** The fix remains
necessary and correct; it simply does not happen to move this number.

**Hypothesis to test first, before anything else.** The full bootstrap may not
be the fixed-family bootstrap Theorem 2 refers to. If FB re-derives cut
locations from each resampled dataset, its candidate family regenerates per
replicate, and the FB/MR gap is then the conditional-versus-marginal estimand
distinction rather than a defect.

**The bootstrap does not freeze the cut grid.** Stated as fact, from the code,
without further investigation:

- `bootstrap_analysis_dofuture.R:377-378` -- "In the consistency/GRF/LASSO path,
  `confounders.candidate` are derived cut-columns that `forestsearch()` rebuilds
  internally on each resample, so they are dropped here and regenerated." The
  derived cut columns are removed from the resampled frame and rebuilt by the
  inner `forestsearch()` call, which re-runs `get_FSdata()` on the resampled
  data. With `cut_type = "default"`, cut locations are sample quantiles, so
  they are re-derived from each resample.
- `bootstrap_analysis_dofuture.R:575-583` -- `grf_res`, `grf_cuts`, `dina_res`
  and `dina_cuts` are all set to `NULL` per replicate, explicitly "so each
  replicate re-fits ... and re-runs its selection from scratch".

The one exception is `conf_force`, whose cut expressions are literal and
therefore identical on every resample.

This supports the hypothesis and should be the first thing checked when the
gap is taken up. It is **not** investigated further here, and no conclusion is
drawn about whether the gap is benign.
