---
bibliography: []
---

# Note: the MR admission set, and what its correction changed

**Status:** fixed. One further observation recorded at the end --
characterisation, not an outstanding defect.
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
only where the comparison required it (the `maxeff` FB/MR observation below,
at 100 replicates -- a detection question, not an estimation one); everything
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

## Observation: the FB/MR gap under maxeff, and why it is expected

**This is an observation to characterise, not a defect awaiting a fix.**
Nothing here should be read as a proposal to change the bootstrap.

Theorem 2 gives first-order agreement between MR and the full bootstrap
**under a fixed family**. Under `sg_focus = "maxeff"` a gap is present, and the
admission fix did not change it:

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

So the admission set was not the cause of this gap. The fix remains necessary
and correct; it simply does not happen to move this number.

### Why a gap is expected here

**The full bootstrap refits everything on each resample by design, and that is
correct.** It re-runs the whole pipeline per replicate -- re-deriving cut
locations, refitting front ends, re-enumerating candidates and re-selecting.
Freezing the cut grid is **not** something to implement, and this note does not
propose it.

Whether a given run has a fixed family is a property of the **configuration**,
not something the code enforces. Two distinct things go by "fixed" and should
not be conflated:

- `family_status = "no-front-end"` (formerly `"fixed"`) records that no model-based front end
  (`use_lasso` / `use_grf` / `use_dina`) shapes the family on the observed
  data. That is what the MR alignment condition checks.
- It does **not** mean the family is invariant under resampling. With
  `cut_type = "default"` the cut locations are sample quantiles, so they move
  with each resample. Only literal cut expressions supplied via `conf_force`
  are resample-invariant.

The validation above was run with all three front ends `FALSE`, so
`family_status = "no-front-end"` -- yet the cut grid still regenerates per replicate,
because the cuts are quantiles. Under a regenerating family FB and MR target
**different estimands**: FB is unconditional over the family-generating
process, MR is conditional on the realized family. That is the
conditional-versus-marginal distinction of Section 5. A gap is therefore the
expected consequence of the configuration, not an anomaly.

What remains genuinely open is only its **size** -- characterising how the gap
behaves as the family-generating step is made more or less stable (for example
`conf_force` cut expressions versus default quantile cuts) -- rather than
whether it should be there at all.

For the record, the code confirming that FB regenerates rather than freezes:
`bootstrap_analysis_dofuture.R:377-378` drops the derived cut columns from the
resampled frame and lets the inner `forestsearch()` rebuild them via
`get_FSdata()`; `:575-583` sets `grf_res`, `grf_cuts`, `dina_res` and
`dina_cuts` to `NULL` per replicate "so each replicate re-fits ... and re-runs
its selection from scratch". This is the intended behaviour of a bootstrap that
replays the whole pipeline.
