# The survival `betaHhat` NA failure on noise-covariate sweeps

Found while answering a different question: whether the disjunction defect
corrected in `betaHhat_truth_glm.R` had ever fired in committed results. It had
not. This had.

Originally measurement only. The recommended fix has since been implemented and
piloted — see **Status: implemented and piloted** at the end. The measurements
and the mechanism below are unchanged and describe the state before that fix.

## Summary

In the two GBSG noise-covariate sweeps, `betaHhat_truth.R` cannot resolve any
realized rule that references a `noise*` variable, because the evaluation frame
has no such column. The failure is silent. It costs **28.4%-77.2%** of
replicates in the `H` block of `C_betaHhat` coverage, and — worse — it puts
**13,777 finite but wrong** values into the `Hc` block, where they pass every
finiteness check and are compared against intervals for a different region.

This is live in committed artifacts, not latent.

## Mechanism

Two facts about `mr_coverage_sweep_h10_knoise3.qmd` that do not fit together.

**The evaluation frame is built once, from the DGM** (`:191`):

```r
eval_df <- build_eval_frame(dgm, analysis_time = analysis_time,
                            cens_adjust = cens_adjust, eval_seed = 20260628L)
```

**The noise covariates are built per replicate, on the trial draw** (`:283-288`),
and appended to the analysis candidate set:

```r
if (k_random_noise > 0) {
  noise_names <- paste0("noise", seq_len(k_random_noise))
  set.seed(seed_base + sim_id + 1000000L)
  for (nm in noise_names) df[[nm]] <- rnorm(nrow(df))
  confs <- c(confs, noise_names)
}
```

So `forestsearch()` may select on `noise1`, and `attach_betaHhat()` then scores
that rule against a frame where `noise1` does not exist.

**The deeper point, which decides the fix.** The noise variables are re-drawn
from `set.seed(seed_base + sim_id + 1000000L)` for *every* replicate. They are
not a population characteristic that the evaluation frame happens to omit —
they have no population definition at all. $\beta(\widehat H)$ is therefore
**undefined** for a rule referencing one, not merely uncomputable. Adding
`noise*` columns to the evaluation frame would not fix this; it would invent a
population realization and score against it.

## Reproduction

Real DGM, real evaluation frame ($n = 100{,}000$, `eval_seed = 20260628`), and
a real NA-target rule taken from
`mr_sweep/m1_h10_knoise3new_s1000/dina_mr_n1000_res.rds`:

```
rule  : {pgr <= 11} & {noise1 >= 0.27207708810350478}
split : '{pgr <= 11}', '{noise1 >= 0.27207708810350478}'

get_dfpred(eval_df, sg):
WARNING: Column 'noise1' not found in data frame
returned a frame; treat.recommend is integer with 32228 NA of 100000

betaHhat_one(rule, eval_df):
 betaHhat_H betaHhat_Hc     nH_eval    nHc_eval
         NA   0.6310372          NA          NA
```

Note what did **not** happen: no error, and `get_dfpred()` did not return
`NULL`. It returned a frame with a partially-`NA` membership vector, so
`betaHhat_one()`'s guard

```r
if (is.null(pred) || is.null(pred$treat.recommend)) return(all-NA)
```

does not fire. The NAs propagate instead:

| step | result |
|---|---|
| `{pgr <= 11}` | TRUE for 32,228 of 100,000 |
| `{noise1 >= ...}` | `NA` everywhere (column absent) |
| `in_harm <- TRUE & NA` | `NA` (32,228) |
| `in_harm <- FALSE & NA` | **`FALSE`** (67,772) |
| `treat.recommend` | 32,228 `NA`, 67,772 ones, **0 zeros** |
| `inH <- treat.recommend == 0L` | 0 TRUE, 67,772 FALSE, 32,228 `NA` |

`eval_df[inH, ]` selects only all-`NA` rows, so the event guard fires and
`betaHhat_H` is `NA` — correct by accident. But `eval_df[!inH, ]` selects the
67,772 real rows plus 32,228 all-`NA` rows, `coxph()` drops the latter, and the
fit lands on `{pgr > 11}`:

```
  HR on {pgr > 11}     = 0.6310372
  reported betaHhat_Hc = 0.6310372
  identical to 7 dp    : TRUE
```

**The reported complement target is the complement of the rule with the noise
clause silently dropped.** It is a finite, plausible hazard ratio for a region
that is not $\widehat H^{c}$.

## Scope

Two grids, both noise-covariate sweeps. All other GBSG sweeps
(`k_random_noise = 0`) are unaffected: **274 of 333 bundles have a zero NA
rate**.

Diagnosis is unambiguous — of rules whose target is NA, 18,149 of 18,149
(100.0%) mention a `noise*` variable; of rules that scored, 0 of 165,783 (0.0%)
do.

Replicates lost from `C_betaHhat` relative to `C_dagger`, block `H`
(identical across the three estimators; 42 cells, 18,064 replicate-losses):

| grid | method | n=500 | 750 | 1000 | 1250 | 1500 | 1750 | 2000 |
|---|---|---|---|---|---|---|---|---|
| knoise3new | consistency | 54.1% | 50.8% | 47.0% | 46.6% | 40.7% | 34.5% | 28.4% |
| knoise3new | dina | 57.6% | 52.3% | 51.7% | 45.1% | 43.6% | 43.8% | 39.4% |
| knoise3new | grf | 59.9% | 53.7% | 53.4% | 48.7% | 47.5% | 45.0% | 37.5% |
| knoise6new | consistency | 74.1% | 69.9% | 67.8% | 64.7% | 60.2% | 51.5% | 47.5% |
| knoise6new | dina | 76.4% | 74.2% | 69.6% | 66.8% | 60.8% | 58.2% | 54.0% |
| knoise6new | grf | 77.2% | 72.7% | 70.2% | 66.4% | 65.2% | 61.8% | 57.6% |

Monotone decreasing in $n$ and worse with six noise covariates than three, both
as expected: a larger trial and a smaller candidate set both make a noise
variable less likely to survive selection.

## Consequence for coverage

`.coverage_meta()` computes `ok <- is.finite(target) & is.finite(lo) &
is.finite(hi)` and `n_eff <- sum(ok)`. There is no error and no warning. So:

**Block `H` — a biased denominator.** Between a quarter and three quarters of
replicates silently leave the denominator. The survivors are exactly the
replicates that selected a noise-free rule, which is the subset closest to the
true subgroup. `C_betaHhat` coverage in these two grids is computed on a
strongly non-random subsample and should not be read as published.

**Block `Hc` — contaminated values, which is worse.** Of the 18,149 NA-`H`
rows, only 4,372 also have `betaHhat_Hc` NA. The other **13,777** carry a
finite hazard ratio for the wrong region. Those are not dropped: they enter
`n_eff` and are compared against `Hc` intervals. That is why the `Hc` block
shows only a 2.9%-25.1% loss while `H` shows 28.4%-77.2% — the difference is
not health, it is contamination.

**The `n_eff` parity guard added in `e6f6024` catches the first and not the
second.** It compares `n_eff` across targets, so a target that is finite-but-wrong
is invisible to it. Recorded here so the guard is not over-trusted.

## Proposed fix

*(The recommended option below is now implemented and piloted; the narrower
options (1)-(3) remain outstanding. See the status section at the end.)*

### Recommended: generate the noise once on `df_super`

**A noise baseline factor is a patient attribute.** Today it is not one. The
current code draws the trial and *then* attaches noise:

```r
df <- tryCatch(                                      #  :273-276
  simulate_from_dgm(dgm, n = n_sample, ...), ...)
...
set.seed(seed_base + sim_id + 1000000L)              #  :286
for (nm in noise_names) df[[nm]] <- rnorm(nrow(df))  #  :287
```

so the same `df_super` subject carries **different** noise in different
replicates. That is replicate-level randomness wearing a covariate's name.
Generating the noise once, on `dgm$df_super`, and letting the sampling carry it
makes it an attribute, and everything else follows.

Four facts establish that this works, all checked against source:

1. **The insertion point exists and is upstream of the eval frame.**
   `dgm$df_super` is assembled at `R/setup_gbsg_dgm.R:107-121` — `df_s <-
   dgm$df_super_rand`, four `.rename()` calls, then `dgm$df_super <- df_s`. It
   is an ordinary data frame. Noise can be attached in the sweep `.qmd`
   immediately after `setup_gbsg_dgm()` and before `build_eval_frame()`, so
   **no `R/` change is required**.
2. **Arbitrary extra columns survive into the trial unchanged.**
   `simulate_from_dgm()` subsets whole rows — `df_sim <- df_super[idx_sample, ]`
   (`R/simulate_from_dgm.R:174`) — and returns `df_sim` with no column
   whitelist. Anything on `df_super` arrives in the trial with the same values.
3. **The evaluation frame gets the identical values by construction.**
   `build_eval_frame()` is `simulate_from_dgm(dgm, n = nrow(dgm$df_super),
   replace = FALSE, ...)`, i.e. every subject exactly once. Trial and eval frame
   then read the same `noise1` for the same subject.
4. **$\beta(\widehat H)$ becomes exact, not approximate.** With the noise a
   fixed population column, a noise-referencing rule resolves on `df_super`
   like any other, and the target is an exact finite computation over the
   population rather than an estimate. This removes the Monte-Carlo-error
   objection to the target entirely — so the separate question of whether the
   eval frame is large enough to read as "population" is moot for these rules.

**One behavioural change to note, correct under this scheme but a change.**
`simulate_from_dgm()` defaults to `replace = TRUE` (`:135`, `:173`), so a
subject drawn twice into one trial appears twice — and under the new scheme
carries **identical** noise both times, because it is the same row copied.
Today each row gets an independent draw. That is the right behaviour for an
attribute (a patient's covariate does not change because they were resampled),
but it means within-trial noise is no longer i.i.d. across rows, and any
analysis that assumed it was should be re-checked.

### The narrower options, still worth doing

1. **Make partial resolution return all-NA.** In `betaHhat_one*()`, after
   `get_dfpred()`, reject any membership vector containing `NA` and return the
   all-NA record. This does not recover any replicate, but it converts 13,777
   silently-wrong complement values into honest NAs, which the parity guard
   then catches. Worth keeping even after the recommended fix, as the guard
   against the *next* unresolvable-rule cause. Note the MD sibling fails
   differently on the same input: `betaHhat_one_md()` **errors** ("missing
   value where TRUE/FALSE needed") because `any(idx)` on an all-NA vector is
   `NA`. Loud beats silent, but an uncaught error aborts a sweep; neither
   module currently returns the all-NA record it should.
2. **Report the exclusion rate.** The fraction of replicates with no computable
   $\beta(\widehat H)$ belongs in the bundle `meta` and the coverage table.
   Under the recommended fix this should be zero, which makes it a useful
   regression check rather than dead weight.
3. **Make `get_dfpred()` loud on a missing column.** It warns and yields `NA`
   membership; an error would have surfaced this the first time it ran. An
   `R/` change with callers well beyond these modules, so it needs its own
   impact review.

**Recommendation: the population-noise scheme, then (1) and (2).** The scheme
is the only option that makes the estimand exist — the others manage the
symptom. It needs no `R/` change, and (1) and (2) then serve as guards rather
than as the fix. (3) remains worth doing separately.

## Blast radius, if the scheme is adopted

**Sweep documents.** 56 `.qmd` files set `k_random_noise > 0` — 34 under
`quarto/simulations/gbsg/`, 22 under `gbsg_redux/` (21 at `k = 3`, 34 at
`k = 5`, 1 at `k = 6`). Only two of them are `mr_coverage_sweep_*` documents
that produce $\beta(\widehat H)$ targets: `mr_coverage_sweep_h10_knoise3.qmd`
and `mr_coverage_sweep_h10_knoise6.qmd`. The rest are batch/combine documents
on other pipelines.

The per-replicate noise block would be **removed, not kept as a fallback**.
Keeping both paths would leave two definitions of what `noise1` means,
selected by whichever branch ran — which is the ambiguity that produced this
defect. `confs <- c(confs, noise_names)` stays; only the `set.seed()` +
`rnorm()` pair goes.

**Committed bundles needing regeneration.** Three directories carry affected
rows:

| directory | rows with a rule | NA targets | share |
|---|---|---|---|
| `gbsg_redux/mr_sweep/m1_h10_knoise3new_s1000` | 15,191 | 7,272 | 47.9% |
| `gbsg_redux/mr_sweep/m1_h10_knoise6new_s1000` | 16,320 | 10,792 | 66.1% |
| `gbsg_redux/results` | 10,632 | 85 | 0.8% |

Totalling 18,149, matching the global count. In every directory the NA count
equals the noise-referencing-rule count exactly.

`gbsg/mr_sweep/m1_h20_knoise5_s1000` is the one other `k > 0` bundle directory
and is **not** affected: its 22 bundles predate the machinery and carry neither
`sg_def` nor `betaHhat_H`.

## A second copy of the survival module

`quarto/simulations/gbsg_redux/betaHhat_truth.R` and
`dev/identifier-alignment/rerun/betaHhat_truth.R` are **byte-identical** — 160
lines, 9057 bytes each. Both carry the split-first disjunction pattern
corrected in the binary module at `e6f6024`, and both would need any fix
applied. Worth resolving whether the `rerun/` copy should exist at all rather
than fixing the same file twice.

Nothing above should be read as a defect in the identification methods. The
sweeps' `detected`, `sens`, `ppv` and the oracle/naive/MR estimates are
untouched — this is confined to the $\beta(\widehat H)$ target column and the
`C_betaHhat` coverage row computed from it.

## Status: implemented and piloted (2026-08-07)

The recommended population-noise scheme above is **implemented**. Two documents
changed, no `R/` change:

```
quarto/simulations/gbsg_redux/mr_coverage_sweep_h10_knoise3.qmd
quarto/simulations/gbsg_redux/mr_coverage_sweep_h10_knoise6.qmd
```

Noise is drawn once on `dgm$df_super` at `noise_seed = 20260807L` — fixed,
distinct from `eval_seed` (20260628L), outside the trial band
`seed_base + sim_id` — immediately after `setup_gbsg_dgm()` and **before**
`build_eval_frame()`. The per-replicate re-draw is removed, not kept as a
fallback. Bundles record `noise_scheme = "population"` and `noise_seed`, so a
reader can tell new bundles from old ones (which are implicitly
`"replicate"`).

### Mechanism, proven

A `.pilot_id` column attached alongside the noise rides the same whole-row
carry, so the id proves the mechanism twice. Joining on it, all three noise
columns and the id are `identical()` across `df_super`, a drawn trial, and the
evaluation frame. The `replace = TRUE` consequence was measured rather than
assumed: 3 subjects were drawn twice into a 500-row trial and carried the same
`noise1` at both rows. An independent redraw at a different seed correctly does
not match.

### Pilot: knoise3, one cell, n = 500, 30 replicates, 102 s

| check | result |
|---|---|
| 1. leak | 11 of 21 detections name a `noise*` variable (52.4%) |
| 2. `betaHhat_H` NA with a rule | **0** — including 0 of the 11 noise-referencing rules, which were **100% NA** under the old scheme |
| 3. partition | exact on all 21 distinct rules, `nH + nHc = 100000` |
| 4. `fs_betaHhat_neff_parity()` | **PASS**, `n_eff` 21 / 21 against `C_dagger` |
| 5. contamination signature | **gone** — `betaHhat_Hc` 0.6768585260 vs the real-part-only complement 0.7211097134, differing by 4.4e-02 where the old bug forced equality to 7 dp; region sizes 17,395 vs 72,059, a factor of four |
| 6. leak comparability | 52.4% vs the committed 54.1% at this cell, **-0.16 SE** — statistically indistinguishable, as expected since noise remains independent of outcome |

Check 5 is the decisive one: under the old bug the noise clause was silently
dropped, so the reported complement was the complement of the real part alone
and the two agreed to 7 dp. They now differ, because the true complement
includes the noise-excluded slice of the real region.

**knoise6 smoke, 2026-08-07** — same settings (one cell, n = 500, 30
replicates, 174 s), same six checks, all pass: leak 18 of 25 detections
(72.0%) against the committed 74.1% at this cell (-0.23 SE); **0** NA
`betaHhat_H` among rows with a rule, including 0 of the 18 noise-referencing
ones; partition exact on all 25 distinct rules; parity PASS at `n_eff` 25/25;
contamination signature gone (0.6804227928 vs 0.6865144474, differing by
6.1e-03, regions 14,995 vs 30,113). Bundle is scratch, not committed.

The pilot bundle is scratch and is **not committed**. Full regeneration of the
noise sweeps is a separate campaign and has not been run.

### Still queued

Sweep re-pointing at `fs_attach_betaHhat()` — which will also let partition
checks read `nH_eval`/`nHc_eval` from the bundle rather than re-deriving them,
as this pilot had to, because the shim strips those columns to preserve the
pre-consolidation frame shape. Then the attach-ITT semantics change, then shim
deletion.
