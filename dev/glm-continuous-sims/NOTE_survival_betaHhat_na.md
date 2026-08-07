# The survival `betaHhat` NA failure on noise-covariate sweeps

Found while answering a different question: whether the disjunction defect
corrected in `betaHhat_truth_glm.R` had ever fired in committed results. It had
not. This had.

Measurement only. Nothing is patched, and the fix proposed at the end is **not**
implemented.

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

## Proposed fix — not implemented

Four options, narrowest first.

1. **Make partial resolution return all-NA.** In `betaHhat_one*()`, after
   `get_dfpred()`, reject any membership vector containing `NA` and return the
   all-NA record. One condition per module. This does not recover any
   replicate, but it converts 13,777 silently-wrong complement values into
   honest NAs, which the parity guard then catches.
2. **Report the exclusion rate.** Whatever else is done, the fraction of
   replicates with no computable $\beta(\widehat H)$ belongs in the bundle
   `meta` and in the coverage table, so a reader sees the denominator.
3. **Make `get_dfpred()` loud on a missing column.** It currently warns and
   yields `NA` membership. An error would have surfaced this the first time it
   ran. This is an `R/` change with callers beyond these modules, so it needs
   its own impact review.
4. **Decide what $\beta(\widehat H)$ means for a per-replicate covariate.** The
   estimand genuinely does not exist for these rules. Either the noise
   covariates become population characteristics defined once on `df_super` and
   subsetted into each trial — which makes the target well-defined and is the
   principled fix — or noise-referencing replicates are excluded by design and
   reported as such under (2). This is a DGM-contract decision, not a bug fix.

**Recommendation: (1) and (2) now, (4) before these two grids are re-run or
cited.** (1) is small and stops wrong numbers from being published; (2) makes
the loss visible; (4) is the real question and should not be settled by whoever
happens to touch the file next. (3) is worth doing but is a separate change
with its own blast radius.

Nothing above should be read as a defect in the identification methods. The
sweeps' `detected`, `sens`, `ppv` and the oracle/naive/MR estimates are
untouched — this is confined to the $\beta(\widehat H)$ target column and the
`C_betaHhat` coverage row computed from it.
