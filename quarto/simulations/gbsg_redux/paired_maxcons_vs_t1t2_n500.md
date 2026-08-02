---
title: "Paired comparison: fs_maxcons (new) vs fs_t1_t2 (legacy), M1 h10 n = 500"
bibliography: []
---

# Bundles

| | legacy | new |
|---|---|---|
| File | `fs_t1_t2_m1_h10_knoise0_n500_combined_1_500.rds` | `fs_maxcons_fb_mr_m1_h10_knoise0_n500_combined_1_500.rds` |
| `sg_focus` | `"eff"` | `"maxcons"` |
| `n_sample` / `n_sims` / `nb_boots` / draws | 500 / 500 / 300 / 5000 | identical |
| `seed_base` | 8316951 | 8316951 |
| `truth` | HR 0.685, marg_H 1.00, marg_Hc 0.657, cde_H 1.00, cde_Hc 0.585 | identical |
| Batches | 4 (1–20, 21–100, 101–200, 201–500) | 3 (1–20, 21–300, 301–500) |

Same design cell, same seed, different batch decomposition. `sim_id` sets are
identical, so all 500 replicates pair one-to-one.

# 1. Detection is bit-identical

|  | new: not detected | new: detected |
|---|---|---|
| **legacy: not detected** | 154 | 0 |
| **legacy: detected** | 0 | 346 |

Zero flips in either direction. True harm $Q$ also matches exactly on all 500
replicates. The detection decision and the underlying trials are untouched by the
revisions.

# 2. Among detections, 81.8% select the identical subgroup

| | count | share |
|---|---|---|
| Identical `label` and `sg_def` | 283 | 81.8% |
| Different rule selected | 63 | 18.2% |

# 3. On the 283 identical-selection replicates, MR is exactly reproduced

Maximum absolute difference across all 283 paired replicates:

| Quantity | max abs diff |
|---|---|
| `mr_H_est` vs `t2_H_est` | 6.7 × 10⁻¹⁶ |
| `mr_H_lo` vs `t2_H_lo` | 5.6 × 10⁻¹⁶ |
| `mr_H_hi` vs `t2_H_hi` | 3.6 × 10⁻¹⁵ |
| `mr_H_se_ij` vs `t2_H_se_ij` | 3.9 × 10⁻¹⁶ |
| `mr_Hc_est` vs `t2_Hc_est` | 1.7 × 10⁻¹⁵ |
| naive, oracle, $\beta(\widehat H)$ | ≤ 8.9 × 10⁻¹⁶ |

These are floating-point noise, not differences. **Given the same selected subgroup
and the same data, the renamed MR path returns numerically the same object as the
legacy gate.** The `debias_gate` → `mr_inference` rename, the `mr_inference_args`
restructuring, and `mr_draws` are confirmed inert on the estimation side.

# 4. FB carries small residual resampling variation

Same 283 replicates, same selected subgroup:

| Quantity | mean abs diff | max abs diff |
|---|---|---|
| `fb_H_est` vs `t1_H_est` | 0.0145 | 0.0438 |
| `fb_H_lo` vs `t1_H_lo` | 0.0138 | 0.0748 |
| `fb_H_hi` vs `t1_H_hi` | 0.0577 | 0.2950 |

Correlation 0.9990 on the point estimate. This is the expected signature of B = 300
bootstrap draws landing differently, not a systematic shift — it is centred, and the
upper bound moves most because it is the least stable end of the interval. Worth
noting for reproducibility purposes: FB is **not** exactly reproducible across these
two runs even when selection is held fixed, whereas MR is.

# 5. The 63 changed replicates select systematically larger regions

| | legacy | new |
|---|---|---|
| mean $|\widehat H|$ (changed reps) | 77.0 | 105.7 |
| median $|\widehat H|$ (changed reps) | 72 | 94 |
| quartiles | 65.5 / 81 | 75.5 / 124.5 |
| direction | — | 44 larger, 7 smaller, 12 same size |
| mean Δ | — | +28.6 (median +24) |
| terms per rule | 2.0 | 2.0 |

Unchanged replicates sit at mean 83.5 in both runs, so the entire marginal shift
(82 → 88) is carried by this 18%.

Covariate overlap between the old and new rule on the changed replicates: 12 share no
covariate, 29 share one, 22 share both (same two variables, different cut or
direction). Rule complexity is unchanged at two terms throughout — the search is not
returning more elaborate rules, it is landing on different, more inclusive ones.

Examples, largest Δ first:

| `sim_id` | legacy rule | $|\widehat H|$ | new rule | $|\widehat H|$ |
|---|---|---|---|---|
| 57 | `q7.1 & q5.0` | 89 | `q7.1 & q15.1` | 220 |
| 395 | `q13.1 & q18.1` | 64 | `q25.0 & q23.1` | 177 |
| 228 | `q12.1 & q17.1` | 68 | `q29.1 & q28.0` | 178 |
| 247 | `q26.1 & q22.0` | 72 | `q27.0 & q26.1` | 169 |
| 315 | `q14.1 & q21.1` | 81 | `q18.0 & q7.1` | 172 |

Classification on the changed subset moves as the larger regions imply: sensitivity
0.496 → 0.578, specificity 0.894 → 0.840, PPV 0.419 → 0.385, NPV 0.926 → 0.933. On
the unchanged subset sensitivity is 0.367 in both.

# 6. Changes are not batch-clustered

Share of detected replicates that changed selection, by `sim_id` segment (segments
chosen to straddle both batch decompositions):

| segment | 1–20 | 21–100 | 101–200 | 201–300 | 301–500 |
|---|---|---|---|---|---|
| changed | 16.7% | 19.6% | 18.9% | 6.9% | 23.4% |
| n changed / n detected | 2/12 | 10/51 | 14/74 | 5/72 | 32/137 |

Broadly uniform at roughly 17–23%, with a dip at 201–300. The differing batch
boundaries (legacy splits at 100/200; new splits at 300) do not line up with where
changes occur, so batch decomposition is not the driver. The 201–300 dip rests on
5 of 72 and I would not build on it without more replicates.

# 7. Coverage of $\beta(\widehat H)$, split by subset

| subset | MR legacy | MR new | FB legacy | FB new |
|---|---|---|---|---|
| unchanged (283) | 96.1% | 96.1% | 78.4% | 78.1% |
| changed (63) | 100% | 100% | 87.3% | 79.4% |

MR coverage is unaffected on both subsets. FB coverage on the changed subset drops
87.3% → 79.4%, but that compares intervals around *different* selected regions
against *different* per-replicate targets, so it is not a like-for-like read; the
FB shortfall against $\beta(\widehat H)$ is present in both runs either way.

# 8. Timings: confounded by batch, no code-level conclusion available

Legacy `t1_secs`/`t2_secs` map to new `fb_secs`/`fit_mr_secs` — same measured
quantities. The means agree closely, but that agreement is coincidental.

| | legacy | new |
|---|---|---|
| FB bootstrap, mean s/rep (detected) | 166.4 | 167.5 |
| fit+MR, mean s/rep (all 500) | 12.78 | 12.66 |

Per-replicate cost is nearly constant *within* a batch and varies more than threefold
*between* batches:

| legacy batch | 1–20 | 21–100 | 101–200 | 201–500 |
|---|---|---|---|---|
| FB mean s (SD) | 113.3 (3.8) | **372.3** (10.4) | 179.8 (6.2) | 114.5 (3.7) |

| new batch | 1–20 | 21–300 | 301–500 |
|---|---|---|---|
| FB mean s (SD) | 187.0 (2.1) | 187.8 (2.3) | **136.5** (8.3) |

Batch identity explains **99.6%** of FB timing variance in the legacy run and **95.4%**
in the new run (70% and 75% for fit+MR). Three checks confirm the replicate itself
contributes almost nothing: FB seconds correlate with $|\widehat H|$ at r = 0.00
(legacy) and 0.03 (new); fit+MR at 0.01 and 0.02; and the paired per-replicate
correlation of fit+MR across the two runs is **−0.23** — same trial, same work,
unrelated timings.

Because the batch decompositions differ (legacy 20/100/200, new 20/300), no two
batches were rendered under comparable conditions, and the confound is large enough to
invert conclusions. On the overlapping 301–500 range, fit+MR reads 15.98 → 5.02 s and
FB reads 114.4 → 136.5 s: naively, one revision made fit+MR 3.2× faster and FB 19%
slower. Both are the batch effect.

What can be said indirectly: §3 shows MR reproduces to machine precision at fixed
selection, so it is doing identical arithmetic and its cost cannot have changed
materially. Any `fit_mr_secs` difference here is environmental.

Two consequences worth noting. The wall-clock projections in both reports (568.6 vs
571.7 CPU-hours at 10,000 sims) extrapolate from these means and inherit the
confound — the `wall (calibrated)` column is empty in combine mode, so nothing
rescales them. And a real timing comparison has to be designed as one: same machine,
same session, identical batch decomposition, back to back. Fifty replicates would
suffice given within-batch SD of 2–10 s. For MR cost specifically,
`fs.est$mr_inference$timing_seconds` isolates it from the enclosing `forestsearch()`
fit, which `fit_mr_secs` conflates.

# 9. Where this leaves attribution

The difference is now localised precisely: **entirely upstream of MR, in which rule
the search returns, on 18% of detected replicates, in the direction of larger
regions.** Everything downstream — detection, MR estimation, MR intervals, MR
coverage — is either bit-identical or unchanged within resampling noise.

That rules out the MR-side renames as a cause and narrows the candidates to those
touching the search itself:

1. **`stop_threshold = NULL` passed explicitly.** Still the leading candidate. The
   setup comment records that `forestsearch()` resets it for the canonical family
   "WITHOUT warning unless supplied explicitly," so explicit supply is a distinct
   branch. Later or absent stopping would let the search continue past the point the
   legacy run halted, and continuing a consistency-driven search generally reaches
   more inclusive regions — which is the observed direction. It also fits the
   partial-penetrance pattern: stopping only binds on replicates where the search
   would otherwise have terminated early, leaving the other 82% untouched.
2. **`effect_neighborhood = 0.10` passed explicitly.** Same explicit-vs-missing
   distinction; documented as inert under `maxcons`.
3. **`sg_focus` alias handling** — whether `"maxcons"` reaches the same canonical key
   and the same tie-break as `"eff"`/`"hr"`. The 22 replicates that keep both
   covariates but change cut/direction are the ones most consistent with a tie-break
   difference.
4. **`future` 1.70.0 → 1.75.0** affecting the consistency-screen resample draws. If
   these were perturbed I would expect changes spread more evenly and in no
   particular size direction; the strong skew toward larger regions (44 vs 7) argues
   against RNG drift being the main mechanism, though it may contribute.

# 10. Suggested next step

The decisive test is now cheap and targeted, because the affected replicates are
identified. Re-run just the 63 changed `sim_id`s under the new pipeline with
`stop_threshold` and `effect_neighborhood` **removed from the call entirely** (not
set to `NULL` / 0.10) and compare `label` and `n_harm` against the legacy bundle:

- all 63 revert → candidates 1–2 confirmed, and the two arguments should be dropped
  from the call rather than pinned;
- a subset reverts → mixed cause, and the residue isolates candidate 3 or 4;
- none revert → attention moves to `.normalize_sg_focus()` and the screen's RNG.

Worth deciding separately, since it is a reporting choice rather than a code question:
whether the intended `maxcons` behaviour is the legacy one or the new one. The new
run finds larger regions with higher sensitivity and lower specificity, which is not
self-evidently worse — but it is a different identifier, and the JASA simulation
tables describe the legacy one.
