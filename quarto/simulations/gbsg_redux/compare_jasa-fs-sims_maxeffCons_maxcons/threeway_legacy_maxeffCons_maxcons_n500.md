---
title: "Three-way comparison: t1_t2_legacy vs maxeffCons vs maxcons (M1, h10, n = 500)"
bibliography: []
---

# Bundles

All three are the same design cell and pair one-to-one on `sim_id`:
`n_sample` 500, `n_sims` 500, `nb_boots` 300, draws 5000, `seed_base` 8316951,
`consistency_method = "resample"`, `forestsearch 0.2.0`, and identical `truth`
(HR 0.685, marg_H 1.000, marg_Hc 0.657, cde_H 1.000, cde_Hc 0.585).

| arm | `sg_focus` | batches |
|---|---|---|
| `t1_t2_legacy` | `"eff"` | 4 (1–20, 21–100, 101–200, 201–500) |
| `maxeffCons` | `"maxeffCons"` | 3 (1–50, 51–300, 301–500) |
| `maxcons` | `"maxcons"` | 3 (1–20, 21–300, 301–500) |

Column sets: legacy uses `t1_*`/`t2_*`/`gate_ok`/`gate_draws`; both new arms use
`fb_*`/`mr_*`/`mr_ok`/`mr_draws` and are structurally identical to each other (52
columns in all three).

# 1. Detection is identical in all three

346 / 500 detections, and the agreement is replicate-by-replicate: legacy vs
maxeffCons, legacy vs maxcons, and maxeffCons vs maxcons all agree on 500 / 500.
Zero flips anywhere. Detection is untouched by either focus rule.

# 2. Selection: the two new arms differ from each other more than from legacy

Share of the 346 detections selecting the identical rule:

| pair | agreement |
|---|---|
| legacy vs `maxcons` | 283 (81.8%) |
| legacy vs `maxeffCons` | 281 (81.2%) |
| **`maxeffCons` vs `maxcons`** | **246 (71.1%)** |

The two new focus settings are genuinely distinct identifiers, not two spellings of
one rule. Their mutual disagreement (28.9%) exceeds either one's disagreement with
legacy.

Subgroup size, detected replicates:

| arm | mean | SD | median | range |
|---|---|---|---|---|
| legacy | 82.3 | 21.8 | 76 | [61, 184] |
| `maxeffCons` | **78.1** | **17.8** | 73 | [61, 167] |
| `maxcons` | **87.5** | **28.1** | 79 | [61, 220] |

They bracket legacy: `maxeffCons` is tighter and slightly more conservative,
`maxcons` is more inclusive with a longer right tail. Classification follows:

| | sens | spec | PPV | NPV |
|---|---|---|---|---|
| legacy | 0.390 | 0.867 | 0.315 | 0.909 |
| `maxeffCons` | 0.368 | 0.874 | 0.309 | 0.907 |
| `maxcons` | 0.405 | 0.858 | 0.309 | 0.910 |

Differences are modest — roughly ±0.02 on sensitivity, ±0.01 on specificity. The mean
conditional target is essentially unchanged (β(Ĥ) = 0.728 / 0.726 / 0.725).

# 3. MR at coincident selection: `maxeffCons` differs as expected, `maxcons` does not

## What MR depends on

The selection map $\mathcal{S}$ is *inside* the correction, not upstream of it. Each
multiplier draw re-selects a winner $\widehat H^{*}_{b}$ under $\mathcal{S}$; that
draw's perturbation $D_{\widehat H^{*}_{b}}(b)$ is what makes the selection-optimism
term positive on average, and it enters both the point correction and the per-draw
residual $r_b$ from which the IJ variance is built. So MR is a function of the
selection rule, not only of the realized $\widehat H$: **two arms with different
$\mathcal{S}$ should produce different MR estimates and intervals even on replicates
where the observed winner happens to coincide.**

## What the data show

Restricting to replicates where the arm selected the identical rule as legacy:

| | identical-selection reps | MR point est, max abs diff | MR upper bound, max abs diff |
|---|---|---|---|
| legacy vs `maxcons` | 283 | **6.7 × 10⁻¹⁶** | 3.6 × 10⁻¹⁵ |
| legacy vs `maxeffCons` | 281 | **2.9 × 10⁻¹** | 9.4 × 10⁻¹ |

`maxeffCons` departs from legacy, systematically:

- **0 of 281** replicates agree even to 10⁻¹⁰; 267 differ by more than 0.01, 179 by
  more than 0.05.
- Mean difference **+0.067** (legacy higher), median +0.064, and **247 of 281** move in
  the same direction — a consistent downward shift in the `maxeffCons` MR estimate.
- The IJ standard error is systematically **larger**: mean 0.386 vs 0.363, mean
  difference +0.024, no replicate agreeing exactly.
- `ij_source` is `"ij"` in both, so this is not a fallback path being taken.

This is the expected signature of a genuinely different $\mathcal{S}$. A rule that
ranks on effect subject to a consistency screen mounts a different competition on every
draw than one ranking on consistency; more optimism subtracted gives a lower corrected
estimate, and a different spread of re-selected winners gives a different IJ variance.
Downstream: MR coverage of $\beta(\widehat H)$ rises to 0.986 (vs 0.968 for legacy and
`maxcons`), and MR harm flags fall to 109 (vs 145 legacy, 134 `maxcons`).

## The result that does need explaining

By the same logic, `maxcons` reproducing legacy MR to **machine precision** on 283
replicates says its selection map is *the same map* as legacy `"eff"` — not a
different rule that happens to agree often. Had $\mathcal{S}$ differed at all, the
per-draw re-selection would have differed and the agreement could not be exact.

That constrains the interpretation of the 18.2% of replicates where `maxcons` and
legacy select different subgroups (documented in
`paired_maxcons_vs_t1t2_n500.md`): with the ranking rule held identical, those
differences must originate in the **candidate family**, not in how candidates are
ranked. The explicit `stop_threshold = NULL` altering where enumeration terminates
remains the natural reading — a different family yields a different argmax on the
affected replicates, while leaving MR's arithmetic identical wherever the family still
delivered the same winner. It also fits the observed pattern there: changes skewed
toward larger regions (44 larger vs 7 smaller) and sparing 82% of replicates, which a
ranking change would not produce.

`mr_ok` is 346 / 346 in both new arms — MR returned an estimate on every detection.

# 4. Estimation and coverage, harm block Ĥ

| arm | estimator | avg | median length | bias vs β(Ĥ) | coverage of β(Ĥ) |
|---|---|---|---|---|---|
| legacy | oracle | 1.135 | 1.441 | +54.9% | 0.835 |
| | naive | 1.629 | 2.082 | +128.0% | 0.145 |
| | FB | 1.222 | 1.659 | +71.0% | 0.801 |
| | MR | 1.003 | 1.473 | +40.0% | 0.968 |
| `maxeffCons` | oracle | 1.135 | 1.441 | +55.2% | 0.838 |
| | naive | 1.643 | 2.157 | +130.9% | 0.153 |
| | FB | 1.229 | 1.749 | +72.6% | 0.806 |
| | MR | **0.946** | 1.536 | **+32.8%** | **0.986** |
| `maxcons` | oracle | 1.135 | 1.441 | +55.6% | 0.827 |
| | naive | 1.611 | 1.964 | +125.9% | 0.130 |
| | FB | 1.223 | 1.648 | +71.3% | 0.783 |
| | MR | 0.992 | 1.395 | +38.7% | 0.968 |

The oracle row is identical across all three (1.135, length 1.441) as it must be — it
refits on the true region and does not depend on the focus rule. Note the oracle's
bias and coverage columns here are computed against β(Ĥ), which is a region mismatch;
they are reported for completeness but are not a bias or a coverage rate for the
oracle. FB under-covers at 0.78–0.81 in every arm; naive at 0.13–0.15.

Complement block Ĥᶜ is stable across arms: MR coverage 1.000 everywhere, MR bias
+0.9% / −0.7% / −0.2%, FB −4.2% / −3.6% / −5.0%.

# 5. Timings

All three arms record per-replicate timings (legacy `t1_secs`/`t2_secs`, new arms
`fb_secs`/`fit_mr_secs`); the rename preserves the measured quantity.

| arm | FB bootstrap, mean s/rep | fit+MR, mean s/rep |
|---|---|---|
| legacy | 166.4 | 12.78 |
| `maxeffCons` | **136.2** | **5.01** |
| `maxcons` | 167.5 | 12.66 |

**These means are not comparable across arms** — per-replicate cost is dominated by
batch render conditions, not by the code. Batch means:

| arm | FB by batch | fit+MR by batch |
|---|---|---|
| legacy | 113.3 / **372.3** / 179.8 / 114.5 | 15.70 / 9.11 / **5.02** / 16.15 |
| `maxeffCons` | 132.8 / 131.7 / 142.9 | 4.85 / 4.92 / 5.16 |
| `maxcons` | 187.0 / 187.8 / **136.5** | 16.57 / 17.85 / **5.02** |

Legacy varies 3.3-fold across its own batches. FB seconds correlate with |Ĥ| at
r = −0.11 (`maxeffCons`) and 0.03 (`maxcons`) — essentially nothing.

Two things can nonetheless be said:

1. **`maxeffCons` was rendered under uniform conditions.** Its batch means vary by
   ~8% (FB) and ~6% (fit+MR), against 3.3-fold and 3.2-fold for legacy. Its numbers
   are the most internally trustworthy of the three.
2. **The fastest batch in each arm converges on the same fit+MR cost**: 5.02 (legacy,
   batch 101–200), 4.85 (`maxeffCons`, batch 1–50), 5.02 (`maxcons`, batch 301–500).
   That agreement across three independently rendered arms is the best available
   estimate of the true uncontended cost — roughly 5 s per replicate — and it is
   **unchanged across focus rules**. The 12–18 s figures are contention, not work.

For FB the equivalent floor is less clean: 113.3 / 131.7 / 136.5. Some of the gap may
be real, but with only one usable batch per arm it cannot be separated from residual
load.

If a timing figure is needed for reporting, use `maxeffCons` — it is the only arm
without a visible batch confound — and quote fit+MR at ~5 s/rep and FB at ~136 s/rep
under 115 workers, noting these are wall-clock under whatever contention obtained.

# 6. Summary

- **Detection**: identical in all three arms, 346/500, zero flips.
- **Selection**: `maxeffCons` and `maxcons` are distinct identifiers (71.1% mutual
  agreement) that bracket legacy — `maxeffCons` tighter (mean |Ĥ| 78.1),
  `maxcons` more inclusive (87.5), legacy between (82.3).
- **MR at coincident selection**: `maxeffCons` differs from legacy systematically
  (−0.067 in the point estimate, +0.024 in the IJ SE, every replicate affected), as a
  different selection map $\mathcal{S}$ implies — $\mathcal{S}$ enters the correction
  through the per-draw re-selected winner. `maxcons` instead reproduces legacy MR to
  machine precision, which identifies its ranking rule as *the same map* as `"eff"`.
- **Coverage**: MR 0.968 (legacy and `maxcons`) vs 0.986 (`maxeffCons`); FB 0.78–0.81
  and naive 0.13–0.15 throughout. Headline JASA conclusions hold in all three arms.
- **Timings**: batch-confounded across arms; uncontended fit+MR cost is ~5 s/rep in
  all three.

`maxeffCons` is a genuinely different identifier, and its MR results differ
accordingly — that is the method working as specified, not a discrepancy to trace.
Reporting it alongside legacy means reporting two identifiers, so the arms should be
labelled as such rather than treated as versions of one procedure.

The question that remains open concerns `maxcons`. Its exact MR agreement establishes
that its ranking rule is unchanged from legacy, which means the 18.2% of replicates
selecting differently must be explained by the candidate family rather than the rule —
most plausibly the explicit `stop_threshold`. That is worth resolving before `maxcons`
numbers are quoted against the JASA tables, since it determines whether the two are
the same procedure on a changed enumeration or something else.
