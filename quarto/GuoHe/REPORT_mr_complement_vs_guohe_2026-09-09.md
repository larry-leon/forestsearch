# REPORT — T1: complement and joint constructions on t7 (Phase A, Mac Studio), 2026-09-09

Governing documents: v1 (`179ae409`) §4, v2 (`669c9ef5`) A2/A3/A4, v4 (`fb65b4a3`) N1–N5, and
Larry's standing authorization of 2026-09-09. Driver:
`quarto/GuoHe/mr_field_complement_vs_guohe_run.R`. Bundles: `mr_field_complement_vs_guohe_t7_beta2_0{0..5}.rds`
(`0fb54669`).

**RESULT: all six cells complete and clean. 2000/2000 replicates per cell, 0 errors,
0 selection mismatches in 12,000, no float column outside `all.equal` 1e-8.**

---

## OPEN ITEMS

Recorded and carried forward, per the standing authorization. None blocks a number in this
report.

1. **`REVIEW_certification_2026-09-09.md` is absent** from the tree (verified again at this
   commit). The B6 certified figures render as the literal marker
   `[certification citation pending sync]`. Fills in when the record syncs (A2).
2. **`theta` and `gamma_s_naive` were reclassified** from the discrete class to the float class
   (`171d2ecc`, approved by Larry). Rationale and measurements in §3 below. The v4 N1 text lists
   "truth lookups keyed by the selection" under the `identical()` class; the reclassification
   narrows that to the *key*, not the returned value.
3. **The N3 rationale quoted a selection-gap minimum of 1.438e-04 from 106 selections. Over the
   full 12,000 production selections the true minimum is 8.463e-08** — three orders smaller.
   The margin against floating point remains ~1.4 × 10⁷ and the observed mismatch count is
   0/12,000, so the conclusion is unchanged, but the 106-selection figure understated the tail
   and should not be requoted. See §2.
4. **Two legacy counters, `naive_mm` and `cur_mm`, are non-zero by construction** and are not
   gates under v4. Both are pre-v4 plain-`identical()` comparisons that include float columns.
   `cur_mm` is 1993–1996 per cell (the cross-platform comparison N1 demotes to a provenance
   measurement); `naive_mm` is 170–262 per cell and is confined to `naive_point`/`naive_lower`
   at the 1e-17 scale. Retained in the bundles for continuity with the committed 16-cell campaign.
5. **The engine returns the complement's 0.95 upper bound as `upper_1s` / `upper_1s_s`**, not
   under a name matching v1's prose. The N5 quotes in §6 pin the actual field names; the driver
   reads them, and nothing is reconstructed.
6. **No Quarto render is attempted in this report.** The T3 render is N10 item 5 and is recorded
   separately.
7. **Line numbers for `R/fs_mr_inference.R` were refreshed on 2026-09-09** against merge
   `f221f75e`, which expanded the `ci_method` roxygen block and shifted everything from line
   298 onward by +8. §6's quotes were re-verified **by content**, not by offset, and are
   unchanged; quoted content unchanged; no result affected. The same merge flipped the
   `ci_method` default from `"ij"` to `"field"`; the T1 driver passes `ci_method` explicitly,
   so the six bundles and every number in this report are unaffected.

---

## 1. Provenance

```
git log -1 --oneline  ->  0fb54669 T1 production: six complement/joint bundles on t7, all six cells clean under the N1/N3 gate
```

Machine: Mac Studio, 14 physical cores, 36 GiB; arm64, R 4.5.2, Accelerate. Run with
`VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1`, 10 forked workers. Stored comparison bundles were
built on `x86_64-pc-linux-gnu`, R 4.6.1, reference BLAS/LAPACK 3.12.0. **No `git fetch`,
`git pull` or `git push` at any point.**

Gate configuration recorded in every bundle (`pair_standard`): rule `v4 N1`, discrete
`identical()`, float `all.equal` at `1e-8`, with both column lists stored. All six bundles carry
the identical definition — the run used `--force` for exactly this reason.

## 2. N3 — selection tally

**Zero selection mismatches in 12,000 replicates.**

| cell | `sel_mismatch` | gap minimum | gap 1st pctile | gap median |
|---|---|---|---|---|
| `t7_beta2_00` | **0 / 2000** | 8.463e-08 | 6.858e-05 | 5.250e-03 |
| `t7_beta2_01` | **0 / 2000** | 2.478e-06 | 6.442e-05 | 5.730e-03 |
| `t7_beta2_02` | **0 / 2000** | 1.125e-05 | 1.027e-04 | 6.035e-03 |
| `t7_beta2_03` | **0 / 2000** | 4.479e-06 | 1.069e-04 | 6.288e-03 |
| `t7_beta2_04` | **0 / 2000** | 9.546e-06 | 1.218e-04 | 6.805e-03 |
| `t7_beta2_05` | **0 / 2000** | 1.257e-05 | 1.305e-04 | 7.364e-03 |
| **total** | **0 / 12000** | **8.463e-08** | — | — |

The gap is top-1 minus top-2 oriented score at selection; `orient = +1` on t7, so the oriented
score is the per-candidate Cox estimate and the selection is `which.max`.

**Honest correction (OPEN ITEM 3).** The N3 rationale carried a minimum of 1.438e-04 measured
over 106 selections. The full 12,000-selection sample finds **8.463e-08**, three orders smaller —
the small sample simply had not reached the tail. Against the worst absolute float deviation
observed anywhere in this run (6.22e-15) that is still a ratio of **1.4 × 10⁷**, and the realized
mismatch count is zero, so the conclusion stands: floating-point deviation at this scale cannot
reach a selection boundary. The 1.438e-04 figure should not be requoted.

## 3. N1 — deviation summary against the stored bundles

**No discrete mismatch and no `all.equal` failure in any cell.** 62,000 float values compared per
cell (372,000 in total).

| cell | `disc_mm` | `ae_fail` | worst absolute | worst relative | bit-identical |
|---|---|---|---|---|---|
| `t7_beta2_00` | **0** | **0** | 2.69e-15 (`fld_upper_2s`, m=298) | 1.03e-11 (`fld_q50`, m=622) | 46.8% |
| `t7_beta2_01` | **0** | **0** | 2.55e-15 (`mr_lower_2s`, m=1896) | 1.19e-11 (`fld_q50`, m=1306) | 46.2% |
| `t7_beta2_02` | **0** | **0** | 1.71e-15 (`fld_lower_se`, m=1419) | 8.97e-11 (`mr_bias_fix`, m=1720) | 45.9% |
| `t7_beta2_03` | **0** | **0** | 4.22e-15 (`mr_upper_2s`, m=724) | 6.87e-12 (`fld_lower_2s`, m=582) | 45.8% |
| `t7_beta2_04` | **0** | **0** | 2.66e-15 (`mr_lower_2s`, m=839) | 7.00e-12 (`mr_bias_fix`, m=1530) | 45.6% |
| `t7_beta2_05` | **0** | **0** | 6.22e-15 (`mr_upper_2s`, m=426) | 6.67e-12 (`fld_lower_2s`, m=1849) | 44.6% |

Worst deviation anywhere: **6.22e-15 absolute, 8.97e-11 relative** — respectively seven and three
orders inside the 1e-8 tolerance.

### The `theta` / `gamma_s_naive` reclassification (OPEN ITEM 2)

**N1's discrete class exists to prove the *selection* is unchanged, not to constrain computed
continuous values.**

- The selection keys — `c_hat`, `c_hat_naive`, `sel`, `n_sel`, `seed_data`, `seed_mr` — are
  verified `identical()` at **24,000/24,000**, and the production run adds `sel_mismatch` 0/12,000
  on top of that.
- The truth lookups are **not labels**. `gh52_truth_at()` is, at
  `quarto/GuoHe/guohe_sec52_truth.R:314`:

```r
  stats::approx(truth$c_grid, y, xout = c_hat, rule = 2)$y
```

  i.e. linear interpolation. Its return is computed arithmetic and belongs to the float class by
  nature, whatever the key.
- Measured residual, obtained by feeding the stored `c_hat` keys straight to the lookup with no
  Cox refits: **24 of 24,000 values differ (0.10%), worst |diff| 5.55e-17**, `all.equal` at 1e-8
  TRUE throughout.

| cell | θ differing | γ_s,naive differing | worst \|diff\| |
|---|---|---|---|
| `t7_beta2_00` | 0/2000 | 0/2000 | 0 |
| `t7_beta2_01` | 3/2000 | 3/2000 | 1.39e-17 |
| `t7_beta2_02` | 2/2000 | 2/2000 | 2.78e-17 |
| `t7_beta2_03` | 1/2000 | 1/2000 | 2.78e-17 |
| `t7_beta2_04` | 3/2000 | 3/2000 | 5.55e-17 |
| `t7_beta2_05` | 3/2000 | 3/2000 | 5.55e-17 |

At β₂ = 0 the truth curve is identically zero, so interpolation returns exact 0 and that cell
shows 0/2000 — which is why the six-replicate Stage 1 probe never surfaced this, and why cell 00
completed clean before the first halt.

### Complement exoneration (N1) — established, not re-run

The Mac isolation test — complement enabled vs disabled, **same machine, same seeds, 29/29 shared
columns `identical()`, both probe cells, all three replicates** — is the proof that enabling the
complement perturbs nothing. It is a stronger proof than the cross-machine comparison A3
specified, and it is cited here with its commit: **`7b2ed976`** (appended to
`REPORT_guohe_supp_stage0_2026-09-09_v2.md`). It matches the engine's own design claim at
`R/fs_mr_inference.R:828-834` and `:897-918`.

## 4. Per-cell results

Complement truth on this design is **exactly 0**: `guohe_sec52_truth.R:75` sets
`b <- ifelse(w <= GH52_C_LO, beta2, 0)` with `GH52_C_LO = 30`, and the cutpoint grid starts at 30,
so every Ĥᶜ = {W > ĉ} lies inside the null region. Complement coverage is `upper >= 0`, with no
dilution term. All values on the oriented log-hazard-ratio scale. Wilson 95% intervals; MCSE at
coverage 0.95 and 2000 replicates is ≈ 0.0049.

### 4a. Complement one-sided 95% upper bound against truth 0

| cell | field-s coverage (Wilson) | field-s mean location | field coverage (Wilson) | field mean location |
|---|---|---|---|---|
| `t7_beta2_00` | 0.928 (0.916, 0.939) | 0.2922 | 0.926 (0.914, 0.937) | 0.2868 |
| `t7_beta2_01` | 0.940 (0.928, 0.949) | 0.2915 | 0.935 (0.923, 0.945) | 0.2897 |
| `t7_beta2_02` | 0.938 (0.927, 0.948) | 0.2823 | 0.939 (0.927, 0.948) | 0.2859 |
| `t7_beta2_03` | 0.941 (0.930, 0.951) | 0.2713 | 0.940 (0.928, 0.949) | 0.2775 |
| `t7_beta2_04` | 0.948 (0.937, 0.957) | 0.2827 | 0.950 (0.940, 0.959) | 0.2892 |
| `t7_beta2_05` | 0.936 (0.925, 0.946) | 0.2753 | 0.939 (0.927, 0.948) | 0.2828 |

Field-s spans **0.928–0.948**, with the single low point at β₂ = 0 (the null cell, where the
selection is loosest). Field-s and field are within 0.005 of each other in every cell and their
Wilson intervals overlap throughout: on this design the studentization neither helps nor hurts
materially, unlike the effMaxSG band regimes where it was worth 0.7–1.8 points.

### 4b. Joint two-subgroup claim (Bonferroni pair, 0.975 each side)

| cell | `joint_s` coverage (Wilson) | `joint` coverage (Wilson) | mean γ | corr(Λ*, Λ*ᶜ) |
|---|---|---|---|---|
| `t7_beta2_00` | 0.933 (0.921, 0.943) | 0.931 (0.919, 0.941) | 0.0252 | +0.018 |
| `t7_beta2_01` | 0.933 (0.921, 0.943) | 0.932 (0.920, 0.942) | 0.0252 | +0.021 |
| `t7_beta2_02` | 0.942 (0.931, 0.951) | 0.941 (0.930, 0.951) | 0.0252 | +0.021 |
| `t7_beta2_03` | 0.946 (0.935, 0.955) | 0.944 (0.933, 0.953) | 0.0252 | +0.018 |
| `t7_beta2_04` | 0.946 (0.935, 0.955) | 0.947 (0.936, 0.956) | 0.0252 | +0.015 |
| `t7_beta2_05` | 0.942 (0.930, 0.951) | 0.943 (0.931, 0.952) | 0.0252 | +0.009 |

Joint coverage spans **0.931–0.947** against a 0.95 target, rising with β₂. The correlation
between the harm and complement fields is +0.009 to +0.021 — small and positive, so the
calibrated pair coincides with Bonferroni and calibration recovers nothing, consistent with every
earlier campaign. 1000 aligned outer draws per replicate.

### 4c. Harm side, and the cross-check against the stored field column

| cell | harm 0.95 (this run) | stored `fld_cover_1s` | vector `identical()` | harm 0.975 |
|---|---|---|---|---|
| `t7_beta2_00` | 0.940 (0.928, 0.949) | 0.940 (0.928, 0.949) | **TRUE** | 0.971 (0.962, 0.977) |
| `t7_beta2_01` | 0.933 (0.922, 0.944) | 0.933 (0.922, 0.944) | **TRUE** | 0.967 (0.958, 0.974) |
| `t7_beta2_02` | 0.938 (0.926, 0.947) | 0.938 (0.926, 0.947) | **TRUE** | 0.970 (0.961, 0.976) |
| `t7_beta2_03` | 0.937 (0.925, 0.947) | 0.937 (0.925, 0.947) | **TRUE** | 0.970 (0.961, 0.976) |
| `t7_beta2_04` | 0.941 (0.929, 0.950) | 0.941 (0.929, 0.950) | **TRUE** | 0.971 (0.962, 0.977) |
| `t7_beta2_05` | 0.951 (0.941, 0.960) | 0.951 (0.941, 0.960) | **TRUE** | 0.976 (0.968, 0.982) |

**The cross-check is exact, not merely aggregate**: the recomputed 0.95 coverage *indicator
vector* is `identical()` to the stored `fld_cover_1s` in all six cells, 2000 replicates each. The
harm-side product is reproduced bit-for-bit as a decision, on a platform where the underlying
continuous values differ in the last bits — the sharpest available demonstration that the
platform difference does not reach any reported conclusion.

Harm coverage spans 0.933–0.951 at 0.95, and 0.967–0.976 at the 0.975 Bonferroni level.

### 4d. p̂(Ĥ), and marginal SD beside error SD

| cell | p̂(Ĥ) median (IQR) | p̂(Ĥ) mean | harm marginal SD | harm error SD | complement marginal SD | complement error SD |
|---|---|---|---|---|---|---|
| `t7_beta2_00` | 0.196 (0.135, 0.276) | 0.221 | 0.1856 | 0.1856 | 0.2014 | 0.2014 |
| `t7_beta2_01` | 0.209 (0.142, 0.299) | 0.234 | 0.1959 | 0.1905 | 0.1907 | 0.1907 |
| `t7_beta2_02` | 0.229 (0.159, 0.324) | 0.256 | 0.1951 | 0.1879 | 0.1828 | 0.1828 |
| `t7_beta2_03` | 0.244 (0.166, 0.347) | 0.273 | 0.2025 | 0.1930 | 0.1802 | 0.1802 |
| `t7_beta2_04` | 0.272 (0.183, 0.393) | 0.301 | 0.2090 | 0.1969 | 0.1776 | 0.1776 |
| `t7_beta2_05` | 0.309 (0.204, 0.442) | 0.335 | 0.2119 | 0.2021 | 0.1802 | 0.1802 |

p̂(Ĥ) — the re-selection probability of the observed winner — rises monotonically with β₂
(median 0.196 → 0.309): a stronger signal makes the pick more stable, as expected. It stays low
in absolute terms throughout, which is the nested family doing what it does — ~151 candidates
differing by one subject apiece share the mass.

On the harm side the marginal SD exceeds the error SD in every cell with β₂ > 0 (e.g. 0.2119 vs
0.2021 at β₂ = 0.5), because γ_ĉ itself moves with the selection; at β₂ = 0 the truth is flat at
zero and the two coincide exactly (0.1856). On the complement side the two are **identical by construction, not by
measurement**: the complement's truth is the constant 0, so the error SD is `sd(x - 0) = sd(x)`.
The column is shown for symmetry with the harm side and carries no independent information. **The complement's marginal SD is the honest scale for its bound, and it is not
inflated relative to the error scale anywhere in this design.**

## 5. Cost against projection

| quantity | Gate 1a projection | realized |
|---|---|---|
| per-replicate (serial) | 1.67 s | ≈ 1.58 s |
| total, 6 × 2000 | 5.6 core-h | ≈ 5.3 core-h |
| wall at 10 workers | 33.4 min | **31.6 min** |
| envelope | ≤ 40 core-h, ≤ 90 min | **WITHIN both** |

Per cell: 5.31, 5.04, 5.10, 5.40, 5.51, 5.24 min. No cell errored; `mr_na`, `fld_na` and `c_na`
are 0 in all six.

## 6. N5 — engine field names for the emitted columns

The probe guards only columns with stored counterparts; the complement and joint columns have
none. These quotes, read from this tree, show the emitted values are **read from the engine, not
reconstructed in the driver**. No `R/` change was made or needed (A2).

**Complement one-sided 95% upper bound** — `R/fs_mr_inference.R:1158-1160`:

```r
    est2 = to_eff(est2_w),
    # Primary: the one-sided 95% UPPER bound (benefit claim "at most U").
    upper_1s = to_eff(bdc - qs[1]),
```

and its studentized companion, `:1171-1172`:

```r
    est2_s = to_eff(est2s_w),
    upper_1s_s = to_eff(bdc - qss[1]), lower_1s_s = to_eff(bdc - qss[5]),
```

**Joint Bonferroni pair and the both-correct indicator** — `.fs_mr_field_joint`, defined at
`R/fs_mr_inference.R:1203`, returning at `:1223-1226`:

```r
       bonf_gamma = alpha / 2,
       bonf_lower_H = to_eff(beta_deb - qh_b), bonf_upper_Hc = to_eff(bdc - qc_b),
       bonf_joint_prob = mean(lh <= qh_b & lc >= qc_b),
       corr = stats::cor(lh, lc), n_joint_draws = n,
```

**Where they attach to the returned object** — `R/fs_mr_inference.R:915-917` and `:1176`:

```r
        field$complement <- fcres$complement
        if (!is.null(fcres$joint)) field$joint <- fcres$joint
        if (!is.null(fcres$joint_s)) field$joint_s <- fcres$joint_s   # field-s (add-beside)
```
```r
  list(complement = complement, joint = joint, joint_s = joint_s)
```

The driver reads `f$complement$upper_1s`, `f$complement$upper_1s_s`, `f$joint$bonf_lower_H`,
`f$joint$bonf_upper_Hc`, `f$joint_s$bonf_upper_Hc` and `f$joint$bonf_joint_prob` from exactly
these, logging each to the oriented scale. The harm-side 0.975 lower bound is taken once: the
harm draws are common to `joint` and `joint_s`, so `bonf_lower_H` is identical in both.

## 7. Reading

- **The complement bound is delivered and honest on this design.** Field-s covers 0.928–0.948
  against a nominal 0.95, its mean location sitting at 0.27–0.29 log-HR above a truth of 0. The
  one cell below 0.93 is the null, where the selection is loosest and ĉ ranges over the whole
  grid.
- **The joint two-subgroup claim holds at 0.931–0.947**, against 0.95, and improves as β₂ grows.
  With corr(Λ*, Λ*ᶜ) ≤ +0.021 the Bonferroni pair is essentially the calibrated one.
- **Nothing here is a comparison with Guo & He.** Their correction is defined through the maximum
  functional over the supplied family; Ĥᶜ is not in that family and is not the argmax of any
  functional, so no analogous bound exists in their framework. This section is a capability
  demonstration, and the supplement says so.
- **The platform question is closed for this deliverable.** 0 selection mismatches in 12,000,
  every float within 1e-8 with the worst at 6.22e-15, and the harm-side coverage *indicator*
  reproduced by `identical()` in all six cells. The arm64/Accelerate build reaches the same
  conclusions as the x86_64/reference-BLAS build on every reported quantity.
