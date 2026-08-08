# STOP A — md/harm MR pilot, re-run after the configuration fix

Deliverable F3 of `dev/glm-continuous-sims/CC_TASK_stopA_fix_and_repilot.md`
(saved verbatim and committed at `49337221`). Baseline tip `2d14b7ae`.

**Transport, not approval.** STOP B is NOT started.

## Header

| | |
|---|---|
| F1 harness fix | `49337221` — five edits, diff vs provenance `169847d5` shows only those plus the saved task doc |
| Re-pilot bundle | `quarto/simulations/actg175/continuous/mr_sweep_md_harm/md_harm_s50_pilot/fs_md_harm_n1000_res.rds`, 9750 B, mtime 2026-08-08 10:57 |
| Package | `forestsearch` 0.2.0, R 4.6.1, installed at `ef504707` |
| Run | 50 replicates, n = 1000, `mr_inference = TRUE`, `mr_draws = 5000`, `consistency_method = "resample"`, 115 workers, PID-watched (PID 473298, `kill -0`) |

The previous pilot bundle (1707 B, detected 0/50) is preserved in history at
`169847d5`; this run overwrites the working-tree file because the cell tag
`md_harm_s50_pilot` is unchanged, as F2 requires.

---

## 1. Detection and MR completion

```
status:      DETECTED  50
detected:    1  -> 50
mr_ok:       1  -> 50
err_msg non-NA : 0
mr_msg  non-NA : 0
betaHhat_status: ok 50
n_true[1] = 341   (pilot bundle: 341 — RNG pin reproduces it)
```

**detected 50/50. mr_ok 50/50.** `mr_ok == detected` on every replicate, so
there are no `err_msg` or MR-skip messages to quote: the F1(c) diagnostics
fired zero times because nothing failed. Against the pre-fix pilot's 0/50 on
the identical seed table, both stated fixes are confirmed effective — the guard
(addendum §D.1) and `consistency_method` (§D.4).

The RNG pin (F1(d)) reproduces the pilot's replicate 1 exactly: `n_true = 341`.

Partition holds on all 50 rows: `nH_eval + nHc_eval = 5000 = N`.

## 2. Per-cell readouts (corrected q-column extraction)

### 2.1 Scale convention — read this before the numbers

The bundle carries **two different scales**, and conflating them is what
CHECK C caught (§5.2):

- `nv_H_est`, `t2_H_est`, `t2_H_lo/hi` are on the **oriented** scale, where
  positive means harm (`adverse_outcome = FALSE` negates Y inside the
  estimator). Range of `nv_H_est`: [51.75, 133.70].
- `betaHhat_H`, `md_trial_H` are on the **raw** `cd4_change` scale, where
  negative means harm. Range of `betaHhat_H`: [-40.00, -26.26].

Every comparison below is made on the **oriented** scale, negating the raw
targets (`bH = -betaHhat_H`). This is stated rather than assumed because the
harness's own readout 1 does not do it — see §5.2.

### 2.2 Bias against the EXACT β(Ĥ)

```
  naive  mean bias :    67.0722   SE  2.3024
  MR     mean bias :    35.4726   SE  2.7009
  reduction        :    31.5996   (47.1% of naive bias removed)

  naive bias  min/q1/med/q3/max : 15.9371 / 57.2312 / 70.1543 / 77.0280 / 102.0428
  MR    bias  min/q1/med/q3/max : -18.9500 / 21.6886 / 41.8487 / 49.4364 /  75.1083
```

MR removes 47.1% of the naive selection bias. **A residual mean bias of
+35.47 (SE 2.70) remains — 13 SEs from zero, so this is not noise.** MR is a
leading-order approximation of the full bootstrap correction
(`R/fs_mr_inference.R:1-29` says so explicitly), and this measures how much of
the optimism the leading order captures at n = 1000: about half.

### 2.3 Coverage of β(Ĥ)

```
  H  coverage : 49/50 = 0.980   n_eff = 50   MC SE 0.0198
  Hc coverage : 50/50 = 1.000
  mean interval width (H) : 121.8052
  mean se_ij              :  31.0733
  ij_source : ij  (50/50)
```

**Coverage is 0.980 against a nominal 0.95**, within 1.5 MC SE. The one
non-covering replicate is id 24:

```
 sim_id                    sg_def  nv_H_est  t2_H_est   t2_H_lo   t2_H_hi  betaHhat_H
     24 {drugs} & !{preanti <= 0}   133.704  106.7696  33.75843  179.7808   -31.66127
```

oriented β(Ĥ) = 31.66, just below the lower limit 33.76.

**Coverage is achieved with wide intervals, not with unbiasedness.** Mean width
121.8 against a residual bias of 35.5: the interval is wide enough to cover
despite being centred ~35 units high. That is the honest reading, and it is the
number STOP B should be judged on at larger n.

### 2.4 Selection

```
  selection_rate : 50/50 = 1.000
  selection effect (naive - exact) : mean 67.0722   SD 16.2802
    min/q1/med/q3/max : 15.9371 / 57.2312 / 70.1543 / 77.0280 / 102.0428
  selected size n_harm min/median/max : 61 / 86 / 354
```

Selection rate 1.0 at this cell is the recorded phenomenon, not a defect
(closed decision, `CC_TASK_md_mr_harness.md`).

### 2.5 Identification vs Q

```
  sens mean 0.1353   min/med/max 0.0000 / 0.0967 / 0.7302
  spec mean 0.9168   min/med/max 0.7671 / 0.9206 / 1.0000
  ppv  mean 0.4196   min/med/max 0.0000 / 0.3984 / 1.0000
  npv  mean 0.6709   min/med/max 0.6061 / 0.6637 / 0.8684
```

**Sensitivity averages 0.135** — the selected region typically captures ~14% of
Q. With median selected size 86 against Q's ~340, that is arithmetic: the
identifier picks small high-effect regions. `sens = 0` occurs on at least one
replicate (a selected region disjoint from Q). NPV is reported per project
convention, not Jaccard.

### 2.6 θ† reference columns — a DISTINCT estimand

```
  theta_dagger_H  = -40.0000000000   (== DGM effect_Q  -40.0000000000)
  theta_dagger_Hc = -26.2552358760   (== DGM effect_Qc -26.2552358760)

  mean exact beta(Hhat) : raw -32.1351   oriented 32.1351
```

θ† is evaluated at the **true flag Q**; β(Ĥ) at the **realized rule Ĥ**. They
are different estimands and are never mixed above. θ† reproduces the DGM's own
effects exactly (that is CHECK B, §5.1).

## 3. Naive vs MR beside the O1 anchors

| quantity | O1 (same seeds) | F2 |
|---|---|---|
| naive mean, in-trial oriented MD of the selected region | 100.8259 | **99.2074** |
| excess over the true region Q | 59.4244 | see note |

The two naive means are close but **not identical, and should not be**: O1 ran
`consistency_method` at its formal default `"split"` with `mr_inference = FALSE`,
F2 runs `"resample"` with MR on. Different consistency machinery selects
different rules, so the realized regions differ (O1 median selected size 73.5;
F2 86). The agreement to ~1.6 units is a consistency check across two different
selection engines on the same seeds, not a reproduction.

The recorder's own in-trial figure agrees with the gate's to machine precision
once put on the same scale: `mean(-md_trial_H) = 99.2074 = mean(nv_H_est)`.

## 4. Per-replicate table (all 50)

```
id st nT nH naive MR lo hi bHhat cov sens ppv secs
1 DETECTED 341 85  105.758   72.839   11.984  133.695  30.599 Y 0.070 0.282 17.7
2 DETECTED 322 72  113.668   83.839   13.750  153.928  30.760 Y 0.068 0.306 19.8
3 DETECTED 335 116   75.546   30.548  -21.315   82.410  31.267 Y 0.116 0.336 17.3
4 DETECTED 352 96  103.087   76.887   20.383  133.392  32.946 Y 0.116 0.427 18.4
5 DETECTED 331 138  108.958   79.470   18.478  140.461  33.261 Y 0.227 0.543 19.1
6 DETECTED 355 103  109.017   84.053   27.485  140.622  31.542 Y 0.110 0.379 19.4
7 DETECTED 347 67  110.134   80.546    2.143  158.949  29.601 Y 0.049 0.254 20.5
8 DETECTED 324 84   98.184   67.666   -1.842  137.174  31.542 Y 0.080 0.310 19.3
9 DETECTED 339 65  120.501   90.274   11.049  169.499  33.487 Y 0.112 0.585 21.0
10 DETECTED 342 81   98.816   78.339   25.116  131.562  26.819 Y 0.020 0.086 21.7
11 DETECTED 379 66  120.374   86.673   18.153  155.194  33.487 Y 0.095 0.545 22.4
12 DETECTED 349 87   89.403   61.703    3.741  119.664  30.079 Y 0.077 0.310 21.8
13 DETECTED 382 200   77.575   39.290   -4.970   83.550  35.365 Y 0.356 0.680 16.2
14 DETECTED 306 111   84.719   46.117   -5.096   97.331  35.258 Y 0.229 0.631 18.4
15 DETECTED 350 85   86.896   46.328  -12.714  105.370  32.976 Y 0.117 0.482 21.0
16 DETECTED 360 61  117.669   86.680    3.824  169.536  35.793 Y 0.100 0.590 19.3
17 DETECTED 358 86  105.551   77.330    4.728  149.933  35.383 Y 0.154 0.640 19.2
18 DETECTED 359 67  111.551   75.733    6.062  145.404  26.255 Y 0.000 0.000 19.9
19 DETECTED 304 63  105.738   80.574   27.048  134.101  31.261 Y 0.092 0.444 22.8
20 DETECTED 336 107   84.542   47.473  -12.518  107.465  32.208 Y 0.125 0.393 20.4
21 DETECTED 341 93  114.164   85.422   31.222  139.622  36.943 Y 0.211 0.774 20.6
22 DETECTED 361 147   88.362   63.370   22.281  104.458  31.567 Y 0.158 0.388 21.8
23 DETECTED 356 185   68.256   35.317   -5.988   76.622  28.516 Y 0.098 0.189 20.0
24 DETECTED 346 76  133.704  106.770   33.758  179.781  31.661 N 0.075 0.342 20.0
25 DETECTED 315 354   51.750   16.863  -15.752   49.477  35.813 Y 0.730 0.650 21.1
26 DETECTED 339 72  104.828   72.213    7.186  137.241  28.379 Y 0.024 0.111 20.3
27 DETECTED 327 88   88.927   46.364  -17.543  110.272  30.388 Y 0.076 0.284 17.0
28 DETECTED 332 141   75.273   38.340  -12.239   88.918  31.721 Y 0.172 0.404 16.8
29 DETECTED 340 63   98.714   62.303   -7.991  132.597  27.058 Y 0.018 0.095 19.1
30 DETECTED 341 103  100.387   72.139   18.287  125.991  32.894 Y 0.164 0.544 21.9
31 DETECTED 355 75   98.298   71.500    3.324  139.676  29.750 Y 0.062 0.293 20.0
32 DETECTED 322 61  103.995   78.026    0.688  155.364  33.348 Y 0.065 0.344 19.8
33 DETECTED 346 106   89.147   57.771    1.858  113.684  30.168 Y 0.064 0.208 21.7
34 DETECTED 327 72  110.092   84.630   24.103  145.157  31.424 Y 0.064 0.292 18.7
35 DETECTED 353 122   83.695   40.126  -11.436   91.688  36.484 Y 0.275 0.795 17.5
36 DETECTED 346 70   87.689   45.732  -18.726  110.189  26.255 Y 0.000 0.000 20.5
37 DETECTED 350 106  103.797   71.002   15.049  126.954  36.641 Y 0.231 0.764 20.4
38 DETECTED 368 94  111.721   83.039   20.983  145.095  36.054 Y 0.190 0.745 19.5
39 DETECTED 339 123   77.798   45.475   -7.020   97.971  32.786 Y 0.159 0.439 20.5
40 DETECTED 342 71  101.621   71.157   11.432  130.882  36.292 Y 0.155 0.746 18.2
41 DETECTED 342 106  103.673   73.455   17.296  129.615  29.211 Y 0.050 0.160 21.6
42 DETECTED 351 70  120.047   89.842   22.918  156.766  30.056 Y 0.085 0.429 22.4
43 DETECTED 359 68  102.382   70.547    4.238  136.857  32.833 Y 0.089 0.471 21.8
44 DETECTED 337 107  100.985   64.846   -2.967  132.658  36.076 Y 0.220 0.692 20.8
45 DETECTED 346 62  111.374   87.323   23.844  150.801  33.148 Y 0.066 0.371 22.0
46 DETECTED 357 63  120.772   82.500    4.009  160.990  27.417 Y 0.020 0.111 19.5
47 DETECTED 359 71  106.628   77.713    8.088  147.338  32.364 Y 0.086 0.437 20.0
48 DETECTED 349 222   83.076   64.584   23.343  105.826  40.000 Y 0.636 1.000 20.0
49 DETECTED 317 115  101.050   73.564   14.420  132.708  26.255 Y 0.000 0.000 20.7
50 DETECTED 358 109   90.476   56.089   -8.904  121.082  35.366 Y 0.207 0.679 20.6
```

---

## 5. Checks A / B / C, with the new guards

### 5.1 CHECK B — PASSED

Ran before the replicates, as before. θ† reproduces the DGM's own effects under
`identical()`:

```
  theta_dagger_H  = -40.0000000000   effect_Q  = -40.0000000000
  theta_dagger_Hc = -26.2552358760   effect_Qc = -26.2552358760
```

### 5.2 CHECK A — PASSED; the F1(e) guard was not needed

The render proceeded past `check-A` to `check-C`, which it could only do by
passing (the chunk `stop()`s otherwise). With 50/50 detections there were
**50 distinct realized rules**, so the empty-rule guard added in F1(e) did not
trigger. It remains in place for the zero-detection case that broke the first
render (overnight O2.4). The guard is therefore untested by this run — it is
exercised only when detection is 0, which no longer happens under the fixed
configuration.

### 5.3 CHECK C — FAILED. Verbatim:

```
Error:
! CHECK C FAILED: 1.984e+02 — the two computations disagree.

Quitting from mr_coverage_sweep_md_harm.qmd:637-661 [check-C]
Execution halted
WARN: Error encountered when rendering files
```

Components:

```
  naive mean bias       (gate)     : 131.3424696483
  mean selection effect (recorder) : -67.0722315241
  difference                       : 1.984e+02
  max per-replicate |gate naive - recorder MD| : 2.674e+02
```

**The render halted here. Per the task, no fix was iterated.** The HTML was not
produced; every number in this report is computed directly from the committed
bundle.

**Diagnosis — a scale mismatch in the check, not a machinery disagreement.**
The two quantities are exact negatives of one another:

```
 id  gate_naive     rec_MD  sum  betaHhat
  1    105.7582  -105.7582    0  -30.5994
  2    113.6683  -113.6683    0  -30.7595
  3     75.5463   -75.5463    0  -31.2669
  4    103.0873  -103.0873    0  -32.9464
  5    108.9576  -108.9576    0  -33.2612

  max |gate_naive + recorder_MD| over 50 replicates = 2.700e-13
```

`nv_H_est` comes from the MR gate on the **oriented** scale (positive = harm);
`md_trial_H` is computed in `record_replicate()` directly from
`df[[outcome_name]]`, un-negated, on the **raw** scale. Both are individually
correct. CHECK C subtracts each from the *raw* `betaHhat_H`, so the gate's term
mixes an oriented estimate with a raw target and the recorder's does not. The
198.4 discrepancy is exactly twice the mean naive effect, as the exact-negative
relation requires.

**CHECK C did its job**: it caught a genuine inconsistency in the harness's own
readout 1, which reports `nv_H_est - betaHhat_H` — an oriented-minus-raw
quantity — as "naive mean bias". That readout is wrong as written; the
corrected value is §2.2's 67.0722, not the gate's 131.3425.

**Stated, not applied** (the task forbids iterating fixes, and this is a sixth
harness edit beyond F1's five): the minimal correction is to compare on one
scale — either negate `betaHhat_H` where the gate's oriented estimates are
used, or have the recorder compute `md_trial_H` on the oriented scale. The
first is preferable: it leaves the stored columns on their natural scales and
puts the conversion at the point of comparison. This affects readout 1 and
CHECK C only; §2.2's bias and §2.3's coverage above are already computed on a
consistent scale and are unaffected.

---

## 6. Wall clock and the PROPOSED STOP B schedule

### 6.1 Measured

```
  per-replicate t2_secs (includes resample consistency + MR 5000 draws):
    min 16.16   q1 19.22   median 20.01   mean 20.01   q3 21.01   max 22.76
  cell wall clock : 99.308 s for 50 replicates on 115 workers
```

**20.0 s median per replicate including MR**, against ~162 s for the flags-off
search measured in the overnight run (`REPORT_overnight_funnel50.md` O1.0,
sequential). The ~8x reduction is `consistency_method = "resample"` replacing
split-consistency's `fs.splits x 2 = 800` refits per evaluation
(`R/forestsearch_main.R:1433-1436`) with a single fit plus a multiplier
representation. MR's 5000 draws are a small addition on top. **This was
measured, not assumed.**

### 6.2 Proposed schedule — STOP B is NOT started

Grid `n in {500, 1000, 2000, 4000}`, `s = 1000` replicates per cell.

Only n = 1000 is measured (20.0 s/replicate). The projection below assumes cost
scales **linearly in n** — the per-candidate fits are O(n) and the candidate
family size is roughly n-independent (overnight O1.3: 1275-1596 combinations).
**This is an assumption, flagged as such**; the n = 500 and n = 2000 cells will
confirm or refute it within the first few minutes of a run.

| cell | projected s/replicate | core-seconds for s = 1000 |
|---|---|---|
| n = 500 | ~10 | 10,000 |
| n = 1000 | 20.0 (measured) | 20,000 |
| n = 2000 | ~40 | 40,000 |
| n = 4000 | ~80 | 80,000 |
| **total** | | **~150,000 core-seconds** |

At 115 workers that is ~1,300 s of wall clock, plus worker startup and
per-cell bundle I/O — call it **30-45 minutes for the full grid**, comfortably
a single session. No reduction in `s` is warranted by cost; the schedule the
task proposes stands as written.

Recommended order: n = 500 first (cheapest, and it tests the linearity
assumption immediately against n = 1000), then 2000, then 4000.

---

## 7. Against the pre-registered context

Stated in the task as context, not acceptance criteria:

| pre-registered | measured | |
|---|---|---|
| detection ≈ 50/50 | **50/50** | met |
| Q floor-clearance ≈ 36/50 | not re-measured in F2 | see note |
| naive optimism ≈ 2.5x | naive 99.21 vs oriented β(Ĥ) mean 32.14 → **3.09x** | higher |
| MR judged on coverage vs exact β(Ĥ) | **0.980** (nominal 0.95), 49/50 | met |

Q floor-clearance is a property of the search configuration, measured at 36/50
in the overnight run under `"split"`; F2 changed the consistency engine to
`"resample"`, so the overnight figure does not transfer and F2 did not
recompute it. Not claimed either way.

The naive optimism ratio is against β(Ĥ) at the realized rule (3.09x), whereas
the overnight anchor 2.5x was against the true region Q's own in-trial effect —
different denominators, both reported rather than reconciled.

---

## 8. Does STOP A show resample-consistency itself misbehaving?

Relevant to the gated `consistency_method` task, which does not flip the
default if it does.

**No.** On the evidence here:

- 50/50 replicates completed with `consistency_method = "resample"`; zero
  errors, zero `err_msg`, zero MR skips.
- `mr_ok == detected` on every replicate — the resample path is what makes MR
  reachable at all on a GLM outcome (`.mr_glm_ok`, `R/forestsearch_main.R:3142`).
- MR coverage of the exact β(Ĥ) is 0.980 against nominal 0.95.
- Median cost 20.0 s/replicate against ~162 s for split-consistency.

The one check failure (§5.3) is a scale-convention defect in the harness's
CHECK C and readout 1, present before this run and independent of the
consistency engine. The residual MR bias (§2.2) is a property of the
leading-order approximation, not of resample-consistency.

---

## Appendix — U3 certification

ONE `rcmdcheck::rcmdcheck(args = "--as-cran", error_on = "never")` run on the
tree carrying U1 (Rd Unicode fixes) and U2 (CLAUDE.md). Final status line,
verbatim:

```
0 errors ✔ | 0 warnings ✔ | 1 note ✖
```

The single NOTE, verbatim — environmental, as anticipated:

```
❯ checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found.
  Please obtain a recent version of HTML Tidy by downloading a binary
  release or compiling the source code from <https://www.html-tidy.org/>.
```

```
── R CMD check results ───────────────────────────────── forestsearch 0.2.0 ────
Duration: 10m 45.2s
* checking tests ...
  Running ‘testthat.R’ [236s/186s]
 [237s/186s] OK
* checking PDF version of manual ... OK
```

**`checking PDF version of manual ... OK`.** The WARNING recorded in
`REPORT_blockers_1-3.md` §1.2 — the LaTeX build failing on U+1D9C and U+2500 —
is resolved by U1, and the `forestsearch-manual.tex` NOTE that was its
byproduct is gone with it. The certification surface is now clean apart from
the missing `tidy` binary, which is a property of this machine, not the
package.

---

## Where this stops

STOP A ends here. STOP B is **not** started: no grid, no additional cells, no
`s = 1000` run. The §5.3 scale fix is **stated, not applied**. No DGM change.
No `R/` change in Task F. Committing this report is transport, not approval.
