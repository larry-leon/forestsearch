# REPORT — overnight 2: maxeffCons MR grid (O1) and FB/MR batches (O2)

Branch `feature/mr-in-replicates`. Session of 2026-08-09, 01:15–06:25 local.
**Transport, not approval.** Every number verbatim. No `R/` changes, no DGM
changes, no G2.

---

## 0. Headline

**O1's pre-registered question is answered, and the answer is yes.** Every trend
STOP B flagged as going the wrong way under the `hr` focus goes the right way
under `maxeffCons`: residual MR bias falls with n instead of rising, bias/SE
falls instead of rising, and coverage holds near 0.98 instead of collapsing to
0.73.

**O2 adds one hard result the FB adjudication was missing.** The two-source
decomposition shows the entire FB correction is the re-selection term:
`fb_src1` (same-rule resampling noise) averages **0.32**, `fb_src2`
(re-selection) averages **−20.11**. The algebra reconciles at 7.105e-15.

**A caution that survives both.** Pooled over 100 replicates, MR's coverage of
the exact β(Ĥ) is **1.0000** — every interval covers — on a mean width of
155.05. That is not calibration, it is width. FB's coverage is 0.7400 on width
131.28. Neither is at nominal 0.95 from the right direction.

---

## 1. O1 — the `maxeffCons` MR-only grid

4 cells × 1000 replicates, `nb_boots = 0L` (bootstrap never called),
`parallel_mode = "sims"`, `mr_draws = 5000`, `consistency_method = "resample"`,
every other knob as committed in the twin. J = 10 grids deliberately NOT applied
(they are Step 1 / G2 only). Bundles:
`quarto/simulations/actg175/continuous/o1_maxeffCons_mr_grid/o1_maxeffCons_mr_n{500,1000,2000,4000}_s1000.rds`

### 1.1 The trend table, beside the `hr` grid

`maxeffCons` (this run):

```
    n    s detected mr_ok naive_bias   mr_bias      mr_se mr_bias_over_se   coverage        mc_se     width
  500 1000      997   997  71.680622 19.426799 0.67484946       28.786862 0.98996991 0.0031558463 152.27912
 1000 1000     1000  1000  69.759645 18.851831 0.64079691       29.419354 0.97900000 0.0045342033 136.55502
 2000 1000     1000  1000  62.733213 16.233971 0.61893098       26.229048 0.98200000 0.0042042835 125.65309
 4000 1000     1000  1000  51.751929 13.177804 0.56043480       23.513536 0.98600000 0.0037153735 110.51511
```

`hr` (committed, `REPORT_stopB_md_harm_grid.md` § 1):

```
    n    s detected mr_ok naive_bias mr_bias   mr_se mr_bias_over_se coverage      mc_se   width
  500 1000      997   997     69.689  29.808 0.82703          36.043  0.95687  0.0064338 137.186
 1000 1000     1000  1000     65.491  33.202 0.71939          46.153  0.94900  0.0069569 121.518
 2000 1000     1000  1000     58.602  37.164 0.65816          56.467  0.88500  0.0100884 106.063
 4000 1000     1000  1000     48.774  37.388 0.53334          70.102  0.73100  0.0140228  86.193
```

| | `maxeffCons` 500 → 4000 | `hr` 500 → 4000 |
|---|---|---|
| `mr_bias` | 19.43 → **13.18** (falls) | 29.81 → **37.39** (rises) |
| `mr_bias_over_se` | 28.79 → **23.51** (falls) | 36.04 → **70.10** (rises) |
| `coverage` | 0.990 → **0.986** (holds) | 0.957 → **0.731** (collapses) |
| `naive_bias` | 71.68 → 51.75 (falls) | 69.69 → 48.77 (falls) |
| `width` | 152.28 → 110.52 | 137.19 → 86.19 |

### 1.2 Full readout, all columns

```
    n    s detected mr_ok naive_bias   mr_bias      mr_se mr_bias_over_se   coverage        mc_se     width        sens       spec        ppv        npv nsel_q10 nsel_q25 nsel_med nsel_q75 nsel_q90 secs_per_rep_med secs_per_rep_mean cell_elapsed_s
  500 1000      997   997  71.680622 19.426799 0.67484946       28.786862 0.98996991 0.0031558463 152.27912 0.166775660 0.86270833 0.39142055 0.66242653     62.0       64       69       78     92.4          10.3790         10.338561        277.079
 1000 1000     1000  1000  69.759645 18.851831 0.64079691       29.419354 0.97900000 0.0045342033 136.55502 0.097650248 0.92585074 0.40607697 0.66087446     62.0       66       74       90    112.0          40.8855         41.402748        637.086
 2000 1000     1000  1000  62.733213 16.233971 0.61893098       26.229048 0.98200000 0.0042042835 125.65309 0.056296263 0.95664546 0.39384738 0.65788742     64.0       70       83      112    138.1         176.4915        160.343487       1769.901
 4000 1000     1000  1000  51.751929 13.177804 0.56043480       23.513536 0.98600000 0.0037153735 110.51511 0.042293086 0.97016137 0.38572633 0.65800459     67.9       81      121      162    233.0         412.4260        375.543706       3681.577
```

### 1.3 Three qualifications on the headline

**(a) Coverage is above nominal at every n**, 0.979–0.990 against 0.95, on
intervals wider than `hr`'s at every n (152 vs 137, 137 vs 122, 126 vs 106,
111 vs 86). Part of the coverage is bought by width — the mechanism STOP B
named. What differs from `hr` is that the bias falls too, so it is not *only*
width.

**(b) bias/SE = 23.5 at n = 4000 is still enormous.** The residual is shrinking,
not gone. On no reading of these numbers is it noise.

**(c) The identification degrades as the inference improves.** `sens` falls
0.1668 → 0.0977 → 0.0563 → 0.0423 while the median selected size grows
69 → 74 → 83 → 121 and `spec` rises 0.8627 → 0.9702. At larger n the search
selects *larger* regions that overlap the true harm region *less*. β(Ĥ) is
estimated better while Ĥ drifts away from Q; PPV is flat near 0.39 throughout.
This is a separate finding from the coverage result and is not explained by it.

---

## 2. O2 — FB/MR batches

### 2.1 What ran

One wave completed: **batch 51–100**, launched 03:26:44, saved 06:21, elapsed
**10475.55 s = 2.91 h**. `parallel_mode = "sims"`, `inner_parallel$plan =
"sequential"`, `nb_boots = 300`, `mr_draws = 5000`.

**No second wave was launched.** The batch finished at 06:22 and a further wave
costs ~2.9 h, landing ~09:20 — well past the 06:30 stop-launching boundary and
past the morning report. One wave, not the 3–4 the task file projected.

**The cost estimate in the task file was wrong by ~3×.** It projected ~1 h per
50-batch wave under `"sims"` on the basis of ~9,000 core-s per replicate; the
measured wave took 2.91 h. For comparison, batch 1–50 under `"boots"` took
3.35 h. The topology change bought ~13%, not ~70%.

### 2.2 The smoke-reproduction check — NOT RUN AS SPECIFIED

The instruction asks that the first batch reproduce the smoke's 5/5 detected,
5/5 mr_ok, 0 fb_err. **It cannot, as written:** the smoke covered sim ids 1–5
and batch 51–100 covers sim ids 51–100, which are disjoint by construction —
that disjointness is the point of the batching scheme. There is no overlap to
compare on.

What the batch does show, reported instead of the check that was asked for:

```
detected 50/50 | mr_ok 50/50 | fb_err 0
```

A strict mode-invariance test — re-running sim ids 1–5 under `"sims"` and
comparing to their `"boots"` values in batch 1–50 — was **not run**. Under
`"sims"` five replicates do not share the 120 workers usefully (each runs its
300 bootstraps sequentially), so it would cost roughly one replicate's serial
FB time, ~2.5 h, for five rows. That was not a good use of the remaining
window. The invariance argument therefore still rests on the seeding basis
recorded in the `parallel_mode` commit, not on a numerical check.

### 2.3 Per-batch and pooled means

Batch 1–50 (`"boots"`, `inner = multisession`):

```
    nv_H_est     or_H_est     fb_H_est     mr_H_est   betaHhat_H 
104.96985740  41.27360063  85.26969673  53.80907790 -31.79814680 
```

Batch 51–100 (`"sims"`, `inner = sequential`):

```
    nv_H_est     or_H_est     fb_H_est     mr_H_est   betaHhat_H 
102.18003586  37.18885508  82.39043494  50.33576856 -31.94130951 
```

Pooled, batches 1–100, 100 detected replicates:

```
    nv_H_est     or_H_est     fb_H_est     mr_H_est   betaHhat_H 
103.57494663  39.23122785  83.83006583  52.07242323 -31.86972815 
```

`betaHhat_H` is **raw**; the oriented target is **+31.86972815**. The other four
columns are oriented.

The two batches agree closely on every column despite different topologies —
consistent with, though not proof of, the invariance claim.

### 2.4 Bias and coverage against the exact β(Ĥ), pooled

```
--- bias against the ORIENTED exact beta(Hhat) ---
  or_H_est  bias    7.361500   se 2.003503   bias/se    3.674
  nv_H_est  bias   71.705218   se 1.710170   bias/se   41.929
  fb_H_est  bias   51.960338   se 1.827123   bias/se   28.438
  mr_H_est  bias   20.202695   se 2.073133   bias/se    9.745

--- coverage of oriented beta(Hhat) ---
  fb coverage 0.7400 (n_eff 100)   mean width 131.2835
  mr coverage 1.0000 (n_eff 100)   mean width 155.0505

--- corrections ---
  FB (nv-fb) mean 19.744881 | MR (nv-mr) mean 51.502523 | ratio 0.383377
```

The FB/MR correction ratio is **0.383** pooled, against 0.385 measured on batch
1–50 alone — stable across the topology change.

**MR covers 1.0000 of 100.** Every interval covers. With a residual bias of
20.20 and a mean width of 155, that is width absorbing bias, not centring. The
`hr` grid's n = 500 coverage of 0.957 was closer to nominal *from a worse
position*. Coverage at these widths is not evidence of calibration.

Note the **oracle is biased too**: 7.36 with bias/SE 3.67. The oracle refits on
the true region Q, so this is not selection bias — it is the n = 500 refit's own
error against a super-population target, and it sets a floor no estimator on
this cell beats.

### 2.5 The two-source decomposition (Step 2)

Batch 51–100, per-replicate, raw scale:

```
    fb_src1           fb_src2          fb_nres   
 Min.   :-5.7584   Min.   :-24.59   Min.   :299  
 1st Qu.:-0.5044   1st Qu.:-21.30   1st Qu.:300  
 Median : 0.2391   Median :-20.28   Median :300  
 Mean   : 0.3236   Mean   :-20.11   Mean   :300  
 3rd Qu.: 1.4610   3rd Qu.:-18.72   3rd Qu.:300  
 Max.   : 4.6998   Max.   :-13.90   Max.   :300  
```

```
ALGEBRA CHECK max|(nv-fb) + (src1+src2)| = 7.105e-15
```

**The result: the FB correction is entirely re-selection.** `fb_src1` — the
same-rule resampling term, `H_star - H_obs` — averages 0.3236 and straddles
zero. `fb_src2` — the re-selection term, `Hstar_star - Hstar_obs` — averages
−20.11 and never changes sign across 50 replicates (range −24.59 to −13.90).
The whole of the 19.74 correction comes from the re-selection term, as the
theory says it should.

`fb_nres` is 300 on 49 of 50 replicates and 299 on one, so essentially no
resamples were lost to non-finite fits.

**Honesty about what the algebra check does and does not verify.** It confirms
`fb_src1 + fb_src2` reconciles with the observed `nv − fb` at machine
precision. Because both terms were *derived* from `H_biasadj_1`,
`H_biasadj_2` and `H0`, their sum reduces to `H0 − mean(H_biasadj_2)`, so what
the check really establishes is that `get_targetEst()` returns the plain mean
with no trimming or extra correction, and that `H0 = H_obs` exactly. **It does
not independently verify the split between `fb_src1` and `fb_src2`** — that
rests on the source algebra quoted in the Step 2 commit, not on a numerical
test. A genuine check of the split would need the four constituents exposed,
which is an `R/` change and out of scope.

### 2.6 Timing — the 147.65 s outlier gets its answer

`fb_secs`, batch 51–100:

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   3182    5496    6402    6556    7812   10443 
```

`fit_mr_secs`, batch 51–100:

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  7.399   9.258   9.899  10.056  10.788  13.241 
```

Against batch 1–50 (`"boots"`):

| | batch 1–50 | batch 51–100 |
|---|---|---|
| `fb_secs` median | 236.99 | 6402 |
| `fit_mr_secs` median | 15.11 | **9.899** |
| `fit_mr_secs` mean | 25.60 | **10.056** |
| `fit_mr_secs` max | **147.654** | **13.241** |

**The 147.65 s outlier was contention, and it does not recur.** Under `"sims"`
the whole `fit_mr_secs` distribution is tighter *and* faster: max 13.24 against
147.65, mean 10.06 against 25.60, and the 90th-percentile-to-max blow-up is
gone. The earlier reading — one replicate ~10× the median, attributed to
contention on a loaded box rather than a slow fit — is confirmed.

**`fb_secs` is NOT comparable across the two rows** and the table above should
not be read as a 27× slowdown. Under `"boots"` one replicate at a time gets 120
workers for its 300 bootstraps; under `"sims"` 50 replicates run concurrently
and each runs its 300 bootstraps serially, so `fb_secs` measures per-replicate
wall clock under 50-way contention. The comparable quantity is total wave time:
3.35 h (`"boots"`) against 2.91 h (`"sims"`).

---

## 3. Wall-clock accounting

| Item | Elapsed |
|---|---|
| Step 1 anchor, J = 10 | 57.5 s |
| Step 1 anchor, default cuts | 18.6 s |
| O1 cell n = 500 | 277.079 s |
| O1 cell n = 1000 | 637.086 s |
| O1 cell n = 2000 | 1769.901 s |
| O1 cell n = 4000 | 3681.577 s |
| **O1 total** | **6365.64 s = 1.77 h** |
| O2 batch 51–100 | **10475.55 s = 2.91 h** |
| **Total compute** | **≈ 4.70 h** |

Session span 01:15–06:25 local. O1 ran 01:40–03:26; O2 ran 03:26–06:22.

---

## 4. Errors and deviations

No R errors. No `fb_err` on any replicate. No failed cells.

Deviations from the instruction set, each recorded where it happened:

1. **`parallel_mode` was not a one-line change** (commit `perf(twin)`). The FB
   call passes `inner_parallel` unconditionally, so flipping the mode alone
   would nest 120 × 120 processes. `inner_parallel` is now conditional.

2. **MR could not be enabled in compare_all's main call** (Step 1 commit). Not
   only because the wrapper cannot scope `...` per arm, but because
   `.validate_mr_configuration()` *aborts* under `consistency` when
   `use_lasso`/`use_grf` is TRUE, which this document sets. Resolved by a
   scoped second call rather than by changing seven pre-existing arms.

3. **The J = 10 family was truncated** (G1 § 8). `max_subgroups_search = 30`
   discarded 98 of the 128 candidate subgroups the grids generate. Left
   unchanged as instructed; flagged in the spec append.

4. **The smoke-reproduction check was not runnable as specified** (§ 2.2), and
   no substitute invariance test was run.

5. **Only one O2 wave completed**, against 3–4 projected (§ 2.1).

## 5. Standing item not addressed

The twin's own introduction (`:34-39`) still asserts:

> For the consistency engine the candidate family is fixed, so FB and MR target
> the same estimand and should agree to leading order

`V_fb_adjudication.md` § V3 falsified this, and the package's own
`.fs_family_status()` roxygen says the same ("Deliberately NOT called 'fixed'
… quantile-derived cuts … are not"). The text was left in place because batch
1–50 was produced under it and editing it mid-series would make the bundles'
provenance inconsistent. It should be corrected before the next series starts.

G2 remains not started.
