# STOP B — md/harm MR coverage grid, n in {500, 1000, 2000, 4000}, s = 1000

Deliverable of the STOP B instruction (supplement items 3–4) and of
`CC_TASK_md_mr_harness.md`'s STOP B phase. Baseline `f0e4caf9`.

**Transport, not approval.**

## Headline

The pre-registered question was whether the ~13-SE residual MR bias measured at
n = 1000 in STOP A shrinks the way a leading-order error should.

**It does not. It grows, and interval coverage collapses with it.** Residual
bias rises 29.81 → 37.39 across the grid while its standard error falls, so the
bias/SE ratio goes 36 → 70, and MR coverage of the exact β(Ĥ) falls
0.957 → 0.731 against a nominal 0.95. Every machinery check passes on every
cell, so this is a property of the estimator, not a defect in the harness.

---

## Header

| | |
|---|---|
| Run | `run_tag = "md_grid_s1000"`, 4 cells x 1000 replicates, `mr_inference = TRUE`, `mr_draws = 5000`, `consistency_method = "resample"` |
| Bundles | `quarto/simulations/actg175/continuous/mr_sweep_md_harm/md_grid_s1000/fs_md_harm_n{500,1000,2000,4000}_res.rds` |
| Package | `forestsearch` 0.2.0 at `f0e4caf9`, R 4.6.1, 115 workers, PID-watched (PID 1446223, `kill -0`) |
| DGM | `cal_target_md = -40`, `n_super = 5000`, `seed_base = 8316951`, `adverse_outcome = FALSE`, `effect_threshold = 30`, `consistency_threshold = 10` |
| Render | completed; HTML produced, so all three checks passed on all four cells |

Truth block, identical across cells:

```
  marg_H        -40.0000000000      effect_Q      -40.0000000000
  marg_Hc       -26.2552358760      effect_Qc     -26.2552358760
  delta         -26.2552358760      beta_inter    -13.7447641240
  prevalence_Q    0.3446000000
```

All quantities below are on the **oriented** scale (positive = harm), per the
orientation bridge committed at `1830ca92`.

---

## 1. The pre-registered trend table

```
    n    s detected mr_ok naive_bias mr_bias   mr_se mr_bias_over_se coverage      mc_se   width
  500 1000      997   997     69.689  29.808 0.82703          36.043  0.95687  0.0064338 137.186
 1000 1000     1000  1000     65.491  33.202 0.71939          46.153  0.94900  0.0069569 121.518
 2000 1000     1000  1000     58.602  37.164 0.65816          56.467  0.88500  0.0100884 106.063
 4000 1000     1000  1000     48.774  37.388 0.53334          70.102  0.73100  0.0140228  86.193
```

### 1.1 What the columns say

**Naive bias falls with n** — 69.69 → 48.77 — as it should. More data, less
selection optimism.

**MR residual bias RISES with n** — 29.81 → 33.20 → 37.16 → 37.39. A
leading-order approximation error should shrink; this does the opposite and
then plateaus near 37. MR removes 57% of the naive bias at n = 500 and only
23% at n = 4000.

**bias/SE rises monotonically** — 36.0 → 46.2 → 56.5 → 70.1. There is no
reading of these numbers on which the residual is noise.

**Coverage degrades** — 0.957 → 0.949 → 0.885 → 0.731, against nominal 0.95
with MC SE ≤ 0.014. The n = 500 and n = 1000 cells are at or near nominal; by
n = 4000 more than a quarter of intervals miss.

**The mechanism is visible in `width`.** Interval width shrinks 137.19 → 86.19
as n grows, while the bias it must absorb stays near 30–37. Coverage at small n
was bought by width, not by centring — exactly what STOP A §2.3 warned when it
recorded 0.980 at n = 1000 with a mean width of 121.8. Widen-and-cover cannot
survive shrinking intervals.

### 1.2 Where the estimates actually sit

```
  n= 500  exact beta(Hhat) 31.7054   naive 101.3946   MR 61.5139
  n=1000  exact beta(Hhat) 31.8933   naive  97.3839   MR 65.0954
  n=2000  exact beta(Hhat) 31.9116   naive  90.5134   MR 69.0759
  n=4000  exact beta(Hhat) 31.8438   naive  80.6173   MR 69.2323
```

The target is essentially constant at ~31.9. The naive estimate converges
toward it from above. **The MR-corrected estimate moves AWAY from it**, from
61.5 to 69.2. MR is not converging on β(Ĥ) on this pathway.

### 1.3 Naive interval coverage, for contrast

```
    n cov_naive naive_bias naive_se
  500    0.2116      69.69   0.7036
 1000    0.1580      65.49   0.6374
 2000    0.1570      58.60   0.6157
 4000    0.1910      48.77   0.5235
```

The naive interval covers 16–21% of the time. MR is a large improvement on
doing nothing at every n; the finding is that it is not a *sufficient*
correction, and that its adequacy degrades as n grows.

---

## 2. Per-cell readouts

### 2.1 Detection and MR completion

| n | detected | mr_ok | status |
|---|---|---|---|
| 500 | 997/1000 | 997/1000 | DETECTED 997, NO-DETECTION 3 |
| 1000 | 1000/1000 | 1000/1000 | DETECTED 1000 |
| 2000 | 1000/1000 | 1000/1000 | DETECTED 1000 |
| 4000 | 1000/1000 | 1000/1000 | DETECTED 1000 |

`mr_ok == detected` on every cell — MR ran wherever a subgroup existed. Zero
`CONFIG-ERROR` rows, zero `err_msg`, zero MR-skip messages across all 4000
replicates.

The three `NO-DETECTION` replicates at n = 500 are the first departure from a
selection rate of 1.0 anywhere in this workstream. They are recorded as
NO-DETECTION, not as errors — the F1(c) status field distinguishing them is
what makes that statement checkable.

### 2.2 Identification vs Q, and selected size

```
    n    sens   spec    ppv    npv sel_n_med
  500 0.20800 0.8351 0.3942 0.6700        76
 1000 0.13085 0.9036 0.4122 0.6656        88
 2000 0.08444 0.9394 0.4137 0.6618       113
 4000 0.06183 0.9581 0.4071 0.6606       167
```

**Sensitivity falls with n** — 0.208 → 0.062. The selected region grows only
76 → 167 while Q grows in proportion to n (roughly 172 → 1378 subjects), so the
identifier captures a steadily smaller share of the true region. PPV is flat
near 0.41 and specificity rises mechanically as the selected region becomes a
smaller fraction of the sample. NPV is reported per project convention.

This is the identification counterpart of §1: the identifier is not converging
on Q, so β(Ĥ) — the effect at the realized rule — stays near 31.9 rather than
approaching θ† = −40 (oriented 40).

### 2.3 Partition and θ†

Partition holds on every row of every cell: `nH_eval + nHc_eval == 5000 == N`,
`betaHhat_status = "ok"` throughout.

θ† is unchanged by n (it is a population quantity at the TRUE flag):
`theta_dagger_H = -40.0000000000`, `theta_dagger_Hc = -26.2552358760`. It is
kept strictly distinct from β(Ĥ) at the realized rule (mean ~31.9 oriented,
i.e. ~-31.9 raw) throughout.

---

## 3. Checks A / B / C — per cell

All three pass on all four cells. The render completed and produced HTML, which
it could not do before the F1(e) guards and the `1830ca92` scale fix.

**CHECK A** — the two-valued CATE identity β(R) = δ + β_inter · P(Q | R) on the
harness's own realized rules:

```
===== CHECK A -- cell fs_md_harm_n500_res.rds  =====  distinct rules = 912;  max |exact - closed form| = 1.421e-14
===== CHECK A -- cell fs_md_harm_n1000_res.rds =====  distinct rules = 876;  max |exact - closed form| = 1.421e-14
===== CHECK A -- cell fs_md_harm_n2000_res.rds =====  distinct rules = 731;  max |exact - closed form| = 1.421e-14
===== CHECK A -- cell fs_md_harm_n4000_res.rds =====  distinct rules = 601;  max |exact - closed form| = 1.421e-14
```

Every cell at floating-point zero, over 601–912 distinct realized rules each.
The empty-rule guard was not triggered anywhere.

**CHECK B** — θ† reproduces the DGM's own effects (runs once, before the cells):

```
CHECK B: identical() = TRUE   max |difference| = 0.000e+00
```

**CHECK C** — the gate's naive bias and the recorder's selection effect, the
same quantity by two paths:

```
===== CHECK C -- cell fs_md_harm_n500_res.rds =====
  naive mean bias        (gate)      : 69.6892111165
  mean selection effect  (recorder)  : 69.6892111165
  difference                         : 0.000e+00
  max per-replicate |gate naive - recorder MD| : 3.268e-13
  PASS

===== CHECK C -- cell fs_md_harm_n1000_res.rds =====
  naive mean bias        (gate)      : 65.4906171178
  mean selection effect  (recorder)  : 65.4906171178
  difference                         : 0.000e+00
  max per-replicate |gate naive - recorder MD| : 3.411e-13
  PASS

===== CHECK C -- cell fs_md_harm_n2000_res.rds =====
  naive mean bias        (gate)      : 58.6018325393
  mean selection effect  (recorder)  : 58.6018325393
  difference                         : -7.105e-15
  max per-replicate |gate naive - recorder MD| : 7.532e-13
  PASS

===== CHECK C -- cell fs_md_harm_n4000_res.rds =====
  naive mean bias        (gate)      : 48.7735051935
  mean selection effect  (recorder)  : 48.7735051935
  difference                         : 0.000e+00
  max per-replicate |gate naive - recorder MD| : 8.029e-13
  PASS
```

CHECK C failed by 1.984e+02 at STOP A and now agrees to **exactly zero** on
three of four cells. That is the sixth edit (`1830ca92`) doing precisely what it
should: the two paths were always computing the same quantity, and once put on
the same scale they agree bit-for-bit. It also means the `naive_bias` column of
§1 is confirmed by two independent computations at every n.

---

## 4. Cost: measured vs the §6.2 linearity projection

```
    n    s elapsed_s s_per_rep_wall med_t2_secs workers linear_proj ratio_meas_over_proj
  500 1000     279.8         0.2798       10.71     115       20.15               0.5314
 1000 1000     641.2         0.6412       40.31     115       40.31               1.0000
 2000 1000    1747.3         1.7473      173.50     115       80.61               2.1522
 4000 1000    3673.4         3.6734      399.91     115      161.23               2.4804
```

Total grid wall clock **6341.7 s = 1 h 45.7 min**, against §6.2's projected
"30–45 minutes". **The linear-in-n assumption is refuted**, and in the expensive
direction: measured cost at n = 4000 is 2.48x the linear projection.

Two distinct errors compounded, and they are worth separating:

1. **The anchor was contaminated.** §6.2 anchored on STOP A's measurement of
   20.01 s per replicate at n = 1000. That cell ran 50 replicates on 115
   workers — under-subscribed, no contention. The same cell under STOP B's
   saturated load costs **40.31 s**, almost exactly 2x. So half the projection
   error is an artefact of measuring a parallel job that was not actually
   parallel-limited.
2. **The true scaling is superlinear.** Within STOP B, where every cell is
   equally saturated, median per-replicate cost goes 10.71 → 40.31 → 173.50 →
   399.91, i.e. factors of 3.76, 4.30, 2.30 per doubling — nearer n^1.8 than
   n^1. Candidate-family size is roughly n-independent, but each candidate's
   fit and its per-subject influence contributions are O(n), and MR's
   influence matrix B is N x S with a 5000-column multiplier matrix, so the
   MR stage carries its own O(n) term on top of the search.

Recorded so that the next schedule anchors on a saturated measurement and
assumes superlinear growth. No cell was reduced and `s = 1000` stood
throughout.

---

## 5. Reading, stated plainly

Within this DGM, this identifier, and this range of n:

- **MR is a large improvement over the naive estimate at every n** (coverage
  0.73–0.96 vs 0.16–0.21) and it is cheap.
- **MR is not consistent for β(Ĥ) here.** The residual bias does not shrink
  with n; it rises and plateaus near 37 while intervals narrow, so coverage
  degrades monotonically and is already at 0.731 by n = 4000. Extrapolating the
  trend, coverage would continue to fall at larger n.
- **The machinery is sound.** Every closed-form check passes at floating-point
  zero on every cell, the partition holds on all 4000 rows, and CHECK C now
  confirms the naive column by two independent paths. The finding is about the
  estimator, not the harness.

This is the honest STOP B result. It is not the pre-registered hope — the
13-SE residual did not behave like a leading-order error — and nothing here is
adjusted to make it look otherwise.

**Not established here**, and deliberately not investigated: whether the
residual reflects MR's leading-order truncation, the known superset gap in the
reconstructed candidate family (`R/forestsearch_main.R:3176-3181` records that
`d0.min`/`d1.min` and `max_subgroups_search` are not replayed, so the family MR
re-selects over is larger than the one the identifier chose among, which
inflates the selection-bias term), or the identifier's own drift away from Q
(§2.2). Those are three different explanations with three different remedies,
and separating them is a design question, not a measurement this grid can
settle.

---

## Gates

No DGM change. No `R/` change. No harness edit beyond the two this instruction
required (per-cell checks; the measured-cost table) and the grid knobs
`run_tag` / `n_grid_run` / `n_sims`, which are the document's own configuration
and which its comments anticipate. Closed decisions untouched: `focus = "harm"`,
`cal_target_md = -40`, `adverse_outcome = FALSE`, thresholds, seed scheme.
`s = 1000` per cell with no reduction. Deferred items (internal-default
hardening; resample-fallback warn-vs-error) not touched.

Committing this report is transport, not approval.
