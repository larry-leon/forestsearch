# REPORT — `dinamr` Block C, the deferred Block B cell, and the GRF cost probes

- **Date:** 2026-09-11. **Machine:** Mac Studio (darwin 25.6.0), R 4.5.2, `forestsearch` 0.3.5.
- **Task:** `dev/tasks/TASK_dinamr_blockC_grfprobe_2026-09-11.md` (committed as the first action of the session).
- **Predecessor:** `dev/tasks/TASK_dinamr_campaign_2026-09-10.md`, `REPORT_dinamr_2026-09-10.md`.
- **Commit only; not pushed. No `R/` change was made or needed.**

## Framing, restated

- **This does not certify DINA.** The fixed-family condition does not hold, so every coverage
  number in this record and in `summary_dinamr.qmd` is coverage of the
  **conditional-on-proposed-family estimand**, and every table says so.
- **Part B is cost measurement for GRF.** It is not a GRF campaign and not an adopted GRF
  configuration. No coverage claim, no acceptance criterion, no recommendation comes out of it.

---

# GATE 0 — Stage 0, verified from source

## Block C's cell definition in the predecessor

`dev/tasks/TASK_dinamr_campaign_2026-09-10.md:36`, the Part C block table:

```
| **C** | both prevalences | HR 1.00 at n = 500, 1000, 1500 — six cells |
```

with the defer order at line 38 ("Block C first (all six, 31% before 12.4%)"), and the grid's
provenance at line 10: "**Grid mirrors FS** (Larry, 2026-09-10): both prevalences, ε = 0.20, the
same HR and n structure as `cert20` / `tier2` / `p12ext`."

## What makes a planted-harm quantity undefined at HR 1.00 — and what does not

`FS_S7_HR` has exactly one meaning in this harness. Quoting the template
(`sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`):

- **line 340** — `target_hr_harm  <- .env_num("FS_S7_HR", 1.0)   # calibrate k_inter to this Cox HR in the harm subgroup`
- **lines 344–346** — the harm rule, in the comment block that documents the prevalence knob:
  ```
  # z1_quantile is the engine's own formal: the harm rule is
  #   {er <= quantile(er, z1_quantile)} & {meno == 0}
  # (the same form at every value).
  ```
- **line 353** — `harm_z1_quantile <- .env_num("FS_S7_Z1Q", 0.25)`
- **line 488** — `dgm_model       <- "alt"`, a literal; there is no `.env_*` read for it anywhere
  in the document.
- **line 661** — `harm_col       <- "flag_harm"`
- **lines 703–708** — `k_inter <- calibrate_k_inter(target_hr_harm = target_hr_harm, ...)` then
  `dgm <- setup_gbsg_dgm(model = dgm_model, k_inter = k_inter, z1_quantile = harm_z1_quantile, ...)`

So the harm region `H` is planted by the **rule**, which is a function of `FS_S7_Z1Q` alone; `HR`
enters only as the calibration **target for the effect inside that region**.
`calibrate_k_inter()` (`R/sim_aft_gbsg.R:1004`) roots `dgm$hr_H_true` on `target_hr_harm`.

The harness does own a genuine global-null path — `model = "null"`, documented at
`R/sim_aft_gbsg.R:93` as "uniform treatment effect", and accepted by
`setup_gbsg_dgm()` at `R/setup_gbsg_dgm.R:89` — but `dgm_model` is pinned to `"alt"` at template
line 488 and is not overridable, so **that path is unreachable from this template.**

### The calibration, run

| `z1q` | target HR | `k_inter` | prevalence | overall causal HR | HR(H) | HR(Hᶜ) |
|---|---|---|---|---|---|---|
| 0.25 (12.4%) | 1.00 | **+0.567917** | 0.1242 | 0.6847 | 1.0005 | 0.6569 |
| 0.25 (12.4%) | 1.50 | +1.111548 | 0.1242 | 0.7041 | 1.5086 | 0.6569 |
| 0.60 (31%)   | 1.00 | **−8.579297** | 0.3065 | 0.7922 | 0.9999 | 0.7206 |
| 0.60 (31%)   | 1.50 | −19.497297 | 0.3065 | 0.8694 | 1.4990 | 0.7206 |

Confirmed on a live 2-replicate render at (12.4%, HR 1.00, n 500): the bundle's `truth` reads
`hr_causal 0.6847, marg_H 1.00045, marg_Hc 0.65689, cde_H 1.0000, cde_Hc 0.58478`.

### Verdict, and the STOP condition

**HR 1.00 is a planted region with HR 1.00 inside it, not a global null.** `k_inter` is not zero
at either prevalence; the planted contrast is HR(H) = 1.00 against HR(Hᶜ) = 0.657 / 0.721. The
treatment effect is heterogeneous at these cells. What is absent is **harm** (HR > 1), not the
region.

**The STOP condition is not met.** The kickoff stops the task only if the predecessor's Block C
definition is *ambiguous* between the two readings. It is not: the predecessor pins the cells
through `FS_S7_HR` (its own `campaign.sh` sets `FS_S7_HR=$H` and nothing else), that knob has a
single documented meaning at template line 340, and the harness offers no way to express a global
null on this path. The reading is forced by the source, not chosen.

## Which quantities are structurally undefined at these cells

**None of the classification set.** Because `flag_harm` is planted by the rule and not by the
effect size, it resolves at HR 1.00 exactly as at HR 1.50 (super-population prevalence 0.1242 and
0.3065, printed by the template's own "Harm rule:" line). Therefore:

| quantity | status at HR 1.00 | how it is reported |
|---|---|---|
| oracle θ̂(H) (`or_*`) | **defined** — a Cox fit on the planted region; its target is log 1.00 = 0 | normally, as at every other cell |
| `sens`, `spec`, `ppv`, `npv` | **defined** — cross-tabulated against the same planted H | normally; finite on every detected replicate of the smoke |
| `betaHhat_H`, `betaHhat_Hc` | **defined** — realized targets of the returned region | normally |
| the nine `fld_recov_*` columns | **defined and populated** (9 of 9 on the smoke) | normally |
| `p_hat_*`, `fld_Hc_scale_ratio` (ρᶜ) | **defined and populated** | normally |
| `n_cons_qual` | present, all-NA | **structural on DINA** (no consistency screen) — @sec-structna, never a failure |
| `band_n` | present, all-NA | **structural on DINA** — never a failure |
| `p_star` | not a recorder column | **structural on DINA** (admission-set term, NULL on this engine) — never a failure |

The three structural columns are structural on **every** DINA cell, at every HR; HR 1.00 adds
nothing to the list. **No quantity is dropped, NA-ed or suppressed on account of HR 1.00.**

## Naming, settled — the rate is a `selection_rate` (Larry, 2026-09-11)

The planted region carries HR 1.00 and DINA's effect floor is log(0.90), **deliberately
sub-null**. A replicate that returns that region is therefore making an **admissible selection**
under the criterion in force; naming the rate for an error would mislabel correct behaviour. The
column is `selection_rate`, it is read as nothing more, and no other reading is carried in the
tables, the captions or this report.

**What speaks to a claim the data do not support is bound location, not selection.**
`summary_dinamr.qmd` @sec-null-location reports, per Block C cell, the **share of field lower
bounds at or above 1.00 and at or above 1.25**, with the rest of the columns of the Block A
location table (`scripts_dinamr/blockA_rest.R:63–76`): the median field lower bound against the
median realized θ(Ĥ) on the HR scale, the field point estimate and the naive estimate, the
bound-to-target gap in difference / ratio / paired-ratio form, and the planted marginal target.
Those columns carry the question. The `NUL` summary table's second share column moves from 0.85
to 1.25 to match.

## Vocabulary in Block C

These cells carry **no planted harm**, so no harm vocabulary is used for them anywhere. The
planted region is described as **differentially null against a benefiting complement** — HR 1.00
inside it, HR 0.657 (12.4%) / 0.721 (31%) outside it — and **every Block C caption states this
beside the conditional-on-proposed-family label**. The shared recorder columns built for the harm
blocks (`harm_*`, `comp_*`) are renamed to `region_*` / `compl_*` in the Block C tables.

---

# GATE 1 — the projection and the go/no-go

## The ceiling and its source

**Ceiling 9 h wall for Part A; hard timeout 12 h.** Source: the kickoff's "Gate 1 — compute
go/no-go" bullet, "**Ceiling 9 h wall for Part A; hard timeout 12 h.**" The 9 h figure originates
in **Amendment 1** of the predecessor (`TASK_dinamr_campaign_2026-09-10.md:68–72`), which replaced
the body's 10 h with 9 h on the grounds that FS projections have run +16% (`cert20`) to −29%
(`p12ext`) against realized and DINA's ~100× family-size skew is less predictable.

## Which original Gate 1 probes covered a Block C cell

**Four of the ten did**, so no new 36-replicate probe was needed. `probe.sh` lines 30–33:

```
run ""    500  1.00 p124_h100_n500
run ""   1500  1.00 p124_h100_n1500
run 0.60  500  1.00 p31_h100_n500
run 0.60 1500  1.00 p31_h100_n1500
```

These include both corners the kickoff would otherwise have had me run — (12.4%, HR 1.00, n 500)
and (31%, HR 1.00, n 1500). All four bundles are on disk.

## The probe corners at HR 1.00, as measured

| prevalence | n | reps | selection rate | K q10 | K med | K q90 | K max | s q10 | s med | s q90 | s max |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 12.4% | 500  | 36 | 0.7778 | 13.5 | 76.0  | 947.9  | 1195 | 0.0545 | 4.421 | 11.886 | 17.533 |
| 12.4% | 1500 | 36 | 0.3611 | 12.8 | 35.0  | 208.4  | 480  | 0.0790 | 0.113 | 6.417  | 9.754  |
| 31%   | 500  | 36 | 0.9444 | 31.9 | 359.5 | 1634.1 | 3053 | 3.2900 | 9.479 | 24.064 | 34.125 |
| 31%   | 1500 | 36 | 0.9444 | 50.7 | 258.5 | 1203.8 | 1699 | 4.2625 | 8.450 | 18.955 | 29.574 |

## Projecting from the distribution, not a mean

`projectC.R` (committed) bootstraps the 2,000-replicate total by resampling each corner's **36
measured per-replicate seconds** 4,000 times, so the projection carries an interval rather than
multiplying a mean. n = 1000 is not probed at HR 1.00 and is interpolated by **pooling the n 500
and n 1500 draws**, which interpolates the distribution rather than only its centre.

### Realized walls come from the bundles' timing columns, not from mtimes

The first pass reconstructed realized walls by differencing result-file mtimes. **The bundles
carry per-replicate timing, so that was never necessary**, and `walls.R` (committed) replaces it.
Verified on the committed bundles before using them:

- `fit_mr_secs` is the **top-level per-replicate worker timer** and is **finite on all 2,000 rows
  of every cell** (1000/1000 per batch), selected or not.
- `fld_H_secs` and `fld_Hc_secs` are **nested inside it** — `fld_H_secs <= fit_mr_secs` on **878 of
  878** rows where both are finite (median ratio 0.829, max 0.940), and
  `fld_H_secs + fld_Hc_secs <= fit_mr_secs` on all of them. They must **not** be added.
- `fb_secs` and `fld_H_uniform_secs` are identically zero on this campaign (`FS_S7_FB=none`).
- `meta$n_workers` = **12** is recorded in every **batch** meta (the combined meta does not carry
  it).

There is no per-replicate *wall* column, so the compute wall is
`sum(fit_mr_secs) / n_workers`. The mtime span is kept beside it **labelled a proxy**, and its
residual over the compute wall is reported rather than absorbed.

| block | HR | n | compute wall (h) | mtime proxy (h) | residual (s) | per render (s) |
|---|---|---|---|---|---|---|
| A | 1.50 | 500  | 0.4154 | 0.5177 | 368.5 | 122.8 |
| A | 1.50 | 1000 | 0.3185 | 0.3888 | 253.1 | 84.4 |
| A | 1.50 | 1500 | 0.2650 | 0.3295 | 232.5 | 77.5 |
| A | 1.75 | 500  | 0.4783 | 0.5723 | 338.5 | 112.8 |
| A | 1.75 | 1000 | 0.3980 | 0.4700 | 259.5 | 86.5 |
| A | 1.75 | 1500 | 0.3652 | 0.4285 | 227.8 | 75.9 |
| B | 1.50 | 500  | 1.2314 | 1.3952 | 589.8 | 196.6 |
| B | 1.50 | 1000 | 1.4453 | 1.5471 | 366.7 | 122.2 |
| B | 1.50 | 1500 | 1.5038 | 1.6268 | 442.5 | 147.5 |
| B | 1.75 | 500  | 1.4410 | 1.5812 | 504.5 | 168.2 |
| B | 1.75 | 1000 | 1.7438 | 1.8361 | 332.3 | 110.8 |

Totals over the eleven cells: **compute 9.606 h**, mtime proxy 10.693 h, residual **1.088 h**.

**The residual is the finding.** `project.R` assumed **30 s** of per-render overhead; the measured
value is a **median of 112.8 s, range 75.9–196.6** — about four times the assumption. It is
render / DGM-build / table cost plus the makespan slack over the `sum/W` bound. Once it is
accounted for, **Gate 1's compute model turns out to have been near-exact**: on the measured
HR 1.50 corners `sum(compute)/sum(gate1) = 0.972`, while the mtime-based multiplier of **1.090**
on those same corners **was overhead, not compute**. (The HR 1.75 cells are excluded from that
reading — they were costed at the HR 1.50 corner, so their ratio absorbs that assumption rather
than measuring anything.)

`projectC.R` therefore **drops the blanket multiplier entirely** and projects
compute-from-distribution **plus the measured per-block overhead**, taken at the largest realized
value per block (12.4% → 122.8 s, 31% → 196.6 s), which is the conservative choice.

## Part A projection

| cell | basis | overhead/render | median (h) | 90% band | at project.R's 30 s |
|---|---|---|---|---|---|
| C 12.4% n 500  | measured | 122.8 s | 0.3462 | 0.3384–0.3543 | 0.2690 |
| C 12.4% n 1000 | pooled draws | 122.8 s | 0.2690 | 0.2616–0.2762 | 0.1918 |
| C 12.4% n 1500 | measured | 122.8 s | 0.1920 | 0.1872–0.1968 | 0.1147 |
| C 31% n 500    | measured | 196.6 s | 0.7084 | 0.6937–0.7228 | 0.5698 |
| C 31% n 1000   | pooled draws | 196.6 s | 0.6716 | 0.6591–0.6842 | 0.5327 |
| C 31% n 1500   | measured | 196.6 s | 0.6347 | 0.6238–0.6455 | 0.4956 |
| **Block C, six cells** | | | **2.822** | 2.764–2.880 | **2.174** |

The final column reproduces the original Gate 1's **2.17 h** for Block C exactly — confirming that
the correction is entirely in the overhead constant, not in the compute model.

**The deferred Block B cell (HR 1.75, n 1500, 31%)** is anchored on realized **compute** walls two
ways, which agree to the second: (a) the n-profile within HR 1.75,
`compute(1.75,1000) × compute(1.50,1500)/compute(1.50,1000)` = 6532 s; (b) the HR-profile within
n = 1500, `compute(1.50,1500) × compute(1.75,1000)/compute(1.50,1000)` = 6532 s. Plus
3 × 196.6 s = **7122 s = 1.978 h**. References: the checkpoint's by-n re-projection 2.225 h
(corrected down), the original Gate 1 1.617 h (corrected up).

| | h |
|---|---|
| deferred Block B cell | 1.978 |
| Block C, six cells | 2.822 |
| **PART A TOTAL** | **4.800** |
| ceiling | 9.000 |
| headroom | 4.200 (47%) |
| room left under the 12 h timeout for Part B's 1.5 h cap | 7.200 |

Under the superseded mtime calibration this read 4.701 h. **Gate 1 is GO either way and no cell is
deferred**, so the running campaign is unaffected by the correction.

**GATE 1: GO — all seven cells run, none deferred.**

---

# Tooling

Everything is in `quarto/simulations/gbsg_020/scripts_dinamr/` and committed; nothing was left in
the session scratchpad. The committed drivers were reused verbatim (`campaign.sh`, `render.sh`,
`gate2.R`, `project.R`, `probe.sh` — `diff` clean against the copies actually executed).

| file | status | role |
|---|---|---|
| `blockC.cells` | **new** | the six Block C cells, in the kickoff's run order |
| `blockB_deferred.cells` | **new** | the one cell the Block A checkpoint deferred |
| `projectC.R` | **new** | Gate 1 for these seven: bootstrap over the per-replicate cost distribution, calibrated on realized walls |
| `grfprobe.sh` | **new** | the five Part B GRF corners, each under `/usr/bin/time -l` with a process-tree RSS sampler |
| `grfprobe.R` | **new** | Gate 3 plus the Part B cost readout |
| `grf_mechanism.R` | **new** | the empty-band diagnostic (below) |
| `gate2.R` | **one-line correction** | `fscomp()` sent Block B's n = 500 cell to `e1stud` unconditionally; `e1stud` was run only at HR 1.50 and 1.75, so at HR 1.00 that named a nonexistent file and the Amendment 3 assertion would have gone unevaluated on one of the six Block C cells. Now `cert20` there, which is on disk and criterion-matched. Regression: `gate2.R A` still reads **204 passes, 0 failures**. |

The `gate2.R` comparator fix is **accepted and recorded** (Larry, 2026-09-11).

`stage1_checks.R` still carries the superseded bonf-vs-raw comparison, labelled as such; every
gate below used `gate2.R`'s corrected identity
(`log(fld_Hc_est2_s) + fld_Hc_lam_mean_s == log(fld_Hc_est2) + fld_Hc_lam_mean`).

**The `SP` export defect, since fixed.** `render.sh` reads `$SP` for its log directory and its own
header states the caller exports what it needs, but `campaign.sh`, `probe.sh` **and** `grfprobe.sh`
all set `SP` as a plain shell variable. In a clean environment `$SP` reached `render.sh` empty and
it died on `mkdir -p /logs` under `set -e`. The Part A and Part B runs in this session exported
`SP` around the scripts; all three callers are now fixed in the tree (`render.sh` unchanged — its
contract puts the export on the caller).

Verified by **execution**, not by reading, in `env -i` with `DINAMR_SCRATCH` unset and a stub
template of the expected filename under `DINAMR_QMD_DIR`:

| | result |
|---|---|
| `campaign.sh` **before** | `mkdir: /logs: Read-only file system` → `CELL FAILED` |
| `campaign.sh` **after** | `batch_1`, `batch_1001`, `combine_1` all `RC=0`, `CELL DONE  wall=4s`, three logs under `<scripts>/logs/` |
| `probe.sh` after | `RC=0`, `PROBE DONE`, `probe_p124_h150_n500.log` written |
| `grfprobe.sh` after | `RC=0`, `GRF PROBE DONE`, plus its `rss_*.txt` sampler and `time_*.txt` sidecars resolved into the same directory |

Nothing was written to `/logs` in any run.

---

# PART A — results

## A note on walls, now that the driver's own figure is available

The campaign driver prints `CELL DONE: <cell>  wall=<n>s` — a **directly measured** wall, better
than either reconstruction. For the seven cells of this task that figure is used, and it also
gives a ground-truth check on the timing-column method used at Gate 1. On the deferred B cell:

| | s | h |
|---|---|---|
| compute wall, `sum(fit_mr_secs)/12` | 6830.0 | 1.8972 |
| residual (3 renders) | 360.0 | |
| **driver wall** | **7190** | **1.9972** |

Compute + residual = 7190 s **exactly**, at **120.0 s per render** — squarely on the 112.8 s median
measured across the eleven earlier cells, and comfortably inside the conservative 196.6 s the
projection used for the 31% block. The timing-column reconstruction is confirmed.

## Cell 1 of 7 — the deferred Block B cell (HR 1.75, n 1500, 31%)

**Wall 7190 s = 1.997 h against a projected 1.978 h — realized/projected 1.010.**

**Gate 2: PASS on every check.** With this cell in place `Rscript gate2.R B` reads **204 passes, 0
failures** over all six Block B cells, none missing.

- completeness: 2,000 rows, `sim_id` 1–2000 no duplicates, no CONFIG-ERROR, two seed-disjoint
  batches, `n_workers` 12 and `forestsearch_version` 0.3.5 recorded in both, seed_base 8316951,
  host Mac-Studio-3.local, R 4.5.2; every campaign knob as set.
- **detection 1.0000 (2000/2000)**; **proposed-family size** min 384, q10 1341.9, med 2355.5,
  q90 3108, max 3726 (mean 2276.6, CV 0.296); realized prevalence 0.30646 trial against 0.30655
  super-population.
- all 38 products finite on detected replicates; nine recovery columns, the p̂ block and ρᶜ all
  present and populated; every interval invariant holds; γ ∈ [0.02500, 0.02700] on both joints.
- **the corrected identity** `log(est2_s) + lam_mean_s == log(est2) + lam_mean`: max |diff|
  **2.78e-16**. The `bonf == raw` identity, gated on the γ-at-floor rows: max |diff| **0**, share
  at floor 0.890 / 0.795.
- classification: sens 0.8215, spec 0.9548, PPV 0.8886, NPV 0.9314, mean |Ĥ| 424.5.
- structurally-NA reported as such: `n_cons_qual` and `band_n` present all-NA, `p_star` not a
  recorder column. **Non-detections: none.**
- **Amendment 3, same draws vs `cert20` (effMaxSG, ε 0.2 — criterion-matched):** `n_true`
  `identical()` on all 2,000 rows **YES**; `truth` `all.equal(tol = 1e-8)` **YES**;
  `identical()` FALSE with max |abs diff| **8.882e-15**, max |rel diff| **4.253e-15** — the
  cross-machine BLAS difference the checker documents, not a mismatch.
- family contrast at this cell: FS enumerated min 1195, med 1297, q90 1396, max 1412, CV 0.0357;
  DINA min 384, med 2355.5, q90 3108, max 3726, CV 0.2958. Both detect 1.0000.

## Cells 2–7 — Block C: all six complete, Gate 2 204 passes / 0 failures

`Rscript gate2.R C` reads **204 passes, 0 failures, 0 cells not on disk**. Every check that passed
on the deferred B cell passes on all six: completeness, every meta knob, all 38 products finite on
selected replicates, the nine recovery columns / p̂ block / ρᶜ present and populated, every interval
invariant, γ inside [0.025, 0.028] on both joints, and the corrected identity
`log(est2_s) + lam_mean_s == log(est2) + lam_mean` at **1.67e-16 to 3.33e-16**. `n_cons_qual` and
`band_n` are reported present-and-all-NA and `p_star` as not-a-recorder-column, structural on DINA
at every HR and never as failures.

### Selection, family, classification

| cell | wall | selection rate | family: med (q90, max), CV | sens / spec / PPV / NPV | mean \|Ĥ\| |
|---|---|---|---|---|---|
| 12.4% n 500  | 1183 s | 0.7135 | 137 (746, 3224), 1.422 | 0.346 / 0.825 / 0.218 / 0.900 | 98.4 |
| 12.4% n 1000 | 712 s  | 0.5245 | 66 (355, 1611), 1.422 | 0.378 / 0.852 / 0.269 / 0.907 | 176.3 |
| 12.4% n 1500 | 445 s  | 0.3435 | 36 (211, 1186), 1.458 | 0.421 / 0.866 / 0.314 / 0.914 | 255.3 |
| 31% n 500    | 2710 s | 0.9300 | 485.5 (1965, 3476), 1.008 | 0.356 / 0.849 / 0.497 / 0.753 | 106.9 |
| 31% n 1000   | 2290 s | 0.9260 | 349 (1348, 3739), 1.083 | 0.457 / 0.870 / 0.593 / 0.791 | 229.9 |
| 31% n 1500   | 1876 s | 0.8880 | 278.5 (896, 3364), 1.078 | 0.578 / 0.871 / 0.656 / 0.833 | 399.6 |

Realized trial prevalence tracks the super-population at every cell (0.12363–0.12401 against
0.12418; 0.30591–0.30646 against 0.30655). Non-detections are 573 / 951 / 1313 and 140 / 148 / 224,
all recorded `NO-DETECTION` with `n_family` NA — never 0, never positive — so the two causes are
not separable from the committed columns, as `gate2.R` states.

### Bound location — the columns that carry the question

| cell | median lower bound | median θ(Ĥ) | bound / θ | paired-ratio med | **share ≥ 1.00** | **share ≥ 1.25** |
|---|---|---|---|---|---|---|
| 12.4% n 500  | 0.4924 | 0.6833 | 0.721 | 0.716 | **0.0189** | **0.0035** |
| 12.4% n 1000 | 0.5721 | 0.7096 | 0.806 | 0.806 | **0.0162** | **0.0010** |
| 12.4% n 1500 | 0.6119 | 0.7267 | 0.842 | 0.860 | **0.0087** | **0.0015** |
| 31% n 500    | 0.5556 | 0.8320 | 0.668 | 0.665 | **0.0301** | **0.0086** |
| 31% n 1000   | 0.6206 | 0.8767 | 0.708 | 0.716 | **0.0157** | **0.0032** |
| 31% n 1500   | 0.6705 | 0.9016 | 0.744 | 0.762 | **0.0096** | **0.0017** |

**This is the Block C headline.** Selection is frequent and admissible — the region clears the
sub-null floor — but the lower bound almost never reaches a level that would assert anything the
data do not carry. The share at or above 1.00 runs **0.9%–3.0%** and the share at or above 1.25
runs **0.10%–0.86%**, and **both fall monotonically with n at both prevalences**. The two
quantities move independently and in opposite directions at 12.4%: the selection rate collapses
0.7135 → 0.5245 → 0.3435 while the bound tightens toward the realized target (0.721 → 0.806 →
0.842); at 31% selection is roughly flat (0.9300 → 0.9260 → 0.8880) while the location shares
still fall by a factor of three. **Selection rate and bound location are not interchangeable
readings, and only the second speaks to a claim the data do not support.**

### The proposed family

DINA's family at these cells is small and extremely volatile — **CV 1.008 to 1.458** — and it
**shrinks with n** (12.4%: 137 → 66 → 36; 31%: 485.5 → 349 → 278.5), with a minimum of 1 at every
cell. FS's family on the **identical draws** is an order of magnitude larger, essentially flat in
n, and two orders of magnitude steadier: median 1223–1300 with **CV 0.036–0.054**. This is the
family-construction difference the confound statement names, measured.

### Coverage of every product, absolute levels, beside the FS comparator

Conditional-on-proposed-family throughout; Block C is differentially null against a benefiting
complement. FS is read as it stands, with its own criterion named — **criterion-matched at 31%,
not at 12.4%**.

| cell | matched | FS criterion | DINA sel / FS sel | DINA field / FS field | DINA IJ 2-sided [Wilson] / FS IJ |
|---|---|---|---|---|---|
| 12.4% n 500  | no  | tier2 / maxeffCons / 0.1  | 0.7135 / 0.6805 | 0.8746 / 0.9625 | 0.9916 [0.9854, 0.9952] / 0.9860 |
| 12.4% n 1000 | no  | p12ext / maxeffCons / 0.1 | 0.5245 / 0.6595 | 0.8494 / 0.9204 | 0.9886 [0.9801, 0.9934] / 0.9803 |
| 12.4% n 1500 | no  | p12ext / maxeffCons / 0.1 | 0.3435 / 0.6240 | 0.8253 / 0.9303 | 0.9898 [0.9791, 0.9951] / 0.9824 |
| 31% n 500    | **yes** | cert20 / effMaxSG / 0.2 | 0.9300 / 0.9205 | 0.9059 / 0.9734 | 0.9919 [0.9867, 0.9951] / 0.9946 |
| 31% n 1000   | **yes** | cert20 / effMaxSG / 0.2 | 0.9260 / 0.9545 | 0.9330 / 0.9560 | 0.9941 [0.9894, 0.9967] / 0.9948 |
| 31% n 1500   | **yes** | cert20 / effMaxSG / 0.2 | 0.8880 / 0.9590 | 0.9471 / 0.9666 | 0.9966 [0.9926, 0.9985] / 0.9964 |

| cell | DINA field-s upper / FS | DINA IJ 2-sided complement / FS | DINA joint_s Bonferroni / FS |
|---|---|---|---|
| 12.4% n 500  | 0.9299 / 0.9398 | 1.0000 / 1.0000 | 0.9075 / 0.9566 |
| 12.4% n 1000 | 0.9466 / 0.9560 | 1.0000 / 1.0000 | 0.8875 / 0.9416 |
| 12.4% n 1500 | 0.9534 / 0.9495 | 1.0000 / 1.0000 | 0.8821 / 0.9431 |
| 31% n 500    | 0.9027 / 0.9115 | 1.0000 / 0.9989 | 0.9027 / 0.9430 |
| 31% n 1000   | 0.9320 / 0.9277 | 0.9989 / 1.0000 | 0.9325 / 0.9371 |
| 31% n 1500   | 0.9358 / 0.9270 | 1.0000 / 1.0000 | 0.9414 / 0.9526 |

The one-sided field lower bound sits below its nominal 0.95 on the region at every Block C cell
(0.825–0.947) and rises with n; the field-s complement upper bound sits at 0.903–0.953; the IJ
two-term two-sided interval is conservative on both blocks (0.989–0.997 on the region, ~1.000 on
the complement); the joint_s Bonferroni pair runs 0.882–0.941. **These are absolute levels,
recorded, not scored** — no acceptance criterion applies and the FS column is a reference line
carrying its own criterion, not a bar.

### The Block C coverage numbers

These are the reason Block C was run. All from the rendered summary's own definitions
(`cov-fns`, `strat-fns`, `strat-xtab`) applied to the committed bundles. **Absolute levels; every
column is coverage of the conditional-on-proposed-family estimand, on selected replicates only.
Block C is differentially null against a benefiting complement — HR 1.00 in the planted region,
0.657 (12.4%) / 0.721 (31%) in its complement.**

`bias_log` is the retained bias on the log scale; `sd_emp` the marginal SD; `sd_err` the error SD;
`se_mean` the mean SE; `cov1` the one-sided 95% coverage on the exposed side with Wilson limits.

**12.4%, β(Ĥ) on the identified region — one-sided lower**

| cell | construction | n | bias_log | sd_emp | sd_err | se_mean | cov1 [Wilson] |
|---|---|---|---|---|---|---|---|
| n 500  | naive | 1427 | 0.6711 | 0.2808 | 0.3016 | 0.3114 | 0.2649 [0.2426, 0.2884] |
| n 500  | field | 1427 | 0.2311 | 0.3016 | 0.3095 | 0.3138 | 0.8746 [0.8564, 0.8907] |
| n 500  | IJ two-term | 1427 | 0.2995 | 0.2851 | 0.2951 | 0.4490 | 0.9650 [0.9541, 0.9733] |
| n 1000 | naive | 1049 | 0.4725 | 0.1899 | 0.2021 | 0.2247 | 0.3041 [0.2770, 0.3326] |
| n 1000 | field | 1049 | 0.1872 | 0.2388 | 0.2330 | 0.2290 | 0.8494 [0.8265, 0.8697] |
| n 1000 | IJ two-term | 1049 | 0.2315 | 0.2148 | 0.2111 | 0.3249 | 0.9647 [0.9518, 0.9743] |
| n 1500 | naive | 687 | 0.3676 | 0.1391 | 0.1519 | 0.1832 | 0.2969 [0.2640, 0.3322] |
| n 1500 | field | 687 | 0.1616 | 0.1976 | 0.1945 | 0.1821 | 0.8253 [0.7951, 0.8519] |
| n 1500 | IJ two-term | 687 | 0.1947 | 0.1731 | 0.1702 | 0.2595 | 0.9738 [0.9590, 0.9834] |

**12.4%, β(Ĥᶜ) on the complement — one-sided upper**

| cell | construction | bias_log | sd_emp | sd_err | se_mean | cov1 [Wilson] |
|---|---|---|---|---|---|---|
| n 500  | naive | −0.09331 | 0.1271 | 0.1254 | 0.1425 | 0.8641 [0.8453, 0.8809] |
| n 500  | field | −0.02388 | 0.1389 | 0.1373 | 0.1399 | 0.9271 [0.9125, 0.9395] |
| n 500  | field-s | −0.02399 | 0.1387 | 0.1371 | 0.1416 | 0.9299 [0.9155, 0.9420] |
| n 500  | IJ two-term | −0.03607 | 0.1358 | 0.1343 | 0.2691 | 1.0000 [0.9973, 1.0000] |
| n 1000 | naive | −0.04397 | 0.0894 | 0.0868 | 0.0993 | 0.9180 [0.8999, 0.9331] |
| n 1000 | field | −0.00755 | 0.0967 | 0.0948 | 0.0978 | 0.9438 [0.9281, 0.9561] |
| n 1000 | field-s | −0.00759 | 0.0966 | 0.0947 | 0.0989 | 0.9466 [0.9313, 0.9587] |
| n 1000 | IJ two-term | −0.01341 | 0.0948 | 0.0928 | 0.1908 | 0.9981 [0.9931, 0.9995] |
| n 1500 | naive | −0.02630 | 0.0725 | 0.0717 | 0.0809 | 0.9374 [0.9168, 0.9532] |
| n 1500 | field | −0.00407 | 0.0769 | 0.0770 | 0.0800 | 0.9520 [0.9333, 0.9656] |
| n 1500 | field-s | −0.00407 | 0.0769 | 0.0770 | 0.0808 | 0.9534 [0.9350, 0.9668] |
| n 1500 | IJ two-term | −0.00755 | 0.0757 | 0.0757 | 0.1570 | 1.0000 [0.9944, 1.0000] |

**31% — field lower on β(Ĥ), field-s upper on β(Ĥᶜ), Bonferroni joint, IJ, naive**

| cell | field lower [Wilson] | field-s upper [Wilson] | joint Bonferroni [Wilson] | joint-s Bonferroni [Wilson] |
|---|---|---|---|---|
| n 500  | 0.9059 [0.8918, 0.9184] | 0.9027 [0.8884, 0.9153] | 0.8957 [0.8810, 0.9088] | 0.9027 [0.8884, 0.9153] |
| n 1000 | 0.9330 [0.9207, 0.9436] | 0.9320 [0.9196, 0.9426] | 0.9255 [0.9126, 0.9366] | 0.9325 [0.9202, 0.9431] |
| n 1500 | 0.9471 [0.9357, 0.9566] | 0.9358 [0.9234, 0.9463] | 0.9381 [0.9259, 0.9484] | 0.9414 [0.9295, 0.9514] |

**Joint and naive, all six cells**

| cell | joint Bonferroni | joint-s Bonferroni | naive β(Ĥ) two-sided | naive β(Ĥᶜ) two-sided |
|---|---|---|---|---|
| 12.4% n 500  | 0.9026 [0.8861, 0.9169] | 0.9075 [0.8914, 0.9215] | 0.4205 [0.3951, 0.4463] | 0.9285 [0.9140, 0.9408] |
| 12.4% n 1000 | 0.8856 [0.8649, 0.9035] | 0.8875 [0.8670, 0.9052] | 0.4681 [0.4380, 0.4983] | 0.9561 [0.9420, 0.9670] |
| 12.4% n 1500 | 0.8821 [0.8558, 0.9041] | 0.8821 [0.8558, 0.9041] | 0.4789 [0.4418, 0.5163] | 0.9651 [0.9485, 0.9764] |
| 31% n 500    | 0.8957 [0.8810, 0.9088] | 0.9027 [0.8884, 0.9153] | 0.4210 [0.3987, 0.4435] | 0.8726 [0.8567, 0.8870] |
| 31% n 1000   | 0.9255 [0.9126, 0.9366] | 0.9325 [0.9202, 0.9431] | 0.5221 [0.4994, 0.5448] | 0.9055 [0.8913, 0.9180] |
| 31% n 1500   | 0.9381 [0.9259, 0.9484] | 0.9414 [0.9295, 0.9514] | 0.6368 [0.6142, 0.6589] | 0.9150 [0.9011, 0.9271] |

**IJ two-term two-sided, miss split by side**

| cell | IJ β(Ĥ) [Wilson] | miss below | miss above | IJ β(Ĥᶜ) [Wilson] | miss below | miss above |
|---|---|---|---|---|---|---|
| 12.4% n 500  | 0.9916 [0.9854, 0.9952] | 0.00771 | 0.00070 | 1.0000 [0.9973, 1.0000] | 0 | 0 |
| 12.4% n 1000 | 0.9886 [0.9801, 0.9934] | 0.01144 | 0 | 1.0000 [0.9964, 1.0000] | 0 | 0 |
| 12.4% n 1500 | 0.9898 [0.9791, 0.9951] | 0.01019 | 0 | 1.0000 [0.9944, 1.0000] | 0 | 0 |
| 31% n 500    | 0.9919 [0.9867, 0.9951] | 0.00807 | 0 | 1.0000 [0.9979, 1.0000] | 0 | 0 |
| 31% n 1000   | 0.9941 [0.9894, 0.9967] | 0.00594 | 0 | 0.9989 [0.9961, 0.9997] | 0 | 0.00108 |
| 31% n 1500   | 0.9966 [0.9926, 0.9985] | 0.00282 | 0.00056 | 1.0000 [0.9978, 1.0000] | 0 | 0 |

**What these say.** The naive interval on the identified region is the outlier: 0.42–0.64 two-sided
against a retained bias of +0.368 to +0.671 log units — the optimism the correction exists to
remove. The field lower bound runs **0.825–0.947**, below nominal at every cell, rising with n at
31% (0.906 → 0.933 → 0.947) and *falling* with n at 12.4% (0.875 → 0.849 → 0.825). The field-s
complement upper bound runs 0.903–0.953. Both Bonferroni joints run 0.882–0.941. The IJ two-term
interval is conservative on the region (0.989–0.997) and essentially saturated on the complement,
its mean SE running 1.5–2.0× the error SD; its misses are almost entirely **below**, i.e. the
interval sits above the realized target.

### Field lower bound by family-size tertile and by p̂ bin

Block C has the most volatile family in the grid (CV 1.008–1.458), so this is where the
stratification matters. Coverage of β(Ĥ) by the field lower bound with Wilson limits; retained bias
is `mean(log(fld_H_est2) − log(betaHhat_H))` within the stratum. K tertiles partition; `K = 1`,
`K <= 5` and `all detected` **overlap** them and must not be summed.

**By `n_family` tertile**

| cell | stratum | n | field cov [Wilson] | retained bias |
|---|---|---|---|---|
| 12.4% n 500 | K T1 [1, 58] | 476 | 0.8824 [0.8503, 0.9083] | 0.1491 |
| | K T2 [59, 268] | 476 | 0.8824 [0.8503, 0.9083] | 0.2427 |
| | K T3 [269, 3224] | 475 | 0.8589 [0.8248, 0.8874] | 0.3018 |
| | K = 1 *(ov)* | 24 | 0.9583 [0.7976, 0.9926] | **−0.0767** |
| | K ≤ 5 *(ov)* | 82 | 0.8902 [0.8044, 0.9412] | 0.0551 |
| 12.4% n 1000 | K T1 [1, 36] | 357 | 0.8908 [0.8541, 0.9191] | 0.1216 |
| | K T2 [37, 125] | 342 | 0.8450 [0.8029, 0.8795] | 0.1971 |
| | K T3 [127, 1611] | 350 | 0.8114 [0.7672, 0.8489] | 0.2443 |
| | K = 1 *(ov)* | 32 | 0.9062 [0.7578, 0.9676] | **−0.0864** |
| 12.4% n 1500 | K T1 [1, 18] | 231 | 0.8225 [0.7681, 0.8664] | 0.1006 |
| | K T2 [19, 73] | 229 | 0.8253 [0.7709, 0.8690] | 0.1779 |
| | K T3 [74, 1186] | 227 | 0.8282 [0.7738, 0.8717] | 0.2073 |
| | K = 1 *(ov)* | 43 | 0.8372 [0.7003, 0.9188] | **−0.0038** |
| 31% n 500 | K T1 [1, 272] | 620 | 0.9403 [0.9188, 0.9564] | 0.0650 |
| | K T2 [273, 872] | 620 | 0.9097 [0.8845, 0.9298] | 0.1744 |
| | K T3 [873, 3476] | 620 | 0.8677 [0.8388, 0.8922] | 0.2536 |
| | K = 1 *(ov)* | 7 | 1.0000 [0.6457, 1.0000] | **−0.2044** |
| 31% n 1000 | K T1 [1, 201] | 619 | 0.9661 [0.9487, 0.9777] | 0.0069 |
| | K T2 [202, 593] | 617 | 0.9481 [0.9277, 0.9630] | 0.0639 |
| | K T3 [594, 3739] | 616 | 0.8847 [0.8571, 0.9076] | 0.1533 |
| 31% n 1500 | K T1 [1, 167] | 593 | 0.9747 [0.9587, 0.9846] | **−0.0124** |
| | K T2 [168, 407] | 591 | 0.9560 [0.9363, 0.9698] | 0.0283 |
| | K T3 [408, 3364] | 592 | 0.9105 [0.8848, 0.9309] | 0.0921 |

**By p̂ tertile**

| cell | stratum | n | field cov [Wilson] | retained bias |
|---|---|---|---|---|
| 12.4% n 500 | p̂ T1 [0.0018, 0.0904] | 477 | 0.9853 [0.9700, 0.9929] | 0.1282 |
| | p̂ T2 [0.091, 0.2284] | 474 | 0.9241 [0.8966, 0.9446] | 0.2207 |
| | p̂ T3 [0.229, 0.9886] | 476 | **0.7143 [0.6721, 0.7530]** | 0.3446 |
| 12.4% n 1000 | p̂ T1 | 350 | 0.9400 [0.9100, 0.9604] | 0.1498 |
| | p̂ T2 | 349 | 0.8539 [0.8130, 0.8871] | 0.1980 |
| | p̂ T3 | 350 | **0.7543 [0.7066, 0.7965]** | 0.2136 |
| 12.4% n 1500 | p̂ T1 | 229 | 0.8908 [0.8438, 0.9250] | 0.1496 |
| | p̂ T2 | 229 | 0.8472 [0.7949, 0.8880] | 0.1610 |
| | p̂ T3 | 229 | **0.7380 [0.6774, 0.7907]** | 0.1742 |
| 31% n 500 | p̂ T1 [2e-04, 0.066] | 622 | 0.9936 [0.9836, 0.9975] | 0.0058 |
| | p̂ T2 | 618 | 0.9612 [0.9429, 0.9738] | 0.1529 |
| | p̂ T3 | 620 | **0.7629 [0.7279, 0.7947]** | 0.3348 |
| 31% n 1000 | p̂ T1 | 618 | 0.9935 [0.9835, 0.9975] | −0.0085 |
| | p̂ T2 | 617 | 0.9708 [0.9544, 0.9815] | 0.0648 |
| | p̂ T3 | 617 | **0.8347 [0.8033, 0.8619]** | 0.1675 |
| 31% n 1500 | p̂ T1 | 595 | 0.9983 [0.9905, 0.9997] | −0.0133 |
| | p̂ T2 | 589 | 0.9626 [0.9441, 0.9752] | 0.0517 |
| | p̂ T3 | 592 | **0.8801 [0.8514, 0.9038]** | 0.0697 |

**The p̂ stratification separates far more sharply than the family-size one.** By K the spread
across tertiles is 2–7 points; by p̂ it is **17–27 points** at every cell, monotone, with the
high-p̂ tertile at 0.714–0.880 and the low-p̂ tertile at 0.891–0.998. Retained bias moves with it
in the same direction. Small families are where the correction over-shoots: at `K = 1` the retained
bias turns **negative** at every cell that has such rows (−0.0038 to −0.2044) and coverage rises
above nominal, on 5–43 replicates.

**Joint count table, K stratum × p̂ bin** — the table that says whether the two stratifications are
measuring the same thing. They are **strongly anti-diagonal**, not diagonal: small families carry
high p̂ and large families carry low p̂.

| cell | K T1: p T1/T2/T3 | K T2: p T1/T2/T3 | K T3: p T1/T2/T3 | K=1 *(ov)* | K≤5 *(ov)* |
|---|---|---|---|---|---|
| 12.4% n 500  | 43 / 142 / 291 | 170 / 180 / 126 | 264 / 152 / 59 | 24 | 82 |
| 12.4% n 1000 | 22 / 92 / 243 | 103 / 157 / 82 | 225 / 100 / 25 | 32 | 105 |
| 12.4% n 1500 | 13 / 51 / 167 | 67 / 109 / 53 | 149 / 69 / 9 | 43 | 113 |
| 31% n 500    | 106 / 194 / 320 | 235 / 225 / 160 | 281 / 199 / 140 | 7 | 36 |
| 31% n 1000   | 85 / 165 / 369 | 241 / 238 / 138 | 292 / 214 / 110 | 5 | 27 |
| 31% n 1500   | 60 / 142 / 391 | 226 / 246 / 119 | 309 / 201 / 82 | 9 | 45 |

### Beside the FS comparator — absolute levels, criterion stated

**At 31% the FS comparator is criterion-matched** (`cert20`, `effMaxSG`, ε 0.20 — DINA's own
criterion exactly). **At 12.4% it is not** (`tier2` at n 500, `p12ext` at n 1000/1500, both
`maxeffCons` at ε 0.10 — a different selection functional at half the band width), so at 12.4% the
criterion is a fourth confounded difference on top of identifier, family construction and detection
set, and a gap there cannot be read as engine behaviour even in part.

| cell | matched | FS criterion | DINA field lower | FS field lower |
|---|---|---|---|---|
| 12.4% n 500  | no | tier2 / maxeffCons / 0.1 | 0.8746 [0.8564, 0.8907] | 0.9625 [0.9511, 0.9714] |
| 12.4% n 1000 | no | p12ext / maxeffCons / 0.1 | 0.8494 [0.8265, 0.8697] | 0.9204 [0.9045, 0.9338] |
| 12.4% n 1500 | no | p12ext / maxeffCons / 0.1 | 0.8253 [0.7951, 0.8519] | 0.9303 [0.9148, 0.9431] |
| **31% n 500**  | **yes** | cert20 / effMaxSG / 0.2 | 0.9059 [0.8918, 0.9184] | 0.9734 [0.9650, 0.9798] |
| **31% n 1000** | **yes** | cert20 / effMaxSG / 0.2 | 0.9330 [0.9207, 0.9436] | 0.9560 [0.9458, 0.9643] |
| **31% n 1500** | **yes** | cert20 / effMaxSG / 0.2 | 0.9471 [0.9357, 0.9566] | 0.9666 [0.9576, 0.9738] |

| cell | DINA field-s upper | FS field-s upper | DINA joint-s | FS joint-s | DINA IJ | FS IJ |
|---|---|---|---|---|---|---|
| 12.4% n 500  | 0.9299 [0.9155, 0.9420] | 0.9398 [0.9258, 0.9512] | 0.9075 | 0.9566 | 0.9916 | 0.9860 |
| 12.4% n 1000 | 0.9466 [0.9313, 0.9587] | 0.9560 [0.9436, 0.9658] | 0.8875 | 0.9416 | 0.9886 | 0.9803 |
| 12.4% n 1500 | 0.9534 [0.9350, 0.9668] | 0.9495 [0.9359, 0.9603] | 0.8821 | 0.9431 | 0.9898 | 0.9824 |
| **31% n 500**  | 0.9027 [0.8884, 0.9153] | 0.9115 [0.8976, 0.9236] | 0.9027 | 0.9430 | 0.9919 | 0.9946 |
| **31% n 1000** | 0.9320 [0.9196, 0.9426] | 0.9277 [0.9152, 0.9385] | 0.9325 | 0.9371 | 0.9941 | 0.9948 |
| **31% n 1500** | 0.9358 [0.9234, 0.9463] | 0.9270 [0.9145, 0.9378] | 0.9414 | 0.9526 | 0.9966 | 0.9964 |

At the **criterion-matched** 31% cells the field lower bound is 6.8, 2.3 and 2.0 points below FS's
and closing with n; the field-s upper bound is within ±0.9 points and crosses over at n 1000; the
IJ interval is indistinguishable. At 12.4% the field gap is 8.8, 7.1 and 10.5 points, but the
criterion differs there, so it is not a like-for-like reading.

### Amendment 3 — same draws, all six cells

`n_true` `identical()` on all 2,000 rows: **YES at every cell.** `truth` `all.equal(tol = 1e-8)`:
**YES at every cell.** Maximum absolute discrepancy **0** (12.4% n 500, where `identical()` is also
TRUE), **2.22e-16** (12.4% n 1000 and n 1500) and **1.443e-15** (all three 31% cells); maximum
relative discrepancy at most **2.199e-15**. These are the cross-machine BLAS differences the
checker documents, not mismatches, and no cell is treated as failing on their account.

## Part A walls — realized against projection

The driver's own per-render and per-cell walls; no reconstruction is used.

| cell | driver wall | compute | batch renders | combine | overhead/batch | projected | realized/projected |
|---|---|---|---|---|---|---|---|
| B HR 1.75 n 1500, 31% | 7190 s | 6830.0 | 7181 | 9 | 175.5 | 1.9780 h | **1.010** |
| C 12.4% n 500  | 1183 s | 923.5  | 1174 | 9 | 125.3 | 0.3462 h | 0.949 |
| C 12.4% n 1000 | 712 s  | 521.6  | 703  | 9 | 90.7  | 0.2690 h | 0.735 |
| C 12.4% n 1500 | 445 s  | 296.8  | 436  | 9 | 69.6  | 0.1920 h | 0.644 |
| C 31% n 500    | 2710 s | 2238.0 | 2701 | 9 | 231.5 | 0.7084 h | 1.063 |
| C 31% n 1000   | 2290 s | 1922.2 | 2281 | 9 | 179.4 | 0.6716 h | 0.947 |
| C 31% n 1500   | 1876 s | 1573.2 | 1867 | 9 | 146.9 | 0.6347 h | 0.821 |

**Part A total: 4.557 h realized against 4.800 h projected — ratio 0.949 — inside a 9 h ceiling,
with all seven cells run and none deferred or dropped.** Start 20:02, finish 00:36.

The overhead model's **shape** is corrected by these numbers, though its total was right: the
combine render is a flat **9 s** on every cell and the whole of the overhead sits in the two batch
renders (the `n_super = 1e5` DGM build), at 69.6–231.5 s each, scaling with prevalence and
inversely with n. `projectC.R` charged 3 × a uniform per-render figure. `partA_accounting.R`
(committed) records this and reconciles exactly: compute + batch overhead + 9 = the driver wall,
to the second, at every cell.

## Stage 3 — the summary

`summary_dinamr.qmd` rendered clean: **"Cells on disk: 18 of 18"**, no chunk skipped, no deferred
cell. The absent-cell guard now skips nothing. The three new sections — `@sec-null-products`,
`@sec-null-location`, `@sec-null-fs` — are present, and the Block C tables carry `region_*` /
`compl_*` column names with the differentially-null statement in every caption.

---

# PART B — GRF cost probes

## GATE 3 — alignment: PASS

Both knobs are template **literals** with no `.env_*` read anywhere in the document, so they
resolve to whatever the source says; `grfprobe.R` reads them back out of the source rather than
asserting them.

| knob | resolves to | template line |
|---|---|---|
| `grf_selection` | `"frontier"` | 503 |
| `grf_select_statistic` | `"effect"` | 504 |
| `dmin.grf` | `0.0` | 506 |

Environment overrides for either knob: **none**.

## `dmin.grf = 0.0` — the recorded rationale, and what tracing the path adds

**The decision (Larry, 2026-09-11), recorded as given.** GRF's DR-scores target RMST for survival
outcomes, whereas FS and DINA both target the Cox hazard ratio. FS's and DINA's floors are
alignable with each other; GRF's is not alignable with either, so there is no GRF value that
reproduces DINA's sub-null log(0.90) floor. 0.0 is the null point on GRF's own scale and
reproduces the setting used in the manuscript's GRF runs.

**Consequence for every GRF record, unchanged:** GRF must **not** be called "FS-analogous". A
GRF-to-FS or GRF-to-DINA comparison differs in identifier, family construction, detection set,
selection criterion **and the scale of the selection criterion**. Unchanged also: the inference
products are computed on β(Ĥ) via the Cox model on the identified region whichever identifier
proposed it, so the estimand is the same kind of object across all three.

**What the path actually does — recorded so the rationale is not left unqualified.** Under this
configuration there are **two floors, on two scales**, and `dmin.grf` is only the first.

1. **`dmin.grf = 0.0` is a DR-score PRE-FILTER on the eligible set.** It is consumed only by the
   *native* frontier select inside the identifier — `R/grf_main.R:291` (survival) and
   `R/grf_subg_harm_glm.R:523` (GLM), both passing `dmin = dmin.grf` into
   `.grf_frontier_select()`, whose eligibility test is
   `elig <- cand[cand$effect >= dmin, , drop = FALSE]` at **`R/grf_subgroup_labels.R:358`**. On
   that call `effect` is the mean DR-score contrast, in RMST units. **The decision's premise is
   exactly right for this filter**: it is on GRF's own scale and not alignable with a log-HR floor.

2. **The binding effect-scale floor on the re-selection path is `hr.threshold = 0.90` — the same
   floor DINA carries, not `dmin.grf`.** With `grf_select_statistic = "effect"` and
   `grf_selection = "frontier"`, `.grf_reselect_on_effect()` re-scores the DR-candidate family on
   the Cox effect MR de-biases and re-selects on that. Its floor comes from the resolved admission
   set, not from `dmin.grf`:
   - **`R/forestsearch_helpers.R:1632–1635`** — `floor_cmp <- admission$effect_floor`, then
     `dmin_eff <- ... } else if (log_scale) exp(floor_cmp) else floor_cmp`.
   - **`R/forestsearch_main.R:2026–2030`** — the admission set is built once from
     `hr.threshold = threshold_config$screening`, which the comment at
     **`R/forestsearch_main.R:2020–2021`** states is *"already on the comparison scale (log for
     ratio measures)"*.
   - **`R/forestsearch_helpers.R:2345`** — `.fs_resolve_admission()` stores it verbatim as
     `effect_floor <- as.numeric(hr.threshold)`.
   - **`R/forestsearch_main.R:2442`** — it reaches the GRF selector as `admission = admission_resolved`.
   - **template line 531** — `hr_threshold <- 0.90`.

   So `dmin_eff = exp(log 0.90) = 0.90`.

**What this does to the recorded claim.** The claim that GRF's floor is **not alignable** with
DINA's is **true of `dmin.grf`, and only of `dmin.grf`**. It is **not** true of the floor that
decides the final selection here: on the effect re-selection path GRF and DINA apply the **same
numeric floor, HR ≥ 0.90**, from the same `hr.threshold`, because that floor is a property of the
admission set rather than of the engine. `dmin.grf = 0.0` therefore does **not** leave GRF
"unfloored" relative to DINA — it makes the DR pre-filter maximally permissive (every candidate
with a non-negative DR contrast survives) and leaves the binding decision to the shared 0.90.
**The decision stands as made; what changes is only what it is a decision about.** Nothing is
acted on here.

## The frontier band cannot come back empty under this configuration

`.compute_inclusion_band()` applies `hr_floor <- (1 - effect_neighborhood) * hr_max` then
`hr_vec >= hr_floor` (**`R/subgroup_consistency_helpers.R:784–785`**). The documented way the band
empties is recorded in-source at **`R/grf_subgroup_labels.R:377–383`**:

> "No empty-band fallback here, deliberately. The band CAN empty: when the maximum effect is
> negative, `(1 - nbhd) * emax` exceeds `emax`, so even the maximum fails its own test."

That requires a **negative maximum over the eligible set**, and the eligible set is
`{effect >= dmin}` (`R/grf_subgroup_labels.R:358`):

- on the **native DR path**, `dmin = dmin.grf = 0.0`, so every eligible effect is ≥ 0 and the
  maximum cannot be negative;
- on the **effect re-selection path** the scored column is a **hazard ratio**, `exp(beta_hat)`,
  strictly positive whatever the floor.

In both applications the maximum passes its own test, so **the band is non-empty whenever the
eligible set is**. What can empty is the **floor**, and on the re-selection path that is recorded
explicitly: `.grf_reselect_on_effect()` sets `grf_res$admitted_n <- 0L` and `sg_def <- NULL`
(**`R/forestsearch_helpers.R:1642–1650`**).

**This bears on Larry's open decision about the frontier-filter asymmetry** — GRF applies the band
frontier-only with no empty-band fallback where DINA uses it as a sort key, and MR's `.inband()`
carries a "never empty" fallback that GRF does not. The finding is that **at `dmin.grf = 0.0` the
asymmetry has no reachable consequence, because the branch it protects cannot be entered.**
**That decision is not resolved here** and nothing is changed on account of it; the frequency is
measured and reported, and the decision remains open.

## Measuring it

`grf_mechanism.R` separates the three no-selection mechanisms, which the bundle cannot: the
template returns its all-NA `NO-DETECTION` row at
`if (!found) { rec$status <- "NO-DETECTION"; return(rec) }` **before** `n_family` is written, so
all three collapse into one indistinguishable record. It re-runs the identifier on every replicate
of each corner and recomputes each cardinality the way `.grf_frontier_select()` does. On a
three-replicate check its recomputed eligible count reproduced the code's own `grf_res$admitted_n`
exactly (**50, 36, 114**), with band cardinalities 3, 11, 1 — never zero — at ~2.2 s per
identification.

## The cost surface

Five 36-replicate probes, 12 workers, tag `grfprobe`. **Whole of Part B ran in 5.7 minutes
(00:36:27 → 00:42:11), against a 1.5 h cap** — the cap never came close to binding, and Part A had
left 7.2 h under the 12 h timeout in any case.

| prevalence | HR | n | wall/replicate: med, p90, max (s) | family: med, q90, max | selection [Wilson 95%] | peak tree RSS |
|---|---|---|---|---|---|---|
| 12.4% | 1.50 | 500  | 14.18, 15.64, 16.70 | 776, 784, 853 | 36/36 = 1.0000 [0.9036, 1.0000] | 16.2 GB |
| 12.4% | 1.50 | 1500 | 16.96, 18.96, 21.48 | 828, 834, 838 | 36/36 = 1.0000 [0.9036, 1.0000] | 19.2 GB |
| 31%   | 1.50 | 500  | 16.00, 17.23, 18.33 | 776, 784, 853 | 36/36 = 1.0000 [0.9036, 1.0000] | 15.8 GB |
| 31%   | 1.50 | 1500 | 19.79, 22.10, 22.87 | 828, 834, 838 | 36/36 = 1.0000 [0.9036, 1.0000] | 19.7 GB |
| 12.4% | 1.00 | 500  | 13.50, 15.27, 16.35 | 776, 784, 853 | 35/36 = 0.9722 [0.8583, 0.9951] | 16.0 GB |

Per-render walls 63, 73, 67, 79, 62 s.

**Wall against family size: essentially flat.** Spearman ρ over the 36 replicates is 0.062,
−0.113, 0.084, −0.171, 0.123 — none of them meaningful at this n. Split at the family median, the
median wall differs by at most 0.2 s (e.g. 14.17 s vs 14.31 s at the first corner). GRF's cost is
driven by the forest fit and the per-candidate Cox re-scoring of an almost-constant family, not by
family size, which is the opposite of DINA's ~100× right-skewed cost profile.

**The proposed family is nearly constant, and identical across prevalences at matched n.** 776 /
784 / 853 (med / q90 / max) at n 500 at *both* prevalences, and 828 / 834 / 838 at n 1500 at both.
That is a consequence of how the family is built: `.grf_dr_candidates(X, dr_scores, n_min)`
enumerates threshold and pair candidates **on the covariate matrix subject to `n.min`**, so the
candidate *set* is a function of `(X, n_min)` and not of the outcome. The two prevalence blocks
share the same GBSG covariates and the same trial seeds and differ only in `k_inter`, i.e. only in
the outcome — hence identical families. For contrast, on the same design DINA's Block C family
runs median 36–485 with CV 1.008–1.458.

**Peak memory.** Two measures, both reported because they mean different things. `/usr/bin/time -l`
maximum resident set size covers the **parent quarto process only**: 1721–1877 MB. The 5-second
sampler sums RSS across the **whole quarto/R process tree**: median 12.7–18.7 GB, peak
**15.8–19.7 GB** against 36 GB physical at 12 workers. The summed figure **over-counts shared
pages** — forked workers share the parent's pages — so the true footprint is lower than 19.7 GB;
it is the right number for judging worker headroom, not for judging absolute footprint.

**Recorder columns on the GRF path: all present and populated.** At every one of the five probes,
p̂ (`p_hat_H`, `p_hat_sum`, `p_hat_top1`) present and populated, ρᶜ (`fld_Hc_scale_ratio`) present
and populated, and **9 of 9 recovery columns present, 9 of 9 populated**.

## The empty-band question, measured

`grf_mechanism.R` re-ran the identifier on **all 36 replicates of all five corners (180 in total)**
and recomputed each filter cardinality the way `.grf_frontier_select()` does.

| corner | band empty | admission empty | band size where admitted: min, med, max | DR pool: min, med, max | admitted (HR ≥ 0.90): min, med, max |
|---|---|---|---|---|---|
| 12.4% HR 1.50 n 500  | **0 / 36** | 0 | 1, 5, 23 | 760, 776, 789 | 4, 92, 247 |
| 12.4% HR 1.50 n 1500 | **0 / 36** | 0 | 1, 5.5, 17 | 946, 950.5, 955 | 18, 93.5, 383 |
| 31% HR 1.50 n 500    | **0 / 36** | 0 | 1, 5, 18 | 760, 776, 789 | 63, 355, 601 |
| 31% HR 1.50 n 1500   | **0 / 36** | 0 | 1, 6.5, 36 | 946, 950.5, 955 | 103, 464.5, 792 |
| 12.4% HR 1.00 n 500  | **0 / 36** | 0 | 1, 4, 18 | 760, 776, 789 | 3, 45.5, 190 |

**FRONTIER BAND EMPTY: 0 of 180 replicates (0.0000).** The source prediction holds exactly. The
band size is **never zero and frequently one** — its minimum is 1 at every corner, which is the
signature of the mechanism: the maximum of a non-negative eligible set always passes its own
`(1 − 0.20) × max` test, so the band retains at least the maximiser and can retain only that.

**The diagnostic reproduces the code's own bookkeeping exactly**: the recomputed eligible count
equals `grf_res$admitted_n` on **36 of 36 rows at every corner — 180 of 180**.

Identification-only re-run wall: median 1.35–2.26 s per replicate.

### What a replicate records when nothing is selected

The template returns its all-NA record at
`if (!found) { rec$status <- "NO-DETECTION"; return(rec) }`, **before** `n_family` is written, so
`n_family` is NA on every no-selection row — never 0, never positive — and an empty DR pool cannot
be told from an empty admission set in the committed columns. **A `forestsearch` error is
distinct**, and the recorder does separate it: `if (is.null(fs.est)) { rec$status <- "CONFIG-ERROR"; return(rec) }`
one branch earlier.

**No error, and no malformed record, on any of the 180 probe replicates.** Every row is either
`DETECTED` or a well-formed all-NA `NO-DETECTION`; there were **zero** `CONFIG-ERROR` rows. There
is nothing here to stop on.

### One pipeline no-selection that the diagnostic could not classify — stated, not resolved

The probes recorded exactly **one** no-selection in 180 replicates: sim_id 21 at
(12.4%, HR 1.00, n 500), `status = NO-DETECTION`, `mr_ok = 0`, `n_family` NA, `sg_def` NA.
**Its status is `NO-DETECTION`, not `CONFIG-ERROR`, so `forestsearch()` returned normally with no
subgroup — nothing errored.**

The diagnostic cannot assign it a mechanism, because **its fits are not bit-identical to the
pipeline's**. On the diagnostic's own fit that replicate selects cleanly (DR pool 771, admitted 68,
band 8, max HR 1.740). Checked while running this down:

- GRF **is** reproducible within a session — three consecutive standalone runs of the same
  replicate returned identical candidate counts, `admitted_n` and selected rule (772 / 50 /
  `{meno <= 0} & {pgr > 61.8}`).
- `seedit` **does** reach the forest: `fit_causal_forest()` passes `seed = seedit` to
  `grf::causal_survival_forest()` (`R/grf_helpers.R:74-85`).
- Yet across the two contexts the per-replicate families differ on **every** detected row —
  0 of 179 match — by up to 75 (n 500) and 131 (n 1500), while the *distributions* are nearly
  identical (median 776 in both at the null corner; pipeline range 729–853 against the
  diagnostic's 760–789).

So the two contexts produce statistically equivalent but not identical GRF fits, and the single
pipeline no-selection sits inside that difference. Candidate causes not chased further: the
`future` worker context the template's replicate loop runs in, and how `n.min = NULL` resolves.
**Chasing it further would risk an `R/` change, which this task forbids, and it is outside Part B's
cost-and-mechanism scope, so it is recorded here and left.**

**What this does and does not affect.** It does not affect the band finding, which is a structural
property read off the source — a non-negative eligible set cannot have a negative maximum — and
which held on 180 of 180 replicates regardless of which context produced the fit. It does mean the
mechanism split is measured on the diagnostic's fits, and that the pipeline's own single
no-selection has no mechanism assigned to it.

### Are GRF's candidate LISTS identical across prevalences, or only their sizes?

**The lists are not stored per replicate.** The probe bundles hold 167 columns, all atomic scalars;
none is a list, `AsIs` or matrix column. The only candidate-related fields are `n_family` (the
*size*), the **selected** rule (`sg_def`, `label`, `covs`) and `p_hat_top_labels` (the top three
by p̂). **A symmetric difference therefore cannot be computed from the committed bundles, and no
recorder change was made to obtain one.**

What the bundles do support is a strong necessary condition, and it holds exactly:

| matched pair | `sim_id` aligned | `n_true` identical | **`n_family` identical on all 36** | max abs difference | selected rule identical | `p_hat_top_labels` identical |
|---|---|---|---|---|---|---|
| n 500, 12.4% vs 31%  | yes | **no** | **yes — 36 of 36** | **0** | no — 8 of 36 | no — 0 of 36 |
| n 1500, 12.4% vs 31% | yes | **no** | **yes — 36 of 36** | **0** | no — 3 of 36 | no — 1 of 36 |

`n_true` differs, so the two blocks really are different data-generating mechanisms — same GBSG
covariates and same trial seeds, different `k_inter`, therefore different outcomes and different
planted regions. Against that, the family size is identical **replicate by replicate**, not merely
in distribution, while the selected rule almost always differs:

| sim | K (12.4% / 31%) | selected at 12.4% | selected at 31% |
|---|---|---|---|
| 1 (n 500) | 779 / 779 | `{age > 45} & {meno <= 0}` | `{age > 45} & {meno <= 0}` |
| 2 (n 500) | 776 / 776 | `{er > 162}` | `{er <= 52.4} & {meno <= 0}` |
| 3 (n 500) | 784 / 784 | `{age <= 45} & {size <= 38.4}` | `{meno <= 0} & {nodes <= 2}` |
| 1 (n 1500) | 830 / 830 | `{er <= 5} & {nodes > 1}` | `{er <= 60} & {meno <= 0}` |
| 2 (n 1500) | 827 / 827 | `{pgr <= 19} & {nodes <= 4}` | `{er <= 58} & {size > 19}` |

**The source settles what the sizes can only be consistent with.** `.grf_dr_candidates()`
(`R/grf_subgroup_labels.R:255-277`) builds each candidate from `X` alone: the cut points are
`stats::quantile(xj, probs = grid_probs)` per covariate column, and admission is
`if (nS < n_min || nS > n - 1L) next`. The DR scores enter **only** the `effect` column,
`mean(ctrl[S]) - mean(trt[S])`, computed *after* the candidate has been admitted — they never
decide which candidates exist. `.grf_dr_candidates_d2()` (`:281` onward) has the same shape. And
`n_min` is outcome-free too: with `n.min = NULL` it resolves to
`max(60L, ceiling(n.min.frac * N_analysis))` (`R/forestsearch_main.R:1368-1377`), a function of the
analysis sample size only.

**Finding: GRF's proposed family does not depend on the outcome.** The candidate *list* is a
deterministic function of `(X, grid_probs, n_min)`; only the per-candidate `effect` attached to it,
and hence the selection, is outcome-driven. That is established from source and is consistent with
identical sizes on 36 of 36 replicates in both matched pairs; it is **not** established by
enumerating the sets, which the bundles do not permit.

**How this bears on the family caveat.** DINA's family is read off a cross-fit surface a bootstrap
would regenerate, which is why the fixed-family condition fails and why every DINA number in this
record is conditional-on-proposed-family. GRF's family is **not** that kind of object: it is fixed
by the covariates and the sample size before any outcome is seen, so on this axis — and only this
axis — GRF is closer to the fixed-family condition than DINA is. Measured against DINA's Block C
families (median 36–485, CV 1.008–1.458, minimum 1) GRF's are 776–830 with the size pinned exactly
across two different DGMs. **This is a characterisation of the family, not a coverage claim, not a
comparison of products, and not a recommendation** — GRF's selection is still outcome-driven, and
everything else separating GRF from FS and DINA (identifier, detection set, selection criterion,
and the scale of that criterion) is unchanged.

### Bearing on the open decision

The frontier-filter asymmetry — GRF applying the band frontier-only with no empty-band fallback
where DINA uses it as a sort key, and MR's `.inband()` carrying a "never empty" fallback that GRF
does not — **has no reachable consequence at `dmin.grf = 0.0`**, because the branch it protects
cannot be entered. That is the evidence, and **the decision is left open.** Nothing was changed on
account of it, and no recommendation follows from it.

**No coverage table, no FS or DINA comparison, no acceptance criterion, no recommendation** appears
anywhere in Part B.
