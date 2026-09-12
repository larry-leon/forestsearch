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

One side issue, **not fixed**: the committed `render.sh` reads `$SP` for its log directory but
`campaign.sh` sets `SP` without exporting it, so a clean-environment run would `mkdir -p /logs`
and fail under `set -e`. The Part A run exported `SP` explicitly rather than editing either
script.

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

## Cells 2–7 — Block C

*(filled in as they land)*

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

*(probe results filled in below)*
