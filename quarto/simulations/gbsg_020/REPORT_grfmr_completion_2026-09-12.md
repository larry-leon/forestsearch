# REPORT — `grfmr` completion: the two deferred harm cells and the six HR 1.00 cells

- **Date:** 2026-09-12
- **Task:** `dev/tasks/TASK_grfmr_completion_2026-09-12.md` (committed as the first action, `6385b42a`)
- **Predecessors:** `REPORT_grfmr_2026-09-11.md`, `REPORT_dinamr_blockC_2026-09-11.md`,
  `REPORT_fs_extraction_2026-09-11.md`, `REPORT_fs_products_reconciliation_2026-09-12.md`
- **Machine:** Mac Studio (`Mac-Studio-3.local`, 14 cores, R 4.5.2, forestsearch 0.3.5). Unattended.
  **Commit only; not pushed.** No `R/` change and no template change.
- **Outcome:** all eight cells completed at 2,000 replicates. **`grfmr` is 18 of 18**, matching the
  DINA grid cell for cell. Nothing deferred, nothing dropped, no replicate count reduced.

**No acceptance criterion and no recommendation appear anywhere in this report.**

## Settled already, and cited rather than re-determined

- **HR 1.00 is not a global null.** The planted region carries HR 1.00 against a benefiting
  complement (HR 0.657 at 12.4%, 0.721 at 31%). This was established from source in the `dinamr`
  Block C Gate 0 on this same template (`REPORT_dinamr_blockC_2026-09-11.md`, "GATE 0"). These
  cells are **differentially null against a benefiting complement**; no harm vocabulary is used for
  them.
- **The rate at these cells is `selection_rate`.** The planted region clears the deliberately
  sub-null log(0.90) floor, so returning it is an admissible selection.
- **Nothing in the classification set is structurally undefined.** sens/spec/PPV/NPV and the
  oracle target (log 1.00 = 0) are defined. The structural columns remain `n_cons_qual`, `band_n`
  and `p_star`, as on every GRF cell.
- **Stratify on `admitted_n`, never on `n_family`**, which on GRF is the outcome-independent
  enumerated pool.

---

# STAGE 1 — 5-replicate smoke at (12.4%, HR 1.00, n 500): PASS

Driver `scripts_dinamr/grfmrC_smoke.sh` (campaign tag `grfmrsmk`), checker
`scripts_dinamr/stage1G.R`. Render wall 26 s. **30 passes, 0 failures.**

| check | result |
|---|---|
| meta: `grf`, HR 1.00, `effMaxSG` ε 0.20, field knobs, `two_term` | PASS |
| 5 rows, no CONFIG-ERROR | PASS |
| `admitted_n` present **and populated on the null path** | PASS — finite on 5 of 5 |
| p̂ block, ρᶜ / scale block, nine recovery columns present and populated | PASS |
| every construction finite on detected replicates | PASS — 45 quantities |
| sens/spec/PPV/NPV finite (classification set defined) | PASS |
| 13 interval / validity invariants | PASS |
| γ (joint, joint-s) in [0.025, 0.05] | PASS — [0.02500, 0.02600] |
| corrected identity `log(est2_s) + lam_mean_s == log(est2) + lam_mean` | PASS — 1.67e-16 |
| realized prevalence | PASS — super-population 0.12418; trial mean 0.1340 (five draws, 0.114–0.152) |

- **Selection rate: 5 / 5.**
- **`admitted_n`:** 39, 114, 24, 41, 48 — min 24, q25 39, **median 41**, q75 48, p90 87.6, max 114.
  `n_family` (the enumerated pool, not the qualified set) 729–784.
- **Non-detections: 0.** So the question "NA or 0 on a non-detection" does not arise on the smoke;
  it is answered on the full cells below.
- Structural: `n_cons_qual`, `band_n` present and all-NA; `p_star` not a recorder column.
- Gate 3 on the smoke batch: 7 passes, 0 failures.

---

# GATE 1 — compute go/no-go: GO for all eight

- **Ceiling 9 h wall. Hard timeout 10 h.** Both from the kickoff; the timeout was enforced by a
  36,000 s watchdog in the driver.
- Projector `scripts_dinamr/projectGC.R`: bootstrap (B = 4,000) over the probes' per-replicate
  `fit_mr_secs`, plus `projectC.R`'s measured overhead **per batch render** (122.8 s at 12.4%,
  196.6 s at 31%, three renders per cell), **calibrated ×0.8962** on `grfmr` Part A's
  realized-over-projected ratio.
- Cost basis, stated per cell: 12.4% HR 1.00 n 500 is **measured** (the null-corner probe,
  13.50 s median). Every other HR 1.00 cell and both HR 1.75 cells use the HR 1.50 pool at their
  own prevalence and n (n 1000 = n 500 and n 1500 draws pooled). That is conservative for the null
  cells: at the one coordinate with both probes, null runs at 0.952 of harm. A null-scaled
  alternative is printed beside it and not used.

| order | cell | basis | projected h (calibrated) | uncalibrated | cumulative |
|---|---|---|---|---|---|
| 1 | 31% HR 1.75 n 1500 | HR 1.50 pool | 0.9758 | 1.0888 | 0.976 |
| 2 | 31% HR 1.75 n 1000 | HR 1.50 pool, pooled n | 0.8941 | 0.9977 | 1.870 |
| 3 | 12.4% HR 1.00 n 500 | **measured** | 0.6477 | 0.7227 | 2.518 |
| 4 | 12.4% HR 1.00 n 1000 | HR 1.50 pool, pooled n | 0.7406 | 0.8264 | 3.258 |
| 5 | 12.4% HR 1.00 n 1500 | HR 1.50 pool | 0.8037 | 0.8968 | 4.062 |
| 6 | 31% HR 1.00 n 500 | HR 1.50 pool | 0.8125 | 0.9066 | 4.875 |
| 7 | 31% HR 1.00 n 1000 | HR 1.50 pool, pooled n | 0.8941 | 0.9977 | 5.769 |
| 8 | 31% HR 1.00 n 1500 | HR 1.50 pool | 0.9758 | 1.0888 | **6.744** |

- **Eight cells: 6.744 h calibrated** (90% band 6.723–6.765 h); 7.526 h uncalibrated. That is
  1.870 h for the two harm cells and 4.874 h for the six null cells, against the kickoff's
  reference of about 1.9, 4.9 and 6.8 h.
- **GATE 1: GO — run 8, defer 0.** Even uncalibrated the eight fit under the ceiling with 1.47 h
  to spare.

---

# The run

Driver `scripts_dinamr/grfmrC.sh` with `scripts_dinamr/grfmrC.cells`. It is `grfmr.sh`'s knob set
unchanged — `FS_S7_METHOD=grf FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20`, the four field knobs,
`FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=grfmr FS_S7_WORKERS=12`, and
`FS_S7_ER_JCUTS` unset. Seeds are 8316951 + sim_id, run as two batches of 1,000 and combined,
with the three thread variables at 1. It adds two things: **Gate 3 after every batch,
stop-on-failure**, and the **10 h watchdog**.

**Started 14:06:42, finished 20:32:13. Wall 6.425 h (23,131 s).** The span equals the sum of the
per-cell walls. The watchdog never fired.

## Walls, realized against the Gate 1 projection

| order | cell | realized | projected (calibrated) | ratio | ratio to uncalibrated |
|---|---|---|---|---|---|
| 1 | 31% HR 1.75 n 1500 | 0.9744 h (3508 s) | 0.9758 h | 0.999 | 0.895 |
| 2 | 31% HR 1.75 n 1000 | 0.8983 h (3234 s) | 0.8941 h | 1.005 | 0.900 |
| 3 | 12.4% HR 1.00 n 500 | 0.6722 h (2420 s) | 0.6477 h | 1.038 | 0.930 |
| 4 | 12.4% HR 1.00 n 1000 | 0.7228 h (2602 s) | 0.7406 h | 0.976 | 0.875 |
| 5 | 12.4% HR 1.00 n 1500 | 0.7428 h (2674 s) | 0.8037 h | 0.924 | 0.828 |
| 6 | 31% HR 1.00 n 500 | 0.7369 h (2653 s) | 0.8125 h | 0.907 | 0.813 |
| 7 | 31% HR 1.00 n 1000 | 0.8075 h (2907 s) | 0.8941 h | 0.903 | 0.809 |
| 8 | 31% HR 1.00 n 1500 | 0.8703 h (3133 s) | 0.9758 h | 0.892 | 0.799 |
| | **TOTAL** | **6.425 h** | **6.744 h** | **0.953** | 0.854 |

- **The two harm cells landed on their projection:** 1.873 h against 1.870 h (ratio 1.002). The
  0.896 calibration, derived from Part A's harm cells, carried over exactly to the remaining
  harm cells.
- **The null cells ran 6.6% under** (4.553 h against 4.874 h, ratio 0.934). The one measured
  basis, 12.4% n 500, ran 3.8% over. The five cells costed from the HR 1.50 pool ran under by
  2.4–10.8%, more so at 31% and larger n. Stated plainly, that is the conservative basis doing
  what it was expected to do. The null probe's 0.952 null-to-harm ratio at n 500 understates the
  gap at larger n.
- Per-render walls are in `logs/grfmrC.driver.log` and `scripts_dinamr/wallsGC.R`
  (`wallsGC.rds`): 1,201–1,775 s per batch and 9–10 s per combine.

## Gate 3 — alignment, per batch: PASS on all sixteen batches

**7 passes, 0 failures on every batch** (`logs/gate3_<cell>_batch_<start>.log`).
`grf_select_statistic` resolved to **`"effect"`**, `grf_selection` to **`"frontier"`** and
`dmin.grf` to **`0`**. Resolution is from source (template lines 503–506), from the batch's own
audit line `method=grf/frontier`, and from the bundle (`admitted_n` finite, which only the
effect/frontier path writes). `dmin.grf` is resolved from source only, as before; it has no
read-back.

Batch-level `admitted_n` finiteness is the first sign of the null cells' non-detections: 987 / 983
(12.4% n 500), 963 / 954 (n 1000), 877 / 866 (n 1500), and 997–1000 at 31%.

## Gate 2 — per cell: PASS on all eight

`Rscript gate2G.R B` (the two new cells among the six 31% harm cells) and `Rscript gate2G.R C`
(the six HR 1.00 cells) each gave **222 passes, 0 failures**, i.e. 37 per cell. That is the
Part A set plus one new check: the **designated comparator resolves** — `tier2` / `p12ext` /
`e1stud` / `cert20` on disk, never a `map1` / `s7` bundle. Output: `logs/gate2G_B.txt`,
`logs/gate2G_C.txt`. The same check on the ten Part A cells gives 370 passes, 0 failures.

- **Completeness and meta:** 2,000 rows, `sim_id` 1..2000, no CONFIG-ERROR, two batches,
  `n_workers` 12, every knob as pinned, `campaign_tag` grfmr.
- **Finiteness:** all 38 gated products finite on every detected replicate. p̂, ρᶜ and the nine
  recovery columns are present and populated.
- **Invariants:** all eleven hold, including `admitted_n >= 1` on every detected replicate.
- **γ** in [0.02500, 0.02700] on both joints, every cell.
- **Corrected identity** `log(fld_Hc_est2_s) + fld_Hc_lam_mean_s == log(fld_Hc_est2) +
  fld_Hc_lam_mean`: **2.22e-16 to 3.33e-16**.
- **Bonferroni identity** on the γ-at-floor rows: **exactly 0**. The share at the floor is
  0.906–0.945 (joint) and 0.802–0.899 (joint-s).
- **Realized prevalence:** 12.4% cells 0.12363 / 0.12388 / 0.12401 against 0.12418; 31% cells
  0.30591–0.30646 against 0.30655.
- **Structural, reported as such:** `n_cons_qual`, `band_n` present and all-NA; `p_star` not a
  recorder column. Structural on every GRF cell, not on account of HR 1.00.

### Designated comparator, resolved per cell before the assertion

| cell | comparator | `sg_focus` / ε | criterion |
|---|---|---|---|
| 31% HR 1.75 n 1000, n 1500 | `cert20` | `effMaxSG` / 0.20 | **matched** |
| 12.4% HR 1.00 n 500 | `tier2` | `maxeffCons` / 0.10 | not matched |
| 12.4% HR 1.00 n 1000, n 1500 | `p12ext` | `maxeffCons` / 0.10 | not matched |
| 31% HR 1.00 n 500, 1000, 1500 | `cert20` | `effMaxSG` / 0.20 | **matched** |

`e1stud` has no HR 1.00 bundle, so 31% HR 1.00 n 500 resolves to `cert20`. No `map1` or `s7`
bundle was substituted.

### Amendment 3 — same draws: holds on all eight

**`n_true` `identical()` on all 2,000 rows of every cell. `truth` `all.equal()` at 1e-8: YES on
every cell.** No DGM-path finding.

| cells | `truth` `identical()` | max abs diff | max rel diff |
|---|---|---|---|
| 12.4% HR 1.00 n 500 (`tier2`) | **TRUE** | **0** | **0** |
| 12.4% HR 1.00 n 1000, n 1500 (`p12ext`) | FALSE | 2.220e-16 | 3.797e-16 |
| 31% HR 1.00, all three n (`cert20`) | FALSE | 1.443e-15 | 2.199e-15 |
| 31% HR 1.75 n 1000, n 1500 (`cert20`) | FALSE | 8.882e-15 | 4.253e-15 |

### Per cell, prominently: selection or detection, `admitted_n`, and the enumerated pool beside it

`n_family` is **the enumerated pool, not the qualified set**. It comes from
`.grf_dr_candidates()`'s quantile enumeration over X subject to `n_min` and does not depend on the
outcome. `admitted_n` is the forest-qualified count and does.

| cell | rate [Wilson] | **`admitted_n`** min / q25 / **median** / q75 / p90 / max | CV | `n_family` (enumerated pool) min / median / max |
|---|---|---|---|---|
| 31% HR 1.75 n 1000 | detection 1.0000 [0.9981, 1.0000] | 116 / 400.75 / **487.5** / 566.25 / 627.1 / 765 | 0.242 | 779 / 830 / 914 |
| 31% HR 1.75 n 1500 | detection 1.0000 [0.9981, 1.0000] | 152 / 409.75 / **493.5** / 562 / 607.1 / 761 | 0.217 | 780 / 830 / 913 |
| 12.4% HR 1.00 n 500 | **selection 0.9840** [0.9775, 0.9886] | 0 / 32 / **61.5** / 113 / 179.1 / 551 | 0.898 | 712 / 776 / 870 |
| 12.4% HR 1.00 n 1000 | **selection 0.9575** [0.9477, 0.9655] | 0 / 17 / **38** / 70 / 109 / 304 | 0.901 | 779 / 830 / 914 |
| 12.4% HR 1.00 n 1500 | **selection 0.8675** [0.8519, 0.8817] | 0 / 11 / **23** / 42 / 70 / 279 | 0.945 | 780 / 829 / 913 |
| 31% HR 1.00 n 500 | **selection 0.9970** [0.9935, 0.9986] | 0 / 95.5 / **166** / 281 / 397.6 / 721 | 0.675 | 712 / 776 / 870 |
| 31% HR 1.00 n 1000 | **selection 0.9980** [0.9949, 0.9992] | 3 / 82 / **141** / 225 / 328.5 / 701 | 0.680 | 779 / 830 / 914 |
| 31% HR 1.00 n 1500 | **selection 0.9935** [0.9889, 0.9962] | 3 / 66 / **115** / 176 / 250 / 649 | 0.684 | 780 / 830 / 913 |

The `admitted_n` quantiles are over every row where it is recorded, which includes the
`admitted_n = 0` non-detections below.

### Non-detections: what `admitted_n` records on them

| cell | non-detections | `admitted_n` NA | `admitted_n` = 0 | status |
|---|---|---|---|---|
| 12.4% HR 1.00 n 500 | 32 (0.0160) | 30 | **2** | NO-DETECTION, `err_msg` none |
| 12.4% HR 1.00 n 1000 | 85 (0.0425) | 83 | **2** | NO-DETECTION |
| 12.4% HR 1.00 n 1500 | 265 (0.1325) | 257 | **8** | NO-DETECTION |
| 31% HR 1.00 n 500 | 6 (0.0030) | 5 | **1** | NO-DETECTION |
| 31% HR 1.00 n 1000 | 4 (0.0020) | 4 | 0 | NO-DETECTION |
| 31% HR 1.00 n 1500 | 13 (0.0065) | 13 | 0 | NO-DETECTION |
| 31% HR 1.75 n 1000, n 1500 | none | — | — | — |

- **Every non-detection is `NO-DETECTION`, never CONFIG-ERROR, with no error message.** They are
  recorded, not chased, and no recorder change was made.
- **The empty-admitted-set path fired for the first time: 13 replicates carry `admitted_n = 0`**,
  all non-detections, 12 of them at 12.4%. It never fired across the 20,000 Part A replicates or
  the 180 probe replicates. Part T2's placement, which writes `admitted_n` before the no-detection
  return, is what makes these separable. The other **392** non-detections carry `admitted_n` NA:
  the re-selection returned before the admission count was written. The two causes are distinct
  in the columns, which DINA's recorder could not show.
- **`admitted_n = 1` now occurs on detected replicates** (5 / 7 / 20 at 12.4% n 500 / 1000 / 1500,
  1 at 31% n 500). It never occurred in Part A.

---

# The question: selection rate and bound location at the null cells

**Two quantities, reported separately.** The selection rate is over all 2,000 replicates. The
location shares are over the detected replicates, restricted to rows with the field lower bound
and θ(Ĥ) finite, the row set `fs_extraction.R` uses for FS. FS numbers are the designated
comparator read as it stands. They match `REPORT_fs_extraction_2026-09-11.md` §2 to every printed
digit.

| cell | engine | criterion | selection rate [Wilson] | n_eval | median lower bound | median θ(Ĥ) | bound − θ | bound / θ | paired ratio | planted | **share ≥ 1.00** [Wilson] | **share ≥ 1.25** [Wilson] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 12.4% n 500 | GRF | `grfmr` `effMaxSG` 0.20 | **0.9840** [0.9775, 0.9886] | 1968 | 0.3943 | 0.6765 | −0.2821 | 0.5830 | 0.5779 | 1.0005 | **0.0061** [0.0035, 0.0106] | **0.0025** [0.0011, 0.0059] |
| 12.4% n 500 | FS | `tier2` `maxeffCons` 0.10 — not matched | 0.6805 [0.6597, 0.7006] | 1361 | 0.4276 | 0.7053 | −0.2777 | 0.6062 | 0.6036 | 1.0005 | 0.0029 [0.0011, 0.0075] | 0.0000 [0.0000, 0.0028] |
| 12.4% n 1000 | GRF | `grfmr` `effMaxSG` 0.20 | **0.9575** [0.9477, 0.9655] | 1915 | 0.4547 | 0.7060 | −0.2513 | 0.6440 | 0.6541 | 1.0005 | **0.0042** [0.0021, 0.0082] | **0.0010** [0.0003, 0.0038] |
| 12.4% n 1000 | FS | `p12ext` `maxeffCons` 0.10 — not matched | 0.6595 [0.6384, 0.6799] | 1319 | 0.5116 | 0.7448 | −0.2332 | 0.6869 | 0.6933 | 1.0005 | 0.0076 [0.0041, 0.0139] | 0.0023 [0.0008, 0.0067] |
| 12.4% n 1500 | GRF | `grfmr` `effMaxSG` 0.20 | **0.8675** [0.8519, 0.8817] | 1735 | 0.5165 | 0.7286 | −0.2122 | 0.7088 | 0.7219 | 1.0005 | **0.0046** [0.0023, 0.0091] | **0.0017** [0.0006, 0.0051] |
| 12.4% n 1500 | FS | `p12ext` `maxeffCons` 0.10 — not matched | 0.6240 [0.6026, 0.6450] | 1248 | 0.5571 | 0.8029 | −0.2459 | 0.6938 | 0.7026 | 1.0005 | 0.0064 [0.0033, 0.0126] | 0.0008 [0.0001, 0.0045] |
| 31% n 500 | GRF | `grfmr` `effMaxSG` 0.20 | **0.9970** [0.9935, 0.9986] | 1994 | 0.4764 | 0.8148 | −0.3383 | 0.5848 | 0.5869 | 0.9999 | **0.0196** [0.0143, 0.0266] | **0.0065** [0.0038, 0.0111] |
| 31% n 500 | FS | `cert20` `effMaxSG` 0.20 — **matched** | 0.9205 [0.9078, 0.9316] | 1841 | 0.4591 | 0.8241 | −0.3650 | 0.5571 | 0.5637 | 0.9999 | 0.0060 [0.0033, 0.0107] | 0.0005 [0.0001, 0.0031] |
| 31% n 1000 | GRF | `grfmr` `effMaxSG` 0.20 | **0.9980** [0.9949, 0.9992] | 1996 | 0.5383 | 0.8492 | −0.3109 | 0.6339 | 0.6427 | 0.9999 | **0.0110** [0.0073, 0.0166] | **0.0025** [0.0011, 0.0059] |
| 31% n 1000 | FS | `cert20` `effMaxSG` 0.20 — **matched** | 0.9545 [0.9445, 0.9628] | 1909 | 0.5512 | 0.8494 | −0.2982 | 0.6489 | 0.6558 | 0.9999 | 0.0105 [0.0068, 0.0161] | 0.0016 [0.0005, 0.0046] |
| 31% n 1500 | GRF | `grfmr` `effMaxSG` 0.20 | **0.9935** [0.9889, 0.9962] | 1987 | 0.5943 | 0.8786 | −0.2843 | 0.6764 | 0.6902 | 0.9999 | **0.0091** [0.0057, 0.0143] | **0.0010** [0.0003, 0.0037] |
| 31% n 1500 | FS | `cert20` `effMaxSG` 0.20 — **matched** | 0.9590 [0.9494, 0.9668] | 1918 | 0.6032 | 0.8709 | −0.2677 | 0.6926 | 0.7058 | 0.9999 | 0.0094 [0.0059, 0.0148] | 0.0010 [0.0003, 0.0038] |

**Confound, with every comparison.** FS and GRF differ in identifier, in family construction and
in detection set. At 12.4% they also differ in the selection criterion (`maxeffCons` ε 0.10 against
`effMaxSG` ε 0.20), so a 12.4% gap cannot be read as engine behaviour even in part. At 31% the
criterion is matched, and the confound is identifier, family construction and detection set.

For reference, DINA on the same draws (`REPORT_dinamr_blockC_2026-09-11.md`, `effMaxSG` ε 0.20 at
both prevalences, with its own identifier, family construction and detection set):

| prevalence | DINA selection rate n 500 / 1000 / 1500 | DINA share ≥ 1.00 | DINA share ≥ 1.25 |
|---|---|---|---|
| 12.4% | 0.7135 / 0.5245 / 0.3435 | 0.0189 / 0.0162 / 0.0087 | 0.0035 / 0.0010 / 0.0015 |
| 31% | 0.9300 / 0.9260 / 0.8880 | 0.0301 / 0.0157 / 0.0096 | 0.0086 / 0.0032 / 0.0017 |

## What the null cells say, descriptively

- **GRF's selection rate at 31% is at 1 within its interval: 0.9970 / 0.9980 / 0.9935.** It does
  not move with n. At these three cells the selection rate carries essentially no information
  about the planted region being differentially null, and **bound location is the whole
  diagnostic.** FS on the same draws at the matched criterion selects 0.9205 / 0.9545 / 0.9590 and
  DINA 0.9300 / 0.9260 / 0.8880.
- **At 12.4% GRF's selection rate is high and falls with n: 0.9840 → 0.9575 → 0.8675.** It falls
  much less than FS's (0.6805 → 0.6595 → 0.6240, under a different criterion) or DINA's (0.7135 →
  0.5245 → 0.3435). The fall comes from non-detections, 32 → 85 → 265. These three are the only
  cells among GRF's 18 where the detection or selection rate moves with n by more than 0.01.
  Elsewhere the largest movement is 0.0045, at 31% HR 1.00. So a 12.4% null n-trend on GRF
  carries a detection-conditioning qualification that no other GRF trend does.
- **The location shares stay small at every null cell.** Share ≥ 1.00 is **0.0042–0.0196** and
  share ≥ 1.25 is **0.0010–0.0065**, across the six GRF null cells. At matched prevalence and n, the harm-cell ≥ 1.00 share
  is **3.8 to 90 times** the null share: from 12.4% HR 1.50 n 500 (0.0230 against 0.0061) to 31%
  HR 1.75 n 1500 (0.8170 against 0.0091). The selection rate at 31% does not separate the harm
  and null cells at all.
  - At 31% the ≥ 1.00 share **falls with n** (0.0196 → 0.0110 → 0.0091), as DINA's does.
  - At 12.4% it is **flat within its Wilson intervals** (0.0061 / 0.0042 / 0.0046).
- **Beside FS at the matched criterion (31%).** At n 500 GRF's ≥ 1.00 share, 0.0196
  [0.0143, 0.0266], sits above FS's 0.0060 [0.0033, 0.0107], and the intervals do not overlap.
  At n 1000 and 1500 the two are indistinguishable (0.0110 against 0.0105; 0.0091 against 0.0094).
  The ≥ 1.25 share follows the same pattern: 0.0065 against 0.0005 at n 500, then 0.0025 against
  0.0016 and 0.0010 against 0.0010. The confound above applies to every one of these.
- **Bound location against the realized target.** The median lower bound sits below median θ(Ĥ)
  at every null cell. The ratio is 0.583–0.709 on GRF and rises with n. The ratio and the paired
  ratio agree within 0.014, so the ordering is not a median artefact. Median θ(Ĥ) at the null is
  0.68–0.73 at 12.4% and 0.81–0.88 at 31%. GRF's returned region is therefore a mix of the planted
  region and benefiting complement, as the classification rates below say directly.

## Coverage and classification at the null, absolute levels

Over detected replicates. The H block at these cells is the identified region; the complement
benefits.

| cell | field lower on β(Ĥ) | field-s upper on β(Ĥᶜ) | IJ two-sided H | IJ two-sided Hᶜ | Bonferroni `joint_s` | sens / spec / PPV / NPV | mean \|Ĥ\| |
|---|---|---|---|---|---|---|---|
| 12.4% n 500 | 0.9395 [0.9281, 0.9492] | 0.9319 [0.9199, 0.9422] | 0.9924 | 1.0000 | 0.9405 [0.9292, 0.9502] | 0.356 / 0.810 / 0.213 / 0.900 | 105.2 |
| 12.4% n 1000 | 0.9473 [0.9363, 0.9564] | 0.9431 [0.9318, 0.9526] | 0.9911 | 1.0000 | 0.9410 [0.9295, 0.9507] | 0.434 / 0.830 / 0.275 / 0.913 | 202.9 |
| 12.4% n 1500 | 0.9568 [0.9462, 0.9654] | 0.9476 [0.9360, 0.9571] | 0.9925 | 1.0000 | 0.9522 [0.9411, 0.9612] | 0.526 / 0.834 / 0.324 / 0.926 | 316.7 |
| 31% n 500 | 0.9418 [0.9307, 0.9513] | 0.9107 [0.8974, 0.9225] | 0.9935 | 1.0000 | 0.9273 [0.9150, 0.9379] | 0.350 / 0.835 / 0.471 / 0.748 | 110.8 |
| 31% n 1000 | 0.9469 [0.9362, 0.9559] | 0.9259 [0.9135, 0.9365] | 0.9905 | 0.9995 | 0.9334 [0.9216, 0.9435] | 0.439 / 0.845 / 0.543 / 0.779 | 242.0 |
| 31% n 1500 | 0.9587 [0.9491, 0.9666] | 0.9361 [0.9245, 0.9460] | 0.9914 | 0.9995 | 0.9467 [0.9359, 0.9557] | 0.568 / 0.833 / 0.603 / 0.821 | 435.2 |

- Field lower on β(Ĥ) **rises with n at both prevalences** (0.9395 → 0.9568; 0.9418 → 0.9587), as
  on the harm cells. Field-s upper on β(Ĥᶜ) rises with n too, and sits lower at 31%
  (0.9107–0.9361) than at 12.4% (0.9319–0.9476).
- **FS products are quoted only with the cell set and construction named**, per
  `REPORT_fs_products_reconciliation_2026-09-12.md`. On the **six HR 1.00 cells**, for the
  **field one-sided lower bound on β(Ĥ)**, the designated FS comparators read 0.9625 / 0.9204 /
  0.9303 at 12.4% (`tier2` / `p12ext`, `maxeffCons` ε 0.10) and 0.9734 / 0.9560 / 0.9666 at 31%
  (`cert20`, `effMaxSG` ε 0.20). Row set: detected, with β(Ĥ), the field lower bound and the field-s
  upper bound finite (`grfmr_numbers.R`). The confound sentence above applies.
- Sensitivity and PPV rise with n; specificity is flat at 0.81–0.85. Mean |Ĥ| grows with n, to
  317 (12.4%) and 435 (31%) at n 1500.

## The `admitted_n` and p̂ strata at the null

From `TABLES_grfmr_completion_2026-09-12.md` §3–§4, extracted by re-executing
`summary_grfmr.qmd`'s own chunks.

- **By `admitted_n` tertile, field lower coverage falls from T1 to T3 at every null cell.** At
  12.4%: 0.9849 → 0.9491 → 0.8841 (n 500), 0.9786 → 0.9652 → 0.8968 (n 1000), 0.9817 → 0.9748
  → 0.9135 (n 1500). At 31%: 0.9851 → 0.9395 → 0.9006, 0.9791 → 0.9671 → 0.8935, 0.9881 →
  0.9589 → 0.9287. Retained bias rises from T1 to T3 (for instance −0.104 → +0.095 → +0.243 at
  12.4% n 500). The T1-to-T3 gap narrows from n 500 to n 1500 at both prevalences: 0.101 → 0.082
  → 0.068 at 12.4%, and 0.085 → 0.086 → 0.059 at 31%, so not monotone there. This is the same
  direction as the harm cells.
- **By p̂ tertile the gradient is sharper.** T1 and T2 cover at 0.994–1.000; T3 covers at
  0.826–0.878, with retained bias +0.113 to +0.369.
- **The joint counts are near-uniform at every null cell.** Base-tertile counts run 172–259
  against 193–222 expected under independence (n_eval / 9), as on the harm cells. On GRF, `admitted_n` and p̂ remain
  separate axes at the null.

---

# STAGE 3 — `summary_grfmr.qmd` extended to 18 cells

**Rendered `summary_grfmr.html` (8.2 MB), render wall 21 s, RC 0.** The render reads
**"Cells on disk: 18 of 18 (12 harm cells, 6 HR 1.00 cells)"**. No guarded chunk skipped: the only
"Skipped" and "NOT ON DISK" strings in the HTML are the echoed source's format strings. All 39 R
chunks parse.

Changes to the document (`0c2de3f3`):

- The six HR 1.00 cells are back on the cell list as `arm = "null"`, so the absent-cell guard now
  applies to them rather than to an out-of-scope statement. Subtitle, grid paragraph and inventory
  line updated.
- `sec-null` **cites** the `dinamr` Block C Gate 0 rather than repeating it. "Nothing in the
  classification set is structurally undefined" is stated; the structural columns are named as
  structural on every GRF cell.
- The detection table and standard-table captions state that the same column is a
  `selection_rate` at HR 1.00. "harm block" became "Hhat block", with the region / complement
  reading stated for the null cells.
- The `null` table carries the `admitted_n` distribution (median, q25, q75, p90) with the
  enumerated pool beside it. Its columns and record lines use region / complement vocabulary.
- **The location table, added** (`null-location`). One GRF row and one FS row per null cell, with
  campaign / `sg_focus` / ε and a `criterion_matched` flag. The **selection rate** sits in its own
  column with Wilson limits. Then median lower bound, median θ(Ĥ), field and naive estimates, the
  gap as difference / ratio / paired ratio, the planted target, and **shares ≥ 1.00 and ≥ 1.25
  with Wilson limits**. The caption states the products, the cells, that rates other than the
  selection rate are over detected replicates, the criterion match per prevalence, and the
  confound. The record chunk echoes the same, with the confound sentence.
- Carried through for the new cells by the existing chunks: the `admitted_n` and p̂ strata and
  their joint count table; the error-SD Gaussian reference beside the marginal one with both
  formulas stated; Wilson limits on every rate; marginal and error SD side by side; absolute
  coverage for every product.
- Two stale citations in the null section were corrected: `scripts_grfmr/blockA_rest.R` →
  `scripts_dinamr/blockA_rest.R`, and `TASK_grfmr_blockC_grfprobe_2026-09-11` →
  `TASK_dinamr_blockC_grfprobe_2026-09-11`. Neither path exists under the old name.

`TABLES_grfmr_completion_2026-09-12.md` is the 18-cell companion to the 10-cell
`TABLES_grfmr_percell_2026-09-12.md`, which is left as committed. The ten shared cells' numbers
are identical in both.

---

# Tooling added or changed (all in `scripts_dinamr/`)

| file | role |
|---|---|
| `grfmrC_smoke.sh`, `stage1G.R` | Stage 1 smoke driver and battery |
| `projectGC.R` → `projectionGC.rds` | Gate 1 projection, calibrated ×0.8962 |
| `grfmrC.sh`, `grfmrC.cells` | the eight-cell driver: Gate 3 per batch, stop-on-failure, 10 h watchdog |
| `wallsGC.R` → `wallsGC.rds` | realized walls against the projection, from the driver's own lines |
| `gate2G.R` | + block `C`; SELECTION RATE label at HR 1.00; designated-comparator check |
| `grfmr_numbers.R` → `grfmr_numbers.rds` | 18 cells; HR 1.00 comparator map; non-detection split; null location table with FS beside GRF |
| `grfmr_tables.R` → `grfmr_tables.rds` | unchanged; re-run over 18 cells |
| `status_inventory.R` | regenerates `current_status.md` §3 from the directory |
| `parse_chunks.R` | parses every R chunk of a `.qmd` (the 39-of-39 check above) |

---

# Side issues, flagged and not fixed

- **`summary_grfmr.qmd` transplant wording outside the sections this task extends.** Line 26 says
  the re-selection floor is "the same one GRF carries", where DINA is meant. Line 41 calls
  `admitted_n` "the analogue of GRF's family-size stratifier", and the strata section says "the
  family-size stratifier the GRF summary uses"; DINA's is meant in both. Left as
  committed under scope discipline.
- **The `strat-miss` caption** labels its first stratification "proposed-family-size" although the
  stratum is `admitted_n` (inherited from the transplant). The numbers are keyed on `admitted_n`;
  only the caption word is stale.
