# REPORT — GRF MR campaign (`grfmr`)

- **Date:** 2026-09-12 (task dated 2026-09-11)
- **Task:** `dev/tasks/TASK_grfmr_campaign_2026-09-11.md` (committed as the campaign's first action)
- **Predecessors:** `dev/tasks/TASK_dinamr_campaign_2026-09-10.md`,
  `dev/tasks/TASK_dinamr_blockC_grfprobe_2026-09-11.md`
- **Machine:** Mac Studio. Unattended. **Commit only; not pushed.**

## Framing, for the record

- **This does not certify GRF.** Every coverage number produced by this campaign is coverage of
  **β(Ĥ) conditional on the proposed family, over selected replicates**, and every table says so.
  Gate 0a's finding below does not change that: no caption moves on account of it.
- **GRF is not "FS-analogous".** A GRF-to-FS or GRF-to-DINA comparison differs in **identifier**,
  in **family construction**, in **detection set**, and — **at the DR pre-filter only, not at
  admission** — in the **scale of the selection criterion**.
- **No acceptance criterion and no recommendation appear anywhere in this campaign.**

---

# GATE 0 — three source determinations, no compute

## 0a. Which candidate set the MR resampling re-evaluates

**Determination: the FULL ENUMERATED POOL, not the forest-admitted subset.** The admission floor
is applied *per draw inside MR*, to the perturbed effects, and not to the family handed to MR.

### The trace, object by object

| # | file:line | what happens |
|---|---|---|
| 1 | `R/grf_subgroup_labels.R:255–277` | `.grf_dr_candidates()` enumerates depth-1 candidates from `stats::quantile(xj, grid_probs)` for each column of X, keeping those with `n_min <= nS <= n-1`. The **enumeration** reads only X and `n_min`. The DR scores enter one column, `effect = mean(ctrl[S]) - mean(trt[S])`. |
| 2 | `R/grf_subgroup_labels.R:281–309` | `.grf_dr_candidates_d2()` does the same for covariate-pair conjunctions. |
| 3 | `R/grf_main.R:269–276` | `.mr_candidates` is the **union of (1) and (2)**, built once. **No `dmin` filter is applied at this point**, and no frontier is taken. |
| 4 | `R/grf_main.R:284–294` | The frontier path builds its **own separate** enumeration `cand`, and `.grf_frontier_select(cand, dmin = config$dmin.grf, …)` selects from it. This object is **not** what reaches MR. |
| 5 | `R/grf_main.R:331` (frontier) and `R/grf_main.R:419` (tree) | `candidates = .mr_candidates` — the unfiltered pool from (3) — is what is attached to the result. |
| 6 | `R/forestsearch_helpers.R:1621` | `.grf_reselect_on_effect()` **adds** a `sel_effect` column to `grf_res$candidates`. It adds; it does not subset. |
| 7 | `R/forestsearch_helpers.R:1625–1627` | `cand_hr <- grf_res$candidates` is a **local copy**, which *is* filtered (to finite effects) and then admitted on `dmin_eff`. It is **never written back** to `grf_res$candidates`. |
| 8 | `R/forestsearch_helpers.R:1654` | `grf_res$admitted_n <- sum(cand_hr$effect >= dmin_eff)` — the admitted count is **recorded**, but the admitted *set* is used only to pick the winner. |
| 9 | `R/forestsearch_helpers.R:1876` | `.forestsearch_grf_select()` returns `candidates = grf_res$candidates` — still the full pool of (3), now carrying `sel_effect`. |
| 10 | `R/forestsearch_main.R:2520` | `.mr_fam <- .fs_mr_family_from_table(.mr_df, gsel$candidates, op_right = ">", n_min = n.min)`. The comment at `R/forestsearch_main.R:2515–2519` states the intent: *"MR's family is GRF's qualifying candidates — the whole of it."* |
| 11 | `R/fs_mr_inference_methods.R:42–69` | `.fs_mr_family_from_table()` drops a candidate **only** if its re-derived membership has fewer than `n_min` rows — the same `n_min` the enumeration at (1) used. No effect filter. |
| 12 | `R/fs_mr_inference.R:590`, `:1046` | `asm <- .fs_mr_assemble(df, candidates, spec)`; the reported `n_family = length(asm$names)` is therefore **the size of the enumerated pool**, less only the candidates whose per-candidate fit failed. |
| 13 | `R/fs_mr_inference.R:661–666` | Inside the draw loop: `pass <- .admit(bs)` applies `admission$effect_floor` to the **perturbed** effects `beta_star[, b]`, then `.fs_mr_select()` re-selects among the passers. |

**Mechanism, stated plainly.** MR re-evaluates the *whole enumerated pool* on every draw, and the
admission floor is a **per-draw filter on perturbed effects** (13), not a fixed subset chosen once
on the observed data. The deciding lines are **(5)**, which puts the unfiltered pool on the result;
**(7)**, where the filtered copy is local and discarded; and **(13)**, where admission happens.

**No conclusion is drawn here about whether GRF satisfies a fixed-family condition, and no caption
in this campaign changes on account of this finding.** Every table stays labelled
conditional-on-proposed-family. The determination is recorded for Larry.

The answer is readable unambiguously from source, so this is not a STOP.

## 0b. Is `admitted_n` a recorded per-replicate column?

**Determination: ABSENT.** `admitted_n` is computed in the package
(`R/forestsearch_helpers.R:1651` for the empty-admitted-set case and `:1654` otherwise) and it
reaches `forestsearch()`'s return as `out$grf_res$admitted_n`
(`R/forestsearch_main.R:2464`, `grf_res = gsel$grf_res`). But the **simulation template's
recorder never read it**: before this campaign, `grep -n "admitted_n"` over
`sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` returned nothing, and `.na_record()` (line 825
onward) declared `n_family`, `n_cons_qual` and `band_n` but no `admitted_n`.

**Part T2 therefore ran.** It is recorded below.

## 0c. The verification standard

**GRF fits are reproducible within a session (3 of 3) but not across contexts (0 of 179 rows
matched between the probe diagnostic and the pipeline).** What that does and does not affect:

- **It does NOT affect Amendment 3's same-draws assertion,** and this is confirmed from source
  rather than assumed. The assertion compares `n_true` and `truth`:
  - `n_true` is written at the recorder's `rec$n_true <- sum(df[[harm_col]] == 1L, na.rm = TRUE)`
    (template line 995), where `df` comes from
    `simulate_from_dgm(dgm, n = n_sample, …, seed = seed_base + sim_id)` (template lines 969–972).
    That is **before the `forestsearch()` call** (template line 1036) and does not read any GRF
    object.
  - `truth` is built in the `build-dgm` chunk (template lines 716–722) from `calibrate_k_inter()`,
    `setup_gbsg_dgm()` and `compute_dgm_cde()`. It reads `dgm_model`, `target_hr_harm`,
    `harm_z1_quantile`, `n_super` and `seed_base`, and **nothing on the identifier path**. It does
    not branch on `subgroup_method`.
  - Both are therefore functions of the DGM and the seed alone. A GRF fit that differs across
    contexts cannot move either.
- **It DOES affect any GRF diagnostic re-run outside the pipeline.** Stated plainly: **a GRF
  diagnostic re-run outside the pipeline measures statistically equivalent but not identical
  fits.** Any such re-run may be read for distributional shape and not for row-level agreement
  with a pipeline bundle. Every GRF quantity in this campaign is read from the campaign's own
  bundles for that reason.
- It also constrains what Gate 3 can check: see below.

---

# PART T2 — recording `admitted_n`

Run because Gate 0b found the column absent. **Add-only, template-level, no `R/` change.**

- **Declaration.** `.na_record()` gains `admitted_n = NA_integer_` beside `n_family`,
  `n_cons_qual` and `band_n`.
- **Population.** `record_replicate()` writes
  `rec$admitted_n <- as.integer(fs.est$grf_res$admitted_n)` on the **GRF path only**, and
  **before the no-detection return**. The placement is deliberate: an empty admitted set is
  written as `0L` (`R/forestsearch_helpers.R:1651`) and is exactly *why* that replicate did not
  detect, so recording it only on detected rows would drop the informative zeros. Every non-GRF
  path leaves the `NA_integer_` in place.

## Gate T2 — PASS

Five replicates at the standing identity cell (effMaxSG ε 0.20, HR 1.50, n 500, `FS_S7_Z1Q=0.60`,
seeds 8316951 + sim_id, sim_id 1–5) with **`FS_S7_METHOD` unset**, pre-change against post-change
on this machine. Driver `scripts_dinamr/t2gate.sh`, checker `scripts_dinamr/t2gate.R`.

| check | result |
|---|---|
| 5 rows each | PASS (pre 5, post 5) |
| the only new column is `admitted_n` | PASS |
| no column removed | PASS |
| pre-existing columns keep their names **and** their order | PASS |
| every non-timing column `identical()` | **PASS — 163 columns, none mismatched** |
| `truth` `identical()` | PASS |
| `admitted_n` NA on every row of the non-GRF run | PASS |

**7 passes, 0 failures.** Timing columns excluded by name: `fit_mr_secs`, `fb_secs`,
`fld_H_secs`, `fld_Hc_secs`. Both bundles are committed
(`results/fs_effMaxSG_…_t2pre_res_1_5.rds` and `…_t2post_res_1_5.rds`), walls 32 s and 31 s.

## The GRF smoke — `admitted_n` is populated, not merely declared

Gate T2 proves the column is `NA` **off** the GRF path; it cannot prove it is populated **on** it.
A 5-replicate GRF smoke (`scripts_dinamr/t2grf_smoke.sh`, 26 s, committed bundle) settles that:

| sim_id | status | `n_family` | `admitted_n` | `n_sel` |
|---|---|---|---|---|
| 1 | DETECTED | 779 | 115 | 100 |
| 2 | DETECTED | 776 | 152 | 98 |
| 3 | DETECTED | 784 | 49 | 85 |
| 4 | DETECTED | 729 | 146 | 139 |
| 5 | DETECTED | 769 | 57 | 82 |

`admitted_n` is `integer`, finite on 5 of 5, and ranges 49–152 against an `n_family` of 729–784.
**That contrast is the campaign's core diagnostic in miniature**: the enumerated pool barely moves,
the admitted set moves by a factor of three.

---

# GATE 3 — GRF alignment

`scripts_dinamr/gate3.R`, run per batch, stop-on-failure. The three values are **template
literals**, not `FS_S7_*` knobs, so the gate resolves them three ways.

| value | required | resolved | how |
|---|---|---|---|
| `grf_select_statistic` | `"effect"` | **`"effect"`** | source, `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd:504` |
| `grf_selection` | `"frontier"` | **`"frontier"`** | source, template line 503; **and** the batch's own audit line, `Run config: method=grf/frontier` |
| `dmin.grf` | `0.0` | **`0`** | source, template line 506 |

Two corroborations beyond the source read:

- The batch's **own audit line** in the rendered document reads `method=grf/frontier`. (It goes to
  the document, not to the render log — `render.sh` captures only quarto's progress output — so
  the gate reads the `.html`.)
- **`admitted_n` finite in the bundle is direct evidence from the code path itself.** It is written
  only by `.grf_reselect_on_effect()`, which `.forestsearch_grf_select()` reaches only when
  `grf_select_statistic == "effect"` **and** `grf_selection == "frontier"`
  (`R/forestsearch_helpers.R:1773–1776`). This is stronger than the audit line.

**`dmin.grf` has no read-back and the gate does not claim one.** It is the DR-score pre-filter
(`R/grf_main.R:291`), consumed before any quantity that reaches the bundle. Source is the only
resolution available for it, which is stated rather than implied. On the smoke: **7 passes, 0
failures.**

---

# GATE 1 — compute go/no-go

- **Ceiling 9 h wall. Hard timeout 10 h.** Both recorded, both from the kickoff; Larry's available
  window is 9–10 h.
- Projector: `scripts_dinamr/projectG.R`, modelled on `projectC.R`.

## Method

- **From the distribution, not a mean.** Each cell's 2,000-replicate total is bootstrap-resampled
  (B = 4,000) from the measured per-replicate `fit_mr_secs` of the relevant probe, so each cell
  carries a 90% band.
- **Overhead charged per batch render**, three per cell (batch 1, batch 1001, combine), at the
  per-block values `projectC.R` measured on the eleven completed `dinamr` cells: **122.8 s** at
  12.4% and **196.6 s** at 31%. Not 3 × uniform.
- **Calibrated on the `dinamr` realized record and reported, not applied.** Block A ran 1.27 over
  projection *before* the overhead correction and **0.949 after**; Block B **0.92**. The corrected
  projection is the one used here, so the residual calibration is within ±8% either way.
- **Interpolation is across n within a prevalence, never across family size.** The probes found
  cost **flat in family size** — Pearson ρ(s, K) between −0.155 and +0.302, Spearman between
  −0.171 and +0.123 — the opposite of DINA's profile. n = 1000 mixes the n 500 and n 1500 draw
  pools 50/50, interpolating the distribution rather than only its mean. HR 1.75 is unprobed and
  is costed from its own prevalence/n HR 1.50 pool, which is stated in the table's `basis` column
  rather than hidden in a constant.

## The probes, as measured

| prevalence | HR | n | detection | K med | K min–max | s q10 | **s med** | s q90 | s max | ρ(s,K) |
|---|---|---|---|---|---|---|---|---|---|---|
| 12.4% | 1.50 | 500 | 1.0000 | 776 | 729–853 | 12.69 | **14.18** | 15.64 | 16.70 | +0.181 |
| 12.4% | 1.50 | 1500 | 1.0000 | 828 | 820–838 | 15.62 | **16.96** | 18.96 | 21.48 | −0.131 |
| 31% | 1.50 | 500 | 1.0000 | 776 | 729–853 | 14.85 | **16.00** | 17.23 | 18.33 | +0.302 |
| 31% | 1.50 | 1500 | 1.0000 | 828 | 820–838 | 18.36 | **19.79** | 22.10 | 22.87 | −0.155 |
| 12.4% | 1.00 | 500 | 0.9722 | 776 | 729–853 | 12.13 | **13.50** | 15.27 | 16.35 | +0.183 |

**Note the K column.** `n_family` is *identical* at 776 (n 500) and 828 (n 1500) across the two
prevalences — the direct measurement of the claim that GRF's enumerated pool does not depend on
the outcome.

## The twelve cells

| order | prevalence | HR | n | basis | h (q05) | **h (median)** | h (q95) | cumulative |
|---|---|---|---|---|---|---|---|---|
| 1 | 12.4% | 1.50 | 500 | measured | 0.7541 | **0.7561** | 0.7580 | 0.756 |
| 2 | 12.4% | 1.50 | 1000 | pooled draws | 0.8231 | **0.8264** | 0.8297 | 1.583 |
| 3 | 12.4% | 1.50 | 1500 | measured | 0.8944 | **0.8968** | 0.8991 | 2.479 |
| 4 | 31% | 1.50 | 500 | measured | 0.9050 | **0.9066** | 0.9082 | 3.386 |
| 5 | 31% | 1.50 | 1000 | pooled draws | 0.9940 | **0.9978** | 1.0017 | 4.384 |
| 6 | 31% | 1.50 | 1500 | measured | 1.0863 | **1.0888** | 1.0913 | 5.473 |
| 7 | 12.4% | 1.75 | 500 | HR 1.50 pool | 0.7542 | **0.7560** | 0.7580 | 6.229 |
| 8 | 12.4% | 1.75 | 1000 | pooled; HR 1.50 pool | 0.8231 | **0.8264** | 0.8299 | 7.055 |
| 9 | 12.4% | 1.75 | 1500 | HR 1.50 pool | 0.8945 | **0.8968** | 0.8992 | 7.952 |
| 10 | 31% | 1.75 | 500 | HR 1.50 pool | 0.9050 | **0.9067** | 0.9082 | **8.858** |
| 11 | 31% | 1.75 | 1000 | pooled; HR 1.50 pool | 0.9939 | **0.9978** | 1.0017 | 9.856 |
| 12 | 31% | 1.75 | 1500 | HR 1.50 pool | 1.0864 | **1.0888** | 1.0912 | 10.945 |

**Twelve cells project to 10.945 h (90% band 10.914–10.976 h) against a 9 h ceiling.** At the
`dinamr` calibration of 0.949 the same twelve would read **10.387 h** — still over. Chat's rough
figure of ≈ 9.3 h compute plus ≈ 1 h overhead ≈ 10.3 h is reproduced to within about 6%, and the
expectation that cells would defer is confirmed.

## The decision

**GATE 1: GO for ten cells; defer two, from the tail of the run order.**

- **Run (cumulative 8.858 h against the 9 h ceiling):** orders 1–10 — the complete HR 1.50
  n-trajectory at both prevalences, the complete HR 1.75 trajectory at 12.4%, and 31% HR 1.75 at
  n 500.
- **Deferred:** order 11 (31%, HR 1.75, n 1000, 0.998 h) and order 12 (31%, HR 1.75, n 1500,
  1.089 h). They are the tail of the stated run order.
- **Replicate count is not reduced.** 2,000 per cell is what makes the grid comparable to the FS
  and DINA grids and was not traded for cell count.
- The HR 1.00 null cells are not in this task.

Cell lists committed as `scripts_dinamr/grfmr.cells` (ten) and
`scripts_dinamr/grfmr_deferred.cells` (two).

---

# PART A — the ten cells

**All ten cells completed. Part A wall 7.939 h (28,579 s) against the 9 h ceiling — headroom
1.061 h, and 2.061 h unused under the 10 h hard timeout.** Started 03:05:42, finished 11:02:01;
the realized span equals the sum of the per-cell walls, so there was no inter-cell slack to
account for. Driver `scripts_dinamr/grfmr.sh` with `scripts_dinamr/grfmr.cells`.

## Walls, realized against the Gate 1 projection

| order | cell | realized | projected | ratio |
|---|---|---|---|---|
| 1 | 12.4% HR 1.50 n 500 | 0.7044 h (2536 s) | 0.7561 h | 0.932 |
| 2 | 12.4% HR 1.50 n 1000 | 0.7511 h (2704 s) | 0.8264 h | 0.909 |
| 3 | 12.4% HR 1.50 n 1500 | 0.8086 h (2911 s) | 0.8968 h | 0.902 |
| 4 | 31% HR 1.50 n 500 | 0.7883 h (2838 s) | 0.9066 h | 0.870 |
| 5 | 31% HR 1.50 n 1000 | 0.8778 h (3160 s) | 0.9978 h | 0.880 |
| 6 | 31% HR 1.50 n 1500 | 0.9403 h (3385 s) | 1.0888 h | 0.864 |
| 7 | 12.4% HR 1.75 n 500 | 0.6986 h (2515 s) | 0.7560 h | 0.924 |
| 8 | 12.4% HR 1.75 n 1000 | 0.7625 h (2745 s) | 0.8264 h | 0.923 |
| 9 | 12.4% HR 1.75 n 1500 | 0.8119 h (2923 s) | 0.8968 h | 0.905 |
| 10 | 31% HR 1.75 n 500 | 0.7950 h (2862 s) | 0.9067 h | 0.877 |
| | **TOTAL** | **7.939 h** | **8.858 h** | **0.8962** |

**The projection ran 10.4% high, and it ran high uniformly** — every cell between 0.864 and 0.932,
no outlier. That sits just below the `dinamr` Block B anchor of 0.92 and below the corrected Block
A anchor of 0.949. The projection method — bootstrap over the per-replicate cost distribution plus
measured per-render overhead — therefore carried over from DINA to GRF without re-tuning, which is
the one thing a projection method has to do.

**The HR 1.75 cells validate their own cost basis.** They were costed from the HR 1.50 pool at
matched prevalence and n because HR 1.75 was unprobed. Their realized ratios — 0.924, 0.923,
0.905, 0.877 — are indistinguishable from the HR 1.50 ratios at the same coordinates (0.932,
0.909, 0.902, 0.870). The assumption was sound and is now measured.

**The deferral was correct.** At the realized ratio the two deferred cells would have cost 1.870 h,
putting all twelve at **9.809 h — 0.809 h over the 9 h ceiling**. Deferring them was not merely
defensible on the projection; it was necessary on the realized cost.

## Deferred, not dropped

| cell | Gate 1 projection | at the realized 0.8962 |
|---|---|---|
| 31% HR 1.75 n 1000 | 0.9978 h | 0.894 h |
| 31% HR 1.75 n 1500 | 1.0888 h | 0.976 h |

Committed as `scripts_dinamr/grfmr_deferred.cells`. **No cell was dropped and no cell's replicate
count was reduced**; every completed cell is 2,000 replicates in two seed-disjoint batches of
1,000, combined.

## Gate 3 — alignment, per batch

**Run on all twenty batches (ten cells x two batches). 7 passes, 0 failures every time.**
`grf_select_statistic` resolved to `"effect"`, `grf_selection` to `"frontier"`, `dmin.grf` to `0`,
on every batch, by all three resolutions described above.

## Gate 2 — per cell

**Every check PASS on all ten cells.** No failure, no warning, nothing deferred to judgement.

- **Completeness**: 2,000 rows, `sim_id` 1..2000 with no duplicates, no CONFIG-ERROR, both batches
  present, `n_workers` 12 and `forestsearch_version` 0.3.5 recorded in every batch meta, host
  `Mac-Studio-3.local`, R 4.5.2, `seed_base` 8316951.
- **Knobs in meta**: `subgroup_method` grf, `sg_focus` effMaxSG, `effect_neighborhood` 0.20,
  `field_complement` / `field_decompose` / `field_recovery` TRUE, `field_scale_complement`
  selected, `ij_residual` two_term, `campaign_tag` grfmr — all ten cells.
- **Finiteness**: all 38 gated quantities finite on every detected replicate, all ten cells. The
  nine recovery columns, the p-hat block and rho-c present and populated throughout.
- **Interval invariants**: all eleven hold on all ten cells, including the added
  `admitted_n >= 1 on every detected replicate`.
- **The corrected identity** `log(fld_Hc_est2_s) + fld_Hc_lam_mean_s == log(fld_Hc_est2) +
  fld_Hc_lam_mean`: **max |diff| 2.22e-16 to 3.33e-16** across the ten cells.
- **The Bonferroni identity**, gated on the gamma-at-floor rows: **max |diff| exactly 0** on every
  cell. Share at the floor 0.876-0.911 (joint) and 0.812-0.866 (joint-s).
- **gamma** within [0.02500, 0.02800] on every cell — inside [0.025, 0.05] throughout.
- **Realized prevalence**: 0.12363-0.12401 against a super-population 0.12418 at 12.4%;
  0.30591-0.30625 against 0.30655 at 31%.
- **Structurally-NA columns, reported as such and never as failures**: `n_cons_qual` and `band_n`
  are present and all-NA on every cell — structural on GRF, which has no consistency screen —
  and `p_star` is not a recorder column at all (an admission-set term, NULL on GRF).

### Amendment 3 — the same-draws assertion, all ten cells

**`n_true` `identical()` on all 2,000 rows of every cell. `truth` `all.equal()` at 1e-8 YES on
every cell.** No cell shows a DGM-path mismatch, so there is no finding to record under that head.

| cells | `truth` `identical()` | max abs diff | max rel diff |
|---|---|---|---|
| 12.4%, HR 1.75 (all three) | **TRUE** | **0** | **0** |
| 12.4%, HR 1.50 (all three) | FALSE | 6.661e-16 | 3.986e-16 |
| 31%, HR 1.50 (all three) | FALSE | 5.551e-15 | 3.249e-15 |
| 31%, HR 1.75 n 500 | FALSE | 8.882e-15 | 4.253e-15 |

`identical()` is reported and not asserted, as the predecessor's gate specifies: the FS comparator
bundles are Linux-produced and cross-machine BLAS moves `truth` at the 1e-16 level. Every
discrepancy above is at least seven orders inside the 1e-8 tolerance. The three 12.4% HR 1.75
cells are bit-identical.

## The per-cell record


### Detection, admitted_n and the enumerated pool

| cell | detection | n_eval | **admitted_n** min / q25 / **med** / q75 / p90 / max | adm CV | adm = 0 | adm NA | n_family min / **med** / max | K CV | med share |
|---|---|---|---|---|---|---|---|---|---|
| 12.4% HR 1.50 n 500 | 0.9985 [0.9956, 0.9995] | 1997 | 2 / 65 / **115** / 179 / 260.4 / 600 | 0.706 | 0 | 3 | 712 / **776** / 870 | 0.0181 | 0.1469 |
| 12.4% HR 1.50 n 1000 | 0.9975 [0.9942, 0.9989] | 1995 | 3 / 57 / **93** / 145 / 198 / 471 | 0.634 | 0 | 5 | 779 / **830** / 914 | 0.0111 | 0.1122 |
| 12.4% HR 1.50 n 1500 | 0.9975 [0.9942, 0.9989] | 1995 | 2 / 47 / **76** / 114 / 156.6 / 369 | 0.618 | 0 | 5 | 780 / **830** / 913 | 0.0081 | 0.0910 |
| 12.4% HR 1.75 n 500 | 0.9995 [0.9972, 0.9999] | 1999 | 3 / 81.5 / **136** / 208 / 293 / 624 | 0.646 | 0 | 1 | 712 / **776** / 870 | 0.0181 | 0.1746 |
| 12.4% HR 1.75 n 1000 | 0.9990 [0.9964, 0.9997] | 1998 | 6 / 78 / **121** / 176 / 234 / 510 | 0.558 | 0 | 2 | 779 / **830** / 914 | 0.0111 | 0.1464 |
| 12.4% HR 1.75 n 1500 | 1.0000 [0.9981, 1.0000] | 2000 | 5 / 69 / **105** / 147 / 191 / 439 | 0.525 | 0 | 0 | 780 / **830** / 913 | 0.0081 | 0.1260 |
| 31% HR 1.50 n 500 | 1.0000 [0.9981, 1.0000] | 2000 | 16 / 269 / **383** / 497 / 578.1 / 752 | 0.382 | 0 | 0 | 712 / **776** / 870 | 0.0181 | 0.4932 |
| 31% HR 1.50 n 1000 | 1.0000 [0.9981, 1.0000] | 2000 | 52 / 310 / **406** / 504 / 586 / 758 | 0.323 | 0 | 0 | 779 / **830** / 914 | 0.0111 | 0.4898 |
| 31% HR 1.50 n 1500 | 1.0000 [0.9981, 1.0000] | 2000 | 96 / 309 / **400.5** / 491 / 555 / 734 | 0.297 | 0 | 0 | 780 / **830** / 913 | 0.0081 | 0.4818 |
| 31% HR 1.75 n 500 | 1.0000 [0.9981, 1.0000] | 2000 | 53 / 341 / **449** / 545.25 / 609.1 / 755 | 0.309 | 0 | 0 | 712 / **776** / 870 | 0.0181 | 0.5782 |

### Coverage of every product, absolute levels [Wilson]

Conditional-on-proposed-family estimand throughout, on detected replicates.

| cell | naive H | field H (1-sided) | field-s Hc (1-sided) | IJ 2-sided H | IJ 2-sided Hc |
|---|---|---|---|---|---|
| 12.4% HR 1.50 n 500 | 0.4982 [0.4763, 0.5202] | 0.9434 [0.9324, 0.9527] | 0.9319 [0.9200, 0.9421] | 0.9890 [0.9834, 0.9927] | 0.9995 [0.9972, 0.9999] |
| 12.4% HR 1.50 n 1000 | 0.6561 [0.6350, 0.6767] | 0.9404 [0.9291, 0.9499] | 0.9519 [0.9416, 0.9604] | 0.9870 [0.9810, 0.9911] | 1.0000 [0.9981, 1.0000] |
| 12.4% HR 1.50 n 1500 | 0.7830 [0.7643, 0.8005] | 0.9579 [0.9482, 0.9659] | 0.9529 [0.9427, 0.9613] | 0.9799 [0.9728, 0.9852] | 0.9995 [0.9972, 0.9999] |
| 12.4% HR 1.75 n 500 | 0.5648 [0.5429, 0.5864] | 0.9415 [0.9303, 0.9509] | 0.9335 [0.9217, 0.9436] | 0.9900 [0.9846, 0.9935] | 0.9995 [0.9972, 0.9999] |
| 12.4% HR 1.75 n 1000 | 0.7297 [0.7098, 0.7487] | 0.9364 [0.9249, 0.9463] | 0.9540 [0.9439, 0.9623] | 0.9775 [0.9700, 0.9831] | 1.0000 [0.9981, 1.0000] |
| 12.4% HR 1.75 n 1500 | 0.8385 [0.8217, 0.8540] | 0.9535 [0.9434, 0.9619] | 0.9560 [0.9461, 0.9641] | 0.9755 [0.9678, 0.9814] | 1.0000 [0.9981, 1.0000] |
| 31% HR 1.50 n 500 | 0.6065 [0.5849, 0.6277] | 0.9295 [0.9174, 0.9399] | 0.9185 [0.9057, 0.9297] | 0.9840 [0.9775, 0.9886] | 1.0000 [0.9981, 1.0000] |
| 31% HR 1.50 n 1000 | 0.7545 [0.7352, 0.7729] | 0.9525 [0.9423, 0.9610] | 0.9460 [0.9352, 0.9551] | 0.9670 [0.9582, 0.9740] | 0.9990 [0.9964, 0.9997] |
| 31% HR 1.50 n 1500 | 0.8570 [0.8410, 0.8717] | 0.9690 [0.9605, 0.9757] | 0.9480 [0.9374, 0.9569] | 0.9730 [0.9649, 0.9792] | 0.9995 [0.9972, 0.9999] |
| 31% HR 1.75 n 500 | 0.6515 [0.6303, 0.6721] | 0.9300 [0.9180, 0.9404] | 0.9210 [0.9084, 0.9320] | 0.9730 [0.9649, 0.9792] | 1.0000 [0.9981, 1.0000] |

### Classification and bound location

| cell | sens | spec | PPV | NPV | mean \|Hhat\| | med bound | med theta(Hhat) | share >= 1.00 | share >= 1.25 |
|---|---|---|---|---|---|---|---|---|---|
| 12.4% HR 1.50 n 500 | 0.5158 | 0.8364 | 0.3187 | 0.9255 | 103.8 | 0.4665 | 0.8373 | 0.0230 [0.0173, 0.0306] | 0.0045 [0.0024, 0.0085] |
| 12.4% HR 1.50 n 1000 | 0.6438 | 0.8855 | 0.4806 | 0.9466 | 180.2 | 0.6053 | 0.9777 | 0.0612 [0.0515, 0.0725] | 0.0160 [0.0114, 0.0226] |
| 12.4% HR 1.50 n 1500 | 0.7448 | 0.8959 | 0.5515 | 0.9617 | 275.5 | 0.6977 | 1.0451 | 0.1083 [0.0954, 0.1227] | 0.0296 [0.0230, 0.0380] |
| 12.4% HR 1.75 n 500 | 0.5746 | 0.8499 | 0.3661 | 0.9349 | 101.6 | 0.5073 | 0.9260 | 0.0420 [0.0341, 0.0517] | 0.0095 [0.0061, 0.0148] |
| 12.4% HR 1.75 n 1000 | 0.6975 | 0.9078 | 0.5634 | 0.9552 | 167.2 | 0.7054 | 1.1497 | 0.1582 [0.1428, 0.1748] | 0.0561 [0.0468, 0.0670] |
| 12.4% HR 1.75 n 1500 | 0.7770 | 0.9207 | 0.6394 | 0.9671 | 248.7 | 0.8301 | 1.2468 | 0.2775 [0.2583, 0.2975] | 0.1080 [0.0951, 0.1224] |
| 31% HR 1.50 n 500 | 0.5219 | 0.8986 | 0.6785 | 0.8162 | 115.0 | 0.6787 | 1.2627 | 0.1365 [0.1221, 0.1522] | 0.0480 [0.0395, 0.0583] |
| 31% HR 1.50 n 1000 | 0.6831 | 0.9305 | 0.8053 | 0.8785 | 257.5 | 0.8314 | 1.3308 | 0.2160 [0.1985, 0.2346] | 0.0605 [0.0509, 0.0718] |
| 31% HR 1.50 n 1500 | 0.8176 | 0.9342 | 0.8520 | 0.9282 | 444.4 | 0.9596 | 1.4254 | 0.4050 [0.3837, 0.4267] | 0.0825 [0.0712, 0.0954] |
| 31% HR 1.75 n 500 | 0.5823 | 0.9217 | 0.7493 | 0.8404 | 116.2 | 0.7997 | 1.4710 | 0.2645 [0.2456, 0.2843] | 0.1100 [0.0970, 0.1245] |

### Beside the FS comparator, with its criterion named

| cell | matched? | FS comparator | GRF det / FS det | GRF field / FS field | GRF IJ2 / FS IJ2 | GRF pool / FS family (med) |
|---|---|---|---|---|---|---|
| 12.4% HR 1.50 n 500 | no | p12ext / maxeffCons / eps 0.1 | 0.9985 / 0.9110 | 0.9434 / 0.9654 | 0.9890 / 0.9769 | 776 / 1223 |
| 12.4% HR 1.50 n 1000 | no | p12ext / maxeffCons / eps 0.1 | 0.9975 / 0.9740 | 0.9404 / 0.9410 | 0.9870 / 0.9322 | 830 / 1299 |
| 12.4% HR 1.50 n 1500 | no | p12ext / maxeffCons / eps 0.1 | 0.9975 / 0.9880 | 0.9579 / 0.9565 | 0.9799 / 0.9008 | 830 / 1297 |
| 12.4% HR 1.75 n 500 | no | tier2 / maxeffCons / eps 0.1 | 0.9995 / 0.9500 | 0.9415 / 0.9663 | 0.9900 / 0.9611 | 776 / 1223 |
| 12.4% HR 1.75 n 1000 | no | tier2 / maxeffCons / eps 0.1 | 0.9990 / 0.9950 | 0.9364 / 0.9437 | 0.9775 / 0.9166 | 830 / 1299 |
| 12.4% HR 1.75 n 1500 | no | tier2 / maxeffCons / eps 0.1 | 1.0000 / 0.9990 | 0.9535 / 0.9600 | 0.9755 / 0.9129 | 830 / 1297 |
| 31% HR 1.50 n 500 | **yes** | e1stud / effMaxSG / eps 0.2 | 1.0000 / 0.9995 | 0.9295 / 0.9745 | 0.9840 / 0.9810 | 776 / 1223 |
| 31% HR 1.50 n 1000 | **yes** | cert20 / effMaxSG / eps 0.2 | 1.0000 / 1.0000 | 0.9525 / 0.9585 | 0.9670 / 0.9710 | 830 / 1299 |
| 31% HR 1.50 n 1500 | **yes** | cert20 / effMaxSG / eps 0.2 | 1.0000 / 1.0000 | 0.9690 / 0.9615 | 0.9730 / 0.9755 | 830 / 1297 |
| 31% HR 1.75 n 500 | **yes** | e1stud / effMaxSG / eps 0.2 | 1.0000 / 0.9995 | 0.9300 / 0.9700 | 0.9730 / 0.9720 | 776 / 1223 |

Confound, with every comparison: FS and GRF differ in identifier, in family construction and
in detection set; at 12.4% they differ in the selection criterion as well (maxeffCons eps 0.10
against effMaxSG eps 0.20), so a 12.4% gap cannot be read as engine behaviour even in part.


## What the Part A numbers say, descriptively

**Detection.** **1.0000 at all four 31% cells and at 12.4% HR 1.75 n 1500; 0.9975-0.9995 at the
other five.** GRF's probe result of 1.0000 holds at full cell size at 31% and is within 0.0025 of
it at 12.4%. **The operative property is that it is flat in n** — 0.9985 / 0.9975 / 0.9975 at
12.4% HR 1.50 and 0.9995 / 0.9990 / 1.0000 at 12.4% HR 1.75, against DINA's 0.7135 / 0.5245 /
0.3435 on the same 12.4% draws. **A detection rate that does not move with n cannot be the reason
a coverage or bias figure moves with n, so the detection-conditioning confound that qualifies
every DINA n-trend does not qualify these.** At 31% it is removed outright. This is a statement
about what a comparison may be read as; it is not a claim about GRF's performance.

**The enumerated pool is outcome-independent, measured.** At each n the `n_family` distribution is
**identical quantile for quantile across both prevalences and both hazard ratios**: 712 / 771 /
776 / 780 / 784 / 870 at n 500 on all four such cells, and 780 / 826 / 830 / 833 / 836 / 913 at n
1500. Its CV is 0.0181 / 0.0111 / 0.0081 at n 500 / 1000 / 1500 and depends on n alone. Meanwhile
`admitted_n` at n 500 runs median 115 (12.4% HR 1.50), 136 (12.4% HR 1.75), 383 (31% HR 1.50) and
449 (31% HR 1.75) — **a factor of 3.9 across cells whose enumerated pool is the same object.**
rho(`admitted_n`, `n_family`) is +0.018 to +0.042 on every cell. **This is the whole justification
for the Part T2 recorder change and the Stage 3 stratifier substitution, and it is measured rather
than argued.**

**`admitted_n` moves differently by prevalence.** It **falls** with n at 12.4% (115 / 93 / 76 at
HR 1.50; 136 / 121 / 105 at HR 1.75) and is **flat** at 31% (383 / 406 / 400 at HR 1.50). Its
median share of the pool falls 0.147 / 0.112 / 0.091 at 12.4% and holds at 0.493 / 0.490 / 0.482
at 31%. Its CV falls with n everywhere (0.71 -> 0.62 at 12.4% HR 1.50; 0.38 -> 0.30 at 31% HR
1.50), so the admitted set does stabilize in relative spread even where its level does not.
**`admitted_n` is never 0 on any of the 20,000 replicates**, so the empty-admitted-set path never
fired; the 21 non-detections across the ten cells all carry `admitted_n` NA, which means the
re-selection returned before reaching the admission step rather than admitting nothing. Part T2
makes that distinction visible, and DINA's columns could not draw it.

**Coverage, absolute levels, conditional on the proposed family throughout.** Field H runs
0.9295-0.9690 over the ten cells and rises with n at both prevalences (0.9434 / 0.9404 / 0.9579 at
12.4% HR 1.50; 0.9295 / 0.9525 / 0.9690 at 31% HR 1.50). Field-s Hc runs 0.9185-0.9560, also
rising with n. IJ two-sided H runs 0.9670-0.9900 and **decays mildly with n at 12.4%** (0.9890 /
0.9870 / 0.9799 at HR 1.50; 0.9900 / 0.9775 / 0.9755 at HR 1.75), which is the same direction FS
shows on those draws. IJ two-sided Hc is 0.9990-1.0000 everywhere. Naive H runs 0.4982-0.8570 and
rises steeply with n, which is the optimism the correction addresses.

**Classification and bound location.** Sensitivity and PPV rise with n and with HR at both
prevalences. The share of field lower bounds at or above 1.00 rises with n and HR at both — 0.0230
/ 0.0612 / 0.1083 at 12.4% HR 1.50 up to 0.4050 at 31% HR 1.50 n 1500 — and the share at or above
1.25 follows at a third to a quarter of the level.

**Beside FS, and the confound.** On the **four criterion-matched 31% rows** (`effMaxSG`, eps 0.20)
GRF's field bound reads 0.9295 / 0.9525 / 0.9690 / 0.9300 against FS's 0.9745 / 0.9585 / 0.9615 /
0.9700, and IJ two-sided 0.9840 / 0.9670 / 0.9730 / 0.9730 against 0.9810 / 0.9710 / 0.9755 /
0.9720. On the **six 12.4% rows the criterion is not matched** (FS is `maxeffCons` at eps 0.10),
and the confound sentence travels with every one of them: FS and GRF differ in identifier, in
family construction and in detection set, and at 12.4% in the selection criterion as well, so a
12.4% gap cannot be read as engine behaviour even in part. GRF's enumerated pool is consistently
smaller than FS's family (776 / 830 against 1223 / 1297-1299) and three to six times steadier
(CV 0.0081-0.0181 against 0.0357-0.0536).

**No acceptance criterion is applied and no recommendation is made.**

---

# PART C — FS extraction

Completed while the cells compute, and never ahead of a gate, a projection or a launch. Written to
**`REPORT_fs_extraction_2026-09-11.md`** beside this report, from committed bundles only — no
re-run, no new simulation, no recorder change. All 18 committed FS cells are covered. See that
document for the tables; its headline is that the field lower bound sits below the realized θ(Ĥ)
in every one of the 18 cells (ratio 0.469–0.694), that the share of bounds at or above 1.00 rises
steeply with n and HR on the harm cells (to 0.798 at 31% HR 1.75 n 1500), and that on every null
cell the same share stays at 0.003–0.011 and does not rise with n.

---

# STAGE 3 — the summary

**Rendered: `summary_grfmr.html`** (7.1 MB, all 79 chunks). `summary_grfmr.qmd`, transplanted
from `summary_dinamr.qmd` by
`scripts_dinamr/transplant_grfmr.py` with every structural edit asserted to match exactly once.
`admitted_n` substitutes for `n_family` as the strata-section stratifier with the reason in every
caption it touches; `n_family` is kept in the descriptive tables with an `admitted_n` table and
plot beside it. Per-cell chunks stay guarded, so the two deferred cells skip rather than render
empty. Every coverage column is labelled as the conditional-on-proposed-family estimand.

The render confirms the guards: **"Cells on disk: 10 of 12 harm cells (the HR 1.00 null cells are
out of scope for this task)"**, the two deferred cells listed under "DEFERRED / NOT ON DISK
(omitted from every table below, not rendered empty)", six guarded chunks emitting a skip note
rather than an empty artifact, and the strata section keyed on `adm T1` / `adm T2` / `adm T3`
rather than on `n_family`.

---

# Artifacts, and one decision left open

**Committed:** the task document; the Part T2 template change with its three Gate T2 / smoke
bundles; the tooling (`projectG.R`, `gate2G.R`, `gate3.R`, `grfmr.sh`, `grfmr.cells`,
`grfmr_deferred.cells`, `t2gate.sh`, `t2gate.R`, `t2grf_smoke.sh`, `grfmr_numbers.R`,
`transplant_grfmr.py`, `fs_extraction.R`); the **31 result bundles** of Part A (31 MB, largest
blob 1.5 MB); `summary_grfmr.qmd`; and both reports. **Commit range `a451230d..6765c372`. Not
pushed.**

**Also committed, on Larry's instruction (2026-09-12):** the **30 batch and combine renders**
(`grfmr_*_batch_*.html`, `grfmr_*_combine_*.html`) and **`summary_grfmr.html`** — 31 files,
138 MB, largest blob 6.7 MB. This completes the `dinamr` convention, under which 55 such renders
and `summary_dinamr.html` are already tracked. No blob approaches the 50 MB limit that
`REPORT_push_size_fix_2026-08-31` documents; the constraint that remains is cumulative repository
size at push, which is Larry's to weigh at push time. **Nothing from this campaign is left
untracked.**
