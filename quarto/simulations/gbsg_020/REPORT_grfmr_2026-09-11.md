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

*(Gate 2 records per cell are appended below as cells complete.)*

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

`summary_grfmr.qmd`, transplanted from `summary_dinamr.qmd` by
`scripts_dinamr/transplant_grfmr.py` with every structural edit asserted to match exactly once.
`admitted_n` substitutes for `n_family` as the strata-section stratifier with the reason in every
caption it touches; `n_family` is kept in the descriptive tables with an `admitted_n` table and
plot beside it. Per-cell chunks stay guarded, so the two deferred cells skip rather than render
empty. Every coverage column is labelled as the conditional-on-proposed-family estimand.
