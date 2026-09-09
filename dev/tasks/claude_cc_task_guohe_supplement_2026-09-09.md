# TASK — Guo & He (2021) comparison for the `fs-post-selection` supplement: Tier 1 (complement + joint on t7) and Tier 2(a) (full Adaptive column on t7)

Date: 2026-09-09. Governing proposal: `claude_proposal_guohe_supplement_2026-09-09_v2.md` (self-contained; decisions resolved below). Protocol: Larry decides; CC executes all repository operations; CC never pushes. Summaries in bullet form. Records live beside results in `quarto/GuoHe/`.

Decisions resolved by Larry (2026-09-09): **D1** scope = t7 only (six cells), one cross-reference sentence to the committed 16-cell record. **D2** Tier 1 = GO, Phase A on the Mac. **D3** Tier 2 = option (a), full 6 × 2000, Phase B on the Linux box after Larry's sync. **D4** B2 includes the grayed MR(IJ) context column and the point-estimate bias block (naive / G&H bias-reduced / est2). Closed lines per the handoff stand; in particular their method is run, never revised.

---

## 1. Session constraint (Larry, 2026-09-09) — in force for this entire task

- **No `git fetch` or `git pull` on the Mac Studio.** A zero-compute session is in flight on the Linux box with unpushed commits; the Mac checkout already contains everything this reading needs (the Mac merge and the certification work are in). If any required input is absent from the tree: **STOP and report — never fetch.** Larry sequences the sync (Linux push → Mac pull → push).
- **Provenance first:** open the record with `git log -1 --oneline` and `git status -sb`, quoted verbatim, so the record states exactly which tree the reading reflects.
- **Seven pre-existing untracked files are out of scope:** three `diag_*` bundles, three `diag_*.html`, and the ACTG175 payload. Never modify, never stage. `git add -A`, `git add .`, and any directory-level add are forbidden; **every add is by explicit named path.** The closing record re-runs `git status -sb` and asserts the identical seven remain untracked and untouched.
- **Commit locally as usual; never push.**

## 2. Stage 0 (Mac) — provenance and source verification. STOP on any failure.

1. Record `git log -1 --oneline` and `git status -sb` verbatim at the top of `quarto/GuoHe/REPORT_guohe_supp_stage0_2026-09-09.md`; enumerate the seven untracked files by name.
2. Input inventory — assert each present in the tree (list each with its byte size in the record):
   - Bundles: `quarto/GuoHe/guohe_repro_t7_beta2_0{0..5}.rds`, `mr_vs_guohe_t7_beta2_0{0..5}.rds`, `mr_field_vs_guohe_t7_beta2_0{0..5}.rds`; truth caches `guohe_sec52_truth_beta2_0{0..5}.rds`.
   - Scripts: `guohe_sec52_sim.R`, `guohe_sec52_truth.R`, `guohe_sec52_run.R`, `guohe_reproduction_run.R`, `guohe_reproduction_sim.R`, `mr_vs_guohe_sim.R`, `mr_field_vs_guohe_run.R`, `mr_field_vs_guohe.qmd`; package sources `R/guohe_algorithm3.R`, `R/guohe_adaptive_r.R`, `R/fs_to_guohe.R`.
   - Records: `NOTE_complement_product_2026-09-08.md`, `REVIEW_certification_2026-09-09.md` (locate; paths as in the Mac tree).
3. Source quotes, with file path and line numbers **as read from this tree** (line numbers from any other tree are not evidence):
   - Q1 — the campaign adapter's disabled-complement lines in `mr_vs_guohe_sim.R` (the `field_complement = FALSE, include_complement = FALSE` signature defaults).
   - Q2 — the current engine signature for `fs_mr_inference` (locate its defining file under `R/`): defaults for `field_complement`, `field_scale_complement`, `return_reselection`, and the interface for requesting one-sided bounds at both 0.95 and 0.975 for the selected subgroup and the complement. If the defaults are not `TRUE` / `"selected"` / `TRUE`, that contradicts the recorded adoption — STOP and report.
   - Q3 — `guohe_adaptive_r()` signature: the `r_grid` default (read at `a9c0d5e` as `c(0.03, 0.10, 0.20, 0.30, 0.40, 0.45)`), the `v` default, and the one-B behavior (one `B` serves the inner CV fits and the final refit).
   - Q4 — the `--adaptive` path in `guohe_reproduction_run.R`: the exact `r_grid`, `v`, and `B` its validated Tables 3–6 Adaptive columns used. Record verbatim; a difference from this task's t7 settings is footnoted in the record, not silently reconciled.
   - Q5 — the certification harness's joint-Bonferroni computation: locate the certification driver(s) behind `REVIEW_certification_2026-09-09.md`, and quote the lines that form the 0.975 one-sided lower (Ĥ) and 0.975 one-sided upper (Ĥᶜ) and the both-correct joint indicator. T1 transplants these lines; they are not re-derived.
   - Q6 — transplant anchors in `mr_field_vs_guohe_run.R`: the `--cells` flag block, `MF_CELLS`, the seed structure (`seed_data = base + m`; `seed_mr = base + m + MV_SEED_MR`; `mv_gh52_base`), the `mf_rep_52` gate-call block, and the row-wise `seed_data` identity assertions.
4. Certification quotes for B6: field-s one-sided upper 0.912–0.960 and joint 0.939–0.963, quoted from `REVIEW_certification_2026-09-09.md` with line numbers.

## 3. T0 — commit this task document first

- Copy this file to `dev/tasks/claude_cc_task_guohe_supplement_2026-09-09.md` and commit it alone (named path), before any other change. If the `~/Downloads` transport failed and this spec was reconstructed from the kickoff, commit the reconstructed text verbatim under the same path as the governing spec.

## 4. T1 — Tier 1: complement + joint on t7 (Phase A, Mac)

**Driver (transplant, never fresh authorship).** Copy `quarto/GuoHe/mr_field_vs_guohe_run.R` → `quarto/GuoHe/mr_field_complement_vs_guohe_run.R`. Named-line changes only:
- Default cells: t7 only (`sprintf("t7_beta2_%02d", 0:5)`); keep the `--cells`/`--force`/skip-if-exists scaffolding.
- Gate call: add, explicitly, `field_complement = TRUE`, `include_complement = TRUE`, `field_scale_complement = "selected"`, `return_reselection = TRUE`; all other arguments unchanged from Q6 (MR B = 5000 centred Poisson; `ci_method = "field"`; `field_R_out/field_R_in = 1000/500`; identical seed derivations, so `seed_data` and `seed_mr` match the stored bundles row-for-row).
- Emit per replicate: the stored-comparable columns (naive; IJ; field est2/lower/margin), the complement upper bounds at 0.95 and 0.975, the harm-side lower bounds at 0.95 and 0.975, the joint both-correct indicator (0.975 lower on Ĥ vs γ_s from `gh52_truth_at()`; 0.975 upper on Ĥᶜ vs truth 0 — transplanted per Q5), complement diagnostics per the E5 join pattern, and p̂(Ĥ).
- Output `mr_field_complement_vs_guohe_<id>.rds` per cell, campaign bundle format (per-replicate rows, seeds, gate tallies, `sessionInfo`).

**Stage 1 identity probe (before any production).** `devtools::install()` first (workers see only the installed package). Run 3 replicates on each of `t7_beta2_00` and `t7_beta2_05` and assert, per replicate, `identical()` of the recomputed naive, IJ, and field columns to the stored `mr_field_vs_guohe_<id>.rds` rows. The engine's complement block precedes the field block; if enabling the complement perturbs **any** stored column, STOP and report — that is an RNG-stream finding for Larry, not a tolerance. If emitting the 0.975 bounds requires any change under `R/`, STOP and route through chat first (no engine edits without authorization).

**Gate 1a (pilot → production, pre-authorized envelope).** Pilot `--pilot` at reps = 20 on `t7_beta2_00`; print the 6 × 2000 projection (core-h and Mac wall-clock) into `REPORT_mr_complement_gate1_2026-09-09.md`. Pre-authorization: proceed to production only if projection ≤ 40 core-h and ≤ 90 min Mac wall; otherwise STOP and report. Production: all six cells; commit the six bundles (named paths).

**Record.** `quarto/GuoHe/REPORT_mr_complement_vs_guohe_2026-09-09.md`, bullet form, quoting per-cell: complement upper coverage vs 0 (Wilson), mean upper-bound location, joint coverage (Wilson), harm-side 0.95 coverage cross-checked against the stored field column (must reproduce it), p̂(Ĥ) summary, marginal SD and error SD side by side, wall/core cost, provenance.

## 5. T2 — Tier 2(a): the Adaptive column on t7 (driver authored and committed now; **executed only in Phase B on the Linux box, after Larry declares the sync complete**)

**Driver (transplant).** New file `quarto/GuoHe/guohe_sec52_adaptive_run.R`: the pilot/flag/skip scaffolding of `guohe_sec52_run.R` plus the `--adaptive` named lines of `guohe_reproduction_run.R`. Per replicate m of cell id:
- Regenerate the data from the stored seed; recompute the naive argmax and naive bound and assert `identical()` to the stored `guohe_repro_t7_<id>.rds` row (the pairing proof).
- Run `guohe_adaptive_r()` with `orient = +1`, `r_grid = c(1/3, 1/12, 1/21, 1/30)` (the published Table-7 grid — pinned; their paper leaves the candidate set to the analyst, and this matches the fixed-r columns), `v = 5` (their tables' value), `B = 200` (`--adaptive-B=200`, the validated reproduction setting; the CV objective is a bootstrap mean), derived seed recorded. Record r̂ and the per-candidate objective values.
- **Primary Adaptive bound = the stored B = 2000 Algorithm-3 bound at r̂**, looked up by (id, m) from `guohe_repro_t7_<id>.rds` — exact resolution parity with the fixed-r columns. **Secondary:** the function's own B = 200 final bound, recorded in the bundle and summarized once in the record. Coverage of both against γ_s from the committed truth caches.
- Output `guohe_adaptive_t7_<id>.rds` per cell, bundle format as above.

**Gate 1b (Phase B, Linux).** Precondition: Larry states the Linux→Mac→push sequencing is complete and gives the Phase-B pointer line. Then `devtools::install()`; pilot reps = 20 on `t7_beta2_00`; print the 6 × 2000 projection. Pre-authorization: proceed only if projection ≤ 2,500 core-h and ≤ 24 h wall at the available cores; otherwise STOP and report. Production: six cells; commit bundles; record `quarto/GuoHe/REPORT_guohe_adaptive_t7_<run-date>.md` (bullet form: Adaptive coverage primary and secondary with Wilson, r̂ distribution per cell, agreement with the best fixed-r column, cost, provenance; carry the two caveats — their Table 6 adaptive caution, and the independent-draws-across-r variance note).

## 6. T3 — assembly qmd (Phase A after T1; Adaptive column joins after Phase B)

- New `quarto/GuoHe/guohe_supp_section.qmd`, transplanted from `mr_field_vs_guohe.qmd` (structure, bundle-reading, Wilson machinery). Never convert formats; render with RStudio's bundled Quarto binary on the executing machine and record the binary path used. HTML untracked per the `GuoHe/*.html` convention.
- **B2** (committed bundles only): per β₂ — coverage of the one-sided 95% lower bound for γ_s (Naive; G&H r ∈ {1/3, 1/12, 1/21, 1/30}; Adaptive when `guohe_adaptive_t7_*.rds` exist, rendered "pending Phase B" otherwise; MR field; MR(IJ) grayed as context), mean margin, and point-estimate bias (naive β̂; G&H bias-reduced; est2). Footnote the ≈ 0.01 replication deficit against print. MCSE ≈ 0.005 stated.
- **B3** (from T1 bundles; caption "capability, not comparison"; visually separate): field-s upper coverage vs 0 and mean location per cell; joint coverage per cell; the two fixed sentences below, verbatim.
- **B5**: distribution of ĉ (note sensitivity ≡ 1 by nesting, S(30) ⊆ S(ĉ); specificity/PPV/NPV as functions of ĉ), p̂(Ĥ), M_eff, marginal SD beside error SD.
- **B6**: the honest-limits paragraph — field low points 0.933–0.937 as the stable-pick under-correction on a design not flattering to us; complement/joint here are new measurements versus the certified 0.912–0.960 and 0.939–0.963 (Q4 quotes); G&H not run on our screened, data-built family (outside its stated scope); one-sided constructions only, hence independent of the open two-sided `ci_method` default.
- Fixed sentences (verbatim; polish is a Larry-side edit later):
  - "Table 8 of Guo and He (2021) reports an upper confidence bound on the hazard-ratio scale; this is the identical one-sided lower-bound construction under their Section 4 convention βᵢ = −log HR, applied to the argmax-selected subgroup, and does not constitute an upper-bound capability for a second, complementary subgroup."
  - "The Guo–He correction is defined through the maximum functional over the supplied candidate family. The complement Ĥᶜ is not a member of that family and is not the argmax of any functional, so no analogous correction exists within their framework; the absence of a complement bound, and hence of a joint two-subgroup claim, is structural rather than an implementation gap."

## 7. Commit plan, order, and STOP list

- Commit order (each by explicit named paths only): (1) T0 task doc; (2) Stage 0 record; (3) T1 driver; (4) Stage 1 probe note appended to the Stage 0 record; (5) Gate 1a record; (6) six T1 bundles + T1 REPORT; (7) T2 driver (authored, not executed); (8) T3 qmd; Phase B commits on Linux follow the same pattern for Gate 1b, six T2 bundles, and the T2 REPORT.
- STOP-and-report, never improvise: any missing input (no fetch); Q2 defaults contradiction; any `identical()` failure in a probe or production pairing check; any needed `R/` change; Gate 1a/1b envelope exceedance; any interaction with the seven untracked files; anything that would touch a closed line.
- Reporting style: all records and the end-of-phase summary to Larry in short bullets, one item per bullet, with verbatim numbers.
