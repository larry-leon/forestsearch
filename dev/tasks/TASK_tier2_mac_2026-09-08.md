# TASK — Tier 2 (Mac Studio): the finalized constructions at 12.5% prevalence, the dominated-complement regime (campaign `tier2`)

Date: 2026-09-08. Author: chat (spec). Executor: Claude Code on the **Mac Studio**, unattended (Larry offline until the morning). Approver: Larry (kickoff paste = compute go). Reviewer: the Linux MR-field chat.
Runs in parallel with the Linux `cert20` campaign; **outputs are disjoint** (different campaign tag, prevalence, focus and stems) and this task makes **no `R/`, template or other code change**, so the two sessions cannot interact. Branch `feature/glm-extension` (no Mac branch: campaign compute adds files only).
Finalized approach: field one-sided lower on β(Ĥ); field-s (R1) one-sided upper on β(Ĥᶜ); IJ two-term two-sided; Bonferroni joint. Caveat on both blocks: stable-pick (high-p̂) under-correction — reported by stratum, never hidden in averages.
Predecessors: `REPORT_mr_field_complement_2026-09-06.md` and `REPORT_complement_refinements_2026-09-06.md` (the s7c / map1c unscaled-field records at this prevalence: complement upper 0.930–0.956, complement SE ≈ naive in the dominated regime), `REPORT_field_studentize_e1_2026-09-08.md`, `dev/notes/HANDOFF_mr_field_linux_2026-09-08.md`.

## Protocol (Mac specifics in bold)

- First action: **`git pull` on `feature/glm-extension`**, then archive stale variants in `~/Downloads` to `~/Downloads/cc_archive/`, copy this file to `dev/tasks/`, commit. Do not push. If the file is missing from `~/Downloads`, reconstruct it from the kickoff and commit it under this name with a reconstruction note.
- **`devtools::install(dependencies = FALSE)` on this machine** (installs do not carry across machines; never `load_all()`); record `git log -1` and the installed version in the record. The Linux tree may be ahead by commits made tonight — expected; state the HEAD used.
- **No `R/`, template, or shared-document change of any kind.** Every knob is set explicitly in the environment, so the Linux defaults work in flight tonight cannot affect these results.
- **Cross-machine note for the record:** DGM draws are seed-determined and identical across machines (`n_true`, `truth`), but fitted quantities differ at BLAS precision, so any comparison to a Linux-produced bundle uses a stated tolerance (~1e-8), **never `identical()`**. Do not assert byte identity against Linux bundles.
- Gates stop per cell; the driver projects from a 5-replicate smoke **at this machine's worker count** and defers, listing, any cell that would cross the ceiling; `.refuse_if_tracked()` live; records beside results.
- Standing conventions: winner-only and winner-floor excluded; bounds by location; Wilson intervals; marginal and error SDs side by side; verify from source. **No recommendation change; report and wait.**

## Stage 0

Quote from the installed package: the `field_scale_complement` argument of `fs_mr_inference()` and the `lam_cs` line of `.fs_mr_field_complement()` (confirming field-s is present at this HEAD); the template's `FS_S7_FIELD_SCALEC` knob. Record the machine's core count and the worker count chosen (leave headroom; the dense-matrix stages are memory-bandwidth-bound — do not oversubscribe).

## Stage 1 — Smoke and projection

5 replicates at HR 1.75 n500 (below): every construction finite; interval invariants on harm, complement, `_s` and joint blocks; γ in range; realized prevalence ≈ the M1 default (state it); ρᶜ recorded. Time per replicate at the chosen worker count and project all four cells. **Gate 1: proceed if the projection is ≤ 9 h wall; hard timeout 11 h;** defer, listing, any cell that would cross the ceiling (defer order: HR 1.75 n1500 first, then HR 1.00 n500).

## Stage 2 — The four cells (campaign `tier2`)

All cells: `FS_S7_FOCUS=maxeffCons` (M1 default prevalence — **do not set `FS_S7_Z1Q`**, so the committed 12.5% configuration of `s7c` / `map1c` is reproduced), `FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=tier2`, J = 10 (default), `return_reselection = TRUE`, seeds `8316951 + sim_id`, sim_id 1–2000, two seed-disjoint batches of 1,000 then combine.

| # | Cell | Comparator on record (unscaled field) | Purpose |
|---|---|---|---|
| 1 | HR 1.75 n500 | `s7c` h175 n500 | field-s where the complement is dominated |
| 2 | HR 1.75 n1000 | `map1c` if present at this cell (else none) | the complement's large-sample behaviour at low prevalence |
| 3 | HR 1.75 n1500 | `map1c` h150 n1500 is a different HR — none exact | largest complement on record |
| 4 | HR 1.00 n500 | `s7c` null | the null at low prevalence |

**Gate 2 per cell:** completeness (2,000 rows, sim_id 1–2000, no duplicates, no CONFIG-ERROR, meta knobs as set, machine and `forestsearch_version` recorded); realized prevalence (distributional criterion); **same DGM draws as the s7c/map1c comparator where one exists — assert `n_true` identical on all rows and `truth` `identical()`** (these are seed-determined, so they hold across machines); every harm / complement / `_s` / joint / β(Ĥ) / β(Ĥᶜ) / p̂ / K / ρᶜ quantity finite on detected replicates; interval invariants; γ ∈ [0.025, 0.05]; achieved joint probability ≥ 0.95 − 2/n_joint; bound↔quantile identities ≤ 1e-12 (within-bundle, so machine-local); p̂ validity; detection / sens / spec / PPV / NPV / |Ĥ| recorded.

## Stage 3 — Record

`summary_tier2.qmd`, **transplanted from the committed `summary_e1stud.qmd`** (change the globs, cell labels and comparator names only — no fresh authorship), rendered; `REPORT_tier2_2026-09-08.md` beside the results, with the Gate 2 record inside and every number verbatim:
1. Standard tables (per-cell constructions; across cells) for the four cells: rows naive / field / field-s / IJ two-term, both blocks; bias (log; marginal-SD and error-SD units), SDs, SE, r, SE/error-SD, one-sided on the exposed side [Wilson], two-sided.
2. **The dominated-regime check:** ρᶜ (mean, q10, q90, share > 1), λ-SDᶜ/naive SE and `se_field_s`/naive SE, and field vs field-s coverage side by side — the pre-registered expectation is ρᶜ ≈ 1 and field-s ≈ field here; state whether it holds, cell by cell, with Wilson overlap.
3. **Against the committed comparators** (unscaled field 0.930–0.956 at this prevalence): the field's coverage in these cells beside the comparator's, noting the comparator's campaign and that the pairing is by DGM draws, not by fitted values.
4. **By p̂ tertile, both blocks, every cell:** the stable-pick caveat at low prevalence — does the T3 shortfall appear here, and does it shrink with n? Gaussian-implied beside observed.
5. Identification and n: detection, |Ĥ|/|H| median and q90, sens / spec / PPV / NPV, β(Ĥ) and β(Ĥᶜ) against the planted values, p̂ mean, complement fits per replicate.
6. Acceptance criteria evaluated as findings (pre-registered): harm field one-sided lower ≥ 0.94 with Wilson support; field-s complement one-sided upper ≥ 0.93 at every cell of this prevalence (the dominated regime — a higher bar than at 31%); IJ two-sided ≥ 0.93; joint Bonferroni ≥ 0.93; the null cell recorded without criterion. Met / not met / not run, per cell, with the number.

## Done means

Stage 0 quotes and machine record; Stage 1 smoke and Gate 1 projection; cells completed / deferred / dropped with walls; Gate 2 records; `summary_tier2` rendered; the report committed; **branch left unpushed** (Larry pushes; the Linux commits go first, then `git pull` here before pushing); one-paragraph closing summary with the criteria table one line per cell and the commit range. Out of scope: any code change, any recommendation change, any comparison asserting byte identity to Linux-produced bundles.
