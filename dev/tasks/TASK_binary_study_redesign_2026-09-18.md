# CC TASK — binary study: feasible design of record, oracle four-cell condition, Stage 0 gate

**Opened:** 2026-09-18 · **Repository:** forestsearch · **Authorized by:** Larry, 2026-09-18.
**Workstream:** admission floors — study-side task. **No edit to `R/`.** No campaign launch: this task
stops before Stage 2. Launch is a separate go/no-go with cell count, replicates and wall clock.
**Where:** `quarto/simulations/actg175/binary_020/` — the template of record
(`sim_fs_mr_field_or_template.qmd`), its helpers, Stage 0 machinery, and the DGM of record.
**Compute:** `fs_dgm_feasibility()` runs (`n_rep = 200`, seconds each) and the existing Stage 0 smoke on the
new design. **Hard abort 20 minutes total.** No fits beyond the Stage 0 smoke.
**Testing:** the template's own gates and smoke; no package tests are touched.
**Prerequisite at HEAD:** `fs_dgm_feasibility()` and the estimator boundary from
`REPORT_estimability_boundary_2026-09-18.md` — verify both exist before starting.

## Larry's design decision

The current design plants the subgroup at ~10% prevalence, which the feasibility check (and the admission
check before it) shows is undeclarable in ~90% of replicates at n = 500 under the strict 60-subject floor.
**Larry directs the planted prevalence be raised into the 12–15% range.** The floors themselves do not move.

## Step 1 — choose the design from the table, not by hand

1. Run `fs_dgm_feasibility()` for every combination of prevalence ∈ {12%, 12.4%, 13%, 14%, 15%}, n ∈ the
   study's n grid (500, 750, 1000, 2000 — confirm from the template), and each OR design point (0.75, 1.00,
   1.50 — confirm), at `tolerance = 0.05`, `n_rep = 200`, with the DGM calibrated **before** any RNG-kind
   switch, exactly as the template does. Record the full table.
2. **Selection rule (Larry's range, applied mechanically):** among prevalences feasible at **every** n in
   the grid, prefer 12.4% (the survival three-identifier grid's prevalence, for cross-design comparability);
   otherwise the smallest feasible. **If no prevalence in the range is feasible at every n, STOP** and
   present the table with options (drop n = 500 from the grid, or go above 15%) — Larry decides.
3. Recalibrate the DGM of record at the chosen prevalence (`calibrate_glm_interaction()`, same seed
   discipline). Commit it as the study's DGM of record, with the feasibility table beside it.

## Step 2 — supersede the old cells, don't delete them

The two committed cells (`orfs_or075_n500`, `orfs_or075_n2000`) were produced under the old design. Mark
them **superseded by design change** in the study's status/record file, with the old prevalence and the
new one named. No bundle is deleted; nothing is re-run in this task.

## Step 3 — the oracle helper

`.logit_or_ci()` (template and sweep helper — find every copy by search): **keep the legacy pooled 5/5
exactly as it is, and add the four-cell condition** — return the NA quadruple if any of control events,
control non-events, treated events, treated non-events is 0. Both conditions apply; either failing gives
NA. The study's existing non-convergence convention already counts NA per estimator. One helper, identical
in every copy — assert the copies match.

## Step 4 — Stage 0 gate

In the template's Stage 0: call `fs_dgm_feasibility()` on the DGM of record at the study's n grid with
`tolerance` passed **explicitly**; print the table into the Stage 0 record; **refuse to proceed to Stage 2
unless `feasible` is TRUE or an explicit override flag is set**, with the override written into the
record. Same pattern as the DINA caps: possible deliberately, never silently.

## Step 5 — Stage 0 smoke on the new design

Run the study's existing Stage 0 smoke (the `smoke_*` machinery) on the new DGM of record with the new
helper, within the compute cap. Gate: smoke green, the assertion check passes, the feasibility table shows
`feasible = TRUE` at every n.

## Report

`quarto/simulations/actg175/binary_020/REPORT_binary_redesign_2026-09-18.md` (beside the study's other
reports): the full feasibility table, the chosen prevalence and why, the superseded-cells note, the helper
diff, the Stage 0 gate, the smoke result, and — for the launch go/no-go that follows — **the projected
campaign size: 18 cells (3 identifiers × 6), replicates per cell, and the wall-clock estimate from the
study's own prior runs.** Findings only beyond that; no task attached. Commit. **Do not push. Do not
launch.**
