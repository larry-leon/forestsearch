# CC TASK — binary campaign: pre-launch, stage-1 launch (FS), stage-2a timing

**Opened:** 2026-09-18 · **Repository:** forestsearch · **Authorized by:** Larry, 2026-09-18.
**Study:** `quarto/simulations/actg175/binary_020/`, design of record from
`REPORT_binary_redesign_2026-09-18.md` (prevalence 14.917%, `sg_quantile = 0.62850`). **No edit to `R/`.**

**Fixed by Larry for this campaign:**
- **Replicates: 1,000 per cell** — replaces the study's earlier 2,000 everywhere it is configured; recorded in
  every bundle's meta.
- **Workers: 63** — the study's configured count (64 physical cores, one left for the parent process). All
  wall-clock estimates below are calibrated at 63. Do not change it.
- **MR only.** The evaluated set is unadjusted, oracle, IJ two-term, the field and field-s bounds, and
  Bonferroni. **No bootstrap and no cross-validation anywhere in the campaign** — asserted in Part 1.

This task has two parts with a **hard stop between them**. Part 1 makes commits and reports; Larry pushes;
Part 2 launches only on a second explicit kickoff.

---

## Part 1 — pre-launch (commits; the only compute is one Stage 0 re-render, cap 5 minutes)

1. **Version bump.** The installed library is HEAD's code labelled 0.3.5 — the same version the superseded
   cells recorded. Bump `DESCRIPTION` to the repository's development-version convention (`0.3.5.9001`
   unless the repo's convention differs — check `NEWS.md`'s header and prior bumps). Commit alone.
2. **Reinstall — authorized by Larry for this launch preparation only:** `devtools::install(quick = TRUE)`,
   then verify `packageVersion("forestsearch")` in a **fresh** R session equals the bumped version. Record
   the build timestamp and version in the study's status file, beside the note that the superseded cells ran
   on the 2026-09-17 build of 0.3.5.
3. **Revert the sweep driver's helper edit** (Larry's disposition: the historical driver stays byte-identical
   to what produced its cells). Restore that file to its content at the commit before `625b10d0`; assert
   byte-identity with `git diff --stat` against that revision; commit alone. The template of record keeps the
   new helper. Update the redesign report's helper line to say one in-scope copy, with the reason.
4. **`smoke_identity.R` `recipe` mode:** record as deferred (retire or re-point after the first new-design cell
   completes) in the status file. No code change.
5. **Replicates → 1,000** in the template of record / runner configuration, wherever the per-cell count
   lives; assert the seed table indexed by global `sim_id` covers 1,000 per cell.
6. **MR-only assertion.** From the template of record and runner configuration, assert that no bootstrap
   argument and no cross-validation is enabled in any cell (every such argument off or absent), and that
   the estimator set is exactly the six named above. Print the assertion result.
7. **Workers assertion.** Confirm the runner is configured at 63 workers; print it.
8. **Stage 0 re-render** once, to confirm the configuration and the feasibility gate still pass.
9. **Report Part 1 in the status file and STOP.** Print the pre-launch checklist: version, build time,
   driver reverted, replicates = 1,000, MR-only assertion, workers = 63, Stage 0 green. **Do not launch.**

---

## Part 2 — stage-1 launch: the six FS cells (compute, on the second kickoff only)

**Go/no-go as approved by Larry:** 6 cells (`orfs`: OR 0.75 / 1.00 / 1.50 × n 500 / 2000) × 1,000
replicates, ≈ **5.1 h at 63 workers** by linear scaling from the measured 2,000-replicate cells
(n = 500 ≈ 1,124 s; n = 2000 ≈ 4,993 s). Unattended.

- **Order:** `orfs_or150_n500` **first** — the configuration that halted twice — then the other n = 500
  cells, then the three n = 2000 cells.
- The study's existing runner, heartbeat, `HALT_*` and Gate 2 conventions apply unchanged. Stop-on-failure
  per cell: a failed cell writes its halt record; the run continues to the next cell only if the runner
  already does so. No new retry logic.
- **Per-cell record must carry:** package version and build time, prevalence and `sg_quantile`, replicates,
  workers, measured wall clock, Gate 2 result, and the count of non-estimable candidates and NA oracle
  replicates (both now reported by the boundary rule and the four-cell helper).
- **Stage-1 report** (`REPORT_binary_stage1_fs_2026-09-18.md`, in the study directory): the six cells' Gate 2
  results, measured wall clock per cell versus the linear-scaling estimate, and the non-estimable / NA counts
  — the first evidence of what the new design and boundary rule produce at scale.

## Part 2a — stage-2 timing, immediately after stage 1 (compute, same kickoff)

Run exactly **two** cells to measure the unmeasured identifiers: `orgrf_or150_n500` and
`ordina_or150_n500` (cheapest n, and the previously-halting design point) at 1,000 replicates, **hard cap
3 h each** — abort and record if exceeded. Report measured wall clock per cell and project the remaining 10
GRF/DINA cells from it. **Then STOP.** The stage-2 launch is Larry's separate go/no-go on those measured
numbers.

## Not in this task

Any change to floors, thresholds, the DGM of record, the worker count, or the package. The remaining 10
GRF/DINA cells. Pushing.
