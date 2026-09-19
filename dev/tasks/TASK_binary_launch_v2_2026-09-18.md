# CC TASK — binary campaign launch, v2 (single authoritative document; idempotent)

**Opened:** 2026-09-18 · **Repository:** forestsearch · **Authorized by:** Larry, 2026-09-18.
**Supersedes every earlier copy of `TASK_binary_launch_stage1_2026-09-18.md`.** Delete any earlier copy from
`dev/tasks/` and commit this one in its place. **This document is safe to run whatever state the repository is
in** — every step first checks whether it is already done and skips if so.
**Study:** `quarto/simulations/actg175/binary_020/`; design of record from
`REPORT_binary_redesign_2026-09-18.md` (prevalence 14.917%, `sg_quantile = 0.62850`). **No edit to `R/`.**
**Unattended. Do not push at any point.**

**Fixed by Larry:** replicates **1,000** per cell · workers **63** (all timings calibrated there; do not change)
· **MR only** — the evaluated set is unadjusted, oracle, IJ two-term, field, field-s, Bonferroni; no bootstrap,
no cross-validation anywhere.

---

## Step 0 — discover the current state (read-only) and print it as a table

Before doing anything, establish and print:

| Item | How to check |
|---|---|
| Is a campaign run active? | running `R`/`Rscript`/`quarto` processes; heartbeat and `HALT_*` files in the study directory |
| Part 1 commits present? | `git log` for the version bump, the driver revert, the replicates change |
| Installed version | `packageVersion("forestsearch")` in a fresh R session, and `DESCRIPTION`'s version |
| Driver reverted? | the sweep driver byte-identical to its content at the commit before `625b10d0` |
| Replicates configured at 1,000? | the template of record / runner configuration |
| MR-only? | no bootstrap or CV argument enabled in any cell; estimator set is the six named |
| Workers = 63? | runner configuration |
| Stage 0 green on the new design? | the most recent Stage 0 record |

**If a stage-1 run is already active: do not launch anything. Monitor it to completion, then continue from
Step 3's report onward.** Otherwise proceed.

## Step 1 — complete whatever pre-launch items are missing (skip each one that Step 0 found done)

1. Version bump in `DESCRIPTION` to the repo's dev convention (`0.3.5.9001` unless prior bumps say otherwise);
   commit alone.
2. Reinstall — **authorized by Larry for this launch preparation only** — `devtools::install(quick = TRUE)`;
   verify `packageVersion()` in a fresh R session; record version and build time in the status file beside the
   note that the superseded cells ran on the 2026-09-17 build of 0.3.5.
3. Revert the sweep driver's helper edit to its pre-`625b10d0` content; assert byte-identity; commit alone.
   The template of record keeps the new helper; update the redesign report's helper line accordingly.
4. Record `smoke_identity.R` `recipe` mode as deferred in the status file (no code change).
5. Replicates → 1,000 wherever configured; assert the global-`sim_id` seed table covers 1,000 per cell.
6. Assert MR-only (no bootstrap, no CV; the six estimators) and workers = 63; print both.
7. Stage 0 re-render once on the new design (cap 5 minutes) — configuration and feasibility gate green.

## Step 2 — the gate

Print the pre-launch checklist: version bumped and installed · driver reverted · replicates 1,000 · MR-only ·
63 workers · Stage 0 green. **Proceed only if every item is green. If any item fails: STOP, leave the failure
in the status file, launch nothing.**

## Step 3 — stage 1: the six FS cells

6 cells (`orfs`: OR 0.75 / 1.00 / 1.50 × n 500 / 2000) × 1,000 replicates, 63 workers ≈ 5.1 h by linear
scaling from the measured 2,000-replicate cells. **Order:** `orfs_or150_n500` first (the configuration that
halted twice), then the other n = 500 cells, then the n = 2000 cells. Existing runner, heartbeat, `HALT_*`
and Gate 2 conventions unchanged; no new retry logic. Every bundle's meta carries version, build time,
prevalence, replicates, workers, wall clock, and the non-estimable / NA-oracle counts.

**Stage-1 report** `REPORT_binary_stage1_fs_2026-09-18.md` in the study directory: Gate 2 per cell, measured
wall clock per cell against the estimate, and the non-estimable / NA counts.

## Step 4 — stage-2a timing, then STOP

`orgrf_or150_n500` and `ordina_or150_n500` only, 1,000 replicates, **3 h cap each** (abort and record if
exceeded). Report measured wall clock per cell and project the remaining 10 GRF/DINA cells. **Then STOP** —
the stage-2 launch is Larry's separate go/no-go on those numbers.

## Not in this task

Any change to floors, thresholds, the DGM of record, the worker count, or the package. The remaining 10
GRF/DINA cells. Pushing.
