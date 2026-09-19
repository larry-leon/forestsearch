# CC TASK — read-only exposure check: the four-cell existence condition against committed work

**Opened:** 2026-09-18 · **Repository:** forestsearch · **Authorized by:** Larry, 2026-09-18.
**Workstream:** admission floors (new register; the thresholds workstream is closed and untouched).
**Type:** read-only. **No edit to `R/`, no `forestsearch()` call, no bootstrap, no fold, no install, no push.**
**Compute:** data regeneration from the committed recipes plus cell counting — no model is fit anywhere.
**Hard abort 30 minutes** total; if exceeded, stop and report what was covered.
**Testing:** none (nothing is changed).

## Why

The next package task makes the estimator boundary return NA-with-reason instead of a finite divergent
estimate when the analysis estimand does not exist on a candidate:
- binary: any of the four cells (control events, control non-events, treated events, treated non-events) = 0;
- survival: a zero-event arm.

Larry's constraint: no committed result may silently move. This check proves invariance before the rule
lands: **if every declared subgroup in every committed campaign already satisfies the condition, the rule
could not have changed any recorded estimate.**

## Method — reuse the admission check's machinery, including its correction

`quarto/simulations/actg175/binary_020/REPORT_binary_admission_check_2026-09-18.md` §4 recomputed declared
memberships from recorded `sg_def` with `forestsearch:::.fs_resolve_membership()` and reproduced `n_sel`
exactly on all 400 rows. Reuse that method. **Reuse its recorded correction too (its finding 6):** regenerate
data exactly as each campaign's template does — the pre-generated seed table by global `sim_id`, the RNG kind
set per replicate, the DGM calibrated *before* any kind switch — and **assert the recomputed membership size
equals the recorded `n_sel` on every row.** A mismatch means the wrong population; stop and report rather
than count cells on it.

## Scope — every committed campaign

For each committed cell, every declared replicate, both blocks (Ĥ and Ĥᶜ where recorded):

| Campaign family | Outcome | Cells to compute | Condition |
|---|---|---|---|
| `orfs_*` (the two committed binary FS cells) | binary | the four cells per arm | all four ≥ 1 |
| survival campaigns behind the three-identifier grid (FS, DINA `dinamr`, GRF `grfmr`, Part A `p12x20`) | survival | events per arm | both ≥ 1 |
| continuous campaigns (`mdsgnb20`, `mdgrf`, `mddina`) | continuous | subjects per arm | report only — no existence condition applies; note any arm < 2 |

Locate every committed campaign by search (`quarto/simulations/**`), not from this list — if one exists that
is not named here, include it and say so. If a campaign's artifacts do not record `sg_def` or enough to
recompute membership, report that cell as **not checkable** rather than skipping it silently.

## Deliverable

`dev/reports/REPORT_four_cell_exposure_2026-09-18.md`, opening with the HEAD SHA, containing:

1. **The exposure table** — one row per campaign × cell × block: replicates checked, the minimum of each
   cell across replicates, and the **count of replicates violating the condition**. Membership-size
   assertion result per row.
2. **The verdict**, in one line: total violating replicates across all committed work. Zero means the
   boundary rule is provably invariant on committed results. Any non-zero row is listed with `sim_id`s.
3. **The residual assumption, stated:** this checks *declared* subgroups. A non-estimable candidate that
   was admitted but not selected is not visible here; under the neighbourhood selection rule a divergent
   effect would have dominated selection and appeared as the declared Ĥ, which is why the declared check is
   the right instrument — but say so, do not imply more.
4. **Not-checkable cells**, if any, with the reason.
5. Wall clock.

Findings only; no fix, no task attached. Commit the report. **Do not push.**
