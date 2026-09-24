# TASK — Section 5 full re-run under the aligned package

**Repo:** `larry-leon/forestsearch`, branch `feature/glm-extension` — run CC from the forestsearch clone
**on the Mac Studio.**

**Purpose: provenance, not correction.** The 18-cell sweep
(`dev/reports/REPORT_mr_alignment_section5_sweep_2026-09-23.md`) established that the MR admission alignment
moves Section 5's numbers in the fourth decimal: declaration identical on every replicate of every cell, the
selected region never changed, the largest single shift anywhere 0.015 on a bound sitting 0.28 from the nearer
threshold. **No stated result, range or conclusion in the manuscript changes.** What this task fixes is that
someone re-running Section 5 from the shipped package should get the published numbers. That is a provenance
defect even when the numbers are immaterial.

**Scope: forest search only.** GRF and DINA carry no consistency screen and their MR output was shown identical
before and after; their committed campaigns stand as published. This task re-runs the FS campaigns —
`p12x20` (9 cells), `cert20` (7), `e1stud` (2) — and nothing else.

**Kind:** simulation. **No `R/` change of any kind.**

**Compute: a separate go.** Committed wall clocks for the 18 cells sum to about 45,000 s, but at mixed worker
counts and on another machine, so that is an order-of-magnitude figure and not a projection. Step 2 replaces it
with a measurement before anything runs at scale. Expect an overnight, unattended run.

---

## 0. First action

1. Copy this file to `dev/tasks/TASK_section5_full_rerun_2026-09-24.md`; `git add` that path; commit
   `docs(tasks): Section 5 full re-run under the aligned package (2026-09-24)`.
2. Record HEAD, R version, platform, core count, and the library path R loads `forestsearch` from.
3. `git status --short`: record pre-existing untracked files; never stage them.

### 0a. Is the clone current? (stop-on-failure)

**Assert that the alignment commits are in this clone** before installing anything:
`git merge-base --is-ancestor 7713942e HEAD` (the MR admission alignment), and the same for `96f84ad8`
(FŴ / settable p-star) and `06ac5391` (`pconsistency.digits` on `forestsearch()`).

If any is not an ancestor, the clone is behind. **Stop and report**, naming which commits are missing. Do not
pull — CC does no remote git operations; Larry pulls via GitHub Desktop and re-runs this task.

Installing a clone that is behind would install a stale build and the assertions in 0b would pass against the
wrong code, so this check comes first.

### 0b. Install, then verify (stop-on-failure)

With `R_LIBS` unset and `R/`, `DESCRIPTION` and `NAMESPACE` clean:

```
Rscript -e 'devtools::install(dependencies = FALSE, upgrade = FALSE)'
```

`dependencies = FALSE, upgrade = FALSE` keeps it to this package with no dependency chain — the same form used
for the 2026-09-23 install and the `p12x20` Stage 0 install. The version stays `0.3.5.9000`; do not bump it.

Then **verify from inside a `multisession` worker, not only in the main process**: `pconsistency.digits` in
`formals(fs_mr_inference)` and `.fs_pcons_eff` in its body, both TRUE. Record the library path, the installed
version and the build timestamp. If either assertion is FALSE, **stop and report** — running on a stale build
is the hazard the 2026-09-23 work existed to remove.

No `git fetch` / `pull` / `push` from CC. Explicit paths on every `git add`.

---

## 1. Cells and structure — from the committed record

The 18 cells, their source campaigns, parameters and committed wall clocks are in
`quarto/simulations/gbsg_020/mrs5sweep/cells.txt`. Take them from there. Seeds are `8316951 + sim_id`,
`sim_id` 1…2000.

**Use each campaign's own driver and knob set unchanged.** Where a driver's batching does not fit, mirror its
knob set exactly and record the deviation.

**Produce bundles in the same structure as the committed campaign** — same batch/combine layout, same file-name
pattern — so that regenerating the manuscript's figures is a change of campaign tag and nothing else.

**Use a new campaign tag.** Nothing committed is overwritten: the published payloads stay as the record of what
was published, and the re-run sits beside them.

---

## 2. Check and timing, then stop for the go

Run **one cell at 200 replicates** — the cheapest, `A7` — with the production worker count.

- **Compare it to the sweep's committed after arm for that cell** (`mrs5probe` A7, 200 replicates, same seeds).
  This is a configuration check: the two runs used the same code and the same seeds, so the comparison
  establishes the knob set, the seed convention and the cell parameters are right before twelve hours of
  compute follow.
  - **Declaration and selected region must match.** A difference on more than a couple of replicates means the
    configuration is wrong: stop and report.
  - **Numeric differences in the estimates and bounds are recorded, not gated.** Report the maximum absolute
    difference per quantity and move on.
- Report wall clock, worker count, the projection to 2,000 replicates, and the projection to all 18 cells.
- **Then stop** for Larry's compute go. Do not start the full run.

---

## 3. The full run (after the go)

All 18 cells, **2,000 replicates each**, forest search only.

- **Cheapest first**, so a partial run is still usable: n = 500, then 1000, then 1500.
- Per-cell overrun guard at 3× its projection; on overrun, record and continue to the next cell.
- Record per cell: wall clock, worker count, replicate count, error count, declaration rate.
- **A partial run is a successful outcome.** Report what completed and what did not, with the reason.

---

## 4. Read-out

One table, cells × quantities, then a short plain-language reading. No prose summary.

Per cell, against the **committed published bundle**:

- declaration rate, before and after;
- mean and median change in the corrected estimate, the field lower bound and the field-s upper bound, with the
  5th / 95th percentiles;
- the count of replicates whose selected region changed;
- the count of bound readings that cross 0.75 or 1.25.

Then state in a few sentences whether any published Section 5 figure or range changes at its reported
precision. The sweep predicts not; say plainly whether the full run bears that out.

---

## 5. Catalogue and commits

Follow `gbsg_020`'s own closeout rule: the curated entry and inventory rules in `status_curated.md` and the
`Updated` line, then the regenerated `current_status.md` committed alone as a child commit, with the pin
asserted equal to HEAD at commit time.

Report to `dev/reports/REPORT_section5_full_rerun_2026-09-24.md`.

Commits, explicit paths, in order: task doc; payloads and logs; report; catalogue entry; regenerated
`current_status.md`.

---

## POST-CONDITIONS (machine-checkable)

1. Step 0a: the three alignment commits asserted as ancestors of HEAD. Step 0b: the install run, and the
   aligned build verified from inside a worker, with the library path, version and build timestamp recorded.
2. Cells, parameters and seeds taken from `mrs5sweep/cells.txt`; seeds are `8316951 + sim_id`.
3. A new campaign tag is used; no committed payload is modified or overwritten.
4. Step 2's comparison reported, with declaration and region matching.
5. Bundle structure matches the committed campaign's, so a figure re-run is a tag change only.
6. Errors: 0 per cell, or every failure enumerated with its replicate index and the denominator stated.
7. Cells not run are listed with the reason; the report states whether the run is complete or partial.
8. No `R/` file modified.
9. Catalogue pin equals HEAD at commit time.

---

## OUT OF SCOPE

No GRF or DINA. No figure or manuscript regeneration — that belongs to the fs-glms-interpretable chats, and
this task produces the payloads they will read. No `R/` change, including the OC predictor's exact-cutoff gate
(`R/fs_oc_predict.R:279`, `R/fs_oc_grid.R:582`), which is recorded and awaiting its own decision. No declcal
consumer fixes. No ACTG 175.
