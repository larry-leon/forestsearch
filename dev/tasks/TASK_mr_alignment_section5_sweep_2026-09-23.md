# TASK — Section 5 sweep: the MR alignment's shift across all 18 cells

**Repo:** `larry-leon/forestsearch`, branch `feature/glm-extension` — run CC from the forestsearch clone.

**Extends** `dev/tasks/TASK_mr_alignment_section5_probe_2026-09-23.md` (`ebfc4548`) from one cell to the whole
survival design. Its OUT OF SCOPE line forbidding a second cell is **superseded by this task**. Everything else
in it — the after-arm-only route, forest search only, the paired-on-seeds design, the read-out shape — carries
over unchanged.

**Purpose.** Measure the alignment's per-replicate shift in every Section 5 cell, so that whether Figures 3 and
4 and Sections 5.2/5.3 need re-running is settled by measurement across the design rather than extrapolated
from one cell.

**Route, already established.** The committed campaign payload **is** the before arm: a clean export of the
fix's parent (`ba595f4b`) reproduced the committed A7 rows exactly, 163 of 163 columns, maximum difference 0.
So only the after arm runs, and nothing committed is re-run.

**Kind:** simulation. **No `R/` change of any kind.**

**Compute (Larry's go, 2026-09-23, given for unattended running while away):** 18 cells, **200 replicates
each**, forest search only, after arm only. Measured basis: A7 at n = 500 took 61 s for 10 replicates on 64
workers, projecting 2–3 min for 200. Cells at larger n and higher declaration rates cost more; the per-cell
timing in §2 replaces the projection with a measurement as the sweep proceeds. **Hard abort at 2.5 h for the
whole task.**

---

## 0. First action

1. Copy this file to `dev/tasks/TASK_mr_alignment_section5_sweep_2026-09-23.md`; `git add` that path; commit
   `docs(tasks): Section 5 sweep for the MR admission alignment (2026-09-23)`.
2. Record HEAD, R version, platform, worker count, and the library path the workers load.
3. `git status --short`: record pre-existing untracked files; never stage them.

No `git fetch` / `pull` / `push`. Explicit paths on every `git add`.

---

## 1. Enumerate the cells from the committed payload

List all 18 cells **from the committed p12x20 record**, not from this document: their identifiers as the
payload names them, prevalence, planted effect, n, committed replicate count, committed declaration rate and
committed wall clock. Do not rely on any cell label or declaration rate quoted in chat — one was wrong.

If A7's after arm was already run under the probe task, reuse it rather than re-running.

---

## 2. Run order, and the budget rule

**Run cheapest first**, so that if the budget runs out the completed cells are still a usable answer:
all n = 500 cells, then n = 1000, then n = 1500; within each, attenuated-benefit before the harm cells.

Before each cell, project its cost from the committed wall clock scaled to 200 replicates and the worker count.
After each cell, record the measured wall clock.

- **If the running total plus the next cell's projection would exceed 2.5 h, stop cleanly**: do not start that
  cell. Report what completed and what did not. **A partial sweep is a successful outcome, not a failure** —
  report it as such.
- If any single cell overruns its projection by more than 3×, stop that cell, record it, and continue with the
  next.

---

## 3. Per cell

200 replicates, seeds identical to the committed rows, paired replicate-for-replicate. Record per replicate:
declaration indicator, region (rule, N), corrected estimate, field lower bound, field-s upper bound,
Bonferroni pair.

Use a campaign tag that cannot overwrite anything committed, as the probe's Gate T run did.

---

## 4. Read-out

**One table, cells × quantities**, followed by a short plain-language reading. No prose summary, no per-cell
narrative.

Columns per cell, over the declaring replicates:

- mean and median per-replicate shift, and the 5th / 95th percentiles, for the corrected estimate, the field
  lower bound and the field-s upper bound;
- the largest single-replicate shift, with its replicate index;
- the share of replicates whose selected region changed (rule or N);
- declaration rate in both arms — **these must be identical**, since declaration is decided by the search and
  MR runs afterwards.

Then, in a few sentences:

- whether the shift varies systematically with n, prevalence or planted effect, and in which direction;
- whether the largest shift anywhere in the design is small relative to the distance between the bounds and
  the thresholds they are read against (a hazard ratio of 0.75 and of 1.25);
- a plain statement of whether re-running Section 5 at full replicates is a refresh or a revision.

Do not estimate coverage. 200 replicates cannot resolve a coverage change, and the shift distribution is what
this task measures.

---

## 5. Gates (stop on failure)

- **Gate 1 — the search is untouched.** In every cell, the declaration indicator matches the committed rows on
  every replicate. Any difference is a failure outright: it would mean the change reached the search.
- **Gate 2 — the build is the right one.** The workers load the HEAD build, not the installed one, which is
  stale (it predates `7713942e` and has no `pconsistency.digits` in `fs_mr_inference`). Assert this per cell
  from inside a worker, and record the library path. **Silently running the stale build would invert this
  task's answer**, so this gate is not a formality.
- **Gate 3 — no `R/` file modified.** `git status` confirms at the end.
- **Gate 4 — errors.** 0 per cell, or every failure enumerated with its replicate index and the denominator
  stated.

---

## 6. Record and commits

Report to `dev/reports/REPORT_mr_alignment_section5_sweep_2026-09-23.md`.

Commits, explicit paths, in order: task doc; the sweep payloads and logs; the report; then the catalogue line
and its regenerated file if the payloads land in a catalogued directory.

---

## POST-CONDITIONS (machine-checkable)

1. Cells enumerated from the committed payload, not from this document.
2. Seeds identical to the committed rows; pairing is replicate-for-replicate.
3. Gate 2 asserted per cell, with the library path recorded.
4. Declaration rate identical to the committed rows in every completed cell.
5. Shift distributions reported with percentiles, not adjectives.
6. Cells not run are listed with the reason (budget or overrun), and the report states the sweep is partial.
7. No `R/` file modified; nothing written outside the simulation directory, `dev/tasks/` and `dev/reports/`.

---

## OUT OF SCOPE

No full 2,000-replicate re-run — that is a decision for after this result. No GRF or DINA: their MR output was
shown identical before and after on GBSG. No coverage estimate. No `R/` change, and no `devtools::install()`:
the installed build stays as it is, and refreshing it is Larry's decision. No ACTG 175. No manuscript edits.
