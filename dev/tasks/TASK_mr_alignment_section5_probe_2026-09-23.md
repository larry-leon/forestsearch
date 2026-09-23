# TASK — Section 5 probe: does the MR alignment move the simulation results?

**Repo:** `larry-leon/forestsearch`, branch `feature/glm-extension` — run CC from the forestsearch clone.

**Purpose.** The MR admission alignment (`57807ec2` and preceding) moved the GBSG application's forest-search
numbers by 0.003 to 0.004 on the hazard-ratio scale. The manuscript's Section 5 reports 18 cells at 2,000
replicates. This task measures the alignment's effect on **one** cell — the one where it has the most room —
so that re-running the rest is a decision with a number behind it rather than a guess.

**Cell:** the **attenuated-benefit** cell at **n = 500, prevalence 12.4%** (narrow). Chosen because every
candidate is null there, so whichever candidate wins is a marginal one, and FS's declaration rate in that cell
is the design's lowest (0.624), which is the regime where the admitted set's composition matters most. A harm
cell would be a weaker probe: true harm candidates sit well clear of the admission threshold and win whatever
the threshold is.

**Identifier: forest search only.** GRF and DINA carry no consistency screen, and their MR output was shown
identical before and after on GBSG (Gate 2 of the alignment task). Running them here would double the compute
for no information.

**Kind:** simulation. **No `R/` change of any kind.**
**Compute: a separate go.** Stop at Gate T and report the measured wall clock, replicate count and cell count
before running anything at scale. **Hard abort at 3 h.**

---

## 0. First action

1. Copy this file to `dev/tasks/TASK_mr_alignment_section5_probe_2026-09-23.md`; `git add` that path; commit
   `docs(tasks): Section 5 probe for the MR admission alignment (2026-09-23)`.
2. Record HEAD, installed forestsearch version and build date, R version, platform, worker count.
3. `git status --short`: record pre-existing untracked files; never stage them.

No `git fetch` / `pull` / `push`. Explicit paths on every `git add`.

---

## 1. Locate the driver and the committed payload — from source

Find, in this repo, the driver that produced Section 5's survival design, and the committed payload for the
target cell. **Reuse the driver unchanged**; do not write a new one. A comparison against a reconstruction of
the pipeline is not a comparison against the published numbers.

Report: the driver path and commit, the payload path and commit, the cell's identifier in the payload, its
replicate count, and its seed convention.

**Then answer the question that decides this task's cost:**

> Does the committed payload retain **per-replicate** rows — corrected estimate, field lower bound, field-s
> upper bound, Bonferroni pair, declaration indicator — or only cell-level summaries?

- **If per-replicate rows are retained:** the "before" arm already exists. Run only the "after" arm, on the
  same seeds, and diff replicate by replicate. This is the preferred route: it halves the compute and does not
  re-run work that is already committed.
- **If only summaries are retained:** both arms must be run. The "before" arm runs on a clean export of the
  pre-change tree (§2).

State which route applies and why.

---

## 2. The pre-change tree (only if §1 says both arms are needed)

Determine the pre-change commit from `git log` — the last commit whose tree still has the exact cutoff in
`R/fs_mr_inference.R`; do not assume which one it is. Export it with `git archive` to a scratch directory
outside the repo and run it with `devtools::load_all()`.

Confirm from source that the exported tree contains `pconsistency.digits` on `forestsearch()` (landed earlier,
`06ac5391`) so the two arms differ **only** in MR's admission threshold.

---

## 3. Gate T — timing, then stop for the compute go

Run **10 replicates** of the cell on the production worker count, on whichever arms §1 requires.

Report: wall clock for 10, the projection to 200, the worker count, and the per-replicate mean. **Then stop.**
Do not proceed to the full probe without Larry's go.

---

## 4. The probe (after the go)

**200 replicates**, paired on identical seeds. Not 2,000: the comparison is paired, so the per-replicate
difference has far smaller variance than either arm alone, and 200 is enough to characterise a shift of the
size GBSG showed. The aim is to measure the shift, not to reproduce the cell's published summaries.

Per replicate, record for both arms: declaration indicator, region (rule, N), corrected estimate, field lower
bound, field-s upper bound, Bonferroni pair.

---

## 5. Read-out

One table, quantities × {before, after, difference}, then a short plain-language reading. No prose summary.

Report, over the declaring replicates:

- **The per-replicate difference** for each quantity: mean, median, and the 5th / 95th percentiles. The
  distribution of the shift is the result — not a change in coverage, which 200 replicates cannot resolve and
  which this task does not attempt to measure.
- **The share of replicates whose selected region changed at all** (rule or N). This is the replicate-level
  counterpart of the alignment task's 26-in-5,000 winner-change count.
- **The declaration rate in both arms.** It must be identical: declaration is decided by the search, and MR
  runs afterwards. A difference means the change reached the search and is a Gate 2 failure.
- **The largest single-replicate shift observed**, with that replicate's index and its region.

Then state plainly whether the shift is small relative to the distance between the bounds and the thresholds
they are read against, and therefore whether re-running Section 5 is a refresh or a revision. If the answer is
"refresh", say so; that is the useful finding, not a disappointing one.

---

## 6. Gates (stop on failure)

- **Gate T (§3):** timing reported, then stopped for the go.
- **Gate 1 — the before arm is the published pipeline.** If §1's route reuses the committed payload, no check
  is needed. If both arms are run, the before arm's first 10 replicates must match the committed payload's
  corresponding rows, where those rows exist. Report the comparison; a mismatch means the driver or the seeds
  are wrong and everything downstream is void.
- **Gate 2 — the search is untouched.** Declaration rate, selected region and N are identical between arms on
  every replicate where the region did not change through re-selection. Any difference in the *declaration*
  indicator is a failure outright.
- **Gate 3 — no `R/` file is modified.** `git status` confirms it at the end.

---

## 7. Record and commits

Report to `dev/reports/REPORT_mr_alignment_section5_probe_2026-09-23.md`.

Commits, explicit paths, in order: task doc; the probe payload and logs; the report. Follow the simulation
directory's own catalogue rule if the payload lands in a catalogued directory.

---

## POST-CONDITIONS (machine-checkable)

1. Driver and payload located in this repo, with paths and commits recorded; the driver is reused unchanged.
2. §1's route stated explicitly, with the reason.
3. Seeds identical across arms and recorded; the pairing is replicate-for-replicate.
4. Gate T reported and the task stopped for the go before any run at scale.
5. Declaration rate identical between arms.
6. Per-replicate difference distributions reported with percentiles, not adjectives.
7. No `R/` file modified; no file written outside the simulation directory, `dev/tasks/` and `dev/reports/`.
8. Errors: 0, or every failure enumerated with its replicate index and the denominator stated.

---

## OUT OF SCOPE

No second cell, no other n, prevalence or planted effect. No GRF or DINA. No coverage estimate — 200 replicates
cannot resolve one, and the shift distribution is what this task measures. No `R/` change. No re-run of the
other 17 cells: that is a decision for after this result. No ACTG 175. No manuscript edits, and no fix to the
consumers that read the removed `pstar_implied` — both separate.
