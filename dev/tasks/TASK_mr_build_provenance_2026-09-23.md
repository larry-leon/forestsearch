# TASK — Build provenance: the MR alignment fits, and the cert20 before arm

**Repo:** `larry-leon/forestsearch`, branch `feature/glm-extension` — run CC from the forestsearch clone.

**Run after** the A7 probe, the Section 5 sweep and the declcal consumers task have completed.

**Purpose.** Two questions about which build produced which numbers.

1. `devtools::load_all()` does not propagate to `multisession` workers — they resolve `forestsearch` from the
   installed library, which is the 02:32 UTC build and predates `06ac5391`, `96f84ad8` and `7713942e`. The six
   before/after fits of `REPORT_mr_admission_alignment_2026-09-23.md` ran on 102 multisession workers, and that
   report does not state which build the workers loaded. If they loaded the stale one, the before/after table
   measures a main-process-aligned / worker-stale mixture rather than the alignment.
2. The `p12x20` committed payload was established as a valid before arm by exporting `ba595f4b` and
   reproducing cell A7 exactly (163/163 columns, max diff 0). The **`cert20`** family (31% prevalence) supplies
   the other nine cells of the sweep and has had no such check.

**Kind:** verification, with a conditional re-run. **No `R/` change. Do not run `devtools::install()`** —
Larry's instruction; the installed build stays as it is.

**Compute:** §1 is read-only. §3 is one ~61 s run. §2's re-run, only if needed, is six fits at ~47 s. **Hard
abort at 45 min.**

---

## 0. First action

1. Copy this file to `dev/tasks/TASK_mr_build_provenance_2026-09-23.md`; `git add` that path; commit
   `docs(tasks): build provenance for the MR alignment fits and cert20 (2026-09-23)`.
2. Record HEAD, R version, platform, and the installed forestsearch version and build date.
3. `git status --short`: record pre-existing untracked files; never stage them.

---

## 1. Where MR's admission actually runs — from source

Determine and **quote the dispatch**: does the MR re-selection — the code that builds the admission floor in
`fs_mr_inference.R` and applies it over the field draws — execute in the calling process, or is it dispatched
to `future` / `doFuture` workers?

The candidate search and the MR field draws may be parallelized differently. Report each separately. The
candidate search is unaffected by the alignment, so stale workers there would be harmless; what matters is the
process that evaluates the admission floor.

---

## 2. The six MR alignment fits

**A hypothesis to test, not to assume.** Stale code has no `pconsistency.digits` argument, so it would ignore
the setting and apply the exact cutoff for every value. If stale code had computed the threshold, the
`digits = 6` supplementary fit and the `digits = 2` after fit would both equal the baseline and there would be
no before/after difference at all. A difference was observed, and `digits = 6` reproduced the baseline. Confirm
or refute that this reasoning holds, from what §1 establishes.

Then:

- **If §1 shows the admission floor is evaluated in the calling process**, and that process ran the modified
  code, the before/after table stands. Record the reasoning and the evidence; no re-run.
- **If it is dispatched to workers**, establish from the run logs or the session's library setup which build
  those workers loaded. If it is not established, **say so plainly** rather than inferring it from the results
  having changed — a mixture would also produce a change.
- **If the workers loaded the stale build**, the before/after table is void. Re-run the six fits with the
  worker build asserted from inside a worker, using the HEAD scratch library via `R_LIBS` as the probe does.
  Report the corrected table and mark the original **superseded** in both reports.

---

## 3. The cert20 before arm

Repeat, for one committed `cert20` cell, the check that validated `p12x20`:

- Export the fix's parent (`ba595f4b`) with `git archive` to a scratch directory; run it with
  `devtools::load_all()`.
- Run 10 replicates of that cell at the committed seeds.
- Compare **every column** against the committed rows and report the count matched and the maximum difference.

- **If identical**, `cert20`'s committed payload is a valid before arm and the sweep's 31% cells stand.
- **If not**, report which columns differ and by how much. The sweep's `cert20` **after**-arm data remains
  valid; only the comparison lacks a before arm, and running one is a separate decision — do not run it here.

---

## 4. Record and commits

Report to `dev/reports/REPORT_mr_build_provenance_2026-09-23.md`.

If §2 requires a re-run, amend `REPORT_mr_admission_alignment_2026-09-23.md` to mark its table superseded and
point to the corrected one. Do not rewrite history; amend forward with a new commit.

Commits, explicit paths, in order: task doc; any re-run payloads; the report; any amendment to the earlier
report.

---

## POST-CONDITIONS (machine-checkable)

1. §1's dispatch is quoted from source, with the candidate search and the MR field draws reported separately.
2. §2 reaches one of the three stated dispositions, and says which — no inference from "the results changed".
3. If a re-run was needed, the worker build is asserted from inside a worker, and both reports record the
   supersession.
4. §3 reports the column count matched and the maximum difference for the `cert20` cell.
5. No `R/` file modified.
6. `devtools::install()` was not run; the installed build is unchanged.

---

## OUT OF SCOPE

No `devtools::install()`. No `R/` change. No before-arm run for `cert20` if §3 fails — that is a separate
decision. No re-run of the sweep or the probe. No manuscript edits.
