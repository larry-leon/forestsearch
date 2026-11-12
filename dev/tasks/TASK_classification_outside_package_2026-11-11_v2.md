# TASK — classification metrics: remaining instances outside `R/`, and the `paper_simulations` re-render, v2, 2026-11-11

Complete on its own. Run from the forestsearch repo root. Authorized by Larry 2026-11-11:
finish all classification instances on the same branch, and re-render the
`paper_simulations` article with the fixes applied. **Larry is away: this task runs
unattended and never stops to ask.**

Supersedes `TASK_classification_outside_package_2026-11-11.md` (v1; do not run it) and
`TASK_paper_sims_scope_2026-11-11.md`.

What v2 changes from v1:

- Step 0 runs the empty-set task first if it has not been run.
- Part R re-renders `paper_simulations` and accounts for what changed.
- The parts are independent: a failed gate stops its own part only.

## Purpose

- The instances are the rows of `dev/reports/REPORT_npv_nofind_instances_2026-11-11.md`,
  §2.3 and §3 (Q1 to Q12 in forestsearch, G1 to G6 in fs-glms-interpretable). Row labels
  below refer to that report.
- The definition (León et al. 2024, Section 3): with no subgroup identified, the estimated
  complement is everyone. Under an alternative the replicate has `sens` = 0, `ppv` = 0,
  `spec` = 1, `npv` = |Hc|/N. Under the null: `spec` = 1, `npv` = 1, `sens` and `ppv` NA.
- Parts:
  - R: re-render `vignettes/articles/paper_simulations.qmd` with the fixed package, and
    compare its numbers with the last recorded render;
  - A: correct document sources in forestsearch `quarto/`;
  - B: correct the stored summary in the tracked payloads (no simulation);
  - C: report what each other affected render needs (no render);
  - D: correct the payload scripts, their CSVs and checks in fs-glms-interpretable.
- **No change under `R/`, `tests/`, `man/` or `vignettes/` in forestsearch. The only
  render is Part R. Nothing under `manuscript/` in fs-glms-interpretable.**

## Rules

- No `git fetch`, `pull`, `push` or `merge`. Every `git add` names explicit paths.
- No install into the user or site library. The only install is Part R's temporary
  library under `/tmp`.
- fs-post-selection: not touched.
- Evidence is the current source. Quote file:line.
- **Gates are stop-on-failure for their own part, never stop-to-ask.** If a gate fails,
  write what failed into the report, leave that part's uncommitted edits out of every
  commit (restore them), and go on to the next part. Step 0 is the exception: a failed
  Step 0 gate stops the task. Part E (the report) always runs.
- A sub-item marked "list as BLOCKED" is not a failure.
- No timing estimates. Report measured wall clocks.

## Step 0 — forestsearch branch and record

1. `git status --porcelain` is empty. **Gate.**
2. Check out `fix/npv-nofind-instances`. **Gate:** it exists.
3. If `dev/reports/REPORT_empty_set_nofind_2026-11-11.md` is not on the branch, the
   empty-set task has not been run. Run it now, completely, with all its own gates:
   `dev/tasks/TASK_empty_set_nofind_2026-11-11.md` if present, else
   `~/Downloads/TASK_empty_set_nofind_2026-11-11.md`. **Gate:** it ends with all its gates
   passed and its report committed. If neither file exists, or it ends on a failed gate,
   stop this task.
4. Copy `~/Downloads/TASK_classification_outside_package_2026-11-11_v2.md` to `dev/tasks/`
   and commit it on the branch (`dev: add classification outside-package task v2`).
5. Record the HEAD SHA before step 4 (the forestsearch start SHA) and the `DESCRIPTION`
   version.

## Part R — re-render `paper_simulations` with the fixed package

Start this part first and let the render run in the background while Parts A to D
proceed. It does not depend on them.

1. **Precedent.** Find how the article was last rendered and recorded: the task document
   and report of the 2026-10-30 articles work (`dev/tasks/`, `dev/reports/`,
   `dev/reports/articles_finish_2026-10-30/`), including
   `results/paper_simulations_output.txt`. Record the method, settings, seeds, the package
   version used then, and any recorded wall clock.
2. **Temporary library.** Export the branch (`git archive HEAD`) to
   `/tmp/fs_src_2026-11-11/` and install it into `/tmp/fs_lib_2026-11-11/`
   (`R CMD INSTALL --library=`). **Gate:** with that library first on the library path,
   `packageVersion("forestsearch")` is the branch's `DESCRIPTION` version and
   `find.package("forestsearch")` is under `/tmp/fs_lib_2026-11-11/`.
3. **Scratch copy.** Copy the article and what it needs to
   `/tmp/paper_sims_rerun_2026-11-11/` and render there. The repository copy of the
   article is not modified and nothing is written under `vignettes/`.
4. **Render** the same way as the precedent (same settings and seeds), with the temporary
   library first on the library path for the R session and its workers. Log to the
   scratch directory. Record start, end and measured wall clock.
   - **Hard cap:** three times the recorded wall clock of the precedent render, and at
     least 1 hour; if no wall clock is on record, 6 hours. On the cap, stop the render
     and its workers, record the wall clock, and do not retry.
   - **Gate:** the render finishes without error and its output holds the classification
     tables.
5. **Outputs.** Put the rendered output and the same kinds of result files the precedent
   kept under `dev/reports/paper_sims_rerun_2026-11-11/`.
6. **Accounting.** Compare every reported number in the new output with the 2026-10-30
   record. State first whether settings and seeds are identical between the two renders.
   One table per block of the article:

   | Block | Metric and cell | 2026-10-30 | 2026-11-11 | Difference |
   |---|---|---|---|---|

   Group the rows:
   - (i) specificity and NPV of the complement (`sens(Hc)`, `ppv(Hc)`, Spec, NPV): these
     are expected to change;
   - (ii) detection rate and PPV: these can change only through the empty-set rule;
   - (iii) everything else. For each number in (iii) that changed, name the forestsearch
     commits between the two package versions that touch the functions the block calls.
     List them; do not explain or judge.

   Then one line per group: how many numbers, how many changed, the largest change.
7. The pkgdown site is not rebuilt.

## Part A — forestsearch document sources (`quarto/`, `.qmd` and `.R` only)

Edit sources only. Do not edit `.html` or `.rmarkdown` files, and nothing under
`dev/reports/` except this task's own outputs.

- **A1 (Q1).** `actg175_binary_methods_template.qmd`, `oc_agg_from_results()`: a no-find
  row stored with NPV missing is filled with the row's `propHc_true` when present, and
  with 1 only when it is absent. The spec fill stays 1.
- **A2 (Q3).** Both `reproduce_table1_m1.qmd` files: the same change at the NPV fill.
- **A3 (Q2, Q8).** Prose that states NPV = 1 for a no-find replicate under H1, or that the
  values match `fs.estimates.out()`: replace with one sentence in the document's own
  notation, saying that when no subgroup is detected the estimated complement is everyone,
  so the replicate contributes sensitivity 0 and PPV 0 for the subgroup, and sensitivity 1
  and PPV |Qc|/N for the complement (León et al. 2024, Section 3). Apply the same
  replacement to every `.qmd` copy. Report the file count and one before/after example per
  cluster of identical text.
- **A4 (Q9).** The four documents whose own extractors return missing values for a no-find
  replicate and then average with `na.rm` (`grf_dmin_calibration.qmd`,
  `validation_glm_simulation_study.qmd`, `archive/actg175_simulations_binary_working.qmd`,
  `dina_subgroup_signature_simulation.qmd`): at each no-find return, store the definition's
  values from the truth vector in scope.
- **A5 (Q5).** `actg175_binary_m1_harm_fs1-fs6_summary.qmd`: it averages bundles produced
  before the A-003 fix. After the bundles are loaded and before any average, re-score the
  no-find rows of an alternative: `spec` = 1 and `npv` = the row's `propHc_true`. Add a
  two-line comment saying why.
- **Not edited, by decision:** Q7 (records of a past render), Q10, Q11 (detected-only
  values, correct for what they are) and Q12.

**Gate:** the files changed by Part A are only `.qmd` and `.R` files under `quarto/` that
belong to rows Q1, Q2, Q3, Q5, Q8 or Q9. Commit them by explicit path
(`fix: classification no-find values in simulation document sources`).

## Part B — stored summary in the tracked payloads (no simulation)

The report (§2.4) found 106 tracked `.rds` payloads that store `$oc$identification` from
`fs_mr_oc_summary()`; in 62 the unconditional NPV holds the old value.

1. Write the script beside the report as
   `dev/reports/REPORT_classification_outside_package_2026-11-11_resummarise.R`. Load the
   package from the forestsearch tree with `pkgload::load_all()`; sequential.
2. For each of the 106 payloads: recompute `fs_mr_oc_summary()` on the payload's own
   stored results with the branch's code, in the way the payload's batch document calls
   it.
3. **Per-file gate:** the recomputed `$oc` is identical to the stored one except the
   unconditional NPV in `$oc$identification`. If anything else differs, do not write that
   file; list it as BLOCKED with the difference.
4. Where only that value differs, replace it in the payload and save the file in place,
   with the same serialization version and compression as the original. Every other
   element of the object must be identical before and after (`identical()`); check it.
   Do not rewrite a file whose value does not change.
5. Write `dev/reports/REPORT_classification_outside_package_2026-11-11_payloads.csv`:
   file, old value, new value.
6. Hard cap: 15 minutes for the whole part. On the cap, stop the part, keep what is
   verified, and list the rest as not done.

**Gate:** the set of changed `.rds` files equals the files listed in the CSV; the report
states the count (expected 62) and explains any difference from it. Commit the payloads,
the script and the CSV by explicit path
(`fix: stored unconditional NPV in MR payloads re-summarised`).

## Part C — what each other affected render needs (read-only; no render)

One table, one row per tracked render that displays classification values produced by the
`run_simulation_analysis()` extractors or by a Part A source: the 17 renders of row Q6,
`actg175_binary_methods_template.html`, and the tracked renders of the Q9 documents.

| Render | Source `.qmd` | What it reads (path; tracked or not; present on this machine) | Does a render run simulations (cells × replicates) or only read stored results? | Do the stored rows hold pre-fix values, so a plain re-render would repeat the old numbers? | Recorded wall clock of a past render (source), or none |
|---|---|---|---|---|---|

## Part D — fs-glms-interpretable (`~/Documents/GitHub/fs-glms-interpretable`)

1. `git status --porcelain` is empty there (**gate**). Record branch and SHA (the
   fs-glms start SHA). **Gate:** branch `fix/npv-nofind-instances` does not exist; create
   it from the current branch.
2. **D1 (G1, G2, G3).** In `dev/payloads/R/identifiers_metrics.R`,
   `identifiers_summary.R` and `idsweep_metrics.R`: a non-declaring replicate contributes
   NPV = 1 − `n_true` / `n_sample` for that replicate, from the quantities each script
   already reads. The fills for `sens`, `spec` and `ppv` are unchanged. If `n_true` is
   missing for a non-declaring replicate, the unconditional NPV of that cell is NA; list
   such cells. Correct the comments in those scripts that state the old value.
3. **D2.** Re-run the three scripts the way the repository normally runs them, to
   regenerate `identifiers_metrics.csv`, `identifiers_summary.csv` and
   `idsweep_metrics.csv` from the stored bundles. If they call forestsearch, they use the
   installed package as they normally do; record its version. Hard cap: 30 minutes in
   total; on the cap, stop this part.
   - **Gate:** in each CSV, only the unconditional NPV columns (`npv_uncond`, and its
     Monte Carlo standard error where present) change. Compare old and new by key; every
     other value is identical.
4. **D3 (G5).** Update the checks that pin the old value
   (`dev/writeup/identifiers/R/reconcile_identifiers_summary.R:91-96`,
   `verify_oc_section.R:772` and `:488-489`, `:1015-1018`, `conv.R`) so that they check
   the new value. Where a check cannot verify the new value from what it reads, remove
   only that assertion and leave a one-line comment saying so. Run each check that does
   not need a render (hard cap 10 minutes each). **Gate:** every check that was run
   passes. List any check not run, with the reason.
5. **D4 (G6).** Correct the statements of the old convention in
   `dev/payloads/identifiers/README.md`, `dev/writeup/identifiers/identifier_evaluation_writeup.qmd`
   and the draft sections under `dev/writeup/identifiers/draft/`. Do not edit decision
   records, requests, reports or status updates: they are history.
6. **D5 — manuscript impact, for Larry's decision.** Nothing under `manuscript/` is
   edited or rendered. Report:
   - for supplement Table `tbl-supp-classification`, every cell of the NPV column: old
     value, new value, difference (from `identifiers_metrics.csv`);
   - for the other two CSVs: the minimum, mean and maximum change in unconditional NPV;
   - every file:line under `manuscript/` that displays or states the old convention.
7. Commit on the fs-glms branch by explicit path
   (`fix: unconditional NPV for a non-declaring replicate is |Hc|/N`).
   **Gate:** `git diff --name-only <fs-glms start>..HEAD` lists nothing under
   `manuscript/`; `git status --porcelain` is empty. If this part failed a gate, restore
   the working tree, leave the branch without the failed edits, and say so.

## Part E — report and close (forestsearch; always runs)

1. Wait for the Part R render to finish or reach its cap. Then write
   `dev/reports/REPORT_classification_outside_package_2026-11-11.md`, in bullets and
   tables:
   - provenance for both repositories, and whether Step 0 ran the empty-set task;
   - **a status table**: one row per part (R, A, B, C, D), `passed / failed gate (which)`;
   - **a summary table for Larry**: one row per instance (Q1 to Q12, G1 to G6), with
     columns `Instance | Done / BLOCKED / not an instance / needs a render or re-run | What changed | What is still needed`;
   - Part R: precedent, library gate, measured wall clock, the accounting tables and the
     three group lines;
   - Part A: files changed per row, with the before/after examples;
   - Part B: counts, and the BLOCKED list if any;
   - the Part C table;
   - Part D: the D5 tables, the checks run and not run, and the fs-glms branch HEAD SHA;
   - OPEN ITEMS.
2. Commit the report and the Part R outputs on the forestsearch branch by explicit path
   (`dev: classification outside-package report and paper_simulations re-render`).
3. **Gates:** `git diff <start> --stat -- R/ tests/ man/ vignettes/` shows only what the
   empty-set task changed, if Step 0 ran it, and otherwise nothing;
   `git status --porcelain` is empty in forestsearch and in fs-glms-interpretable;
   fs-post-selection is unchanged.
4. Copy the report to `~/Downloads/REPORT_classification_outside_package_2026-11-11.md`.
5. Print the status table, the Part R group lines, the gate results and both branch HEAD
   SHAs. Do not push or merge in either repository.
