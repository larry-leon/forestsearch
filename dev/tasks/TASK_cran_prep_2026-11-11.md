# TASK — CRAN preparation: confirmed codebase fixes, documentation, and a check of the remaining audit findings, 2026-11-11

Complete on its own. Run from the forestsearch repo root. Authorized by Larry 2026-11-11:
the priority is a robust codebase for a CRAN submission. This task runs unattended and
never stops to ask.

## Purpose

- **Part 1** fixes defects already confirmed at `master` by earlier reports. Each is a
  **change of behaviour, not of method**, except F7 (documentation only).
- **Part 2** corrects one article source and re-renders that article.
- **Part 3** is read-only: it checks the remaining audit findings at `master` and proposes
  the smallest fix for each. It changes nothing.
- Out of scope: every simulation document re-run, both manuscript repositories, and the
  two items held for Larry's decision (the GBSG generator's censoring scale, and the fixed
  seed inside `subgroup.consistency()`).

Evidence for Part 1 (pointers; verify each from the current source before editing):

- `dev/reports/REPORT_audit_revisit_2026-11-11.md` (A-009, A-001);
- `dev/reports/REPORT_empty_set_nofind_2026-11-11.md` (open items 2 and 3);
- `dev/reports/REPORT_npv_nofind_instances_2026-11-11.md` (open item 6);
- `dev/reports/REPORT_classification_outside_package_2026-11-11.md` (open item 6);
- `dev/reports/REPORT_code_audit_2026-11-09.md` (A-004, A-059, INV-2, INV-3, C-001).

## Rules

- No `git fetch`, `pull`, `push` or `merge`. Every `git add` names explicit paths.
- No install into the user or site library. The only install is Part 2's temporary
  library under `/tmp`.
- fs-post-selection and fs-glms-interpretable: not touched.
- Evidence is the current source. Quote file:line. If a defect is not as described, make
  no change for that item and say so.
- **Each fix is its own commit.** A fix whose own gate fails is reverted completely and
  reported as FAILED; the task goes on to the next item. Never stop to ask.
- A sub-item marked "list as BLOCKED" is not a failure.
- Smallest edit that satisfies the rule. No refactor, no renaming, no other fix.
- Tests: hand-built inputs where possible. A test that has to run `forestsearch()` uses
  small synthetic data and a sequential plan, and carries `skip_on_cran()`. Hard cap:
  each test file runs in under 2 minutes.
- No timing estimates. Report measured wall clocks.

## Step 0 — record and branch

1. On `master`, `git status --porcelain` is empty. **Gate.**
2. **Gate:** `DESCRIPTION` version is 0.4.1.9006 (the classification fixes are merged).
3. Copy `~/Downloads/TASK_cran_prep_2026-11-11.md` to `dev/tasks/` and commit it on
   `master` (`dev: add CRAN preparation task`).
4. **Gate:** branch `cran-prep` does not exist. Create it from `master`. Record the start
   SHA.

## Part 1 — fixes

### F1 — GLM bootstrap uses the parent fit's estimator settings (A-009)

- Site: `R/bootstrap_dofuture_main.R`, the `make_effect_estimator()` call that builds the
  bootstrap's estimator (near `:302-310`). It omits `adverse_outcome` and
  `ps_adjust_method`, so the estimates behind the bias correction use the defaults.
- **Rule.** The bootstrap's estimator is built with the parent fit's resolved
  `adverse_outcome`, taken from the parent's stored arguments.
- The parent's own rebuild of its estimator under active propensity adjustment
  (`R/forestsearch_main.R` near `:2703-2712`) omits `adverse_outcome` too: pass it there.
- **Propensity adjustment.** Pass the parent's resolved `ps_adjust_method` only if the
  data the bootstrap estimates on carry what the weighted estimator needs. If they do not,
  leave this half unchanged and list it as BLOCKED with exactly what is missing.
- **Invariant to test:** for a parent fit with `adverse_outcome` FALSE and for one with
  TRUE, the bootstrap's estimate of the parent's subgroup on the original data equals the
  parent fit's own estimate for that subgroup, sign included. With TRUE, every bootstrap
  output is unchanged from the start SHA.
- **Consumers.** Tracked `quarto/` documents that re-orient the bootstrap output by
  `adverse_outcome` (the `.orient` lines in
  `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_fb_mr_*.qmd`): remove the
  re-orientation, so the recorded values stay in the parent's orientation. Source edit
  only, no re-run. List the files. List, without editing, anything else that assumes the
  old orientation (e.g. captions that say "RAW scale").

### F2 — continuous DINA floor comes from the resolved threshold (A-001)

- Sites: the gaussian-family harm floor in `.forestsearch_dina_select()`
  (`R/forestsearch_helpers.R` near `:1525`, called from `R/forestsearch_main.R:2766`) and
  `m_diff_sel` (`R/forestsearch_main.R` near `:3237`). With both thresholds unset the fit
  uses the raw 1.25 while the returned object and the replays use the resolved 0.
- **Rule.** At both sites the gaussian floor is the resolved screening threshold that the
  returned object records. The package documentation already says so, and the code's own
  comment calls 1.25 meaningless for a mean difference.
- With either threshold supplied, and for every non-gaussian family, the floor is
  unchanged.
- **Test:** both unset → the fit's floor, the recorded threshold and a replay's floor are
  equal; either threshold supplied → the floor equals the supplied value at both sites.

### F3 — survival GRF covariate matrix (A-059)

- Site: `R/grf_main.R` near `:226`, `apply(data[, confounders.name], 2, as.numeric)`.
  When a confounder is a factor or character column, labels become NA and numeric columns
  are formatted to 7 significant digits, while the selected cut is evaluated on the raw
  data (`R/grf_subgroup_labels.R` near `:225`).
- **Rule.** Build the matrix so that numeric columns keep their exact values and
  non-numeric columns are coded the way the GLM path codes them
  (`R/grf_subg_harm_glm.R`, `.build_grf_X`). The matrix GRF sees and the data the cut is
  evaluated on must agree.
- **Gate (test):** for all-numeric confounders the matrix is `identical()` to the one the
  start-SHA code builds. With a factor confounder present, the numeric columns equal the
  raw values exactly.

### F4 — identified set is everyone: NPV

- Sites: the found-case `npv` in `R/run_simulation_analysis.R` (both extractors) and the
  fill for a detected replicate with no stored NPV in `R/fs_mr_oc_summary.R`.
- When the identified set is everyone, the estimated complement is empty and NPV is 0/0.
  The extractors store NA; `fs_mr_oc_summary()` fills 1.
- **Rule F.** When the estimated complement is empty, NPV is 0, by the same 0/0
  convention that gives PPV = 0 for an empty identified set. In `fs_mr_oc_summary()` apply
  it only where the replicate's identified set is everyone (its size equals the sample
  size); any other missing value keeps its present treatment.
- Nothing else changes. Stored payloads are not regenerated.
- **Test:** both extractors and `fs_mr_oc_summary()` give NPV = 0 for such a replicate;
  a replicate with a non-empty complement is unchanged.

### F5 — `generate_gbsg_bootstrap_general()` (A-004)

- Site: `R/synthetic_data_perturbation.R` near `:366`. `data("gbsg")` loads nothing under
  current `survival`, so the function fails on every call. Confirm by calling it.
- **Rule.** Take the data from `survival::gbsg`. **Test:** the function returns without
  error on a small call.

### F6 — the caller's `future` plan is restored (INV-2, plan half)

- Sites: wherever package code sets a `future` plan (e.g. `R/subgroup_search.R` near
  `:174`). The caller's plan is replaced and not put back.
- **Rule.** Each exported function that sets a plan restores the caller's plan on exit
  (`on.exit(..., add = TRUE)`). A function whose documented purpose is to set the plan is
  left alone; list it.
- The fixed seed in `subgroup.consistency()` is **not** changed here.
- **Test:** the plan after a call is the plan before it. Results of the call are unchanged.

### F7 — documentation matches the code (INV-3, C-001; no behaviour change)

- For every exported function, where the roxygen states a default that differs from the
  formal's default, correct the text to the formal. Start from the audit's INV-3 list
  (12 sites, plus `field_complement`, `ci_method`, `batch_size`) and sweep for others.
- Correct the documented values of `est.scale` to the values the code accepts (C-001).
- `NEWS.md`: the development-section heading states the current development version.
- Run `devtools::document()` once. **Gate:** no file under `R/` changes outside roxygen
  comment lines in this commit; list every line changed.

### Common to Part 1

- `DESCRIPTION`: 0.4.1.9006 → 0.4.1.9007, once, in the last fix commit.
- `NEWS.md`: one entry per fix F1 to F6, each saying what changes and for whom.
- For each fix, the report gives: the site at the start SHA, the diff, the test, and every
  output whose value changes (name, old, new) on a small example.

## Part 2 — `paper_simulations` article

1. `vignettes/articles/paper_simulations.qmd`, key-metrics chunk (near `:1340-1355`): the
   loop variable has the same name as a column, so the FS and GRF lines print the same
   pooled values. Rename the loop variable so each line is computed for its own analysis.
   Change nothing else in the article. Commit.
2. Re-render the article exactly as
   `dev/reports/REPORT_classification_outside_package_2026-11-11.md` §4.2 did, with the
   branch (after Part 1) installed into a fresh temporary library under `/tmp`, in a
   scratch copy. Hard cap 1 hour. Outputs to `dev/reports/paper_sims_rerun_2026-11-11b/`.
3. Accounting against `dev/reports/paper_sims_rerun_2026-11-11/`: every number that
   changed, by block, with the fix (F1 to F6 or the chunk) that accounts for it. Numbers
   that are run-dependent (dates, timings, memory) are counted, not listed.
4. The pkgdown site is not rebuilt.

## Part 3 — remaining audit findings (read-only)

For each item, from the current source at the branch head:

| ID | Site at master (file:line) | Holds? (YES / NO / PARTLY) | Reachable from an exported function with ordinary inputs? (how; default path or not) | Smallest fix, described | Kind (behaviour / method / documentation) | Would a default-path result change? |
|---|---|---|---|---|---|---|

Items (claims as the audit states them; they are pointers, not evidence):

- A-024: DINA ignores `offset.name` for count outcomes (`forestsearch_helpers.R`).
- A-042: for GLM outcomes LASSO returns the intercept as "selected", so when it keeps
  nothing every confounder is dropped (`get_FSdata_helpers.R`).
- A-047: GRF/DINA cuts dropped as duplicates when one name is a substring of another
  (`get_fsdata.R`).
- A-055: the RD standard error is near 0 or very large when an arm is at 0% or 100%, yet
  the fit reports converged (`glm_effect_estimators.R`).
- A-062: `dmin.grf` acts as an all-or-nothing gate (`grf_main.R`).
- A-064: count-outcome GRF is not oriented (`grf_subg_harm_glm.R`).
- A-103: the CV "Exact" metric counts a two-factor superset as exact
  (`forestsearch_cross_validation.R`).
- A-107: `fpr_calibration()` assumes a binary outcome when none is given.
- A-125: the MRCT simulation's default confounder set includes internal generator columns
  and the region indicator (`mrct_simulation.R`).
- A-147: summary bins count a value exactly on a boundary one bin low
  (`summarize_bootstrap_subgroups.R`).
- A-150, A-162: categorical and ordinal columns not perturbed
  (`synthetic_data_perturbation.R`).
- A-151: continuous covariates left out of the HTE test (`test_hte_crump.R`).
- A-156: `k_inter` calibrated on 5,000 rows while the generator is built on 25,000
  (`fs_fdr_report.R`).
- INV-5: sequential and multisession results differ under `split`, `IRD` or the
  closed-form fallback (`subgroup_consistency_main.R`).
- The `dina_subgroup()` call on the `use_dina` route omits `tau_sign`
  (`forestsearch_main.R` near `:3237-3251`).
- A non-empty FS set with NA members stops with "missing value where TRUE/FALSE needed"
  (`run_simulation_analysis.R` near `:1133`); CV fold assignments {1, NA} or {0, NA} stop
  with "subscript out of bounds" (`forestsearch_cross_validation.R` near `:1412`).

Optional reproducers as in earlier tasks: synthetic data, sequential, `pkgload::load_all()`,
2 minutes each, 15 minutes in total. No change to any file for Part 3.

## Part 4 — tests and check

1. After each fix, run that fix's own test files. **Gate (per fix):** 0 failures.
2. After Parts 1 and 2, once: the full suite and `devtools::check()`, on the check route
   of `dev/reports/REPORT_empty_set_nofind_2026-11-11.md` §6.
   - Baseline, from that report (do not re-run it): suite 613 blocks, 81 files, 8,270
     expectations, 0 failures, 3 skips, 38 warnings; check 0 errors / 0 warnings /
     0 notes.
   - Hard caps: suite 30 minutes, `check()` 60 minutes.
   - Pre-set decision: a failing expectation that pins a value one of F1 to F6 changes on
     purpose is updated and listed (file:line, old, new, which fix).
   - Any other failure, or a new check ERROR, WARNING or NOTE: find the fix commit
     responsible, revert it, record it as FAILED, and run the suite and `check()` once
     more. If they still fail, stop and report.

## Part 5 — report and close

1. Write `dev/reports/REPORT_cran_prep_2026-11-11.md`, in bullets and tables:
   - provenance and the commit list;
   - **a status table for Larry**: one row per item F1 to F7, the article chunk and the
     render, with columns `Item | Fixed / BLOCKED / FAILED / not as described | What changes for a user | Outputs changed on a small example`;
   - per fix: site, diff, test, changed outputs;
   - Part 2: wall clock and the accounting;
   - the Part 3 table;
   - suite and check results beside the baseline, with measured wall clocks;
   - OPEN ITEMS.
2. Commit the report and the Part 2 outputs on `cran-prep` by explicit path.
3. **Gates:** `git status --porcelain` is empty; the two sibling repositories are
   unchanged; `git diff <start> --stat` lists only files the report names.
4. Copy the report to `~/Downloads/REPORT_cran_prep_2026-11-11.md`.
5. Print the status table, the gate results and the branch HEAD SHA. Do not push or merge.
