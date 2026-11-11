# TASK — audit revisit: A-009, A-001, A-003, A-163 (read-only), 2026-11-11

Complete on its own. Run from the forestsearch repo root on `master`.

## Purpose

- Four findings from the 2026-11-09 code audit were passed on with no work. Larry wants to
  know, for each one, two facts:
  1. **Does the claim hold in the current source at `master`?**
  2. **Is any committed analysis exposed to it?**
- The output is one table plus the evidence behind each row. **No fixes, no proposed diffs,
  no recommendations.** Larry picks what, if anything, goes further.
- A verdict of "does not hold" or "no analysis exposed" is a complete answer. Do not stretch
  a claim to keep it alive.

## Rules

- **Read-only.** In forestsearch the only files this task adds are this task document and
  the report. No change to `R/`, `tests/`, `man/`, `vignettes/`, `inst/`, `quarto/`,
  `DESCRIPTION` or `NAMESPACE`. In every other repository: no write of any kind, no branch
  change.
- No `git fetch`, `pull`, `push` or `merge`. Every `git add` names explicit paths.
- No package install into any library. No test suite, `check()`, `document()`, `build()`,
  render, simulation or application re-run. No parallel workers.
- **Evidence is the current source.** The audit report, comments, roxygen, NEWS and earlier
  task documents are pointers only. Trace the executing code path and quote it.
- A numerical gap between two different methods (e.g. full bootstrap vs multiplier
  resampling) is not evidence. Evidence is code against its own specification.
- Gates are stop-on-failure, never stop-to-ask. A missing repository, file or line number is
  not a failure: record it under OPEN ITEMS as `CANNOT DETERMINE, because …` and continue.
- Line numbers below are the audit's, taken at base `5b994177`. They may have moved. Locate
  each site by function and content, and report the line at `master`.

## Step 0 — record

1. `git status --porcelain` is empty and the branch is `master`. **Gate.**
2. Copy `~/Downloads/TASK_audit_revisit_2026-11-11.md` to `dev/tasks/` and commit it on
   `master` (`dev: add audit revisit task`).
3. Record for the report: hostname, `master` SHA, `DESCRIPTION` version. The handoff expects
   `02f88dac`, 0.4.1.9003. If HEAD differs, record `git log --oneline 02f88dac..HEAD` and
   continue.
4. Locate the audit report `dev/reports/REPORT_code_audit_2026-11-09.md` (at `master`, else
   `git show fs-code-audit:dev/reports/REPORT_code_audit_2026-11-09.md`). Record where it
   was found. It is a pointer to sites and reproducers only. If it is not found, continue;
   this document carries the claims.
5. Locate local clones of `fs-post-selection` and `fs-glms-interpretable` (try
   `~/Documents/GitHub` and `~/GitHub`). For each, record path, checked-out branch, SHA and
   whether the tree is clean. Read what is checked out; do not change it.

## What "committed analysis" means (for every exposure question)

Tracked files only (`git grep` / `git ls-files`), in:

- forestsearch: `quarto/` (applications and simulation campaign templates), `vignettes/`;
- fs-post-selection: `manuscript_rev1/` and `manuscript/`, with their generators;
- fs-glms-interpretable: `quarto/` and the manuscript directory.

`dev/` scratch is out of scope. For each exposed call site report repo, file:line, and the
argument values as written in the call (or the default they resolve to).

## Step 1 — A-009: GLM full bootstrap and the parent's estimator settings

- **Audit's claim** (`R/bootstrap_dofuture_main.R:302`): on GLM outcomes,
  `forestsearch_bootstrap_dofuture()` rebuilds its effect estimator without the parent
  fit's `adverse_outcome` or `ps_adjust_method`.
- **Holds at master?** Trace how the parent fit's estimator settings reach a bootstrap
  replicate. Answer separately for:
  1. the replicate's own search (the re-run of `forestsearch()` on the resample);
  2. the estimates the bias correction is built from (the bootstrap-selected subgroup on
     the bootstrap sample and on the original data, and the parent's subgroup on the
     bootstrap sample).
  For each, say whether `adverse_outcome` and `ps_adjust_method` carry the parent's value
  or fall back to a default, and give the default. Use the argument names found in the
  source; if the audit's names are not the source's, say so.
- **Exposure.** Every committed analysis that calls `forestsearch_bootstrap_dofuture()` on
  a non-survival fit. Exposed = the parent fit set either argument away from its default.

## Step 2 — A-001: continuous outcome with DINA, thresholds unset

- **Audit's claim** (`R/forestsearch_main.R:2766`, `:3237`): with a continuous outcome,
  the DINA identifier and both screening thresholds left unset, the fit applies a harm
  floor of 1.25 while the returned object records 0 and bootstrap/CV replays resolve 0.
- **Holds at master?** Trace the DINA floor for the gaussian family at each site where it
  is derived, and the value recorded in the returned object and used by a replay. Give the
  floor used by the fit and by a replay in three cases: both thresholds unset; only the
  first set; only the second set. Name the two threshold arguments as the source has them.
- Say in one line whether this is the same site as the "`use_dina` unoriented floor"
  that an earlier review placed near `forestsearch_main.R:3054-3063`, or a different one.
- **Do not judge which floor is correct.** That is a method choice and is Larry's.
- **Exposure.** Every committed continuous-outcome DINA run. Exposed = both thresholds
  unset in the call. For an exposed run, name which of its reported products come from a
  replay (bootstrap, CV) and which do not.

## Step 3 — A-003: classification metric when no subgroup is found

- **Audit's claim** (`R/run_simulation_analysis.R:1025` and three further sites): under an
  alternative, a replicate that identifies no subgroup stores specificity as |Hc|/N, not 1.
- **Holds at master?** For each site, quote the lines and write the formula stored for
  every classification column in two cases: a subgroup is found; none is found. Then give
  one verdict per column:
  - `DEFECT`: the column's own found-case formula, applied to an empty identified set,
    gives a different value from what is stored;
  - `LABEL`: the stored value is correct for a different defined metric (say which, e.g.
    the positive predictive value of the complement);
  - `AS DESIGNED`: stored value and definition agree.
  Take the definition of each column from where it is reported: quote the defining sentence
  or caption from the manuscript source or the function's documentation, with file:line.
- **Exposure.** Every reported table or figure that shows the column or a summary of it.
  For each, say whether the summary runs over all replicates or only over replicates with a
  subgroup found. Exposed = summarised over all replicates under an alternative.
- **Magnitude, from stored payloads only.** For each exposed table whose committed payload
  holds per-replicate values, `readRDS()` it (read-only) and report, per cell: the detection
  rate, the reported mean, and the mean with no-find replicates set to 1. No re-run. If
  per-replicate values are not stored, record `CANNOT DETERMINE`.

## Step 4 — A-163: censoring model in the GBSG AFT generator

- **Audit's claim, marked "suspected" only** (`R/sim_aft_gbsg.R:638`): censoring-model
  coefficients are back-transformed with the event model's scale, not the censoring
  model's.
- **Holds at master?** Read the generator. Write out, as equations with the source's
  variable names, how an event time and a censoring time are generated: which fitted model
  each coefficient vector comes from, and which scale parameter is applied to which. Quote
  the lines. Verdict: `CONSISTENT` (each model's coefficients use its own scale, or the
  transformation is correct as written for another stated reason), `INCONSISTENT`, or
  `CANNOT DETERMINE`.
- `CONSISTENT` closes the row with no exposure section.
- Only if `INCONSISTENT`: list the generator entry points that reach this code. Nothing
  more.

## Optional reproducers

- Allowed only where reading leaves a verdict uncertain. Synthetic data, sequential plan,
  `Rscript --vanilla` from `/tmp/audit_revisit_2026-11-11/`, package loaded with
  `pkgload::load_all()` from the repo tree (no install). The audit report's reproducer
  appendix may be transplanted.
- Hard caps: 2 minutes per reproducer, 10 minutes in total. On a cap, stop that reproducer,
  record the wall clock, and continue with the reading verdict.
- Nothing is written into any repository.

## Step 5 — report

Write `dev/reports/REPORT_audit_revisit_2026-11-11.md`:

1. **Provenance**: the Step 0 records.
2. **Summary table**, one row per finding, exactly these columns:

   | ID | Site at master (file:line) | Claim (one line) | Holds at master | Exposed committed analyses | Would alter a reported number |
   |---|---|---|---|---|---|

   - "Holds at master": `YES`, `NO`, `PARTLY` (say which part), or `CANNOT DETERMINE`.
     For A-003 use the per-column verdicts; for A-163 use its verdict words.
   - "Exposed": a count and the names, or `none`.
   - "Would alter a reported number": `yes`, `no` or `unknown`, with the table or figure
     named when `yes`.
3. **One section per finding**: quoted code (at most 30 lines per quote), the trace, the
   exposure table, any reproducer output verbatim.
4. **OPEN ITEMS**: every `CANNOT DETERMINE`, each with its reason.

No fixes, diffs or recommendations anywhere in the report.

## Step 6 — close

1. `git add dev/reports/REPORT_audit_revisit_2026-11-11.md` and commit on `master`
   (`dev: audit revisit report`).
2. **Gates:** `git diff --name-status HEAD~2 HEAD` lists exactly two added files (the task
   document and the report); `git status --porcelain` is empty; `git status --porcelain` in
   each other repository is unchanged from Step 0.
3. Copy the report to `~/Downloads/REPORT_audit_revisit_2026-11-11.md`.
4. Print the summary table and the gate results as the closing output. Do not push.
