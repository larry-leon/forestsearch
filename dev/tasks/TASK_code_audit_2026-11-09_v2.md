# TASK — forestsearch code audit (read-only), 2026-11-09, v2

Supersedes `TASK_code_audit_2026-11-09.md`. This document is complete on its own; do not
read the earlier file for instructions.

## Goal

- Assess whether any functions in forestsearch can be improved.
- The audit **finds and ranks; it changes no code**. No file under `R/`, `tests/`, `man/`,
  `vignettes/`, `inst/`, `DESCRIPTION` or `NAMESPACE` is modified, created or deleted.
- Deliverables: one report and one findings CSV, on branch `fs-code-audit`.

## Rules for this task

- Gates are **stop-on-failure, never stop-to-ask**. Every decision is pre-set below. If a
  gate fails, write what failed into the report (or, before the report exists, print it),
  and stop.
- Verify from the current source and git only. Do not rely on descriptions, memory, earlier
  reports or any knowledge-base document for what the code does.
- Do **not** run the full test suite, `devtools::check()`, `R CMD check`, `document()`,
  `build()`, or any pkgdown/vignette build. Do not profile.
- Do not push. Larry pushes via GitHub Desktop.
- Hard caps: any single reproducer ≤ 120 s wall clock; total reproducer wall clock ≤ 20 min.
  A reproducer that hits its cap is abandoned and the finding is classed `traced` or
  `suspected` instead. Record the wall clock of each reproducer run.

## Step 0 — session, task record, pre-conditions

1. Name this session `forestsearch-code-audit-2026-11-09`.
2. On `master`: `git status --porcelain` must be empty. **Gate:** not empty → stop.
3. Copy `~/Downloads/TASK_code_audit_2026-11-09_v2.md` to `dev/tasks/` and commit it on
   `master` (message: `dev: add code-audit task v2`). If the v1 file is already in
   `dev/tasks/`, leave it untouched. If v1 is not there, do nothing about it.
4. Record for the report header: `master` HEAD SHA after step 3, `DESCRIPTION` version,
   date, counts of files in `R/` and `tests/testthat/`.
5. Record (report only, take no action): whether `fs-articles-wrappers` is merged into
   `master` (`git branch --merged master`, and `git log` evidence).
6. **Gate:** if a branch `fs-code-audit` already exists, or
   `dev/reports/REPORT_code_audit_2026-11-09.md` already exists on any local branch → stop
   and print the branch state (`git log --oneline master..fs-code-audit`).
7. Create and check out `fs-code-audit` from `master`.

## Step 1 — scope

**Tier 1 — deep read (line by line):** `forestsearch()` and everything it calls:

- candidate cut construction;
- the subgroup search;
- the consistency screen;
- final selection (all `sg_focus` settings);
- the GRF and DINA identifiers;
- all outcome types (survival, continuous, binary);
- the MR call and what it invokes.

Start by building the call graph from `forestsearch()` by reading source, and state it in
the report (function → file). That call graph, together with the bootstrap / bias-correction
and MR inference functions an analysis calls on a `forestsearch()` result, is the
**conducted-analysis path** for the `on_path` field. State the path definition you used.

**Tier 2 — static sweep:** every other function in `R/` (exported and internal), read for
the dimensions below at lower depth.

## Step 2 — dimensions

| Code | Dimension | What to look for |
|---|---|---|
| A | Correctness and robustness | wrong results; unhandled inputs (zero events, single-level factors, missing values, tiny subgroups, ties); silent coercion; arguments accepted but ignored; error paths that return partial objects; partial matching; `1:n` on empty; `sapply` type instability; scoping / `<<-` hazards |
| B | Dead code and duplication | unreachable branches; functions with no caller and no export; near-duplicate helpers |
| C | Documentation drift | roxygen `@param` / `@return` / defaults / examples that disagree with the code; undocumented arguments; stale `@seealso` |
| D | API consistency | argument naming and ordering, default conventions, return-object structure and class across sibling functions |
| E | Test gaps (static reference only) | exported functions and Tier 1 internals never referenced under `tests/testthat/` (grep; do not run tests) |
| F | CRAN hygiene | `T`/`F`; `cat`/`print` instead of `message`; `set.seed` or `options`/`par`/`plan`/`setwd` changed without `on.exit` restore; writes outside `tempdir()`; `:::`; undeclared or unused dependencies; `requireNamespace` guards for Suggests; cores > 2 in examples/tests; `installed.packages()` |
| G | Maintainability | the **ten worst only**: by length, nesting, argument count, mixed responsibilities. One row each. |

**Out of scope — do not raise as findings:**

- Efficiency or speed of any kind.
- Methodology: no alternative estimands, selection rules, thresholds or variance
  constructions.

**Closed — do not raise, do not call dead code:**

- Winner-only and winner-floor IJ variance code paths (retained research options).
- `field_uniform` (documented research option).
- The existence of the deprecated arguments `vi.grf.min`, `max_n_confounders`,
  `max.minutes` (their *inertness* is invariant 4 below and is in scope).
- The removal of `fit_causal_forest_glm()`.
- `mr_estimates_table()`'s attribute and the case study's MR print block reading stored
  limits; a `labels` argument on the KM plot; articles with untracked HTML;
  `mrct_record()` reading CSV only.

If one of these is encountered, list it once under "Closed items encountered" with
file:line and no severity.

## Step 3 — the five invariants (each gets an explicit verdict)

For each: verdict `holds` / `violated` / `not established`, the evidence class, and the
file:line trail.

1. With `use_grf = FALSE` and `use_lasso = FALSE`, no model-based pre-processing occurs
   (no GRF, lasso or other fitted model influences cuts or candidate factors).
2. `future::plan`, `options()` and the RNG state (`.Random.seed`, RNG kind) are restored on
   exit, including on error exit, by every function that changes them.
3. Every documented default (roxygen) equals the formal default.
4. Deprecated arguments are inert: supplying them changes nothing but a deprecation
   message.
5. Sequential and multisession execution give the same selection for the same seed.

Invariants 1–4 are established by trace. Invariant 5: trace the seeding / `future.seed`
handling; one small reproducer (2 workers, smallest dataset that selects a subgroup,
within the caps) is permitted. A `violated` verdict is a finding and is ranked with the
others.

## Step 4 — recording findings

One row per finding, in both the report table and the CSV:

| Field | Values |
|---|---|
| `id` | `A-01`, `B-03`, … (dimension letter + sequence); invariant violations `INV-1`… |
| `tier` | 1 / 2 |
| `dimension` | A–G or INV |
| `file`, `line`, `function` | exact; line numbers at the audited SHA |
| `severity` | S1–S4 (below) |
| `on_path` | on / off the conducted-analysis path |
| `evidence` | reproduced / traced / suspected |
| `summary` | one sentence: what is wrong |
| `alters_numbers` | yes / no / unknown — would fixing it change the numbers of an analysis already run with default or commonly used settings |
| `fix_summary` | one sentence: what a fix would be (describe; do not write the patch) |
| `fix_class` | add-only / moves code / changes behaviour / method-adjacent |

**Severity:**

- **S1** — a wrong numerical result, wrong selection, or a failure in a supported,
  documented configuration.
- **S2** — wrong or misleading behaviour in an edge or non-default configuration; a
  silent failure where an error is warranted; an invariant violated without wrong numbers.
- **S3** — a quality defect with no wrong result: documentation drift, duplication, API
  inconsistency, test gap, CRAN-hygiene item.
- **S4** — cosmetic.

**Evidence:**

- `reproduced` — a minimal script was run and showed it. Run reproducers from
  `tempdir()` with `pkgload::load_all()` on the checked-out tree; write nothing into the
  repo; paste the script, the output and the wall clock into the report appendix.
- `traced` — established by reading: give the file:line chain from input to effect.
- `suspected` — plausible, not established.

**Rules:**

- An **S1 cannot be `suspected`**. An S1 without a reproducer or a complete trace is
  recorded as S2 `suspected` with a note "candidate S1".
- Do not inflate: a pattern repeated in many places is one finding with a location list.
- Do not propose patches, and do not attach a follow-up task to any finding.

## Step 5 — the report

Write `dev/reports/REPORT_code_audit_2026-11-09.md`:

1. **Header** — SHA, version, date, file counts, `fs-articles-wrappers` merge status, total
   reproducer wall clock.
2. **Gate checklist** — every gate in this document, pass/fail.
3. **Summary counts** — table: severity × on/off path; and dimension × severity.
4. **Invariants** — table of the five verdicts.
5. **Ranked findings** — **one table**, all findings, ordered S1→S4, then on-path before
   off-path, then reproduced → traced → suspected. Columns as in Step 4.
6. **Notes per finding** — only for S1 and S2: the trace or reproducer reference, in
   bullets.
7. **Call graph and path definition** (Tier 1).
8. **Test-gap list** (dimension E) — exported functions and Tier 1 internals with no
   reference under `tests/testthat/`.
9. **Maintainability — ten worst.**
10. **Closed items encountered.**
11. **Coverage statement** — the files read at Tier 1, the files swept at Tier 2, and
    anything in `R/` not read, with the reason.
12. **Appendix** — reproducer scripts, output, wall clocks.

Format: plain markdown; tables for results; bullets for everything else; no narrative
prose sections.

Write `dev/reports/FINDINGS_code_audit_2026-11-09.csv` with the Step 4 fields as columns,
one row per finding, same order as the ranked table.

## Step 6 — final gates and close

1. Commit the two files on `fs-code-audit` (message: `dev: code audit report and findings`).
2. **Gate:** `git diff --name-status master..fs-code-audit` lists exactly two lines, both
   `A`: the report and the CSV. Anything else → revert the stray change, re-check, and
   record the incident in the report's gate checklist.
3. **Gate:** `git status --porcelain` is empty (no untracked reproducer output, no
   `.Rhistory`, no compiled artefacts).
4. **Gate:** row count of the CSV equals the row count of the ranked table; no S1 row has
   `evidence = suspected`.
5. Do not merge, do not push, do not delete the branch.
6. Copy to `~/Downloads`:
   - `dev/reports/REPORT_code_audit_2026-11-09.md`
   - `dev/reports/FINDINGS_code_audit_2026-11-09.csv`
7. Print: the SHA audited, the commit SHA on `fs-code-audit`, the summary counts table, and
   the gate checklist.
