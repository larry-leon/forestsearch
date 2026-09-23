# TASK — `R CMD check --as-cran` cleanup

**Repo:** `larry-leon/forestsearch`, branch `feature/glm-extension` — run CC from the forestsearch clone.

**Run only after** `TASK_pconsistency_digits_argument_2026-09-23.md` has committed in full. Do not start while
any check from that task is still running.

**Purpose.** The `--as-cran` baseline captured on 2026-09-23 showed 2 ERRORs, 1 WARNING and 3 NOTEs. The ERRORs
block submission outright. This task diagnoses each item and fixes the ones that can be fixed without touching
the logic of an exported function.

**`R/` CALLOUT.** Expected to **move existing code and correct documentation, encoding, examples and tests**.
It must **not change the behaviour or method of any exported function**. §3's boundary rule enforces that: if an
item's root cause is in function logic, it is recorded and left alone.

**Time.** Two `R CMD check` runs at roughly 10 minutes each, plus diagnosis and fixes. **Hard abort at 2 h.**
Any single item estimated at more than 30 minutes is skipped, recorded with its estimate, and the rest continue
— do not stop the task for it.

---

## 0. First action

1. Copy this file to `dev/tasks/TASK_cran_check_cleanup_2026-09-23.md`; `git add` that path; commit
   `docs(tasks): add R CMD check cleanup task (2026-09-23)`.
2. Record HEAD, installed forestsearch version and build date, R version, platform.
3. `git status --short`: record pre-existing untracked files; never stage them.

No `git fetch` / `pull` / `push`. Explicit paths on every `git add`.

---

## 1. Baseline

The baseline for this task is HEAD **after** the digits task committed — not the pre-digits baseline.

If that task's post-change `R CMD check` output was saved, reuse it and say so. Otherwise run
`rcmdcheck::rcmdcheck(args = "--as-cran")` once on a clean export of current HEAD and save it. Record the full
item list verbatim.

---

## 2. Diagnosis — report the verbatim text for every item

For each of the six known items, quote the actual check output and state the root cause from source:

1. **ERROR** — the `fs_dgm_feasibility` example.
2. **ERROR** — tests: the `test-fs-dgm-feasibility.R:103` file read, and the `effect_measure` guard at
   `test-fs-family-report.R:115` (`effect_measure` appears both classified and on the out-of-scope list).
3. **WARNING** — non-ASCII characters in `R/fs_bias_coverage.R`. Report **which** characters, at which lines,
   and whether each sits in a comment, in roxygen, or inside a string literal that reaches user-visible output.
   Also report whether `DESCRIPTION` already declares `Encoding: UTF-8`.
4. **NOTE** — CRAN incoming feasibility. Quote it; several of its sub-items (new submission, maintainer field,
   days since last update) are not fixable in the package. Say which are and which are not.
5. **NOTE** — missing global-variable declarations in `fs_plot_bias_coverage`. List the variables.
6. **NOTE** — HTML tidy. **Out of scope**: this reflects `tidy` not being installed on this machine, not a
   package defect. Confirm that reading and do not chase it.

Also diagnose the `devtools::document()` inline-code error at `R/fs_bias_coverage.R:17`, which is in the same
file as item 3 and may share a cause.

---

## 3. Fixes — the boundary rule

**Permitted:** encoding, roxygen, examples, tests, `NAMESPACE`, `DESCRIPTION`, and global-variable
declarations.

**Not permitted:** any change to the logic of an exported function. If an item's root cause is there, record
the cause and the proposed fix, and move on. That is a finding, not a failure.

Guidance per item, subject to what §2 finds:

- **Non-ASCII.** Replace with ASCII where the character is in a comment or roxygen. Where it is in a string
  literal that reaches output, use a `\uXXXX` escape rather than changing the rendered text. Report every
  substitution. If the rendered text would change for a user, stop on that instance and report it.
- **Global variables.** Follow whatever convention the package already uses — check for an existing
  `utils::globalVariables()` call or `.data$` pronoun usage before introducing either. Do not impose a new
  pattern.
- **`effect_measure` guard.** Determine from source whether `fs_family_report()` actually classifies it. Fix
  whichever list is wrong — remove it from the out-of-scope list if it is classified, or classify it if it is
  not. One line either way.
- **Example and test ERRORs.** Fix if the cause is in the example or the test. If the cause is in function
  logic, the boundary rule applies.

---

## 4. Verification

- Re-run `R CMD check --as-cran` on the modified tree.
- **Report a before/after table, item by item**, against the §1 baseline. Every item is resolved, unchanged, or
  newly introduced.
- `devtools::test()`: report the suite line before and after.
- **Any new ERROR, WARNING or NOTE is a failure — stop and report.** Items that remain unresolved because of
  the boundary rule or the 30-minute cap are not failures; they are recorded.

---

## 5. NEWS and commits

- NEWS.md under the current development version, only for user-visible fixes. Encoding and documentation
  corrections are worth one line; test-only changes are not.
- Commits, explicit paths, in order: task doc; encoding and roxygen fixes; global-variable declarations; example
  and test fixes; regenerated `man/`; NEWS.md; the verification record beside the existing `REPORT_*` files.

---

## POST-CONDITIONS (machine-checkable)

1. Every one of the six items has verbatim check text and a stated root cause.
2. No new ERROR, WARNING or NOTE relative to the §1 baseline.
3. The before/after table lists every item with its disposition: resolved, unchanged (with reason), or skipped
   under the 30-minute cap (with estimate).
4. `devtools::test()` suite line reported before and after; no test passing before now fails.
5. No exported function's logic changed. Any item whose cause lies there is recorded, not fixed.
6. Every non-ASCII substitution is listed with its line and its character.
7. `man/` regenerated wherever roxygen changed.
8. Pre-existing untracked files never staged.

---

## OUT OF SCOPE

The HTML tidy NOTE (environment, not package). Any change to `pconsistency.digits`, the rounding design, or
`fs_declaration_calibration.R`. Any simulation. Any new functionality. Sub-items of the CRAN incoming NOTE that
are not fixable in the package — record them and move on.
