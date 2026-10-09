# TASK — CRAN release, Stage 2: remote checks (forestsearch 0.4.0)

- Date: 2026-10-09
- Repo: `larry-leon/forestsearch`, branch `feature/glm-extension`, at or after the Stage 1 report commit.
- Authorization (Larry's kickoff of this file): upload the package to **win-builder R-devel, win-builder R-release and mac-builder** — one upload each. Nothing else leaves the machine.
- Not authorized: CRAN submission (`devtools::submit_cran()`, `devtools::release()`, the web form), R-hub, `git fetch`/`pull`/`push`, tags, installs into the user library, any edit to `R/`, `man/`, `DESCRIPTION`, `NAMESPACE`, `tests/`, `vignettes/`.
- R/ callout: none.
- Gates are stop-on-failure.

## First action

1. Copy this file to `dev/tasks/TASK_cran_release_remote_checks_2026-10-09.md` and commit it alone.

## Preconditions (stop on failure)

- P1. `git status --porcelain` shows no tracked-file modifications.
- P2. `DESCRIPTION` reads `Version: 0.4.0`; HEAD contains `dev/reports/REPORT_cran_release_remediation_2026-10-09.md`. Record HEAD SHA.
- P3. `pandoc`, `quarto`, `pdflatex` on `PATH`; `devtools` loadable (`requireNamespace("devtools")`). Do not install anything.

## Steps

1. Build once in a scratch directory (outside the repo): `R CMD build .` with vignettes (≤ 45 min). Record tarball name and size.
2. Upload that built tarball (not a rebuild) — each call ≤ 10 min:
   - `devtools::check_win_devel(pkg = <tarball path>? )` — if the devtools version requires a package directory rather than a tarball, call it from the repo root with `pkg = "."` and record that it rebuilt.
   - `devtools::check_win_release(...)` (same rule).
   - `devtools::check_mac_release(...)` (same rule).
   - Answer any interactive confirmation only with "yes, upload"; if a call prompts for anything else (e.g. email change, submission), abort that call and record it.
   - Record each call's console output verbatim, including the mac-builder results URL.
3. Write `dev/reports/REPORT_cran_release_remote_checks_2026-10-09.md` (bullets): HEAD SHA, tarball size, the three upload outcomes verbatim, the mac-builder URL, and the note that win-builder results arrive by email to the maintainer address `larry.leon.05@post.harvard.edu` (~30–90 min).
4. Commit the report (explicit path). Copy the report to `~/Downloads`.
5. Delete the scratch build directory and tarball.

## Post-conditions

- PC1. `git diff --name-only <task-doc-commit>^ HEAD` lists only this task doc and the report.
- PC2. `git status --porcelain` clean; no tarball inside the repo.
- PC3. The report records an outcome line for each of the three uploads.

## Out of scope

- Updating `cran-comments.md` with remote results (follow-up after the emails arrive).
- The CRAN submission itself (Larry).
