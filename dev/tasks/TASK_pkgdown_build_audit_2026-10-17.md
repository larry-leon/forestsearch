# TASK: Stage A — read-only audit of the full pkgdown build against CRAN 0.4.0

- Date: 2026-10-17
- Repo: `larry-leon/forestsearch`, branch `master`. The starting HEAD must be `377dfb68` or a descendant of it.
- Purpose: assess whether the GitHub Actions site build will run cleanly. That build installs CRAN forestsearch (`cran::forestsearch`) and renders every article, vignette and reference example.
- Nature: **READ-ONLY. No compute: no chunk is evaluated, no example is run, nothing is rendered.** The only commit is this task doc.
- Report: `~/Downloads/REPORT_pkgdown_build_audit_2026-10-17.md`

## Standing rules

- Do not run `git fetch`, `git pull` or `git push`. Run each git command separately. Never stage pre-existing untracked files.
- Do NOT run `R CMD check`, `devtools::check()`, testthat, any `pkgdown::build_*()`, `quarto render`, `rmarkdown::render()`, `example()` or `devtools::run_examples()`. Do not install or update packages.
- Loading the installed namespace with `asNamespace("forestsearch")` is allowed. Base R `parse()`/`getParseData()` only for source files.
- Gates are stop-on-failure, not stop-to-ask.
- Hard cap: 30 minutes for the whole task.

## Step 0 — Record the task

1. Copy this file from `~/Downloads/` to `dev/tasks/TASK_pkgdown_build_audit_2026-10-17.md`.
2. `git add` that path.
3. Commit with the message `dev/tasks: add pkgdown build audit task (2026-10-17)`.

## Step 1 — Baseline (report section 1)

Record:
- `git rev-parse HEAD`
- `packageVersion("forestsearch")`, its `packageDescription()$Repository`, and the CRAN version from `available.packages(repos = "https://cloud.r-project.org")`
- the `development:` block of `_pkgdown.yml`

**Release commit.** Find the most recent commit at which `DESCRIPTION` had `Version: 0.4.0` exactly, using `git log -p --follow -- DESCRIPTION` on the history reachable from HEAD. Call it `<rel>`, and record its hash, date and subject.

Gate: the installed forestsearch is 0.4.0 from CRAN.

## Step 2 — Divergence since the release (report section 2)

- `git diff --stat <rel> HEAD -- R man NAMESPACE DESCRIPTION`
- the full `git diff <rel> HEAD -- NAMESPACE`
- the list of `man/*.Rd` files changed since `<rel>`. For each one, state whether it contains an `\examples` section and whether that section changed (`git diff <rel> HEAD -- <file>` restricted to `\examples`).
- every exported function name in HEAD's `NAMESPACE` that is not in `getNamespaceExports("forestsearch")` of the installed CRAN release

## Step 3 — Inventory of documents the site build renders (report section 3)

1. List every `.Rmd`/`.qmd` under `vignettes/` (recursively), excluding `vignettes/articles/extreme_subgroups.qmd`, which is already verified.
2. Quote the `articles:` section of `_pkgdown.yml`.
3. For each document, record:
   - **a.** Its size, line count, chunk count, and whether all chunks `parse()`.
   - **b.** Its global and per-chunk `eval`, `cache` and `freeze` settings, and how many chunks have `eval = FALSE`.
   - **c. Call audit against the installed CRAN release.** List every call to a forestsearch function, whether exported or via `:::`, with:
     - its named arguments
     - a `MISSING` flag (not in the installed namespace)
     - an `INTERNAL` flag (not exported, and not called with `:::`)
     - an `UNKNOWN` flag for any argument not among the formals of a function without `...`

     Restrict this to chunks that evaluate.
   - **d. Package dependencies.** List every package attached (`library()`/`require()`) or called with `pkg::`. For each, record whether `DESCRIPTION` lists it in Imports/Suggests. Packages outside both may not be installed on the runner.
   - **e. External reads.** Record every `readRDS`, `load`, `read.*`, `fread`, `system.file`, `here::here` or `file.path` call, and every literal path. Say whether each target exists in the repo, tracked or untracked.
   - **f. Network.** Record any `download.file`, `url(`, `httr`/`curl` call, or `https?://` literal inside an evaluated chunk.
   - **g. Compute footprint.** Record every assignment whose name matches `n_sim|nsim|n_sims|B|nb|n_boot|boots|reps|n_rep` (case-insensitive) with a literal value, every `for`/`foreach`/`lapply`/`future_lapply` whose bound is one of those, and every call to `forestsearch()`, `forestsearch_bootstrap_dofuture()`, `run_simulation_analysis()` or another function with `boot` in its name. Do not estimate timings.
   - **h.** Any citation key (`[@...]`) and whether the `bibliography:` field resolves to a file that exists.

## Step 4 — Reference examples (report section 4; static only)

For every `man/*.Rd` file with an `\examples` section:
- Parse the examples, and record the forestsearch calls that fall outside `\dontrun{}`. Treat `\donttest{}` as run, since pkgdown runs it.
- Flag any `MISSING`/`UNKNOWN` against the installed CRAN namespace, as in 3c.

Give totals, then list only the flagged files.

## Step 5 — Verdict table (report section 5)

Give one row per document, plus one row for "reference examples":

| Document | MISSING | INTERNAL | UNKNOWN | missing pkgs | missing files | network | heavy compute (quote) | verdict |
|---|---|---|---|---|---|---|---|---|

The verdict is one of:
- `OK`: no flags.
- `RISK`: a flag that may fail on the runner.
- `FAIL`: certain to fail, for example a `MISSING` function in an evaluated chunk, or a missing file.

## Post-conditions (report section 6)

1. Exactly one new commit since the starting HEAD (the task doc).
2. `git status --porcelain` shows no modified tracked files.
3. Nothing under `docs/` was modified.

## Report

Write `~/Downloads/REPORT_pkgdown_build_audit_2026-10-17.md` with sections 1–6.

In the terminal, print one line only: the report path and the count of `OK`/`RISK`/`FAIL` verdicts.
