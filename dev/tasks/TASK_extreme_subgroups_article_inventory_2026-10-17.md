# TASK: Extreme-subgroups article — Stage 0 inventory against forestsearch 0.4.0

- Date: 2026-10-17
- Repo: `larry-leon/forestsearch`, branch `master`
- Nature: **READ-ONLY inventory.** No edits to any tracked file other than committing this task doc. Nothing in `R/` is touched.
- Output: one bundle file, `~/Downloads/BUNDLE_extreme_subgroups_article_inventory_2026-10-17.md`

## Standing rules (apply throughout)

- Do not run `git fetch`, `git pull` or `git push`.
- Every `git add` names explicit paths. Never stage pre-existing untracked files.
- Do NOT run `R CMD check`, `rcmdcheck`, `devtools::check()`, `devtools::build_vignettes()`, `pkgdown::build_*()`, `quarto render`, `rmarkdown::render()` or the testthat suite.
- Do NOT execute any article/vignette chunk or any forestsearch analysis function.
- Gates are stop-on-failure, not stop-to-ask. If a gate fails, write what you have to the bundle with a `FAILED GATE` section and stop.
- Hard time cap: each step 5 minutes; the whole task 30 minutes.

## Step 0 — Record the task

1. Copy this file from `~/Downloads/` to `dev/tasks/TASK_extreme_subgroups_article_inventory_2026-10-17.md`.
2. `git add dev/tasks/TASK_extreme_subgroups_article_inventory_2026-10-17.md`
3. Commit with the message `dev/tasks: add extreme-subgroups article inventory task (2026-10-17)`.
4. Gate: `git status --porcelain` shows no staged or modified tracked files after the commit.

## Step 1 — Provenance (bundle section 1)

- Output of `git rev-parse HEAD` and `git branch --show-current`.
- `Version:` line from `DESCRIPTION`.
- Gate: Version is `0.4.0` or later. If not, record it under FAILED GATE and stop.

## Step 2 — Locate the article(s) (bundle section 2)

- List every file under `vignettes/`, recursively, with its size.
- Run `grep -ril "extreme"` over `vignettes/` and `_pkgdown.yml`. List the hits.
- Identify the extreme-subgroups article file(s). For each one, record:
  - the path
  - the format (`.qmd` / `.Rmd`) and its `VignetteEngine` line, if any
  - whether it is a pkgdown-only article (under `vignettes/articles/`, or matched by `.Rbuildignore`) or a CRAN vignette
  - its entry in `_pkgdown.yml` (quote the full `articles:` section verbatim)
- List any other article or vignette that links to it.
- Gate: at least one file is identified. If none is found, or the result is ambiguous, list the candidates and stop.

## Step 3 — Article text (bundle section 3)

- The full text of each identified article, verbatim, in a fenced block.

## Step 4 — Call audit against 0.4.0 source (bundle section 4)

Use `Rscript` with base R only (`parse()`, `getParseData()`). Do not load or install forestsearch.

1. From `NAMESPACE`, list the exported functions.
2. Parse every file in `R/` and build a map of top-level function name → formal argument names (does it have `...`: yes/no).
3. Parse the R code chunks of each article (extract them as text only; do not evaluate). For every call to a function that is defined in `R/` or prefixed `forestsearch::`, record:
   - the function and the chunk label/line
   - whether the function is exported, internal-only (and called with `:::`), or missing from `R/`
   - each named argument passed, flagged `UNKNOWN` if it is not among the formals and the function has no `...`
4. Present the result as a table: function | exported? | arguments flagged | locations.
5. Deprecations: grep `R/` for `lifecycle::deprecate_`, `.Deprecated`, `.Defunct` and `deprecat` (case-insensitive). List each function or argument affected, and mark any that the article uses.

## Step 5 — Execution posture (bundle section 5)

For each article, record (from the text only):

- the setup/global chunk options (`eval`, `cache`, `echo`)
- every `readRDS`, `load`, `system.file`, `here::here` or file path it reads, and whether that file exists in the repo (tracked or untracked: note which)
- every chunk whose code contains bootstrap or resampling arguments (for example `B =`, `nb =`, `boots`, `parallel`, `workers`), quoted with its chunk label. Do not run or estimate timings.

## Step 6 — Citations and API notes (bundle section 6)

- Every citation, DOI or reference in the article. Flag the presence or absence of `10.1002/sim.10163`, `arXiv`, `2609.38361`.
- `NEWS.md`: verbatim, all entries from 0.1.0 onward.

## Post-conditions (machine-checkable; record each as PASS/FAIL in bundle section 7)

1. Exactly one new commit since the starting HEAD. It adds only `dev/tasks/TASK_extreme_subgroups_article_inventory_2026-10-17.md` (check with `git show --stat HEAD`).
2. `git status --porcelain` shows no modified tracked files.
3. `~/Downloads/BUNDLE_extreme_subgroups_article_inventory_2026-10-17.md` exists and contains sections 1–7.
4. No file under `vignettes/`, `R/`, `docs/` or `_pkgdown.yml` has a modification time later than the start of this task.

## Final report (to the terminal)

- Five bullets at most: the article path(s), the counts of UNKNOWN arguments, missing functions and deprecated items used, the citation status, and the PASS/FAIL list of post-conditions.
