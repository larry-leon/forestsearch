# TASK: Replace the extreme-subgroups article with the fixed-baseline Quarto version

- Date: 2026-10-17
- Repo: `larry-leon/forestsearch`, branch `master`. The starting HEAD must be `01e297f3` or a descendant of it.
- Input file: `~/Downloads/extreme_subgroups.qmd`. This is Larry's fixed-baseline document (`extreme_subgroups_fixed.qmd`), adapted for pkgdown in chat.
- Change:
  - Add `vignettes/articles/extreme_subgroups.qmd`.
  - Remove `vignettes/articles/extreme_subgroups.Rmd`. The edits made to it in `082a5489` are superseded.
  - The article name stays `extreme_subgroups`, so `_pkgdown.yml` is unchanged.
- **No `R/` changes. No edits to `DESCRIPTION`, `NAMESPACE`, `_pkgdown.yml` or any other file.**
- Report: `~/Downloads/REPORT_extreme_subgroups_qmd_replace_2026-10-17.md`

## Standing rules

- Do not run `git fetch`, `git pull` or `git push`. Every `git add`/`git rm` names explicit paths. Never stage pre-existing untracked files.
- Do not delete, clean or revert anything under `docs/` or `~/Downloads`. Do not commit any rendered output.
- Do NOT run `R CMD check`, `rcmdcheck`, `devtools::check()`, testthat, `pkgdown::build_site()`, `build_articles()` or `clean_*()`. Do NOT install the package.
- Gates are stop-on-failure, not stop-to-ask. On failure, write the report with a `FAILED GATE` section and stop. Nothing is committed after a failed gate, apart from the Step 0 commit.
- Hard caps: each non-render step 5 minutes; the render 30 minutes.

## Step 0 — Record the task

1. Copy this file from `~/Downloads/` to `dev/tasks/TASK_extreme_subgroups_qmd_replace_2026-10-17.md`.
2. `git add` that path.
3. Commit with the message `dev/tasks: add extreme-subgroups qmd replacement task (2026-10-17)`.

## Step 1 — Provenance (report section 1)

Record:
- the output of `git rev-parse HEAD`
- the `Version:` line from `DESCRIPTION`
- `quarto --version`
- `packageVersion("pkgdown")`
- `packageVersion("forestsearch")` (the installed version)
- the line count and `shasum -a 256` of `~/Downloads/extreme_subgroups.qmd`

## Step 2 — Compatibility gate against the repo source (report section 2)

Use `Rscript --vanilla` with base R only (`parse()`, `getParseData()`). Do not load forestsearch, and do not evaluate any chunk.

1. Extract every `{r ...}` chunk of `~/Downloads/extreme_subgroups.qmd` as text.
   - Gate: every chunk `parse()`s.
2. Using `NAMESPACE` and the top-level functions in `R/`, audit every call in those chunks to a function defined in `R/`, as in the Stage 0 inventory. For each call, record:
   - whether the function is exported
   - its named arguments, flagged `UNKNOWN` if an argument is not among the formals and the function has no `...`
3. Record whether the following exist in `R/` (grep plus formals):
   - `calibrate_k_treat` is exported, and has a `use_ahr` formal
   - `simulate_from_dgm` has a `baseline` formal
   - `generate_aft_dgm_flex` has a `k_treat` formal
   - `df_source` is assigned in the DGM constructor
   - `AHR` is assigned inside `hazard_ratios`
4. **Gate:** 0 missing functions, 0 `UNKNOWN` arguments, and every item in 3 present. If the gate fails, stop: the remedy is a decision for Larry, not part of this task.
5. **Package availability (record only, no gate).** List every package the chunks attach with `library()` or call with `pkg::`. For each one, record whether `DESCRIPTION` lists it under Depends, Imports or Suggests.

## Step 3 — Compare with the repo's copy (report section 3; record only)

If `quarto/extreme_subgroups/extreme_subgroups_fixed.qmd` exists, record:
- `git log -1 --format='%h %ad %s' -- <that path>`
- `diff <that path> ~/Downloads/extreme_subgroups.qmd | head -400`

There is no gate on this step.

## Step 4 — Place the article (report section 4)

1. Copy `~/Downloads/extreme_subgroups.qmd` to `vignettes/articles/extreme_subgroups.qmd`.
   - Gate: the copy is byte-identical (`cmp`).
2. `git rm vignettes/articles/extreme_subgroups.Rmd`
3. Gate: `ls vignettes/articles/extreme_subgroups.*` lists only the `.qmd`.
4. Gate: these fixed-string counts on the placed `.qmd`:

   | String | Expected |
   |---|---|
   | `n_sims_uniform <- 150L` | 1 |
   | `install.packages("forestsearch")` | 1 |
   | `2609.38361` | ≥ 1 |
   | `10.1002/sim.10163` | ≥ 1 |
   | `significan` (case-insensitive) | 0 |
   | `glm-extension` | 0 |
   | `install_github` | 0 |
   | `BRAF` | 0 |
   | `Section 1.3` | 0 |
   | `Section 3.3.1` | 0 |
   | `self-contained` | 0 |

5. `git add vignettes/articles/extreme_subgroups.qmd`
6. Commit with this message:
   ```
   articles/extreme_subgroups: replace with fixed-baseline Quarto version (n_sims 150, CRAN install, references)
   ```
7. Gate: `git show --stat HEAD` lists exactly those two paths (the `.qmd` added and the `.Rmd` deleted).

## Step 5 — Preserve the earlier render, then render (report section 5)

1. **Preserve the earlier render.** If `docs/dev/articles/extreme_subgroups.html` exists, **copy** it and `docs/dev/articles/extreme_subgroups_files/` to `~/Downloads/extreme_subgroups_render_rmd_2026-10-17/`.
   - Gate: `cmp` on the HTML.
2. Touch a marker file outside the repo.
3. Render from the repo root, capturing all output to a file outside the repo:
   ```r
   t0 <- Sys.time()
   callr::r(
     function() pkgdown::build_article("articles/extreme_subgroups", new_process = FALSE),
     timeout = 1800, show = TRUE
   )
   Sys.time() - t0
   ```
   - Record the wall-clock time.
   - On an error or timeout, record the last 80 lines of output under FAILED GATE and stop.
4. Find the output with `find docs -name 'extreme_subgroups.html' -newer <marker>`.
   - Gate: exactly one match. Call its directory `<adir>`.

## Step 6 — Save a local copy and inspect (report section 6)

1. **Copy** `<adir>/extreme_subgroups.html`, and `<adir>/extreme_subgroups_files/` if it exists, to `~/Downloads/extreme_subgroups_render_qmd_2026-10-17/`.
   - Gate: `cmp` on the HTML.
   - Gate: the original files under `<adir>` still exist after the copy.
2. Strip tags from the rendered HTML (Python or R; no browser) and record these lines verbatim:
   - lines containing `Simulations completed`
   - lines containing `AHR (target 0.70)`
   - lines containing `Fixed-baseline panel`
   - lines containing `Parallel plan`
   - lines containing `Throughput`
   - the `forestsearch_` line in the session info
3. Also record:
   - the count of lines containing `Error`. If nonzero, list them.
   - the count of lines containing `Warning`. If nonzero, list up to 10.
   - the number of `<img` tags
   - whether the section headings carry numbers (quote the first three `<h1>`/`<h2>` texts)

## Step 7 — Publishing configuration (report section 7; read-only, record only)

1. List `.github/workflows/`.
2. For each workflow file mentioning `pkgdown`, quote verbatim:
   - its `on:` block
   - any `branches:` lines
   - its deploy step
3. Quote the `development:` block of `_pkgdown.yml`.

## Step 8 — Post-conditions (report section 8, PASS/FAIL)

1. Exactly two new commits since the starting HEAD (the task doc and the article). `git show --stat HEAD` shows only the two article paths.
2. `git status --porcelain -- R vignettes man DESCRIPTION NAMESPACE _pkgdown.yml` is empty.
3. `git diff --cached --name-only` is empty.
4. `~/Downloads/extreme_subgroups_render_qmd_2026-10-17/extreme_subgroups.html` exists and is byte-identical to `<adir>/extreme_subgroups.html`.

## Report

Write `~/Downloads/REPORT_extreme_subgroups_qmd_replace_2026-10-17.md` with sections 1–8 and any FAILED GATE detail.

In the terminal, print one line only: the report path and overall PASS/FAIL.
