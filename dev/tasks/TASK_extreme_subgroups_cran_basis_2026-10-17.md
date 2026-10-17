# TASK: Extreme-subgroups article — ease of execution, CRAN basis, main-page publishing

- Date: 2026-10-17
- Repo: `larry-leon/forestsearch`, branch `master`. The starting HEAD must be `1e613340` or a descendant of it.
- Input file: `~/Downloads/extreme_subgroups_v2.qmd`. This replaces `vignettes/articles/extreme_subgroups.qmd`.

## Standing directives behind this task (Larry, 2026-10-17)

- **Ease of execution.** Any user should be able to run the article as simply as possible.
- **Main page.** The site's main page must update as each article is revised.
- **CRAN basis.** Every vignette/article is based on the latest CRAN release of forestsearch (now 0.4.0). That applies both locally and on the GitHub Actions site build.

## What changes

1. `vignettes/articles/extreme_subgroups.qmd`:
   - A visible **User settings** chunk with `n_sims_uniform` and `n_workers_user`, plus a check that forestsearch is at least 0.4.0.
   - A step-by-step **Running this article yourself** section, which replaces the Installation section.
   - The parallel setup now honours `n_workers_user`.
2. `_pkgdown.yml`: `development: mode: auto` → `mode: release`, so a push publishes to the main page instead of `/dev/`.
3. `.github/workflows/pkgdown.yaml`: the site build installs forestsearch **from CRAN** (`cran::forestsearch`) instead of from the repo source (`local::.`). The articles and reference examples then run against the CRAN release.

**No `R/`, `DESCRIPTION`, `NAMESPACE` or `man/` changes.**

- Report: `~/Downloads/REPORT_extreme_subgroups_cran_basis_2026-10-17.md`

## Standing rules

- Do not run `git fetch`, `git pull` or `git push`. Every `git add` names explicit paths. Never stage pre-existing untracked files.
- Do not delete, clean or revert anything under `docs/` or `~/Downloads`. Do not commit rendered output.
- Do NOT run `R CMD check`, `rcmdcheck`, `devtools::check()`, testthat, `pkgdown::build_site()`, `build_articles()` or `clean_*()`. Do NOT install or update any package.
- Gates are stop-on-failure, not stop-to-ask. On failure, write the report with a `FAILED GATE` section and stop.
- Hard caps: each non-render step 5 minutes; the render 30 minutes.

## Step 0 — Record the task

1. Copy this file from `~/Downloads/` to `dev/tasks/TASK_extreme_subgroups_cran_basis_2026-10-17.md`.
2. `git add` that path.
3. Commit with the message `dev/tasks: add extreme-subgroups CRAN-basis task (2026-10-17)`.

## Step 1 — Preconditions (read-only; report section 1)

1. **Installed forestsearch is the latest CRAN release.** Record:
   - `packageVersion("forestsearch")`
   - `packageDescription("forestsearch")$Repository`
   - `available.packages(repos = "https://cloud.r-project.org")["forestsearch", "Version"]`

   Gate: Repository is `"CRAN"`, and the installed version equals the CRAN version.
2. **Workflow.** Quote `.github/workflows/pkgdown.yaml` verbatim in full.
   - Gate: it uses `r-lib/actions/setup-r-dependencies@v2`.
   - Gate: its `extra-packages:` value contains `local::.` exactly once.
   - Gate: it does not set `install-quarto: false`. Version `@v2` installs Quarto automatically when the repo contains a `.qmd` file.
3. **`_pkgdown.yml`.** Gate: it contains the two lines `development:` / `  mode: auto` exactly once.
4. **pkgdown behaviour (record only).** Record:
   - `deparse(formals(pkgdown::build_site)$devel)`
   - every line of `deparse(body(pkgdown::build_site_github_pages))` that mentions `devel` or `install`

## Step 2 — Audit the article against the installed CRAN release (report section 2)

Do not evaluate any chunk. Loading the installed namespace with `asNamespace("forestsearch")` is allowed.

1. Extract the `{r}` chunks of `~/Downloads/extreme_subgroups_v2.qmd`.
   - Gate: all chunks parse.
2. For every call to a function in the installed forestsearch namespace, check two things:
   - the function is in `getNamespaceExports("forestsearch")`
   - its named arguments are in `formals()`, or the function has `...`

   Gate: 0 missing or non-exported functions, and 0 `UNKNOWN` arguments.
3. Gate: in the installed CRAN package:
   - `calibrate_k_treat` has `use_ahr`
   - `simulate_from_dgm` has `baseline`
   - `generate_aft_dgm_flex` has `k_treat`
   - deparsing the namespace's functions finds both `df_source` and `AHR`
4. Record only: for each of survival, ggplot2, dplyr, patchwork, future, doFuture, foreach and weightedsurv, check whether `packageDescription("forestsearch")$Imports` lists it. The article text says installing forestsearch installs all of them.

## Step 3 — Place the article (report section 3)

1. Copy `~/Downloads/extreme_subgroups_v2.qmd` to `vignettes/articles/extreme_subgroups.qmd`.
   - Gate: `cmp` byte-identical.
2. Gate: these fixed-string counts on the placed file:

   | String | Expected |
   |---|---|
   | `n_sims_uniform <- 150L` | 1 |
   | `n_workers_user <- NULL` | 1 |
   | `## Running this article yourself` | 1 |
   | `raw.githubusercontent.com/larry-leon/forestsearch/master/vignettes/articles/extreme_subgroups.qmd` | 1 |
   | `utils::packageVersion("forestsearch") < "0.4.0"` | 1 |
   | `pak::pak` | 0 |
   | `## Installation` | 0 |
   | `self-contained:` | 0 |

3. Gate: `git diff --stat` touches only `vignettes/articles/extreme_subgroups.qmd`.
4. `git add vignettes/articles/extreme_subgroups.qmd`
5. Commit with the message `articles/extreme_subgroups: user settings chunk, run-it-yourself steps, CRAN 0.4.0 check`.

## Step 4 — Publishing configuration (report section 4)

1. In `_pkgdown.yml`, replace `  mode: auto` (the line under `development:`) with `  mode: release`.
2. In `.github/workflows/pkgdown.yaml`, within the `extra-packages:` value, replace `local::.` with `cran::forestsearch`. Leave every other line unchanged.
3. Gate: `git diff` shows exactly those two one-line changes.
4. `git add _pkgdown.yml .github/workflows/pkgdown.yaml`
5. Commit with this message:
   ```
   site: publish to the main page (mode: release); build against CRAN forestsearch
   ```

## Step 5 — Render locally against the installed CRAN release (report section 5)

1. Touch a marker file outside the repo.
2. Render exactly as in the previous task: `callr::r(function() pkgdown::build_article("articles/extreme_subgroups", new_process = FALSE), timeout = 1800, show = TRUE)`.
   - Record the wall-clock time.
   - On an error or timeout, record the last 80 lines of output under FAILED GATE and stop.
3. Find the output with `find docs -name 'extreme_subgroups.html' -newer <marker>`.
   - Gate: exactly one match. Call its directory `<adir>`. With `mode: release`, expect `docs/articles`.
4. **Copy** the HTML and `extreme_subgroups_files/` to `~/Downloads/extreme_subgroups_render_qmd_v2_2026-10-17/`.
   - Gate: `cmp` on the HTML.
   - Gate: the originals still exist after the copy.
5. Strip tags from the HTML and record verbatim the lines containing:
   - `n_workers_user`
   - `Parallel plan`
   - `Simulations completed`
   - `AHR (target 0.70)`
   - `Throughput`
   - `forestsearch_`
6. Also record:
   - counts of lines containing `Error` and `Warning`
   - the number of `<img` tags
   - whether the text `Running this article yourself` appears

## Step 6 — Post-conditions (report section 6, PASS/FAIL)

1. Exactly three new commits since the starting HEAD: the task doc, the article, and the site configuration.
2. `git status --porcelain -- R vignettes man DESCRIPTION NAMESPACE _pkgdown.yml .github` is empty.
3. `git diff --cached --name-only` is empty.
4. The `~/Downloads` render copy is byte-identical to `<adir>/extreme_subgroups.html`.
5. `session info` in the render shows `forestsearch_0.4.0`.

## Report

Write `~/Downloads/REPORT_extreme_subgroups_cran_basis_2026-10-17.md` with sections 1–6 and any FAILED GATE detail.

In the terminal, print one line only: the report path and overall PASS/FAIL.
