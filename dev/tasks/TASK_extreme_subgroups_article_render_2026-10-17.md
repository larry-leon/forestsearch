# TASK: Extreme-subgroups article — Stage 3 render (keep HTML locally)

- Date: 2026-10-17
- Repo: `larry-leon/forestsearch`, branch `master`. The starting HEAD must be `082a5489` or a descendant of it.
- Nature: render one pkgdown article. **No source edits, no `R/` changes.** The only commit is this task doc.
- Report: `~/Downloads/REPORT_extreme_subgroups_article_render_2026-10-17.md`

## Standing rules

- Do not run `git fetch`, `git pull` or `git push`. Every `git add` names explicit paths. Never stage pre-existing untracked files.
- **Do not delete or revert any HTML or site file.** Forbidden:
  - `pkgdown::build_site()`, `build_articles()`, `clean_site()`, `clean_cache()`
  - `git clean`, `git checkout`/`git restore`/`git stash` on anything under `docs/`
  - `rm`/`unlink` of any file under `docs/` or `~/Downloads`
- **Do not commit any rendered output.** That means nothing under `docs/` and no `.html` file.
- Do NOT run `R CMD check`, `rcmdcheck`, `devtools::check()` or testthat. Do NOT install the package.
- Gates are stop-on-failure, not stop-to-ask. On failure, write the report with a `FAILED GATE` section and stop.
- Hard caps: each non-render step 5 minutes; the render 30 minutes (Step 3).

## Step 0 — Record the task

1. Copy this file from `~/Downloads/` to `dev/tasks/TASK_extreme_subgroups_article_render_2026-10-17.md`.
2. `git add` that path.
3. Commit with the message `dev/tasks: add extreme-subgroups article render task (2026-10-17)`.

## Step 1 — Pre-render state (report section 1)

Record:
- the output of `git rev-parse HEAD` and `git status --porcelain`
- the `destination:` setting in `_pkgdown.yml`, or "default (docs)" if there is none. Call the resulting folder `<dest>` below.
- whether `<dest>` is tracked: the output of `git ls-files <dest> | head -5` and `git ls-files <dest> | wc -l`
- the output of `git check-ignore -v <dest>/articles/extreme_subgroups.html`, or "not ignored"
- whether `<dest>/articles/extreme_subgroups.html` exists, with its size and modification time

## Step 2 — Back up the existing HTML

If `<dest>/articles/extreme_subgroups.html` exists, **copy** it (do not move it) to `~/Downloads/extreme_subgroups_before_2026-10-17.html`.

Gate: the copy exists and is byte-identical to the original (`cmp`).

## Step 3 — Render (hard cap 30 minutes)

macOS has no coreutils `timeout`, so run the render from the repo root inside `callr::r()` with a timeout:

```r
t0 <- Sys.time()
callr::r(
  function() pkgdown::build_article("articles/extreme_subgroups", new_process = FALSE),
  timeout = 1800, show = TRUE
)
Sys.time() - t0
```

- Capture the full console output to a temporary file outside the repo.
- Record the measured wall-clock time.
- If the render errors or times out, record the last 60 lines of output under FAILED GATE and stop. Nothing is deleted in either case.

## Step 4 — Save a local copy and inspect (report sections 2–3)

1. **Copy** `<dest>/articles/extreme_subgroups.html` to `~/Downloads/extreme_subgroups_2026-10-17.html`.
   - Gate: the copy is byte-identical (`cmp`).
2. Gate: `<dest>/articles/extreme_subgroups.html` still exists after the copy.
3. From the rendered HTML (strip tags with R or Python; no browser), record:
   - the line containing `Null DGM` and `true overall HR`, verbatim
   - the line containing `Simulations completed`, verbatim
   - the line(s) in the session-info output that mention `forestsearch` and its version
   - the count of lines containing `## Error` or `Error in`. If nonzero, list them; record this but do not fail.
   - the count of lines containing `Warning`. If nonzero, list up to 10.
   - the number of `<img` tags (the figures)
4. List every file created or modified under `<dest>` by the render (with `find <dest> -newer <a marker touched just before Step 3>`).

## Step 5 — Post-conditions (report section 4, PASS/FAIL)

1. Exactly one new commit since the starting HEAD. It adds only the task doc (check with `git show --stat HEAD`).
2. No file under `R/`, `vignettes/`, `man/`, `DESCRIPTION`, `NAMESPACE` or `_pkgdown.yml` is modified (`git status --porcelain` filtered to those paths is empty).
3. `<dest>/articles/extreme_subgroups.html` exists.
4. `~/Downloads/extreme_subgroups_2026-10-17.html` exists and is byte-identical to it.
5. Nothing under `<dest>` is staged (`git diff --cached --name-only` is empty).

Also record the full `git status --porcelain` after the render. Leave any modified or untracked files under `<dest>` exactly as they are.

## Report

Write `~/Downloads/REPORT_extreme_subgroups_article_render_2026-10-17.md` with sections 1–4 above, the render wall-clock time, and any FAILED GATE detail.

In the terminal, print one line only: the report path and overall PASS/FAIL.
