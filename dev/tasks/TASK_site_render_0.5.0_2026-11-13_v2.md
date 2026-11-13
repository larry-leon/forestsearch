# TASK — render every pkgdown article under forestsearch 0.5.0, for Larry's review (2026-11-13), v2

Supersedes `TASK_site_render_0.5.0_2026-11-13.md`. The only change from v1: Step 5b (compare
each new render with the rendered `.html` already on disk, where one exists) and the matching
lines in Step 6. Everything else is the same.

Session name: `fs-site-render-0.5.0`

Host `pop-os`, repository `~/Documents/GitHub/forestsearch`. Larry wants to see every pkgdown
vignette and article rendered with 0.5.0 before he submits to CRAN. Read this whole file first.

## Purpose

Build the pkgdown site, with every article, from the current code in a **scratch copy**, and
put the result where Larry can open it in a browser. Nothing is published, and the repository's
own files and any rendered pages already on disk are left exactly as they are. Where an
article already has a rendered `.html` on disk, compare the new render with it.

## Rules

- **Lean run.** The main session does everything. No agents. One short report.
- **No `R CMD check`, no test suite.**
- **The repository is read-only**, apart from two path-limited commits on `cran-prep`: this
  task file and the report. No article source is edited. No file is written into the working
  tree by a render.
- **Packages.** No install into the user or site library (it holds forestsearch 0.4.1.9002,
  and article workers load whatever is installed). The 0.5.0 build goes into
  `/tmp/site050/lib` only, first on `R_LIBS` for every R process.
- **Scratch:** `/tmp/site050`. **Site output for Larry:** `~/fs050_site/` (outside every
  repository; kept).
- **One article's failure does not stop the others.** Record it and go on.
- **No push, pull, merge or tag. Nothing is deployed.**
- fs-post-selection and fs-glms-interpretable are not touched.

## Step 0 — start

1. Copy this file into `dev/tasks/` and commit it on `cran-prep`.
2. **Gate (the only stop):** branch `cran-prep`; `DESCRIPTION` says `Version: 0.5.0`.
3. Record: HEAD of `cran-prep` and of `master`; whether `R/`, `man/`, `vignettes/`, `inst/`,
   `DESCRIPTION`, `NAMESPACE`, `_pkgdown.yml`, `README.md` and `NEWS.md` are identical between
   the two; `git status --porcelain`.

## Step 1 — inventory (read-only)

1. List every page the site builds as an article: the CRAN vignette(s) under `vignettes/` and
   everything under `vignettes/articles/`, with format (`.qmd` / `.Rmd`) and whether
   `_pkgdown.yml` lists it.
2. Find how the site was last built here (the most recent pkgdown or article task and report
   under `dev/tasks/` and `dev/reports/`) and use the same route and tool versions. Record the
   route in one or two lines.
3. For each article, record the last render wall clock if a report gives one.

## Step 2 — scratch copy and library

1. Copy the working tree to `/tmp/site050/src`, leaving out `.git/`, `quarto/`,
   `legacy_vignette/` and any existing built-site directory. Untracked local files that the
   articles read are copied, so that each render behaves as it would on this machine.
2. `git archive` HEAD's package paths and `R CMD INSTALL` them into `/tmp/site050/lib`.
3. **Confirm** (record, and stop only if it fails): with `R_LIBS=/tmp/site050/lib`, the main
   process and one `multisession` worker load forestsearch 0.5.0 from `/tmp/site050/lib`.

## Step 3 — home, reference and news pages

Build these once into `/tmp/site050/site`, with the scratch library and without re-installing
the package. Record the wall clock and any error or warning.

## Step 4 — the articles

- Render **each article separately**, so that one failure or time-out does not stop the rest.
- Order: shortest recorded render time first; articles with no recorded time last.
- Each article runs as written: its own seeds, simulation counts and worker counts.
- **Time limits (safety stops):** 60 minutes per article; 6 hours in total from the start of
  this step. On the per-article limit, stop that article and go on. On the total limit, stop
  and report what is left as NOT RUN.
- Record for each article: status (`OK`, `ERROR`, `TIMEOUT`, `NOT RUN`), wall clock, the error
  text and the call that raised it, and the number of warnings in the render log.
- If an article needs a file that is not in the copy, record the file and the article.

## Step 5 — checks on the rendered pages (records only)

For each rendered article page:

1. the number of R error outputs and of R warning outputs visible on the page (say how they
   were detected), with the first one of each quoted;
2. whether the page or its source names any of the four functions removed in part 3 (their
   names are in `dev/reports/REPORT_cran_prep_docs_legacy_2026-11-12.md`);
3. the forestsearch version the page prints, if it prints one.

## Step 5b — compare with the renders already on disk (records only)

Some articles have an earlier rendered `.html` in `vignettes/` or `vignettes/articles/`
(tracked or untracked); some have none. Never modify, move or delete those files.

1. List every `.html` there. Pair each with its article source by name. List any `.html`
   with no matching source, and every article with no earlier render (nothing to compare).
2. For each pair, record the earlier file's modification time, the forestsearch version it
   prints (if any), and the date of the last commit to the article's source. If the source
   was committed after the earlier render was made, flag the pair `source edited since`:
   differences may then come from the text, not from the package.
3. The two pages differ in layout (site template, time stamps, embedded images), so compare
   **content**, extracted the same way from both pages:
   - every table, cell by cell (tables paired in page order, and by caption where there is
     one);
   - every block of printed R output;
   - the section headings;
   - the number of figures.
4. For each pair, report: tables compared and tables that could not be paired; cells changed
   of cells compared; the largest numeric change; output blocks changed of blocks compared;
   whether the headings and the figure count agree. Mark cells that are plainly timings or
   host details (seconds, minutes, cores, dates) and count them separately.
5. Classify each pair as exactly one of:
   - `SAME` — every compared table cell and output block is identical, timings aside;
   - `NUMBERS MOVED` — the structure agrees and some values differ;
   - `STRUCTURE DIFFERS` — tables or sections do not line up, so the pages are not comparable
     cell by cell (say what differs);
   - `NEW RENDER FAILED` — there is nothing new to compare.
6. For each pair that is not `SAME`, write the full list of differences (earlier value, new
   value, where on the page) to `~/fs050_site/_compare/<article>.txt`. The report quotes at
   most the first ten per article. State facts only: no view on which render is right.

## Step 6 — hand over and close

1. Copy `/tmp/site050/site` to `~/fs050_site/` (replace only that directory's own contents).
   Put a **copy** of each earlier render in `~/fs050_site/_previous/`, so that Larry can open
   the earlier and the new page side by side.
2. Write `dev/reports/REPORT_site_render_0.5.0_2026-11-13.md`. Short:
   - **first lines:** the two paths Larry opens, `~/fs050_site/index.html` and
     `~/fs050_site/articles/index.html`, and a one-line count (articles rendered, failed,
     timed out, not run);
   - **top table:** one row per article — source file, status, wall clock, warnings in the
     log, visible errors and warnings on the page, removed-function names found, the direct
     path of the rendered page;
   - **comparison table (Step 5b):** one row per article that has an earlier render — earlier
     render date, `source edited since` or not, cells changed of cells compared (timings
     counted separately), largest numeric change, output blocks changed, classification, and
     the paths of the earlier copy and the new page; then the articles with no earlier
     render;
   - the pages most worth looking at first (any with an error, a visible warning, a removed
     name, or a classification other than `SAME`), in bullets;
   - Step 0's record, the build route, and the Step 3 result;
   - every error text. Facts only, no recommendation.
3. Commit the report (path-limited) and copy it to `~/Downloads`.
4. **Final checks** (report the result): `git status --porcelain` is as in Step 0; this
   session's two commits touch only this task file and the report; nothing was pushed, merged
   or tagged.
5. Delete `/tmp/site050/lib`. Keep `/tmp/site050/site` and `~/fs050_site/` (with `_previous/`
   and `_compare/`).
6. Print the first lines and the top table as the closing message.
