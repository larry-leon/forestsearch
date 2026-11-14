# TASK — website: build the site with the repository's forestsearch, and dry-run it (2026-11-14)

Session name: `fs-site-build-from-repo`

Host `pop-os`, repository `~/Documents/GitHub/forestsearch`, branch **`master`**, version 0.5.0.
Larry asked for this on 2026-11-14 ("I want to have the public site current asap"). Read this
whole file first.

## Why

The site workflow `.github/workflows/pkgdown.yaml` installs forestsearch from CRAN
(`cran::forestsearch`, line 36) and builds with `install = FALSE`. CRAN still has 0.4.0, the
articles now need 0.5.0, and every site build since 2026-11-13 has failed in the batch Quarto
render. 0.5.0 is with CRAN for a manual check that may take days. Larry wants the public site
current now, so the site is to build with the repository's own code.

This does not touch the package under review: `.github/` is not in the CRAN tarball.

## Rules

- **Only `.github/workflows/pkgdown.yaml` changes** (plus this task file and the report). No
  article, no `R/`, no test, no `_pkgdown.yml` change.
- **Lean run.** No agents. One short report. No `R CMD check`, no test suite.
- Work on `master`. Commits are path-limited. **No push, pull, merge or tag.** Nothing is
  deployed from this machine.
- No install into the user or site library. Scratch: `/tmp/site060`, scratch build first on
  `R_LIBS=/tmp/site060/lib`.
- Only Step 0 has a stop gate. After it, record a failure and go on.

## Step 0 — start

1. Copy this file into `dev/tasks/` and commit it. Record the start SHA.
2. **Gate (the only stop):** branch `master`; no modified tracked files; line 36 of
   `.github/workflows/pkgdown.yaml` contains `cran::forestsearch` exactly once in the file.

## Step 1 — the workflow

1. In `.github/workflows/pkgdown.yaml`, in the `extra-packages:` line, replace
   `cran::forestsearch` with `local::.`, so that the site build installs the checked-out
   package. Add one comment line above it:
   `# 2026-11-14: build the site with the repository's forestsearch (0.5.0 is not on CRAN yet); see dev/tasks/TASK_site_build_from_repo_2026-11-14.md`
2. **Show the error if the article render fails.** The failed runs print only "Rerun with
   `quiet = FALSE` to see the full error message". If the installed pkgdown (2.2.1) has a
   documented way to make `build_site_github_pages()` print Quarto's output (read its help
   pages), use it in the "Build site" step and quote the documentation line. If it has none,
   leave the step as it is and say so.
3. Confirm the file is still valid YAML. One commit.

## Step 2 — dry run, as the GitHub build will see it

The GitHub runner has only the tracked files. Reproduce that:

1. `git archive HEAD` into `/tmp/site060/src`, leaving out `quarto/` (tracked files only; no
   untracked local files).
2. `R CMD INSTALL` the package from that export into `/tmp/site060/lib`; confirm that the main
   process and one `multisession` worker load forestsearch 0.5.0 from there.
3. In `/tmp/site060/src`, with `CI=true` and `R_LIBS=/tmp/site060/lib`, run what the workflow's
   "Build site" step runs, with the output written to `/tmp/site060/site` (not `docs/`).
   Time limit: 90 minutes. For scale: on 2026-11-13 the reference pages took 2 min 36 s and the
   13 articles 26 minutes on this machine.
4. Record: exit status; wall clock; every warning and error verbatim; the article pages
   written. Confirm `articles/gbsg-table1.html`, `articles/mrct-consistency.html` and
   `articles/mrct-consistency-simulations.html` exist, carry the pkgdown marker and print
   forestsearch 0.5.0.
5. **If the build fails:** identify the failing article, render that article alone with full
   output, and quote the error. Change nothing; this is for the report.
6. List the pages the build makes from top-level `.md` files other than `README.md` and
   `NEWS.md` (for example `CLAUDE.html`, `HANDOFF_*.html`). List only.

## Step 3 — report and close

1. Write `dev/reports/REPORT_site_build_from_repo_2026-11-14.md`: first line SITE BUILD PASSES
   LOCALLY or SITE BUILD FAILS LOCALLY; the workflow diff; Step 1.2's finding; Step 2's record;
   Step 2.6's list; this session's commits. Facts only.
2. Commit the report. Copy it to `~/Downloads`.
3. **Final checks** (report the result): this session's commits touch only
   `.github/workflows/pkgdown.yaml`, this task file and the report; no modified tracked files;
   nothing was pushed, merged, tagged or deployed.
4. Delete `/tmp/site060/lib`. Keep `/tmp/site060/site`.
5. Print the first line and Step 2's key lines as the closing message.
