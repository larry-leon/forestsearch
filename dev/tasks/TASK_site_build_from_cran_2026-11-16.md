# TASK — website: build the site with the CRAN release again (2026-11-16)

Session name: `fs-site-build-from-cran`

Host `pop-os`, repository `~/Documents/GitHub/forestsearch`, branch **`master`**.
Larry asked for this on 2026-11-16. Read this whole file first.

## Why

forestsearch 0.5.0 is on CRAN (confirmed by Larry on 2026-11-16 with `available.packages()`),
and the tag `v0.5.0` marks the accepted commit `58a1c27f`. Larry's standing rule is that every
vignette and article on the site is based on the latest CRAN release. On 2026-11-14 the site
workflow was switched to the repository's code (`local::.`, commit `b98bd543`) only because
0.5.0 was not on CRAN yet. This task switches it back.

## Rules

- **Only `.github/workflows/pkgdown.yaml` changes** (plus this task file and the report). No
  article, no `R/`, no test, no `_pkgdown.yml` change.
- **Lean run.** No agents. One short report. No `R CMD check`, no test suite, **no site build
  and no article render**: the full site build passed on 2026-11-14 with this package code
  (`dev/reports/REPORT_site_build_from_repo_2026-11-14.md`); Step 1 confirms that still applies.
- Work on `master`. Commits are path-limited. **No push, pull, merge or tag.** Nothing is
  deployed from this machine.
- No install into any library. Reading a repository index over the network is allowed.
- Gates stop on failure; they never stop to ask. On a stop, change nothing further, write the
  report, and close.

## Step 0 — start

1. Copy this file into `dev/tasks/` and commit it. Record the start SHA.
2. **Gate:** branch `master`; no modified tracked files (untracked files are fine); `local::.`
   occurs exactly once in `.github/workflows/pkgdown.yaml`, in the `extra-packages:` line.

## Step 1 — two checks that replace a new dry run

1. **The runner will get 0.5.0.** Read the workflow to see where the runner resolves CRAN
   packages from (for example `use-public-rspm: true` on `r-lib/actions/setup-r` means Posit
   Public Package Manager). Query that repository's current index, and CRAN's, for the
   forestsearch version, for example:
   `available.packages(repos = "https://packagemanager.posit.co/cran/latest")["forestsearch", "Version"]`
   and the same with `https://cloud.r-project.org`. Record both answers and the repository the
   runner uses.
   **Gate:** the repository the runner uses serves forestsearch 0.5.0. If it still serves an
   older version, stop: the report's first line is `STOPPED: runner repository still serves
   <version>`. (The switch then waits until that repository has caught up.)
2. **The CRAN release is the code the site was built with.** Compare the tag with the current
   tree over the package sources: `git diff --stat v0.5.0 HEAD` limited to `R/`, `DESCRIPTION`,
   `NAMESPACE`, and `data/`, `inst/`, `man/` where they exist.
   **Gate:** no difference. If there is one, stop: first line `STOPPED: package sources differ
   from v0.5.0`, with the `--stat` output.
3. For the record only (no gate): `git diff --stat b98bd543 HEAD` limited to `vignettes/` and
   `_pkgdown.yml`. Say whether any article or the site configuration changed since the
   2026-11-14 dry run.

## Step 2 — the workflow

1. In `.github/workflows/pkgdown.yaml`, in the `extra-packages:` line, replace `local::.` with
   `cran::forestsearch`. Replace the 2026-11-14 comment line above it with:
   `# 2026-11-16: forestsearch 0.5.0 is on CRAN; the site builds with the CRAN release again; see dev/tasks/TASK_site_build_from_cran_2026-11-16.md`
2. Leave the "Build site" step as it is, including `quiet = FALSE`.
3. Confirm the file is still valid YAML and that nothing else in it changed. One commit.

## Step 3 — report and close

1. Write `dev/reports/REPORT_site_build_from_cran_2026-11-16.md`: first line `SITE WORKFLOW
   BACK ON CRAN` or the `STOPPED: ...` line; Step 1's three results; the workflow diff; this
   session's commits. Facts only.
2. Commit the report. Copy it to `~/Downloads`.
3. **Final checks** (report the result): this session's commits touch only
   `.github/workflows/pkgdown.yaml`, this task file and the report; no modified tracked files;
   nothing was pushed, pulled, merged, tagged or deployed.
4. Print the first line and Step 1's results as the closing message.
