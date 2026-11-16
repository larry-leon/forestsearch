SITE WORKFLOW BACK ON CRAN

# Report: site build from the CRAN release again (2026-11-16)

Task: `dev/tasks/TASK_site_build_from_cran_2026-11-16.md`. Host pop-os, branch `master`.
Start SHA `7d277a0f`. R 4.6.2. No site build, article render, `R CMD check`, test suite or
install was run.

## Step 0

Task file committed as `b1dde7f3` (`cmp` identical to `~/Downloads`). Gate passed: branch
`master`. No modified tracked files, and no untracked files. `local::.` occurs once in
`.github/workflows/pkgdown.yaml`, on line 37, in the `extra-packages:` line.

## Step 1

**1.1 The runner will get 0.5.0. Gate passed.** The workflow runs `r-lib/actions/setup-r@v2`
with `use-public-rspm: true`. So the runner resolves CRAN packages from Posit Public Package
Manager, which serves Linux binaries under `__linux__/<Ubuntu codename>/latest`. Queried with
`available.packages(repos = ...)["forestsearch", "Version"]` at 2026-11-16 18:04:06 UTC:

| Repository | URL | forestsearch |
|---|---|---|
| Posit Public Package Manager, source index | `https://packagemanager.posit.co/cran/latest` | 0.5.0 |
| Posit Public Package Manager, Ubuntu 24.04 (noble) binaries | `https://packagemanager.posit.co/cran/__linux__/noble/latest` | 0.5.0 |
| CRAN | `https://cloud.r-project.org` | 0.5.0 |

**1.2 The CRAN release is the code the site was built with. Gate passed.**
`git diff --stat v0.5.0 HEAD -- R DESCRIPTION NAMESPACE inst man` printed nothing.
`--no-renames --name-status` also listed nothing. `data/` exists at neither commit.
`v0.5.0` (tag object `6271c032` → `58a1c27f`) is an ancestor of `HEAD`. The 2026-11-14 dry run
built from `b98bd543`, and `git diff --stat v0.5.0 b98bd543` over the same paths is also empty.

**1.3 For the record.** `git diff --stat b98bd543 HEAD -- vignettes _pkgdown.yml` printed
nothing. No article and no site configuration changed since the 2026-11-14 dry run.

## Step 2: the workflow (commit `2bf79833`)

```diff
@@ -33,8 +33,8 @@ jobs:
 
       - uses: r-lib/actions/setup-r-dependencies@v2
         with:
-          # 2026-11-14: build the site with the repository's forestsearch (0.5.0 is not on CRAN yet); see dev/tasks/TASK_site_build_from_repo_2026-11-14.md
-          extra-packages: any::pkgdown, any::katex, local::.
+          # 2026-11-16: forestsearch 0.5.0 is on CRAN; the site builds with the CRAN release again; see dev/tasks/TASK_site_build_from_cran_2026-11-16.md
+          extra-packages: any::pkgdown, any::katex, cran::forestsearch
           needs: website
 
       - name: Build site
```

Checks:

- Line 36 is byte-identical to the comment line the task gives.
- `local::.` now occurs 0 times and `cran::forestsearch` once.
- The file parses with R `yaml` 2.3.12 and with Python PyYAML 6.0.1. In both, the parsed
  workflow equals `HEAD`'s, apart from the `extra-packages` value.
- The "Build site" step is unchanged:
  `pkgdown::build_site_github_pages(new_process = FALSE, install = FALSE, quiet = FALSE)`.
- The file is still 50 lines, with LF endings and a final newline. The diff is 2 insertions and
  2 deletions.

## This session's commits

```
b1dde7f3 task: site-build-from-cran 2026-11-16 (as received)
2bf79833 ci: build the pkgdown site with the CRAN release (cran::forestsearch) again
(this report's commit)
```

## Final checks (run with this report staged, before its commit)

- This session's commits and the staged report touch only three paths:
  `.github/workflows/pkgdown.yaml`, `dev/tasks/TASK_site_build_from_cran_2026-11-16.md` and
  `dev/reports/REPORT_site_build_from_cran_2026-11-16.md`.
- No modified tracked files besides the staged report.
- This session pushed, pulled, merged, tagged and deployed nothing.
- `refs/remotes/origin/master` moved to `b1dde7f3` ("update by push") at 10:04:18 -0800. That
  was 37 s after the task commit, with a `FETCH_HEAD` update at 10:04:19. This session made
  neither change; GitHub Desktop was running on the host. A push to `master` triggers
  `pkgdown.yaml` on GitHub. That run uses `b1dde7f3`'s workflow, which still has `local::.`.
  `2bf79833` and this report's commit were not on `origin/master` when this was checked.
- The `v0.5.0` tag ref was written at 09:46:18 -0800, before the task file arrived in
  `~/Downloads` at 10:03:03.
