SITE BUILD PASSES LOCALLY

# Report: site build from the repository's forestsearch (2026-11-14)

Task: `dev/tasks/TASK_site_build_from_repo_2026-11-14.md`. Host pop-os, branch `master`.
Start SHA `9fe0cd12`. R 4.6.2, pkgdown 2.2.1, quarto (R) 1.5.1, Quarto CLI 1.10.18,
rmarkdown 2.32, knitr 1.52.

## Step 0

Task file committed as `c9cfeffd`. Gate passed: branch `master`. No modified tracked files; the
only untracked file is `CRAN-SUBMISSION`. `cran::forestsearch` occurs once in
`.github/workflows/pkgdown.yaml`, on line 36.

## Step 1: the workflow (commit `b98bd543`)

```diff
@@ -33,11 +33,12 @@ jobs:
 
       - uses: r-lib/actions/setup-r-dependencies@v2
         with:
-          extra-packages: any::pkgdown, any::katex, cran::forestsearch
+          # 2026-11-14: build the site with the repository's forestsearch (0.5.0 is not on CRAN yet); see dev/tasks/TASK_site_build_from_repo_2026-11-14.md
+          extra-packages: any::pkgdown, any::katex, local::.
           needs: website
 
       - name: Build site
-        run: pkgdown::build_site_github_pages(new_process = FALSE, install = FALSE)
+        run: pkgdown::build_site_github_pages(new_process = FALSE, install = FALSE, quiet = FALSE)
         shell: Rscript {0}
```

The file is still valid YAML. `yaml::read_yaml()` and Python `yaml.safe_load()` both parse it.
The parsed workflow is identical to the previous one apart from the two changed values.

### Step 1.2: showing Quarto's output

pkgdown 2.2.1 documents a way, and the "Build site" step now uses it (`quiet = FALSE`). The
help pages chain as follows:

- `?build_site_github_pages`, argument `...`: "Additional arguments passed to `build_site()`."
- `?build_site`, argument `quiet`: "If `FALSE`, generate build messages for build functions
  that take `quiet` arguments."
- `?build_articles`, argument `quiet`: "Set to `FALSE` to display output of knitr and pandoc.
  This is useful when debugging."

The source follows the same path. `build_site()` passes `quiet` to `build_articles()`, then to
`build_quarto_articles()`, then to pkgdown's internal `quarto_render()`, and finally to
`quarto::quarto_render(quiet = quiet)`. The hint "Rerun with `quiet = FALSE` to see the full
error message" comes from the quarto R package (`quarto:::wrap_quarto_error()`). It is added
only when `--quiet` was passed to the Quarto CLI.

In the dry run, `quiet = FALSE` printed Quarto's per-file and per-chunk progress, for example
`[ 1/13] articles/actg175-glm-outcomes.qmd`, `processing file: actg175-glm-outcomes.qmd`,
`2/41 [setup]`.

## Step 2: dry run

**Export.** `git archive HEAD` at `b98bd543`, all top-level entries except `quarto/`, into
`/tmp/site060/src`: 3,303 files, equal to `git ls-files -- ':!quarto'`. No tracked
`.gitattributes` exists, so no export-ignore applies.

**Install.** `R CMD INSTALL --no-test-load --library=/tmp/site060/lib /tmp/site060/src`: exit
0, 16.6 s. Nothing was installed into the user or site library. The load check below stands in
for the skipped install-time load test.

**Load check** (in `/tmp/site060/src`, `CI=true R_LIBS=/tmp/site060/lib`):

```
main: forestsearch 0.5.0 from /tmp/site060/lib/forestsearch
worker pid 1405370 (main pid 1405305 ): forestsearch 0.5.0 from /tmp/site060/lib/forestsearch ; .libPaths()[1] = /tmp/site060/lib
nbrOfWorkers: 2
LOAD CHECK PASS
```

The worker ran under `plan(multisession, workers = 2)`. The user library also holds
forestsearch 0.5.0, but `/tmp/site060/lib` is first on `.libPaths()`, and both processes loaded
from it.

**Build.** In `/tmp/site060/src`, `CI=true R_LIBS=/tmp/site060/lib Rscript build_site.R`, which
runs
`pkgdown::build_site_github_pages(new_process = FALSE, install = FALSE, quiet = FALSE, dest_dir = "/tmp/site060/site")`.
The call is wrapped in a `withCallingHandlers()` that only records R warnings to a file and does
not muffle them. The export's `.Rprofile` (`options(pkgdown.internet = FALSE)`) was read, as it
would be on the runner.

- Exit status: **0**.
- Wall clock: **1,635 s (27 min 15 s)**, from 13:37:17 to 14:04:32 PST.
- Phases, from the mtimes of the site files:

  | Phase | Time |
  |---|---|
  | Clean, init, home and reference | 2 min 06 s (to 13:39:23) |
  | `forestsearch-2024.Rmd` | 25 s |
  | Batch `quarto render` of the 13 `.qmd` files | 23 min 37 s (13:39:48 to 14:03:25) |
  | Quarto pages written | to 14:03:53 |
  | News, redirects, search index, checks, `.nojekyll` | 39 s |

**Warnings and errors, verbatim:**

- R warnings signalled in the main process: **0**. The captured-warnings file is empty.
- I searched the whole log, including Quarto's output (`quiet = FALSE` was on), for `warn`,
  `error`, `fail`, `skipping` and `abort` (case-insensitive), `✖` and `! `. The only hits are
  log lines 359 to 361. They come from running the `\donttest` example of `collect_results()`,
  which stops task 3 of 5 on purpose with `stop("simulated failure")`:
  ```
  *** demo: 1 of 5 tasks failed ***
    First error message(s):
       simulated failure 
  ```
- The "Checking for problems" step printed:
  ```
  Missing images in 'README.md': 'reference/figures/logo.png'
  ℹ pkgdown can only use images in 'man/figures' and 'vignettes'
  ```
  `README.md` line 2 references `man/figures/logo.png`, and no `man/figures/` is tracked.
- Inside a page, not in the log: `articles/mrct-consistency-simulations.html` has one
  `cell-output-stderr` block with three lines, `Warning: Ran out of iterations and did not
  converge`. The chunk re-signals the DGM build's warnings on purpose (`# kept by the build,
  shown here`). No other article page has a warning or error output block.

**Article pages written** (15 files: 14 articles and `index.html`). Every page has the pkgdown
marker on line 2 and navbar `title="Released version">0.5.0`:

```
actg175-glm-outcomes.html  biomarker_effects.html  causal_effects_brief_review.html
extreme_subgroups.html  forestsearch-2024.html  forestsearch.html  gbsg-table1.html
guohe-nested-continuum.html  index.html  methodology.html  mrct-consistency.html
mrct-consistency-simulations.html  paper_simulations.html  survival-post-selection.html
treatment_effect_definitions.html
```

**The three named pages.** All exist, have `<!-- Generated by pkgdown: do not edit by hand -->`
on line 2, and print forestsearch 0.5.0. None mentions any 0.4.x version.

| Page | Size (bytes) | Version printed in the body |
|---|---|---|
| `gbsg-table1.html` | 144,415 | session info `forestsearch_0.5.0` |
| `mrct-consistency.html` | 273,240 | `c(forestsearch = as.character(packageVersion("forestsearch")))` prints `"0.5.0"` |
| `mrct-consistency-simulations.html` | 487,212 | the same print, `"0.5.0"`, plus "Passed under forestsearch 0.5.0" |

The site is 31 MB, 1,227 files, kept at `/tmp/site060/site`. The full log is at
`/tmp/site060/build.log` (1,995 lines).

### Step 2.5

Not applicable: the build passed.

### Step 2.6: pages from top-level `.md` files other than `README.md` and `NEWS.md`

```
CLAUDE.html
HANDOFF_efficiency_robustness_release_2026-10-29.html
HANDOFF_vignettes_release_0.4.1_2026-10-27.html
HANDOFF_vignettes_release_0.4.1_2026-10-28.html
HANDOFF_vignettes_release_0.4.1_2026-10-28_v2.html
```

Each also has a `.md` copy beside it in the site root, as `index.md` and `authors.md` do.
`cran-comments.md` is not rendered, because it is on pkgdown's no-render list.

## This session's commits

```
c9cfeffd task: site-build-from-repo 2026-11-14 (as received)
b98bd543 ci: build the pkgdown site with the repository's forestsearch (local::.), quiet = FALSE
(this report's commit)
```

Nothing was pushed, pulled, merged, tagged or deployed.
