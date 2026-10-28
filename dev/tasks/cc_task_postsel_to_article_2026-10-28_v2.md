# CC task v2 — finish moving `survival-post-selection` to the pkgdown articles (2026-10-28)

**Supersedes** `dev/tasks/cc_task_postsel_to_article_2026-10-28.md` (v1) from its G13 stop onward. v1 stopped at G13 because `pkgdown::check_pkgdown()` ran without `RSTUDIO_PANDOC` set ("Pandoc not available"); the fault was in the task document, not in the edits. This document is the complete instruction set for the rest of the task.

**Decision (Larry, 2026-10-28), unchanged.** The post-selection analysis becomes a pkgdown article. `vignettes/forestsearch.qmd` stays the only CRAN vignette and runs its analyses in full. Larry's checks follow this task; they are not part of it.

**Starting state (from the v1 report)**

- Branch `feature/vignettes-0.4.1`. Commits since `b2e43397`: `0c07d261` (v1 task document) and `119b58b4` (rename, R100).
- The v1 §5 edits sit uncommitted in exactly eight files: `.Rbuildignore`, `NEWS.md`, `_pkgdown.yml`, `vignettes/articles/forestsearch-2024.Rmd`, `vignettes/articles/methodology.Rmd`, `vignettes/articles/paper_simulations.Rmd`, `vignettes/articles/survival-post-selection.qmd`, `vignettes/forestsearch.qmd` (the "eight files" below).
- Untracked from v1: the record directory `dev/reports/postsel_to_article_2026-10-28/` and `dev/reports/REPORT_postsel_to_article_2026-10-28.md`. Pre-existing untracked: `vignettes/articles/methodology.html`, `vignettes/articles/paper_simulations.html`.
- v1 gates passed: G1–G12, G14, G15.

**Rules**

- Gates stop on failure — a failed gate is a stop even when the reason is clear. On a stop: update the report with the failing gate and a proposed fix, do the close-out copy (§6), and end.
- Every `git add` names explicit paths. Never stage rendered HTML. Never push, fetch or pull.
- Do not run `R CMD check`, `rcmdcheck`, `devtools::check()` or the test suite.
- No changes under `R/`, `man/`, `DESCRIPTION`, `NAMESPACE`, `inst/`, `tests/`, `quarto/`, `dev/` (apart from this task's documents, record and report), or to any `HANDOFF_*.md`. Tracked `.html` files are never edited.
- **New in v2:** set `RSTUDIO_PANDOC=/usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64` in every shell invocation that runs R, Rscript, Quarto or pandoc (environment variables do not carry over between invocations).
- **New in v2:** never overwrite a v1 record file; new outputs get new names (suffix `_v2` where a v1 file of that name exists).

**Fixed values**

- Quarto: `/usr/lib/rstudio/resources/app/bin/quarto/bin/quarto`.
- Text checks: `$RSTUDIO_PANDOC/pandoc --wrap=none -t plain <file>.html`; count occurrences with `grep -o … | wc -l`, never `grep -c`.
- Article address, written `URL` below: `https://larry-leon.github.io/forestsearch/articles/survival-post-selection.html`
- Record directory: `dev/reports/postsel_to_article_2026-10-28/`. Report: `dev/reports/REPORT_postsel_to_article_2026-10-28.md`. Close-out bundle: `~/Downloads/postsel_to_article_2026-10-28/`.

---

## 1. Commit this task document

Copy `~/Downloads/cc_task_postsel_to_article_2026-10-28_v2.md` to `dev/tasks/` and commit it alone. Stage only that path, and assert that `git diff --cached --name-only` is exactly that one file before committing:

`dev: task v2 — post-selection to pkgdown article, resume from G13 (2026-10-28)`

## 2. Preconditions

- **V1** `git rev-parse HEAD~1` begins `119b58b4` and `git rev-parse HEAD~2` begins `0c07d261`.
- **V2** `git diff --name-only` lists exactly the eight files; `git diff --cached --name-only` is empty.
- **V3** `git status --porcelain` shows only: the eight modified files, `?? dev/reports/postsel_to_article_2026-10-28/`, `?? dev/reports/REPORT_postsel_to_article_2026-10-28.md`, `?? vignettes/articles/methodology.html`, `?? vignettes/articles/paper_simulations.html`.
- **V4** `git diff --quiet 936b6734 HEAD -- R/ DESCRIPTION NAMESPACE inst/` exits 0.

Record this run's start time (`date -Is`).

## 3. One more edit: restore a dropped NEWS sentence

v1's NEWS replacement removed a sentence that describes a different change. In `NEWS.md`, directly after the line that begins `* The post-selection analysis is now an article on the package website`, insert this line:

```
* The methodology and simulation articles on the package website carry a note placing them as the 2024 workflow, with the deprecated arguments removed from their listings.
```

- **V5** The development section (from `# forestsearch (development version)` to the next `# forestsearch ` heading) contains that line exactly once, directly after the post-selection bullet; nothing else in `NEWS.md` changes in this step.

## 4. Verification (all before the §5 commit)

- **G10 purl identity.** For each of the five edited `.qmd`/`.Rmd` files: `knitr::purl(<file>, output = <tmp>, documentation = 0)` is byte-identical to the purl of its HEAD version (`git show HEAD:<path>` into a temp file with the same extension).
- **G11 references.** In `vignettes/**` (the article and all `.html` excluded) and `README.md`: zero `vignette("survival-post-selection`, zero `post-selection vignette`, zero `companion vignette` (case-insensitive), and the count of `survival-post-selection` equals the count of `URL`.
- **G12 NEWS.** In the development section, `survival-post-selection` occurs once and only inside `URL`; `git diff HEAD -- NEWS.md` changes no line outside that section.
- **G13 pkgdown (amended).** `timeout 120 Rscript -e 'pkgdown::check_pkgdown()' > <record dir>/check_pkgdown_v2.txt 2>&1`, with `RSTUDIO_PANDOC` set. Exit status 0 and the output contains `No problems found`.
- **G14 timing sentence.** In `vignettes/forestsearch.qmd`: `about a minute and a half` 0 times; `about a minute on one core` once.
- **G15 scope.** `git diff --name-only HEAD` lists exactly the eight files.
- **G16 article render** (cap 6 min; measured 2:07 on two workers). From `vignettes/articles/`: `timeout 360 <quarto> render survival-post-selection.qmd > <record dir>/render_article.log 2>&1`. Exit status 0; `survival-post-selection.html` is newer than this run's start time; no `survival-post-selection_files/` directory is left; the log has no line matching `citeproc` or `not found` (case-insensitive). Record the wall clock.
- **G17 article text** (plain text of the new HTML):
  - `In forestsearch 0.4.0 this was a CRAN vignette` occurs once;
  - each of `1.44 (0.62, 3.36)`, `0.64 (0.40, 1.03)`, `1.11 (0.45, 2.75)`, `0.66 (0.40, 1.09)`, `1.05 (0.47, 2.34)`, `0.66 (0.41, 1.06)` occurs at least twice;
  - each of `1.88`, `0.83`, `1.29` occurs at least once; `Efron` occurs at least once.
  - Recorded, not gated: the plain-text diff against `~/Downloads/postsel_to_article_2026-10-28/pre_move_render/survival-post-selection.html`.
- **G18 primary render** (cap 4 min; measured 1:14). From `vignettes/`: `timeout 240 <quarto> render forestsearch.qmd > <record dir>/render_primary.log 2>&1`. Exit status 0; the plain text has `about a minute on one core` once, no `about a minute and a half`, and `post-selection article` at least once. Record the wall clock.
- **G19 what CRAN would receive** (cap 5 min). In a temporary directory outside the repo: `R CMD build --no-build-vignettes --no-manual --no-resave-data <repo path>`, then `tar tzf` the tarball. The tarball holds no `HANDOFF_` file, nothing under `dev/` or `vignettes/articles/`, and nothing named `survival-post-selection`. Record the tarball size, its top-level entries and its `vignettes/` entries; flag (record only, not a gate) any top-level entry other than `DESCRIPTION`, `NAMESPACE`, `NEWS.md`, `LICENSE`, `LICENSE.md`, `README.md`, `R`, `man`, `tests`, `vignettes`, `inst`, `data`, `build`. Remove the temporary directory afterwards; nothing in the repo changes.

## 5. Commits and post-conditions

1. Commit the edits — explicit paths, exactly the eight files; before committing, assert that `git diff --cached --name-only` equals the eight files and contains no `.html`:

   `docs: post-selection analysis as a pkgdown article — note, cross-references, _pkgdown.yml, NEWS; primary vignette timing sentence`

2. Update the report (§6) and commit the record directory and the report (explicit paths, no `.html`):

   `dev: report — post-selection to pkgdown article (2026-10-28)`

- **G20** `git log --oneline b2e43397..HEAD` shows exactly five commits: v1 task document, rename, v2 task document, edits, report.
- **G21** `git diff --name-only b2e43397 HEAD -- R man DESCRIPTION NAMESPACE inst tests quarto` is empty, and no `HANDOFF_*.md` changed.
- **G22** `git status --porcelain` lists only `??` entries for `vignettes/articles/methodology.html`, `vignettes/articles/paper_simulations.html` and `vignettes/articles/survival-post-selection.html`.
- **G23** Nothing at the top level of `vignettes/` begins `survival-post-selection`.

## 6. Report and close-out

Keep the v1 report's content as the record of the stop. Change its status line to the final status, add this run's results (V1–V5, G10–G23, marked "v2") to the gate table, and append a section "Resumed under v2", in bullets:

- the restored NEWS line;
- render wall clocks;
- G19: tarball size, flagged top-level entries, `vignettes/` contents.

Close-out (always, after a stop too): copy the report, the record directory, `vignettes/articles/survival-post-selection.html` and `vignettes/forestsearch.html` (whichever exist) into `~/Downloads/postsel_to_article_2026-10-28/`, replacing the earlier copies there; `pre_move_render/` stays.
