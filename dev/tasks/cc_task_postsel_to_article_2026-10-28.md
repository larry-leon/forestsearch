# CC task — `survival-post-selection` moves to the pkgdown articles (2026-10-28)

**Decision (Larry, 2026-10-28).** To cut the 0.4.1 CRAN check time, the post-selection analysis leaves the CRAN vignettes and becomes a pkgdown article. `vignettes/forestsearch.qmd` stays the only CRAN vignette and keeps running its analyses in full (no evaluation switch). Larry's two checks (`devtools::check()` and `devtools::check_win_devel()`) follow this task; they are not part of it.

**Branch:** `feature/vignettes-0.4.1`, starting at `b2e43397` ("dev: handoff — vignettes and 0.4.1 release (2026-10-28)").

**Rules**

- Gates stop on failure — a failed gate is a stop even when the reason is clear. On a stop: write the report with the failing gate and a proposed fix, do the close-out copy (§8), and end.
- Every `git add` names explicit paths. Never stage rendered HTML. Never push, fetch or pull.
- Do not run `R CMD check`, `rcmdcheck`, `devtools::check()` or the test suite.
- No changes under `R/`, `man/`, `DESCRIPTION`, `NAMESPACE`, `inst/`, `tests/`, `quarto/`, `dev/` (apart from this task's own document, record and report), or to any `HANDOFF_*.md`. References found there are recorded for the release task, not edited. Tracked `.html` files are never edited.

**Fixed values**

- Quarto: `/usr/lib/rstudio/resources/app/bin/quarto/bin/quarto`, run with `RSTUDIO_PANDOC=/usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64`.
- Text checks: `$RSTUDIO_PANDOC/pandoc --wrap=none -t plain <file>.html`; count occurrences with `grep -o … | wc -l`, never `grep -c`.
- Article address, written `URL` below: `https://larry-leon.github.io/forestsearch/articles/survival-post-selection.html`
- Record directory: `dev/reports/postsel_to_article_2026-10-28/`. Report: `dev/reports/REPORT_postsel_to_article_2026-10-28.md`. Close-out bundle: `~/Downloads/postsel_to_article_2026-10-28/`.

---

## 1. Commit this task document

Copy `~/Downloads/cc_task_postsel_to_article_2026-10-28.md` to `dev/tasks/` and commit it alone:

`dev: task — post-selection vignette to pkgdown article (2026-10-28)`

## 2. Preconditions

- **G1** On `feature/vignettes-0.4.1`; `git rev-parse HEAD~1` begins `b2e43397`.
- **G2** `git status --porcelain` lists only `?? vignettes/articles/methodology.html` and `?? vignettes/articles/paper_simulations.html`.
- **G3** The caption fix is committed: `git show HEAD:vignettes/survival-post-selection.qmd | md5sum` gives `660b2f87bd26edd0a5e1fdd9801fd09e`.
- **G4** The installed build is still valid: `git diff --quiet 936b6734 HEAD -- R/ DESCRIPTION NAMESPACE inst/` exits 0.

Record the start time (`date -Is`) and `Rscript -e 'cat(future::availableCores())'`.

## 3. Inventory (read-only; write `inventory.md` in the record directory)

1. **References.** `git grep -n -i -e 'survival-post-selection' -e 'post-selection vignette' -e 'companion vignette' -- ':(exclude)*.html'`. For each hit give the file, the line, the whole sentence (or code line), whether it sits in prose or in code, and its class:
   - *edit* — anything under `vignettes/` except `vignettes/survival-post-selection.qmd` itself, plus `README.Rmd` and `README.md`;
   - *pkgdown* — `_pkgdown.yml`;
   - *news* — `NEWS.md` between the `# forestsearch (development version)` heading and the next `# forestsearch ` heading;
   - *record-only* — everything else (`R/`, `man/`, released sections of `NEWS.md`, `dev/`, `HANDOFF_*.md`, `quarto/`, …).
   Also list, as record-only, tracked `.html` files that contain `survival-post-selection` (names only).
2. **Relative paths** in `vignettes/survival-post-selection.qmd`: YAML header (`bibliography`, `csl`, `css`, `include-*`, …), markdown images and links, and file paths in code (`include_graphics()`, `readRDS()`, `source()`, …). Mark each one that resolves from `vignettes/` but not from `vignettes/articles/`.
3. **`n_vig`**: whole-word, case-insensitive matches of `vignette` in the prose of `vignettes/survival-post-selection.qmd` — outside the YAML header, fenced code chunks and inline code.
4. **Rendered leftovers**: entries in `vignettes/` whose names begin `survival-post-selection`, other than the `.qmd`, with `git check-ignore -v` for each.
5. **`.Rbuildignore` coverage** of the paths `vignettes/articles`, `HANDOFF_vignettes_release_0.4.1_2026-10-28.md` and `dev`: a path is excluded when any non-empty, non-`#` line of `.Rbuildignore` matches it with `grepl(pattern, path, perl = TRUE, ignore.case = TRUE)`.

## 4. Move

1. Move each leftover from §3.4 to `~/Downloads/postsel_to_article_2026-10-28/pre_move_render/` (moved, not deleted; never staged).
2. `git mv vignettes/survival-post-selection.qmd vignettes/articles/survival-post-selection.qmd` and commit it alone:

   `docs: move survival-post-selection.qmd to vignettes/articles/ (pkgdown article)`

   - **G5** `git show --name-status --format= HEAD` is exactly one line: `R100`, the old path, the new path.

## 5. Edits

1. **Article note.** In `vignettes/articles/survival-post-selection.qmd`, directly after the line that closes the YAML header, insert a blank line, the note below as a single line, and a blank line:

   ```
   > **A website article.** In forestsearch 0.4.0 this was a CRAN vignette; from 0.4.1 it is an article on the package website, because the full analysis takes a few minutes to run, too long for CRAN's package checks.
   ```

   If `n_vig` ≥ 1, append to the same line (leading space included):

   ```
    The text is otherwise unchanged; where it says "vignette", read "article".
   ```

2. **Article paths.** Rewrite each path marked in §3.2 so it resolves from `vignettes/articles/` (e.g. `references.bib` → `../references.bib`). No other change to the article.
3. **Cross-references (class *edit*).** Each prose hit becomes a link with text `post-selection article` and target `URL`, changing only that sentence and only as much as it needs to read correctly. Examples: "see `vignette("survival-post-selection")`" → "see the [post-selection article](URL)"; "the post-selection vignette" → "the [post-selection article](URL)". For a hit inside a fenced code block of `README.Rmd` or `README.md`: delete that code line and put this line directly after the block:

   ```
   The post-selection analysis is an article on the package website: [post-selection article](URL).
   ```

   - **G6** No hit of class *edit* lies inside a code chunk of a `.qmd`/`.Rmd` under `vignettes/`.
4. **`_pkgdown.yml`.** In the `articles:` index, the entry `survival-post-selection` becomes `articles/survival-post-selection`, same section and position. Navbar and home links to `articles/survival-post-selection.html` stay as they are (the address does not change).
   - **G7** `_pkgdown.yml` has no other occurrence of the slug.
5. **`NEWS.md`, development section only.** Replace the bullet (with its continuation lines) that describes the post-selection vignette by this bullet, keeping the file's bullet marker:

   ```
   - The post-selection analysis is now an article on the package website rather than a CRAN vignette, because the full analysis takes a few minutes to run: <https://larry-leon.github.io/forestsearch/articles/survival-post-selection.html>. It reproduces the GBSG analysis of León and Anderson (2026a), illustrates the full bootstrap at a small number of replicates, and shows `mr_estimates_table()`.
   ```

   - **G8** The development section has exactly one such bullet before the edit.
6. **`vignettes/forestsearch.qmd`.** Replace `about a minute and a half on one core` with `about a minute on one core` (handoff §6 item 11; the render measures 1:14).
   - **G9** The old string occurs exactly once before the edit.
7. **`.Rbuildignore`.** For each path in §3.5 not yet excluded, append its pattern, and only these: `^vignettes/articles$`, `^HANDOFF_.*\.md$`, `^dev$`.

## 6. Verification (all before the §7 commit)

- **G10 purl identity.** For every edited `.qmd`/`.Rmd`: `knitr::purl(<file>, output = <tmp>, documentation = 0)` is byte-identical to the purl of its HEAD version (`git show HEAD:<path>` into a temp file with the same extension). The only lines allowed to differ are §5.2 path rewrites; list any.
- **G11 references.** In `vignettes/**` (the article and all `.html` excluded), `README.Rmd` and `README.md`: zero `vignette("survival-post-selection`, zero `post-selection vignette`, zero `companion vignette` (case-insensitive), and the count of `survival-post-selection` equals the count of `URL`.
- **G12 NEWS.** In the development section, `survival-post-selection` occurs once and only inside `URL`; `git diff HEAD -- NEWS.md` changes no line outside that section.
- **G13 pkgdown.** `Rscript -e 'pkgdown::check_pkgdown()'` (cap 2 min) prints `No problems found`; save the output as `check_pkgdown.txt` in the record directory.
- **G14 timing sentence.** In `vignettes/forestsearch.qmd`: `about a minute and a half` 0 times; `about a minute on one core` once.
- **G15 scope.** `git diff --name-only HEAD` lists only files among: the article, `vignettes/forestsearch.qmd`, `_pkgdown.yml`, `NEWS.md`, `.Rbuildignore`, `README.Rmd`, `README.md`, and the *edit*-class files of §3.1.
- **G16 article render** (cap 6 min; measured 2:07 on two workers). From `vignettes/articles/`: `timeout 360 <quarto> render survival-post-selection.qmd > <record dir>/render_article.log 2>&1`. Exit status 0; `survival-post-selection.html` is newer than the start time; no `survival-post-selection_files/` directory is left; the log has no line matching `citeproc` or `not found` (case-insensitive). Record the wall clock.
- **G17 article text** (plain text of the new HTML):
  - `In forestsearch 0.4.0 this was a CRAN vignette` occurs once;
  - each of `1.44 (0.62, 3.36)`, `0.64 (0.40, 1.03)`, `1.11 (0.45, 2.75)`, `0.66 (0.40, 1.09)`, `1.05 (0.47, 2.34)`, `0.66 (0.41, 1.06)` occurs at least twice;
  - each of `1.88`, `0.83`, `1.29` occurs at least once; `Efron` occurs at least once.
  - Recorded, not gated: the plain-text diff against the pre-move render, if §4.1 moved one.
- **G18 primary render** (cap 4 min; measured 1:14). From `vignettes/`: `timeout 240 <quarto> render forestsearch.qmd > <record dir>/render_primary.log 2>&1`. Exit status 0; the plain text has `about a minute on one core` once and no `about a minute and a half`; if §5.3 edited this file, `post-selection article` is present. Record the wall clock.
- **G19 what CRAN would receive** (cap 5 min). In a temporary directory outside the repo: `R CMD build --no-build-vignettes --no-manual --no-resave-data <repo path>`, then `tar tzf` the tarball. The tarball holds no `HANDOFF_` file, nothing under `dev/` or `vignettes/articles/`, and nothing named `survival-post-selection`. Record the tarball size, its top-level entries and its `vignettes/` entries; flag (record only, not a gate) any top-level entry other than `DESCRIPTION`, `NAMESPACE`, `NEWS.md`, `LICENSE`, `LICENSE.md`, `README.md`, `R`, `man`, `tests`, `vignettes`, `inst`, `data`, `build`. Remove the temporary directory afterwards; nothing in the repo changes.

## 7. Commits and post-conditions

1. Commit the edits — explicit paths, exactly the G15 list; before committing, assert that `git diff --cached --name-only` equals that list and contains no `.html`:

   `docs: post-selection analysis as a pkgdown article — note, cross-references, _pkgdown.yml, NEWS; primary vignette timing sentence`

2. Write the report (§8) and commit the record directory and the report (explicit paths, no `.html`):

   `dev: report — post-selection to pkgdown article (2026-10-28)`

- **G20** `git log --oneline b2e43397..HEAD` shows exactly four commits: task document, rename, edits, report.
- **G21** `git diff --name-only b2e43397 HEAD -- R man DESCRIPTION NAMESPACE inst tests quarto` is empty, and no `HANDOFF_*.md` changed.
- **G22** `git status --porcelain` lists only `??` entries for `vignettes/articles/methodology.html`, `vignettes/articles/paper_simulations.html` and `vignettes/articles/survival-post-selection.html`.
- **G23** Nothing at the top level of `vignettes/` begins `survival-post-selection`.

## 8. Report and close-out

The report, in bullets:

- commits (hash, subject) and a gate table (G1–G23, PASS/FAIL);
- before → after, verbatim: every cross-reference sentence (§5.3), the NEWS bullet, the `_pkgdown.yml` line, any path rewrites, any `.Rbuildignore` lines added;
- record-only references, grouped by file (these go to the release task);
- render wall clocks and `availableCores()`;
- G19: tarball size, flagged top-level entries, `vignettes/` contents.

Close-out (always, after a stop too): copy the report, the record directory, and `vignettes/articles/survival-post-selection.html` and `vignettes/forestsearch.html` (whichever exist) into `~/Downloads/postsel_to_article_2026-10-28/`, which already holds `pre_move_render/`.
