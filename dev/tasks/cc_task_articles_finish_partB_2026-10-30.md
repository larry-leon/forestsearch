# CC task — finish the articles, part B: the ACTG 175 article, the `?forestsearch` examples, the vignette's next steps, NEWS/README/cran-comments, and the release re-check

**File:** `dev/tasks/cc_task_articles_finish_partB_2026-10-30.md` · **Issued:** 2026-10-30 by chat, commissioned by Larry ("Let's finish anything outstanding with the vignettes"; his answers of 10-30: reopen the 0.4.1 tarball for `README.md` and `NEWS.md`; the ACTG 175 article covers 2026b's two ACTG 175 analyses and 2026a's Table S3; "incorporate the 'Two decisions' on the affirmative": fix `?forestsearch` Example 1 and add the Guo–He and ACTG 175 bullets to the vignette's "Where to go next"; Gate C of part A's stage 0 stands at the 2-dp reading).
**Repo:** `forestsearch` · **Branch:** `feature/articles-finish` at `da17be3c` (part A's closing block). Same branch; no new branch.
**Transport:** `~/Downloads` → `dev/tasks/`: this document, and `articles_finish_partB_inputs_2026-10-30.zip` (md5 `dc861f543191ecb2b756c1f730bff53e`), unzipped to `dev/tasks/articles_finish_partB_2026-10-30/`. Filenames may arrive with hyphens or underscores altered; locate by pattern and restore the names.

**This document is the complete, authoritative task.** A gate failure stops the task at that gate: record it, leave the uncommitted work uncommitted, write the report, and stop. Never ask, never work around, except under the one deviation clause of §5.4.

**What chat could not do.** Chat cannot run R or Quarto. The article, the three gate scripts and the edit set were written by reading the package source at `da17be3c` and part A's stage-0 record (`dev/reports/actg175_article_stage0_2026-10-30/`: `calls_2026b/`, `calls_2026a/`, `settings.md`, `package_facts.md`, `ladder.csv`, the compare logs). The renders of §5 are the article's first execution.

---

## ⚠ CATEGORY

**Documentation only; no package code changes behaviour.** Files written:

| file | in the tarball | what changes |
|---|---|---|
| `vignettes/articles/actg175-glm-outcomes.qmd` (new) | no (`^vignettes/articles$`) | the ACTG 175 article |
| `_pkgdown.yml` | no | navbar entry and articles-index entry for it |
| `R/forestsearch_main.R` | yes | **roxygen `@examples` comment and example lines only, inside `\dontrun{}`**: Example 1 calls `mr_estimates_table()` with the bootstrap result it requires; Example 1's "Sec. 5" → "Sec. 6"; Example 2's header names the analysis it runs (2026a, Table S3) |
| `man/forestsearch.Rd` | yes | regenerated from the above by `devtools::document()` |
| `vignettes/forestsearch.qmd` | yes | "Where to go next": two bullets added (Guo–He; ACTG 175), the `?forestsearch` bullet reworded |
| `NEWS.md`, `README.md` | yes | the ACTG 175 article; the examples |
| `cran-comments.md` | no | the Documentation bullet |

Larry reopened the tarball for exactly these: `README.md`, `NEWS.md` (10-30), and the example fix and vignette bullets (10-30, "on the affirmative"). §2.3 checks that the `R/` change touches roxygen lines only. Also written: this document and the inputs under `dev/tasks/`; the record `dev/reports/articles_finish_partB_2026-10-30/`; the report `dev/reports/REPORT_articles_finish_partB_2026-10-30.md`. The article's local HTML render stays on disk, ignored by git (`vignettes/articles/*.html`, part A §1.2). `docs/` (git-ignored) is written by §5.

**Compute, its own line.** Two renders of the new article (pkgdown, then local Quarto), each about 6 min: five MR fits at 5,000 draws on one core (stage-0 measured 64 / 99 / 52 / 3 / 18 s), two declaration calibrations (< 1 s), and bootstraps of B = 25 / 25 / 10 on at most two workers (stage-0: 3.1 / 1.3 / 4.0 s per replicate on one core); cap 30 min per render. `?forestsearch` Example 1 once (a GBSG fit with 5,000 draws on 4 workers and a 20-replicate bootstrap; a few minutes; cap 15 min). `devtools::document()`. `R CMD build` (rebuilds the vignette, about 75 s) and `R CMD check --as-cran --run-donttest` on the built tarball (cap 45 min). `urlchecker` and `spelling`. About an hour in all. Approved by Larry with the scope above.

---

## 0. Provenance and first commit — GATE

1. No `git fetch`. `git rev-parse --abbrev-ref HEAD` must be `feature/articles-finish` and `git rev-parse HEAD` must start `da17be3c`; `git status --porcelain | grep -v '^?? '` must print nothing. Record `git status --porcelain | grep '^?? '` verbatim (**UNTRACKED_0**; untracked files are never staged). *GATE:* otherwise stop.
2. `md5sum` of the zip must equal the header value. Unzip into `dev/tasks/` (it creates `dev/tasks/articles_finish_partB_2026-10-30/`); `(cd dev/tasks/articles_finish_partB_2026-10-30 && md5sum -c MANIFEST.md5)` must report all six files `OK`. *GATE:* otherwise stop.
3. Copy this document to `dev/tasks/cc_task_articles_finish_partB_2026-10-30.md`. `git add` this document and `dev/tasks/articles_finish_partB_2026-10-30/` by path, and commit them alone: `task: finish the articles, part B — task document and inputs`.
4. Record `R --version | head -1`; `packageVersion()` of forestsearch, pkgdown, quarto (the R package), knitr, rmarkdown, roxygen2, devtools, xml2, gt, speff2trial, future; the Quarto CLI and version (as part A §0.5: `quarto` on `PATH`, else the RStudio binary with `QUARTO_PATH` exported for every R call); `pandoc --version | head -1`; whether `urlchecker` and `spelling` are installed. Record `git diff --name-only c2967362 HEAD -- R/ DESCRIPTION NAMESPACE inst/`. *GATE:* forestsearch is `0.4.1`; that diff prints nothing; `speff2trial`, `xml2` and `gt` are installed; **roxygen2 is `8.1.1`** (`DESCRIPTION` carries `Config/roxygen2/version: 8.1.1`, so `document()` then rewrites nothing but the edited page). Otherwise stop — no installs.

Create the record directory `dev/reports/articles_finish_partB_2026-10-30/` with `logs/` and `results/`. Every log below goes to `logs/`, every gate output to `results/`.

## 1. Pre-state — GATE

`md5sum` must print: `R/forestsearch_main.R` `e0f48d94d712765c6c1659dd1d6e2b22`; `man/forestsearch.Rd` `1f9a6708b578c01c23ece464709faa97`; `vignettes/forestsearch.qmd` `2d1c763d1fcd4419afc1a222673cc103`; `NEWS.md` `4b1815355976269f6851f4488c659bc0`; `README.md` `2d3ff9f92188b80092a786eb4993da97`; `cran-comments.md` `29a206ea1b763ff019423eb7c53e44fd`; `_pkgdown.yml` `8b94b6bebb229bd95b79f41367c41870`. `vignettes/articles/actg175-glm-outcomes.qmd` must not exist. *GATE:* otherwise stop.

## 2. Place the article and apply the edit set — GATE

1. `cp dev/tasks/articles_finish_partB_2026-10-30/actg175-glm-outcomes.qmd vignettes/articles/`. *GATE:* md5 `6782b084ee7ceb26e73ba5caf8bbe021`.
2. From the repo root:
   ```
   python3 -I dev/tasks/apply_alignment_edits.py dev/tasks/articles_finish_partB_2026-10-30/articles_finish_partB_edits_2026-10-30.json . --check
   python3 -I dev/tasks/apply_alignment_edits.py dev/tasks/articles_finish_partB_2026-10-30/articles_finish_partB_edits_2026-10-30.json .
   ```
   (the applier committed at `794863f3`, md5 `1056eed3b82f1d889500c662649f6210`; the file keys are repo-relative). *GATE:* both end `RESULT: OK`, with transitions `R/forestsearch_main.R` `e0f48d94→d640bdbf`, `vignettes/forestsearch.qmd` `2d1c763d→037aad7c`, `NEWS.md` `4b181535→d86f61ac`, `README.md` `2d3ff9f9→2d344e97`, `cran-comments.md` `29a206ea→a5cc1963`, `_pkgdown.yml` `8b94b6be→ba8d8365`. The edits, by id (texts in the JSON; `build_partB_edits.py` regenerates it from the `da17be3c` files):
   - **R1** Example 1 header: "2026a, Sec. 5" → "Sec. 6" (Section 6 of 2026a is Applications).
   - **R2** Example 1: after `print(fs_surv)`, a 20-replicate `forestsearch_bootstrap_dofuture()` (4 workers, as the example's fit) and `mr_estimates_table(fs_surv, boot_results = fs_boot, est.scale = "hr")`, with a four-line comment; replaces `mr_estimates_table(fs_surv)`, which has no `boot_results` and stops ("boot_results must contain FSsg_tab and SG_CIs").
   - **R3** Example 2 header: it runs 2026a's forest-search call (Table S3: `effect_neighborhood = 0.10`, `fs.splits = 500`; stage-0 `settings.md` §§1, 3), not 2026b's (0.20), and says that 2026b's analysis selects the same region. Comment only.
   - **V1, V2** the vignette's "Where to go next": a Guo–He bullet (live URL) and an ACTG 175 bullet (it links the **articles index**, because the article's own page is not on the site until the merge); the `?forestsearch` bullet no longer says the examples are 2026b's analyses.
   - **N1** NEWS 0.4.1 Documentation: the article (no URL, for the same reason) and the examples. **M1** README Documentation: the same. **C1** cran-comments: the Documentation bullet. **P1, P2** `_pkgdown.yml`: navbar "Binary and Continuous Outcomes (ACTG 175)" after the Guo–He entry; `articles/actg175-glm-outcomes` in the "Post-Selection Inference" index section.
3. *GATE — roxygen lines only:* `git diff -U0 -- R/forestsearch_main.R | grep -E '^[+-]' | grep -vE '^(\+\+\+|---) ' | grep -vE "^[+-]#'"` prints nothing, and `git diff --numstat -- R/forestsearch_main.R` prints `17	5	R/forestsearch_main.R`.

## 3. The Rd page — GATE

`Rscript -e 'devtools::document()'` → `logs/document.log`. *GATE:* `git status --porcelain -- DESCRIPTION NAMESPACE man/` prints exactly ` M man/forestsearch.Rd`; `git diff --numstat -- man/forestsearch.Rd` prints `17	5	man/forestsearch.Rd` (each roxygen example line is one Rd line); in `man/forestsearch.Rd`, `grep -c 'Anderson (2026a, Sec. 6)'` = 1, `grep -c 'mr_estimates_table(fs_surv, boot_results = fs_boot, est.scale = "hr")'` = 1, `grep -c 'mr_estimates_table(fs_surv)$'` = 0, `grep -c '2026a, Table S3'` = 1. Otherwise stop (do not restore anything).

## 4. `?forestsearch` Example 1 runs — GATE

```
Rscript dev/tasks/articles_finish_partB_2026-10-30/example1_extract.R man/forestsearch.Rd dev/reports/articles_finish_partB_2026-10-30
cd dev/reports/articles_finish_partB_2026-10-30 && timeout 15m Rscript example1.R > logs/example1.log 2>&1; echo "exit $?" >> logs/example1.log
```
`example1_extract.R` writes every example as shipped (`tools::Rd2ex(..., commentDontrun = FALSE)`) to `forestsearch_examples.R`, and Example 1 alone to `example1.R`, with `library(forestsearch)` before it (R CMD check attaches the package; a script must) and recording lines after it. *GATE:* exit 0; the log has `EXAMPLE1_OK`, a `REGION` line, an `MR corrected` line and a `BOOT` line; `example1_mr_table.html` exists and is not empty. Record the three lines and the wall-clock. (Examples 2 and 3 are unchanged code and are not run here.)

## 5. The article — renders — GATE

Record `git status --porcelain` → `logs/status_before_render.txt`.

1. **Site configuration.** `Rscript -e 'pkgdown::check_pkgdown()'` → `logs/check_pkgdown.txt`. *GATE:* `No problems found`.
2. **The site build.** `Rscript -e 'pkgdown::init_site()'`; then `timeout 30m Rscript -e 'pkgdown::build_article("articles/actg175-glm-outcomes", lazy = FALSE)'` → `logs/render_pkgdown.log`, wall-clock recorded; then `Rscript -e 'pkgdown::build_articles_index()'`. *GATE:* `docs/articles/actg175-glm-outcomes.html` exists; `grep -c 'actg175-glm-outcomes.html' docs/articles/index.html` ≥ 1; `grep -c 'Binary and Continuous Outcomes (ACTG 175)' docs/articles/actg175-glm-outcomes.html` ≥ 1 (the navbar).
3. **The local copy, after the site build** (`build_article()` removes a Quarto article's local HTML — part A's finding). `timeout 30m quarto render vignettes/articles/actg175-glm-outcomes.qmd` → `logs/render_local.log`, wall-clock recorded. *GATE:* `vignettes/articles/actg175-glm-outcomes.html` exists and `git check-ignore -v` names the `vignettes/articles/*.html` line of `.gitignore`.
4. **Deviation clause — the article's assembly code only.** The article's analysis calls are fixed: the five `forestsearch()` calls, the two `fs_declaration_calibration()` calls, the bootstrap call, `B`, `n_workers`, the data preparation and the `stopifnot()` checks on the regions. If a render fails, or §6's table gate fails on a **formatting or assembly** defect (a helper that reads an element under another name than the one `fit$mr_inference` carries, a `sprintf` format, a `kable` column count, an inline expression), make the smallest edit to the article's assembly code that restores the intended output, record it verbatim (old line, new line, reason, the error) under "Article deviations", and re-run §5.2–5.3. At most three such edits in all. Never change an analysis call, a published value typed into the article, the gate scripts, or a sentence's meaning. A failure in an analysis call or in a `stopifnot()` on a region, a fourth defect, or a number that differs from its target in §6 stops the task.
5. After the renders, `git status --porcelain` → `logs/status_after_render.txt`. It must equal the before-file (the HTML and `docs/` are ignored); delete any new untracked path a render leaves (e.g. `vignettes/articles/.quarto/`), say so, and compare again.

## 6. The article — output gates — GATE

On each of the two renders:
```
Rscript dev/tasks/articles_finish_partB_2026-10-30/actg175_article_gate.R docs/articles/actg175-glm-outcomes.html dev/reports/articles_finish_partB_2026-10-30/results/gate_pkgdown.txt
Rscript dev/tasks/articles_finish_partB_2026-10-30/actg175_article_gate.R vignettes/articles/actg175-glm-outcomes.html dev/reports/articles_finish_partB_2026-10-30/results/gate_local.txt
```
The script reads only what the page printed. It checks: seven tables; **Table S16** as run equals the published two columns in every cell but the Region row (whose cells must be the package's labels, `cd40 <= 380` and `wtkg <= 86`; `cd40 <= 507` and `gender`); **Table S17** as run equals the published two columns in every cell; **Table S3** as run equals the published table in Method, SG, n (%), Naive and in every Multiplier cell but FS Ĥ, which must read `1.79 (0.63, 5.09)` against the published `1.79 (0.63, 5.08)` (Gate C of part A, Larry 10-30: "2-dp reading is what I meant"); the bootstrap column six finite positive numbers; the cost table three rows of finite numbers; twelve sentences of prose and printed output that carry computed values (`FS and GRF classify 1060 of 1083 patients alike`, `brings it to 1.773`, `field lower bound of 0.686`, … `about half the naive 3.58`), each exactly once; `NaN` and `?@` zero times. It records the count of `Warning`. *GATE:* both files end `RESULT: OK`. Also record, from `logs/render_pkgdown.log` or the page, each fit's printed region line, the bootstrap seconds, and the cost table, in `results/keyvalues.txt`; and the two render wall-clocks. If either render exceeded 10 min, record it as an OPEN ITEM (the article says "about six minutes").

## 7. URLs and spelling

1. If `urlchecker` is installed: `Rscript -e 'urlchecker::url_check()'` → `logs/url_check.txt`. *GATE:* no failing URL in `README.md`, `NEWS.md`, `vignettes/forestsearch.qmd` or `man/*.Rd` (others recorded, not gated). If not installed, record that; §8's `--as-cran` check covers the incoming URL check.
2. If `spelling` is installed: `Rscript -e 'spelling::spell_check_package()'` → `logs/spelling.txt`, and the words it flags on lines this task added (recorded, not gated).

## 8. The release re-check — GATE

1. Build from the working tree's tracked files only, outside the repo, so that no untracked file can enter the tarball:
   ```
   rm -rf /tmp/fs_partB && mkdir -p /tmp/fs_partB/src
   git ls-files --cached | grep -vE '^(dev|quarto|docs|vignettes/articles)/' | tr '\n' '\0' | xargs -0 cp --parents -t /tmp/fs_partB/src
   cd /tmp/fs_partB && R CMD build src > build.log 2>&1
   timeout 45m R CMD check --as-cran --run-donttest forestsearch_0.4.1.tar.gz > check.log 2>&1
   ```
   (`git ls-files --cached` copies the working-tree version of each tracked file, so the edits of §2–3 are in it; the excluded paths are all `.Rbuildignore`d.) Copy `build.log`, `check.log` and `forestsearch.Rcheck/00check.log` to `logs/`.
2. *GATE:* `00check.log` ends `Status: OK` (0 errors | 0 warnings | 0 notes). Record the vignette re-build and examples timings the check prints, beside `cran-comments.md`'s "about 75 s" and "about 40 s"; if either differs by more than half, record it as an OPEN ITEM (Larry updates the sentence; no edit here).
3. The vignette as shipped: `tar -xzf /tmp/fs_partB/forestsearch_0.4.1.tar.gz -C /tmp/fs_partB forestsearch/inst/doc/forestsearch.html`, then `Rscript dev/tasks/articles_finish_partB_2026-10-30/vignette_gate.R /tmp/fs_partB/forestsearch/inst/doc/forestsearch.html dev/reports/articles_finish_partB_2026-10-30/results/gate_vignette.txt`. *GATE:* `RESULT: OK` (the two new bullets and their two links, the reworded `?forestsearch` bullet, and no "examples for this analysis").
4. `tar -tzf` the tarball → `logs/tarball_files.txt`. *GATE:* no path under `forestsearch/vignettes/articles/`, `forestsearch/dev/` or `forestsearch/docs/`.

## 9. Commits — GATE

Before each commit: `git restore --staged .` (removes nothing from the working tree), then `git add` only that commit's paths, then `git diff --cached --name-only` must list only them. Never `git commit -a`. In this order:
1. (made in §0) the task document and inputs;
2. `vignettes/articles/actg175-glm-outcomes.qmd`, `_pkgdown.yml` — `docs: article — the ACTG 175 analyses (2026b S3, 2026a Table S3)`;
3. `R/forestsearch_main.R`, `man/forestsearch.Rd` — `docs: ?forestsearch examples — mr_estimates_table() with its bootstrap template; Example 2's source; a section reference`;
4. `vignettes/forestsearch.qmd`, `NEWS.md`, `README.md`, `cran-comments.md` — `docs: vignette next steps, NEWS, README and cran-comments for the articles`;
5. `dev/reports/articles_finish_partB_2026-10-30/` and `dev/reports/REPORT_articles_finish_partB_2026-10-30.md` — `dev: report — articles finish, part B (2026-10-30)`.

Post-conditions, printed in the report: `git diff --name-only 0bdc2036 HEAD | grep -vE '^(dev/|vignettes/articles/|\.gitignore$)'` prints exactly `NEWS.md`, `R/forestsearch_main.R`, `README.md`, `_pkgdown.yml`, `cran-comments.md`, `man/forestsearch.Rd`, `vignettes/forestsearch.qmd`; of these, the tarball carries `NEWS.md`, `R/forestsearch_main.R`, `README.md`, `man/forestsearch.Rd`, `vignettes/forestsearch.qmd`. `git ls-files 'vignettes/articles/*.html'` prints nothing.

**Report — `dev/reports/REPORT_articles_finish_partB_2026-10-30.md`:** provenance (§0); §1–3 outcomes; Example 1's recorded lines and wall-clock; per render: wall-clock, the gate file; the key values of §6; article deviations (or "none"); URL and spelling outputs; the check's status line and timings; the tarball check of §8.4; the commits; OPEN ITEMS (findings only); a ten-line verdict.

## 10. Closing — GATE

`git log --oneline da17be3c..HEAD` (five commits); `git status --porcelain | grep -v '^?? '` prints nothing; `git status --porcelain | grep '^?? '` equals UNTRACKED_0. Append a closing block to this document under `dev/tasks/` with those outputs and the outcome, and commit it alone: `task: articles finish part B — closing block`. Copy the report, `results/`, `logs/00check.log`, `example1_mr_table.html`, the local render `vignettes/articles/actg175-glm-outcomes.html` and `/tmp/fs_partB/forestsearch_0.4.1.tar.gz` to `~/Downloads/articles_finish_partB_2026-10-30/`. No push, no fetch, no merge.

## 11. After this task (Larry; not CC)

Review; merge `feature/articles-finish` into `master` and push — the pkgdown workflow then renders the new article on the site (under `cran::forestsearch`, 0.4.0, as it renders the other articles); check the page and that the navbar entry resolves; win-builder release and devel on the tarball (rebuilt from `master` if it moved); fill the two win-builder lines of `cran-comments.md`; submit on 10 November.

## 12. Out of scope

`print()`/`summary()` labelling the GLM consistency column "HR" (stage-0 finding; `R/forestsearch_methods.R:277`, `:499`) — a code change, not made here. `mr_estimates_table()` on a binary fit (untested; the article does not call it). Examples 2 and 3's code. Every other article. Pushing and merging.

---

## Closing block (CC, 2026-10-30)

`git log --oneline da17be3c..HEAD`:
```
5e28fa7a dev: report — articles finish, part B (2026-10-30)
cd153dfe docs: vignette next steps, NEWS, README and cran-comments for the articles
75285f8d docs: ?forestsearch examples — mr_estimates_table() with its bootstrap template; Example 2's source; a section reference
17f0c718 docs: article — the ACTG 175 analyses (2026b S3, 2026a Table S3)
2150ba5c task: finish the articles, part B — task document and inputs
```
`git status --porcelain | grep -v '^?? '`: (empty). `git status --porcelain | grep '^?? '`: (empty) = UNTRACKED_0 (empty).

**Outcome:** every gate (§0–§10) passed, with no article deviations. Article gates: `RESULT: OK` on both renders (349 s pkgdown, 328 s local).
Example 1 ran in 57 s. `R CMD check --as-cran --run-donttest`: `Status: OK` (vignette 74 s, examples 44 s).
Execution note: `git ls-files` quoting of two non-ASCII paths under `dev/` and `quarto/` made §8.1's `cp` fail on them. The copied tree was verified identical to the intended set (627 files); see the report.
Report: `dev/reports/REPORT_articles_finish_partB_2026-10-30.md` (four OPEN ITEMS, none blocking). No push, no fetch, no merge.
