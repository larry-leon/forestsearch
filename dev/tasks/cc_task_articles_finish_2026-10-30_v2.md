# CC task — finish the articles, part A: Quarto conversion, local HTML, README and NEWS, DINA note; ACTG 175 article stage 0 — v2

**File:** `dev/tasks/cc_task_articles_finish_2026-10-30_v2.md` · **Issued:** 2026-10-30 by chat, commissioned by Larry (14:19: "Let's finish anything outstanding with the vignettes. I don't know why we would postpone anything for after the release"; 14:2x, his answers: reopen the 0.4.1 tarball for `README.md` and `NEWS.md`; article HTML kept on disk and ignored by git; convert all five aligned articles to `.qmd`; the ACTG 175 article covers 2026b's two ACTG 175 analyses **and** 2026a's Table S3).
**Repo:** `forestsearch` · **Branch:** `feature/articles-finish`, created from `master` at `0bdc2036` (the merge of both article branches). **All remaining article work uses this one branch** until Larry merges it; the ACTG 175 article task (part B, next) continues on it.
**Transport:** `~/Downloads` → `dev/tasks/`: this document, and `articles_finish_inputs_2026-10-30.zip` (md5 `25f6119061f6cb0586df2c81836ca4ee`), unzipped to `dev/tasks/articles_finish_2026-10-30/`. Filenames may arrive with hyphens or underscores altered; locate by pattern and restore the names.

**v2 supersedes v1** (`dev/tasks/cc_task_articles_finish_2026-10-30.md`, committed with the inputs in §0). Both parts of v1 stopped at a gate, with all their work done and checked on disk:
- **Part 1 stopped at §1.9.** Commit 1 (`8cde6811`) also took the five `.Rmd → .qmd` renames, because §1.3's `git mv` stages a rename at once and v1 did not say to unstage it. That was a defect of v1, not of CC's work. Larry's disposition (2026-10-30): `git reset --soft HEAD~1`, then make commits 1–4 exactly as §1.9 lists them. §1.9 below now says how.
- **Part 2 stopped at Gate C.** Under the installed 0.4.1, 2026a's forest-search call gives the multiplier interval (0.63, **5.09**) where Table S3 prints (0.63, **5.08**). Table S3's payload was built under forestsearch 0.2.0, and v1's stop report says it was rendered from text committed in `43b051b6`, not from the commit the manifest names. Gates A and B and Table S17 reproduced in full. v2 adds a check that the transplant matches what was rendered (§2.2.4), restates Gate C with a tolerance on the interval limits only (§2.4), and runs the remaining rungs.

**This document is the complete, authoritative task.** What differs from v1: this paragraph, §1.3's note, §1.9, §2.2.4, Gate C in §2.4, and the resume rule in §3. Nothing that passed under v1 is re-run: not the renders, the text gates, the README example, Gates A and B, or Table S17.

**Two independent parts.** Part 1 edits the articles, README and NEWS. Part 2 is a read-only inventory and measurement for the ACTG 175 article; it writes only under `dev/`. A gate failure stops **that part**: record it, leave that part's uncommitted work uncommitted, and run the other part. Never ask, never work around.

**What chat could not do.** Chat cannot run R or Quarto. The five `.qmd` files and the README were written by reading. The renders in §1.6 are their first execution, and the deviation clause there covers render defects of the conversion only.

---

## ⚠ CATEGORY

**No package code is modified:** nothing under `R/`, `man/`, `tests/`, `inst/`, `DESCRIPTION`, `NAMESPACE`, `vignettes/forestsearch.qmd`, `_pkgdown.yml`, `quarto/`. **Two tarball files change, `README.md` and `NEWS.md`** (Larry reopened the 0.4.1 tarball for them, 2026-10-30). The release re-check (`R CMD check --as-cran` on the built tarball, `cran-comments.md`) is **not** run here: it runs once at the end of part B, after part B's NEWS line, on the final state. `.gitignore` also changes; `.Rbuildignore` excludes it from the tarball.

**Writes:** this document and the unzipped inputs under `dev/tasks/`; `.gitignore`; the five articles renamed `.Rmd` → `.qmd` under `vignettes/articles/`; `vignettes/articles/survival-post-selection.qmd` (one sentence); `README.md`; `NEWS.md`; the records `dev/reports/articles_finish_2026-10-30/` and `dev/reports/actg175_article_stage0_2026-10-30/`; the reports `dev/reports/REPORT_articles_finish_2026-10-30.md` and `dev/reports/REPORT_actg175_article_stage0_2026-10-30.md`. **Article HTML is rendered in place and stays on disk, ignored by git** (§1.2) — never staged again. `docs/` (pkgdown's local build, git-ignored) is written by §1.6 and stays local.

**Compute, its own line.** Part 1: six Quarto renders in place and six `pkgdown::build_article()` builds — five conversions (methodology 4 s, paper_simulations 47 s and the rest under 3 s under R Markdown on 10-30) and `survival-post-selection` (minutes; it runs three fits and a small bootstrap) — plus one ~45-second run of the README's example; about 15–20 min in all, cap 20 min per render. Part 2: the ladder of §2.4 — five MR fits, two declaration calibrations and 15 bootstrap replicates, all `plan = "sequential"`; unmeasured for these settings (the GBSG fits ran 37 / 8 / 25 s); cap 20 min per rung, 2 h for the ladder. Approved by Larry with the scope above.

---

## 0. Provenance and first commit — GATE

1. No `git fetch`. `git rev-parse --abbrev-ref HEAD` must be `master`; `git merge-base --is-ancestor 0bdc2036 HEAD` must succeed; `git status --porcelain | grep -v '^?? '` must print nothing. Record HEAD (**BASE_SHA**) and `git status --porcelain | grep '^?? '` verbatim (untracked files are never staged). *GATE:* otherwise stop the whole task.
2. `git switch -c feature/articles-finish`.
3. `md5sum` of the zip must equal the header value. Unzip into `dev/tasks/` (it creates `dev/tasks/articles_finish_2026-10-30/`); `(cd dev/tasks/articles_finish_2026-10-30 && md5sum -c MANIFEST.md5)` must report every file `OK`. *GATE:* otherwise stop.
4. Copy this document to `dev/tasks/cc_task_articles_finish_2026-10-30.md`. `git add` this document and `dev/tasks/articles_finish_2026-10-30/` by path, and commit them alone: `task: finish the articles, part A — task document and inputs`.
5. Record `R --version | head -1`; `packageVersion()` of forestsearch, pkgdown, quarto (the R package), rmarkdown, knitr; the Quarto CLI and its version (`quarto --version` on `PATH`, else `/usr/lib/rstudio/resources/app/bin/quarto/bin/quarto --version` — say which; if `quarto` is not on `PATH`, export `QUARTO_PATH` to that binary for every R call below); `pandoc --version | head -1`; `python3 --version`. Record `git diff --name-only c2967362 HEAD -- R/ DESCRIPTION NAMESPACE inst/` (must print nothing: the installed 0.4.1 is the release's code). *GATE:* forestsearch is `0.4.1` and that diff is empty; otherwise stop.

---

## Part 1 — the articles, README and NEWS

### 1.1 Pre-state — GATE

`md5sum` must print: `vignettes/articles/methodology.Rmd` `5d79611e4a538e2b54aadaa151fd27d1`, `paper_simulations.Rmd` `6154d7ea3b2fd8df2b831b1ae5524443`, `treatment_effect_definitions.Rmd` `2c46b87fccf0114413d602c4b74e618e`, `biomarker_effects.Rmd` `42a79e8cbc708576fd5950929a95c4b3`, `causal_effects_brief_review.Rmd` `94caca1f6166919499fce0f1d6b79aaa` (the aligned texts of `a88f8ea4`); `README.md` `915aa3504518e3a6c41e2200f7fd88ff`; `NEWS.md` `29aa87839a6aebab4c67479302d52d95`; `vignettes/articles/survival-post-selection.qmd` `90b5128fa0320bbcf7e001fb49e54bd1`. Record `git ls-files 'vignettes/articles/*.html'` (expected: `forestsearch-2024.html`, `guohe-nested-continuum.html`, `methodology.html`, `paper_simulations.html`, `survival-post-selection.html`) and `md5sum vignettes/articles/*.html` → `logs/html_before.txt`. *GATE:* any md5 differs → stop Part 1.

### 1.2 Article HTML: on disk, ignored by git — GATE

1. Append one line to the root `.gitignore`: `vignettes/articles/*.html` (with a comment line above it: `# rendered article HTML: kept on disk, never tracked (Larry, 2026-10-30)`).
2. `git rm --cached` each tracked `vignettes/articles/*.html` listed in §1.1. Nothing is deleted from disk.
3. *GATE:* `git ls-files 'vignettes/articles/*.html'` prints nothing; every file of `logs/html_before.txt` still exists with the same md5; `git check-ignore -v vignettes/articles/<each>.html` names the new `.gitignore` line; `git status --porcelain` lists `.gitignore` as modified and the five HTML files as deleted from the index (`D `), and **no** `??` line for any of them.

### 1.3 The Quarto conversion — GATE

For each `<name>` in `methodology`, `paper_simulations`, `treatment_effect_definitions`, `biomarker_effects`, `causal_effects_brief_review`: `git mv vignettes/articles/<name>.Rmd vignettes/articles/<name>.qmd`, then `cp dev/tasks/articles_finish_2026-10-30/qmd/<name>.qmd vignettes/articles/<name>.qmd`. *GATE:* each placed file's md5 equals the `qmd/` entry of `MANIFEST.md5`, and `ls vignettes/articles/*.Rmd` lists only `forestsearch-2024.Rmd`.

Note (v2): `git mv` stages the rename at once. The renames therefore belong to commit 2 of §1.9, and must not be in the index when commit 1 is made (§1.9 says how).

What the conversion is (for the report; nothing to do): the YAML header becomes the Quarto `format: html` form the package's other articles use (`embed-resources`, `toc`, `toc-depth: 3`, `number-sections`, `code-fold` where the R Markdown version had `code_folding: hide`, the vignette engine `quarto::html`), with title, subtitle, author, bibliography, `nocite` and `VignetteIndexEntry` unchanged; the body is unchanged except (a) the three `DiagrammeR::grViz()` diagrams (two in methodology, one in paper_simulations) become native Quarto `{dot}` cells with the same DOT source, (b) methodology's power-approximation equation gets the label `{#eq-power}` that its existing cross-reference `@eq-power` points to (R Markdown rendered it as "(eq-power?)"), and (c) treatment_effect_definitions' bookdown reference `Section \@ref(theoretical-relationships)` becomes a link. `qmd_grviz/` holds methodology and paper_simulations with (a) undone — the fallback of §1.6.

### 1.4 README and NEWS — GATE

1. `cp dev/tasks/articles_finish_2026-10-30/README.md README.md`. *GATE:* md5 `2d3ff9f92188b80092a786eb4993da97`.
2. From the repo root:
   ```
   python3 -I dev/tasks/apply_alignment_edits.py dev/tasks/articles_finish_2026-10-30/articles_finish_edits_2026-10-30.json . --check
   python3 -I dev/tasks/apply_alignment_edits.py dev/tasks/articles_finish_2026-10-30/articles_finish_edits_2026-10-30.json .
   ```
   (`dev/tasks/apply_alignment_edits.py` is the applier committed at `794863f3`, md5 `1056eed3b82f1d889500c662649f6210`; the edit set's file keys are repo-relative paths.) *GATE:* both end `RESULT: OK`; transitions `NEWS.md` `29aa8783→4b181535` and `vignettes/articles/survival-post-selection.qmd` `90b5128f→5b840137`. The NEWS edit replaces the 0.4.1 Documentation bullet on the methodology and simulation articles with three bullets (the aligned articles and their conversion; the Guo–He article; the README). The article edit adds, after DINA's floor of 0.90, that this is the floor in the document that produced Table S2 and that Section 6 of 2026a gives a hazard ratio of 1.0 (open finding 9 of the 10-28 handoff).

### 1.5 The README's example runs — GATE

Extract the R block that starts with `library(forestsearch)` from `README.md` to `dev/reports/articles_finish_2026-10-30/readme_example.R` (verbatim), append the three lines below, and run it with `Rscript` (cap 10 min) → `logs/readme_example.log`:
```r
cat("REGION", paste(fit$sg.harm, collapse = " & "), "\n")
cat(sprintf("CHECK corrected %.3f field_lower %.3f complement %.3f field_s %.3f\n",
            fit$mr_inference$debiased$est, fit$mr_inference$field$lower_1s,
            fit$mr_inference$complement$debiased$est, fit$mr_inference$field$complement$upper_1s_s))
```
*GATE:* exit 0; the `CHECK` line prints `corrected 1.276 field_lower 0.634 complement 0.645 field_s 0.804` (Table 1 of 2026b, as `vignette("forestsearch")` reproduces it); `print(fit)`, `summary(fit)` and `fs_declaration_calibration()` produce output without error. A mismatch is a stop of Part 1 (the README's call is not the vignette's).

### 1.6 Render — GATE

Record `git status --porcelain` before (`logs/status_before_render.txt`).

1. **Local copies (in place).** For each of the six articles changed by this part (the five `.qmd` and `survival-post-selection.qmd`), from the repo root: `quarto render vignettes/articles/<name>.qmd` under `timeout 20m` → `logs/render_local_<name>.log`, wall-clock recorded. It writes `vignettes/articles/<name>.html` beside the source — the copy Larry keeps; it is ignored.
2. **The site build.** `Rscript -e 'pkgdown::init_site()'` once, then for each of the six `Rscript -e 'pkgdown::build_article("articles/<name>", lazy = FALSE)'` under `timeout 20m` → `logs/render_pkgdown_<name>.log`. It writes `docs/articles/<name>.html` (git-ignored). This is the rendering the site will use (pkgdown renders Quarto articles with its own template and `minimal: true`).
3. **The diagrams — GATE, with a pre-set fallback.** On the pkgdown builds: `grep -c '>Step 1: Construct<' docs/articles/methodology.html` ≥ 1 and `grep -c '>Average Bias<' docs/articles/methodology.html` ≥ 1; `grep -c '>1. Create DGM<' docs/articles/paper_simulations.html` ≥ 1 (Graphviz draws each label line as an SVG `<text>` element, so these match only a drawn diagram). If an article fails: replace it with `dev/tasks/articles_finish_2026-10-30/qmd_grviz/<name>.qmd` (md5 from `MANIFEST.md5`), re-run 1–2 for it, and check instead that its pkgdown build carries the widget (`grep -c 'class="grViz html-widget' ≥ 1` per diagram) and that every `<script src="…">` in the page resolves to an existing file under `docs/`. If the fallback also fails → stop Part 1. Record which variant each article uses; the local render's diagrams are recorded, not gated.
4. **Deviation clause (conversion defects only).** If a render fails on a Quarto error in a construct the conversion produced (YAML key, cell option, cross-reference label), make the smallest edit that restores the intended rendering, record it verbatim (file, old line, new line, reason) under "Conversion deviations", and re-render that article once. Never change a number, a citation, a section reference or a sentence's meaning. A failure in text the conversion did not touch, or a second failure, stops Part 1.
5. After all renders, `git status --porcelain` again (`logs/status_after_render.txt`). It must equal the before-file except for ignored paths. If a new untracked path appeared (for example `vignettes/.quarto/` or `vignettes/_quarto.yaml` left behind), delete it, say so, and confirm the two files then agree.

### 1.7 Text gates — GATE (output only)

On each pkgdown build `docs/articles/<name>.html`: `results/<name>.txt` = `pandoc --wrap=none -t plain` of the page; `results/<name>_output.txt` = the same after every `<pre …>…</pre>` element is removed (dot-all; the strip script of `dev/reports/articles_alignment_2026-10-30/` may be reused — keep the one used in the record). Counts with `grep -o … | wc -l`, fixed strings:

| # | article | string | expected |
|---|---|---|---|
| 1 | methodology | `2026 note` | 13 |
| 2 | methodology | `(eq-power?)` · `?@eq-power` | 0 · 0 |
| 3 | methodology | `Equation` | ≥ 1 (the resolved cross-reference) |
| 4 | methodology | `These are the 2024 guidelines` | 1 |
| 5 | paper_simulations | `44 simulated trials under the alternative and 44 under the null` | 1 |
| 6 | paper_simulations | `Identifying Patient Subgroups` | 0 |
| 7 | treatment_effect_definitions | `@ref(` | 0 |
| 8 | treatment_effect_definitions | `The Post-Selection Target: the Standard Analysis of the Selected Region` | ≥ 1 |
| 9 | biomarker_effects | `generate_aft_dgm_flex` | ≥ 2 |
| 10 | causal_effects_brief_review | `What forestsearch Reports, Read Against This Review` · `not a claim of causal optimality` | ≥ 1 · 2 |
| 11 | survival-post-selection | `which this call reproduces` | 1 |
| 12 | all six | `NaN` · `?@` | 0 · 0 (a Quarto citation or cross-reference that failed to resolve prints as `?@key`) |

Headings are counted "≥ 1" because pkgdown builds its contents list in JavaScript. Record every count in `results/keyvalues.txt`. *GATE:* any count outside its value → stop Part 1.

### 1.8 Site configuration and URLs — GATE

1. `Rscript -e 'pkgdown::check_pkgdown()'` → `logs/check_pkgdown.txt`. *GATE:* `No problems found` (the five articles keep their names, so `_pkgdown.yml` is unchanged; `git diff --name-only -- _pkgdown.yml` prints nothing).
2. If `urlchecker` is installed: `Rscript -e 'urlchecker::url_check()'` → `logs/url_check.txt`. *GATE:* no failing URL in `README.md` or `NEWS.md` (failures elsewhere are recorded, not gated). If it is not installed, record that; the as-cran check of part B covers the incoming URL check.
3. If `spelling` is installed: `spelling::spell_check_files(c("README.md", "NEWS.md"))` → `logs/spelling.txt` (recorded, not gated).

### 1.9 Commits — GATE

**Resume (v2).** If `8cde6811` is HEAD, run `git reset --soft HEAD~1`, as Larry instructed (the branch is not pushed). If Larry's instruction has already been carried out, with commits 1–4 made as listed below, do not repeat them: check that they match the list, and go to §2. Before **each** commit, set the index to exactly that commit's paths: `git restore --staged .` (this removes nothing from the working tree), then `git add` (or `git rm --cached`, for the index removals) only that commit's paths, then check `git diff --cached --name-status`. For commit 2, re-stage the renames with `git add -A vignettes/articles/<name>.Rmd vignettes/articles/<name>.qmd` for each of the five, so that git records each as a rename (`R`). Never use `git commit -a`.

By explicit path, in this order, each with `git diff --cached --name-only` listing only its paths:
1. `.gitignore` and the five index removals — `chore: rendered article HTML stays on disk and is ignored by git`;
2. the five renames with their new content — `docs: the five aligned articles as Quarto documents`;
3. `README.md`, `NEWS.md`, `vignettes/articles/survival-post-selection.qmd` — `docs: README for the current method; NEWS for the articles; DINA floor note`;
4. `dev/reports/articles_finish_2026-10-30/` and `dev/reports/REPORT_articles_finish_2026-10-30.md`.

Post-conditions, printed in the report: `git diff --name-only 0bdc2036 HEAD | grep -vE '^(dev/|vignettes/articles/|\.gitignore$)'` prints exactly `NEWS.md` and `README.md`; `git ls-files 'vignettes/articles/*.html'` prints nothing; `ls vignettes/articles/*.html` lists the five HTML of §1.1 plus `biomarker_effects.html`, `causal_effects_brief_review.html`, `treatment_effect_definitions.html` (new local renders).

**Report — `dev/reports/REPORT_articles_finish_2026-10-30.md`:** provenance; §1.1–1.4 outcomes; the README example's output (region, `CHECK` line); per article: local and pkgdown wall-clock, diagram variant used; the text-gate table; `check_pkgdown()`, URL and spelling outputs; conversion deviations (or "none"); OPEN ITEMS (findings only — e.g. anything in the pkgdown builds that renders differently from the R Markdown version: section numbers, figure captions, code folding); a ten-line verdict.

---

## Part 2 — ACTG 175 article, stage 0: provenance, reproduction, timing (read-only)

Sibling repos are **read only**: nothing edited, built, rendered or committed in `fs-glms-interpretable` or `fs-post-selection`. Record and transplants go to `dev/reports/actg175_article_stage0_2026-10-30/`. The ACTG 175 data come from `speff2trial::ACTG175` (Suggests); if `speff2trial` is not installed, stop Part 2 and say so — no installs.

**The comparison targets.**

*2026b, Table S16, the two ACTG 175 columns (forest search; printed precision 3 dp for the OR column, 1 dp for the MD column):*

| row | binary, OR | continuous, MD (CD4 cells) |
|---|---|---|
| Whole trial N (events); estimate, 95% CI | 1,083 (430); 0.665 (0.520, 0.850) | 1,083; −27.6 (−43.1, −12.1) |
| Region | weight > 86 kg and CD4 > 380 | CD4 > 507 and male |
| N (events) | 72 (37) | 79 |
| Unadjusted, 95% CI | 3.584 (1.315, 9.774) | 71.1 (−5.1, 147.3) |
| Corrected (two-term), IJ interval | 1.773 (0.604, 5.206) | 22.1 (−77.5, 121.6) |
| Field, lower bound | 0.686 | −38.3 |
| Complement N (events) | 1,011 (393) | 1,004 |
| Complement unadjusted, 95% CI | 0.587 (0.454, 0.759) | −36.2 (−51.4, −20.9) |
| Complement corrected, IJ interval | 0.620 (0.372, 1.033) | −32.7 (−62.2, −3.2) |
| Field-s, upper bound | 0.786 | −17.6 |
| Bonferroni pair | 0.597 / 0.817 | −51.2 / −15.8 |
| Selection bias; re-selection family | 0.704; 2,794 | 47.9; 4,935 |

*2026b, Table S17, the two ACTG 175 columns (α = 0.10):* protected level OR 0.75 / MD −10 cells; T of the region / largest T 2.494 / 2.494 and 1.701 / 2.043; calibrated cutoff 2.889 and 2.952; implied level 0.996 and 0.997; diagnostic at the default screen 0.838 and 0.909; cutoff at the claim threshold 3.611 and 3.625; diagnostic there 0.989 and 0.997; the calibrated rule admits nothing in both. Grids: OR 0.70–0.85; MD 5–15 cells.

*2026a, Table S3 (binary, OR; 2 dp):*

| Method | SG | n (%) | Naive | Full bootstrap | Multiplier |
|---|---|---|---|---|---|
| FS | Ĥ | 72 (7) | 3.58 (1.31, 9.77) | 2.43 (0.80, 7.33) | 1.79 (0.63, 5.08) |
| FS | Ĥᶜ | 1011 (93) | 0.59 (0.45, 0.76) | 0.61 (0.41, 0.89) | 0.62 (0.37, 1.03) |
| DINA | Ĥ | 98 (9) | 1.94 (0.84, 4.47) | 0.85 (0.23, 3.11) | 0.95 (0.29, 3.10) |
| DINA | Ĥᶜ | 985 (91) | 0.60 (0.46, 0.78) | 0.63 (0.40, 1.02) | 0.63 (0.38, 1.06) |
| GRF | Ĥ | 95 (9) | 2.25 (0.97, 5.20) | 1.09 (0.40, 3.01) | 1.10 (0.43, 2.80) |
| GRF | Ĥᶜ | 988 (91) | 0.58 (0.45, 0.76) | 0.63 (0.38, 1.04) | 0.62 (0.37, 1.04) |

Regions: FS `wtkg > 86 & cd40 > 380`; DINA `preanti >= 849.4 & cd40 >= 338`; GRF `wtkg > 84.37 & cd40 > 368.2`.

### 2.1 Inventory A — where 2026b's ACTG 175 numbers came from

1. Find `fs-glms-interpretable` (`~/Documents/GitHub/`; record HEAD and `git status --short`). Locate the manuscript source of Tables S16 and S17: grep its `.qmd`/`.tex`/`.md` for `The five analyses side by side`, `The three declarations read against`, `3.584`, `1.773`, `0.686`, `71.1`, `22.1`. Record the fragment, the "applications extract" it reads, and how that extract is built. Copy fragment and extract builder to `manuscript_source_2026b/`.
2. The producing documents — expected `quarto/actg175/analysis_actg175_binary_mr.qmd` and `analysis_actg175_continuous_mr.qmd` with their `<doc>_payload.rds` (the 10-27 stage 0 listed them: rendered 2026-09-23 under 0.3.5.9000). Find the document behind Table S17's declaration calibration (it may be these, or a separate one). *GATE:* the payloads (or the extract) carry Table S16's region, N, unadjusted, corrected, field and field-s cells for both ACTG 175 columns at printed precision (the whole-trial row may be computed elsewhere — record where); record render date, `forestsearch_version`, commit.
3. Transplants into `calls_2026b/`, verbatim with a header (source path, chunk names and fence lines, sibling HEAD): data preparation; the forest-search call for each endpoint; the declaration-calibration calls; the code that assembles the Table S16/S17 cells.

### 2.2 Inventory B — where 2026a's Table S3 came from

1. Find `fs-post-selection` (record HEAD and status). Locate the Table S3 fragment (grep `De-biased subgroup odds ratios for the ACTG175`, `3.58`, `2.43`, `1.79`), the payload it reads, and its manifest entry (`manuscript/payload_manifest.qmd`). Copy them to `manuscript_source_2026a/`.
2. Follow the manifest to the producing document and payload — expected in forestsearch `quarto/applications/actg175/` (an `analysis_actg175_binary_multimethod*.qmd`; record the commit as of the payload's build time, as the 10-28 Table S2 stage 0 did). Record whether the B = 1,000 bootstrap ran in the same render. *GATE:* the document and a payload carrying Table S3's FS row (3.58 / 2.43 / 1.79 at 2 dp, n 72) are located.
3. Transplants into `calls_2026a/`: data preparation; the three identifier calls; the bootstrap calls with their arguments; the table assembly.
4. **(v2) The rendered text.** v1 found that the payload was rendered from the document as committed in `43b051b6`, not from the commit the manifest names. For every chunk transplanted in step 3, diff that commit's version against the transplant (`git show 43b051b6:<path>`) → `logs/s3_transplant_vs_43b051b6.diff`. If any chunk differs, re-transplant it from `43b051b6`, mark the header to say so, and re-run M-S3-FS on the new transplant before Gate C is evaluated. Record which.

### 2.3 Settings and package facts

1. `settings.md`: one table per analysis listing every argument that differs from the GBSG Table 1 calls (`dev/reports/vignettes_stage0_2026-10-27/table1_calls/`) and from `?forestsearch` Examples 2 and 3 (`R/forestsearch_main.R` ≈ lines 1473–1560), with values side by side.
2. `package_facts.md`, quoting file and line: the element paths of every Table S16 cell on the GLM paths (OR and MD) in `fit$mr_inference`; whether `print()`/`summary()` show them; `fs_declaration_calibration()` on the MD path (a negative `c0` on the difference scale); `forestsearch_bootstrap_dofuture()` and `mr_estimates_table()` on a binary fit (`est.scale`, column names).

### 2.4 The ladder — compute (approved)

Installed package, single process (`parallel_args = list(plan = "sequential")`), the transplanted calls verbatim plus `quiet = TRUE`, `details = FALSE`, `plot.sg = FALSE` (`plot.grf = FALSE` where the formal exists). Driver `timeout 1200 /usr/bin/time -v Rscript ladder.R <rung>`, per-rung logs, `ladder.csv`, per-rung `.rds` of args and fit.

| rung | what | records |
|---|---|---|
| M-bin | 2026b binary forest search as its document ran it | every Table S16 OR cell at full precision; wall-clock; peak RSS |
| M-cont | 2026b continuous forest search | every Table S16 MD cell; wall-clock; peak RSS |
| D-bin, D-cont | the Table S17 declaration calibrations on M-bin / M-cont | every Table S17 cell; wall-clock |
| M-S3-FS, M-S3-DINA, M-S3-GRF | 2026a's three identifier calls | region, n, naive and MR estimate with IJ interval, both regions; wall-clock |
| B-S3-FS, -DINA, -GRF | the document's bootstrap call on each M-S3 fit with `nb_boots = 5` (its seed) | per-replicate wall-clock (total / 5 and the first separately); whether it runs for DINA and GRF |

**GATE A (M-bin) and GATE B (M-cont):** every Table S16 cell of that column reproduces at printed precision — the region by membership and N (label order may differ; record it), and the whole-trial row from the transplanted data preparation with the analysis's own working model (logistic regression for the OR, linear regression for the MD). **GATE C (M-S3-FS), restated in v2:** the region and n 72, and both naive rows (3.58 (1.31, 9.77) and 0.59 (0.45, 0.76)), at 2 dp exactly. The multiplier rows — 1.79 (0.63, 5.08) and 0.62 (0.37, 1.03) — must match at 2 dp in the point estimates, and within 0.01 in each interval limit. Table S3 was computed under forestsearch 0.2.0, and the 10-28 stage 0 found the forest-search multiplier cells of Table S2 differing from 0.2.0 in the third decimal (open finding 5 of the 10-28 handoff), so a one-unit difference in the second decimal of a limit is within that drift. Record every multiplier cell at full precision under 0.4.1, beside the payload's value and their difference. The same rule applies to DINA and GRF, which are recorded, not gated. A failed gate stops Part 2 with the full `sg.harm`, the counters and the settings in force. Recorded, not gated: Table S17 (D-bin, D-cont); DINA and GRF against Table S3; every quantity against the located payloads at full precision (`all.equal`, tolerance 1e-8, and the maximum absolute difference).

**Budget arithmetic in the report:** each fit's wall-clock; per-replicate bootstrap cost per identifier, and from it the B each identifier can afford in 60 s, 120 s and 300 s of single-core time. No recommendation — part B sets B on these numbers.

### 2.5 Report and commit

`dev/reports/REPORT_actg175_article_stage0_2026-10-30.md`: provenance; §2.1–2.2 (fragments, extracts, manifests, producing documents and payloads, transplant list); §2.3 facts; the ladder table with both comparisons and the budget arithmetic; deviations verbatim; OPEN ITEMS (findings only). Commit the record directory and the report by path, alone: `dev: report — ACTG 175 article stage 0 (2026-10-30)`. Copy the report, `ladder.csv`, `settings.md`, `package_facts.md`, `calls_2026b/`, `calls_2026a/`, `manuscript_source_2026b/`, `manuscript_source_2026a/` and the producing documents' `.qmd` to `~/Downloads/actg175_article_stage0_2026-10-30/`; `ls -la` into the report's last section.

---

## 3. Closing — GATE

**Resume rule (v2).** Copy this document to `dev/tasks/cc_task_articles_finish_2026-10-30_v2.md` and commit it alone, before anything else is committed. Then run Part 1 from §1.9 only (the resume there), and Part 2 from §2.2.4: the diff, M-S3-FS again only if a chunk differed, Gate C as restated, then M-S3-DINA, M-S3-GRF and the three bootstrap rungs, and the budget arithmetic. The v1 rungs that completed (M-bin, M-cont, D-bin, D-cont) and their records stand as recorded. Their report sections are kept, and the report gains a "v2" section. The closing block v1 appended to the v1 document, if any, is left uncommitted: remove it with `git checkout -- dev/tasks/cc_task_articles_finish_2026-10-30.md`, and put the closing block on this document instead.

After both parts: `git log --oneline 0bdc2036..HEAD` (the v1 task commit, the v2 task commit, Part 1's four, Part 2's one — fewer if a part stopped); `git status --porcelain | grep -v '^?? '` prints nothing; `git status --porcelain | grep '^?? '` equals the §0 list. Append a closing block to this document under `dev/tasks/` with those outputs and each part's outcome, and commit it alone. Copy `REPORT_articles_finish_2026-10-30.md`, its `results/` and the six local HTML renders to `~/Downloads/articles_finish_2026-10-30/`. No push, no fetch, no merge.

## 4. Out of scope

The ACTG 175 article itself, `cran-comments.md`, the `R CMD check --as-cran` re-check and anything in `R/`, `man/`, `vignettes/forestsearch.qmd` (all part B, or Larry's call); the Guo–He article; `forestsearch-2024.Rmd` (stays R Markdown, as decided 2026-10-27); `extreme_subgroups.qmd`; merging.
