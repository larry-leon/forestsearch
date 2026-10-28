# CC task — the primary vignette: Table 1 of León & Anderson (2026b) end to end, the 2024 vignette moved to the articles, site navigation (2026-10-27)

**Repo:** `larry-leon/forestsearch` — run CC from this repo's root, on branch `feature/vignettes-0.4.1` (it exists; stage 0 is on it).
**Kind:** documentation. Touches `vignettes/` (one new file; one `git mv` into `vignettes/articles/` with a small edit; one two-line edit in `survival-post-selection.qmd`), `_pkgdown.yml`, `NEWS.md`. **No edit to `R/`, `man/`, `DESCRIPTION`, `tests/`, `.Rbuildignore`.** If a gate can only be passed by changing something under `R/`, STOP and report — that is Larry's decision, not this task's.
**Inputs from `~/Downloads`:** this file and `forestsearch.qmd` (the vignette, written by chat). Stage 0's record is in the repo: `dev/reports/vignettes_stage0_2026-10-27/` (`package_facts.md`, `ladder_quantities.csv`, `fits/fit_<rung>.rds`, `table1_calls/`).
**Compute:** the vignette's three fits (measured at stage 0: 75 s on one core) run twice — once in the render, once in the identity gate — plus one render of the moved article. Hard cap 20 min on each render and on the gate script. No campaign, nothing parallel.
**Unattended.** Gates stop the task; never ask, never work around. No `git fetch` / `pull` / `push`. Explicit paths on every `git add`; never stage pre-existing untracked files or any rendered `.html`. No `R CMD check`, no `devtools::check()`, no test suite, no `pkgdown::build_site()`, no `load_all()` — the installed build (0.4.0.9000, HEAD's tree) is what the vignette runs against.

Standing decisions this task executes (Larry, 2026-10-27): D1 the current `vignettes/forestsearch.Rmd` moves as-is to `vignettes/articles/forestsearch-2024.Rmd` with a pointer note, no `.Rmd` → `.qmd` conversion; D2 the new primary vignette is `.qmd`; D3 all three identifiers; D4 the vignette runs the analyses at the manuscript's settings outright (stage 0 measured 75 s and exact reproduction of Table 1 — no reduced draws, no `params` switch). D5–D7 are later tasks.

---

## 0. First action — branch, task doc, inputs — GATE

1. `git status --short`: record pre-existing untracked or modified paths; never stage them. *GATE:* current branch is `feature/vignettes-0.4.1` and `git merge-base --is-ancestor 30aee979 HEAD` succeeds (stage 0's report commit is in history). Record HEAD.
2. Copy this file from `~/Downloads` to `dev/tasks/cc_task_vignette_primary_2026-10-27.md` (the downloaded name may have its hyphens stripped; the committed name is this one). `git add` that path only; commit **`dev: task — primary vignette (2026-10-27)`**.
3. *GATE:* `~/Downloads/forestsearch.qmd` exists (or `forestsearch.qmd` with the same stem under a stripped name), and content-grep verification passes: the file contains the line `  %\VignetteIndexEntry{Interpretable subgroup identification with forestsearch}`, exactly three occurrences of `<- forestsearch(`, the string `declaration_c0 = c(0.70, 0.75, 0.80, 0.85)`, and the chunk labels `fit-fs`, `fit-dina`, `fit-grf`, `table1-helper`, `table1`, `declaration`. Record `wc -l` and `md5sum`.
4. Provenance to `dev/reports/vignette_primary_2026-10-27/provenance.txt`: hostname, HEAD, branch, R version, installed forestsearch version and `Built`, quarto CLI version, pkgdown version. *GATE:* installed forestsearch is `0.4.0.9000` and its installed `NEWS.md` is byte-identical to the repo's (as at stage 0).

---

## 1. Move the 2024 vignette to the articles — commit 2

1. `git mv vignettes/forestsearch.Rmd vignettes/articles/forestsearch-2024.Rmd`.
2. In the moved file, change only these lines:
   - YAML `title:` → `"Exploratory subgroup identification with forestsearch (the 2024 workflow)"` (keep the `subtitle:` line as it is).
   - `%\VignetteIndexEntry{Exploratory Subgroup Identification with ForestSearch}` → `%\VignetteIndexEntry{Exploratory subgroup identification with forestsearch (the 2024 workflow)}`.
   - Insert, as the first content after the YAML block (before the first heading or setup chunk), this blockquote verbatim (the file is an `.Rmd` rendered by rmarkdown under pkgdown, so a plain blockquote, not a Quarto callout):

   ```
   > **About this article (added October 2026).** This is the package's original vignette, written for
   > forestsearch as described in León, Jemielita, Guo, Marceau West and Anderson (2024), *Statistics in
   > Medicine*, <https://doi.org/10.1002/sim.10163>. It is kept as a reference for that publication: the full
   > argument listing, the bootstrap bias correction and the cross-validation workflow it documents are still
   > available in the package. The current methodology — the search as closed-form stability selection with
   > its declaration calibration (León and Anderson 2026b, arXiv:2610.28379) and the selection-adjusted
   > inference reported with `mr_inference = TRUE` (León and Anderson 2026a, arXiv:2609.38361) — is documented
   > in `vignette("forestsearch")` and `vignette("survival-post-selection")`. The code below runs as written;
   > its settings are those of the 2024 workflow and may differ from the current defaults and from the
   > manuscript analyses.
   ```

   No other change to the file.
3. *GATE:* `git diff --cached -M` after staging shows a rename with similarity ≥ 95% and a diff confined to the title line, the index-entry line and the inserted block.
4. Render check, cap 20 min: `Rscript -e 'rmarkdown::render("vignettes/articles/forestsearch-2024.Rmd")'` (the installed package; the article runs its bootstrap and cross-validation at the small settings it already has). *GATE:* the render succeeds. Leave the HTML on disk, untracked and unstaged (Larry's standing: rendered article HTML stays locally and is not deleted). Record wall-clock.
5. Commit by explicit path (`vignettes/articles/forestsearch-2024.Rmd`): **`docs: move the 2024 vignette to the pkgdown articles as forestsearch-2024, with a pointer note`**.

---

## 2. The new vignette — commit 3

### 2.1 Placement

`cp ~/Downloads/forestsearch.qmd vignettes/forestsearch.qmd`; `cmp` the two (identical). The file takes the `forestsearch` stem freed by §1, so `vignette("forestsearch")` and pkgdown's "Getting Started" both resolve to it; `VignetteBuilder: knitr, quarto` and `Suggests: quarto` are already in `DESCRIPTION` (stage 0, §2.4).

### 2.2 The only edits CC makes to the draft (each recorded in the report, old → new)

- **The `fs-candidates` chunk.** The draft indexes `fit_fs$grp.consistency$out_sg$result` by column names `N`, `hr`, `Pcons` and takes column 1 as the candidate definition. Read `names(readRDS("dev/reports/vignettes_stage0_2026-10-27/fits/fit_S0.rds")$grp.consistency$out_sg$result)` first and adjust the names and the column index to what exists; keep the table's four columns (candidate, N, HR, p_cons), the ordering by size, the band computation and the `stopifnot`. Then delete the two `# CC (task §2.2) …` comment lines. If the result has no per-candidate consistency rate or no size column, STOP and report.
- **A missing element path.** If any element the draft reads does not exist in the result objects (`package_facts.md` §2.3 is the recorded map), replace it with the recorded path. Nothing else.
- **A rendering option.** If a fit chunk leaks progress output into the HTML despite `quiet = TRUE` and `results: hide`, or a warning appears that is not a defect, add the minimal chunk option and record it.

**No edit to prose.** If a sentence would become false under what the fit returns — a different region, a different count, a sign — STOP and report the sentence and the value; chat rewrites, not CC.

### 2.3 Render — GATE

`( time quarto render vignettes/forestsearch.qmd ) > dev/reports/vignette_primary_2026-10-27/render_vignette.log 2>&1`, cap 20 min (`/usr/bin/time -v` around it if available, for peak memory). *GATE:* exit 0; `vignettes/forestsearch.html` exists; in the HTML, `grep -c` of `Error in`, `Warning`, `NaN`, `NA (` is 0 each; the HTML size is recorded (expect under 2 MB; it is self-contained by `embed-resources: true`). Do not stage the HTML or any `forestsearch_files/` directory.

### 2.4 Identity with stage 0 — GATE

Write `dev/reports/vignette_primary_2026-10-27/gate_identity.R` and run it in a fresh `Rscript` (cap 20 min), output to `gate_identity_output.txt`:

1. `knitr::purl("vignettes/forestsearch.qmd", output = <tempfile>)`; `source()` it (this re-runs the three fits, ≈ 75 s).
2. For each of `fit_fs` / `fit_dina` / `fit_grf` against stage 0's `fits/fit_S3.rds` / `fit_D3.rds` / `fit_G3.rds`: compare, with `identical()` or `abs(diff) == 0`, the region (`sg.harm`), `g$n_selected`, `g$complement$n`, `g$naive$est/lower/upper`, `g$complement$naive$est/lower/upper`, `g$debiased$est`, `g$complement$debiased$est`, `g$field$lower_1s`, `g$field$complement$upper_1s_s`, `g$field$joint_s$bonf_lower_H`, `g$field$joint_s$bonf_upper_Hc`, `g$n_family`, and events in Ĥ / Ĥᶜ computed from `df.est` membership. Also compare the same quantities with `ladder_quantities.csv` rows S3 / D3 / G3 at the CSV's 6 dp.
3. *GATE:* every comparison equal. A difference means a dropped or trimmed argument in the vignette's calls was not a no-op (the draft drops `max_n_confounders`, `fs.splits`, `cut_type`, `cont.cutoff`, `collapse_cuts`, `m1.threshold`, `minp`, `use_lasso/use_grf/use_dina` and the default-valued `mr_inference_args` elements of the transplanted calls): STOP and report which quantity and which fit differ — do not restore arguments yourself.
4. Record alongside, not gated: `funnel` equals 2,211 / 1,410 / 121 / 9 / 1,744; the candidates table has 9 rows, band `>= 2.03`, selected N 75, max N 82, max-β N 61; the three-way partition equals Table S15 (35/19, 40/22, 78/42, 7/4, 526/212).

### 2.5 The declaration section — GATE

From the same sourced session: `tb$kappa_hat[tb$c0 == 0.75]` = 2.726, implied level 0.994, `tb$fw_size[tb$c0 == 0.75]` = 0.667, `T_H` = 2.438, `T_max` = 2.564, `cal_c2$kappa_hat` = 3.503, `cal_c2$fw_size` = 0.974, every `admits` count 0 — all at 3 dp, the GBSG column of Supplementary Table S17. *GATE:* all equal; otherwise STOP and report the values (the vignette's prose is live but the section is published as reproducing the paper, and a difference must be understood first).

### 2.6 Self-containment — GATE

`grep -nE 'dev/|readRDS|read\.csv|load\(|\.\./' vignettes/forestsearch.qmd` returns nothing.

### 2.7 Commit

By explicit path (`vignettes/forestsearch.qmd`): **`docs: new getting-started vignette — Table 1 of León & Anderson (2026b) end to end`**.

---

## 3. Site navigation, NEWS, cross-reference — commit 4

### 3.1 `_pkgdown.yml` — exact edits

1. In `navbar: components: articles: menu:`, after the item
   ```
           - text: Extreme Subgroup Effects
             href: articles/extreme_subgroups.html
   ```
   append
   ```
           - text: "---"
           - text: The 2024 workflow (León et al. 2024)
             href: articles/forestsearch-2024.html
   ```
2. In `articles:`, after the `Simulation Studies` block, append
   ```
     - title: Reference
       desc: The package's original vignette, kept as a reference for León et al. (2024)
       contents:
         - "articles/forestsearch-2024"
   ```
3. `home: description:` → `Exploratory subgroup identification in clinical trials: interpretable harm-subgroup search with selection-adjusted inference for survival and generalized-linear outcomes.` (one line; the current text names only the bootstrap bias correction).
4. `Getting Started` (`contents: - forestsearch`) and the `Introduction to forestsearch` menu item stay as they are: they now resolve to the `.qmd`.
5. *GATE:* `Rscript -e 'pkgdown::check_pkgdown()'` passes (every vignette and article listed, no missing topic).

### 3.2 `NEWS.md`

Under `# forestsearch (development version)`, in a `## Documentation` subsection (append to it if one exists, otherwise create it after the existing subsections), add:

```
* New getting-started vignette, `vignette("forestsearch")`: the German Breast Cancer Study Group analysis of
  León and Anderson (2026b), Table 1, run end to end — the three identifiers (forest search, DINA, GRF) and
  the selection-adjusted estimates and bounds of León and Anderson (2026a) — at the manuscript's settings,
  with the declaration calibration at a protected level.
* The previous getting-started vignette, which documents the workflow of León et al. (2024), is now the
  pkgdown article "Exploratory subgroup identification with forestsearch (the 2024 workflow)" and is no
  longer built on CRAN.
```

### 3.3 `vignettes/survival-post-selection.qmd` — one cross-reference

Replace exactly
```
`fs$mr_inference`; `field_recovery = TRUE` adds the membership-agreement
diagnostics used below. Bootstrap and cross-validation are not run here; see
`vignette("forestsearch")`.
```
with
```
`fs$mr_inference`; `field_recovery = TRUE` adds the membership-agreement
diagnostics used below. Bootstrap and cross-validation are not run here; see
the article *Exploratory subgroup identification with forestsearch (the 2024
workflow)* on the package website.
```
Nothing else in that file changes in this task (its alignment to p\* 0.90 and 2026b is the next task).

### 3.4 Commit

By explicit paths (`_pkgdown.yml`, `NEWS.md`, `vignettes/survival-post-selection.qmd`): **`docs: pkgdown navigation, NEWS and cross-reference for the vignette change`**.

---

## 4. One read-only trace for the record (no edit)

Stage 0 found 121 candidates clearing the screening threshold (`filter_counts$n_passed_hr`) but 120 evaluated by the consistency stage (`grp.consistency$n_candidates_evaluated`). Using `fits/fit_S0.rds` and the source, say in one or two sentences which candidate is dropped and by which code path (file:line), cap 30 minutes; if not found in that time, say so. The vignette prints the 121 and the 9, not the 120.

---

## 5. Report — `dev/reports/REPORT_vignette_primary_2026-10-27.md`

Bullets, one item each: provenance; §1's diff summary and render wall-clock; §2.2's edits (old → new, verbatim); §2.3 render wall-clock, peak memory, HTML size; §2.4 gate output (the comparison table: quantity × fit × vignette value × stage-0 value × equal); §2.4's recorded checks; §2.5's values beside Table S17; §2.6; §3's `check_pkgdown()` output; §4's finding; deviations verbatim; OPEN ITEMS (findings only — nothing is fixed or proposed as a task).

---

## 6. Close-out — commit 5

1. Commit by explicit path: `dev/reports/vignette_primary_2026-10-27/` (provenance, `gate_identity.R`, `gate_identity_output.txt`, `render_vignette.log`, the article render log) and the report. **`dev: report — primary vignette (2026-10-27)`**. No push.
2. Copy to `~/Downloads/vignette_primary_2026-10-27/` (create it): the report; `vignettes/forestsearch.html` (the render); `vignettes/forestsearch.qmd` as committed; `vignettes/articles/forestsearch-2024.html` (the render); `_pkgdown.yml`; `NEWS.md`; `gate_identity_output.txt`. `ls -la` that directory into the report's last section. Larry attaches these to the chat; he reads the HTML before the next task.

---

## 7. Out of scope

- No edit to `R/`, `man/`, `DESCRIPTION`, `tests/`, `.Rbuildignore`, `README.md`; no edit to `methodology`, `paper_simulations` or the other articles (next task); no further edit to `survival-post-selection.qmd` beyond §3.3.
- No `R CMD check` / `devtools::check()`, no test suite, no `pkgdown::build_site()`, no `devtools::install()`.
- No change to any stage-0 record; no edit in fs-glms-interpretable; no push.
