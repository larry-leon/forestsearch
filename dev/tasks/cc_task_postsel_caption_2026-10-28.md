# CC task — post-selection vignette: the cost caption and the recovery sentence made exact for the two-worker bootstrap (2026-10-28)

Closes the two open items of the math pass (report section "Math, presentation and two-worker bootstrap pass", commit `25318e8c`): the cost-table caption still said "Single-core timings" after the bootstrap moved to two workers, and the recovery sentence rendered as "retain 0% or less".

**Repo:** `larry-leon/forestsearch`, branch `feature/vignettes-0.4.1`. **Kind:** documentation — `vignettes/survival-post-selection.qmd` ← v6. Changes from the committed v5: the caption of the cost table is built from `n_workers` (a `plan_label` line and a `sprintf()` caption in the `cost` chunk); the paragraph after it says the bootstrap column is wall-clock on that many workers; the recovery sentence chooses its wording by `rc$q10` (inline R, one line); "uses two, the most a CRAN build may use" → "uses at most two, the limit a CRAN check imposes" in two places. **No edit to `R/`, `man/`, `DESCRIPTION`, `NEWS.md`, `tests/`, `_pkgdown.yml`, the articles, or `vignettes/forestsearch.qmd`.**
**Input from `~/Downloads`:** this file and `survival-post-selection_v6.qmd` (md5 `660b2f87bd26edd0a5e1fdd9801fd09e`, 499 lines; line 447 is a 214-character inline-R expression, deliberately on one line).
**Compute:** one render (2:07 measured on two workers). Cap 20 min. No identity run: the fit, bootstrap and three-columns chunks are byte-identical to the committed file and the purl gate proves it.
**Unattended.** Gates stop the task — a failed gate is a stop even when the reason is clear; CC reports with a proposed fix. No `git fetch` / `pull` / `push`; explicit `git add` paths; never stage a rendered `.html`; no `R CMD check`, test suite, site build, install or `load_all()`.

HTML text gates run on `pandoc --wrap=none -t plain vignettes/survival-post-selection.html` (not on the raw HTML, whose MathJax polyfill contains `NaN`; `--wrap=none` so that a phrase is never split across lines), counting occurrences with `grep -o … | wc -l`.

---

## 0. First action — GATE

1. `git status --short`; never stage pre-existing untracked paths (the two article `.html` files). *GATE:* branch `feature/vignettes-0.4.1`; `git merge-base --is-ancestor 25318e8c HEAD`. Record HEAD.
2. Copy this file to `dev/tasks/cc_task_postsel_caption_2026-10-28.md`; commit it alone: **`dev: task — post-selection caption fix (2026-10-28)`**.
3. *GATE (input and base):* `~/Downloads/survival-post-selection_v6.qmd` has md5 `660b2f87bd26edd0a5e1fdd9801fd09e`; `md5sum vignettes/survival-post-selection.qmd` at HEAD is `8c56077e1fdbd0b1c13fbf4c5d082e15` (the committed v5 — if it differs, STOP and report the diff against the input rather than placing it).
4. Installed-build check: `git diff 936b6734 HEAD -- R/ DESCRIPTION NAMESPACE inst/` empty. Record `Rscript -e 'cat(future::availableCores())'` (the vignette's `n_workers` resolves to 2 here).

---

## 1. `vignettes/survival-post-selection.qmd` ← v6 — commit 2

1. `cp ~/Downloads/survival-post-selection_v6.qmd vignettes/survival-post-selection.qmd`.
2. *GATE (code change confined):* `knitr::purl()` of HEAD's file and of the new file differ only inside the `cost` chunk — one added line beginning `plan_label <-` and the `caption =` argument (now a `sprintf(paste(…), plan_label)` spanning four lines). `diff` the purled scripts and assert every differing line lies between the `## ----cost` marker and the next chunk marker. The `fit-fs`, `fit-dina`, `fit-grf`, `bootstrap` and `three-columns` chunks are byte-identical, so the identity result of `e7f4fa6d` stands for the fits and the bootstrap values of `ce9cfb38` are expected unchanged. (Inline R is not purled; the recovery-sentence change is prose to `purl()`.)
3. *GATE (markdown):* `grep -nE '^\s*[0-9]+[.)]( |$)' vignettes/survival-post-selection.qmd` returns nothing; `grep -c 'Single-core' vignettes/survival-post-selection.qmd` is 0.
4. Render under `/usr/bin/time -v`, cap 20 min: `quarto render vignettes/survival-post-selection.qmd`. Then `pandoc --wrap=none -t plain vignettes/survival-post-selection.html > /tmp/postsel_plain.txt`. *GATE:* exit 0; in the plain text, occurrence counts of `Error in`, `NaN`, `Warning` are 0 and `grep -cP '\bNA \('` is 0; `grep -c 'gt_table' vignettes/survival-post-selection.html` ≥ 2; `grep -ci mathjax vignettes/survival-post-selection.html` ≥ 1; the six published multiplier strings each occur at least twice in the plain text (`1.44 (0.62, 3.36)`, `0.64 (0.40, 1.03)`, `1.11 (0.45, 2.75)`, `0.66 (0.40, 1.09)`, `1.05 (0.47, 2.34)`, `0.66 (0.41, 1.06)`); in the raw HTML, `the bootstrap on 2 workers` occurs exactly once (the caption) and `Single-core` not at all; in the plain text, `0% or less` occurs not at all and `at least a tenth of the draws retain none of them` exactly once; the side-by-side table's "Full bootstrap estimate" column reads `1.88`, `0.83`, `1.29` for the three Ĥ rows (the values of the two previous renders — a difference is a stop). Record wall-clock, size, the printed `boot_secs` beside the previous run's (FS 32.6 / DINA 23.6 / GRF 21.8 s).
5. Recorded, not gated: the rendered cost table with its caption, and the paragraph beginning "Read together" (both from the plain text), pasted into the report.
6. Commit `vignettes/survival-post-selection.qmd`: **`docs: post-selection vignette — cost caption and recovery sentence made exact for the two-worker bootstrap`**.

---

## 2. Report and close-out — commit 3

Append `## Caption and recovery-sentence fix (2026-10-28)` to `dev/reports/REPORT_postsel_vignette_2026-10-28.md` (gates, wall-clock, size, `boot_secs`, the two recorded excerpts, deviations verbatim, OPEN ITEMS — the two items of the previous section are closed by this commit; say so). Commit it with the render log under `dev/reports/postsel_vignette_2026-10-28/` (a new file name; do not overwrite the previous logs): **`dev: report — post-selection caption fix (2026-10-28)`**. Copy the rendered HTML and the report to `~/Downloads/postsel_caption_2026-10-28/`. No push.

---

## 3. Out of scope

No edit to `R/`, `man/`, `DESCRIPTION`, `NEWS.md`, `tests/`, `_pkgdown.yml`, the articles or `vignettes/forestsearch.qmd`; no change to any chunk of `survival-post-selection.qmd` other than `cost`; no identity run; no check, test suite, site build or install; no push.
