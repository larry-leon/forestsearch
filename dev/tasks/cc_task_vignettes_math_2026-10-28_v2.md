# CC task v2 — mathematical typesetting in both CRAN vignettes; post-selection presentation fixes; bootstrap on two workers (2026-10-28)

**Supersedes** `cc_task_vignettes_math_2026-10-28.md` (v1, never executed; do not execute it). Changes from v1: the post-selection input is v5 (bootstrap on `min(2, availableCores())` workers through a visible `n_workers` knob; one sentence on the fits' plan), the §2 purl gate admits the `bootstrap` chunk, the previous run's bootstrap values are a recorded comparison rather than a gate, and a sequential-versus-two-workers check is recorded.

**Repo:** `larry-leon/forestsearch`, branch `feature/vignettes-0.4.1`. **Kind:** documentation — `vignettes/forestsearch.qmd` ← v6 (prose only: Unicode math replaced by LaTeX math, two display equations), `vignettes/survival-post-selection.qmd` ← v5 (the same, plus: the side-by-side bootstrap column as a point estimate only; the cost sentence saying what each number measures; the membership-recovery sentence made live; the bootstrap run on two workers behind an `n_workers` knob). **No edit to `R/`, `man/`, `DESCRIPTION`, `NEWS.md`, `tests/`, `_pkgdown.yml`, the articles.**
**Inputs from `~/Downloads`:** this file, `forestsearch_v6.qmd` (md5 `8ed0d81d373a4427a1a74a0b0dc9817f`, 662 lines), `survival-post-selection_v5.qmd` (md5 `8c56077e1fdbd0b1c13fbf4c5d082e15`, 494 lines).
**Compute:** one render of each vignette (74 s; and 2:51 measured sequential, expected less on two workers) and one FS bootstrap at B = 25 sequential for the comparison (≈ 55 s). Caps 20 min per render, 10 min for the comparison. No identity run of the fits: the fit chunks are byte-identical to the committed files and the gates prove it.
**Unattended.** Gates stop the task — a failed gate is a stop even when the reason is clear; CC reports with a proposed fix. No `git fetch` / `pull` / `push`; explicit `git add` paths; never stage a rendered `.html`; no `R CMD check`, test suite, site build, install or `load_all()`.

Math rendering: Quarto's default `html-math-method` (MathJax from its CDN at viewing time; nothing downloaded at build) — the arrangement of `rmarkdown::html_vignette` and the package's articles. Nothing is set in the YAML for it.

---

## 0. First action — GATE

1. `git status --short`; never stage pre-existing untracked paths. *GATE:* branch `feature/vignettes-0.4.1`; `git merge-base --is-ancestor e7f4fa6d HEAD`. Record HEAD.
2. Copy this file to `dev/tasks/cc_task_vignettes_math_2026-10-28_v2.md`; commit it alone: **`dev: task — vignette math v2 (2026-10-28)`**.
3. *GATE (inputs and base):* the two inputs have the md5s above; `md5sum vignettes/forestsearch.qmd` at HEAD is `234322275f3df0ef4fee47c5509b81e5` (the file v6 was built from — if it differs, STOP and report the diff against the input rather than placing it); `md5sum vignettes/survival-post-selection.qmd` at HEAD is `2f955e6926d231e4eae1bd191eeb2fa9`.
4. Installed-build check: `git diff 936b6734 HEAD -- R/ DESCRIPTION NAMESPACE inst/` empty. Record `Rscript -e 'cat(future::availableCores())'` on this machine (the vignette's `n_workers` will resolve to 2 here).

---

## 1. `vignettes/forestsearch.qmd` ← v6 — commit 2

1. `cp ~/Downloads/forestsearch_v6.qmd vignettes/forestsearch.qmd`.
2. *GATE (prose only):* `knitr::purl()` of HEAD's file and of the new file are byte-identical.
3. *GATE (markdown):* `grep -nE '^\s*[0-9]+[.)]( |$)'` returns nothing; outside code chunks the file contains none of `β`, `Ĥ`, `Φ`, `κ̂` (strip the chunks with `awk` first; code and output strings never had them).
4. Render under `/usr/bin/time -v`, cap 20 min. *GATE:* exit 0; HTML greps of `Error in`, `Warning`, `NaN`, `NA (` all 0; `grep -c '<ol'` 0; `grep -ci 'mathjax'` ≥ 1; `grep -c 'class="math inline"'` ≥ 20 and `grep -c 'class="math display"'` = 2; the HTML contains `1.276`, `1.076`, `1.111`. Record wall-clock, size.
5. Recorded, not gated: `pandoc -t plain` of the HTML around the first display equation and around "T(g)", pasted into the report (pandoc shows the TeX source; the browser renders it).
6. Commit `vignettes/forestsearch.qmd`: **`docs: vignette — LaTeX math for the thresholds, the closed-form rate, the selection rule and the declaration statistic`**.

---

## 2. `vignettes/survival-post-selection.qmd` ← v5 — commit 3

1. `cp ~/Downloads/survival-post-selection_v5.qmd vignettes/survival-post-selection.qmd`.
2. *GATE (code change confined):* `knitr::purl()` of HEAD's file and of the new file differ only inside the `three-columns` chunk (the `bootstrap =` element, the `boot_label` assignment, the `caption =` argument) and the `bootstrap` chunk (the two `n_workers` / `boot_plan` lines and `parallel_args = boot_plan`). `diff` the purled scripts and assert every differing line lies in one of those two chunks. The `fit-fs`, `fit-dina` and `fit-grf` chunks are byte-identical, so the identity result of `e7f4fa6d` stands for the fits.
3. *GATE (markdown):* as §1.3.
4. Render under `/usr/bin/time -v`, cap 20 min. *GATE:* exit 0; HTML greps of `Error in`, `NaN`, `Warning` 0; `grep -c 'gt_table'` ≥ 2; `grep -ci mathjax` ≥ 1; the six published multiplier strings each occur at least twice (`1.44 (0.62, 3.36)`, `0.64 (0.40, 1.03)`, `1.11 (0.45, 2.75)`, `0.66 (0.40, 1.09)`, `1.05 (0.47, 2.34)`, `0.66 (0.41, 1.06)`); in the pandoc plain text of the side-by-side table the "Full bootstrap estimate" column holds bare numbers (no parenthesis). Record wall-clock, size, and the printed `boot_secs` beside the previous run's (53 / 30.5 / 36 s sequential).
5. **Recorded, not gated — plan invariance.** In a fresh `Rscript` (cap 10 min): purl the vignette, source it up to and including the `fit-fs` chunk only (stop before the bootstrap — e.g. evaluate the purled script's expressions until `fit_fs` exists), then run `forestsearch_bootstrap_dofuture(fs.est = fit_fs, nb_boots = 25L, seed = 8316951L, mr_in_replicates = FALSE, parallel_args = list(plan = "sequential"))` and the same with `list(plan = "multisession", workers = 2L)`. Compare `H_estimates$H2`, `Hc_estimates$H2` and `results` (the per-replicate table) between the two: identical or not, and the wall-clocks. The report states which; the vignette's text claims nothing about it.
6. Recorded: the rendered recovery sentence (pandoc plain text, the paragraph beginning "Read together") and the cost table; the side-by-side bootstrap estimates beside the previous run's (FS 1.88, DINA 0.83, GRF 1.29).
7. Commit `vignettes/survival-post-selection.qmd`: **`docs: post-selection vignette — LaTeX math; bootstrap on two workers behind an n_workers knob; bootstrap column as point estimate at small B; cost and recovery sentences made exact`**.

---

## 3. Report and close-out — commit 4

Append `## Math, presentation and two-worker bootstrap pass (2026-10-28)` to `dev/reports/REPORT_postsel_vignette_2026-10-28.md` (gates, wall-clocks, sizes, the recorded excerpts and the plan-invariance result, deviations verbatim, OPEN ITEMS); commit it with the record additions under `dev/reports/postsel_vignette_2026-10-28/`: **`dev: report — vignette math pass (2026-10-28)`**. Copy both rendered HTMLs and the report to `~/Downloads/vignettes_math_2026-10-28/`. No push.

---

## 4. Out of scope

No edit to `R/`, `man/`, `DESCRIPTION`, `NEWS.md`, `tests/`, `_pkgdown.yml`, the articles; no change to any code chunk of `forestsearch.qmd` or to the fit chunks of `survival-post-selection.qmd`; no identity run of the fits; no check, test suite, site build or install; no push.
