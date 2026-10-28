# CC task — the post-selection vignette rewritten: the GBSG analysis of León & Anderson (2026a) with the full bootstrap illustrated (2026-10-28)

**Repo:** `larry-leon/forestsearch`, branch `feature/vignettes-0.4.1`. **Kind:** documentation — `vignettes/survival-post-selection.qmd` replaced by the draft written by chat; one attribution paragraph in `vignettes/forestsearch.qmd`; one NEWS bullet. **No edit to `R/`, `man/`, `DESCRIPTION`, `tests/`, `_pkgdown.yml`, the articles.** If a gate can only be passed by changing something under `R/`, STOP and report.
**Inputs from `~/Downloads`:** this file and `survival-post-selection_v2.qmd` (md5 `2f955e6926d231e4eae1bd191eeb2fa9`, 469 lines). Stage 0's record `dev/reports/postsel_stage0_2026-10-28/` supplies the comparison fits.
**Compute:** the vignette's three fits (stage 0: 31.8 + 1.7 + 11.3 s single-core), its bootstraps at B = 25 / 25 / 10 (stage 0's per-replicate costs project to about 44 + 11 + 37 s), one render of each CRAN vignette, and one identity run that repeats the three fits. Caps: 20 min per render, 20 min for the identity script. No campaign.
**Unattended.** Gates stop the task — a failed gate is a stop even when the reason is clear; CC reports with a proposed fix. No `git fetch` / `pull` / `push`; explicit `git add` paths; never stage a rendered `.html`; no `R CMD check`, test suite, `pkgdown::build_site()`, `devtools::install()` or `load_all()`.

What the new vignette is (Larry, 2026-10-28): the primary analysis of León & Anderson (2026a) — its Table S2 — run for the three identifiers at the manuscript's settings, with multiplier resampling in full (5,000 draws) and the full bootstrap implemented and run at a small B for illustration, `summarize_bootstrap_results()` and `mr_estimates_table()` shown, Table S2 printed as published for comparison; the p̂ and membership-recovery diagnostics and the reading-by-location kept. Exact reproduction of the bootstrap column is not sought.

---

## 0. First action — GATE

1. `git status --short`; never stage pre-existing untracked paths. *GATE:* branch `feature/vignettes-0.4.1`; `git merge-base --is-ancestor db5056f9 HEAD`. Record HEAD.
2. Copy this file to `dev/tasks/cc_task_postsel_vignette_2026-10-28.md`; commit it alone: **`dev: task — post-selection vignette (2026-10-28)`**.
3. *GATE:* `~/Downloads/survival-post-selection_v2.qmd` has md5 `2f955e6926d231e4eae1bd191eeb2fa9`; installed forestsearch `0.4.0.9000` and `git diff 936b6734 HEAD -- R/ DESCRIPTION NAMESPACE inst/` empty.

---

## 1. `vignettes/survival-post-selection.qmd` ← v2 — commit 2

1. `cp ~/Downloads/survival-post-selection_v2.qmd vignettes/survival-post-selection.qmd` (a replacement; the old file's content is in history at `5a1af95d`).
2. *GATE (markdown):* `grep -nE '^\s*[0-9]+[.)]( |$)' vignettes/survival-post-selection.qmd` returns nothing; `grep -nE 'dev/|readRDS|read\.csv|load\(|\.\./'` returns nothing.
3. **Render** under `/usr/bin/time -v`, cap 20 min: `quarto render vignettes/survival-post-selection.qmd`. *GATE:* exit 0; `grep -c` of `Error in`, `NaN` in the HTML is 0; `Warning` is 0 (the bootstrap chunk suppresses the package's small-B warning by `warning: false`; any other warning is a stop). Record wall-clock, peak RSS, HTML size. The gt tables (`summ_fs$table`, `mr_estimates_table(...)`) must appear in the HTML: `grep -c 'gt_table'` ≥ 2.
4. **Identity with stage 0 — GATE.** `dev/reports/postsel_vignette_2026-10-28/gate_identity.R`: `knitr::purl()` the vignette, source it in a fresh `Rscript` (this reruns the fits and the bootstraps), then compare with stage 0's saved M fits (`dev/reports/postsel_stage0_2026-10-28/fits/` or wherever `ladder.R` saved them — read the record's layout first):
   - for each of `fit_fs` / `fit_dina` / `fit_grf` against M-FS / M-DINA / M-GRF: `sg.harm`, `g$n_selected`, `g$complement$n`, `g$naive$est/lower/upper`, `g$complement$naive$est/lower/upper`, `g$debiased$est/lower/upper`, `g$complement$debiased$est/lower/upper`, `g$n_family`, `g$selection_bias` — all `identical()` or `abs(diff) == 0`. The vignette's calls drop the transplanted calls' inert and deprecated arguments and set `ci_method` explicitly to what each call resolved to; a difference means a dropped argument was not inert — STOP and report which.
   - for `fit_fs`, recorded not gated: `g$field$lower_1s`, `g$field$complement$upper_1s_s`, `g$field$joint_s$bonf_lower_H / bonf_upper_Hc` against M-FS's (the vignette adds `field_recovery = TRUE`; expected identical, say so either way), and the recovery values printed (`sens_H`, `ppv_H`, q10/q50/q90, `share_equal_1`).
   - *GATE (Table S2, naive and multiplier):* the `side_by_side` data frame's `naive` and `multiplier` cells equal `table_s2`'s for all six rows (string equality after the vignette's own 2-dp formatting). Its `bootstrap` cells are recorded, not gated.
   - Record the bootstrap objects' `H_estimates$H2` / `Hc_estimates$H2` for the three identifiers and `boot_secs`, and the `cost` table.
5. Commit `vignettes/survival-post-selection.qmd`: **`docs: post-selection vignette — the GBSG analysis of León & Anderson (2026a): three identifiers, multiplier resampling in full, the full bootstrap at a small B, mr_estimates_table()`**.

---

## 2. `vignettes/forestsearch.qmd` — one attribution paragraph — commit 3

Replace exactly
```
**Inference** (`mr_inference = TRUE`) follows León and Anderson (2026a): from
5,000 multiplier draws over every size-eligible candidate the identifier
proposes, the corrected (two-term) estimate of the effect in the identified
region Ĥ and in its complement Ĥᶜ; the nominal 95% one-sided *field* lower
bound on β(Ĥ); the *field-s* upper bound on β(Ĥᶜ), the studentized complement
bound; and the Bonferroni pair, a joint 95% claim on both. The field's own
```
with
```
**Inference** (`mr_inference = TRUE`): from 5,000 multiplier draws over every
size-eligible candidate the identifier proposes, the corrected (two-term)
estimate of the effect in the identified region Ĥ and in its complement Ĥᶜ
(León and Anderson 2026a); the nominal 95% one-sided *field* lower bound on
β(Ĥ), the *field-s* upper bound on β(Ĥᶜ), the studentized complement bound,
and the Bonferroni pair, a joint 95% claim on both (León and Anderson 2026b).
The field's own
```
*GATE:* `knitr::purl()` of the file before and after is byte-identical (prose only). Render, cap 20 min; *GATE:* exit 0; the HTML contains `1.276`, `1.076` and `1.111` (the three corrected estimates of Table 1) and `grep -c '<ol'` is 0. Commit: **`docs: vignette — field bounds and the two-subgroup claim attributed to 2026b, the correction to 2026a`**.

---

## 3. `NEWS.md` — commit 4

Replace exactly
```
* `vignette("survival-post-selection")` now runs at the manuscript's settings (consistency level 0.90,
  5,000 multiplier draws, the closed-form consistency screen) and cites both 2026 papers; the
  methodology and simulation articles on the package website carry a note placing them as the 2024
  workflow, with the deprecated arguments removed from their listings.
```
with
```
* `vignette("survival-post-selection")` now reproduces the GBSG analysis of León and Anderson (2026a):
  the three identifiers at the manuscript's settings, multiplier resampling in full and the full
  bootstrap at a small B for illustration, with `summarize_bootstrap_results()` and
  `mr_estimates_table()`. The methodology and simulation articles on the package website carry a
  note placing them as the 2024 workflow, with the deprecated arguments removed from their listings.
```
Commit `NEWS.md`: **`docs: NEWS — post-selection vignette rewritten`**.

---

## 4. Report and close-out — commit 5

`dev/reports/REPORT_postsel_vignette_2026-10-28.md`, record `dev/reports/postsel_vignette_2026-10-28/` (gate script and output, render logs): bullets — provenance; §1 render numbers; the identity comparison table (quantity × identifier × vignette × stage 0 × equal); the Table S2 naive/multiplier check; the bootstrap values and timings at B = 25 / 25 / 10 with the per-replicate cost beside stage 0's; §2's checks; deviations verbatim; OPEN ITEMS. Commit: **`dev: report — post-selection vignette (2026-10-28)`**. No push.

Copy to `~/Downloads/postsel_vignette_2026-10-28/`: the report, both rendered HTMLs, `vignettes/survival-post-selection.qmd` as committed. `ls -la` into the report's last section.

---

## 5. Out of scope

No edit to `R/`, `man/`, `DESCRIPTION`, `tests/`, `_pkgdown.yml`, the articles, or any other part of `vignettes/forestsearch.qmd`; no `R CMD check` / `devtools::check()`, test suite or site build; no change to the stage-0 record; no push.
