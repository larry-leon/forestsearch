# CC task — vignette rendering polish (2026-10-28)

**Repo:** `larry-leon/forestsearch`, branch `feature/vignettes-0.4.1`. **Kind:** documentation — `vignettes/forestsearch.qmd` (v5: layout and strings only) and two YAML lines in `vignettes/survival-post-selection.qmd`. **No edit to `R/`, `man/`, `DESCRIPTION`, `tests/`, `_pkgdown.yml`, `NEWS.md`, the articles.**
**Inputs from `~/Downloads`:** this file and `forestsearch_v5.qmd` (md5 `64f4a604da9b761520e1cb4fb958dde7`, 648 lines).
**Compute:** two renders (74 s and 37 s measured) and one identity run (≈ 75 s). Caps 20 min each.
**Unattended.** Gates stop the task. No `git fetch` / `pull` / `push`; explicit `git add` paths; never stage a rendered `.html`; no `R CMD check`, test suite, `pkgdown::build_site()` or `devtools::install()`.

What v5 changes, for the record: YAML gains `code-overflow: wrap` and `fontsize: 15px`; a raw-HTML `<style>` block after the setup chunk sets tables to content width; one introduction bullet is rewrapped so no source line begins with a number followed by `)` (Pandoc read "2023)" as a numbered-list marker and broke the bullet); long code lines are rewrapped; the declaration table has short column names with the explanation in a caption and its row label `HR 1.00 (claim threshold)`; two `cat()` lines are split. No computation changes.

## 0. First action — GATE

1. `git status --short`; never stage pre-existing untracked paths. *GATE:* branch `feature/vignettes-0.4.1`; `git merge-base --is-ancestor c1fb38d3 HEAD`. Record HEAD.
2. Copy this file to `dev/tasks/cc_task_vignette_render_2026-10-28.md`; commit it alone: **`dev: task — vignette rendering polish (2026-10-28)`**.
3. *GATE:* `~/Downloads/forestsearch_v5.qmd` has md5 `64f4a604da9b761520e1cb4fb958dde7`. Installed-build check as before: `git diff 936b6734 HEAD -- R/ DESCRIPTION NAMESPACE inst/` is empty (no NEWS identity).

## 1. `vignettes/forestsearch.qmd` ← v5 — GATES

1. `cp ~/Downloads/forestsearch_v5.qmd vignettes/forestsearch.qmd`.
2. *GATE (markdown):* `grep -nE '^\s*[0-9]+[.)]( |$)' vignettes/forestsearch.qmd` returns nothing.
3. Render under `/usr/bin/time -v`, cap 20 min: `quarto render vignettes/forestsearch.qmd`. *GATE:* exit 0; HTML greps of `Error in`, `Warning`, `NaN`, `NA (` all 0. Record wall-clock, size.
4. *GATE (identity, re-run because the chunks changed):* the task-1 gate script `dev/reports/vignette_primary_2026-10-27/gate_identity.R`, pointed at the new file — purl, source in a fresh `Rscript`, the 54 comparisons against `dev/reports/vignettes_stage0_2026-10-27/fits/fit_{S3,D3,G3}.rds` and `ladder_quantities.csv`, all equal; the declaration values at 3 dp equal Table S17's GBSG column (2.726 / 0.994 / 0.667 / 2.438 / 2.564 / 3.503 / 0.974 / all admits 0). Output to `dev/reports/vignette_render_2026-10-28/gate_identity_output.txt`.
5. *GATE (rendered content):* the Table 1 block of the new HTML (pandoc plain text, "Table 1" heading to "Declaration at a protected level") is identical to the previous render's except for whitespace (`diff -w`); the rendered introduction contains no `<ol` element (no numbered list) — `grep -c '<ol' vignettes/forestsearch.html` is 0.
6. Recorded, not gated: the rendered declaration table's header row text; whether any `<pre>` block still overflows (quarto's `code-overflow: wrap` sets `white-space: pre-wrap` on `pre > code` — confirm by grepping the HTML's inline CSS for `pre-wrap`).
7. Commit `vignettes/forestsearch.qmd`: **`docs: vignette rendering — content-width tables, wrapped code, 15px base, short declaration headers; one bullet rewrapped`**.

## 2. `vignettes/survival-post-selection.qmd` — two YAML lines

Replace exactly
```
    toc-depth: 3
    code-fold: false
vignette: >
```
with
```
    toc-depth: 3
    code-fold: false
    code-overflow: wrap
    fontsize: 15px
vignette: >
```
Render, cap 20 min; *GATE:* exit 0, no `Error in` / `Warning` / `NaN` in the HTML; the printed quantities equal the previous render's (pandoc plain text `diff -w` of the `#>` output lines). Commit: **`docs: post-selection vignette — code wrapping and 15px base, matching the primary vignette`**.

## 3. Report and close-out

Append a section `## Rendering polish (2026-10-28)` to `dev/reports/REPORT_vignettes_align_2026-10-28.md` (gates, wall-clocks, sizes, the recorded items, deviations verbatim); commit it with `dev/reports/vignette_render_2026-10-28/`: **`dev: report — vignette rendering polish (2026-10-28)`**. Copy the two rendered HTMLs and the report to `~/Downloads/vignette_render_2026-10-28/`. No push.
