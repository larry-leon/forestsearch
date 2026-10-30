# CC task — the Guo & He (S4) article: place, gate against the manuscript's record, render, list on the site — v2

**File:** `dev/tasks/cc_task_guohe_article_2026-10-30_v2.md` · **Issued:** 2026-10-30 by chat, commissioned by Larry (all recommendations of 2026-10-30 approved: default run 2 cells × 25 replicates at B = 200, MR at the manuscript's settings, truth at n_big = 2e5; complement and pair included; published S4 results displayed as literal data gated against the committed payload; file `vignettes/articles/guohe-nested-continuum.qmd`; version basis "0.4.0 or later"; merge when green).
**Repo:** `forestsearch` · **Branch:** `feature/articles-guohe` (continues from Stage 0; HEAD `f2ddb5c0`) · **Transport:** `~/Downloads` → `dev/tasks/` (this document, committed first and alone); `~/Downloads/guohe-nested-continuum.qmd` → `vignettes/articles/` (§2). Filenames may arrive with hyphens stripped; locate by pattern.

**v2 supersedes v1** (`dev/tasks/cc_task_guohe_article_2026-10-30.md`, committed `4ea458bb`). v1 stopped at its §5: its text gates counted the article's printed code — `code-fold: false` echoes every chunk, so each `sprintf()` caption and the `-Inf` in `argmax_candidate()` appeared once in the code and once in the output — and four counts came out one too high; gate 3 passed only because its caption is a template string. Larry's disposition (2026-10-30): count the rendered output only, no second render, the article unchanged. **This document is the complete, authoritative task.** Only §5 (output-only counting) and §7 (the resume rule) differ from v1; the render of v1's §4 stands.

**Inputs.** `guohe-nested-continuum.qmd`, md5 `18091f48356588306e784c8580e388d4`, 939 lines. The article is chat's; it carries its own data-generating mechanism, nested family, truth curve, bounds, Monte Carlo loop and summaries, and reads nothing from `quarto/` or `dev/`. The two package functions it calls, `fs_mr_inference()` and `guohe_algorithm3()`, are exported in CRAN 0.4.0 with code identical to `master` (`git diff v0.4.0 master -- R/fs_mr_inference.R` is empty; `R/guohe_algorithm3.R` differs in roxygen only), and the pkgdown workflow renders the site against `cran::forestsearch`, so the article is true on the site now.

**What this task produces.** The article in `vignettes/articles/`, its reference tables machine-checked against the manuscript's S4 payload, one local render, the `_pkgdown.yml` entries, and a report. Nothing else.

---

## ⚠ CATEGORY

**No package code is modified.** No edit under `R/`, `tests/`, `man/`, `inst/`, `DESCRIPTION`, `NAMESPACE`, `NEWS.md` (in the tarball — the article's NEWS line waits for the 0.4.1 release; OPEN ITEM), `vignettes/forestsearch.qmd`, or anything under `quarto/`. No `R CMD check`, no test suite, no install. Chat could not run R; the article's code was checked by reading, so §4's deviation clause exists for defects of the article's own code, never for its method.

**Writes:** this document under `dev/tasks/`; `vignettes/articles/guohe-nested-continuum.qmd`; two entries in `_pkgdown.yml`; the record directory `dev/reports/guohe_article_2026-10-30/` (gate script, results, logs); the report `dev/reports/REPORT_guohe_article_2026-10-30.md`. **The rendered HTML is never staged**: it stays untracked in `vignettes/articles/` (Larry decides tracking), and a copy goes to `~/Downloads`.

**Compute (its own line):** one render of the article at its default settings on `n_workers = min(2, availableCores()) = 2`: two truth curves at n_big = 2e5 (≈ 6 s each), one replicate at the manuscript's settings (Guo–He at B = 2000: ≈ 75 s from the Stage 0 rate), and 2 cells × 25 replicates at B = 200 (≈ 9.7 s each single core; ≈ 4 min on two workers). **≈ 6–8 min; cap 40 min** on the render. One extra render only if §4's deviation clause is invoked. No full-size run.

**Unattended.** Gates stop the task; never ask, never work around. Pre-existing untracked files (the Stage 0 folder is committed; the only expected untracked path after the render is the article's own HTML) are never staged.

---

## 1. Provenance and first commit — GATE

1. `git rev-parse --abbrev-ref HEAD` must be `feature/articles-guohe`; `git merge-base --is-ancestor f2ddb5c0 HEAD` must succeed. Record HEAD as **BASE_SHA** and `git status --porcelain` verbatim (expected empty). *GATE:* branch or ancestry differ → stop.
2. `md5sum ~/Downloads/guohe-nested-continuum.qmd` must print `18091f48356588306e784c8580e388d4` (locate by pattern if the name arrived altered). *GATE:* mismatch or absent → stop and report.
3. Copy this document to `dev/tasks/cc_task_guohe_article_2026-10-30.md`, `git add` by path, commit it alone. Report both filenames.
4. Record `Rscript -e 'cat(as.character(packageVersion("forestsearch")), future::availableCores(), "\n")'` (expect `0.4.1 128`), the quarto binary and version used in the previous article tasks (`/usr/lib/rstudio/resources/app/bin/quarto/bin/quarto --version`; if that path is absent, the `quarto` on `PATH`, and say which), and `pandoc --version | head -1`.

## 2. Place the article

`cp ~/Downloads/guohe-nested-continuum.qmd vignettes/articles/guohe-nested-continuum.qmd`; `md5sum` of the placed file must equal the input's. Do not commit yet — §3 may correct literals.

## 3. Reference-literal gate against the manuscript's S4 payload — GATE

The article's chunk `reference-data` holds the manuscript's Tables S20–S25 and S28 as literal data frames (`ref_lower`, `ref_select`, `ref_location`, `ref_upper`, `ref_joint`, `ref_16`) and Guo and He's published Table 7 (`ref_gh_published`). Their source of truth is the committed payload the manuscript's section was drafted from, `quarto/GuoHe/_payloads/guohe_supp_section/guohe_supp_section_payload.rds` (read-only; elements `bound_t7`, `ident_t7`, `comp_t7`, `joint_t7`, `pre16`, `meta`; written by `quarto/GuoHe/guohe_supp_section.qmd`, chunk `export-payload`), and `GH52_TABLE_7` in `quarto/GuoHe/guohe_sec52_sim.R` (≈ line 96) for the published Table 7 values.

Script `s1_reference_gate.R` in `dev/reports/guohe_article_2026-10-30/`, run from the repo root:

1. Extract from the placed article the text of chunk `summaries` lines defining `LOWER_METHODS` and `UPPER_METHODS`, and the whole text of chunk `reference-data` (the lines between the chunk's opening fence and its closing fence, dropping `#|` option lines); evaluate both in a fresh environment. *GATE:* evaluation error → stop (quote it).
2. Read the payload and source `GH52_TABLE_7` (evaluate only that assignment from the file's text, or `source()` the file with `sys.nframe()` guarding its executable block — it already does).
3. Compare every literal with its payload value, at the **printed precision**: a literal with `d` decimals matches when `abs(literal - payload) <= 0.5 * 10^-d + 1e-9`. The map (article column ↔ payload column; article method name ↔ payload method name):
   - `ref_lower` ↔ `bound_t7`, methods `Unadjusted`↔`Naive`, `Guo-He r = 1/12`↔`G&H r=1/12`, `MR (field)`↔`MR (field)`, `MR (IJ two-term)`↔`MR (IJ)`; `coverage`↔`cover` (3), `wilson_lo`↔`cover_lo` (3), `wilson_hi`↔`cover_hi` (3), `margin`↔`mean_margin` (3), `bias`↔`bias` (3), `lower_hr`↔`mean_lower_hr` (2).
   - `ref_select` ↔ `ident_t7`: `c_hat_median`↔`chat_med` (1), `c_hat_q25`↔`chat_q25` (1), `c_hat_q75`↔`chat_q75` (1), `c_hat_min`↔`chat_min` (1), `c_hat_max`↔`chat_max` (1), `specificity`↔`spec_mean` (3; both `NA` at β₂ = 0), `ppv`↔`ppv_mean` (3; same), `n_sel_mean`↔`n_sel_mean` (1), `n_sel_median`↔`n_sel_med` (0), `gamma_mean`↔`gamma_mean` (3), `dilution`↔`dilution` (3), `p_hat_median`↔`phat_med` (3), `p_hat_q25`↔`phat_q25` (3), `p_hat_q75`↔`phat_q75` (3), `share_at_30`↔`p_chat_30` (3), `M_eff_median`↔`M_eff_med` (2), `sd_marginal`↔`sd_marg` (3), `sd_error`↔`sd_err` (3).
   - `ref_location`: `theta_H`↔`exp(beta2)` (2), `hr_in_H_hat`↔`exp(ident_t7$gamma_mean)` (2).
   - `ref_upper` ↔ `comp_t7`, methods `Unadjusted`↔`naive`, `MR (IJ two-term)`↔`ij`, `MR (field-s)`↔`field_s`; `coverage`↔`cover` (3), `wilson_lo`↔`cover_lo` (3), `wilson_hi`↔`cover_hi` (3), `upper_hr`↔`mean_upper_hr` (2).
   - `ref_joint` ↔ `joint_t7`: `joint_coverage`↔`joint_s` (3), `wilson_lo`↔`joint_s_lo` (3), `wilson_hi`↔`joint_s_hi` (3), `joint_unscaled`↔`joint` (3), `unscaled_lo`↔`joint_lo` (3), `unscaled_hi`↔`joint_hi` (3), `corr_field_s`↔`corr_s_mean` (3), `corr_unscaled`↔`corr_mean` (3), `pair_lower_H_hr`↔`pair_lower_H_hr` (2), `pair_upper_Hc_hr`↔`pair_upper_Hc_s_hr` (2).
   - `ref_16` ↔ `pre16`, methods `G&H r=1/12` → `gh_*` and `MR (field)` → `field_*`: `*_coverage`↔`cover` (3), `*_margin`↔`margin` (3), `*_bias`↔`bias` (3); cells in payload order (`t35_beta2_00..05`, `t6_k02, k06, k10, k12`, `t7_beta2_00..05`) against the article's rows in the same order.
   - `ref_gh_published` ↔ `GH52_TABLE_7`: `gh_r_1_12`↔`r_1_12` (3), `unadjusted`↔`naive` (3).
4. Write every comparison to `results/reference_gate.csv` (table, row, column, literal, payload, abs_diff, pass) and the totals to `results/keyvalues.txt`.
5. **If every comparison passes:** the gate is green; go to §4.
   **If any fails:** the payload is authoritative. Replace the failing literal in `vignettes/articles/guohe-nested-continuum.qmd` with the payload value formatted at that column's printed precision (`sprintf("%.3f", …)`, `"%.2f"`, `"%.1f"`, `"%.0f"`), touching only that number, re-run the comparison, and list every correction in the report (table, row, column, article value → payload value). A literal that cannot be located for replacement, or a payload element or column that is absent, is a stop. *GATE after corrections:* every comparison passes.
6. Record the payload's `meta$nrep`, `meta$draws`, `meta$field_R_out`, `meta$field_R_in`, `meta$head_at_render`, `meta$forestsearch_version` in `keyvalues.txt`, and the article's stated manuscript settings beside them (2000 replicates, 5000 draws, 1000/500) — they must agree (record; a disagreement is a stop).

## 4. Render — GATE

From the repo root, with the quarto binary of §1.4:

```
<quarto> render vignettes/articles/guohe-nested-continuum.qmd 2>&1 | tee dev/reports/guohe_article_2026-10-30/logs/render.log
```

under `timeout 40m`. Record wall-clock (`/usr/bin/time -v` or `date` before and after) and `future::availableCores()` as the article sees it. *GATE:* non-zero exit → stop, with the log's last 40 lines in the report — **except** under the deviation clause below. The article's own `stopifnot()` (selection agreement across the three identifications) is a stop like any other; if it fires, report the replicate and do not relax it.

**Deviation clause (defects of the article's code only).** If the render fails on an R error inside one of the article's chunks that is a defect of the article's code — a misspelled element name, a wrong column name, a missing argument, a formatting call with the wrong number of arguments — make the smallest edit that restores the stated computation, re-render once, and list the edit verbatim (chunk, old line, new line, reason) in the report's "Article deviations". Never change what is computed, compared or seeded, which bound is which, or any literal in `reference-data`; if that would be needed, stop and report with a proposed fix. A second failure is a stop.

The render writes `vignettes/articles/guohe-nested-continuum.html` (and possibly a `_files/` directory — with `embed-resources: true` it should not; if one appears, leave it untracked and say so). Neither is staged.

## 5. Text gates on the render — GATE (v2: output only)

The article prints its code (`code-fold: false`), so the rendered HTML carries every chunk's source as well as its output. **Counts are taken on the output only.** Produce two files in `results/`:

- `article.txt` — `pandoc --wrap=none -t plain vignettes/articles/guohe-nested-continuum.html` (the whole page, kept for the record);
- `article_output.txt` — the same conversion with the echoed code removed: delete every `<pre …>…</pre>` element from the HTML before converting (or, equivalently, drop the code blocks from the plain text — pandoc renders them as indented blocks); state the method in the report and keep the script that does it in the record directory.

All counts below are on `article_output.txt`, with `grep -o … | wc -l`:

1. `Replicate 1, beta_2 = 0.3, at the manuscript.s settings` (regex; the apostrophe may be typographic) — exactly 1 (the walkthrough ran at B = 2000).
2. `This article.s run: one-sided 95% lower bounds` — exactly 1; `This article.s run: what the argmax selects` — exactly 1; `This article.s run: the Bonferroni pair` — exactly 1 (regexes, as above).
3. `Wall-clock of this article.s run on 2 workers` — exactly 1 (the plan resolved to two workers on this machine).
4. The published values are in the reference tables: `0.952 (0.941, 0.960)` — at least 2; `0.928 (0.916, 0.939)` — at least 1; `0.945 / 0.302 / +0.010` — exactly 1 (the first row of Table S23).
5. `NaN` — 0; `grep -cP '\bNA \('` — 0; `grep -oE '(^|[^[:alpha:]])-?Inf([^[:alpha:]]|$)'` — 0 (whole-word Inf only, on the output text; the `-Inf` in the article's code is not output).
6. Selection agreement: `grep -c 'Guo and He.s argmax and the multiplier resampling both carry the same candidate' results/article_output.txt` — 1 (the walkthrough's reading; the Monte Carlo's agreement is enforced by the article's `stopifnot`).
7. Numbered-list hazard, on the source: every line matching `^\s*[0-9]+[.)]( |$)` in `vignettes/articles/guohe-nested-continuum.qmd` must lie between the `## Limits` heading and the `## Running the manuscript's design` heading (the nine intended items). Any other such line is a stop.
8. The two figures rendered: `grep -c 'Bias of the point estimate against beta(H-hat)' results/article_output.txt` ≥ 1 and `grep -c 'The complement and the joint claim, by cell' results/article_output.txt` ≥ 1 (the captions).

Record each count in `keyvalues.txt`, alongside v1's whole-page counts for comparison (they differ from the output-only counts by exactly the echoed occurrences: +1 on gates 1, 2b, 2c and 5c, and nowhere else — record that difference). *GATE:* any output-only count outside its stated value → stop and report the count.

## 6. `_pkgdown.yml` — GATE

Two additions, nothing else:

1. Navbar, in `navbar: components: articles: menu:`, directly after the entry `- text: Post-Selection Inference` / `href: articles/survival-post-selection.html`:
   ```
           - text: Comparison with Guo & He (2021)
             href: articles/guohe-nested-continuum.html
   ```
   (same indentation as its neighbours).
2. Articles index, in the `articles:` list, section `title: Post-Selection Inference`, append to its `contents:`:
   ```
         - articles/guohe-nested-continuum
   ```
   after `- articles/survival-post-selection`.

Then `Rscript -e 'pkgdown::check_pkgdown()'` → `logs/check_pkgdown.txt`. *GATE:* must print `No problems found`. `git diff --stat _pkgdown.yml` must show one file changed, 3 insertions (±1 for a blank line), 0 deletions; record it.

## 7. Commit, report, post-conditions — GATE

**Resume rule.** If §5–§7 were already completed under Larry's option-1 instruction before this document arrived (the article and `_pkgdown.yml` committed, the report written): re-run nothing. Copy this document to `dev/tasks/cc_task_guohe_article_2026-10-30_v2.md`, commit it alone as the last commit, check that the counts recorded in `results/keyvalues.txt` satisfy this §5 (output-only values) and that `results/article_output.txt` exists (produce it from the existing render if the re-run used another route, and say so), and append one line to the report's provenance naming this document and the instruction it records. Otherwise: copy and commit this document alone now, then execute §5–§7 as written here on the existing render.

Commits by explicit path, in this order, no push, no fetch, no merge:
1. (done, v1 §1.3) the v1 task document alone (`4ea458bb`); this document alone, per the resume rule;
2. `vignettes/articles/guohe-nested-continuum.qmd` and `_pkgdown.yml` together — message `docs: article — post-selection bounds on the Guo & He (2021) nested continuum (S4 of 2026b)`;
3. `dev/reports/guohe_article_2026-10-30/` (gate script, `results/`, `logs/`) and `dev/reports/REPORT_guohe_article_2026-10-30.md`.

Before each commit `git diff --cached --name-only` must list only the paths named for it. After the last commit, print in the report's closing block:

```
git diff --name-only c2967362 HEAD | grep -vE '^(dev/|vignettes/articles/guohe-nested-continuum\.qmd$|_pkgdown\.yml$)'   # must print nothing
git status --porcelain | grep -v '^?? '                                                                          # must print nothing
git status --porcelain | grep '^?? '                                                                             # only vignettes/articles/guohe-nested-continuum.html (and a _files/ dir if one appeared)
```

If the first line prints anything, do not fix it — stop and report what it prints. Copy the report, `results/` and the rendered HTML to `~/Downloads/guohe_article_2026-10-30/`.

**Report — `dev/reports/REPORT_guohe_article_2026-10-30.md`:**
1. Provenance: BASE_SHA, HEAD after each commit, input md5, quarto/pandoc versions, cores, wall-clock of the render, the plan the article resolved (`2 workers`); v1's stop at its §5 and this document's output-only rule, with the method used to strip the echoed code.
2. Reference gate: totals (comparisons, passes), the corrections list (or "none"), the payload meta beside the article's stated settings.
3. Render: exit status, wall-clock, the article's own timing table transplanted from `results/article.txt` (the "Wall-clock of this article's run" table and the truth-curve seconds), and the one-replicate table transplanted (so the walkthrough's values are on record).
4. Text gates: each count.
5. `_pkgdown.yml`: the diff, `check_pkgdown()` output.
6. Article deviations (verbatim edits, or "none"); Script notes.
7. OPEN ITEMS (never a stop): the NEWS line after the 0.4.1 release; tracking of the rendered HTML (Larry's decision); anything the render surfaced.
8. Ten-line verdict.

## 8. Out of scope

- No `R/`, `man/`, `NEWS.md`, `DESCRIPTION`, `vignettes/forestsearch.qmd` edit; no install, check or tests.
- No full-size run; no change to the article's default settings; no change to its seeds, bounds or literals beyond §3's payload corrections and §4's deviation clause.
- No edit to any other article (the alignment task is separate). No push, no merge.
