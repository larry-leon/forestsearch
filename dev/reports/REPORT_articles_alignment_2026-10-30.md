# REPORT — align five pkgdown articles with León & Anderson 2026a/2026b (2026-10-30)

Task: `dev/tasks/cc_task_articles_alignment_2026-10-30_v2.md` (v2, authoritative; supersedes v1
`dev/tasks/cc_task_articles_alignment_2026-10-30.md`). Branch `feature/articles-alignment`.
Record directory: `dev/reports/articles_alignment_2026-10-30/` (`render_all.sh`,
`strip_and_gate.sh`, `results/`, `logs/`).

## 0. Guo–He HTML (v2 §0)

`vignettes/articles/guohe-nested-continuum.html` (the Guo–He task's untracked render) was
already absent on `feature/articles-guohe` before v1 began: `git status --porcelain` there
was empty at the start of v1.

What the listings showed:
- **Stash** (`logs/guohe_html_stash.txt`): two entries, both GitHub Desktop stashes from other
  branches (`stash@{0}` on `feature/sg-focus-transparency`, `stash@{1}` on `feature/glm-extension`);
  neither has a third (untracked) parent, and neither lists a `guohe` path. **No stash holds the file.**
- **Trash** (`logs/guohe_html_trash.txt`): older Guo–He working files (July–September: `_08_guohe_supplement.qmd`,
  `guohe_algorithm3.R*`, `test_guohe_algorithm3*.qmd`, …); **no `guohe-nested-continuum.html`.**
- **Downloads** (`logs/guohe_html_downloads.txt`): `~/Downloads/guohe_article_2026-10-30/guohe-nested-continuum.html`,
  1,774,977 bytes, the Guo–He task's §7 copy (same size as the render recorded there).

So the file was neither stashed nor moved to the trash by a GitHub Desktop "discard"; how it
left the working tree is not recoverable from these records.

**Source restored from: (b), the `~/Downloads` copy.** Verification:
`diff <(pandoc --wrap=none -t plain vignettes/articles/guohe-nested-continuum.html) <(git show feature/articles-guohe:dev/reports/guohe_article_2026-10-30/results/article.txt)`
printed nothing (`logs/guohe_html_verify_diff.txt`, 0 lines). `guohe_html_md5 = e6173419e0f3c07d64514f6a06e2b1f1`.
Untracked list at that point: exactly `?? dev/reports/articles_alignment_2026-10-30/` and
`?? vignettes/articles/guohe-nested-continuum.html`. Not staged; not re-rendered.

Order note: v2's resume rule lists §0 before the v2 commit; Larry's instruction was to commit
v2 alone "before anything else resumes", so the v2 commit came first. §0 writes only an
untracked file, so the order does not change any commit's content.

## 1. Provenance

| item | value |
|---|---|
| BASE_SHA | `c2967362e7ee9a6361cb1b07fdf83c0098e712db` (`master`, clean, ancestry confirmed; no fetch) |
| commit 1 (v1 task files) | `794863f3fb79805554b25937ff9425238e2daff1` |
| commit 1b (v2 task doc alone) | `bdeb75f8caebb32f39e2b0c19cce49156c1f03ac` |
| commit 2 (five articles) | `a88f8ea4df9dbacebcddbba85a6c76660aec9057` |
| commit 3 (record + report) | see closing block in the v2 task document |
| commit 4 (closing block, v2 doc) | see closing block |
| `alignment_edits_2026-10-30.json` | `d3859a8c37c2710edd2cef9b4966a664` (matched) |
| `apply_alignment_edits.py` | `1056eed3b82f1d889500c662649f6210` (matched) |
| v1 task doc (not gated) | `a41bb71b1611a3fe2c883152a4a4e75e` |
| v2 task doc (recorded) | `3e1034ebcd302004cfd3a3b9ad96b606` |
| forestsearch / rmarkdown / pkgdown | 0.4.1 / 2.32 / 2.2.1 |
| pandoc / python | 3.1.3 / Python 3.12.3 |

v1's §1 untracked list was **empty** (see §0). v1 stopped at §4 gate 12a (count 1 against an
expected 2); v2 re-states 12a as exactly 1 (`toc_float` builds its contents in JavaScript, so the
saved HTML has only an empty `<div id="TOC" class="tocify">`).

## 2. Apply (v1 §2)

```
apply_check.log: RESULT: OK   (check only; all five OK)
apply.log:       RESULT: OK
biomarker_effects.Rmd:            applied 3 edits;  md5 bdcd672a -> 42a79e8c
causal_effects_brief_review.Rmd:  applied 15 edits; md5 787e60f6 -> 94caca1f
methodology.Rmd:                  applied 14 edits; md5 2a311304 -> 5d79611e
paper_simulations.Rmd:            applied 7 edits;  md5 6f0282e5 -> 6154d7ea
treatment_effect_definitions.Rmd: applied 5 edits;  md5 4d2a5e6d -> 2c46b87f
```

`git diff --stat` (`logs/diff_stat.txt`):

```
 vignettes/articles/biomarker_effects.Rmd           |  42 ++++--
 vignettes/articles/causal_effects_brief_review.Rmd | 128 +++++++++++++------
 vignettes/articles/methodology.Rmd                 | 142 ++++++++++++++++++++-
 vignettes/articles/paper_simulations.Rmd           |  48 ++++++-
 .../articles/treatment_effect_definitions.Rmd      |  63 ++++++++-
 5 files changed, 365 insertions(+), 58 deletions(-)
```

No path outside `vignettes/articles/`. No hand edits. The md5s were re-confirmed unchanged at
the v2 resume (`5d79611e 6154d7ea 2c46b87f 42a79e8c 94caca1f`).

## 3. Render (v1 §3; not repeated under v2)

`rmarkdown::render(..., output_dir = "/tmp/fs_align_render_2026-10-30", intermediates_dir = ".../intermediates_<name>")`
via `render_all.sh`, each under `timeout` at its cap:

| article | exit | cap | wall-clock |
|---|---|---|---|
| methodology | 0 | 10 min | 4 s |
| treatment_effect_definitions | 0 | 10 min | 2 s |
| biomarker_effects | 0 | 10 min | 0 s |
| causal_effects_brief_review | 0 | 10 min | 2 s |
| paper_simulations | 0 | 45 min | 47 s (44 + 44 trials; gate 7a confirms the inline counts evaluated) |

`logs/status_before_render.txt` and `logs/status_after_render.txt` identical (five ` M` articles,
the record directory untracked). No `*_cache/` or `*_files/` directory appeared; nothing
deleted. `git diff --name-only -- 'vignettes/articles/*.html'` printed nothing.

## 4. Text gates (output only)

**Strip method** (`strip_and_gate.sh`, same as `dev/reports/guohe_article_2026-10-30/s2_text_gates.sh`):
every `<pre …>…</pre>` element deleted from the rendered HTML with R
`gsub("(?s)<pre[ >].*?</pre>", "", …, perl = TRUE)` (guarded by a check that none remain), then
`pandoc --wrap=none -f html -t plain` → `results/<name>_output.txt`; `<name>.txt` is the whole page.
`<pre>` elements removed: methodology 9, paper_simulations 251, treatment_effect_definitions 8,
biomarker_effects 0, causal_effects_brief_review 0.

| gate | article | string | count | expected | result |
|---|---|---|---|---|---|
| 1 | methodology | `2026 note` | 13 | = 13 | PASS |
| 2 | methodology | `consistency_method` | 2 | ≥ 1 | PASS |
| 2 | methodology | `fs_declaration_calibration` | 2 | ≥ 1 | PASS |
| 2 | methodology | `forestsearch_bootstrap_dofuture` | 1 | ≥ 1 | PASS |
| 2 | methodology | `hrMaxSG` | 2 | ≥ 1 | PASS |
| 3 | methodology | `1.276` / `0.634` / `0.645` / `0.804` | 1 / 1 / 1 / 1 | ≥ 1 each | PASS |
| 4 | methodology | `These are the 2024 guidelines` | 1 | = 1 | PASS |
| 5 | paper_simulations | `Identifying Patient Subgroups` | 0 | = 0 | PASS |
| 5 | paper_simulations | `Exploratory subgroup identification in the heterogeneous Cox model` | 1 | = 1 | PASS |
| 6 | paper_simulations | `arXiv:2609.38361` / `arXiv:2610.28379` | 2 / 2 | ≥ 2 each | PASS |
| 7 | paper_simulations | `44 simulated trials under the alternative and 44 under the null` | 1 | = 1 | PASS |
| 7 | paper_simulations | `n_sims = 1000` | 0 | = 0 | PASS |
| 8 | paper_simulations | `random-split consistency` | 0 | = 0 | PASS |
| 8 | paper_simulations | `consistency_method` | 1 | ≥ 1 | PASS |
| 9 | treatment_effect_definitions | `The Post-Selection Target: the Standard Analysis of the Selected Region` | 2 | = 2 | PASS |
| 9 | treatment_effect_definitions | `identification factor` | 1 | = 1 | PASS |
| 10 | treatment_effect_definitions | `default estimand for ForestSearch simulation` | 0 | = 0 | PASS |
| 10 | treatment_effect_definitions | `arXiv:2609.38361` / `arXiv:2610.28379` | 1 / 1 | ≥ 1 each | PASS |
| 11 | biomarker_effects | `create_gbsg_dgm` | 1 | ≥ 1 | PASS |
| 11 | biomarker_effects | `generate_aft_dgm_flex` | 2 | ≥ 2 | PASS |
| 11 | biomarker_effects | `arXiv:2609.38361` / `arXiv:2610.28379` | 1 / 1 | ≥ 1 each | PASS |
| **12a** | causal_effects_brief_review | `What forestsearch Reports, Read Against This Review` | **1** | v1: = 2 (FAIL); **v2: = 1** | **PASS (v2)** |
| 12b | causal_effects_brief_review | `not a claim of causal optimality` | 2 | = 2 | PASS |
| 13 | causal_effects_brief_review | `GRF/LASSO` / `CHR alongside AHR` / `biased comparator` / `GRF baseline covariate conditioning` | 0 / 0 / 0 / 0 | = 0 each | PASS |
| 14 | causal_effects_brief_review | `arXiv:2609.38361` / `arXiv:2610.28379` | 1 / 1 | ≥ 1 each | PASS |
| 15 | all five | `NaN` | 0 in each | = 0 | PASS (no master render needed) |
| 16 | all five (source) | numbered lines not on master | 0 in each, after the three allowed `paper_simulations` reference items `1.` (B6), `4.`, `5.` (B7) | 0 | PASS |

Gate 12a: v1 value 1 against an expected 2 (stop); v2 evaluation 1 against an expected 1, on
the existing `causal_effects_brief_review_output.txt` (no re-render, no re-conversion). The other
51 counts are v1's, not recomputed. Records: `results/keyvalues.txt`, `logs/text_gates.log`.

## 5. `pkgdown::check_pkgdown()`

```
✔ No problems found.
```

`git diff --name-only -- _pkgdown.yml` printed nothing (`logs/check_pkgdown.txt`).

## 6. Article deviations; script notes

**Article deviations: none.** The applier was the only writer of the five articles.

Script notes:
- `render_all.sh` renders the five in the stated order, each under `timeout` at its cap, stopping
  at the first non-zero exit; it writes `logs/render_<name>.log`, `logs/time_<name>.txt` and
  `logs/render_summary.txt`.
- `strip_and_gate.sh` writes `results/<name>.txt` / `<name>_output.txt` and prints every §4 count;
  gate 16 compares each numbered source line with `git show c2967362:vignettes/articles/<name>.Rmd`
  (exact-line match). It was run once under v1; under v2 only 12a was re-evaluated, inline.
- `keyvalues.txt` keeps v1's block (including its RED line) followed by the v2 lines
  (`guohe_html_md5`, `g12a_v2`, the v2 GREEN line).

## 7. OPEN ITEMS (not stops)

1. **`methodology.Rmd` citeproc warning, in text this task did not touch:** the render log has
   `[WARNING] Citeproc: citation eq-power not found`. Line 916 reads "The power approximation
   from @eq-power …" — a Quarto-style equation cross-reference, which R Markdown's citeproc reads
   as a citation key. It is on `master` too (`c2967362`), so it predates this task; the rendered
   page reads "from (eq-power?) was" in place of a link.
2. **Tracking of `vignettes/articles/guohe-nested-continuum.html`** on `feature/articles-guohe`
   (Larry's call). As an untracked file it has already disappeared once without a stash or
   trash trace; tracked, it is stale against its `.qmd` on every build (the reasoning of
   `c2967362`). The same staleness applies to the four tracked article HTMLs
   (`forestsearch-2024`, `methodology`, `paper_simulations`, `survival-post-selection`), now
   behind the aligned `methodology.Rmd` and `paper_simulations.Rmd`.
3. **`toc_float` articles:** a heading counts once in output-only text; future gates should
   expect 1.
4. The v1 estimate for `paper_simulations` (cap 45 min) was conservative: 47 s.

## 8. Verdict

1. BASE_SHA `c2967362` on a clean `master`; branch `feature/articles-alignment`; no fetch.
2. Task inputs' md5s matched; the v1 files and the v2 document were each committed alone.
3. The applier applied 44 edits over five articles, `RESULT: OK` twice, all md5 transitions as stated.
4. Five out-of-tree renders, all exit 0, under a minute in total; the repo was unchanged by them.
5. 52 output-only counts pass: 51 under v1, gate 12a under its v2 re-statement (= 1).
6. `check_pkgdown()` clean; `_pkgdown.yml` untouched.
7. No article deviations; no package code, no tracked HTML touched.
8. The Guo–He HTML was restored from the `~/Downloads` copy and verified byte-for-text (empty diff).
9. One untouched-text issue surfaced (`@eq-power` in `methodology.Rmd`, pre-existing).
10. No push, no fetch, no merge; ready for Larry's merge decision.
