# CC task — align five pkgdown articles with León & Anderson 2026a/2026b: apply the anchored edit set, render, gate, commit — v2

**File:** `dev/tasks/cc_task_articles_alignment_2026-10-30_v2.md` · **Issued:** 2026-10-30 by chat, commissioned by Larry ("I would like the alignment with the manuscripts to be comprehensive"; frame approved 2026-10-30: the two 2024 articles stay the historical record with every divergence stated explicitly in place, the three conceptual articles have their conflicting guidance rewritten; 2026-10-30 13:37: the Guo–He article's rendered HTML "needs to be restored").
**Repo:** `forestsearch` · **Branch:** `feature/articles-alignment` (created from `master` at `c2967362` = BASE_SHA; HEAD `794863f3`) · **Transport:** `~/Downloads` → `dev/tasks/` for this document, committed alone before anything else resumes. The other two task files (`alignment_edits_2026-10-30.json`, `apply_alignment_edits.py`) are already committed at `794863f3`.

**v2 supersedes v1** (`dev/tasks/cc_task_articles_alignment_2026-10-30.md`, committed `794863f3`). v1 stopped at its §4 gate 12a: it expected the new closing heading of `causal_effects_brief_review.Rmd` to count 2 (heading plus table-of-contents entry), but that article's YAML sets `toc_float`, whose contents list is built by JavaScript at page load and is absent from the saved HTML, so the count is 1. The render is correct; the gate was wrong. The other 51 counts passed. CC's stop was right. **This document is the complete, authoritative task.** What differs from v1: the new §0 (restore the Guo–He article's HTML, which was found missing before v1 began), gate 12a in §4, the resume rule and the untracked-list expectation in §6, and the report's §0 item. Nothing is re-rendered.

**Inputs (md5; verified under v1 §1).** `alignment_edits_2026-10-30.json` `d3859a8c37c2710edd2cef9b4966a664` (44 edits over 5 files); `apply_alignment_edits.py` `1056eed3b82f1d889500c662649f6210`. The edit set was built by chat from the five articles exactly as committed on `master` at `c2967362` (`vignettes/articles/`: `methodology.Rmd` `2a311304…`, `paper_simulations.Rmd` `6f0282e5…`, `treatment_effect_definitions.Rmd` `4d2a5e6d…`, `biomarker_effects.Rmd` `bdcd672a…`, `causal_effects_brief_review.Rmd` `787e60f6…`) and verified against them: every anchor occurs exactly once, every replacement zero times, the applier is idempotent, and it refuses a file that is neither the master text nor the already-edited text. Every statement in the edits was checked by chat against the manuscripts' text (2026a arXiv:2609.38361 §1.2, §2.1–2.5, §3, §5, §6, §7, §8, S1, S3, S7, S8, S9; 2026b arXiv:2610.28379 §2.2, §2.3, §3, §4, §5.1, §6.1, §7, S1.8.4, S2.3, S3.4) and against the package source on `master` (`consistency_method` default and fallback, `sg_focus` spellings, `selection_rule`, `effect_neighborhood`, `field_recovery`, the `forestsearch_bootstrap_dofuture` export, the columns `run_simulation_analysis()` records, the 0.4.1 NEWS entry on the variable-importance screen). §A lists every edit with its basis.

**State at v2 (from v1's stop report).** BASE_SHA `c2967362`; `794863f3` = the three v1 task files; the five articles edited by the applier (md5 transitions as §2) and the record directory `dev/reports/articles_alignment_2026-10-30/` (`render_all.sh`, `strip_and_gate.sh`, `results/`, `logs/`) present and uncommitted; all five renders in `/tmp/fs_align_render_2026-10-30/` exited 0 (methodology 4 s, treatment_effect_definitions 2 s, biomarker_effects 0 s, causal_effects_brief_review 2 s, paper_simulations 47 s); `git status` identical before and after the renders; no deviations. The untracked list recorded in v1 §1 was **empty**: `vignettes/articles/guohe-nested-continuum.html` (the Guo–He task's render, left untracked by design) was already absent on `feature/articles-guohe` before v1 began; a copy was placed in `~/Downloads/guohe_article_2026-10-30/` by the Guo–He task's §7.

**What this task produces.** The five edited articles in `vignettes/articles/`, one local render of each with output-only text gates, `pkgdown::check_pkgdown()`, a record directory and a report; and the Guo–He HTML restored to `vignettes/articles/` (untracked, as the Guo–He task left it). Nothing else.

---

## ⚠ CATEGORY

**No package code is modified.** No edit under `R/`, `tests/`, `man/`, `inst/`, `DESCRIPTION`, `NAMESPACE`, `NEWS.md`, `vignettes/forestsearch.qmd`, `_pkgdown.yml` (no navbar change: the five articles are already listed), or anything under `quarto/`. No `R CMD check`, no test suite, no install. **No hand edits to the five articles**: the only writer is `apply_alignment_edits.py`; if it fails, the task stops (§2). Chat could not run R, so the render (§3) is the first execution of the edited articles; the deviation clause in §3 exists for a render defect of the inserted text only (a markdown or math construct pandoc rejects), never for content.

**Writes:** this document under `dev/tasks/` (the other two task files are committed); the five `.Rmd` files under `vignettes/articles/`; the record directory `dev/reports/articles_alignment_2026-10-30/` (`results/`, `logs/`, the scripts); the report `dev/reports/REPORT_articles_alignment_2026-10-30.md`; and, by §0, the untracked file `vignettes/articles/guohe-nested-continuum.html` restored byte-for-byte from a verified copy. **Rendered HTML of the five articles is never staged and never written into the repo**: the renders are in `/tmp/fs_align_render_2026-10-30/` and a copy goes to `~/Downloads/articles_alignment_2026-10-30/`. The tracked `vignettes/articles/*.html` files are not touched. The restored Guo–He HTML is not staged either (tracking it is Larry's decision, §6 OPEN ITEMS).

**Compute (its own line):** none in v2. The five renders of v1 §3 stand (total under one minute; v1's 45-min cap for `paper_simulations.Rmd` was conservative — it took 47 s with its 44 + 44 trials, confirmed by gate 7a). Nothing is re-rendered; the Guo–He article is **not** re-rendered (a re-render would change its wall-clock table and is a separate go/no-go).

**Unattended.** Gates stop the task; never ask, never work around. Pre-existing untracked files are never staged. No push, no fetch, no merge.

---

## 0. Restore the Guo–He article's HTML — GATE (v2; runs first, before the resume)

The file to restore is `vignettes/articles/guohe-nested-continuum.html`, the render the Guo–He task (`dev/tasks/cc_task_guohe_article_2026-10-30_v2.md` on `feature/articles-guohe`) left untracked in `vignettes/articles/`. Untracked files are not protected by git: a GitHub Desktop branch switch that "leaves changes" stashes them, and a "discard changes" moves them to the trash. Record first, then restore from the first source that holds it; the restored bytes are verified against the committed plain-text conversion of the original render.

1. **Record (no stop):** `git stash list` verbatim; for each stash, `git stash show --include-untracked --name-only stash@{n}` (or `git show --stat stash@{n}^3` where the third parent exists) → `logs/guohe_html_stash.txt`. `ls -la ~/.local/share/Trash/files/ 2>/dev/null | grep -i guohe` → `logs/guohe_html_trash.txt`. `ls -la ~/Downloads/guohe_article_2026-10-30/` → `logs/guohe_html_downloads.txt`. These say what happened; they are recorded in the report, not acted on beyond step 2.
2. **Restore, in this order of preference, from the first that exists:** (a) a stash entry holding the file — `git checkout stash@{n}^3 -- vignettes/articles/guohe-nested-continuum.html`, then `git restore --staged vignettes/articles/guohe-nested-continuum.html` so it is untracked again (checking out from a stash's untracked tree stages the file); do **not** drop or apply the stash; (b) `cp ~/Downloads/guohe_article_2026-10-30/guohe-nested-continuum.html vignettes/articles/`; (c) the trash copy. Say which source was used. If none holds it: stop and report — a re-render is Larry's go/no-go.
3. **Verify — GATE:** `diff <(pandoc --wrap=none -t plain vignettes/articles/guohe-nested-continuum.html) <(git show feature/articles-guohe:dev/reports/guohe_article_2026-10-30/results/article.txt)` must print nothing (same pandoc, same bytes, same text). Record `md5sum vignettes/articles/guohe-nested-continuum.html` in `results/keyvalues.txt` as `guohe_html_md5`, and `git status --porcelain | grep '^?? '` — must list exactly two entries at this point: `?? dev/reports/articles_alignment_2026-10-30/` (uncommitted until §6) and `?? vignettes/articles/guohe-nested-continuum.html`. A non-empty diff → stop (the copy is not the render). The file is never staged.

---

## 1. Provenance and first commit — GATE (done under v1; recorded here)

Completed under v1 and not repeated: BASE_SHA `c2967362` on `master` (clean, ancestry confirmed); `feature/articles-alignment` created; the JSON and applier md5s matched; the three v1 files committed alone as `794863f3`; versions forestsearch 0.4.1, rmarkdown 2.32, pkgdown 2.2.1, pandoc 3.1.3, Python 3.12.3. The untracked list was empty (see §0).

**v2 first commit.** `md5sum ~/Downloads/cc_task_articles_alignment_2026-10-30_v2.md` is recorded (not gated). Copy it to `dev/tasks/cc_task_articles_alignment_2026-10-30_v2.md`, `git add` by path, commit it alone: `task: article alignment — v2 (gate 12a re-stated for toc_float; Guo–He HTML restore)`. `git diff --cached --name-only` must list only that path. The five edited articles and the record directory stay uncommitted until §6.

## 2. Apply the edit set — GATE (done under v1; recorded here)

Completed under v1: the check run and the real run both ended `RESULT: OK` with the five md5 transitions below; the diff is 5 files, 365 insertions, 58 deletions, all under `vignettes/articles/`; nothing edited by hand. Not repeated. For the record, the commands were, from the repo root:

```
python3 -I dev/tasks/apply_alignment_edits.py dev/tasks/alignment_edits_2026-10-30.json vignettes/articles --check 2>&1 | tee dev/reports/articles_alignment_2026-10-30/logs/apply_check.log
python3 -I dev/tasks/apply_alignment_edits.py dev/tasks/alignment_edits_2026-10-30.json vignettes/articles         2>&1 | tee dev/reports/articles_alignment_2026-10-30/logs/apply.log
```

*GATE:* both runs end in `RESULT: OK`; the second prints `applied` for all five files with the md5 transitions `2a311304→5d79611e`, `6f0282e5→6154d7ea`, `4d2a5e6d→2c46b87f`, `bdcd672a→42a79e8c`, `787e60f6→94caca1f`. Any `FAIL` line → stop and report it verbatim (a `FAIL md5 … neither master nor applied` means the article on `master` is not the text the edits were cut against; do not edit by hand).

Then `git diff --stat` → `logs/diff_stat.txt` (five files, insertions only except the lines §A marks as rewrites), and `git diff --name-only | grep -v '^vignettes/articles/'` must print nothing. Do not commit yet.

## 3. Render — GATE (done under v1; recorded here)

Completed under v1: all five renders exited 0 into `/tmp/fs_align_render_2026-10-30/` (times in the header), `git status` identical before and after, no cache directories, tracked HTML unchanged, no deviations. **Not repeated**; if `/tmp/fs_align_render_2026-10-30/` no longer exists (a reboot), stop and report — do not re-render under v2. For the record, the procedure was: `git status --porcelain` before rendering (`logs/status_before_render.txt`); for each of the five, from the repo root, under `timeout` at its cap:

```r
rmarkdown::render("vignettes/articles/<name>.Rmd",
                  output_dir = "/tmp/fs_align_render_2026-10-30",
                  intermediates_dir = "/tmp/fs_align_render_2026-10-30/intermediates_<name>",
                  quiet = FALSE)
```

with the console to `logs/render_<name>.log` and wall-clock recorded (`date` before and after). Order: `methodology`, `treatment_effect_definitions`, `biomarker_effects`, `causal_effects_brief_review`, `paper_simulations`. The input directory stays the knit root, so `references.bib` resolves. *GATE:* non-zero exit → stop, with the log's last 40 lines in the report — **except** under the deviation clause below. A timeout is a stop.

After all five, `git status --porcelain` again (`logs/status_after_render.txt`). The render must leave the repo as §2 left it: if a `*_cache/` or `*_files/` directory appeared under `vignettes/articles/`, delete it, say so in the report, and confirm the two status files then agree. The tracked `vignettes/articles/*.html` must be unmodified (`git diff --name-only -- 'vignettes/articles/*.html'` prints nothing).

**Deviation clause (render defects of the inserted text only).** If a render fails on a pandoc or knitr error located in text this task inserted — an unbalanced `$`, a blockquote line pandoc rejects, an inline `` `r …` `` that does not evaluate — make the smallest edit that restores the intended rendering, record it verbatim (file, old line, new line, reason) under "Article deviations", and re-render that article once. Never change a number, a citation, a section reference or a sentence's meaning; if that would be needed, stop and report with a proposed fix. A failure in text this task did not insert, or a second failure, is a stop.

## 4. Text gates on the renders — GATE (output only; v2: gate 12a re-stated)

Under v1, 51 of the 52 counts passed on the existing `results/<name>_output.txt` files; only gate 12a failed, for the reason given in the header. **Under v2, re-evaluate gate 12a alone on the existing `results/causal_effects_brief_review_output.txt`** (no re-render, no re-conversion), record it in `results/keyvalues.txt` beside v1's value, and leave the other 51 counts as recorded. The full gate set is restated below as the authoritative list.

The articles echo some code and fold others, so counts are taken on the **output only**. For each rendered HTML, `results/` holds:

- `<name>.txt` — `pandoc --wrap=none -t plain <html>` (whole page, for the record);
- `<name>_output.txt` — the same conversion after every `<pre …>…</pre>` element is removed from the HTML (dot-all regex; the same method as `dev/reports/guohe_article_2026-10-30/`, whose strip script may be copied — keep the script used in the record directory and state it in the report).

All counts on `<name>_output.txt` with `grep -o … | wc -l` (fixed strings unless marked regex):

**methodology**
1. `2026 note` — exactly 13 (the banner mention and the twelve blockquote notes).
2. `consistency_method` — at least 1; `fs_declaration_calibration` — at least 1; `forestsearch_bootstrap_dofuture` — at least 1; `hrMaxSG` — at least 1.
3. The GBSG figures of 2026b Table 1 in the Applications note: `1.276`, `0.634`, `0.645`, `0.804` — at least 1 each (the article's own computed tables may repeat a value by coincidence, so no upper bound).
4. `These are the 2024 guidelines` — exactly 1.

**paper_simulations**
5. `Identifying Patient Subgroups` — 0 (the wrong 2024 citation is gone); `Exploratory subgroup identification in the heterogeneous Cox model` — exactly 1.
6. `arXiv:2609.38361` — at least 2; `arXiv:2610.28379` — at least 2 (banner and references).
7. `44 simulated trials under the alternative and 44 under the null` — exactly 1 (the inline `r nsims_alt` / `r nsims_null` evaluated); `n_sims = 1000` — 0.
8. `random-split consistency` — 0; `consistency_method` — at least 1.

**treatment_effect_definitions**
9. `The Post-Selection Target: the Standard Analysis of the Selected Region` — exactly 2 (the new section heading, repeated once by the table of contents that `toc: true` adds); `identification factor` — exactly 1.
10. `default estimand for ForestSearch simulation` — 0; `arXiv:2609.38361` — at least 1; `arXiv:2610.28379` — at least 1.

**biomarker_effects**
11. `create_gbsg_dgm` — at least 1; `generate_aft_dgm_flex` — at least 2; `arXiv:2609.38361` — at least 1; `arXiv:2610.28379` — at least 1.

**causal_effects_brief_review**
12. `What forestsearch Reports, Read Against This Review` — **exactly 1** (v2). The heading appears once in the body; this article's YAML sets `toc_float`, so its contents list is built by JavaScript at page load and the saved HTML carries only an empty `<div id="TOC">` — there is no static entry to count. (v1 said 2 and was wrong; `treatment_effect_definitions.Rmd` has a static `toc: true` list, which is why its gate 9a counts 2.) `not a claim of causal optimality` — exactly 2 (the Estimand Selection Guide and the closing section).
13. `GRF/LASSO` — 0; `CHR alongside AHR` — 0; `biased comparator` — 0; `GRF baseline covariate conditioning` — 0.
14. `arXiv:2609.38361` — at least 1; `arXiv:2610.28379` — at least 1.

**all five**
15. `NaN` — expected 0 in every article. A non-zero count is a stop unless the same count appears in a render of that article's `master` text (`git show c2967362:vignettes/articles/<name>.Rmd` rendered to `/tmp/fs_align_render_2026-10-30/master/`): equal counts pass and are recorded, a larger count is a stop.
16. Numbered-list hazard, on the source: every line of each edited `.Rmd` matching `^\s*[0-9]+[.)]( |$)` must exist verbatim on the `master` version of the file (`git show c2967362:vignettes/articles/<name>.Rmd`), except in `paper_simulations.Rmd` the three reference items `1.` (rewritten by B6), `4.` and `5.` (added by B7). Any other such line is a stop (an inserted sentence that begins a line with a number would render as a list).

Record each count in `results/keyvalues.txt`. *GATE:* any count outside its stated value → stop and report the count. (Under v2 the only count evaluated afresh is 12a; the 51 recorded under v1 are not recomputed.)

## 5. `pkgdown::check_pkgdown()` — GATE

`Rscript -e 'pkgdown::check_pkgdown()'` → `logs/check_pkgdown.txt`. *GATE:* must print `No problems found`. `_pkgdown.yml` is not edited by this task; `git diff --name-only -- _pkgdown.yml` must print nothing.

## 6. Commit, report, post-conditions — GATE

**Resume rule (v2).** Order of execution under v2: §0 (restore and verify the Guo–He HTML) → §1's v2 first commit (this document alone) → §4 gate 12a on the existing output file → §5 → the commits below. Nothing of v1's §2–§3 is re-run; if any v1 artefact needed by a step is missing (the edited articles, `results/*_output.txt`, the render directory), stop and report rather than regenerate.

Commits by explicit path, in this order, no push, no fetch, no merge:
1. (done, v1 §1.4) the three v1 task files alone (`794863f3`); (done, §1 above) this document alone;
2. the five articles together — `docs: align methodology, paper_simulations, treatment_effect_definitions, biomarker_effects and causal_effects_brief_review with León & Anderson 2026a/2026b`;
3. `dev/reports/articles_alignment_2026-10-30/` (`results/`, `logs/`, `render_all.sh`, `strip_and_gate.sh`) and `dev/reports/REPORT_articles_alignment_2026-10-30.md`;
4. the closing block appended to this document (v2) under `dev/tasks/`, committed alone. The v1 document is left as committed.

Before each commit `git diff --cached --name-only` must list only the paths named for it. After the last commit, print in the report's closing block:

```
git diff --name-only c2967362 HEAD | grep -vE '^(dev/|vignettes/articles/(methodology|paper_simulations|treatment_effect_definitions|biomarker_effects|causal_effects_brief_review)\.Rmd$)'   # must print nothing
git status --porcelain | grep -v '^?? '                                                                                                                              # must print nothing
git status --porcelain | grep '^?? '                                                                                                                                 # must print exactly: ?? vignettes/articles/guohe-nested-continuum.html (restored by §0, never staged)
git log --oneline c2967362..HEAD                                                                                                                                      # the five commits of this task: 794863f3, this document, the articles, the record, the closing block
```

If the first line prints anything, do not fix it — stop and report what it prints. Copy the report, `results/` and the five rendered HTML files to `~/Downloads/articles_alignment_2026-10-30/`.

**Report — `dev/reports/REPORT_articles_alignment_2026-10-30.md`:**
0. Guo–He HTML (§0): what the stash, trash and Downloads listings showed (so Larry knows what happened), the source restored from, `guohe_html_md5`, and the empty verification diff.
1. Provenance: BASE_SHA, HEAD after each commit (v1's `794863f3`, then v2's four), the three input md5s and this document's, package/rmarkdown/pkgdown/pandoc versions; v1's stop at gate 12a and v2's re-statement.
2. Apply: both applier logs' `RESULT` lines and the five md5 transitions; `git diff --stat`.
3. Render: per article exit status and wall-clock; the status-before/after comparison; any cache directory deleted.
4. Text gates: every count of §4, in order, with the strip method; gate 12a's v1 value (1, against an expected 2) and its v2 evaluation (1, expected 1).
5. `check_pkgdown()` output.
6. Article deviations (verbatim edits, or "none"); Script notes.
7. OPEN ITEMS (never a stop): anything the renders surfaced in text this task did not touch (a broken link, a stale number); whether `vignettes/articles/guohe-nested-continuum.html` should be tracked on `feature/articles-guohe` (Larry's call — as an untracked file it can vanish again on the next stash or discard; the four tracked article HTMLs are stale relative to their `.Rmd` on every build, which is the reasoning of `c2967362` for removing the vignette's pre-build render); for future gates on `toc_float` articles, a heading counts once.
8. Ten-line verdict.

## 7. Out of scope

Everything not named: `forestsearch.qmd`, `survival-post-selection.qmd`, `forestsearch-2024.Rmd`, `extreme_subgroups.qmd`, `guohe-nested-continuum.qmd`, the site build, the DINA floor check (`hr.threshold = 0.90` vs 2026a §6's "floors of 0 and log 1.0"; standing OPEN ITEM before any edit of the post-selection article), the NEWS line for the Guo–He article, and merging either feature branch.

---

## A. The edit set — what each edit changes and its basis

IDs match `alignment_edits_2026-10-30.json` (`edits[].id`); the JSON carries the exact anchor and replacement text. "Note" = a `> **2026 note.**` blockquote inserted after the anchored passage, leaving the 2024 text in place. "Rewrite" = the anchored passage is replaced.

### methodology.Rmd — the 2024 method, kept as record (14 edits, all notes except A1 and A14)

| ID | Where | What the article said | What the manuscripts say | Edit |
|---|---|---|---|---|
| A1 | banner | — | — | one sentence announcing the in-place notes |
| A2 | Type-1 error, scenario (ii) | per-candidate Type-1 framing | 2026b §4: the family-wise declaration rate under a uniform benefit $c_0$; diagnostic and calibrated cutoff; `fs_declaration_calibration()`; run in `vignette("forestsearch")` | note |
| A3 | Marginal vs CDE | "targets marginal hazard ratios" | 2026a §1.2, eq. 15, Thm 2, §7/S8; 2026b §1: the target is $\beta(\widehat H)$, the least-false within-region coefficient, family held fixed; data-adaptive; what the MR interval covers | note |
| A4 | Step 3c consistency rate | Monte Carlo over R splits | 2026b §2.2 and 2026a S3: closed form $\max\{0, 2\Phi((\hat\beta - c_{\text{cons}})/\sigma_{D,g}) - 1\}$; `consistency_method = "resample"` default, `fs.splits` only on the fallback (verified in `R/forestsearch_main.R` roxygen); complementary-pairs stability selection; the $(1+p^\star)/2$ lower-bound reading | note |
| A5 | selection table | 2024 criteria | 2026a §6 (10% neighbourhood), 2026b §6.1 (max N(ε) at 0.20), §2.3 (an operating-characteristic choice, no optimality claimed); spellings `hrMaxSG`/`effMaxSG`, `selection_rule`, `effect_neighborhood`, `maxeffCons`, `maxeff` (verified in source); default `"hr"` | note |
| A6 | threshold rationale table | per-candidate Type-1/power | 2026b §4.1–4.2: law of the maximum of $T_g$, $\hat\kappa_\alpha(c_0)$, the diagnostic; thresholds of the applications (HR 1.0, $p^\star$ 0.90: 2026a §6, 2026b §6.1) vs package defaults | note |
| A7 | bootstrap CI paragraph | bootstrap-only inference | 2026a §2.5, §3, Thm 1 (MR reproduces the correction to first order), §6 (seconds vs minutes), §5 (fixed vs model-generated family); 2026b §3 (field bounds and pair from the same draws); `forestsearch_bootstrap_dofuture()` exported | note |
| A8 | cross-validation | CV as the stability tool | the post-selection article's $\hat p(\widehat H)$ and `field_recovery = TRUE` (verified in `R/fs_mr_inference.R`), read as flags | note |
| A9 | Key results: estimation | 2024 OC table | 2026b §5/S2 and 2026a §7/S8 report the current inference's OCs; their designs (2026b §5.1: 2,000 replicates, conventional thresholds, no strict-null cell; 2026a S8: $\log 0.9$, $\log 0.8$) are not comparable with the table | note |
| A10 | Applications: GBSG | 2024 GBSG result (0.64 vs 0.69, "statistically significant") | 2026b Table 1 / §6.1–6.2: region {er ≤ 0} ∩ {pgr ≤ 26} (75), corrected 1.276, field lower 0.634; complement 0.645, field-s upper 0.804; whole trial 0.695; read by location, not as a test at HR = 1 | note |
| A11 | Applications: ACTG-175 | 2024 time-to-event benefit search | 2026a S7/Table S3 and 2026b S3.4: binary (OR of no CD4 improvement) and continuous (mean CD4 change) analyses in the harm orientation; no counterpart of the 2024 analysis | note |
| A12 | variable-selection table | LASSO/GRF recommendations | 2026a §5 and S1 (fully enumerated family recommended; a prognostic lasso prefilter shapes the family, conditional target); NEWS 0.4.1 (VI screen removed); 2026a S8.7–S8.8 (noise covariates) | note |
| A13 | timing table | 2024 timings | the vignette's timed fit (under a minute, MR and field included); 2026a §6 (FB 72 s–10.1 min vs MR 0.11–0.38 s) | note |
| A14 | Interpretation Guidelines | 2024 list | one paragraph before the list: the current products (corrected estimate, field lower / field-s upper, Bonferroni pair), read by location (2026b §3, §6.1) | paragraph inserted; list kept |

### paper_simulations.Rmd — the 2024 simulation framework (7 edits)

| ID | Where | What the article said | What the manuscripts say | Edit |
|---|---|---|---|---|
| B1 | banner | — | 2026b §5.1, §6.1; 2026a §7: the settings here (1.25, LASSO+GRF, `fs.splits = 400`; verified in the article's chunks) are not the 2026 studies' | sentence appended to the banner |
| B2 | after the capability list | "Estimate Type I error under null" | 2026a §8 reporting guidance (estimate the declaration rate by simulation calibrated to one's own trial, report it beside the de-biased effect); 2026b §4 diagnostic; S1.8.4 | paragraph inserted |
| B3 | Step 2 | "Because ForestSearch uses random-split consistency evaluation…" | 2026b §2.2; `consistency_method = "resample"` default, not overridden by `run_simulation_analysis()` (verified) | rewrite |
| B4 | Step 4, `build_estimation_table()` | bias-corrected rows | `run_simulation_analysis()` records `hr.H.hat`/`hr.Hc.hat`; `hr.H.bc` rows are the bootstrap correction when present (`R/simulation_tables.R`); no MR column; OCs of MR in 2026b §5/S2, 2026a §7/S8 | sentences appended |
| B5 | comparing with benchmarks | "`n_sims = 1000`" | the chunk sets `nsims_alt <- 44L`, `nsims_null <- 44L` | rewrite to inline `r nsims_alt` / `r nsims_null` |
| B6 | references | "Leon LF, Marceau-West CT, He W, et al. (2024). Identifying Patient Subgroups…" (wrong) | the published record: León, Jemielita, Guo, Marceau West, Anderson (2024) SIM 43(20):3921–3942, doi 10.1002/sim.10163 | rewrite |
| B7 | references | — | 2026a, 2026b | items 4–5 appended |

### treatment_effect_definitions.Rmd — estimands of the DGM (5 edits)

| ID | Where | What the article said | What the manuscripts say | Edit |
|---|---|---|---|---|
| C1 | introduction | three measures as "the" estimands | 2026a §1.2: design quantities vs the reported estimand $\beta(\hat H)$ | paragraph inserted |
| C2 | new section before "Guidance for Choosing an Estimand" | — | 2026a §1.2, eq. 15, Thm 2, §7, S8, S8.1 (random across replicates; sits below the planted marginal HR; $\kappa = \beta(H)/\beta(\hat H)$); 2026b §3, S2.3 | section "The Post-Selection Target: the Standard Analysis of the Selected Region" |
| C3 | guidance: AHR | "Use AHR as the default estimand…" | design quantity; the reported inference targets $\beta(\hat H)$ (2026a §1.2) | rewrite |
| C4 | Role in OCs | columns table | these are `run_simulation_analysis()`'s columns; MR not among them (verified) | sentences appended |
| C5 | references | 2024 only | + 2026a, 2026b | items appended |

### biomarker_effects.Rmd — the Weibull AFT/Cox DGM (3 edits)

| ID | Where | What the article said | What the manuscripts say | Edit |
|---|---|---|---|---|
| D1 | introduction | — | 2026a §7 and S9 (`create_gbsg_dgm()`, `generate_aft_dgm_flex()` named by the manuscript); 2026b §5.1 (GBSG-built design) | paragraph inserted |
| D2 | "Within the ForestSearch framework" bullets | "evaluation uses the AHR and CDE"; "Bootstrap bias correction (`hr.H.bc`)" | 2026a §1.2, §3, §7, §8, S8, Thm 1: reported estimand $\beta(\hat H)$ by institutional choice; evaluation against both target families; bootstrap reproduced to first order by MR | rewrite of the four bullets |
| D3 | references | 2024 (title truncated) | full 2024 record + 2026a, 2026b | rewrite + items appended |

### causal_effects_brief_review.Rmd — causal reading of hazard ratios (15 edits)

| ID | Where | What the article said | What the manuscripts say | Edit |
|---|---|---|---|---|
| E1 | purpose blockquote | — | signpost to the closing section | sentence appended |
| E2–E6 | mapping table (`paste0()` strings, rows 1, 2, 3, 6, 8) | "observed Cox HR treated as biased comparator"; "CHR reported as primary summary"; "report CHR alongside AHR"; "GRF/LASSO … conditions on baseline X_i"; "GRF baseline covariate conditioning is justified" | 2026a §1.2 (reported estimand), §2.1 (family from covariate cuts), §5 (enumerated family recommended); 2026b §3 (bounds); VI screen removed 0.4.1 | rewrite of the five strings (R syntax preserved) |
| E7 | guide: primary simulation estimand | "Primary simulation estimand: AHR" | design quantity; scored also against $\beta(\hat H)$ (2026a §7) | rewrite |
| E8 | guide: applied reporting | "CHR alongside AHR" | 2026a §1.2 (institutional choice, "not a claim of causal optimality"), §8 (interpretive limits, non-collapsibility); 2026b §3; CHR/AHR may accompany | rewrite |
| E9 | guide: GRF conditioning | "GRF baseline covariate conditioning is theoretically justified" | conditioning through the subgroup definition; 2026a §5; VI screen removed 0.4.1 | rewrite |
| E10 | guide: per-replicate Cox HR | — | the selection part is what MR removes (2026a §3, Thm 2) | sentences appended |
| E11 | Key Takeaways item 6 | "(GRF/LASSO in ForestSearch)" | through the covariate cuts that define a subgroup | rewrite |
| E12 | Key Takeaways item 10 | — | 2026a §1.2, §8 | sentence appended |
| E15 | summary table (gt), Fay & Li vaccine row | "(GRF/LASSO)" | as E11 | rewrite of the string |
| E13 | new closing section before References | — | 2026a §1.2, §8 (two short verbatim phrases, checked); 2026b §6.1 (reading by location) | section "What forestsearch Reports, Read Against This Review" |
| E14 | references | — | + 2026a, 2026b after Knudsen (2025) | items inserted |

**Deliberately unchanged.** The 2024 text of `methodology.Rmd` and `paper_simulations.Rmd` (figures, tables, thresholds, the bootstrap and CV sections, the 2024 application numbers); every passage of the three conceptual articles that is about the DGM, the AHR/CDE theory or the reviewed papers rather than about what the package reports; `references.bib` (the 2026 papers are cited inline, as the banners already do). `methodology.Rmd` says "Forest Search targets marginal hazard ratios" and keeps saying it — A3 states beside it what the 2026 papers target instead.

---

## Closing block (appended by CC, 2026-10-30)

Checks run after commit 3 (`356bba4f`), before this block was written; this block is commit 4, committed alone.

```
$ git diff --name-only c2967362 HEAD | grep -vE '^(dev/|vignettes/articles/(methodology|paper_simulations|treatment_effect_definitions|biomarker_effects|causal_effects_brief_review)\.Rmd$)'
(nothing)
$ git status --porcelain | grep -v '^?? '
(nothing)
$ git status --porcelain | grep '^?? '
?? vignettes/articles/guohe-nested-continuum.html
$ git log --oneline c2967362..HEAD
356bba4f Report: article alignment with 2026a/2026b (applier OK, five renders, output-only gates, check_pkgdown)
a88f8ea4 docs: align methodology, paper_simulations, treatment_effect_definitions, biomarker_effects and causal_effects_brief_review with León & Anderson 2026a/2026b
bdeb75f8 task: article alignment — v2 (gate 12a re-stated for toc_float; Guo–He HTML restore)
794863f3 task: article alignment with 2026a/2026b — task document, edit set, applier
(plus this closing-block commit)
```

Report: `dev/reports/REPORT_articles_alignment_2026-10-30.md`. Copied to `~/Downloads/articles_alignment_2026-10-30/`: the report, `results/`, and the five rendered HTML files.
