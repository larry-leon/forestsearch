# REPORT — primary vignette (2026-10-27)

Task: `dev/tasks/cc_task_vignette_primary_2026-10-27.md`, branch `feature/vignettes-0.4.1`. Run in two sessions: the first stopped at the §2.3 render gate (inline-code parse error at line 167 of the draft); the second resumed at §2.1 with `~/Downloads/forestsearch_v2.qmd` and ran through §6. All gates pass. Not pushed.

## Commits

- `0908a3c1` dev: task — primary vignette (2026-10-27)
- `225cac4d` docs: move the 2024 vignette to the pkgdown articles as forestsearch-2024, with a pointer note
- `c810dbd8` docs: new getting-started vignette — Table 1 of León & Anderson (2026b) end to end
- `0de9a8e5` docs: pkgdown navigation, NEWS and cross-reference for the vignette change
- (this commit) dev: report — primary vignette (2026-10-27)

## Provenance

- `dev/reports/vignette_primary_2026-10-27/provenance.txt`: host pop-os; HEAD at task start `30aee979`; R 4.6.1; installed forestsearch 0.4.0.9000, Built `R 4.6.1; ; 2026-10-27 21:01:28 UTC; unix`; quarto CLI 1.10.18 (`/usr/lib/rstudio/resources/app/bin/quarto/bin/quarto`, not on PATH); pkgdown 2.2.1.
- §0 gates: branch correct; `30aee979` is an ancestor of HEAD; installed `NEWS.md` byte-identical to the repo's; pre-existing untracked/modified paths: none (tree clean).
- Input draft `forestsearch.qmd`: 531 lines, md5 `2353159e4232ad9ba1fe82c0bd36ba32`; content-grep checks all pass (index entry line 1×, `<- forestsearch(` 3×, `declaration_c0` string 1×, each of the six chunk labels 1×).
- Input used for the commit, `forestsearch_v2.qmd`: 530 lines, md5 `9aa23ea0609ff423a100d39813162661` (= committed `vignettes/forestsearch.qmd`).

## §1 — the 2024 vignette moved

- `git mv vignettes/forestsearch.Rmd vignettes/articles/forestsearch-2024.Rmd`; `git diff --cached -M`: rename, similarity 97%; diff confined to the `title:` line, the `%\VignetteIndexEntry{}` line and the inserted pointer blockquote (13 insertions, 2 deletions). The blockquote is the first content after the closing YAML `---`, before the `setup` chunk.
- Render (`rmarkdown::render`, installed package, `RSTUDIO_PANDOC` set): exit 0, wall-clock 19 s; HTML 4,428,091 bytes, left untracked on disk. Log: `render_article_2024.log`.

## §2.2 — edits to the draft (old → new)

- `fs-candidates` chunk. The result `fit$grp.consistency$out_sg$result` has columns `Pcons hr N E g m K M.1 M.2`: `N`, `hr`, `Pcons` exist as named; the definition is split over `M.1` / `M.2` (`M.2` is `""` for the one single-cut candidate), so column 1 alone would drop the second cut.
  - old:
    ```
    # CC (task §2.2): confirm the column names used below against names(cand) — the draft
    # assumes a subgroup-definition column, `N`, `hr` and `Pcons`; adjust and record.
    cand <- cand[order(-cand$N), ]
    knitr::kable(cand[, c(1, match(c("N", "hr", "Pcons"), names(cand)))], digits = 3, row.names = FALSE,
    ```
  - new:
    ```
    cand$candidate <- ifelse(cand$M.2 == "", cand$M.1, paste(cand$M.1, cand$M.2, sep = " & "))
    cand <- cand[order(-cand$N), ]
    knitr::kable(cand[, c("candidate", "N", "hr", "Pcons")], digits = 3, row.names = FALSE,
    ```
- Missing element paths: none. Every path the draft reads on `fit_fs` exists in stage 0's `fit_S3.rds`; the paths absent from `fit_D3` / `fit_G3` (`filter_counts`, `grp.consistency$n_passed`, `declaration_field`) are read only on `fit_fs`.
- Rendering options: none added. The two `progress|\|=+` matches in the HTML are Bootstrap CSS (`progress-bar` rules), not leaked output.
- Line 167 (chat's fix, carried in `forestsearch_v2.qmd`; the only line differing from the first-session working copy):
  - old: ``The identified region is `r sprintf("`%s`", paste(fit_fs$sg.harm, collapse = " & "))`:``
  - new: ``The identified region is `r sprintf("\x60%s\x60", paste(fit_fs$sg.harm, collapse = " & "))`:``

## §2.3 — render

- `quarto render vignettes/forestsearch.qmd` under `/usr/bin/time -v`: exit 0; wall-clock 73.8 s (`real 1m13.847s`); peak RSS 1,108,816 kB; `vignettes/forestsearch.html` 1,277,105 bytes (self-contained; no `forestsearch_files/`). Greps of `Error in`, `Warning`, `NaN`, `NA (` in the HTML: 0 / 0 / 0 / 0. HTML is git-ignored (`.gitignore:15 vignettes/*.html`), not staged. Log: `render_vignette.log` (the first session's failing render log was overwritten by this run).

## §2.4 — identity with stage 0: PASS

`gate_identity.R` (purl + source of the vignette, fresh `Rscript`, sourced in 23.9 s), output `gate_identity_output.txt`. 54 comparisons (18 quantities × 3 fits), each `identical()` against the stage-0 RDS and string-equal against `ladder_quantities.csv` at 6 dp: 54/54 and 54/54 equal.

| quantity | FS (S3) | DINA (D3) | GRF (G3) | equal (rds / csv) |
|---|---|---|---|---|
| region | {er <= 0} & {pgr <= 26} | {grade3 >= 1} & {pgr <= 20} | {er <= 0} | all TRUE / TRUE |
| n_selected / Ĥᶜ n | 75 / 611 | 113 / 573 | 82 / 604 | all TRUE / TRUE |
| events Ĥ / Ĥᶜ | 41 / 258 | 61 / 238 | 45 / 254 | all TRUE / TRUE |
| naive Ĥ (lo, hi) | 2.221839 (1.169437, 4.221320) | 1.314631 (0.745910, 2.316974) | 1.951393 (1.036610, 3.673451) | all TRUE / TRUE |
| naive Ĥᶜ (lo, hi) | 0.605388 (0.465249, 0.787738) | 0.640768 (0.489297, 0.839130) | 0.614995 (0.472315, 0.800778) | all TRUE / TRUE |
| corrected Ĥ / Ĥᶜ | 1.276181 / 0.644653 | 1.076093 / 0.658316 | 1.110994 / 0.663287 | all TRUE / TRUE |
| field lower Ĥ | 0.633561 | 0.637786 | 0.510471 | all TRUE / TRUE |
| field-s upper Ĥᶜ | 0.803570 | 0.842063 | 0.846291 | all TRUE / TRUE |
| Bonferroni pair | 0.508684 / 0.839594 | 0.572467 / 0.886822 | 0.403679 / 0.874401 | all TRUE / TRUE |
| n_family | 1744 | 30 | 858 | all TRUE / TRUE |

The vignette value equals the stage-0 value in every cell; the full per-row table is in `gate_identity_output.txt`. The arguments the draft drops from the transplanted calls are no-ops on all 54 quantities.

## §2.4 — recorded checks (not gated): all as expected

- funnel: 2,211 / 1,410 / 121 / 9 / 1,744.
- candidates table: 9 rows; band `>= 2.03`; selected N 75; max N 82; max-β N 61.
- three-way partition (N/events): FS∧DINA∧GRF 35/19; FS∧GRF only 40/22; DINA only 78/42; GRF only 7/4; none 526/212 — equals Table S15.

## §2.5 — declaration: PASS

| quantity | vignette | Table S17 (GBSG) |
|---|---|---|
| κ̂ at c₀ 0.75 | 2.726 | 2.726 |
| implied level at c₀ 0.75 | 0.994 | 0.994 |
| fw_size at c₀ 0.75 | 0.667 | 0.667 |
| T_H | 2.438 | 2.438 |
| T_max | 2.564 | 2.564 |
| κ̂ at c₂ | 3.503 | 3.503 |
| fw_size at c₂ | 0.974 | 0.974 |
| admits (c₀ 0.70 / 0.75 / 0.80 / 0.85 / c₂) | 0 / 0 / 0 / 0 / 0 | 0 |

Full grid (recorded): c₀ 0.70 κ̂ 2.578, implied 0.990, fw 0.558; 0.80 κ̂ 2.891, 0.996, 0.769; 0.85 κ̂ 3.039, 0.998, 0.844.

## §2.6 — self-containment: PASS

- `grep -nE 'dev/|readRDS|read\.csv|load\(|\.\./' vignettes/forestsearch.qmd`: no output (exit 1), on both the first draft and v2.

## §3 — `check_pkgdown()`: PASS

- Output (`check_pkgdown.txt`): `✔ No problems found.`
- Navbar: `---` separator and "The 2024 workflow (León et al. 2024)" → `articles/forestsearch-2024.html` appended after Extreme Subgroup Effects; `articles:` gains the `Reference` block; `home: description:` replaced (see deviations); "Getting Started" and "Introduction to forestsearch" unchanged.
- `NEWS.md`: two bullets under a new `### Documentation` after `### Bug fixes` in the development section.
- `vignettes/survival-post-selection.qmd`: the one cross-reference replaced as specified; no other change.

## §4 — 121 screened vs 120 evaluated

- The dropped candidate is `{er <= 3} & {pgr <= 77}` (row 72 of `find.grps$out.found$hr.subgroups`, grp 409). It has the same membership as `{er <= 3} & {pgr <= 45}` (grp 407): no patient has er ≤ 3 and 45 < pgr ≤ 77, so both are the same 122 patients with identical n / E / d1 / m1 / m0 / HR / CI. `remove_near_duplicate_subgroups()` drops the second via `!duplicated(dup_key)` on the rounded K…U(HR) columns (`R/subgroup_consistency_helpers.R:1133`). It is called from `R/subgroup_consistency_main.R:618` for every `sg_focus` other than `"maxeff"` (the vignette's `effMaxSG` normalizes to `hrMaxSG`).

## Deviations (verbatim)

- First session: the §2.3 render failed at draft line 167 (`Error: Failed to parse the inline R code`). Task stopped; resumed with chat's `forestsearch_v2.qmd`.
- `provenance.txt` `HEAD:` is the post-task-doc commit `0908a3c1`; the task-start HEAD `30aee979` is recorded on the next line. A line for the v2 input was appended on resume.
- `gate_identity.R`: the first run failed in the script itself (`Error: unexpected 'else'`, a top-level `if … else` split across lines in my `fmt` helper) after the vignette had been sourced. I wrapped the helper in braces and reran; the reported output is from the rerun.
- `NEWS.md`: the subsection is `### Documentation`, not `## Documentation` as the task text gives. Every existing subsection under `# forestsearch (development version)` uses `###`.
- `_pkgdown.yml` `home: description:`: the specified text, on one line, is double-quoted. Unquoted, it is invalid YAML (`Scanner error: mapping values are not allowed in this context at line 23, column 70`, from the `: ` after "clinical trials").
- Commit 2 was made with both the new and the old path (`vignettes/articles/forestsearch-2024.Rmd vignettes/forestsearch.Rmd`) so that the rename's deletion side was included.
- `check_pkgdown.txt` is an extra file in the report directory (the §3.1 gate output).

## OPEN ITEMS (findings only)

- The gate's purl-and-source of the three fits took 23.9 s; the full render took 73.8 s; stage 0's fit times for S3 + D3 + G3 sum to 69.3 s. The vignette text says "about a minute and a half on one core."
- `remove_near_duplicate_subgroups()` keys on summary statistics, so it removes exact-membership duplicates such as the §4 pair. It would also collapse two distinct subgroups whose ten statistics agree to 0.001; the comment at `R/subgroup_consistency_main.R:601-607` makes the same point for `maxeff`.
- Rendered HTML left on disk: `vignettes/forestsearch.html` (git-ignored) and `vignettes/articles/forestsearch-2024.html` (untracked, not ignored).

## §6 — copies in `~/Downloads/vignette_primary_2026-10-27/`


```
total 5768
drwxrwxr-x 2 larryleon larryleon    4096 Oct 27 21:03 .
drwxr-xr-x 6 larryleon larryleon   20480 Oct 27 21:03 ..
-rw-rw-r-- 1 larryleon larryleon 4428091 Oct 27 21:03 forestsearch-2024.html
-rw-rw-r-- 1 larryleon larryleon 1277105 Oct 27 21:03 forestsearch.html
-rw-r--r-- 1 larryleon larryleon   27783 Oct 27 21:03 forestsearch.qmd
-rw-rw-r-- 1 larryleon larryleon    8762 Oct 27 21:03 gate_identity_output.txt
-rw-r--r-- 1 larryleon larryleon   92420 Oct 27 21:03 NEWS.md
-rw-r--r-- 1 larryleon larryleon   18738 Oct 27 21:03 _pkgdown.yml
-rw-rw-r-- 1 larryleon larryleon   10293 Oct 27 21:03 REPORT_vignette_primary_2026-10-27.md
```

(The report's size above is before this listing was appended; the copy is refreshed after.)
