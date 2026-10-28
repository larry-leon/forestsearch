# Inventory — post-selection to pkgdown article (2026-10-28)

Start time: 2026-10-28T13:31:58-07:00. `future::availableCores()` = 128.
HEAD at inventory: `0c07d261` (task document), parent `b2e43397`.

## 1. References

Command: `git grep -n -i -e 'survival-post-selection' -e 'post-selection vignette' -e 'companion vignette' -- ':(exclude)*.html'`

### Class *edit* (4 hits, all prose — blockquotes or bullet)

| File:line | Context | Sentence |
|---|---|---|
| `vignettes/articles/forestsearch-2024.Rmd:29` | prose (blockquote) | "The current methodology — the search as closed-form stability selection with its declaration calibration (León and Anderson 2026b, arXiv:2610.28379) and the selection-adjusted inference reported with `mr_inference = TRUE` (León and Anderson 2026a, arXiv:2609.38361) — is documented in `vignette("forestsearch")` and `vignette("survival-post-selection")`." |
| `vignettes/articles/methodology.Rmd:64` | prose (blockquote) | "The current methodology is documented in `vignette("forestsearch")` and `vignette("survival-post-selection")`; the argument listings below are kept for the 2024 workflow and are not the complete current signature of `forestsearch()` — see `?forestsearch`." |
| `vignettes/articles/paper_simulations.Rmd:59` | prose (blockquote) | "The operating characteristics of the current methodology — … — are reported in León and Anderson (2026a, arXiv:2609.38361) and (2026b, arXiv:2610.28379); the methodology is documented in `vignette("forestsearch")` and `vignette("survival-post-selection")`." |
| `vignettes/forestsearch.qmd:611` | prose (bullet under "Where to go next") | "- `vignette("survival-post-selection")` reads the bounds in more depth on the same trial: the re-selection frequency p̂ as a diagnostic, membership agreement, and what each construction does and does not say." |

No hits in `README.Rmd` or `README.md`.

### Class *pkgdown* (2 hits)

- `_pkgdown.yml:45` — navbar `href: articles/survival-post-selection.html` (stays).
- `_pkgdown.yml:546` — articles index `- survival-post-selection` (to become `articles/survival-post-selection`).

### Class *news* (1 hit; development section is lines 1–76)

- `NEWS.md:71` — bullet lines 71–75: "`vignette("survival-post-selection")` now reproduces the GBSG analysis of León and Anderson (2026a): the three identifiers at the manuscript's settings, multiplier resampling in full and the full bootstrap at a small B for illustration, with `summarize_bootstrap_results()` and `mr_estimates_table()`. The methodology and simulation articles on the package website carry a note placing them as the 2024 workflow, with the deprecated arguments removed from their listings."

### Class *record-only*

- `vignettes/survival-post-selection.qmd:113` (prose) — "The companion vignette uses …"; `:127` (prose) — "… companion vignette and the membership-recovery diagnostics used below." (the article itself; text unchanged by rule.)
- `R/forestsearch_main.R:1435` (roxygen `@examples`, `\dontrun{}` comment) — "## León and Anderson (2026a, Sec. 5) and the survival-post-selection ## vignette."
- `man/forestsearch.Rd:1309` — generated copy of the same comment.
- `HANDOFF_vignettes_release_0.4.1_2026-10-27.md`: lines 58, 74, 77, 81, 89, 94, 102.
- `HANDOFF_vignettes_release_0.4.1_2026-10-28.md`: lines 25, 35, 38, 39, 43, 61, 74, 99, 163.
- `quarto/` (all "companion vignette(s)" referring to other documents, plus two gbsg_020 reports naming `vignettes/survival-post-selection.qmd`):
  `quarto/extreme_subgroups/gbsg_design_comparison.qmd:37,1251`; `quarto/extreme_subgroups/gbsg_survival_fixed.qmd:57,515,767`; `quarto/extreme_subgroups/gbsg_survival_resample.qmd:57,456,692`; `quarto/guides/extreme_subgroups.qmd:384,620`; `quarto/guides/forestsearch_diagnostic_scoping.qmd:52`; `quarto/presentations/my use of claude/extreme_subgroups.qmd:384,620`; `quarto/simulations/README.qmd:456`; `quarto/simulations/actg175/README.qmd:456`; `quarto/simulations/gbsg_020/REPORT_field_recovery_2026-09-09.md:446,484`; `quarto/simulations/gbsg_020/REPORT_print_vignette_2026-09-09.md:1,227`.
- `dev/` — hit counts per file:

```
      1 dev/glm-continuous-sims/REPORT_bootstrap_stage0_2026-10-26.md
      1 dev/glm-continuous-sims/REPORT_bootstrap_stage0b_2026-10-26.md
      1 dev/glm-continuous-sims/stage0b_glm_2026-10-26/logs/b0_provenance.log
      1 dev/glm-continuous-sims/stage0b_glm_2026-10-26/logs/b3_summary.log
      1 dev/glm-continuous-sims/stage0b_glm_2026-10-26/results/keyvalues.txt
      1 dev/glm-continuous-sims/stage0b_glm_2026-10-26/results/stage0b_summary.md
      1 dev/glm-continuous-sims/stage0b_glm_2026-10-26/stage0b_common.R
      1 dev/glm-continuous-sims/stage0_bootstrap_2026-10-26/logs/s1_profile.log
      1 dev/glm-continuous-sims/stage0_bootstrap_2026-10-26/logs/s4_summary.log
      1 dev/glm-continuous-sims/stage0_bootstrap_2026-10-26/results/keyvalues.txt
      1 dev/glm-continuous-sims/stage0_bootstrap_2026-10-26/results/stage0_summary.md
      2 dev/glm-continuous-sims/stage0_bootstrap_2026-10-26/stage0_common.R
      1 dev/reports/postsel_stage0_2026-10-28/ladder.R
      1 dev/reports/postsel_stage0_2026-10-28/postprocess.R
      1 dev/reports/postsel_stage0_2026-10-28/run_ladder.sh
      1 dev/reports/postsel_vignette_2026-10-28/gate_identity_output.txt
      4 dev/reports/postsel_vignette_2026-10-28/gate_identity.R
      5 dev/reports/postsel_vignette_2026-10-28/render_postsel_caption.log
      5 dev/reports/postsel_vignette_2026-10-28/render_postsel.log
      5 dev/reports/postsel_vignette_2026-10-28/render_postsel_math.log
      3 dev/reports/REPORT_cran_release_audit_2026-10-09.md
      6 dev/reports/REPORT_docs_refs_maxminutes_2026-10-27.md
      1 dev/reports/REPORT_extreme_subgroups_inventory.md
      4 dev/reports/REPORT_messaging_inventory_2026-10-09.md
      5 dev/reports/REPORT_messaging_revision_2026-10-09.md
      3 dev/reports/REPORT_postsel_stage0_2026-10-28.md
     17 dev/reports/REPORT_postsel_vignette_2026-10-28.md
      2 dev/reports/REPORT_rerender_compare_2026-10-26.md
      1 dev/reports/REPORT_vignette_primary_2026-10-27.md
      4 dev/reports/REPORT_vignettes_align_2026-10-28.md
      2 dev/reports/REPORT_vignettes_stage0_2026-10-27.md
      3 dev/reports/rerender_2026-10-26/inventory.txt
      4 dev/reports/rerender_2026-10-26/render_V1.log
      1 dev/reports/vi_deprecation_2026-10-26/vi_common.R
      5 dev/reports/vignette_render_2026-10-28/render_post_selection.log
      1 dev/reports/vignettes_align_2026-10-28/logs/render_methodology.log
      5 dev/reports/vignettes_align_2026-10-28/logs/render_survival-post-selection.log
      2 dev/reports/vignettes_align_2026-10-28/postsel_full_precision.R
      2 dev/reports/vignettes_stage0_2026-10-27/package_facts.md
      3 dev/tasks/cc_task_bootstrap_stage0_2026-10-26.md
      1 dev/tasks/cc_task_bootstrap_stage0b_2026-10-26.md
      1 dev/tasks/cc_task_docs_refs_maxminutes_2026-10-27.md
     11 dev/tasks/cc_task_postsel_caption_2026-10-28.md
      5 dev/tasks/cc_task_postsel_stage0_2026-10-28.md
     24 dev/tasks/cc_task_postsel_to_article_2026-10-28.md
     15 dev/tasks/cc_task_postsel_vignette_2026-10-28.md
      2 dev/tasks/cc_task_rerender_compare_2026-10-26.md
      5 dev/tasks/cc_task_vignette_primary_2026-10-27.md
      3 dev/tasks/cc_task_vignette_render_2026-10-28.md
     13 dev/tasks/cc_task_vignettes_align_2026-10-28.md
      7 dev/tasks/cc_task_vignettes_math_2026-10-28_v2.md
      2 dev/tasks/cc_task_vignettes_stage0_2026-10-27.md
      1 dev/tasks/TASK_messaging_inventory_2026-10-09.md
      4 dev/tasks/TASK_messaging_revision_2026-10-09.md
      3 dev/tasks/TASK_print_vignette_2026-09-09.md
```

- Tracked `.html` containing `survival-post-selection`: `vignettes/articles/forestsearch-2024.html`.

## 2. Relative paths in `vignettes/survival-post-selection.qmd`

YAML header (lines 1–18): no `bibliography`, `csl`, `css`, `include-*` or other file keys (`format: html` with `embed-resources: true`; the CSS is an inline `{=html}` `<style>` block). No markdown images or links to local files; no `include_graphics()`, `readRDS()`, `source()`, `load()`, `file.path()` or `{{< include >}}`. **No path needs rewriting** from `vignettes/articles/`.

## 3. `n_vig`

Whole-word, case-insensitive `vignette` in prose (YAML, fenced chunks and inline code removed): **6**, on lines 70 ("This vignette"), 78 ("this vignette"), 113 ("The companion vignette"), 117 ("that vignette's"), 127 ("companion vignette"), 253 ("This vignette uses at most two"). n_vig ≥ 1 → the note gets the appended sentence.

## 4. Rendered leftovers

- `vignettes/survival-post-selection.html` — `.gitignore:15:vignettes/*.html	vignettes/survival-post-selection.html`

## 5. `.Rbuildignore` coverage

- `vignettes/articles` — excluded by `^vignettes/articles$`
- `HANDOFF_vignettes_release_0.4.1_2026-10-28.md` — **not excluded**
- `dev` — excluded by `^dev$`
