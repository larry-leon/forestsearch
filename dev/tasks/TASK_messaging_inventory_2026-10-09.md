# TASK — Inventory of user-facing interpretive text (read-only)

- Date: 2026-10-09
- Repo: `larry-leon/forestsearch`, branch `feature/glm-extension`, current HEAD.
- Why: before the 0.4.0 CRAN submission, Larry wants the interpretive language that reaches users toned down and the methods paper (arXiv; under review at JRSS-B) cited. The `print()`/`summary()` output of an MR fit prints simulation-coverage caveats and points users to `dev/notes/NOTE_survival_products_2026-09-09.md`, which is not shipped; the `survival-post-selection` vignette cites `dev/` and `REPORT_*` files throughout. This task only finds and lists such text; nothing is edited.
- Type: read-only. No edits anywhere except the report. No installs, no checks, no uploads.
- R/ callout: none.

## First action

1. Copy this file to `dev/tasks/TASK_messaging_inventory_2026-10-09.md` and commit it alone.

## Steps (≤ 30 min total)

1. In `R/` (code and roxygen), `vignettes/*.Rmd`, `vignettes/*.qmd` (not `vignettes/articles/`), `README.md`, `inst/`, list every user-facing occurrence of:
   - a. any path or file name under `dev/`, or `REPORT_`, `NOTE_`, `TASK_`, `HANDOFF_`;
   - b. "certified" / "not certified" / "certification" / "source of truth";
   - c. simulation operating-characteristic figures quoted as text (coverage, bias, prevalence %, cell counts, e.g. "0.913-0.917", "12.4% prevalence", "n >= 1000");
   - d. phrases such as "unstable-pick pole", "stable-pick", "exposed one", "over-correction", "under-corrects".
   "User-facing" = text that a user sees: strings passed to `cat`/`message`/`warning`/`stop`/`print`/`format`/`sprintf` that reach output, roxygen text that becomes `man/`, and vignette prose or printed output. Skip comments that never reach output.
2. For each hit record: file:line; the enclosing function (for R/) or section heading (for vignettes); whether it reaches output and by which path (e.g. `print.forestsearch()` → `.print_mr_block()`); and the full text block verbatim (the whole paragraph, not only the matching line).
3. Group hits by the paragraph/message they belong to, so each distinct block appears once with all its line numbers. For each block, state which outcome paths print it (survival / binary / continuous / count) and under which conditions (e.g. only when `mr_inference = TRUE`, only when p-hat < 0.20).
4. Separately list, in `vignettes/forestsearch.Rmd`, the lines for these known issues (verbatim, file:line): the GBSG treatment description ("versus chemotherapy"), the "included in forestsearch package" comment, the Key Findings `exp(fs$grp.consistency$out_sg$result$hr[1])` line, and "95%95%".

## Deliverable

`dev/reports/REPORT_messaging_inventory_2026-10-09.md`: one section per distinct block (verbatim text, all file:lines, output path, outcome paths, conditions), then the `forestsearch.Rmd` list, then a count summary. Commit it, then copy it to `~/Downloads`.

## Post-conditions

- PC1. `git diff --name-only <task-doc-commit>^ HEAD` lists only this task doc and the report.
- PC2. `git status --porcelain` clean.
