# CC task — read-only inventory of `quarto/applications/` ahead of the 0.3.1 re-run

**Date:** 2026-08-31 · **Repo:** `~/Documents/GitHub/forestsearch` · **Branch:** `feature/glm-extension` · **Installed package expected:** forestsearch 0.3.1
**Type:** read-only inventory. **No renders. No cache files touched. No `R/` changes. Nothing modified or deleted anywhere.** The only new files from this task are this document (copied to `dev/tasks/` and committed on arrival, per kickoff) and the report below.
**Report:** `dev/tasks/REPORT_applications_inventory_2026-08-31.md`
**Compute:** negligible (grep, `git log`, `readRDS`); no go/no-go applies.

## Context — complete; do not consult notes outside the repository

The applications under `quarto/applications/` (gbsg, actg175) are to be re-rendered under 0.3.1. Before any render is specced, Larry needs a from-source inventory: what each application document actually passes to `forestsearch()`, where payloads and leave-one-out (LOO) caches sit on disk, what the 08-27 reproduction check against 0.2.0 actually did (documents, method, baseline, cache handling, wall-clock), and what a re-render would cost. Two facts drive the design: (i) since 0.3.0 the default `vi.grf.min` is `NULL`, and the applications are believed to rely on the default, so an as-is render takes the `NULL` path for the first time; (ii) the LOO cache key (`cv_out_<doc>_<focus>_<rule>.rds`) carries no package version, so a stale cache would silently report 0.2.x LOO results under a 0.3.1 render. Expectations from the project record are stated below only so deviations can be flagged; every claim in the report must trace to repository source or command output gathered in this task.

## Gates — stop on failure

On any gate failure: write the report with everything gathered so far, put `STOPPED AT GATE <n>: <reason>` as the first line, commit it, and end. Do not ask questions.

- **Gate 1 — repo state.** `git rev-parse --abbrev-ref HEAD` must be `feature/glm-extension`. Record `git log -1 --oneline` and `git status --porcelain` (a dirty tree is a finding to record, not a failure).
- **Gate 2 — targets exist.** `quarto/applications/` exists and contains at least one `.qmd`.
- **Gate 3 — the 08-27 record is locatable.** Locate `REPORT_repro_check_vs_0.2.0_2026-08-27.md` and `HANDOFF_continuous_2026-08-27_v5.md` (expected under `dev/tasks/` and/or `dev/glm-continuous-sims/`). If not at the expected names, run `find dev -iname '*repro*' -o -iname '*2026-08-27*'` before declaring failure. Fail only if neither file can be found.

## Steps

### 0. Environment
Record in the report header: hostname; branch; `git log -1 --oneline`; `git status --porcelain | wc -l` (list the lines if ≤ 20); `Rscript -e 'cat(as.character(packageVersion("forestsearch")), "|", R.version.string, "\n")'`.

### 1. Document census
- `git ls-files 'quarto/applications/**'` filtered to `.qmd`, and `find quarto/applications -name '*.qmd' | sort`. Report both; flag untracked files and anything under archive-like paths.
- For every `.qmd`: `git log -1 --format='%h %ad %s' --date=short -- <path>`.

### 2. Per-document extraction — every `.qmd` from step 1
Read each document's source and record:
- Whether it calls `forestsearch()`. If yes, **quote every `forestsearch()` call verbatim** in a fenced block with line numbers. These quotes are the primary record; the summary-table fields are read off them, never paraphrased from memory.
- `vi.grf.min`: explicit value, or **absent** (⇒ takes the 0.3.x default `NULL` on a re-render).
- `sg_focus`; selection rule under either spelling (`selection_rule` / `fs_selection_rule`); `consistency_method`; `parallel_args`; seed handling (`set.seed(...)` and any seed-bearing arguments); any occurrence of `max_subgroups_search` (expected: provenance assertions and prose only — quote every hit).
- Bootstrap section: present? Quote its gating flag/condition and the replicates argument as written.
- MR section: present? Quote its gating condition and the number of draws as written.
- LOO section: present? Quote its gating condition and the cache path construction as written (`gh_dir`, `cv_out_*` naming).
- Payload write: quote the path construction as written, and whether `results_dir` / `dirout` override parameters exist in the document (quote the params lines). Expected per the record: `<qmd dir>/_payloads/<stem>/<stem>_payload.rds` with overrides defaulting `NULL` — flag deviations in §9.
- Documents that do **not** call `forestsearch()` (summaries, cross-analysis): quote their input lines (`readRDS` etc.) instead.

### 3. Payloads on disk
- `find quarto/applications -type d -name '_payloads*'` and `find quarto/applications -name '*_payload.rds' -printf '%p\t%TY-%Tm-%Td %TH:%TM\t%s bytes\n'`. Also note any version-tagged payload directories already present.
- For every payload found:
  `Rscript -e 'p <- readRDS("<path>"); cat(as.character(p$forestsearch_version), "|", format(p$built_at), "|", p$est_scale, "|", paste(names(p), collapse=","), "|", paste(names(p$extras), collapse=","), "\n")'`
- Deep dive on exactly two payloads — one gbsg (expected `est_scale = "hr"`), one actg175 binary (expected `"or"`): `str(p, max.level = 2)` into a fenced block, and **name the field(s) where the selected subgroup's definition lives** (needed later for a membership-based comparison). If any other payload's top-level names differ from these two, add one `str()` per differing shape. Read-only throughout.

### 4. LOO caches on disk
- Resolve each document's `gh_dir` from the step 2 quotes; list `cv_out_*.rds` there with mtime and size. Also run repo-wide catch-all: `find . -name 'cv_out_*.rds' -not -path './.git/*' -printf '%p\t%TY-%Tm-%Td %TH:%TM\t%s bytes\n'`.
- Table: cache file | mtime | size | parsed `<doc>_<focus>_<rule>` key | which document(s) from step 2 would consume it on a re-render.

### 5. The 08-27 record, quoted
From the located report, its task document (if present under `dev/tasks/`), and `HANDOFF_continuous_2026-08-27_v5.md`, quote verbatim — with file path and line numbers — the passages stating:
1. the exact list of documents rendered; resolve "and their siblings" to names if any source names them;
2. the package version and commit rendered under;
3. the comparison method — what was compared (payload fields, memberships, definition strings, rendered output);
4. how LOO caches were handled;
5. the baseline: which 0.2.0 payloads, produced on which machine at which commit;
6. wall-clock: the ≈19-minute total and any per-document timings.

Where a source is silent on an item, write "not stated in `<file>`". Do not infer.

### 6. Crossanalysis state
- In the crossanalysis document (expected `actg175/analysis_actg175_crossanalysis_summary.qmd`; locate by census if moved): quote the lines reading `selected_subgroups_binary.rds` and `comparison_binary.rds` and their `file.exists()` guards.
- Repo-wide: `find . \( -name 'selected_subgroups_binary.rds' -o -name 'comparison_binary.rds' \) -not -path './.git/*'`, and `grep -rn` tracked files for any writer of those filenames.

### 7. Render-cost estimates — no timing runs
Build a per-document estimate table using, in order of preference: (i) per-document timings quoted in step 5; else (ii) the 19-minute total apportioned by structural weight (bootstrap B, MR draws, LOO on/off), labeled as an apportionment; and separately (iii) the LOO-recompute cost as N × single-search time, with N confirmed from each document's data-loading code (record: gbsg N = 686 — confirm) and single-search time taken from a step 5 source if stated, else marked "unknown — not measured in this task". **Label every number with its basis.** Run nothing to obtain a timing.

### 8. Report and commit
Assemble `dev/tasks/REPORT_applications_inventory_2026-08-31.md` using the skeleton below. `git add` the report and commit with message `applications: read-only inventory ahead of 0.3.1 re-run (report)`. Do not push. Print the report's path and stop.

## Report skeleton — use exactly these sections

```
# REPORT — applications inventory (read-only) — 2026-08-31
[STOPPED banner here if any gate failed]
## 0 Environment
## 1 Document census
## 2 Per-document extraction   (2.1, 2.2, … one subsection per document; verbatim call quotes first, then the field summary)
## 3 Summary table   (doc | forestsearch()? | vi.grf.min | sg_focus | rule | consistency | parallel_args | boot B | MR draws | LOO | results_dir override | payload path)
## 4 Payloads on disk   (per-payload version/built_at/est_scale table; the two str() deep dives; where the sg definition lives)
## 5 LOO caches
## 6 The 08-27 record, quoted   (6.1 documents · 6.2 version+commit · 6.3 method · 6.4 cache handling · 6.5 baseline provenance · 6.6 wall-clock)
## 7 Crossanalysis state
## 8 Render-cost estimates   (every number labeled with its basis)
## 9 Deviations from the record's expectations   (payload path/shape, cache key form, vi.grf.min reliance, max_subgroups_search occurrences, est_scale values — anything differing from the Context above)
```

## Hard rules
- Read-only: no renders; no modifications outside `dev/tasks/`; no deletions; no cache files created, moved, or removed; nothing under `R/`.
- Verify from source: every claim in the report traces to a quoted line (with path and line number) or a command output captured in this task.
- No follow-up work on findings: findings land in §9 and stop there.
