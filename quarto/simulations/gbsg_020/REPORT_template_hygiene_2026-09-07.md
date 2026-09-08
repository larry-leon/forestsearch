# RECORD — Template hygiene (Part C): pooled-meta knob fields and the campaign-tag guard

**Task:** `dev/tasks/TASK_complement_variance_banddial_2026-09-07.md` (4b245516), Part C. Decision T-3 at default (accept underscores). Document-level, add-only edits to `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`; no `R/` change; no committed bundle rewritten.
**Date:** 2026-09-07. Executor: Claude Code, unattended.

## C1 — Pooled-meta assembly

**Finding on the tree before the edit.** The combine-mode poolability key vector already gated all four fields — `"harm_z1_quantile"`, `"sg_focus"`, `"effect_neighborhood"`, `"er_jcuts"` (nb20 Stage 1 added the last two) — and the pooled `meta = list(...)` already carried `sg_focus` and `harm_z1_quantile`, so the omission was the two knob fields in the pooled meta only (nb20 Gate 2 note 1; REPORT_nb20 side issue 1). **Edit:** two lines added to the pooled meta list, `effect_neighborhood = bundles[[1]]$meta$effect_neighborhood %||% NA_real_` and `er_jcuts = bundles[[1]]$meta$er_jcuts %||% NA_integer_`, carried from the first batch as `seed_base`, `sg_focus` and `harm_z1_quantile` are (the gate above has already required every batch to agree). Nothing else in the combine path changed.

**Identity (re-combining the committed nb20 batch bundles).** Two checks:

1. *The template's own combine code, verbatim.* `knitr::purl()` of the edited template, cut before the first report chunk (`mr-settings-readout`), with one line redirected — `combined_path` to a scratch path (`…_combined_recheck.rds`) — because the template's `.refuse_if_tracked()` guard forbids writing over a committed pooled bundle and `FS_S7_SAVE_COMBINED=FALSE` writes nothing to compare. Run under the nb20 env (`FS_S7_MODE=combine`, `FS_S7_FOCUS=effMaxSG`, `FS_S7_Z1Q=0.60`, `FS_S7_NBHD=0.20`, `FS_S7_FIELD_COMPLEMENT=TRUE`, `FS_S7_IJ_RESIDUAL=two_term`, `FS_S7_FB=none`; arm B adds `FS_S7_ER_JCUTS=20`) for each of the seven cells:

| Cell | `results` identical to the committed pooled bundle | `truth` identical | pooled meta fields (was → now) | carries `effect_neighborhood` / `er_jcuts` / `harm_z1_quantile` / `sg_focus` | = batch metas | other shared meta fields identical |
|---|---|---|---|---|---|---|
| A HR 1.50 n500 | TRUE (2000 × 136) | TRUE | 24 → 26 | 0.20 / 10 / 0.60 / effMaxSG | yes | TRUE |
| A HR 1.75 n500 | TRUE (2000 × 136) | TRUE | 24 → 26 | 0.20 / 10 / 0.60 / effMaxSG | yes | TRUE |
| B HR 1.00 n500 | TRUE (2000 × 136) | TRUE | 24 → 26 | 0.20 / 20 / 0.60 / effMaxSG | yes | TRUE |
| B HR 1.50 n500 | TRUE (2000 × 136) | TRUE | 24 → 26 | 0.20 / 20 / 0.60 / effMaxSG | yes | TRUE |
| B HR 1.75 n500 | TRUE (2000 × 136) | TRUE | 24 → 26 | 0.20 / 20 / 0.60 / effMaxSG | yes | TRUE |
| B HR 1.50 n1000 | TRUE (2000 × 136) | TRUE | 24 → 26 | 0.20 / 20 / 0.60 / effMaxSG | yes | TRUE |
| B HR 1.00 n1000 | TRUE (2000 × 136) | TRUE | 24 → 26 | 0.20 / 20 / 0.60 / effMaxSG | yes | TRUE |

`identical()` on the pooled `results` data frame (row order, all 136 columns) and on `truth`; "other shared meta fields" = every field of the committed pooled meta except `n_workers` and `built_at`.

2. *The full document in combine mode with `FS_S7_SAVE_COMBINED=FALSE`* (arm A HR 1.50 n500, output to the session scratchpad): renders clean (rc 0, 19 s), knob echo `focus=effMaxSG nbhd=0.20 er_jcuts=10 … campaign=p30sgnb20 run_mode=combine`, "Combined 2 batch file(s) -> 2000 rows, sim_id 1-2000", every table and the re-selection callout rendered; **no pooled bundle written** (no "Saved pooled bundle" line; `results/` unchanged, `git status` clean of new bundles).

**The committed nb20 pooled bundles are not rewritten** (read-only by convention); they carry the two knobs in their batch metas, which `summary_complement_variance.qmd` and the next campaign's Gate 2 read through `meta$source_files`. Bundles pooled from here on carry all four fields in the pooled meta.

## C2 — Campaign-tag guard

**Edit:** `stopifnot(grepl("^[A-Za-z0-9]+$", campaign_tag))` → `stopifnot(grepl("^[A-Za-z0-9_]+$", campaign_tag))`, with the reason recorded in the template comment: the tag is the stem's final token and `combine_glob` is `<stem>_res_*`, so an underscore inside the tag cannot make one campaign's glob match another's. Checked: `banddial`, `p30sg_nb20`, `p30sgnb20` pass; `bad-tag`, `a b`, `""` fail. Every committed tag is alphanumeric and every committed stem is unchanged; the nb20 record's tags (`p30sgnb20`, `p30sgnb20j20`) stand as written.

## Diff

One file, +14 / −1: the two-line meta addition with its comment (C1) and the one-line guard with its comment (C2). Scripts: `recombine/template_to_combine.R`, `recombine/run_identity.sh` (session scratchpad).

Nothing blocked. Side issue (not fixed): `combined_path` is a literal `NULL` in the template rather than an `FS_S7_*` knob, so a verification combine that wants to *keep* its pooled bundle must edit the document (as the purl here did); an `FS_S7_COMBINED_PATH` knob would make the identity check a pure env-driven render. Not part of this task.
