# TASK — O-1: make the DINA and GRF branches operational for the current product set (argument forwarding, add-only) and verify what a campaign needs

Date: 2026-09-10. Author: chat (spec). Executor: Claude Code (Linux). Approver: Larry (O-1 authorized 2026-09-10). Reviewer: the Linux MR-field chat.
Predecessors: `REPORT_dinamr_stage0_2026-09-10.md` (the STOP this resolves), `dev/notes/NOTE_survival_products_2026-09-09.md`, `REPORT_cimethod_flip_2026-09-09.md` (the gate pattern), manuscript §4 (alignment), §5 (the fixed-family condition), Supplementary S1/S2.

**Why.** `.fs_apply_mr()` forwards 15 of `fs_mr_inference()`'s 25 arguments. Two field products (`field_decompose`, `field_recovery`) are therefore **unreachable** on the DINA and GRF branches, and three more (`field_complement`, `field_scale_complement`, `ij_residual`) are **inert while their defaults coincide with the intended values** — so a campaign's meta would record them as set while they controlled nothing, and flipping one to test the alternative would silently do nothing. That is a provenance defect on two engines, worth closing independently of any campaign.

**What this task does NOT do.** It does **not** change what any branch defaults to. Forwarding uses `fs_mr_inference()`'s own defaults, so every existing DINA and GRF call stays byte-identical. Whether DINA/GRF should *default* to the field constructions is a separate decision and is not taken here.

**Scope note on the criterion.** Verified from source: DINA's `effMaxSG` reuses the inclusion-band logic shared with `forestsearch()`, honors `effect_neighborhood`, and its semantics match. The FS-analogous criterion for these engines is therefore **`effMaxSG` with the effect floor and no consistency floor** — not `maxeff`. Item V1 below confirms the same for GRF.

## Protocol

- First action: archive stale variants in `~/Downloads` to `~/Downloads/cc_archive/` — **do not archive `HANDOFF_guohe_comparison_2026-09-09.md`**. Copy this file to `dev/tasks/` and commit. **Commit only; do not push.**
- **One `R/` change, confined to `R/fs_mr_inference_methods.R`**, classified **adds code; byte-identical defaults**. No other `R/` file. No template change in this task.
- Compute: verification renders and pilot fits only (≤ 20 replicates). No campaign.
- Gates stop on failure; **on failure revert the touched file to HEAD, re-install, record the failure, and stop.**
- Standing conventions: verify from source; transplant-first; winner-only and winner-floor excluded. Leave the seven pre-existing untracked files alone.

## Part F — The forwarding fix

**F1.** In `.fs_apply_mr()`, forward the ten arguments it currently drops — `return_reselection`, `field_R_out`, `field_R_in`, `field_uniform`, `field_M_cap`, `field_complement`, `field_decompose`, `field_scale_complement`, `ij_residual`, `field_recovery` — each through the existing `.g(mr_inference_args$<name>, <default>)` idiom, **with `<default>` taken from `fs_mr_inference()`'s own formals rather than re-stated by hand**. Quote both the before and after argument lists in the record, and quote the formals you read the defaults from, so the equality is verifiable rather than asserted.

**F2.** Roxygen note on `.fs_apply_mr()` (or its enclosing documented function): the wrapper now forwards the full MR argument set, so the DINA and GRF branches are controllable through `mr_inference_args`; defaults are unchanged. `devtools::document()` if any exported documentation changes; `NEWS.md` one bullet under the development header.

**F3.** `devtools::install(dependencies = FALSE)`; `deparse()` of the installed `.fs_apply_mr()` equals source; read the installed wrapper's call node back and quote the resolved defaults.

**Gate F** (DINA and GRF both; M1 DGM, HR 1.75, n = 500, seeds `8316951 + sim_id`, sim_id 1–5):

- **Fa — byte-identical when nothing is asked for.** With `mr_inference_args` **omitted entirely**, the returned MR object on each engine is `identical()` to the same call on the pre-change installed package (build both captures explicitly; exclude timing elements). **This is the criterion that makes the change add-only.**
- **Fb — byte-identical against committed work.** Re-run one of the seven committed DINA/M1 drivers (HR 1.00, n 500) at its own settings and confirm its recorded columns are `identical()` to the committed bundle, timing excluded. If the driver cannot be re-run as committed, say why and substitute the closest check you can defend.
- **Fc — now controllable.** With `field_decompose = TRUE, field_recovery = TRUE` passed through `mr_inference_args`, `field$complement$decomp_fields` and `field$recovery` are **PRESENT** on both engines and their contents satisfy the standing invariants (`0 ≤ sens_H, ppv_H ≤ 1`, `scale_ratio_c > 0`). Report the present/absent table for both engines, before and after, in the format of the Stage 0 report §3.3.
- **Fd — the inert three are now live.** Setting `field_complement = FALSE` and `field_scale_complement = "none"` explicitly now **changes the output** on both engines (the complement and `_s` blocks absent / present as asked) — proving they control rather than coincide. This is the provenance half of the fix.
- **Fe — the consistency branch is untouched.** One consistency-engine fit at the standing identity cell is `identical()` to the committed `e1stud` rows 1–5 on all pre-existing non-timing columns.

Full test suite and `devtools::check()` re-run; report both tallies against the previous task's.

## Part V — Verification for a future campaign (no code)

**V1. GRF's band.** Does GRF's `effMaxSG` honor `effect_neighborhood` and share `forestsearch()`'s inclusion-band logic, as DINA's does? Quote the source. **If it does not, say so plainly** — the FS-analogous criterion would then be unavailable on GRF and the campaign design would have to change.
**V2. Admission per engine.** Extend the Stage 0 report's `.fs_resolve_admission` table to GRF, and to `effMaxSG` on both engines at ε = 0.20, so the effect floor and the absence of a consistency term are on record for the exact campaign setting.
**V3. Stem collision.** Under a DINA (and GRF) transplant of the field template at `sg_focus = "effMaxSG"`, derive the output stem and confirm it cannot collide with any committed pre-field bundle. Quote both stems.
**V4. Recorder meaning.** List, measured on a pilot fit per engine, which of the 158 recorder columns are populated, which are structurally `NA` without a consistency screen (`n_cons_qual`, `band_n`, `p_star`, and any others), and which field/`_s`/joint/p̂/recovery columns are populated **after** the fix.
**V5. Engine knobs.** State whether the template's existing `dina_select_statistic`, `grf_select_statistic`, `grf_selection`, `grf_depth`, `dmin.grf` defaults are the campaign-appropriate values (alignment requires ranking on the inferential effect — confirm both `*_select_statistic` resolve to `"effect"`), or name what a campaign must set.
**V6. Cost.** A 20-replicate pilot per engine at `effMaxSG` ε 0.20, HR 1.75, n = 500: detection rate, proposed-family size distribution (mean, median, 10th/90th percentile, min, max), and per-replicate wall, so a future Gate 1 can project from the family-size distribution rather than from FS walls. **No projection is offered here** — the data are the deliverable.

**Output:** `REPORT_o1_forwarding_2026-09-10.md` beside the Stage 0 report, with Part F's gate results and Part V's six answers, every number verbatim.

## Done means

Part F committed on PASS with Fa–Fe concrete values (or the failure recorded and the file reverted); Part V's six items answered with source quotes and measured values; test suite and `check()` tallies reported; one-paragraph closing summary naming what a DINA/GRF campaign still needs and what it no longer needs. **Out of scope:** any default change on any branch, any template change, any campaign, the prevalence and ε decisions (Larry's).
