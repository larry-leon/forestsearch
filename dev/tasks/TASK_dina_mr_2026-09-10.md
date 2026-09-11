# TASK — MR campaign for the DINA interpretable identifier (campaign `dinamr`): the current product set on a model-generated family

Date: 2026-09-10. Author: chat (spec). Executor: Claude Code (Linux). Approver: Larry (2026-09-10). Reviewer: the Linux MR-field chat.
Predecessors: `REPORT_cert20_2026-09-08.md` and `REPORT_tier2_2026-09-08.md` (the FS certification this is analogous to), `dev/notes/NOTE_survival_products_2026-09-09.md`, `REVIEW_certification_2026-09-09.md`, manuscript Section 4 (alignment), Section 5 (the fixed-family condition), Supplementary S2 (the DINA interpretable construction) and S8.3 with Table S5 (the existing DINA evidence).

## Framing — to be stated in the record's opening paragraph, not buried

- **This campaign does not certify DINA and must not be described as doing so.** The manuscript's fixed-family condition (§2.1) requires a family that "depends on no outcome, fitted effect, or learned surface." DINA's candidates are read off a cross-fit surface that a genuine bootstrap would regenerate, so — after the alignment repair (re-ranking qualifiers on the inferential coefficient β̂(g)) — the de-biased estimate and interval are "first-order exact, **conditional on the proposed family**," with the omitted family-generation component as the gap to the unconditional target.
- **The target is therefore the conditional-on-proposed-family estimand**, and every coverage number in the record is coverage of it. State this beside each coverage table.
- **What is new.** The existing DINA evidence (Supplementary Table S5) predates the current product set: it carries no field-s complement bound, no Bonferroni joint pair, no p̂, and no recovery diagnostics. Those have **never** been measured on a model-generated family.
- **Correction of a common premise, to be recorded:** DINA is not `maxeff`-only. Per `forestsearch_main.R` (the per-engine `sg_focus` table), for DINA and GRF `maxeff` and `maxeffCons` are **synonyms** — neither engine has a consistency floor — while `effMaxSG`, `effMinSG`, `maxSG` and `minSG` remain available. Quote that table at Stage 0.

## Protocol

- First action: archive stale variants in `~/Downloads` to `~/Downloads/cc_archive/` — **do not archive `HANDOFF_guohe_comparison_2026-09-09.md`**. Copy this file to `dev/tasks/` and commit. **Commit only; do not push.**
- **No `R/` change.** If the campaign appears to need one, STOP and report.
- Gates stop per cell; `.refuse_if_tracked()` live; `devtools::install(dependencies = FALSE)` if installed does not match HEAD; never `load_all()`. Leave the seven pre-existing untracked files alone.
- Standing conventions: winner-only and winner-floor excluded; bounds read **by location**; Wilson intervals; marginal and error SDs side by side; NPV beside sensitivity/specificity/PPV; **no recommendation change, no repair proposal; report and wait.**

## Stage 0 — Discovery, with a hard STOP

1. **Does a DINA simulation document exist?** Search `quarto/simulations/` and the repo for any committed document that drives `subgroup_method = "dina"` under the M1 DGM. Report what exists.
   - **If one exists:** name it, quote its knobs, and use it. Proceed.
   - **If none exists:** **STOP.** Do not author one. Report (a) which committed document is the closest transplant source — expected to be `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`; (b) the exact list of changes a transplant would require (the `forestsearch()` call's `subgroup_method`, `dina_*` arguments, which `FS_S7_*` knobs become inert, which recorder columns lose meaning without a consistency screen — e.g. `n_cons_qual`, `p_star`); (c) whether the template's recorder and the `.fs_apply_mr()` path record the field/field-s/joint/p̂/recovery columns on the DINA branch **at all**. That report is the deliverable; the campaign becomes a separate task once Larry has read it.
2. Quote the per-engine `sg_focus` resolution table and the DINA branch's `.fs_mr_family_from_table()` call, establishing what family MR receives and that it is frozen before MR runs.
3. **Quote `.fs_apply_mr()`'s argument list** and state plainly which of `ci_method`, `field_complement`, `field_scale_complement`, `field_recovery`, `field_decompose` it forwards. **It is known not to forward several.** If the field products cannot reach the DINA branch without an `R/` change, **STOP and report** — that is decision O-1, which is Larry's and is not taken here.
4. Quote the M1 DGM's DINA-relevant settings and the typical proposed-family size (how many candidates DINA yields per replicate) from any committed DINA run or a single pilot fit — this determines whether the MR field is even non-degenerate.

## Stage 1 — Smoke and projection (only if Stage 0 clears)

5-replicate smoke at the primary cell: every construction finite; interval invariants on harm, complement, `_s`, joint and IJ; γ in range; **the proposed-family size per replicate recorded (mean, min, max)** — a family of one makes selection deterministic and the correction vacuous, which is itself a finding to report; realized prevalence; p̂ and the recovery columns if they reach the branch. Time per replicate at 100 workers (DINA fits a model per replicate, so the cost profile is unlike FS — project from measurement, never from the FS walls). **Gate 1: proceed if the total projection is ≤ 6 h wall; hard timeout 8 h.**

## Stage 2 — Cells (campaign `dinamr`)

Primary setting `sg_focus = "maxeff"` (≡ `maxeffCons` on this engine — state the synonymy in the record), M1 DGM, J = 10 where applicable, `FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_RECOV=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=dinamr`, seeds `8316951 + sim_id`, sim_id 1–2000 in two seed-disjoint batches then combine, 100 workers.

| # | Cell | Purpose |
|---|---|---|
| 1 | HR 1.75 n500 | the anchor; pairs with Supplementary Table S5's design where it matches |
| 2 | HR 1.75 n1000 | does conditional coverage recover with n, as S8.3 reports for the data-driven identifiers? |
| 3 | HR 1.75 n1500 | same |
| 4 | HR 1.50 n500 | effect-size contrast |
| 5 | HR 1.50 n1000 | " |
| 6 | HR 1.00 n500 (null) | declaration rate and the level dimension on a model-generated family |

Defer order if the ceiling threatens: cell 5, then cell 3.

**Gate 2 per cell:** completeness (2,000 rows, sim_id 1–2000, no duplicates, no CONFIG-ERROR, meta knobs as set, `n_workers` and `forestsearch_version` recorded); **detection / declaration rate recorded prominently** (DINA's differs sharply from FS's and every summary conditions on it); proposed-family size distribution; every harm / complement / `_s` / joint / IJ / β(Ĥ) / β(Ĥᶜ) / p̂ quantity finite on detected replicates where the column exists; interval invariants; γ ∈ [0.025, 0.05]; bound↔quantile identities ≤ 1e-12; **columns that do not reach the DINA branch reported as absent, not silently skipped.**

## Stage 3 — Report

`summary_dinamr.qmd`, **transplanted from the committed `summary_cert20.qmd`** (globs, labels, comparator names only); `REPORT_dinamr_2026-09-10.md` beside the results with the Gate 2 record inside and every number verbatim:

1. **Standard tables** (per-cell constructions; across cells), both blocks, rows naive / field / field-s / IJ two-term where available: bias (log; marginal-SD and error-SD units), SDs, SE, r, SE/error-SD, one-sided on the exposed side [Wilson], two-sided [Wilson]. **Every coverage column labelled as coverage of the conditional-on-proposed-family estimand.**
2. **Classification:** detection/declaration rate, sens, spec, PPV, NPV, mean |Ĥ|, |Ĥ|/|H|, β(Ĥ) and β(Ĥᶜ) against the planted values, per cell and against n.
3. **The proposed family:** size distribution per cell and against n — does it stabilize, as S8.3's explanation of recovering coverage requires?
4. **Against the FS reference:** the corresponding FS numbers from `cert20` / `tier2` placed beside, **with the confound stated** — FS and DINA differ in identifier *and* in family construction, and each summary conditions on a different detection set, so the comparison is descriptive, not a contest. Quote Supplementary S8.3's own caution to that effect.
5. **Against Supplementary Table S5** where the design matches: the record's numbers beside the published ones, noting replicate counts and any design differences. **If they disagree materially, report the disagreement — do not reconcile it.**
6. **p̂ and the recovery diagnostics** on a model-generated family, if they reach the branch: distributions, and whether the low-p̂/high-p̂ bias structure found for FS appears here. Describe; do not interpret beyond the record.
7. **No acceptance criteria are pre-registered.** This is exploratory. Report the FS certified ranges as reference lines only and make no recommendation.

## Done means

Stage 0 report (**or the STOP, which is a complete deliverable in itself**); if cleared, Stage 1 and Gate 1; cells completed/deferred/dropped with walls; Gate 2 records; `summary_dinamr` rendered; the report committed; one-paragraph closing summary. **Out of scope:** any `R/` change, decision O-1, any recommendation change, GRF (a separate task if Larry wants it).
