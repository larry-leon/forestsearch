# TASK — Studentized complement field, E1: the R1 "field-s" construction and the six-cell evaluation campaign (report-and-wait)

Date: 2026-09-08. Author: chat (spec). Executor: Claude Code (Linux). Approver: Larry (kickoff paste = the compute go per P-5). Reviewer: the Linux MR-field chat.
Governing proposal: `dev/tasks/PROPOSAL_complement_field_scale_2026-09-08_v2.md` §5–§6. Predecessors in this line: `dev/tasks/TASK_field_studentize_stage1_e0_2026-09-08.md` (1bb8bdf1), `REPORT_field_studentize_stage1_2026-09-08.md` (7245e898), `REPORT_field_studentize_e0_2026-09-08.md` (…9d71f736). **Variant decided by the Linux chat's E0 review, on record: R1** (canonical studentized form; E0 checks i–iv passed; winner-scale CV 0.054 median removed the R2 fallback condition). Standing constraint: the effMaxSG band is capped at ε ≤ 0.25 for adoption (Larry, 2026-09-08); ε 0.30 appears below strictly as a mapped stress comparator.

## Protocol

- First action: archive stale variants in `~/Downloads` to `~/Downloads/cc_archive/`, then copy this file to `dev/tasks/` and commit. Do not push.
- **Scope: this session implements the field-s construction, gates it, runs the six-cell campaign, and reports. No adoption decision, no headline-analysis change, no E2.** The record reports; the chat and Larry decide.
- **`R/` changes authorized (both named up front): `R/fs_mr_inference.R` and the one-line forwarding pass-through in `R/forestsearch_main.R`** — classified **adds code; byte-identical defaults; new outputs only on the enabled path**. Nothing else under `R/`. Gates stop on failure without asking; fail-fast per cell; `.refuse_if_tracked()` respected; `devtools::install()` before any parallel run; committed work never re-run — per-cell pairing identity to the committed bundles is the proof.
- Winner-only and winner-floor excluded from every table, figure, and line. Bounds read by location; Wilson intervals; the marginal-vs-error SD caveat (Part A, A0) applies to any SD-unit column.

## Stage 0 — Discovery (quote from HEAD; STOP on material difference)

Quote with current line numbers: the Stage-1 decomposition block and the `sel` ensure-fit as committed at 7245e898; the `lam_c` loop; the bound assembly and `.fs_mr_field_joint()` signature and its call; the `forestsearch_main.R` `.g_mr` forwarding block; the template's `FS_S7_FIELD_DECOMP` knob, the recorder rows for `fld_Hc_se` and the four scale columns, and the existing complement bound and joint recorder columns (names verbatim — the `_s` companions mirror them). Confirm the four committed comparator bundle sets and metas: `…_nb20_p30sgnb20` (ε 0.20, HR 1.50/1.75 n500), `…_nb30_banddial` (ε 0.30, both HRs), `…_maxSG_…_banddial` and `…_minSG_…_banddial` (HR 1.75), seeds `8316951 + sim_id`, 136 columns.

## Stage 1 — The R1 construction (add-beside; defaults byte-identical)

**S1a. `R/fs_mr_inference.R`.**

1. `fs_mr_inference()` gains `field_scale_complement = c("none", "selected")` (`match.arg`; roxygen `@param`), forwarded at the complement call site. `.fs_mr_field_complement()` gains the same argument.
2. Inside `.fs_mr_field_complement()`, let `scale_on <- identical(field_scale_complement, "selected")`. When `scale_on || isTRUE(field_decompose)`: move the existing ensure-`sel`-fit (the Stage-1 block's lazy fit, unchanged) to immediately after the original lazy-fit loop, then compute `s <- sqrt(colSums(Bc * Bc))` there; the decompose block reuses `s` instead of recomputing. Output-invariance argument to record in a comment: `Zo_c[sel, ·]` / `Zi_c[sel, ·]` are only ever indexed when `sel` wins a draw, in which case `sel ∈ need` and was already fitted by the original loop — so populating `Bc[, sel]` before `crossprod` changes nothing the existing path reads. Gates G2a–G2c prove it.
3. In the existing `lam_c` loop, add (existing assignments untouched): `lam_cs <- rep(NA_real_, R_out)` before the loop, and after the `lam_c[r]` line, guarded by `scale_on`:
   `lam_cs[r] <- (s[sel]/s[G]) * Zo_c[G, r] - mean((s[sel]/s[wi[ok_in]]) * Zi_c[cbind(wi[ok_in], ok_in)])`
   — the same draws, winners, and drop logic as `lam_c`; per-draw candidate-wise rescale exactly as the proposal's R1.
4. After the existing quantile/`sd_c` lines, when `scale_on`: `lfs <- lam_cs[ok_c]`; quantiles of `lfs` (same probs, type 7); append to the `complement` list the `_s` companions — `est2_s`, `upper_1s_s`, `lower_1s_s`, `lower_2s_s`, `upper_2s_s`, `se_field_s = sd(lfs)`, `lambda_mean_s = mean(lfs)` — each the exact analogue of its unscaled twin (inverted around the same `bdc`). Compute `field$joint_s` by the identical `.fs_mr_field_joint(lam_H[ok_c], lfs, beta_deb, bdc, …)` call, attached beside `joint`. **No existing field changes under any setting; with `field_scale_complement = "none"` the function is byte-identical to 7245e898.**
5. `R/forestsearch_main.R`: `field_scale_complement = .g_mr(mr_inference_args$field_scale_complement, "none")` — the same one-line pattern as `field_decompose`.
6. `NEWS.md`: one development bullet. `devtools::document()`; `devtools::install(dependencies = FALSE)`; `deparse()` of installed bodies against source.

**S1b. Template (document-level, add-only):** knob `FS_S7_FIELD_SCALEC` (default `none`) → `mr_inference_args$field_scale_complement`; echoed in the knob audit line and MR-settings readout; recorder gains the `_s` companion of every existing `fld_Hc_*` bound/estimate/SE column plus the `joint_s` mirrors of the existing joint columns, filled `%||% NA_real_`.

**S1c. Identity gates** (template-driven, seeds `8316951 + sim_id`; compare all pre-existing non-timing columns and `truth` by `identical()`):

- **G2a — both knobs off** (nb20-A HR 1.50 config, sim_id 1–5): every pre-existing column identical to the committed `p30sgnb20` rows; all new columns (four scale + all `_s`) present and NA.
- **G2b — `FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_SCALEC=selected`** (same rows): every pre-existing column identical; `_s` columns finite on detected rows; invariants — `lower_2s_s ≤ lower_1s_s ≤ est2_s ≤ upper_1s_s ≤ upper_2s_s` on the working scale, `se_field_s > 0`, `joint_s` γ in range; report the ratio `se_field_s / (fld_Hc_scale_ratio · fld_Hc_se)` per row (R1 vs the global rescale; expected ≈ 1 at CV ≈ 0.05, recorded not gated).
- **G2c — regression on Stage 1** (e0stud config, decompose on, scale off, sim_id 1–5): the four scale columns identical to the committed `e0stud` rows — the ensure-fit/`s` relocation changed nothing.

STOP on any mismatch. Commit Stage 1 with `REPORT_field_studentize_e1_stage1_2026-09-08.md` (diff summary; gate results with concrete values).

## Gate 1 — Compute go (pre-authorized by the kickoff)

Project the six cells from a 5-replicate timing at 100 workers. Proceed only if projection ≤ **3.5 h wall**; hard timeout **5 h**; cells beyond the ceiling deferred and listed. Reference walls: band cells 31–32 min, size-rule cells 23 min (≈ 2 h 50 m total).

## Stage 2 — The six cells (campaign `e1stud`; sim_id 1–2000, two seed-disjoint batches of 1,000 then combine; 100 workers; fail-fast per cell)

All cells: `FS_S7_Z1Q=0.60`, `FS_S7_N=500`, `FS_S7_FIELD_COMPLEMENT=TRUE`, `FS_S7_IJ_RESIDUAL=two_term`, `FS_S7_FB=none`, `FS_S7_FIELD_DECOMP=TRUE`, `FS_S7_FIELD_SCALEC=selected`, `FS_S7_CAMPAIGN=e1stud`, `return_reselection = TRUE`, seeds `8316951 + sim_id`, J = 10.

| # | Setting | Cells | Committed comparator |
|---|---|---|---|
| 1–2 | `effMaxSG` ε 0.20 | HR 1.50, HR 1.75 | `p30sgnb20` |
| 3–4 | `effMaxSG` ε 0.30 (stress comparator only, per the ε ≤ 0.25 cap) | HR 1.50, HR 1.75 | `banddial` |
| 5 | `maxSG` | HR 1.75 | `banddial` |
| 6 | `minSG` | HR 1.75 | `banddial` |

**Gate 2, per cell:** completeness (2,000 rows, no duplicates, meta knobs as set including `field_scale_complement = selected`); **pairing identity — every pre-existing column `identical()` to the corresponding committed bundle** (the add-beside proof: harm block, unscaled complement, joint, identification untouched); `_s` columns finite on detected rows; interval invariants on both the unscaled and `_s` blocks; `joint_s` achieved probability within its tolerance; sens/spec/|Ĥ|/p̂ recorded. Stop-on-failure per cell.

## Stage 3 — Evaluation (report-and-wait)

On identical replicates per cell, produce the two standard summary tables (per-cell constructions; across-cells) for **Ĥ and Ĥᶜ with rows: naive, field, field-s, IJ two-term** — bias, SD, SE, SE/SD, two-sided, one-sided on the exposed side; marginal and error-scale SDs beside each other per the Part A caveat. Then the pre-registered criteria, stated as findings only:

1. **Coverage:** field-s Ĥᶜ one-sided-upper coverage in the four band cells, against ≥ 0.92 with Wilson support (committed field: 0.897–0.912); the harm block and unscaled complement unchanged (structural, asserted by pairing).
2. **Shape test:** field-s Ĥᶜ upper coverage by |Ĥ|/|H| tertile — flattened vs the field's falling pattern.
3. **No-regression ends:** field-s ≈ field at `minSG` and `maxSG` (coverage within Wilson overlap; mean ρᶜ ≈ 1 confirmed there).
4. **Bound locations:** field-s Ĥᶜ upper mean and shares < 0.85 / < 0.80 beside field and IJ; `joint_s` Bonferroni/calibrated beside `joint`; margins.
5. The `se_field_s / (ρᶜ · se_field)` distribution per cell (R1 vs global, informational).

Output: cross-cell `summary_e1stud.qmd` transplanted from the committed `summary_banddial.qmd` (copy, change the named globs/rows — no fresh authorship), rendered; `REPORT_field_studentize_e1_2026-09-08.md` with Gate 2 records beside the results. **No adoption recommendation; no change to any documented rule; report and wait.**

## Decisions (defaults in brackets)

- E-1 Variant: R1 [decided on record by the E0 review; not revisited here].
- E-2 Cells: the six listed [default].
- E-3 Knob name `FS_S7_FIELD_SCALEC`, values `none` / `selected` [default].
- E-4 Ceiling 3.5 h / timeout 5 h at 100 workers [default].
- E-5 Recorder: `_s` companions of the existing `fld_Hc_*` columns plus `joint_s` mirrors [default].

## Done means

Stage 0 quotes; Stage 1 committed with G2a/G2b/G2c PASS and concrete values; six cells with Gate 2 records; `summary_e1stud` rendered; the E1 report committed; branch left unpushed for Larry; one-paragraph closing summary listing cells completed/deferred/dropped and the commit range. Out of scope: adoption, headline analyses, any further campaign.
