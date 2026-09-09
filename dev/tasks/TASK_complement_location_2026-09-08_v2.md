# TASK (v2) — Post-merge identity gate, the field-s adoption NOTE, and localization of the complement's residual bias by selection stratum (no simulation compute)

Date: 2026-09-08. Author: chat (spec). Executor: Claude Code (Linux). Approver: Larry (D-1 and D-2 approved 2026-09-08). Reviewer: the Linux MR-field chat.
Version note: v2 supersedes `TASK_complement_location_2026-09-08.md` (committed as received, 0337a6f3). The v2 file never reached `~/Downloads`; this copy is the spec as executed, reconstructed at Stage 0 from the non-v2 text plus the invocation's Stage 0 delta (record fix, one commit, after 64814a7a).
Predecessors: `REPORT_field_studentize_e1_2026-09-08.md` (…57e00f69; Finding 2 — the p̂-tertile coverage 0.877–0.887 under both field and field-s; Finding 1 — the Gaussian reference), `REPORT_complement_variance_2026-09-07.md` (Part A: variances by stratum, no means), `summary_complement_variance.qmd` (reads any bundle set via `FS_SUMCV_GLOBS`).

## Protocol

- First action: archive stale variants in `~/Downloads` to `~/Downloads/cc_archive/`; copy this file to `dev/tasks/` and commit. Do not push. Proceed unattended.
- **Precondition (merge / install).** The Mac merge (`origin/feature/glm-extension-mac`, ACTG175 continuous intervals) must be in HEAD; list the `R/` files it changed; the installed package must match HEAD by `deparse()` over the namespace before any render.
- **No simulation compute; one verification render (Stage 0, <= 5 replicates) is authorized. No `R/` change. No campaign bundle written** (assert at the end: `git status` shows no new campaign `.rds` under `results/`; the Stage 0 gate bundle is committed beside the campaign's gate bundles). Analysis on the committed `e1stud` bundles only.
- Transplant-first: the new section is built by copying the existing A2 stratification chunk of `summary_complement_variance.qmd` and changing the statistic from variances to means — no fresh authorship of the loading, filtering, or tertile machinery.
- Winner-only and winner-floor excluded; bounds by location; Wilson intervals; marginal and error SDs beside each other where an SD appears. **No repair proposal in the record; report and wait.**

## Stage 0 — Post-merge identity gate (G0)

Confirm the Mac merge is in HEAD and list the `R/` files it changed; confirm the installed package matches HEAD by `deparse()`. Run the G0 identity render: the committed template driven by env only, nb20-A HR 1.50 n500 config (`FS_S7_FOCUS=effMaxSG FS_S7_Z1Q=0.60 FS_S7_NBHD=0.20 FS_S7_N=500 FS_S7_HR=1.50 FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none`), 5 replicates (`FS_S7_NSIMS=5 FS_S7_START=1`), both knobs off (`FS_S7_FIELD_DECOMP` and `FS_S7_FIELD_SCALEC` at their defaults), tag `FS_S7_CAMPAIGN=postmerge`. Gate: **131 / 131 pre-existing non-timing columns `identical()` to the committed `..._nb20_p30sgnb20_res_1_1000.rds` rows 1–5** (timing columns excluded as always), `truth` identical, the post-nb20 columns present and NA — **or STOP.** Commit the gate bundle and render beside the campaign's.

## Part N — The adoption NOTE (D-1, decided)

Write and commit `dev/notes/NOTE_complement_product_2026-09-08.md` with exactly this content (a rule statement, not analysis):

> **Complement product (decided 2026-09-08, Larry).** The documented one-sided upper bound on β(Ĥᶜ) is **field-s** (`field_scale_complement = "selected"`; the studentized complement field, R1). Its shortfall is stated: ≈ 3 points below 0.95 at n = 500 in the effMaxSG band regimes (0.912–0.921 across the four band cells of `e1stud`), concentrated in the high-p̂ (stable-pick) stratum; the analysis-time caution flag for the complement is therefore a **high** p̂(Ĥ). `field` is retained as the unscaled comparator; the IJ two-term bound is the conservative option; the naive bound is never reported as a product. The harm-side product is unchanged: the **field** one-sided lower bound on β(Ĥ), with IJ two-term as the conservative option and `ci_method = "ij"` as the two-sided reporting default — field-s is not defined on the harm side because the harm estimate's error is not naive-SE-scaled (naive SE / error SD 0.83–0.89 in the band cells vs 0.99–1.00 for the complement). Campaign convention from here: `FS_S7_FIELD_SCALEC=selected`. The package default stays `"none"` so committed bundles remain byte-reproducible; flipping the default is a separate decision. Caveat on record: under a pure-size pick (maxSG) the naive SE itself mis-calibrates (0.876 of the error SD) and field-s inherits it; ε > 0.25 is not adoptable (Larry, 2026-09-08).

## Part L — Location by stratum (the analysis)

**Data.** The six committed `e1stud` pooled bundles (2,000 rows each; `_s`, `joint_s`, scale-diagnostic and `p_hat_H` columns present). Detected replicates with every input finite (state n per cell; expect the detections).

**Per detected replicate, log scale** (the Part A definitions, plus the two field corrections):
- `a = log(nv_Hc_est) − log(betaHhat_Hc)` (naive error); `c = log(nv_Hc_est) − log(mr_Hc_est)` (two-term correction); `e = a − c` (de-biased error).
- `lam = fld_Hc_lam_mean` (the field's additional correction on the working scale, mean of Λ*ᶜ); `lam_s = fld_Hc_lam_mean_s`; the field estimates' errors `ef = log(fld_Hc_est2) − log(betaHhat_Hc)` and `ef_s = log(fld_Hc_est2_s) − log(betaHhat_Hc)`. Verify from source (Stage 0 quote) that `fld_Hc_est2` = `to_eff(bdc − mean(lf))` so that `ef = e − lam` holds on the working scale, and assert it numerically (≤ 1e-12) as the identity check of this record.

**L1 — Means by p̂ tertile, per cell** (tertiles within the cell as in A2; report n, and the mean with its SE for each of): `a`, `c`, `e`, `lam`, `ef`, `lam_s`, `ef_s`; also mean `nv_Hc_se`, mean `fld_Hc_se_s`, mean ρᶜ, mean |Ĥ|/|H|. Beside them the **observed** field-s one-sided upper coverage in the stratum (must reproduce E1 Finding 2's 0.936–0.954 / 0.917–0.937 / 0.877–0.887 values — quote both) and the **Gaussian-implied** coverage from the stratum's own (mean `ef_s`, SD `ef_s`, mean `fld_Hc_se_s`).

**L2 — The same by |Ĥ|/|H| tertile**, per cell.

**L3 — Regime sequence across the six cells** (ε 0.20 → 0.30 → maxSG / minSG): the three stratum means of `ef_s` and of the correction shortfall `a − c − lam_s` versus the stratum's naive optimism `a`, to show whether the residual is a stable fraction of the optimism or concentrated by stratum.

**L4 — Reading (in the record, not a task).** State, per cell, which of the following the numbers support: (i) in the high-p̂ tertile the naive optimism persists (mean `a` ≈ its cell-level value) while `c` and `lam_s` shrink toward 0, leaving mean `ef_s` clearly negative (the "correction collapses with the stable pick" hypothesis; prediction: mean `ef_s` ≈ −0.05 to −0.09, Gaussian-implied ≈ 0.88–0.90); (ii) mean `a` itself → 0 in the high-p̂ tertile, so the corrections are not the issue; (iii) something else (say what). Then whether the Gaussian-implied coverage tracks the observed per stratum (location + scale explain the shape) or not (a shape effect remains). No recommendation.

**Document.** Add section "A5 — Location by stratum" to `summary_complement_variance.qmd` (copy of the A2 chunk pattern with means in place of variances; the new columns guarded so the section is skipped for bundles lacking `_s` columns), render with `FS_SUMCV_GLOBS` pointed at the six `e1stud` bundles; output `REPORT_complement_location_2026-09-08.md` beside the `e1stud` records with L1–L4, every number verbatim from the rendered document.

## Done means

Stage 0 PASS (merge in HEAD, installed = HEAD by `deparse()`, G0 131 / 131; gate bundle and render committed); NOTE committed (Part N); `summary_complement_variance.qmd` extended and rendered; the location record committed; no campaign bundle written (`git status` confirms); branch left unpushed; one-paragraph closing summary with the commit range. Out of scope: any repair, any simulation compute beyond the Stage 0 render, any change under `R/`.
