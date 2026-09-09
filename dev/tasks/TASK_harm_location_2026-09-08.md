# TASK — Harm-side location by selection stratum on `e1stud` (zero compute): does the stable-pick under-correction sit on both blocks?

Date: 2026-09-08. Author: chat (spec). Executor: Claude Code (Linux). Approver: Larry (F-3 research thread opened 2026-09-08). Reviewer: the Linux MR-field chat.
Predecessors: `REPORT_complement_location_2026-09-08.md` (complement: uncorrected fraction ≈ 0 / 0.13–0.22 / 0.41–0.44 across p̂ tertiles, stable across ε and HR), `summary_complement_variance.qmd` §A5 (the means-by-stratum machinery), `REPORT_field_studentize_e1_2026-09-08.md` (constructions tables, harm block), `dev/notes/HANDOFF_mr_field_linux_2026-09-08.md`.

## Protocol

- First action: archive stale variants in `~/Downloads` to `~/Downloads/cc_archive/`; copy this file to `dev/tasks/` and commit. Do not push. If the file is missing from `~/Downloads`, reconstruct it from the kickoff text, commit it under this name with a note that it was reconstructed, and proceed.
- **No compute. No `R/` change. No bundle written** (assert at the end: `git status` shows no new `.rds` under `results/`). Analysis on the six committed `e1stud` pooled bundles only.
- Transplant-first: the new section A6 is the A5 chunk with the harm-block columns substituted — no fresh authorship of loading, filtering, or tertile machinery.
- Winner-only and winner-floor excluded; bounds by location; Wilson intervals; marginal and error SDs side by side. **No repair proposal; report and wait.**

## Stage 0 — Quotes for the record (from HEAD)

- The field's center: quote `R/fs_mr_inference.R` lines 776–782 of the current tree (the "Shrunk field w = beta_hat with the winner's entry replaced by the two-term de-biased estimate" comment and the `v_r`, `G_r`, `m-hat`, `Lambda*` definitions) and the line where `w` is assembled. Quote `.fs_mr_select()`'s `effMaxSG` branch.
- The template's harm-block recorder columns (names verbatim): the naive, MR two-term, and field estimates/SEs/bounds for Ĥ, the harm-side `betaHhat_H` truth column, the harm-side field `lambda_mean` if recorded, and `p_hat_H`. State which of the harm-side analogues of the complement quantities exist as columns; if the harm field's `lambda_mean` is not recorded, say so and derive the field's residual from `fld_H_est2` alone.

## A6 — Harm-side location by stratum

**Per detected replicate, log scale** (harm-block analogues of A5): `a_H = log(nv_H_est) − log(betaHhat_H)` (naive error; positive = harm over-estimated), `c_H = log(nv_H_est) − log(mr_H_est)` (two-term correction), `e_H = a_H − c_H`, `lam_H` (the harm field's own correction, if recorded), `ef_H = log(fld_H_est2) − log(betaHhat_H)`. Identity check if `lam_H` is recorded: `ef_H = e_H − lam_H` on the working scale to ≤ 1e-12.

- **H1 — Means by p̂ tertile, per cell** (tertiles within the cell as in A5; n; mean with SE of `a_H`, `c_H`, `e_H`, `lam_H`, `ef_H`; mean naive SE and field SE on the harm block; mean |Ĥ|/|H|), beside the **observed** harm-side field one-sided *lower* coverage and the IJ two-term one-sided lower coverage in the stratum, and the **Gaussian-implied** coverage from the stratum's own (mean `ef_H`, SD `ef_H`, mean field SE).
- **H2 — The same by |Ĥ|/|H| tertile.**
- **H3 — The paradox test, both blocks:** per cell, corr(p̂, `a_H`) and corr(p̂, `a` [complement, from A5]) with Spearman and Pearson; and the uncorrected fraction `ef_H / a_H` by p̂ tertile beside the complement's `ef_s / a` from A5, in one table.
- **H4 — Reading (in the record, not a task).** State, per cell, which the numbers support: (A) *competition collapse on both blocks* — `a_H` largest in the high-p̂ tertile while `c_H` (and `lam_H`) shrink, leaving `ef_H` clearly positive there with harm-side coverage lowest in that stratum; (B) *complement-only* — `ef_H` ≈ 0 and harm-side coverage flat across p̂ tertiles while the complement's residual persists; (C) something else (say what). Then whether the Gaussian-implied harm coverage tracks the observed per stratum.

**Document.** Add section "A6 — Harm-side location by stratum" to `summary_complement_variance.qmd` (copy of A5 with the harm-block columns; guarded on their presence), render with `FS_SUMCV_GLOBS` on the six `e1stud` bundles; output `REPORT_harm_location_2026-09-08.md` beside the `e1stud` records with the Stage 0 quotes and H1–H4, every number verbatim from the rendered document.

## Done means

Stage 0 quotes; A6 added and rendered; the record committed; no bundle written; branch unpushed; one-paragraph closing summary with the commit range and the H4 reading per cell in one line each. Out of scope: any repair, any compute, any `R/` change.
