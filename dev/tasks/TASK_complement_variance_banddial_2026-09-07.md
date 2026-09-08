# TASK — Complement variance decomposition, the identifier band dial, and template hygiene

Date: 2026-09-07. Author: chat (spec). Executor: Claude Code (Linux). Approver: Larry. Reviewer: the Linux MR-field chat (from `HANDOFF_mr_field_linux_2026-09-07.md`).
Predecessors: campaign `nb20` (`REPORT_nb20_2026-09-07.md`; arms `p30sgnb20` J = 10 and `p30sgnb20j20` J = 20; commits 095d8c4d..e8158f4c), `p30sg` (ε = 0.10), `p30` (maxeffCons). The band statement: `REPORT_p30sg_stage0_2026-09-06.md`; the `maxSG` gate map was confirmed aligned at nb20's Stage 0.

## Protocol

- First action: copy this file to `dev/tasks/` and commit. Do not push.
- Part A is analysis on committed bundles (no compute, no renders beyond one summary document). Part B is a campaign with Gate 1 as the compute go/no-go (may be pre-authorized). Part C is two document-level, add-only edits to the combined template with a knob-inert identity.
- No change under `R/`. Standing conventions: winner-only and winner-floor excluded; bounds read by location; Wilson intervals; bias in SD units; verify from source; records beside the results.

## Part A — Decomposition of the complement's error variance (no compute)

Data: the nb20 arm-B n = 500 cells (HR 1.50, HR 1.75) as primary; the nb20 arm-A cells, `p30sg` and `s7c` h175 as the regime sequence. Per detected replicate, on the log scale: `a = log(nv_Hc_est) − log(betaHhat_Hc)` (naive error), `c = log(nv_Hc_est) − log(mr_Hc_est)` (the two-term correction), `e = a − c` (de-biased error), and, if recorded separately, the complement's `selection_bias` and `fixed_bias` components of `c`.

A1. Per cell: Var(a), Var(c), Cov(a, c), Var(e) [= Var(a) + Var(c) − 2Cov(a,c), checked]; mean naive SE², mean λ-SDᶜ², mean IJ se_ij²; the ratios Var(e)/naive SE² (the regime diagnostic squared), λ-SDᶜ²/Var(e), and λ-SDᶜ²/Var(a).
A2. The same stratified by p̂(Ĥ) tertile (nb20 cells only, where p̂ is recorded) and by |Ĥ|/|H| tertile.
A3. Across the regime sequence (s7c h175 → p30 h175 → p30sg h175 → nb20 A/B h175): the three components and λ-SDᶜ², to show which component grows as the complement's variability grows.
A4. Reading (in the record, not a task): (i) if Var(a) ≈ naive SE² and the excess in Var(e) is Var(c) − 2Cov(a,c), the field under-simulates the correction's data-to-data variability — name the analysis-time objects a repair would need (the field's own uncertainty rather than conditioning on the observed field); (ii) if Var(a) exceeds naive SE², the excess is the complement's identity varying with the winner, which the field simulates and underweights — say by how much λ-SDᶜ² falls short of Var(a); (iii) either way, whether the shortfall is a stable fraction of the regime diagnostic across the sequence. The Linux chat decides from A4 whether a repair proposal is written or the documented rule stands (field bound with its shortfall stated against SD(β̃ᶜ)/naive SE and p̂; two-term IJ as the conservative option when the complement is not dominated).

Output: `summary_complement_variance.qmd` (rendered; reads the committed bundles by campaign glob) and `REPORT_complement_variance_<date>.md` beside the nb20 records.

## Part B — The identifier band dial (campaign `banddial`)

Cells: J = 10, 31% prevalence (`FS_S7_Z1Q=0.60`), `effMaxSG`, seeds `8316951 + sim_id`, `sim_id` 1–2,000, HR 1.50 n500 and HR 1.75 n500, at two settings: `FS_S7_NBHD=0.30`, and `maxSG` (`FS_S7_FOCUS=maxSG` — the aligned pure-size map; Stage 0 re-quotes it and STOPs if the alignment finding of nb20 does not hold). Four cells, `ci_method = "field"`, `include_complement = TRUE`, `field_complement = TRUE`, `return_reselection = TRUE`, FB none, 100 workers.

Stage 0: quote the `maxSG` map on both sides and the ε = 0.30 forwarding; confirm the tag guard's handling of the campaign name. Stage 1: knob-inert identity against nb20 arm A (5 replicates, all non-p̂ columns); smoke at both settings (band admits at least as many as ε = 0.20; |Ĥ| never smaller; fields finite; γ in range; alignment check on one replicate); projection at 100 workers. Gate 1: compute go (pre-authorizable; expected ≈ 2 h; ceiling 3 h, timeout 4 h). Stage 2: the four cells; Gate 2 per cell as nb20. Stage 3: **the dial table** on identical replicates, ε ∈ {0.10 (p30sg), 0.20 (nb20 A), 0.30, maxSG}: detection, mean |Ĥ|, |Ĥ|/|H| median and 90th percentile, sensitivity, specificity, PPV, NPV, share of replicates at or above the planted size, naive optimism on Ĥ (SD units), β(Ĥ) against the planted effect, β(Ĥᶜ) against 0.721, complement SD(β̃ᶜ)/naive SE and p̂(Ĥ); beside it the constructions' one-sided coverage on both blocks (field, IJ two-term) at each setting. The band is chosen by Larry on the capture/specificity trade-off; the record does not recommend one.

Output: `REPORT_banddial_<date>.md` plus rendered documents.

## Part C — Template hygiene (document-level, add-only)

C1. Pooled-meta assembly: carry `harm_z1_quantile`, `effect_neighborhood`, `er_jcuts` and `sg_focus` from the batch metas into the combined bundle's meta and include them in the combine-mode poolability keys. Identity: re-combining the committed nb20 batch bundles (with `FS_S7_SAVE_COMBINED=FALSE`) reproduces the committed pooled `results` exactly and the meta now carries the four fields.
C2. Campaign-tag guard: either accept `_` in `campaign_tag` (the stem already namespaces by tag) or leave the guard and record in the template comment that tags are alphanumeric; the nb20 record already names the tags as written (`p30sgnb20`, `p30sgnb20j20`). Larry's choice (decision T-3); default: accept underscores.

## Decisions (defaults in brackets)

- T-1 Part A cells: as listed [default].
- T-2 Part B settings: ε = 0.30 and `maxSG` on the two harm cells [default]; ε = 0.40 as a third setting only if `maxSG` and 0.30 differ materially [default: no].
- T-3 Tag guard: accept underscores [default].
- T-4 Order: A first (no compute; its record informs the reading of B), then C, then B under pre-authorization.
- T-5 Compute (Part B): pre-authorization ≤ 3 h wall at 100 workers, hard timeout 4 h, cells beyond the ceiling deferred.

## Done means

Part A record and summary document committed; Part C edits committed with the re-combine identity recorded; Part B Stage 3 report, Gate 2 records and rendered documents committed; branch left unpushed for Larry.
