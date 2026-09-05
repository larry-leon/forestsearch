# REPORT — MR de-biased gate vs Guo & He (2021) on the replicated Guo–He simulations

Date: 2026-09-04. Task: `dev/tasks/TASK_mr_vs_guohe_2026-09-04.md` + addendum A.
Full tables, figures and provenance: `quarto/GuoHe/mr_vs_guohe.qmd` (rendered `mr_vs_guohe.html`, untracked per the `GuoHe/*.html` gitignore convention; re-render reproduces it from the committed bundles). Gate 2 record: `REPORT_mr_vs_guohe_gate2_2026-09-04.md` (GREEN — 32,000/32,000 paired replicates, naive column bit-identical throughout, zero MR NAs). Identity measurements (1b + addendum A1): Stage 1 report.

Design in one line: on the identical replicates of the committed §5.1/§5.2 replication (16 cells × 2000), compare naive, Guo–He-as-replicated (stored per-replicate values, no re-run), and MR (`fs_mr_inference`, maxeff, unrestricted admission, B = 5000 centred Poisson, IJ one-sided 95% bound), all targeting θ_{k̂} on the oriented log-HR scale.

## Headline numbers (G&H shown at r = 1/30; all four r columns in the qmd)

Bias against θ_{k̂} / empirical SD / SE-to-SD / one-sided coverage, per method:

| Cell | Naive bias | G&H bias | MR bias | SD (N/G&H/MR) | SE/SD (N/G&H/MR) | Coverage (N/G&H/MR) |
|---|---|---|---|---|---|---|
| t35 β₂=0.0 | +0.102 | +0.001 | +0.028 | 0.157/0.158/0.173 | 1.19/1.02/1.76 | 0.901/0.949/0.999 |
| t35 β₂=0.1 | +0.091 | −0.008 | +0.020 | 0.156/0.157/0.174 | 1.19/1.02/1.75 | 0.912/0.956/0.999 |
| t35 β₂=0.2 | +0.077 | −0.020 | +0.016 | 0.162/0.164/0.184 | 1.15/0.97/1.69 | 0.915/0.950/0.999 |
| t35 β₂=0.3 | +0.055 | −0.038 | +0.005 | 0.177/0.179/0.201 | 1.05/0.88/1.59 | 0.910/0.947/1.000 |
| t35 β₂=0.4 | +0.030 | −0.061 | −0.009 | 0.184/0.187/0.207 | 1.01/0.84/1.58 | 0.920/0.953/0.997 |
| t35 β₂=0.5 | +0.015 | −0.072 | −0.013 | 0.185/0.188/0.206 | 1.01/0.84/1.64 | 0.933/0.954/0.997 |
| t6 k=2  | +0.101 | +0.000 | +0.026 | 0.157/0.158/0.174 | 1.19/1.02/1.75 | 0.897/0.950/0.998 |
| t6 k=6  | +0.237 | +0.004 | +0.072 | 0.125/0.126/0.155 | 1.50/1.05/1.56 | 0.747/0.948/0.987 |
| t6 k=10 | +0.289 | +0.004 | +0.088 | 0.109/0.111/0.142 | 1.72/1.11/1.57 | 0.599/0.953/0.982 |
| t6 k=12 | +0.305 | +0.002 | +0.093 | 0.108/0.109/0.143 | 1.74/1.10/1.53 | 0.555/0.949/0.972 |
| t7 β₂=0.0 | +0.122 | +0.011 | +0.045 | 0.180/0.176/0.183 | 1.03/1.05/1.91 | 0.866/0.954/0.998 |
| t7 β₂=0.1 | +0.116 | +0.005 | +0.040 | 0.186/0.184/0.192 | 1.01/1.01/1.84 | 0.870/0.954/0.998 |
| t7 β₂=0.2 | +0.112 | −0.000 | +0.039 | 0.185/0.182/0.192 | 1.04/1.02/1.87 | 0.882/0.958/0.999 |
| t7 β₂=0.3 | +0.103 | −0.010 | +0.033 | 0.189/0.185/0.198 | 1.04/1.00/1.84 | 0.895/0.962/0.999 |
| t7 β₂=0.4 | +0.094 | −0.021 | +0.028 | 0.195/0.191/0.204 | 1.02/0.97/1.79 | 0.901/0.967/0.999 |
| t7 β₂=0.5 | +0.076 | −0.041 | +0.015 | 0.196/0.192/0.207 | 1.02/0.97/1.79 | 0.918/0.973/0.999 |

(Naive coverage reproduces its stored replication values by construction — the pairing is exact. G&H SE is the implied one-sided margin /z₀.₉₅; MR SE is se_ij. MCSE on coverage ≈ 0.005; on bias 0.002–0.005.)

## What the numbers show

1. **Bias.** MR removes most of the naive optimism everywhere but, unlike G&H, never over-corrects: across the 12 β₂-graded cells MR's bias declines monotonically toward zero and stays within ±0.013 of it at β₂ ≥ 0.3 in t35, while the replicated G&H at small r goes increasingly negative (−0.038 to −0.072 in t35, −0.010 to −0.041 in t7) — the over-shrinkage of true effects its own paper's Table 5 documents. In the tie/null cells the ordering reverses: G&H is nearly exactly unbiased (|bias| ≤ 0.011) while MR retains a positive residual that grows with k (+0.026 at k=2 to +0.093 at k=12, +0.045 nested).
2. **Tie theory (addendum A1.4).** The retained fraction of naive optimism in tie scenarios is 0.275/0.261 (k=2, MCSE ≈ 0.03), 0.301/0.305/0.305 (k=6/10/12, MCSE ≈ 0.01) against the two-candidate tie constant 0.293 — consistent at every k for the disjoint families. The per-replicate implied tie residual (1−2^{−1/2})·m̂(0)·σ̂_D tracks the observed MR bias closely there: 0.0306 vs 0.0280±0.0039 (t35₀), 0.0690 vs 0.0715±0.0035 (k6), 0.0840 vs 0.0884±0.0032 (k10), 0.0889 vs 0.0931±0.0032 (k12). On the nested family the prediction under-shoots: implied 0.0327 vs observed 0.0449±0.0041 (retained fraction 0.369±0.022) — the one place the numbers depart from the supplement's tie arithmetic, consistent with the heavy overlap that M_eff ≈ 2.1 (of K ≈ 151) is summarizing.
3. **Variance and intervals.** MR's IJ SE is conservative throughout: SE/SD 1.5–1.9, one-sided coverage 0.97–1.00 with correspondingly longer margins (qmd length columns; mean one-sided margin ≈ 0.55–0.63 vs G&H ≈ 0.30–0.33 at r=1/30). The replicated G&H interval sits near-nominal (0.947–0.973; §5.2's known ≈0.01 conservatism deficit against print is footnoted in the qmd and not re-litigated here). MR's empirical SD is 5–30% larger than naive/G&H — the price of the two-term correction's own Monte Carlo/selection variability.
4. **Agreement.** Per-replicate MR and G&H debiased estimates are nearly collinear: corr 0.990–0.997 in every cell, mean |MR − G&H| 0.027–0.091, largest where the corrections themselves are largest (k=12).
5. **Diagnostics.** (A6) retained mass is structurally satisfied everywhere (per-scenario minima of A6_mass_std ≥ 1 in the qmd table; §5.1 off-diagonals exactly zero). M_eff reads as designed: ≈ k on the disjoint families (2.00/5.99/9.98/11.94), ≈ 2.13–2.24 on the ~151-candidate nested family, drifting up with β₂; p̂_Ĥ rises with separation (0.74→0.92 in t35). `ij_source = "ij"` on every replicate of every cell.
6. **Cost.** MR per replicate: 0.09 s (k=2) to ~1.4 s (151 candidates) single-core at B=5000 — versus ≈ 58 to ≈ 690 core-s for the replicated G&H at B=2000 — i.e. the MR pass is two to three orders of magnitude cheaper than the bootstrap it approximates on these designs. The addendum-A record (M_eff Monte Carlo) costs more than MR itself on the nested family (≈ 4.3 s single-core; and it is memory-bandwidth-bound, which is what inflated the Stage 2 wall-clock — actual 89.5 min / ≈179 core-h vs the 10 min / 20 core-h projection; deviation detailed in the Gate 2 record).

Interpretation stops here; the observations above are restatements of the tables. No task-blocking finding: nothing in this run proposes a change to the MR method or to the replication.

## Provenance

- Comparison bundles `mr_vs_guohe_*.rds` (16, tracked) — per-replicate records, seeds, gate tallies, sessionInfo. Driver `mr_vs_guohe_run.R`; adapter `mr_vs_guohe_sim.R` (Stage 1, identities green).
- Guo & He columns are arithmetic on the committed replication bundles; the naive column is recomputed and `identical()` to stored on all 32,000 replicates (the pairing proof).
- MR engine: `forestsearch:::fs_mr_inference` at commit `736088e3` (D6 add-only argument; suite 0 fail).
- Decisions: D1 all 16 cells; D2 pure argmax, no second MR row; D3 B=5000 centred Poisson; D4 no full-bootstrap row; D5 moot; D6 as directed.
