# REPORT — Localization of the complement's residual bias by selection stratum (zero compute), with the post-merge identity gate and the field-s adoption NOTE

Date: 2026-09-08. Task: `dev/tasks/TASK_complement_location_2026-09-08_v2.md` (the spec as executed: the non-v2 text plus the Stage 0 post-merge gate, the amended compute line, the merge/install precondition and the amended Done-means; the v2 file never reached `~/Downloads`, so it was reconstructed at Stage 0 from the invocation and committed as a record fix after this report — the non-v2 file, committed as received at 0337a6f3, is superseded). Executor: Claude Code (Linux). **Report-and-wait: no repair proposal, no recommendation.** No compute beyond the 5-replicate identity gate; no `R/` change; no campaign bundle written.

Predecessors: `REPORT_field_studentize_e1_2026-09-08.md` (57e00f69; Finding 1 cell-level field-s coverage 0.912–0.921 in the four band cells, Finding 2 the p̂-tertile shape 0.936–0.954 / 0.917–0.937 / 0.877–0.887), `REPORT_complement_variance_2026-09-07.md` (Part A: variances by stratum, no means), `summary_complement_variance.qmd` (A2 stratification; extended here with A5).

## Stage 0 — post-merge identity gate (G0): PASS

- **Merge in HEAD.** `c0f48a7c` = Merge `origin/feature/glm-extension-mac` (`540c16e8`, ACTG175 continuous intervals Stage 2) into `feature/glm-extension` (first parent `57e00f69`, the E1 Stages 2–3 commit). The only `R/` file the merge changed is `R/fs_bias_coverage.R` (add-only `scale = c("log", "identity")` on `fs_sim_bias_coverage()`, default `"log"` reproducing the previous output; `man/fs_sim_bias_coverage.Rd` regenerated). Everything else the merge brought is under `dev/tasks/`, `quarto/applications/actg175/` and `quarto/simulations/actg175/continuous/`.
- **Installed package = HEAD.** forestsearch 0.3.5 at `~/R/x86_64-pc-linux-gnu-library/4.6/forestsearch` (Packaged 2026-09-08 23:53 UTC): every function in `R/` (661) has a namesake in the installed namespace (661), none only on one side, and `deparse()` is `identical()` on 661 / 661; `fs_sim_bias_coverage` carries the merged `scale` formal.
- **G0 render.** Template `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` (committed, unchanged), env only: `FS_S7_FOCUS=effMaxSG FS_S7_Z1Q=0.60 FS_S7_NBHD=0.20 FS_S7_N=500 FS_S7_HR=1.50 FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_NSIMS=5 FS_S7_START=1 FS_S7_WORKERS=5 FS_S7_CAMPAIGN=postmerge` (both knobs at their defaults: `field_decompose = FALSE`, `field_scale_complement = "none"`; J 10; seeds `8316951 + sim_id`), wall 68 s, 5 / 5 detected, mean fit+MR 36.1 s per replicate (e1inert 35.7 s). Output `results/fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_nb20_postmerge_res_1_5.rds` (158 columns) and `..._postmerge_batch_1_5.html`, committed beside the campaign's gate bundles.
- **Result: 131 / 131** pre-existing non-timing columns `identical()` to the committed `..._nb20_p30sgnb20_res_1_1000.rds` rows 1–5 (the five timing columns excluded as always); `truth` identical; the 22 post-nb20 columns (4 scale, 9 `_s`, 9 `joint_s`) present and all NA; and 153 / 153 non-timing columns identical to the committed `e1inert` gate bundle (the pre-merge G2a render of the same config). The merge changed nothing the simulation pipeline reads. It did change `R/fs_bias_coverage.R`, the **reporting** path (`fs_sim_bias_coverage()` is what `e1stud_findings.R` and the constructions tables call for coverage, Wilson limits and the Gaussian reference): the post-merge regression check on that file is Part L's exact reproduction of E1 Finding 2 (every observed coverage and every tertile count, both stratifications), computed post-merge from the same bundles.

## Part N — the adoption NOTE (D-1): committed

`dev/notes/NOTE_complement_product_2026-09-08.md`, text exactly as given in the task (commit `30c98635`). A rule statement, not analysis; it is not restated here.

## Part L — location by stratum

**Data.** The six committed `e1stud` pooled bundles (`results/*_e1stud_combined_1_2000.rds`, 2,000 rows each, 158 columns; `_s`, `joint_s`, scale-diagnostic and `p_hat_H` columns present). Detected replicates with every input finite: **1999 per cell** (the six detections, as expected). Winner-only and winner-floor excluded. Read by `summary_complement_variance.qmd` with `FS_SUMCV_GLOBS="results/*_e1stud_combined_1_2000.rds"`; the rendered `summary_complement_variance.html` (committed) is the document of record and every number below is its A5 record section verbatim.

**Definitions (per replicate, log scale).** `a = log(nv_Hc_est) − log(betaHhat_Hc)`, `c = log(nv_Hc_est) − log(mr_Hc_est)`, `e = a − c` (Part A); `lam = fld_Hc_lam_mean`, `lam_s = fld_Hc_lam_mean_s`; `ef = log(fld_Hc_est2) − log(betaHhat_Hc)`, `ef_s = log(fld_Hc_est2_s) − log(betaHhat_Hc)`. **Source quote** (`R/fs_mr_inference.R`, complement field block): `est2_w <- bdc - mean(lf)` … `est2 = to_eff(est2_w)`; under `scale_on`, `est2s_w <- bdc - mean(lfs)` … `est2_s = to_eff(est2s_w)`; with `mr_Hc_est = to_eff(bdc)` and `to_eff = exp` on the log-HR scale, `ef = e − lam` and `ef_s = e − lam_s` hold on the working scale. Asserted in the document (`stopifnot(... <= 1e-12)`): the realized maxima are 2.2e-16 to 4.4e-16 (table below). It follows that L3's correction shortfall `a − c − lam_s` **is** `ef_s`.

**Observed coverage** is the field-s one-sided upper bound's, `mean(betaHhat_Hc <= fld_Hc_up1s_s)`, with the unscaled field's (`fld_Hc_up1s`) beside it; Wilson 95% limits. **Gaussian-implied coverage** is Φ(z₀.₉₅ · r + b) with b = mean(ef_s) / SD(ef_s) and r = mean(fld_Hc_se_s) / SD(ef_s) from the stratum's own numbers (the **error** SD; the constructions tables' Gaussian reference uses the marginal SD, so the cell-level values differ slightly from E1 Finding 1's). Tertiles within the cell, as in A2 (`quantile(v, c(0, 1/3, 2/3, 1))`, `cut(..., include.lowest = TRUE)`).

**Reproduction of E1 Finding 2** (observed field-s upper coverage, p̂ tertiles T1 / T2 / T3; E1 report row vs this record's L1 rows):

| cell | E1 Finding 2 (field → field-s) | A5 L1 field | A5 L1 field-s |
|---|---|---|---|
| effMaxSG eps 0.20, HR 1.50 | 0.912 / 0.907 / 0.872 → 0.936 / 0.917 / 0.884 | 0.912 / 0.907 / 0.872 | 0.936 / 0.917 / 0.884 |
| effMaxSG eps 0.20, HR 1.75 | 0.943 / 0.916 / 0.878 → 0.954 / 0.928 / 0.877 | 0.943 / 0.916 / 0.878 | 0.954 / 0.928 / 0.877 |
| effMaxSG eps 0.30, HR 1.50 | 0.920 / 0.917 / 0.875 → 0.945 / 0.932 / 0.887 | 0.920 / 0.917 / 0.875 | 0.945 / 0.932 / 0.887 |
| effMaxSG eps 0.30, HR 1.75 | 0.918 / 0.920 / 0.872 → 0.936 / 0.937 / 0.884 | 0.918 / 0.920 / 0.872 | 0.936 / 0.937 / 0.884 |
| maxSG, HR 1.75 | 0.948 / 0.874 / 0.953 → 0.922 / 0.893 / 0.964 | 0.948 / 0.874 / 0.953 | 0.922 / 0.893 / 0.964 |
| minSG, HR 1.75 | 0.955 / 0.935 / 0.898 → 0.957 / 0.932 / 0.898 | 0.955 / 0.935 / 0.898 | 0.957 / 0.932 / 0.898 |

Every value and every tertile count (669 / 664 / 666, 667 / 666 / 666, 671 / 662 / 666, 667 / 666 / 666, 667 / 666 / 666, 669 / 664 / 666) reproduces; the |Ĥ|/|H| stratification (L2) likewise reproduces E1 Finding 2's field-s values 0.916 / 0.904 / 0.917, 0.915 / 0.908 / 0.935, 0.910 / 0.931 / 0.923, 0.913 / 0.929 / 0.914, 0.888 / 0.937 / 0.956, 0.944 / 0.931 / 0.911 and counts (667 / 678 / 654 at maxSG, 679 / 654 / 666 at minSG). The Wilson limits reproduce as well.

### Tables (verbatim from the rendered `summary_complement_variance.html`, section "Numbers for the record (A5)")

Cells are labelled as the document labels them: `HR x n500 | e1stud (focus [eps], J 10, prev 31%)`. Order: effMaxSG ε 0.20 (HR 1.50, 1.75) → ε 0.30 (HR 1.50, 1.75) → maxSG → minSG. Mean (SE) on the log scale; "all" is the whole cell.

Identity check (max |ef - (e - lam)|, max |ef_s - (e - lam_s)|; n detected / n analysed):

| cell | regime | detected | n | max abs ef | max abs ef_s |
|---|---|---|---|---|---|
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | 1999 | 1999 | 2.8e-16 | 2.8e-16 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | 1999 | 1999 | 2.8e-16 | 2.8e-16 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | 1999 | 1999 | 2.8e-16 | 2.8e-16 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | 1999 | 1999 | 2.8e-16 | 2.8e-16 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | 1999 | 1999 | 4.4e-16 | 3.3e-16 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | 1999 | 1999 | 2.2e-16 | 2.2e-16 |

L1 – location by p-hat tertile (mean (SE) on the log scale):

| cell | regime | stratum | n | mean p-hat | mean Hhat/H | mean rho^c | a | c | e | lam | ef | lam_s | ef_s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat all | 1999 | 0.103 | 0.834 | 1.054 | -0.1234 (0.0032) | -0.0790 (0.0004) | -0.0443 (0.0032) | -0.0137 (0.0002) | -0.0306 (0.0033) | -0.0135 (0.0002) | -0.0308 (0.0033) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T1 [0.000, 0.043] | 669 | 0.022 | 1.009 | 1.103 | -0.1111 (0.0057) | -0.0839 (0.0007) | -0.0272 (0.0057) | -0.0167 (0.0004) | -0.0105 (0.0058) | -0.0173 (0.0004) | -0.0099 (0.0058) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T2 [0.044, 0.109] | 664 | 0.071 | 0.814 | 1.047 | -0.1228 (0.0054) | -0.0816 (0.0006) | -0.0412 (0.0054) | -0.0144 (0.0004) | -0.0269 (0.0055) | -0.0140 (0.0004) | -0.0272 (0.0055) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T3 [0.109, 0.818] | 666 | 0.215 | 0.677 | 1.012 | -0.1362 (0.0055) | -0.0716 (0.0006) | -0.0647 (0.0056) | -0.0100 (0.0004) | -0.0546 (0.0056) | -0.0092 (0.0004) | -0.0555 (0.0056) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat all | 1999 | 0.119 | 0.846 | 1.058 | -0.1021 (0.0032) | -0.0723 (0.0004) | -0.0299 (0.0033) | -0.0108 (0.0002) | -0.0191 (0.0033) | -0.0105 (0.0003) | -0.0193 (0.0033) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T1 [0.000, 0.051] | 667 | 0.026 | 1.009 | 1.107 | -0.0768 (0.0057) | -0.0783 (0.0007) | 0.0015 (0.0057) | -0.0142 (0.0004) | 0.0157 (0.0058) | -0.0147 (0.0004) | 0.0162 (0.0058) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T2 [0.051, 0.129] | 666 | 0.085 | 0.828 | 1.049 | -0.1091 (0.0054) | -0.0751 (0.0006) | -0.0340 (0.0055) | -0.0120 (0.0004) | -0.0220 (0.0056) | -0.0115 (0.0004) | -0.0225 (0.0056) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T3 [0.129, 0.928] | 666 | 0.245 | 0.700 | 1.017 | -0.1205 (0.0055) | -0.0634 (0.0006) | -0.0571 (0.0056) | -0.0061 (0.0004) | -0.0510 (0.0056) | -0.0053 (0.0004) | -0.0518 (0.0056) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat all | 1999 | 0.064 | 1.107 | 1.087 | -0.1367 (0.0034) | -0.0894 (0.0004) | -0.0473 (0.0034) | -0.0162 (0.0003) | -0.0311 (0.0035) | -0.0157 (0.0003) | -0.0316 (0.0035) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T1 [0.000, 0.027] | 671 | 0.013 | 1.315 | 1.151 | -0.1168 (0.0060) | -0.0952 (0.0008) | -0.0217 (0.0060) | -0.0202 (0.0005) | -0.0015 (0.0061) | -0.0208 (0.0005) | -0.0009 (0.0061) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T2 [0.027, 0.065] | 662 | 0.043 | 1.078 | 1.073 | -0.1365 (0.0060) | -0.0920 (0.0007) | -0.0446 (0.0060) | -0.0170 (0.0004) | -0.0275 (0.0060) | -0.0161 (0.0005) | -0.0285 (0.0060) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T3 [0.066, 0.687] | 666 | 0.135 | 0.926 | 1.037 | -0.1568 (0.0057) | -0.0810 (0.0007) | -0.0758 (0.0057) | -0.0114 (0.0004) | -0.0644 (0.0058) | -0.0101 (0.0004) | -0.0657 (0.0058) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat all | 1999 | 0.080 | 1.094 | 1.083 | -0.1135 (0.0034) | -0.0804 (0.0004) | -0.0331 (0.0035) | -0.0121 (0.0003) | -0.0209 (0.0035) | -0.0114 (0.0003) | -0.0217 (0.0035) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T1 [0.000, 0.035] | 667 | 0.017 | 1.291 | 1.143 | -0.0968 (0.0062) | -0.0881 (0.0008) | -0.0086 (0.0062) | -0.0168 (0.0004) | 0.0081 (0.0063) | -0.0169 (0.0005) | 0.0083 (0.0063) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T2 [0.035, 0.083] | 666 | 0.056 | 1.063 | 1.069 | -0.1094 (0.0058) | -0.0833 (0.0007) | -0.0261 (0.0059) | -0.0133 (0.0004) | -0.0128 (0.0060) | -0.0121 (0.0004) | -0.0140 (0.0060) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T3 [0.084, 0.872] | 666 | 0.167 | 0.927 | 1.037 | -0.1343 (0.0057) | -0.0698 (0.0007) | -0.0646 (0.0058) | -0.0063 (0.0004) | -0.0582 (0.0058) | -0.0051 (0.0004) | -0.0595 (0.0058) |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | p-hat all | 1999 | 0.415 | 2.633 | 0.997 | -0.1914 (0.0091) | -0.1277 (0.0013) | -0.0637 (0.0088) | -0.0229 (0.0006) | -0.0408 (0.0088) | -0.0245 (0.0008) | -0.0392 (0.0088) |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | p-hat T1 [0.005, 0.123] | 667 | 0.061 | 2.197 | 0.824 | -0.3061 (0.0074) | -0.1821 (0.0010) | -0.1240 (0.0076) | -0.0390 (0.0008) | -0.0849 (0.0079) | -0.0410 (0.0007) | -0.0829 (0.0078) |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | p-hat T2 [0.123, 0.644] | 666 | 0.365 | 2.759 | 1.086 | -0.2549 (0.0173) | -0.1419 (0.0012) | -0.1130 (0.0173) | -0.0304 (0.0009) | -0.0826 (0.0174) | -0.0386 (0.0012) | -0.0744 (0.0174) |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | p-hat T3 [0.644, 1.000] | 666 | 0.821 | 2.943 | 1.081 | -0.0131 (0.0177) | -0.0589 (0.0013) | 0.0459 (0.0178) | 0.0007 (0.0010) | 0.0451 (0.0180) | 0.0061 (0.0012) | 0.0398 (0.0180) |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | p-hat all | 1999 | 0.086 | 0.406 | 1.006 | -0.0488 (0.0030) | -0.0297 (0.0003) | -0.0190 (0.0030) | -0.0061 (0.0002) | -0.0129 (0.0030) | -0.0061 (0.0002) | -0.0129 (0.0030) |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | p-hat T1 [0.000, 0.014] | 669 | 0.005 | 0.405 | 1.005 | -0.0027 (0.0050) | -0.0219 (0.0004) | 0.0192 (0.0050) | -0.0029 (0.0003) | 0.0221 (0.0050) | -0.0029 (0.0003) | 0.0221 (0.0050) |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | p-hat T2 [0.014, 0.078] | 664 | 0.037 | 0.407 | 1.006 | -0.0533 (0.0049) | -0.0317 (0.0004) | -0.0216 (0.0049) | -0.0068 (0.0003) | -0.0148 (0.0049) | -0.0068 (0.0003) | -0.0148 (0.0049) |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | p-hat T3 [0.078, 0.867] | 666 | 0.216 | 0.406 | 1.007 | -0.0904 (0.0053) | -0.0355 (0.0004) | -0.0549 (0.0052) | -0.0087 (0.0003) | -0.0461 (0.0052) | -0.0087 (0.0003) | -0.0462 (0.0052) |

L1 – coverage by p-hat tertile (field-s upper bound; Gaussian-implied from the stratum’s own mean ef_s, SD ef_s, mean fld_Hc_se_s):

| cell | regime | stratum | n | mean naive SE | mean se_field_s | SD ef_s | mean ef_s | b | r | field cov | field-s cov [Wilson] | Gaussian-implied |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat all | 1999 | 0.1426 | 0.1423 | 0.1466 | -0.0308 | -0.210 | 0.970 | 0.897 | 0.912 [0.899, 0.924] | 0.917 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T1 [0.000, 0.043] | 669 | 0.1484 | 0.1474 | 0.1494 | -0.0099 | -0.066 | 0.986 | 0.912 | 0.936 [0.915, 0.952] | 0.940 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T2 [0.044, 0.109] | 664 | 0.1418 | 0.1412 | 0.1418 | -0.0272 | -0.192 | 0.996 | 0.907 | 0.917 [0.894, 0.936] | 0.926 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T3 [0.109, 0.818] | 666 | 0.1377 | 0.1382 | 0.1451 | -0.0555 | -0.382 | 0.952 | 0.872 | 0.884 [0.858, 0.906] | 0.882 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat all | 1999 | 0.1431 | 0.1432 | 0.1490 | -0.0193 | -0.130 | 0.961 | 0.912 | 0.919 [0.907, 0.931] | 0.927 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T1 [0.000, 0.051] | 667 | 0.1484 | 0.1480 | 0.1501 | 0.0162 | 0.108 | 0.986 | 0.943 | 0.954 [0.935, 0.967] | 0.958 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T2 [0.051, 0.129] | 666 | 0.1421 | 0.1419 | 0.1442 | -0.0225 | -0.156 | 0.984 | 0.916 | 0.928 [0.906, 0.945] | 0.928 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | p-hat T3 [0.129, 0.928] | 666 | 0.1386 | 0.1396 | 0.1449 | -0.0518 | -0.358 | 0.964 | 0.878 | 0.877 [0.850, 0.900] | 0.890 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat all | 1999 | 0.1529 | 0.1525 | 0.1565 | -0.0316 | -0.202 | 0.974 | 0.904 | 0.921 [0.909, 0.932] | 0.919 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T1 [0.000, 0.027] | 671 | 0.1610 | 0.1597 | 0.1582 | -0.0009 | -0.006 | 1.010 | 0.920 | 0.945 [0.925, 0.960] | 0.951 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T2 [0.027, 0.065] | 662 | 0.1513 | 0.1506 | 0.1548 | -0.0285 | -0.184 | 0.973 | 0.917 | 0.932 [0.910, 0.949] | 0.922 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T3 [0.066, 0.687] | 666 | 0.1465 | 0.1470 | 0.1498 | -0.0657 | -0.439 | 0.982 | 0.875 | 0.887 [0.861, 0.909] | 0.880 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat all | 1999 | 0.1523 | 0.1526 | 0.1581 | -0.0217 | -0.137 | 0.965 | 0.903 | 0.919 [0.906, 0.930] | 0.927 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T1 [0.000, 0.035] | 667 | 0.1596 | 0.1589 | 0.1620 | 0.0083 | 0.051 | 0.981 | 0.918 | 0.936 [0.914, 0.952] | 0.952 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T2 [0.035, 0.083] | 666 | 0.1508 | 0.1509 | 0.1541 | -0.0140 | -0.091 | 0.979 | 0.920 | 0.937 [0.916, 0.953] | 0.936 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | p-hat T3 [0.084, 0.872] | 666 | 0.1466 | 0.1480 | 0.1504 | -0.0595 | -0.396 | 0.984 | 0.872 | 0.884 [0.858, 0.906] | 0.889 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | p-hat all | 1999 | 0.3554 | 0.3584 | 0.3954 | -0.0392 | -0.099 | 0.906 | 0.925 | 0.926 [0.914, 0.937] | 0.918 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | p-hat T1 [0.005, 0.123] | 667 | 0.2343 | 0.2292 | 0.2025 | -0.0829 | -0.410 | 1.132 | 0.948 | 0.922 [0.899, 0.940] | 0.927 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | p-hat T2 [0.123, 0.644] | 666 | 0.3741 | 0.3763 | 0.4497 | -0.0744 | -0.165 | 0.837 | 0.874 | 0.893 [0.868, 0.915] | 0.887 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | p-hat T3 [0.644, 1.000] | 666 | 0.4579 | 0.4697 | 0.4658 | 0.0398 | 0.085 | 1.008 | 0.953 | 0.964 [0.947, 0.976] | 0.959 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | p-hat all | 1999 | 0.1298 | 0.1298 | 0.1335 | -0.0129 | -0.097 | 0.972 | 0.929 | 0.929 [0.917, 0.939] | 0.933 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | p-hat T1 [0.000, 0.014] | 669 | 0.1297 | 0.1299 | 0.1298 | 0.0221 | 0.170 | 1.001 | 0.955 | 0.957 [0.938, 0.970] | 0.965 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | p-hat T2 [0.014, 0.078] | 664 | 0.1300 | 0.1303 | 0.1266 | -0.0148 | -0.117 | 1.029 | 0.935 | 0.932 [0.911, 0.949] | 0.942 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | p-hat T3 [0.078, 0.867] | 666 | 0.1298 | 0.1292 | 0.1355 | -0.0462 | -0.341 | 0.954 | 0.898 | 0.898 [0.873, 0.919] | 0.890 |

L2 – location by |Hhat|/|H| tertile (mean (SE) on the log scale):

| cell | regime | stratum | n | mean p-hat | mean Hhat/H | mean rho^c | a | c | e | lam | ef | lam_s | ef_s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H all | 1999 | 0.103 | 0.834 | 1.054 | -0.1234 (0.0032) | -0.0790 (0.0004) | -0.0443 (0.0032) | -0.0137 (0.0002) | -0.0306 (0.0033) | -0.0135 (0.0002) | -0.0308 (0.0033) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T1 [0.372, 0.690] | 667 | 0.162 | 0.569 | 1.001 | -0.1136 (0.0053) | -0.0741 (0.0006) | -0.0395 (0.0053) | -0.0121 (0.0004) | -0.0274 (0.0054) | -0.0114 (0.0004) | -0.0281 (0.0054) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T2 [0.690, 0.946] | 666 | 0.093 | 0.806 | 1.041 | -0.1251 (0.0055) | -0.0769 (0.0006) | -0.0482 (0.0056) | -0.0129 (0.0004) | -0.0353 (0.0056) | -0.0125 (0.0004) | -0.0357 (0.0056) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T3 [0.947, 2.085] | 666 | 0.052 | 1.127 | 1.120 | -0.1314 (0.0058) | -0.0861 (0.0006) | -0.0454 (0.0059) | -0.0161 (0.0004) | -0.0292 (0.0060) | -0.0167 (0.0005) | -0.0287 (0.0060) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H all | 1999 | 0.119 | 0.846 | 1.058 | -0.1021 (0.0032) | -0.0723 (0.0004) | -0.0299 (0.0033) | -0.0108 (0.0002) | -0.0191 (0.0033) | -0.0105 (0.0003) | -0.0193 (0.0033) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T1 [0.350, 0.705] | 667 | 0.184 | 0.578 | 1.002 | -0.1043 (0.0051) | -0.0666 (0.0006) | -0.0377 (0.0052) | -0.0091 (0.0004) | -0.0286 (0.0053) | -0.0085 (0.0004) | -0.0292 (0.0053) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T2 [0.705, 0.972] | 666 | 0.110 | 0.833 | 1.047 | -0.1059 (0.0056) | -0.0710 (0.0006) | -0.0349 (0.0057) | -0.0099 (0.0004) | -0.0250 (0.0058) | -0.0095 (0.0004) | -0.0255 (0.0058) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T3 [0.972, 1.929] | 666 | 0.062 | 1.127 | 1.125 | -0.0961 (0.0060) | -0.0792 (0.0006) | -0.0169 (0.0060) | -0.0132 (0.0004) | -0.0037 (0.0061) | -0.0136 (0.0005) | -0.0033 (0.0061) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H all | 1999 | 0.064 | 1.107 | 1.087 | -0.1367 (0.0034) | -0.0894 (0.0004) | -0.0473 (0.0034) | -0.0162 (0.0003) | -0.0311 (0.0035) | -0.0157 (0.0003) | -0.0316 (0.0035) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T1 [0.372, 0.968] | 667 | 0.096 | 0.772 | 1.017 | -0.1254 (0.0056) | -0.0818 (0.0007) | -0.0436 (0.0057) | -0.0131 (0.0004) | -0.0305 (0.0057) | -0.0120 (0.0004) | -0.0316 (0.0057) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T2 [0.969, 1.235] | 666 | 0.066 | 1.098 | 1.082 | -0.1234 (0.0059) | -0.0873 (0.0007) | -0.0361 (0.0060) | -0.0151 (0.0004) | -0.0211 (0.0061) | -0.0143 (0.0005) | -0.0218 (0.0061) |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T3 [1.236, 2.405] | 666 | 0.029 | 1.451 | 1.162 | -0.1613 (0.0062) | -0.0992 (0.0007) | -0.0621 (0.0063) | -0.0204 (0.0005) | -0.0417 (0.0063) | -0.0206 (0.0006) | -0.0415 (0.0064) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H all | 1999 | 0.080 | 1.094 | 1.083 | -0.1135 (0.0034) | -0.0804 (0.0004) | -0.0331 (0.0035) | -0.0121 (0.0003) | -0.0209 (0.0035) | -0.0114 (0.0003) | -0.0217 (0.0035) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T1 [0.393, 0.982] | 667 | 0.116 | 0.783 | 1.020 | -0.1052 (0.0056) | -0.0734 (0.0007) | -0.0318 (0.0057) | -0.0097 (0.0004) | -0.0221 (0.0058) | -0.0087 (0.0004) | -0.0231 (0.0058) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T2 [0.982, 1.217] | 666 | 0.087 | 1.094 | 1.080 | -0.0998 (0.0058) | -0.0784 (0.0007) | -0.0214 (0.0060) | -0.0111 (0.0005) | -0.0104 (0.0062) | -0.0101 (0.0005) | -0.0113 (0.0062) |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T3 [1.218, 2.354] | 666 | 0.036 | 1.404 | 1.150 | -0.1354 (0.0062) | -0.0894 (0.0007) | -0.0460 (0.0063) | -0.0156 (0.0005) | -0.0304 (0.0064) | -0.0153 (0.0005) | -0.0307 (0.0064) |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | Hhat/H all | 1999 | 0.415 | 2.633 | 0.997 | -0.1914 (0.0091) | -0.1277 (0.0013) | -0.0637 (0.0088) | -0.0229 (0.0006) | -0.0408 (0.0088) | -0.0245 (0.0008) | -0.0392 (0.0088) |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | Hhat/H T1 [0.454, 2.519] | 667 | 0.096 | 2.085 | 0.828 | -0.3378 (0.0071) | -0.1777 (0.0012) | -0.1601 (0.0075) | -0.0353 (0.0009) | -0.1248 (0.0078) | -0.0367 (0.0008) | -0.1234 (0.0077) |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | Hhat/H T2 [2.520, 2.892] | 678 | 0.489 | 2.733 | 1.033 | -0.1685 (0.0168) | -0.1137 (0.0023) | -0.0549 (0.0165) | -0.0209 (0.0012) | -0.0339 (0.0165) | -0.0225 (0.0015) | -0.0324 (0.0165) |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | Hhat/H T3 [2.893, 3.640] | 654 | 0.665 | 3.088 | 1.131 | -0.0658 (0.0188) | -0.0912 (0.0018) | 0.0254 (0.0187) | -0.0123 (0.0011) | 0.0377 (0.0188) | -0.0142 (0.0015) | 0.0397 (0.0188) |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | Hhat/H all | 1999 | 0.086 | 0.406 | 1.006 | -0.0488 (0.0030) | -0.0297 (0.0003) | -0.0190 (0.0030) | -0.0061 (0.0002) | -0.0129 (0.0030) | -0.0061 (0.0002) | -0.0129 (0.0030) |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | Hhat/H T1 [0.330, 0.391] | 679 | 0.083 | 0.375 | 1.004 | -0.0204 (0.0052) | -0.0274 (0.0004) | 0.0070 (0.0052) | -0.0051 (0.0003) | 0.0121 (0.0052) | -0.0051 (0.0003) | 0.0121 (0.0052) |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | Hhat/H T2 [0.391, 0.416] | 654 | 0.088 | 0.404 | 1.007 | -0.0513 (0.0052) | -0.0298 (0.0005) | -0.0215 (0.0051) | -0.0062 (0.0003) | -0.0152 (0.0051) | -0.0062 (0.0003) | -0.0152 (0.0051) |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | Hhat/H T3 [0.416, 1.042] | 666 | 0.087 | 0.440 | 1.007 | -0.0752 (0.0052) | -0.0321 (0.0005) | -0.0432 (0.0051) | -0.0070 (0.0003) | -0.0362 (0.0051) | -0.0070 (0.0003) | -0.0362 (0.0051) |

L2 – coverage by |Hhat|/|H| tertile (field-s upper bound; Gaussian-implied from the stratum’s own mean ef_s, SD ef_s, mean fld_Hc_se_s):

| cell | regime | stratum | n | mean naive SE | mean se_field_s | SD ef_s | mean ef_s | b | r | field cov | field-s cov [Wilson] | Gaussian-implied |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H all | 1999 | 0.1426 | 0.1423 | 0.1466 | -0.0308 | -0.210 | 0.970 | 0.897 | 0.912 [0.899, 0.924] | 0.917 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T1 [0.372, 0.690] | 667 | 0.1344 | 0.1344 | 0.1393 | -0.0281 | -0.202 | 0.965 | 0.918 | 0.916 [0.893, 0.935] | 0.917 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T2 [0.690, 0.946] | 666 | 0.1404 | 0.1404 | 0.1457 | -0.0357 | -0.245 | 0.964 | 0.892 | 0.904 [0.879, 0.924] | 0.910 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T3 [0.947, 2.085] | 666 | 0.1530 | 0.1520 | 0.1546 | -0.0287 | -0.185 | 0.983 | 0.881 | 0.917 [0.894, 0.936] | 0.924 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H all | 1999 | 0.1431 | 0.1432 | 0.1490 | -0.0193 | -0.130 | 0.961 | 0.912 | 0.919 [0.907, 0.931] | 0.927 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T1 [0.350, 0.705] | 667 | 0.1344 | 0.1347 | 0.1359 | -0.0292 | -0.215 | 0.991 | 0.922 | 0.915 [0.891, 0.933] | 0.921 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T2 [0.705, 0.972] | 666 | 0.1413 | 0.1417 | 0.1506 | -0.0255 | -0.169 | 0.941 | 0.905 | 0.908 [0.884, 0.928] | 0.916 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | Hhat/H T3 [0.972, 1.929] | 666 | 0.1535 | 0.1531 | 0.1585 | -0.0033 | -0.021 | 0.966 | 0.910 | 0.935 [0.914, 0.952] | 0.942 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H all | 1999 | 0.1529 | 0.1525 | 0.1565 | -0.0316 | -0.202 | 0.974 | 0.904 | 0.921 [0.909, 0.932] | 0.919 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T1 [0.372, 0.968] | 667 | 0.1405 | 0.1404 | 0.1483 | -0.0316 | -0.213 | 0.947 | 0.909 | 0.910 [0.886, 0.929] | 0.911 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T2 [0.969, 1.235] | 666 | 0.1523 | 0.1523 | 0.1565 | -0.0218 | -0.139 | 0.973 | 0.916 | 0.931 [0.909, 0.948] | 0.928 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T3 [1.236, 2.405] | 666 | 0.1660 | 0.1648 | 0.1640 | -0.0415 | -0.253 | 1.005 | 0.887 | 0.923 [0.901, 0.941] | 0.919 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H all | 1999 | 0.1523 | 0.1526 | 0.1581 | -0.0217 | -0.137 | 0.965 | 0.903 | 0.919 [0.906, 0.930] | 0.927 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T1 [0.393, 0.982] | 667 | 0.1410 | 0.1414 | 0.1501 | -0.0231 | -0.154 | 0.942 | 0.913 | 0.913 [0.889, 0.932] | 0.919 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T2 [0.982, 1.217] | 666 | 0.1522 | 0.1527 | 0.1593 | -0.0113 | -0.071 | 0.959 | 0.919 | 0.929 [0.907, 0.947] | 0.934 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | Hhat/H T3 [1.218, 2.354] | 666 | 0.1639 | 0.1637 | 0.1641 | -0.0307 | -0.187 | 0.997 | 0.878 | 0.914 [0.891, 0.933] | 0.927 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | Hhat/H all | 1999 | 0.3554 | 0.3584 | 0.3954 | -0.0392 | -0.099 | 0.906 | 0.925 | 0.926 [0.914, 0.937] | 0.918 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | Hhat/H T1 [0.454, 2.519] | 667 | 0.2257 | 0.2201 | 0.1998 | -0.1234 | -0.618 | 1.101 | 0.916 | 0.888 [0.861, 0.909] | 0.884 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | Hhat/H T2 [2.520, 2.892] | 678 | 0.3913 | 0.3963 | 0.4306 | -0.0324 | -0.075 | 0.920 | 0.925 | 0.937 [0.916, 0.953] | 0.925 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | Hhat/H T3 [2.893, 3.640] | 654 | 0.4503 | 0.4601 | 0.4815 | 0.0397 | 0.082 | 0.956 | 0.934 | 0.956 [0.937, 0.969] | 0.951 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | Hhat/H all | 1999 | 0.1298 | 0.1298 | 0.1335 | -0.0129 | -0.097 | 0.972 | 0.929 | 0.929 [0.917, 0.939] | 0.933 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | Hhat/H T1 [0.330, 0.391] | 679 | 0.1292 | 0.1290 | 0.1346 | 0.0121 | 0.090 | 0.958 | 0.946 | 0.944 [0.924, 0.959] | 0.952 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | Hhat/H T2 [0.391, 0.416] | 654 | 0.1300 | 0.1301 | 0.1306 | -0.0152 | -0.116 | 0.996 | 0.927 | 0.931 [0.909, 0.948] | 0.936 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | Hhat/H T3 [0.416, 1.042] | 666 | 0.1303 | 0.1304 | 0.1310 | -0.0362 | -0.276 | 0.995 | 0.916 | 0.911 [0.887, 0.931] | 0.913 |

L3 – regime sequence, p-hat tertiles (T1 low … T3 high p-hat):

| cell | regime | quantity | all | T1 | T2 | T3 |
|---|---|---|---|---|---|---|
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | mean a | -0.1234 | -0.1111 | -0.1228 | -0.1362 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | mean c | -0.0790 | -0.0839 | -0.0816 | -0.0716 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | mean lam_s | -0.0135 | -0.0173 | -0.0140 | -0.0092 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | mean ef_s = a - c - lam_s | -0.0308 | -0.0099 | -0.0272 | -0.0555 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | ef_s / a | 0.2498 | 0.0887 | 0.2215 | 0.4072 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | mean a | -0.1021 | -0.0768 | -0.1091 | -0.1205 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | mean c | -0.0723 | -0.0783 | -0.0751 | -0.0634 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | mean lam_s | -0.0105 | -0.0147 | -0.0115 | -0.0053 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | mean ef_s = a - c - lam_s | -0.0193 | 0.0162 | -0.0225 | -0.0518 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.20, J 10, prev 31%) | ef_s / a | 0.1894 | -0.2110 | 0.2061 | 0.4300 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | mean a | -0.1367 | -0.1168 | -0.1365 | -0.1568 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | mean c | -0.0894 | -0.0952 | -0.0920 | -0.0810 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | mean lam_s | -0.0157 | -0.0208 | -0.0161 | -0.0101 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | mean ef_s = a - c - lam_s | -0.0316 | -0.0009 | -0.0285 | -0.0657 |
| HR 1.50 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | ef_s / a | 0.2313 | 0.0075 | 0.2088 | 0.4189 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | mean a | -0.1135 | -0.0968 | -0.1094 | -0.1343 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | mean c | -0.0804 | -0.0881 | -0.0833 | -0.0698 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | mean lam_s | -0.0114 | -0.0169 | -0.0121 | -0.0051 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | mean ef_s = a - c - lam_s | -0.0217 | 0.0083 | -0.0140 | -0.0595 |
| HR 1.75 n500 | e1stud (effMaxSG eps 0.30, J 10, prev 31%) | ef_s / a | 0.1913 | -0.0859 | 0.1277 | 0.4430 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | mean a | -0.1914 | -0.3061 | -0.2549 | -0.0131 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | mean c | -0.1277 | -0.1821 | -0.1419 | -0.0589 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | mean lam_s | -0.0245 | -0.0410 | -0.0386 | 0.0061 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | mean ef_s = a - c - lam_s | -0.0392 | -0.0829 | -0.0744 | 0.0398 |
| HR 1.75 n500 | e1stud (maxSG, J 10, prev 31%) | ef_s / a | 0.2048 | 0.2709 | 0.2918 | -3.0399 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | mean a | -0.0488 | -0.0027 | -0.0533 | -0.0904 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | mean c | -0.0297 | -0.0219 | -0.0317 | -0.0355 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | mean lam_s | -0.0061 | -0.0029 | -0.0068 | -0.0087 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | mean ef_s = a - c - lam_s | -0.0129 | 0.0221 | -0.0148 | -0.0462 |
| HR 1.75 n500 | e1stud (minSG, J 10, prev 31%) | ef_s / a | 0.2647 | -8.0868 | 0.2779 | 0.5106 |

## L4 — Reading (in the record, not a task; no recommendation)

Per cell, which of (i) "the correction collapses with the stable pick" (high-p̂ tertile: mean `a` ≈ its cell-level value while `c` and `lam_s` shrink toward 0, mean `ef_s` clearly negative, prediction −0.05 to −0.09 and Gaussian-implied 0.88–0.90), (ii) "mean `a` itself → 0 in the high-p̂ tertile", or (iii) something else the numbers support; then whether the Gaussian-implied coverage tracks the observed per stratum.

- **effMaxSG ε 0.20, HR 1.50.** (i), with one amendment: the naive optimism does not merely persist in T3, it is larger there (`a` −0.1362 vs cell −0.1234 vs T1 −0.1111), while `c` shrinks from −0.0839 (T1) to −0.0716 (T3) and `lam_s` from −0.0173 to −0.0092. Residual `ef_s` −0.0099 / −0.0272 / −0.0555 (T1 / T2 / T3): inside the predicted −0.05 to −0.09 in T3. Of the −0.0247 by which T3's residual exceeds the cell's, −0.0128 is the larger optimism and −0.0119 the smaller corrections (`c` +0.0074, `lam_s` +0.0043) — about half and half. Not (ii). Gaussian-implied 0.940 / 0.926 / 0.882 vs observed 0.936 / 0.917 / 0.884: tracks (max gap 0.009).
- **effMaxSG ε 0.20, HR 1.75.** Same shape: `a` −0.0768 / −0.1091 / −0.1205 (cell −0.1021), `c` −0.0783 / −0.0751 / −0.0634, `lam_s` −0.0147 / −0.0115 / −0.0053; `ef_s` +0.0162 / −0.0225 / −0.0518. T1 is over-corrected (the corrections, −0.093 together, exceed the optimism −0.077; observed 0.954). (i) in T3 with the same amendment (optimism larger, corrections smaller); Gaussian-implied 0.958 / 0.928 / 0.890 vs observed 0.954 / 0.928 / 0.877: tracks (max gap 0.013, inside the T3 Wilson interval [0.850, 0.900]).
- **effMaxSG ε 0.30, HR 1.50.** `a` −0.1168 / −0.1365 / −0.1568 (cell −0.1367), `c` −0.0952 / −0.0920 / −0.0810, `lam_s` −0.0208 / −0.0161 / −0.0101; `ef_s` −0.0009 / −0.0285 / −0.0657 — T1 fully corrected, T3 the largest residual of the four band cells (still inside −0.05 to −0.09). (i) with the amendment; not (ii). Gaussian-implied 0.951 / 0.922 / 0.880 vs observed 0.945 / 0.932 / 0.887: tracks (max gap 0.010).
- **effMaxSG ε 0.30, HR 1.75.** `a` −0.0968 / −0.1094 / −0.1343 (cell −0.1135), `c` −0.0881 / −0.0833 / −0.0698, `lam_s` −0.0169 / −0.0121 / −0.0051; `ef_s` +0.0083 / −0.0140 / −0.0595. (i) with the amendment; not (ii). Gaussian-implied 0.952 / 0.936 / 0.889 vs observed 0.936 / 0.937 / 0.884: tracks (max gap 0.016 in T1, inside its Wilson interval [0.914, 0.952]).
- **The four band cells together (L3).** The residual is not a stable fraction of the optimism: `ef_s / a` is ≈ 0 in T1 (−0.21 to +0.09), 0.13–0.22 in T2 and **0.41–0.44 in T3**, and the T3 fraction is nearly constant across ε and HR while the T3 optimism itself moves (−0.1205 to −0.1568). Where p̂ is low the two corrections together (`c + lam_s` = −0.093 to −0.116) match or exceed the optimism; where p̂ is high they cover only 56–59% of it. The corrections scale with the pick's instability (ρᶜ 1.10–1.15 in T1 vs 1.01–1.04 in T3, hence field-s's gain being confined to T1–T2, E1 Finding 2), and the part of the naive optimism that survives a stable pick is not what either correction is measuring. A small scale component adds to the location one in T3: the naive SE falls with p̂ (0.148–0.161 in T1 to 0.138–0.147 in T3) faster than SD(ef_s) does (0.149–0.162 to 0.145–0.150), so r falls from 0.98–1.01 to 0.95–0.98.
- **maxSG, HR 1.75.** (ii) in T3 and (iii) below it. In T3 (p̂ 0.64–1.00, mean 0.821) the naive optimism vanishes (`a` −0.0131 (0.0177)) while `c` −0.0589 is still applied (`lam_s` +0.0061), so `ef_s` is +0.0398: over-corrected, observed 0.964 (Gaussian-implied 0.959). The residual sits in T1 and T2 (`ef_s` −0.0829 / −0.0744, fractions 0.27 / 0.29 of an optimism of −0.306 / −0.255); T2 is the cell's shortfall (observed 0.893, Wilson [0.868, 0.915]; Gaussian-implied 0.887) and there the naive SE is 0.837 of the error SD (the NOTE's caveat) — the location and the scale mis-calibration act together. T1's residual is masked by r = 1.132 (ρᶜ 0.824 < 1: the selected complement is smaller-scaled than the winners), observed 0.922 vs Gaussian-implied 0.927. Tracks (max gap 0.008).
- **minSG, HR 1.75.** (iii): the optimism grows with p̂ (`a` −0.0027 / −0.0533 / −0.0904, cell −0.0488) and the corrections grow with it rather than shrinking (`c` −0.0219 / −0.0317 / −0.0355, `lam_s` −0.0029 / −0.0068 / −0.0087; ρᶜ ≈ 1.005 everywhere, so field-s = field), but more slowly, so the uncorrected fraction rises 0 (T1 over-corrected, `ef_s` +0.0221, observed 0.957) → 0.28 → 0.51 (`ef_s` −0.0462 in T3, observed 0.898, Gaussian-implied 0.890). Neither the "collapse" of (i) nor the vanishing optimism of (ii). Tracks (max gap 0.010).
- **Gaussian-implied vs observed, all strata.** Across the 24 p̂ rows (six cells × all / T1 / T2 / T3) the largest |Gaussian-implied − observed| is 0.016 and across the 24 |Ĥ|/|H| rows 0.013 — every gap inside the stratum's Wilson half-width (0.017–0.025). Location (mean `ef_s`) and scale (SD `ef_s` against mean `fld_Hc_se_s`) explain the per-stratum coverage; no shape effect remains at this resolution.
- **By |Ĥ|/|H| (L2).** In the band cells the location is flat or non-monotone across |Ĥ|/|H| tertiles (`ef_s` −0.028 / −0.036 / −0.029 at ε 0.20 HR 1.50; −0.029 / −0.026 / −0.003 at ε 0.20 HR 1.75), consistent with field-s having flattened this shape (E1 Finding 2) — the residual is organised by p̂, not by region size. At maxSG the T1 |Ĥ|/|H| stratum (small regions, ρᶜ 0.828) carries the location (`ef_s` −0.1234, b −0.618) and is the 0.888 stratum.

## Side issues (flagged, not acted on)

1. The `_v2` task file named in the invocation never reached `~/Downloads`; the non-v2 file (the only copy) was committed as received and archived, Stage 0 followed the invocation text, and the v2 was reconstructed and committed as a record fix (see the header).
2. First render of the extended document: the A5 record chunk's markdown tables for L2 broke on the pipe characters in the A2-style label `|Hhat|/|H| T1 [...]`; the label is written pipe-free (`Hhat/H`) in the record chunk only (the kbl tables keep the A2 label), and the committed HTML is the second render.
3. With `FS_SUMCV_GLOBS` pointed at the e1stud set, the pre-existing A3 chunk's HR 1.00 / n 1000 tables render empty (no such cells) and the A3 plot's regime axis is unlabelled (the e1stud tag is not in `order_key`); left as is — no change to the existing chunks.
4. The three untracked `diag_*` renders and `_diag_res_*.rds` files under `gbsg_020/` and the untracked `actg175_extreme_sims_continuous_fixed_10000_payload.rds` predate this task and are untouched.

## Done means

Stage 0 PASS (gate bundle + render committed); NOTE committed (Part N); `summary_complement_variance.qmd` extended with A5 and rendered on the six e1stud bundles (HTML committed); this record committed; no campaign bundle written; branch left unpushed.
