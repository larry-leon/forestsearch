# REPORT — c₀ campaign record: six reads for Supplementary S1.8 (read-only)

Task: `dev/tasks/TASK_declcal_c0_campaign_reads_2026-09-24_v2.md` (committed as received, `9224dfe9`). It supersedes v1 (`dev/tasks/TASK_declcal_c0_campaign_reads_2026-09-24.md`, `097ea7ab`). Branch `feature/glm-extension`.

**Revision.** This file first answered v1's five questions (`ab25e9d6`). Under v2 it adds Q6; Q1–Q5 are unchanged. v2's Step 3 names the same report path, so here the report is a **modification** of a committed file rather than an added one. P1's "exactly two added files" is therefore met as: one added file (the v2 task document) plus one modified file (this report).

**No simulation was re-run and no rate was computed.** R was used once, for Q6 only (P2): an `Rscript --vanilla` session that read one payload's `names()`, column classes and the value of `pconsistency_digits`, and wrote nothing. The payload's md5 was the same before and after.

## 0. Gates and sources

- **G0.1.** No tracked file was modified. `git status` listed five untracked paths that were already there at session start: `actg175/binary_020/...` (two directories and two `.html` files) and `gbsg_020/scripts_dinamr/logs/nullmr_findings.err`. They are unrelated to this task, which neither touched nor committed them. I treated the tree as clean on that basis, the same scope as P1, which checks tracked files only.
- **G0.2.** Repository `forestsearch`, branch `feature/glm-extension`: pass.
- **G0.3.** All three files exist at HEAD (`5b225003` at task start). The `10b56b80` fallback was not needed.

**Files read.** Paths are relative to `quarto/simulations/gbsg_020/` unless they start with `dev/` or `R/`.

| short name | path | used for |
|---|---|---|
| REPORT | `dev/reports/REPORT_declcal_c0_campaign_2026-09-22.md` | Q2–Q5 |
| FINDINGS | `scripts_dinamr/logs/declcalc0_findings.txt` | Q2–Q5 |
| STATUS | `current_status.md` | Q1 |
| HANDOFF | `~/Downloads/HANDOFF_declaration_calibration_c0_2026-09-22.md` | the wording each question checks against |

**Where the uniform-benefit DGPs are defined:**

- `scripts_dinamr/declcalc0_inull.cells`, lines 3–8. These are the cell → (DGM, HR, n) rows.
- `scripts_dinamr/declcalc0_run.R`, at the campaign commit `a46bf9b7`:
  - lines 154–161: `uniroot` on `k_treat` so that the super-population `hr_causal` equals `FS_S7_HR`;
  - lines 209 and 212: the design-point gate's printout.
- `R/sim_aft_gbsg.R`:
  - line 236: `.create_gbsg_dgm_()`;
  - line 417: `gamma["treat"] <- k_treat * gamma["treat"]`;
  - lines 517–520: `hr_causal`, a Cox fit on the stacked potential outcomes, which is the marginal Cox HR.
- `scripts_dinamr/logs/declcalc0_inull_B{1..6}.log`, lines 2, 7 and 8: the design point as it was executed.

**Extra read, outside Step 1's list.** For Q3's fixed-p\* row I also read `scripts_dinamr/logs/declcal_fixedk_practical.txt`. REPORT line 268 names it as the source of that row. It was only read. Q6 also cites it for its header lines, 1–4.

**Q6 sources.**

- `scripts_dinamr/declcalc0_run.R` at `a46bf9b7`: its record schema (lines 221–237) and per-replicate definitions (lines 330–421).
- The committed payload `quarto/simulations/gbsg_020/results/declcalc0_inull_B4_res_1_2000.rds`: its structure only (`names()` and column classes).

---

## Q1 — Scale of the uniform-benefit designs

**Answer: two designs, and both labels are on the marginal Cox scale.** The campaign ran two uniform-benefit designs:

| design | cells | target, marginal Cox HR (`hr_causal`) | executed marginal Cox HR | uniform patient-level HR, exp(b0[treat]) | `k_treat` |
|---|---|---|---|---|---|
| "HR 0.657" | B1, B2, B3 (n 500 / 1000 / 1500) | 0.657 | 0.657000011553 | 0.582907533499 | 1.272 |
| "HR 0.721" | B4, B5, B6 (n 500 / 1000 / 1500) | 0.721 | 0.720999999973 | 0.656561772814 | 0.9918 |

**Sources:**

- `scripts_dinamr/declcalc0_inull.cells` lines 3–8 give `B1 null 0.657 500` … `B6 null 0.721 1500`. The HR is `FS_S7_HR`.
- `declcalc0_run.R@a46bf9b7` lines 159–160 solve `.hr_at(kt) - target_hr_harm = 0`, where `.hr_at` returns `setup_gbsg_dgm(model = "null", ...)$hr_causal`. `R/sim_aft_gbsg.R` lines 517–520 define `hr_causal` as the Cox HR on the stacked potential outcomes, so the target is a marginal Cox HR.
- `logs/declcalc0_inull_B1.log` lines 7–8 hold the 0.657 values, and `B4.log` lines 7–8 hold the 0.721 values. B2/B3 and B5/B6 are identical.
- STATUS line 19: "The null design points are fixed by their **marginal Cox HR** (0.657, 0.721) … **0.582908** at the 0.657 point and **0.656562** at the 0.721 point … Never call 0.657 'the uniform patient-level HR'." STATUS line 153 also has "marginal Cox HR **0.657**, **0.721**".

**Disagreement with the HANDOFF wording.**

- Both HANDOFF phrasings are consistent with the source:
  - line 200: "uniform HR 0.657 / 0.721";
  - line 222: "their patient-level HR is 0.657 at the 0.721 marginal point".
- What makes them look like one design on two scales is a numerical coincidence. The 0.721 design's patient-level HR, 0.656562, rounds to 0.657, which is also the marginal Cox HR of the other design.
- Line 200 omits the scale. Line 222's "0.657" is a patient-level value, while the design *labelled* 0.657 is a marginal value.
- STATUS line 19 forbids calling 0.657 "the uniform patient-level HR". HANDOFF line 222 uses 0.657 as a patient-level HR, although of the other design.

---

## Q2 — The conventional screen's uniform-benefit rates, cell by cell

The `p* = 0.90` screen as executed is `declared_conv`, the rounded rule on the post-reduction family (REPORT line 98). The table is copied from FINDINGS lines 31–41, the "p\* = 0.90 as executed" column. The same column appears again at lines 57–67 and in REPORT lines 105–115 and 136–146. Replicate counts come from FINDINGS lines 3–8 ("2000 rows ; ok 2000").

| cell | design (marginal Cox HR) | n | replicates | p\* 0.90 as executed [Wilson 95%] | FINDINGS line |
|---|---|---|---|---|---|
| B1 | uniform benefit HR 0.657 | 500 | 2000 | 0.0950 [0.0829, 0.1086] | 31 |
| B2 | uniform benefit HR 0.657 | 1000 | 2000 | 0.0325 [0.0256, 0.0412] | 33 |
| B3 | uniform benefit HR 0.657 | 1500 | 2000 | 0.0070 [0.0042, 0.0117] | 35 |
| B4 | uniform benefit HR 0.721 | 500 | 2000 | 0.2325 [0.2145, 0.2515] | 37 |
| B5 | uniform benefit HR 0.721 | 1000 | 2000 | 0.1110 [0.0980, 0.1255] | 39 |
| B6 | uniform benefit HR 0.721 | 1500 | 2000 | 0.0325 [0.0256, 0.0412] | 41 |

**Disagreement with the HANDOFF wording.** None. HANDOFF line 146 gives the range "0.7% to 23%", which matches the endpoints B3 0.0070 and B4 0.2325. HANDOFF line 141 quotes a different range, "realized rates of 0.007–0.25", with 0.25 as the upper end. 0.25 is B4's *pre-family exact* rate (`max_T_pre >= 1.6449`, FINDINGS line 233), not its executed rate. That is outside Q2 and is recorded here only.

---

## Q3 — Which cell is worst

Each design below is written as (marginal Cox HR, n). "Tie" means the record shows equal worst values to its 4 printed decimals. The record's own "worst B (cell)" column names a single cell even where there is a tie.

### HANDOFF §4.1 (α = 0.10)

| HANDOFF row | value | worst cell | tie? | source |
|---|---|---|---|---|
| p\* = 0.90 as executed | 0.2325 | HR 0.721, n 500 (B4) | no; next is B5 at 0.1110 | FINDINGS 37, 39; REPORT 266 |
| tuned fixed p\* 0.9545 (k 2.0) | 0.0725 | HR 0.721, n 500 (B4) | no; next is B5 at 0.0305. This is the post-reduction rate | `declcal_fixedk_practical.txt` lines 93, 111; REPORT 264 |
| c₀ 0.70 | 0.0705 | HR 0.721, n 500 (B4) | no; B5 0.0645 | FINDINGS 63, 65; REPORT 255 |
| c₀ 0.75 | 0.0350 | HR 0.721, n 500 (B4) | no; B5 0.0255 | FINDINGS 63, 65; REPORT 257 |
| c₀ 0.80 | 0.0165 | HR 0.721, n 500 (B4) | no; B5 0.0115 | FINDINGS 63, 65; REPORT 259 |
| c₀ 0.85 | 0.0055 | HR 0.721, n 500 (B4) | no; B5 0.0035 | FINDINGS 63, 65; REPORT 261 |
| c₀ = c2 = 1.0 | 0.0005 | **tie:** HR 0.657, n 500 (B1) and HR 0.721, n 500 (B4); the other B cells are 0.0000 | **yes**; the record names B1 | FINDINGS 57, 63 (last column); REPORT 263 |

### HANDOFF §4.2 (α = 0.05)

| HANDOFF row | value | worst cell | tie? | source |
|---|---|---|---|---|
| c₀ 0.70 | 0.0315 | HR 0.721, n 500 (B4) | no; B5 0.0285 | FINDINGS 37, 39; REPORT 238 |
| c₀ 0.75 | 0.0130 | HR 0.721, n 500 (B4) | no; B5 0.0100 | FINDINGS 37, 39; REPORT 240 |
| c₀ 0.80 | 0.0040 | **tie:** HR 0.721, n 500 (B4) and HR 0.721, n 1000 (B5) | **yes**; the record names B4 | FINDINGS 37, 39; REPORT 242 |
| c₀ 0.85 | 0.0010 | **tie:** HR 0.721 at n 500 (B4), n 1000 (B5) and n 1500 (B6) | **yes**; the record names B4 | FINDINGS 37, 39, 41; REPORT 244 |
| c₀ = c2 = 1.0 | 0.0005 | **tie:** HR 0.657, n 500 (B1) and HR 0.721, n 500 (B4) | **yes**; the record names B1 | FINDINGS 31, 37 (last column); REPORT 246 |

**Is it the same cell in every calibrated row?**

- **α 0.10, c₀ ∈ {0.70, 0.75, 0.80, 0.85}:** yes. HR 0.721, n 500 is the unique worst cell in all four rows.
- **α 0.05, c₀ grid:** HR 0.721, n 500 attains the worst rate in all four rows. It is unique only at c₀ 0.70 and 0.75, and ties at 0.80 (with n 1000) and 0.85 (with n 1000 and n 1500).
- **c₀ = c2 row, both α:** the record labels the worst cell B1 (HR 0.657, n 500). The value is tied with B4 (HR 0.721, n 500).

**Disagreement with the record's wording.** REPORT line 269 says "It is B4 (HR 0.721, n 500) in every calibrated exact row". STATUS line 197 says "(always B4)". Both statements are about the c₀ rows. At α 0.05 they hold only in the weak sense (B4 attains the worst rate), not the unique sense. The c₀ = c2 rows are labelled B1, although B4 has the same value.

---

## Q4 — The cutoff the campaign's Eq. 8 used

**Answer: both cutoffs.** The campaign computed Eq. 8 twice per replicate and per c₀: once against the exact cutoff 1.6449 (`fw_1645_<c0>`) and once against the executed cutoff 1.621 (`fw_1621_<c0>`). The record reports both side by side.

**Campaign runner, `scripts_dinamr/declcalc0_run.R` at `a46bf9b7`** (the version every payload ran under; REPORT lines 17–21 and 29):

```
230  z_exact <- qnorm((1 + 0.90) / 2)          # 1.644854
231  z_round <- qnorm((1 + 0.895) / 2)         # 1.621 -- effective cutoff of the rounded rule
394      r[[paste0("fw_1645_", sx)]] <- mean(Mj > z_exact)
395      r[[paste0("fw_1621_", sx)]] <- mean(Mj > z_round)
```

**Record:**

- REPORT line 168: "**fw_1645, fw_1621:** means of `mean(Mstar_c0 > cutoff)`". Both columns appear in REPORT §4, lines 171–222, and in FINDINGS lines 81–132.
- FINDINGS line 1: "Values computed at the exact cutoff, before the 2026-09-23 alignment (96f84ad8, 7713942e). See dev/tasks/TASK_declcal_consumers_2026-09-24_v3.md."

**Which one HANDOFF §4.4 quotes.** HANDOFF line 221 says "on the same pre-reduction family". In the record, that pairing is `fw_1645_c0` against `max_T_pre >= 1.6449`, the exact cutoff (FINDINGS line 219 header; REPORT lines 370 and 373). The executed-cutoff pairing, `fw_1621_c0` against `declared_conv`, is on the post-reduction family.

**Note on HEAD.** The runner at HEAD (`3dd61162`, lines 230–233) no longer records `fw_1645` or `declared_conv_exact`, and takes the cutoff from `fs_declaration_calibration(mr)$z_pstar`. That change came after the campaign. FINDINGS line 1 marks the committed values as computed before it.

---

## Q5 — The free check's margins

**Answer: confirmed, with one qualification.** Each margin is a within-cell, paired difference. Per replicate, the findings script computes the replicate's `fw_c0` minus that same replicate's declaration indicator, within one cell. It reports the mean, with SE = sd/√2000:

- `scripts_dinamr/declcalc0_findings.R` at `5ace4093` (the version that wrote FINDINGS), line 195: `pre_ex <- as.integer(r$max_T_pre >= qnorm(0.95))`;
- line 198: `d1 <- f45 - pre_ex; d2 <- f21 - r$declared_conv`;
- `r` is one cell's replicate rows.

The margins, all at c₀ = 0.70, in the HR 0.721 cells:

| cell | design | c₀ | fw_1645 − pre-family exact rate (paired SE) | fw_1621 − executed rate (paired SE) | FINDINGS line | REPORT line |
|---|---|---|---|---|---|---|
| B4 | HR 0.721, n 500 | 0.70 | 0.3113 − 0.2500 = **+0.0613** (0.0098) | 0.3248 − 0.2325 = **+0.0923** (0.0096) | 233 | 387 |
| B5 | HR 0.721, n 1000 | 0.70 | 0.1571 − 0.1050 = **+0.0521** (0.0069) | 0.1659 − 0.1110 = **+0.0549** (0.0071) | 237 | 391 |
| B6 | HR 0.721, n 1500 | 0.70 | 0.0703 − 0.0300 = **+0.0403** (0.0038) | 0.0751 − 0.0325 = **+0.0426** (0.0040) | 241 | 395 |

The record states the range at REPORT line 400: "The smallest gaps are at c0 0.70 in B4–B6 (+0.04 to +0.09)".

**Disagreement between the record and the HANDOFF wording.** The range "+0.04 to +0.09" spans **both** pairings:

- the pre-family exact-cutoff pairing (`fw_1645`) runs from +0.0403 (B6) to +0.0613 (B4);
- the executed pairing (`fw_1621`, post-reduction family) runs from +0.0426 (B6) to +0.0923 (B4).

REPORT line 358 gives the `fw_1645` pairing alone as "0.04–0.06". HANDOFF line 221 attaches "+0.04 to +0.09" to "the same pre-reduction family". The +0.09 end, however, is the post-reduction executed pairing (B4, `fw_1621`). On the pre-reduction family alone, the record's range is +0.04 to +0.06.

---

## Q6 — Is the executed rate at a screening level other than 0.90 recoverable without a re-run?

**Answer: yes. The column is `max_T_post`, not `max_T_pre`.**

**What `max_T_pre` is.** It is the maximum of T(g) over the **pre-reduction** family: the covariate-measurable enumeration, before the event floors and the near-duplicate reduction. It is **not** the family that the executed search admits from.

- `declcalc0_run.R@a46bf9b7` line 330: `fam <- .enum_family(cap$Z, maxk, aux$nmin_fit)`.
- Line 364: `r$max_T_pre <- max(T_pre)`, with `T_pre` computed over that family's field.
- REPORT line 460: "the pre-reduction family, which the calibrated rule reads, but not in the post-reduction family that the executed screen evaluated".
- `declcal_fixedk_practical.txt` line 3: "max_T_pre: pre-reduction family".

**Which committed column holds the maximum over the executed family.** `max_T_post` does. It is the maximum of T(g) over the post-reduction family, which is the family the executed screen evaluated. `declared_conv`, the "p\* = 0.90 as executed" rate (0.2325 at B4), is computed on that same family, from the same `Tp`.

- `declcalc0_run.R@a46bf9b7` line 406: `# ---- post-reduction family: the one the executed screen evaluated`.
- Lines 407 and 412: `red <- ns$.fs_decl_reduction(fit, fld$family_id)` and `Tp <- T_pre[post]`.
- Line 414: `adm_round <- round(rate, digits) >= p_star`.
- Line 415: `r$max_T_post <- max(Tp)`.
- Line 416: `r$declared_conv <- as.integer(any(adm_round))`.
- Line 417: `r$declared_conv_exact <- as.integer(any(Tp >= z_exact))`.
- Line 421: when the post-reduction family is empty, `max_T_post` is `NA_real_` and `declared_conv` is `0L`. The record counts empty families in the first campaign's payloads. B4 (HR 0.721, n 500) has 33 of 2000 (`declcal_fixedk_practical.txt` line 14).
- `declcal_fixedk_practical.txt` line 3: "max_T_post: executed post-reduction family (after d0/d1 floors)".

**Evidence that the column ties to the executed search:**

- **Fidelity gate.** `declared_conv` equals the fitted search's own declaration indicator on 2000 of 2000 replicates in every cell (REPORT lines 71–82; FINDINGS lines 15–26).
- **Identity gate.** `max_T_post` is one of the thirteen columns that match the first campaign's committed `declcal` payloads exactly, on 20,000 of 20,000 replicates (REPORT line 84).
- **Prior use.** The record has already used this column for other cutoffs. It re-evaluated `max_T_post` at fixed k ∈ {1.6449, 1.8, 2.0, 2.2, 2.4, 2.6} on the first campaign's payloads (`declcal_fixedk_practical.txt` lines 1–4 and 24–82). Its k 2.0 post-reduction column "is reproduced from the committed payloads in all 10 of 10 cells" (REPORT line 268).
- **NA convention.** That log states it: "max_T_post is NA when the post-reduction family is empty; counted as not declared (matches declared_conv_exact == (max_T_post >= qnorm(0.95)) in all 13 cells)" (line 4).

**What the payload carries** (`declcalc0_inull_B4_res_1_2000.rds`, structure read only):

- **Top level:** a list with `results`, `aux` and `meta`.
- **`results`:** a data.frame, 2000 × 79. Its columns include `G_pre`, `G_post`, `max_T_pre`, `max_T_post` (numeric), `declared_conv`, `declared_conv_exact`, `n_admitted_conv`, `n_band` and `pconsistency_digits` (integer, value 2). It also has, per c₀ ∈ {c070, c075, c080, c085}: `kappa_hat_05_*`, `kappa_hat_10_*`, `pstar_implied_05_*`, `declared_cal05_*`, `declared_cal10_*`, `n_admitted_cal05_*`, `fw_1645_*`, `fw_1621_*` and `Mstar_c0_q90/q95/q99_*`.
- **`meta`:** includes `z_exact`, `z_round`, `p_star`, `floors` and `family`.
- The schema is defined at `declcalc0_run.R@a46bf9b7` lines 221–237. STATUS line 197 describes it as "the `declcal` schema plus per-c0 columns".

**What is and is not recoverable, as the code and structure state it:**

- **Executed screen at another level.** Recoverable from `max_T_post`, with an empty post-reduction family (NA) counted as not declared.
  - The unrounded rule `any(Tp >= z)` is exactly `max_T_post >= z`. The record uses this at line 4 of the fixed-k log.
  - The executed rule also rounds (line 414, `round(rate, digits) >= p_star`, with `pconsistency_digits` = 2). The record gives the rounded p\* 0.90 rule as the cutoff 1.621 (`z_round`, runner line 231), not 1.6449. The record does not say which cutoff the rounded rule has at other p\* values. That point is **not established from source**. What the record does establish is `max_T_post` as the column that holds the executed family's maximum.
- **Pre-reduction screen at another level.** Recoverable from `max_T_pre`. Unlike `max_T_post`, it is not the executed rate.
- **FŴ(c₀) at a cutoff other than 1.6449 or 1.621.** **Not recoverable.** The per-draw shifted field `Mstar_c0` is not in the payload; only its per-replicate summaries are (`fw_1645_*`, `fw_1621_*`, `Mstar_c0_q90/q95/q99_*`). Runner line 350 passes `keep_field_matrix = FALSE`. Lines 381–398 reduce `fld$Mstar_c0` to those summaries.
- **Calibrated cutoff κ̂(c₀) at an α other than 0.10, 0.05 or 0.01.** **Not recoverable**, for the same reason. Only the 0.90, 0.95 and 0.99 quantiles of `Mstar_c0` are stored.

**Disagreement between the HANDOFF wording and the record.**

- HANDOFF lines 329–331 list the payloads as carrying "`max_T_pre`, `Mstar_c0`, `Mstar_q90/95/99`, per-replicate κ̂ and implied p\* — any cutoff, α or c₀ on the grid re-evaluable without re-running". The committed `results` frame has no `Mstar_c0` field. It has `Mstar_c0_q90/q95/q99_<c0>` summaries, so "any α" and "any cutoff" do not hold for κ̂(c₀) or FŴ(c₀).
- The same HANDOFF lines omit `max_T_post`. That is the column holding the executed-family maximum, and so the one that gives the executed rate. `max_T_pre`, which the handoff does list, gives the pre-reduction rate.
