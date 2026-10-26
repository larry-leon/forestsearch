# REPORT — re-render under the lean-fit build and compare with the legacy records (2026-10-26)

Task: `dev/tasks/cc_task_rerender_compare_2026-10-26.md` (committed alone as `b3ad07e0`). Measurement only: no `R/` change, no edit to any document source, no install, no check. Record: `dev/reports/rerender_2026-10-26/`.

## 1. Provenance

| item | value |
|---|---|
| branch | `feature/vi-deprecation-lean-fits` (the merge into `master` had **not** happened; feature branch accepted per §1) |
| start HEAD | `ee8ee08e89a186e2fcc51cbfd4c566928b22f463` |
| GATE `merge-base --is-ancestor 2e2b7437 HEAD` | passes |
| `git status --porcelain` at start | `?? quarto/extreme_subgroups/fixed_random/` (pre-existing; never staged) |
| installed forestsearch | 0.4.0.9000, Built `R 4.6.1; ; 2026-10-26 17:45:59 UTC; unix` — **GATE passes** (the A+F build) |
| quarto | 1.10.18 (`/usr/local/bin/quarto`; the RStudio-bundled copy reports the same) |
| R / survival | R 4.6.1 (2026-06-24) / survival 3.8.12 |
| speff2trial | available (`requireNamespace` TRUE); xml2 TRUE (rvest not installed; not needed) |
| machine | Linux, 128 physical cores; load average 0.43 before the renders |
| legacy vignette MD5s (**GATE**) | `CRAN_0.4.0_survival-post-selection.html` 69effc528fdda3c2100074243f216878 ✓; `CRAN_0.4.0_forestsearch.html` 1fc55176f5506245a3c7e91220f06508 ✓ (no download needed) |
| legacy application HTMLs | from `git show HEAD:<path>` — both tracked. C last committed `f3f188f9` 2026-09-16; B last committed `43b051b6` 2026-08-16 |

## 2. Inventory summary (full listing: `rerender_2026-10-26/inventory.txt`)

| id | document | bootstrap | MR draws | splits | workers | eval flags | `vi.grf.min` lines | files written | render wall-clock | legacy render platform |
|---|---|---|---|---|---|---|---|---|---|---|
| V1 | `vignettes/survival-post-selection.qmd` (copy) | none | 500 (field_recovery) | 200 | sequential | all evaluated | none | none | 32 s | Mac, R 4.5.2, survival 3.8-9, forestsearch 0.4.0 (CRAN) |
| V2 | `vignettes/forestsearch.Rmd` (copy) | NB = 2 | — | 100 | multisession, 2 | `oob-cv` eval=FALSE; bootstrap / K-fold eval=TRUE | none (no longer passes it) | none | 19 s | Mac, R 4.5.2, forestsearch 0.4.0 (CRAN), glmnet 5.0 (now 5.1), grf 2.6.1 both |
| C | `quarto/applications/actg175/analysis_actg175_continuous_intervals.qmd` (in place) | none | 5000 (ci_method field, R_out/R_in 1000/500) | 500 | sequential | all evaluated | none | `_payloads/analysis_actg175_continuous_oc_intervals/analysis_actg175_continuous_oc_intervals_effmaxsg_payload.rds` (**tracked** — overwritten) | 95 s | this Linux machine (payload path `/home/larryleon/…`), **forestsearch 0.3.5** (document's own provenance line), 2026-09-16 |
| B | `quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd` (in place) | NB = 1000 × 3 fits (FS, FS-dina, FS-grf) | 5000 per fit | 500 | multisession, floor(0.8 × 128) = 102 | `run_cv`/`run_loo` FALSE → CV/LOO chunks not evaluated; gate chunks conditional on a subgroup (all ran) | L577 (`forestsearch-main`), L1040 (`fs-dina-screen`) — both `-0.2`, left as committed | `_payloads/analysis_actg175_binary_multimethod_fixed_family/…_payload.rds` (untracked, new) | 400 s | **Mac (14 cores), R 4.5.2, survival 3.8-9; sessionInfo lists forestsearch_0.2.0**, committed 2026-08-16 |

All four renders exited 0, well inside the caps (20 min / 3 h). Renders ran sequentially 15:03:14–15:12:46 PDT.

## 3. Per-document comparison

Line-diff counts (code-output lines after noise filtering; `diff_<id>.txt`, `outputs_<id>_{legacy,new}.txt`, `diff_counts.csv`):

| id | lines legacy | lines new | removed | added | key quantities | identical |
|---|---|---|---|---|---|---|
| V1 | 105 | 105 | 0 | 0 | 27 | 27 |
| V2 | 260 | 260 | 5 | 5 | 31 | 31 |
| C | 60 | 60 | 19 | 19 | 28 | 12 |
| B | 770 | 1619 | 101 | 950 | 73 | 51 (+3 "same clauses, reordered") |

"Identical" on a key quantity means string-equal at the printed precision (column `printed_precision` in each CSV gives the decimals printed).

### What changed and why

**(c) other — the findings of this task, first.**

1. **C — the MR/field layer moved; the selection did not.** Ĥ = `!{cd40 <= 507} & {gender}`, n = 79, T_obs 71.095513, Pcons 0.91, the band (4 of 8), the naive rows for Ĥ and Ĥᶜ, p̂(Ĥ) = 0.064 and the family size 4935 are identical. Everything built on the multiplier layer moved by more than a last digit: β̃ 21.79 → 22.05, IJ SE 50.60 → 50.78, field one-sided lower −38.59 → −38.31, Bonferroni Ĥ lower −51.45 → −51.16, Ĥᶜ field-s upper −17.66 → −17.64 (unstudentized −17.37 → −17.38), field outer draws used 983 → 985, complement fits 1149 → 1147, Λ*-mean +5.049 → +4.894, selection rate 0.972 → 0.974. The payload agrees (68 of 194 leaves differ, all under `extras.intervals`; `payload_C_compare.txt`). **This is not attributable to the lean fits or the deprecation from this record**: the legacy C record was rendered under **forestsearch 0.3.5** on 2026-09-16 (its own provenance line), not under the pre-deprecation 0.4.0.9000 build, so every MR/field change committed since then lies between the two renders. Same machine, so not platform. The Stage 0b continuous fixture ran `mr_inference = FALSE`, so it does not cover this layer either. Isolating it would need a render under a pre-lean build of this tree (out of scope).
2. **B — the legacy record is older than assumed.** It was rendered on the Mac (14 cores, R 4.5.2) and committed 2026-08-16; its sessionInfo lists `forestsearch_0.2.0`. Differences in B therefore combine platform, worker count (14 → 102), and roughly two months of package history, besides the deprecation.
3. **B — FS-dina Tier-1 bootstrap moved, although that fit passes no `vi.grf.min`.** Success 978 (97.8%) → 977 (97.7%); bias-corrected H OR 0.85 (0.23, 3.11) → 0.87 (0.25, 2.97). The selected DINA subgroup (`{preanti >= 849.4} & {cd40 >= 338}`, n = 98, naive OR 1.94) and its Tier-2 gate (0.955 (0.294, 3.098)) are identical. Not attributable from this record (see 2).
4. **B — Pareto frontier row of the selected FS subgroup.** Its "Split ~" and "FSBC ~" columns changed, (0.50, 25.61) / 3.30 (0.46, 23.60) → (1.23, 10.48) / 3.42 (1.17, 10.02); the other three frontier rows are identical. It is a diagnostic column only, and the cause is not established here.
5. **B — the output is twice as long.** 950 added lines are mostly a much longer DINA per-covariate frontier listing in the FS-DINA-screen print (with a "max_print" display-cap note), plus parallel-plan messages and an `sg_focus 'effMaxSG' resolves to canonical rule 'hrMaxSG'` message. These are display/verbosity changes from intervening package history, not results.

**(b) deprecation-driven family change.**

- **B, FS main fit:** passed effect screen 145 → 160; near-duplicates removed 9 → 23; consistency pool 136 → 137 (`Evaluated:`); passed consistency 5 = 5; the pre-consistency preview's "On frontier" 10 → 11. `Candidate family: 2794 (full <= maxk space); screened survivors = 145 → 160`. These are exactly the expected numbers. Selected subgroup `!{wtkg <= 86} & !{cd40 <= 380}` → `!{cd40 <= 380} & !{wtkg <= 86}`: **same clauses, reordered**, with identical N 72, naive OR 3.58 (1.35, 9.50), Pcons 0.99, so membership is identical (`membership_B.txt`). **Tier-1 bootstrap is identical** — 990/1000 successful, H 2.43 (0.80, 7.33), Hc 0.61 (0.41, 0.89), exact-match 31 (3.1%). The Tier-2 gate moved in the third decimal: de-biased OR 1.794 (0.633, 5.082) → 1.795 (0.632, 5.095); selection bias 0.6915 → 0.6922; IJ SE 0.531 → 0.532; complement 0.616 (0.370, 1.026) → 0.616 (0.370, 1.027). That is consistent with the larger screened family.
- **B, FS-DINA-screen fit** (also passes `vi.grf.min = -0.2`): passed effect screen 168 → 184; consistency pool 165 → 166. Its selected subgroup and size are identical.
- **B, FS-grf:** every key quantity is identical (selection, size 95, bootstrap 1000/1000, gate 1.097 (0.429, 2.804) / 1.091 (0.396, 3.010), candidate family 1289).
- **V2 (dropped `vi.grf.min = -0.2`):** no result moved. All 31 key quantities are identical: GRF cut `er <= 0`, LASSO 4 of 7, 406 combinations, 7 pass the effect screen, selected `{er <= 0}` n = 82 Pcons 0.96, bootstrap rows, and every CV metric. The only label change is the clause order of the third candidate, `{grade3} & {pgr <= 7}` → `{pgr <= 7} & {grade3}` (same N 72, HR 1.710, Pcons 0.860).

**(a) platform / last digit.**

- V1: nothing — 105/105 lines identical, all 27 key quantities identical (HR 2.222, field lower 0.609, field-s upper 0.812, Bonferroni 0.500/0.848, IJ (0.549, 2.900), p̂ 0.006, n_family 1744, sens_H 0.559, ppv_H 0.563, containment 0.000/0.547/1.000, share 0.244 over 998 draws), Mac R 4.5.2 vs Linux R 4.6.1 notwithstanding.
- V2: one GRF leaf `control.mean` 3.85 → 3.86 (leaf 5 of the depth-2 tree; grf 2.6.1 on both, so platform floating point), and the bootstrap "Projected for 1000 boots" time (timing noise the filter did not catch).
- B: MR harm confirmation OR 1.794 → 1.795 follows the gate (b), not platform.


#### V1 key quantities (27; identical 27)

| quantity | legacy | new | result | printed precision |
|---|---|---|---|---|
| selected subgroup | {er <= 0} & {pgr <= 26} | {er <= 0} & {pgr <= 26} | identical | integer/text |
| N | 75 | 75 | identical | integer/text |
| HR | 2.222 | 2.222 | identical | 3 |
| Pcons | 0.99 | 0.99 | identical | 2 |
| candidates evaluated | 120 | 120 | identical | integer/text |
| candidates passed | 16 | 16 | identical | integer/text |
| candidate family (family_status) | no-front-end | no-front-end | identical | integer/text |
| admission set | effect floor 0; consistency floor 0 at p* = 0.8 | effect floor 0; consistency floor 0 at p* = 0.8 | identical | 1 |
| max subgroup HR | 2.537 | 2.537 | identical | 3 |
| naive HR | 2.222 | 2.222 | identical | 3 |
| de-biased HR | 1.262 | 1.262 | identical | 3 |
| field one-sided lower (H) | 0.609 | 0.609 | identical | 3 |
| field-s one-sided upper (Hc) | 0.812 | 0.812 | identical | 3 |
| Bonferroni H lower | 0.500 | 0.500 | identical | 3 |
| Bonferroni Hc upper | 0.848 | 0.848 | identical | 3 |
| IJ two-sided (H) | (0.549, 2.900) | (0.549, 2.900) | identical | 3 |
| se_field | 0.411 | 0.411 | identical | 3 |
| se_field_s | 0.135 | 0.135 | identical | 3 |
| se_ij | 0.425 | 0.425 | identical | 3 |
| p-hat(H) | 0.006 | 0.006 | identical | 3 |
| top-3 re-selection mass | q1.1 & q18.1 0.182 \| q1.1 & q27.1 0.156 \| q10.0 & q30.0 0.056 | q1.1 & q18.1 0.182 \| q1.1 & q27.1 0.156 \| q10.0 & q30.0 0.056 | identical | 1/3 |
| n_family | 1744 | 1744 | identical | integer/text |
| sens_H | 0.559 | 0.559 | identical | 3 |
| ppv_H | 0.563 | 0.563 | identical | 3 |
| containment q10/q50/q90 | 0.000 / 0.547 / 1.000 | 0.000 / 0.547 / 1.000 | identical | 3 |
| share containing all of Hhat | 0.244 | 0.244 | identical | 3 |
| recovery draws | 998 | 998 | identical | integer/text |

#### V2 key quantities (31; identical 31)

| quantity | legacy | new | result | printed precision |
|---|---|---|---|---|
| GRF subgroup | er <= 0 | er <= 0 | identical | integer/text |
| Cox-LASSO selected | 4 of 7 candidate factors | 4 of 7 candidate factors | identical | integer/text |
| LASSO omitted | age, meno, er | age, meno, er | identical | integer/text |
| candidate factors | 14 | 14 | identical | integer/text |
| combinations evaluated | 406 | 406 | identical | integer/text |
| passed effect threshold | 7 | 7 | identical | integer/text |
| candidates found | 7 | 7 | identical | integer/text |
| consistency evaluated | 3 | 3 | identical | integer/text |
| consistency passed | 3 | 3 | identical | integer/text |
| Pcons (selected) | 0.96 | 0.96 | identical | 2 |
| selected subgroup | {er <= 0} | {er <= 0} | identical | integer/text |
| subgroup size | 82 (12.0% of ITT) | 82 (12.0% of ITT) | identical | 1 |
| max HR estimate | 2.54 | 2.54 | identical | 2 |
| ITT row | 686 (100.0%) \| 246 (35.9%) \| 299 (43.6%) \| 66.3 \| 50.2 \| 7.8 \| 0.69 (0.54, 0.89) | 686 (100.0%) \| 246 (35.9%) \| 299 (43.6%) \| 66.3 \| 50.2 \| 7.8 \| 0.69 (0.54, 0.89) | identical | 1/2 |
| Questionable (H) row | 82 (12.0%) \| 26 (31.7%) \| 45 (54.9%) \| 22.9 \| 43.7 \| -14 \| 1.95 (1.04, 3.67) | 82 (12.0%) \| 26 (31.7%) \| 45 (54.9%) \| 22.9 \| 43.7 \| -14 \| 1.95 (1.04, 3.67) | identical | 1/2 |
| Recommend (Hc) row | 604 (88.0%) \| 220 (36.4%) \| 254 (42.1%) \| 66.7 \| 52.6 \| 9.3 \| 0.61 (0.47, 0.80) | 604 (88.0%) \| 220 (36.4%) \| 254 (42.1%) \| 66.7 \| 52.6 \| 9.3 \| 0.61 (0.47, 0.80) | identical | 1/2 |
| bootstrap success | 2 (100.0%) | 2 (100.0%) | identical | 1 |
| SG_CIs H row (bias-corrected) | 82 (12.0%) \| 26 (31.7%) \| 45 (54.9%) \| 22.9 \| 43.7 \| -14 \| 1.95 (1.04, 3.67) \| 3.0206 (1e-04,148775.2913) | 82 (12.0%) \| 26 (31.7%) \| 45 (54.9%) \| 22.9 \| 43.7 \| -14 \| 1.95 (1.04, 3.67) \| 3.0206 (1e-04,148775.2913) | identical | 1/2/4 |
| SG_CIs Hc row (bias-corrected) | 604 (88.0%) \| 220 (36.4%) \| 254 (42.1%) \| 66.7 \| 52.6 \| 9.3 \| 0.61 (0.47, 0.80) \| 0.75 (0.19,2.86) | 604 (88.0%) \| 220 (36.4%) \| 254 (42.1%) \| 66.7 \| 52.6 \| 9.3 \| 0.61 (0.47, 0.80) \| 0.75 (0.19,2.86) | identical | 1/2 |
| H bias-corrected HR | 3.021 | 3.021 | identical | 3 |
| Hc bias-corrected HR | 0.746 | 0.746 | identical | 3 |
| H bootstrap CV | 38.5% | 38.5% | identical | 1 |
| Hc bootstrap CV | 18.1% | 18.1% | identical | 1 |
| subgroup identified in bootstrap | 100.0 | 100.0 | identical | 1 |
| CV Sensitivity (H) | 19.5 | 19.5 | identical | 1 |
| CV Sensitivity (Hc) | 76.2 | 76.2 | identical | 1 |
| CV PPV (H) | 10.0 | 10.0 | identical | 1 |
| CV PPV (Hc) | 87.5 | 87.5 | identical | 1 |
| CV Any Found | 100.0 | 100.0 | identical | 1 |
| CV Exact Match | 0.0 | 0.0 | identical | 1 |
| CV At Least 1 | 0.0 | 0.0 | identical | 1 |

#### C key quantities (28; identical 12)

| quantity | legacy | new | result | printed precision |
|---|---|---|---|---|
| N | 1083 | 1083 | identical | integer/text |
| ITT MD | -27.591 (SE 7.889) | -27.591 (SE 7.889) | identical | 3 |
| selected subgroup | !{cd40 <= 507} & {gender} | !{cd40 <= 507} & {gender} | identical | integer/text |
| n(H-hat) | 79 | 79 | identical | integer/text |
| T_obs (MD, oriented) | 71.095513 | 71.095513 | identical | 6 |
| Pcons | 0.91000000 | 0.91000000 | identical | 8 |
| band line | 87.916667; band floor (1 − ε) × max = 70.333333; 4 of 8 qualifying candidates in band; effMaxSG selects the largest in band | 87.916667; band floor (1 − ε) × max = 70.333333; 4 of 8 qualifying candidates in band; effMaxSG selects the largest in band | identical | 6 |
| H Naive | 71.10 \| 38.86 (robust) \| (-5.07, 147.26) \| 7.17 \| -71.10 \| (-147.26, 5.07) \| -7.17 | 71.10 \| 38.86 (robust) \| (-5.07, 147.26) \| 7.17 \| -71.10 \| (-147.26, 5.07) \| -7.17 | identical | 2 |
| H MR (IJ two-term) | 21.79 (β̃) \| 50.60 (IJ) \| (-77.38, 120.96) \| -61.44 \| -21.79 \| (-120.96, 77.38) \| 61.44 | 22.05 (β̃) \| 50.78 (IJ) \| (-77.48, 121.58) \| -61.47 \| -22.05 \| (-121.58, 77.48) \| 61.47 | changed | 2 |
| H MR (field) | 16.74 (est₂); β̃ 21.79 \| 32.31 (λ-SD) \| (-51.45, 73.46) \| -38.59 \| -16.74 \| (-73.46, 51.45) \| 38.59 | 17.16 (est₂); β̃ 22.05 \| 32.47 (λ-SD) \| (-51.16, 74.78) \| -38.31 \| -17.16 \| (-74.78, 51.16) \| 38.31 | changed | 2 |
| Hc Naive | -36.15 \| 7.78 (robust) \| (-51.41, -20.90) \| -23.35 \| 36.15 \| (20.90, 51.41) \| 23.35 | -36.15 \| 7.78 (robust) \| (-51.41, -20.90) \| -23.35 \| 36.15 \| (20.90, 51.41) \| 23.35 | identical | 2 |
| Hc MR (IJ two-term) | -32.74 (β̃ᶜ) \| 15.06 (IJ) \| (-62.25, -3.23) \| -7.97 \| 32.74 \| (3.23, 62.25) \| 7.97 | -32.72 (β̃ᶜ) \| 15.05 (IJ) \| (-62.23, -3.21) \| -7.96 \| 32.72 \| (3.21, 62.23) \| 7.96 | changed | 2 |
| Hc MR (field, unstudentized) | -31.10 (est₂); β̃ᶜ -32.74 \| 8.09 (λ-SDᶜ) \| (-46.69, -15.38) \| -17.37 \| 31.10 \| (15.38, 46.69) \| 17.37 | -31.09 (est₂); β̃ᶜ -32.72 \| 8.09 (λ-SDᶜ) \| (-46.67, -15.37) \| -17.38 \| 31.09 \| (15.37, 46.67) \| 17.38 | changed | 2 |
| Hc MR (field-s) | -31.12 (est₂ˢ); β̃ᶜ -32.74 \| 7.89 (λ-SDᶜ, studentized) \| (-46.32, -15.82) \| -17.66 \| 31.12 \| (15.82, 46.32) \| 17.66 | -31.11 (est₂ˢ); β̃ᶜ -32.72 \| 7.89 (λ-SDᶜ, studentized) \| (-46.29, -15.82) \| -17.64 \| 31.11 \| (15.82, 46.29) \| 17.64 | changed | 2 |
| regime diagnostic | λ-SDᶜ / naive SE = 1.039 unstudentized, 1.014 studentized (the IJ SE / naive SE = 1.935). | λ-SDᶜ / naive SE = 1.039 unstudentized, 1.013 studentized (the IJ SE / naive SE = 1.934). | changed | 3 |
| joint subtitle | 0.025 on the 0.025–0.050 grid, corr(Λ*, Λ*ᶜ) = -0.038 over 983 aligned draws; field-s: γ = 0.025, corr = -0.035 over 983 | 0.025 on the 0.025–0.050 grid, corr(Λ*, Λ*ᶜ) = -0.042 over 985 aligned draws; field-s: γ = 0.025, corr = -0.039 over 985 | changed | 3 |
| separate one-sided pair | -38.59 \| -17.37 \| 38.59 \| 17.37 \| < 0.95 by construction | -38.31 \| -17.38 \| 38.31 \| 17.38 \| < 0.95 by construction | changed | 2 |
| Bonferroni pair (unstudentized) | -51.45 \| -15.38 \| 51.45 \| 15.38 \| 0.949 | -51.16 \| -15.37 \| 51.16 \| 15.37 \| 0.949 | changed | 2/3 |
| Bonferroni pair (field-s) | -51.45 \| -15.82 \| 51.45 \| 15.82 \| 0.949 | -51.16 \| -15.82 \| 51.16 \| 15.82 \| 0.949 | changed | 2/3 |
| p-hat(H) | 0.064 | 0.064 | identical | 3 |
| top re-selected | !{wtkg <= 81} & !{cd40 <= 413} \| 0.077 | !{wtkg <= 81} & !{cd40 <= 413} \| 0.078 | changed | 3 |
| 2nd re-selected | {age <= 37} & !{cd40 <= 507} \| 0.077 | {age <= 37} & !{cd40 <= 507} \| 0.077 | identical | 3 |
| 3rd re-selected | {age <= 40} & !{cd40 <= 507} \| 0.066 | {age <= 40} & !{cd40 <= 507} \| 0.066 | identical | 3 |
| selection rate; n_family | 0.972; 4935 candidates | 0.974; 4935 candidates | changed | 3 |
| field Lambda*-mean | +5.049 | +4.894 | changed | 3 |
| field lambda-SD ratios | 32.31 / 50.60 = 0.64; 32.31 / 38.86 = 0.83 | 32.47 / 50.78 = 0.64; 32.47 / 38.86 = 0.84 | changed | 2 |
| field draw usage | 983 / 1000; 499 / 500 | 985 / 1000; 499 / 500 | changed | integer/text |
| complement fits; lambda-SDc ratio | 1149; 8.09 / 7.78 = 1.039 | 1147; 8.09 / 7.78 = 1.039 | changed | 2/3 |

#### B key quantities (73; identical 51)

| quantity | legacy | new | result | printed precision |
|---|---|---|---|---|
| ITT OR line | 0.665  (95% CI: 0.520, 0.850) | 0.665  (95% CI: 0.520, 0.850) | identical | 3 |
| FS: combinations evaluated | 3240 | 3240 | identical | integer/text |
| FS: passed effect screen | 145 | 160 | changed | integer/text |
| FS: near-duplicates removed | 9 | 23 | changed | integer/text |
| FS: consistency pool (Evaluated) | 136 | 137 | changed | integer/text |
| FS: consistency passed | 5 | 5 | identical | integer/text |
| FS: on frontier | 10 | 11 | changed | integer/text |
| FS: selected subgroup (search print) | !{wtkg <= 86} & !{cd40 <= 380} | !{cd40 <= 380} & !{wtkg <= 86} | same clauses, reordered | integer/text |
| FS: MR harm confirmation | OR = 1.794 (point rule vs 1.00) -> consistent with harm | OR = 1.795 (point rule vs 1.00) -> consistent with harm | changed | 3/2 |
| FS: candidate subgroups found (table) | 145 | 160 | changed | integer/text |
| FS: max OR estimate | 3.58 | 3.58 | identical | 2 |
| FS: ITT row | 1083 (100.0%) \| 522 (48.2%) \| 44.39% \| 34.67% \| -9.71% \| 0.67 (0.52, 0.85) | 1083 (100.0%) \| 522 (48.2%) \| 44.39% \| 34.67% \| -9.71% \| 0.67 (0.52, 0.85) | identical | 1/2 |
| FS: H row (naive) | 72 (6.6%) \| 34 (47.2%) \| 36.84% \| 67.65% \| +30.80% \| 3.58 (1.35, 9.50) | 72 (6.6%) \| 34 (47.2%) \| 36.84% \| 67.65% \| +30.80% \| 3.58 (1.35, 9.50) | identical | 1/2 |
| FS: Hc row (naive) | 1011 (93.4%) \| 488 (48.3%) \| 44.93% \| 32.38% \| -12.56% \| 0.59 (0.45, 0.76) | 1011 (93.4%) \| 488 (48.3%) \| 44.93% \| 32.38% \| -12.56% \| 0.59 (0.45, 0.76) | identical | 1/2 |
| FS: identified subgroup (H) | !{wtkg <= 86} & !{cd40 <= 380} | !{cd40 <= 380} & !{wtkg <= 86} | same clauses, reordered | integer/text |
| FS: subgroup size | 72 (6.6% of ITT) | 72 (6.6% of ITT) | identical | 1 |
| FS: Pareto selected row | !{wtkg <= 86} \| !{cd40 <= 380} \| 72 \| 37 \| 3.58 \| 0.99 \| 2 \| (1.35, 9.50) \| (0.50, 25.61) \| 3.30 (0.46, 23.60) | !{cd40 <= 380} \| !{wtkg <= 86} \| 72 \| 37 \| 3.58 \| 0.99 \| 2 \| (1.35, 9.50) \| (1.23, 10.48) \| 3.42 (1.17, 10.02) | changed | 2 |
| FS vs GRF agreement | 1015 / 1083 (93.7%) | 1015 / 1083 (93.7%) | identical | 1 |
| FS vs GRF kappa | 0.632 | 0.632 | identical | 3 |
| FS: bootstrap success | 990 (99.0%) | 990 (99.0%) | identical | 1 |
| FS: SG_CIs H (bootstrap) | 72 (6.6%) \| 34 (47.2%) \| 36.84% \| 67.65% \| +30.80% \| 3.58 (1.35, 9.50) \| 2.43 (0.80,7.33) | 72 (6.6%) \| 34 (47.2%) \| 36.84% \| 67.65% \| +30.80% \| 3.58 (1.35, 9.50) \| 2.43 (0.80,7.33) | identical | 1/2 |
| FS: SG_CIs Hc (bootstrap) | 1011 (93.4%) \| 488 (48.3%) \| 44.93% \| 32.38% \| -12.56% \| 0.59 (0.45, 0.76) \| 0.61 (0.41,0.89) | 1011 (93.4%) \| 488 (48.3%) \| 44.93% \| 32.38% \| -12.56% \| 0.59 (0.45, 0.76) \| 0.61 (0.41,0.89) | identical | 1/2 |
| FS: SG_CIs H (MR gate) | 72 (6.6%) \| 34 (47.2%) \| 36.84% \| 67.65% \| +30.80% \| 3.58 (1.35, 9.50) \| 1.79 (0.63,5.08) | 72 (6.6%) \| 34 (47.2%) \| 36.84% \| 67.65% \| +30.80% \| 3.58 (1.35, 9.50) \| 1.79 (0.63,5.09) | changed | 1/2 |
| FS: SG_CIs Hc (MR gate) | 1011 (93.4%) \| 488 (48.3%) \| 44.93% \| 32.38% \| -12.56% \| 0.59 (0.45, 0.76) \| 0.62 (0.37,1.03) | 1011 (93.4%) \| 488 (48.3%) \| 44.93% \| 32.38% \| -12.56% \| 0.59 (0.45, 0.76) \| 0.62 (0.37,1.03) | identical | 1/2 |
| FS: H bias-corrected OR | 2.427 | 2.427 | identical | 3 |
| FS: Hc bias-corrected OR | 0.606 | 0.606 | identical | 3 |
| FS: H bootstrap CV | 79.9% | 79.9% | identical | 1 |
| FS: boot exact match with original | 31 (3.1%) | 31 (3.1%) | identical | 1 |
| FS: boot partial match | 210 (21.2%) | 210 (21.2%) | identical | 1 |
| FS: gate de-biased OR (T2 \| T1) | 1.794 (0.633, 5.082) 2.427 (0.803, 7.333) | 1.795 (0.632, 5.095) 2.427 (0.803, 7.333) | changed | 3 |
| FS: gate complement de-biased OR (T2 \| T1) | 0.616 (0.370, 1.026) 0.606 (0.411, 0.893) | 0.616 (0.370, 1.027) 0.606 (0.411, 0.893) | changed | 3 |
| FS: candidate family (family_status) | 2794 (full <= maxk space); screened survivors = 145. | 2794 (full <= maxk space); screened survivors = 160. | changed | integer/text |
| FS: scale cross-check | PASS (3.584 vs 3.584). | PASS (3.584 vs 3.584). | identical | 3 |
| FS: gate harm confirmed | TRUE | TRUE | identical | integer/text |
| FS: gate selection bias (log-OR) | 0.6915 | 0.6922 | changed | 4 |
| FS: IJ SE vs robust SE | IJ SE 0.531 vs robust SE 0.512 | IJ SE 0.532 vs robust SE 0.512 | changed | 3 |
| FS: gate complement | 1011): naive OR 0.587 -> de-biased 0.616 (selection bias -0.0499, log-OR). | 1011): naive OR 0.587 -> de-biased 0.616 (selection bias -0.0501, log-OR). | changed | 3/4 |
| DINA-screen: GRF cut | (wtkg <= 83 & preanti > 1292) \| (wtkg > 83 & cd40 > 390) | (wtkg <= 83 & preanti > 1292) \| (wtkg > 83 & cd40 > 390) | identical | integer/text |
| DINA-screen: combinations evaluated | 3081 | 3081 | identical | integer/text |
| DINA-screen: passed effect screen | 168 | 184 | changed | integer/text |
| DINA-screen: consistency pool | 165 | 166 | changed | integer/text |
| DINA-screen: consistency passed | 17 | 17 | identical | integer/text |
| DINA-screen: selected subgroup | !{wtkg <= 83.4624} & !{cd40 <= 413} | !{wtkg <= 83.4624} & !{cd40 <= 413} | identical | 4 |
| DINA-screen: subgroup size | 71 (6.6% of ITT) | 71 (6.6% of ITT) | identical | 1 |
| standalone DINA: selection 1 | cd40 >= 500 | cd40 >= 500 | identical | integer/text |
| standalone DINA: selection 2 | drugs >=   0 | drugs >=   0 | identical | integer/text |
| standalone DINA: direction/N 1 | right \| N = 100 | right \| N = 100 | identical | integer/text |
| standalone DINA: direction/N 2 | left \| N = 100 | left \| N = 100 | identical | integer/text |
| FS-dina: candidates searched | 10284 | 10284 | identical | integer/text |
| FS-dina: candidates qualifying | 147 | 147 | identical | integer/text |
| FS-dina: SELECTED | {preanti >= 849.4} & {cd40 >= 338.0}  (n = 98, mean tau-hat = 0.0098) | {preanti >= 849.4} & {cd40 >= 338.0}  (n = 98, mean tau-hat = 0.0098) | identical | 1/4 |
| FS-dina: selected subgroup (H) | {preanti >= 849.4} & {cd40 >= 338} | {preanti >= 849.4} & {cd40 >= 338} | identical | 1 |
| FS-dina: subgroup size | 98 (9.0% of ITT) | 98 (9.0% of ITT) | identical | 1 |
| FS-dina: bootstrap success | 978 (97.8%) | 977 (97.7%) | changed | 1 |
| FS-dina: SG_CIs H (bootstrap) | 98 (9.0%) \| 43 (43.9%) \| 49.09% \| 65.12% \| +16.03% \| 1.94 (0.85, 4.40) \| 0.85 (0.23,3.11) | 98 (9.0%) \| 43 (43.9%) \| 49.09% \| 65.12% \| +16.03% \| 1.94 (0.85, 4.40) \| 0.87 (0.25,2.97) | changed | 1/2 |
| FS-dina: SG_CIs Hc (bootstrap) | 985 (91.0%) \| 479 (48.6%) \| 43.87% \| 31.94% \| -11.93% \| 0.60 (0.46, 0.78) \| 0.63 (0.40,1.02) | 985 (91.0%) \| 479 (48.6%) \| 43.87% \| 31.94% \| -11.93% \| 0.60 (0.46, 0.78) \| 0.63 (0.40,1.02) | identical | 1/2 |
| FS-dina: gate de-biased OR (T2 \| T1) | 0.955 (0.294, 3.098) 0.851 (0.233, 3.106) | 0.955 (0.294, 3.098) 0.867 (0.253, 2.969) | changed | 3 |
| FS-dina: gate complement de-biased OR (T2 \| T1) | 0.632 (0.376, 1.062) 0.634 (0.395, 1.017) | 0.632 (0.376, 1.062) 0.634 (0.395, 1.017) | identical | 3 |
| FS-dina: gate candidate family | 147 | 147 | identical | integer/text |
| FS-dina: gate selection bias (log-OR) | 0.7081 | 0.7081 | identical | 4 |
| FS-grf: SELECTED | wtkg > 84.369600000000005 & cd40 > 368.20000000000005 | wtkg > 84.369600000000005 & cd40 > 368.20000000000005 | identical | 15/14 |
| FS-grf: selected subgroup (H) | {wtkg > 84.37} & {cd40 > 368.2} | {wtkg > 84.37} & {cd40 > 368.2} | identical | 2/1 |
| FS-grf: subgroup size | 95 (8.8% of ITT) | 95 (8.8% of ITT) | identical | 1 |
| FS-grf: bootstrap success | 1000 (100.0%) | 1000 (100.0%) | identical | 1 |
| FS-grf: SG_CIs H (bootstrap) | 95 (8.8%) \| 47 (49.5%) \| 39.58% \| 59.57% \| +19.99% \| 2.25 (0.99, 5.11) \| 1.09 (0.40,3.01) | 95 (8.8%) \| 47 (49.5%) \| 39.58% \| 59.57% \| +19.99% \| 2.25 (0.99, 5.11) \| 1.09 (0.40,3.01) | identical | 1/2 |
| FS-grf: SG_CIs Hc (bootstrap) | 988 (91.2%) \| 475 (48.1%) \| 44.83% \| 32.21% \| -12.62% \| 0.58 (0.45, 0.76) \| 0.63 (0.38,1.04) | 988 (91.2%) \| 475 (48.1%) \| 44.83% \| 32.21% \| -12.62% \| 0.58 (0.45, 0.76) \| 0.63 (0.38,1.04) | identical | 1/2 |
| FS-grf: gate de-biased OR (T2 \| T1) | 1.097 (0.429, 2.804) 1.091 (0.396, 3.010) | 1.097 (0.429, 2.804) 1.091 (0.396, 3.010) | identical | 3 |
| FS-grf: gate complement de-biased OR (T2 \| T1) | 0.622 (0.370, 1.044) 0.627 (0.378, 1.039) | 0.622 (0.370, 1.044) 0.627 (0.378, 1.039) | identical | 3 |
| FS-grf: gate candidate family | 1289 | 1289 | identical | integer/text |
| FS-grf: gate selection bias (log-OR) | 0.7261 | 0.7261 | identical | 4 |
| summary: definition | !{wtkg <= 86} & !{cd40 <= 380} | !{cd40 <= 380} & !{wtkg <= 86} | same clauses, reordered | integer/text |
| summary: H size | 72 (6.6%) | 72 (6.6%) | identical | 1 |
| summary: Hc size | 1011 (93.4%) | 1011 (93.4%) | identical | 1 |


## 4. Deprecation warnings (`deprecation_counts.csv`)

| id | legacy `grep -c deprecated` | new `grep -c deprecated` | new `vi.grf.min` deprecation warnings |
|---|---|---|---|
| V1 | 0 | 0 | 0 |
| V2 | 0 | 0 | 0 |
| C | 0 | 0 | 0 |
| B | 1 | 6 | 2 — one per `forestsearch()` call that passes it (L577, L1040), as expected |

The other "deprecated" hits in B are forestploter notices (`refline_*`, `title/arrow_lab/footnote`, `xlim/ticks_at`, `theme` "will be deprecated"): 1 in the legacy record, 4 now.

## 5. Working-tree state (left for Larry)

Modified tracked files (in-place renders):

| path | size |
|---|---|
| `quarto/applications/actg175/analysis_actg175_continuous_intervals.html` | 2,210,477 B (legacy 2,212,626) |
| `quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.html` | 3,269,521 B (legacy 3,316,211) |
| `quarto/applications/actg175/_payloads/analysis_actg175_continuous_oc_intervals/analysis_actg175_continuous_oc_intervals_effmaxsg_payload.rds` | 2,401 B (legacy 2,395) — written by C's `export-payload` chunk |

New untracked file (by-product of B): `quarto/applications/actg175/_payloads/analysis_actg175_binary_multimethod_fixed_family/analysis_actg175_binary_multimethod_fixed_family_payload.rds` (1,923 B). C left no `_files/` directory.

**For Larry, in GitHub Desktop:** commit the two HTMLs and C's payload `.rds` (plus B's new payload, if wanted) to accept them as the new records; discard them to keep the legacy records. Keep C's HTML and its payload together either way.

## 6. Deviations (verbatim)

1. The merge into `master` had not happened; ran on `feature/vi-deprecation-lean-fits` at `ee8ee08e` (accepted by §1).
2. C's render overwrote a **tracked non-HTML** file — its payload `.rds` — so the third §6 status line also lists that file, not only the two HTMLs. The legacy payload was saved from git as `legacy_C_payload.rds` and compared leaf-by-leaf (`payload_C_compare.txt`). No document source was edited to prevent it.
3. "Nothing else running": while B rendered, CC ran a few short single-process `Rscript` HTML extractions on finished files (seconds each, on 26 idle cores). They do not touch B's computation; B's own timings may carry a negligible effect.
4. Output extraction, beyond the §4.1 recipe: (a) V2's knitr output uses the default `## ` prefix, not `#>`, so both prefixes are accepted. (b) Tables inside quarto `cell-output-display` blocks and rmarkdown sections (gt/kable) are included as `<<table>>` rows (`cell | cell`), because C and B report most quantities in tables. "Computational Timing" tables are excluded as wall-clock noise. (c) Additional noise filters: remaining sessionInfo cruft (`Matrix products`, `locale:`, `LC_*=`, `[1] stats graphics…`), `gc()` memory tables (`Ncells`/`Vcells`), `Average per iteration`, `System max cores`, and absolute `payload to: /…` paths (Mac `/Users/…` vs Linux `/home/…`). V2's "Projected for 1000 boots" timing line is not filtered and shows as a diff line. All filters are in `extract_lib.R`.
5. Key quantities from B's tier comparisons are printed as `results: asis` markdown (not code output), so `compare.R` also reads those prose/table lines (`source = prose` in `keyvalues_B.csv`).
6. No document prints a field named `family_status`; the `Candidate family:` line is used (V1: `no-front-end`; B: `2794 (full <= maxk space); screened survivors = N`). C and V2 print none.
7. `FS: on frontier` in `keyvalues_B.csv` is read from the pre-consistency preview line (10 → 11); the post-consistency summary's "On frontier: 4" is identical.
8. B membership `setequal()` was not possible: no committed or on-disk legacy payload exists. The judgment rests on identical N, naive OR and complement rows (`membership_B.txt`).
9. `new_C.html` was copied to `~/Downloads` right after its render (§3), and again at the end with the others.
10. `legacy_C_payload.rds` (2.4 KB, binary) is committed with the record, as a non-HTML file under the record directory.

## Verdict

1. Gates: provenance ✓, A+F build ✓, CRAN MD5s ✓, all four renders exit 0 inside their caps (32 s / 19 s / 95 s / 400 s).
2. **Finding (c):** C's MR/field quantities moved (β̃ 21.79 → 22.05, Ĥ field lower −38.59 → −38.31, Ĥᶜ field-s upper −17.66 → −17.64) with identical selection. The legacy C record is a **0.3.5** render, so this measures the change since 2026-09-16, not the lean fits.
3. **Finding (c):** the legacy B record is an August Mac render (sessionInfo `forestsearch_0.2.0`, 14 cores), so B compares across platform, workers and two months of history. The FS-dina Tier-1 bootstrap moved (0.85 → 0.87; 978 → 977 successes) and the selected FS row's Pareto Split/FSBC CIs moved; neither is attributable from this record.
4. V1 (survival + MR + field recovery) is identical in every code-output line and all 27 key quantities.
5. V2 is identical in all 31 key quantities despite dropping `vi.grf.min = -0.2`. Its only differences are one GRF leaf mean in the last digit, one reordered label and a timing projection.
6. B FS main: the pool moved exactly as expected (145 → 160 screened, 136 → 137 pool). The selected subgroup has the same clauses, reordered, and identical membership (N 72, OR 3.58).
7. B FS Tier-1 bootstrap is identical (2.43 (0.80, 7.33)). The Tier-2 gate moved in the third decimal (1.794 → 1.795, upper CI 5.082 → 5.095).
8. B FS-grf is identical throughout. FS-DINA-screen's pool moved (168 → 184; 165 → 166) with the same selection.
9. Deprecation warnings: exactly 2 in B (one per call that passes `vi.grf.min`), none elsewhere.
10. To isolate lean-fit/deprecation effects from older history in C and B, the clean comparator is a re-render of C and B under a pre-lean build of this tree. That is Larry's call; it was not done here.

## 7. Closing block

(§6 post-conditions, appended after the commit.)
