# S2 example trial: definitions (TASK_s2-example-km_2026-10-06)

- Route: A (a committed replicate reproduced). Cell S2_rmlog4 of the campaign at r = −log 4 (b = 0.55, h = 2, n = 1000, 10000 replicates); records `results/ident_rlog4_full/S2_rmlog4` at commit `fca4c1e`; build `build_full("S2")` of `sims/ident_rlog4_full_design.R` (its build record equals the committed `analytic/build.csv`); forestsearch 0.3.5.9002; repo HEAD at the run `e997ff15a47a6f5b336bd743ed68896d1d9e960d`.
- Selection: of the 10000 replicates, 2398 reach the trigger (ITT ub < 1.0 and region ITT > 1.0), 2397 of those declare Ĥ, and 1228 also have region(Ĝ) < 1.0. The first by index is replicate 5; replicates 1–5 were examined (`s2_example_trials.csv`).
- Chosen replicate: 5, seed 613800005 (call seed 613800000 + 5), redrawn by `simulate_from_dgm()` under L'Ecuyer-CMRG with the campaign's arguments and the identifier run on its 806 non-region patients with the campaign's arguments (identifier run 10.1 s, 0 warnings).
- Reproduction against the committed record (P5): ITT HR 0.808414 vs 0.808414, region HR 1.036400 vs 1.036400, region(Ĝ) HR 0.843608 vs 0.843608, Ĥ `{z_cd40 <= 236} & !{z_preanti <= 196}` vs `{z_cd40 <= 236} & !{z_preanti <= 196}`; largest absolute difference over the four sets' HR and limits 4.9e-15.

## Ĥ, in the identifier's notation

`{z_cd40 <= 236} & !{z_preanti <= 196}` (`fs_result$sg.harm` joined by ` & `, as `mrct_region_sims()` records it; `!{…}` is the complement). H* = {z_cd40 ≤ 230} (`flag_harm`). G* = region minus H*; Ĝ = region minus Ĥ.

## The four sets (unweighted two-sample Cox fit of treatment on the set; Wald 95% CI)

| Panel | Set | n | experimental n / events | control n / events | HR | 95% CI |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | All patients, ITT | 1000 | 500 / 187 | 500 / 218 | 0.808 | 0.665, 0.983 |
| 2 | Region, ITT | 194 | 104 / 22 | 90 / 18 | 1.036 | 0.556, 1.932 |
| 3 | Region minus H* | 159 | 83 / 12 | 76 / 14 | 0.736 | 0.340, 1.592 |
| 4 | Region minus H-hat | 174 | 92 / 16 | 82 / 16 | 0.844 | 0.422, 1.688 |

## Sizes and the 2 × 2 of Ĥ against H*

- Region size 194; |Ĥ ∩ region| = 20; |Ĝ| = 174 (= 194 − 20); |G*| = 159; |H* ∩ region| = 35. Non-region size 806; |Ĥ ∩ non-region| = 73; |H* ∩ non-region| = 130.
- Region: TP 18, FP 2, FN 17, TN 157; Sens 0.514, Spec 0.987.
- Non-region (where the identifier ran): TP 65, FP 8, FN 65, TN 668; Sens 0.500, Spec 0.988.
- All counts and the 2 × 2 equal the committed `estimators.csv` and `hhat_confusion.csv` rows of replicate 5.

## Figure annotations (`km_s2_example.pdf` / `.png`)

- Panel 1, All patients, ITT: HR 0.808 (95% CI 0.665, 0.983); experimental n = 500, events 187; control n = 500, events 218
- Panel 2, Region, ITT: HR 1.036 (95% CI 0.556, 1.932); experimental n = 104, events 22; control n = 90, events 18
- Panel 3, Region minus H*: HR 0.736 (95% CI 0.340, 1.592); experimental n = 83, events 12; control n = 76, events 14
- Panel 4, Region minus H-hat: HR 0.844 (95% CI 0.422, 1.688); experimental n = 92, events 16; control n = 82, events 16

## Trigger and post-conditions

- P1: ITT ub 0.983 < 1.0, region ITT 1.036 > 1.0, Ĥ declared: TRUE.
- P3: panel 3 uses H* = {z_cd40 ≤ 230} (flag_harm equals z_cd40 ≤ 230 for all 1000 patients); panel 4 uses Ĝ = region minus Ĥ; 174 + 20 = 194: TRUE.
- P4: specification block identical to `sims/121_ident_rlog4_full_call.R` (`spec_diff.txt` empty): TRUE; equal to the stored call meta: TRUE.
- P5: hr_itt_3dp TRUE, hr_region_3dp TRUE, hr_region_Ghat_3dp TRUE, Hhat_definition TRUE, declared TRUE.
