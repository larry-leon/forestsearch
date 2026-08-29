# REPORT — the realized selection distribution: `sg_def` tabulated, evaluated on the population, and compared against the analytic

**Task:** `dev/tasks/cc_task_oc_wrapper_sgdef_2026-08-29.md` (commit `69ccc84d`)
**Repo / branch:** `~/Documents/GitHub/forestsearch` (pop-os), `feature/glm-extension`
**Follows:** `HANDOFF_oc_wrapper_2026-08-29_v1.md` §5 (b), (c) and §9; `REPORT_oc_wrapper_null_2026-08-29.md`.
**Category:** read-only. No `R/` change, no export, no edit to any package file or document. Artefacts: `sgdef_selection_2026-08-29.R`, its `.log`, and `sgdef_selection_2026-08-29.rds` (per-rule and per-replicate tables), all under `dev/glm-continuous-sims/`. No `fs_oc_*` re-runs; the analytic side is read from `oc_wrapper_grid_2026-08-29.rds` (alt) and `oc_wrapper_null_2026-08-29.rds` (null).

---

## 1. Provenance (§1 — GATE passed)

```
pop-os · /home/larryleon/Documents/GitHub/forestsearch · feature/glm-extension
HEAD at start 648c2af8 (descends from itself) · git status --porcelain: empty
git log -4: 648c2af8 e3725912 6574beb5 399ccabe
packageVersion forestsearch 0.2.6
```

`~/Downloads/cc_task_oc_wrapper_sgdef*`: one match, `cc_task_oc_wrapper_sgdef_20260829.md` (hyphens stripped) → `dev/tasks/cc_task_oc_wrapper_sgdef_2026-08-29.md`, committed **`69ccc84d`**.

DGMs rebuilt deterministically (alt: `k_inter = truth$beta_inter`; null: `model = "null"`), gated against the payloads' committed `truth` (alt500 / alt700 `effect_Q`, `prevalence_Q`; null `Q` empty, `effect_Q` NA, `effect_Qc`). Frames from `fs_build_eval_frame(dgm, "continuous")` (= `df_super`). Scale: `tauQc = 26.2552`, `bint = 13.7448`; null common effect `26.2552`.

**One fact about the analytic family that this task surfaced and that bears on §5:** the driver ran with `include_str2 <- TRUE` (`sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_1000.qmd` L267–269, appending `str2` to `analysis_binary_vars`), i.e. **13 confounders**; the wrapper's stored families were enumerated on the 12-variable list without `str2` (`fs_args$confounders.name` in `oc_wrapper_grid_2026-08-29.R`). The realized rules use `str2` in 3–5% of winners. The payload's `meta` does not record the confounder list; the `sg_def` column is what reveals it.

---

## 2. §2 — the realized selection distribution

`results$sg_def` over detected replicates (alt500: 1000/1000; alt700: 999/1000; null500: 998/1000). The search's syntax is `{var <= v}` / `!{var <= v}` for continuous cuts and `{var}` / `!{var}` for binaries, `" & "`-joined.

**(a) Verbatim strings** — heavy fragmentation, as expected from sample-quantile cuts:

| cell | detected | distinct `sg_def` | most frequent string (count) |
|---|---:|---:|---|
| alt500 | 1000 | **922** | `!{age <= 35} & {age <= 39}` (4), `!{age <= 42} & {homo}` (4), `{drugs}` (4) |
| alt700 | 999 | **879** | `!{age <= 37} & {age <= 39}` (7), `!{age <= 33} & {age <= 35}` (6), `{race} & !{age <= 35}` (6) |
| null500 | 998 | **912** | `!{age <= 42} & {homo}` (4), `!{homo} & {karnof <= 90}` (4), `{drugs}` (4) |

No string reaches 1% of replicates; top-15 lists are in the log.

**(b) Covariate signatures** (variables and directions, thresholds stripped; `age>` = a negated `age <=` clause, `homo=0` = `!{homo}`):

| cell | distinct signatures | top 8 (count) |
|---|---:|---|
| alt500 | **158** | `age<= & age>` 36 · `age> & preanti<=` 36 · `age> & cd40<=` 33 · `age> & cd80>` 28 · `age> & cd80<=` 25 · `age> & wtkg<=` 24 · `age> & cd40>` 23 · `preanti> & wtkg>` 20 |
| alt700 | **153** | `age<= & age>` 65 · `preanti<= & preanti>` 57 · `age> & preanti<=` 35 · `age> & wtkg<=` 25 · `age> & cd80>` 23 · `age> & cd80<=` 22 · `age<= & preanti>` 19 · `age<= & wtkg>` 19 |
| null500 | **160** | `age<= & age>` 27 · `age<= & preanti>` 27 · `cd40> & preanti>` 24 · `age> & cd80<=` 23 · `preanti> & wtkg>` 23 · `age> & cd40<=` 21 · `age> & preanti<=` 21 · `cd80> & preanti>` 21 |

**K** (factors in the winner): alt500 1: 4 / 2: 996; alt700 1: 13 / 2: 986; null500 1: 4 / 2: 994 — the winner is a two-factor rule in ≥98.7% of replicates in every cell.

**Fraction of winners using each covariate:**

| cell | age | preanti | wtkg | karnof | cd40 | cd80 | hemo | homo | drugs | race | gender | symptom | str2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| alt500 | 0.440 | 0.310 | 0.238 | 0.068 | 0.259 | 0.281 | 0.000 | 0.095 | 0.015 | 0.089 | 0.037 | 0.059 | 0.044 |
| alt700 | 0.453 | 0.295 | 0.216 | 0.063 | 0.237 | 0.235 | 0.002 | 0.070 | 0.049 | 0.066 | 0.067 | 0.074 | 0.031 |
| null500 | 0.404 | 0.328 | 0.240 | 0.068 | 0.267 | 0.283 | 0.002 | 0.103 | 0.014 | 0.088 | 0.041 | 0.058 | 0.050 |

**Does the realized distribution live on the analytic family's high-`p_sel` covariates?** Yes. The prediction document's two axes (`age`, `preanti`) appear in **67.0% / 67.4% / 65.2%** of winners (alt500 / alt700 / null), and the same signatures head both lists (§5). The inventory's quoted winner, `!{cd40 <= 415} & !{cd80 <= 1022}`, uses neither axis; **that is not typical** — a third of winners use neither, two thirds use at least one, and the single most frequent signature in every cell is an age-band (`age<= & age>`) or an age × prior-ART rule. The null cell's covariate usage is indistinguishable from the alternative's (the candidate set is the same 1601 rules and the winner is the noise maximiser in both).

---

## 3. §3 — every realized rule evaluated on the population

Resolver: the package's own `.fs_resolve_membership(frame, rule)` (`R/betaHhat_truth.R` L75 — the eval-frame path `fs_betaHhat_table()` uses, dispatching to `get_dfpred()` exactly as `fs_build_eval_frame()`'s scoring does; no parser was written). **All distinct rules resolved: 922 / 879 / 912 `ok`, 0 `empty`, 0 `unresolved`** in the three cells — every referenced column is on `df_super`, including `str2`.

Per distinct rule, on `df_super` (5000 rows): `Pg_pop`, `nPg_pop = n·Pg_pop`, `PQg_pop = P(g∩Q)/P(g)`, `sens_pop = P(g∩Q)/P(Q)`, `spec_pop = 1 − P(g∩Qc)/P(Qc)`, `beta_pop = tauQc + bint·PQg_pop` (null: `PQg = 0`, `sens = NA`, `beta = 26.2552`). Saved per rule (`$cells[[cell]]$rules`) and joined onto replicates (`$reps`) in the `.rds`.

One structural fact that the evaluation makes explicit: the payload's `betaHhat_H` is **itself a population functional of the realized rule** (the driver scores `fs_betaHhat` on the eval frame), so `beta_real − beta_pop` is 0 to machine precision in every replicate (max |diff| ≈ 4×10⁻¹⁵). The "effect" row of §4 is therefore an identity check on the resolver, not a comparison; it passes.

---

## 4. §4 — the within-rule comparison (the discriminating test)

Per covariate signature with ≥ 20 detected replicates: realized (from `results`) minus population (from §3, over the same replicates' own rules), paired difference ± SE.

**alt500** (8 signatures, 225 replicates):

| signature | n | size real / pop | Δ size ± SE | PPV real / pop | Δ PPV ± SE | Δ sens ± SE | Δ spec ± SE | β real = pop | naive bias real ± SE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| age<= & age> | 36 | 83.6 / 84.7 | −1.13 ± 1.50 | 0.410 / 0.410 | −0.000 ± 0.005 | −0.004 ± 0.004 | +0.001 ± 0.003 | 31.89 | 68.2 ± 2.9 |
| age> & preanti<= | 36 | 72.1 / 73.1 | −1.02 ± 1.02 | 0.970 / 0.968 | +0.002 ± 0.004 | +0.003 ± 0.005 | +0.001 ± 0.001 | 39.56 | 66.3 ± 2.8 |
| age> & cd40<= | 33 | 71.0 / 70.6 | +0.39 ± 1.36 | 0.677 / 0.672 | +0.005 ± 0.008 | +0.007 ± 0.006 | +0.001 ± 0.002 | 35.48 | 72.0 ± 3.0 |
| age> & cd80> | 28 | 69.3 / 68.9 | +0.41 ± 1.45 | 0.666 / 0.675 | −0.009 ± 0.012 | +0.004 ± 0.006 | −0.002 ± 0.003 | 35.53 | 76.9 ± 3.6 |
| age> & cd80<= | 25 | 69.9 / 70.6 | −0.69 ± 1.42 | 0.701 / 0.698 | +0.003 ± 0.011 | −0.002 ± 0.006 | +0.001 ± 0.003 | 35.85 | 76.9 ± 4.0 |
| age> & wtkg<= | 24 | 72.5 / 71.3 | +1.14 ± 1.66 | 0.575 / 0.579 | −0.004 ± 0.010 | +0.007 ± 0.006 | −0.001 ± 0.003 | 34.21 | 71.6 ± 3.6 |
| age> & cd40> | 23 | 67.8 / 69.8 | −1.98 ± 1.46 | 0.696 / 0.695 | +0.001 ± 0.014 | +0.000 ± 0.006 | +0.003 ± 0.003 | 35.81 | 70.6 ± 4.2 |
| preanti> & wtkg> | 20 | 67.6 / 68.5 | −0.96 ± 1.26 | 0.188 / 0.170 | +0.018 ± 0.009 | +0.005 ± 0.003 | +0.004 ± 0.004 | 28.59 | 80.0 ± 3.6 |
| **weighted, 8 sigs** | 225 | | **−0.48** | | **+0.0015** | **+0.0024** | **+0.0010** | 0 | |
| **all 1000, paired** | 1000 | 72.34 / 72.90 | **−0.56 ± 0.24** | 0.4055 / 0.4008 | **+0.0047 ± 0.0016** | **+0.0013 ± 0.0008** | **+0.0021 ± 0.0006** | 0 | 74.27 ± 0.58 |

**alt700** (6 signatures, 227 replicates): weighted Δ size **+0.55**, Δ PPV +0.0006, Δ sens +0.0029, Δ spec +0.0000; all 999 paired: Δ size **−0.02 ± 0.26**, Δ PPV +0.0008 ± 0.0016, Δ sens +0.0004 ± 0.0007, Δ spec +0.0003 ± 0.0004. The largest single signature, `age<= & age>` (65), has Δ size +2.6 ± 1.0; `age> & preanti<=` (35) −1.5 ± 1.4; the rest within ±1.3 ± 1.5.

**null500** (9 signatures, 207 replicates): weighted Δ size **+0.01**; all 998 paired: Δ size **−0.50 ± 0.24**, Δ spec +0.0010 ± 0.0005; PPV 0 = 0; sensitivity NA (0/0) on both sides.

**The decomposition of the handoff's §5b gaps** — analytic (`"resample"`, stored) → population value of the rules the search actually picked → realized:

| cell | quantity | analytic | pop. of realized rules | realized | **between-rule** (pop − analytic) | **within-rule** (realized − pop) | total |
|---|---|---:|---:|---:|---:|---:|---:|
| alt500 | E\\|Ĥ\\| | 70.88 | 72.90 | 72.34 | **+2.02** | −0.56 | +1.45 |
| alt500 | E[sens] | 0.1638 | 0.1693 | 0.1706 | **+0.0055** | +0.0013 | +0.0068 |
| alt500 | E[PPV] | 0.3994 | 0.4008 | 0.4055 | +0.0014 | **+0.0047** | +0.0061 |
| alt700 | E\\|Ĥ\\| | 73.03 | 73.65 | 73.63 | **+0.62** | −0.02 | +0.60 |
| alt700 | E[sens] | 0.1217 | 0.1224 | 0.1229 | +0.0007 | +0.0004 | +0.0012 |
| alt700 | E[PPV] | 0.4027 | 0.4007 | 0.4016 | −0.0020 | +0.0008 | −0.0011 |
| null500 | E\\|Ĥ\\| | 71.00 | 72.55 | 72.05 | **+1.55** | −0.50 | +1.05 |

### The verdict, in §4's terms

**Realized and population agree within signatures for size, and the size gap is between-rule.** Within every signature with ≥ 20 replicates, in all three cells, the realized `n_harm` sits within ±2 of its rule's population size and the pooled paired difference is **zero or slightly negative** (−0.56 ± 0.24 at n = 500, −0.02 ± 0.26 at n = 700, −0.50 ± 0.24 at the null) — the search does *not* pick rules whose realized sample size runs above their population size; if anything the opposite. The E|Ĥ| gap (+1.45 / +0.60 / +1.05) is therefore **not a within-rule sample-realization effect**. It is **between-rule**: the population size of the rules the search actually selects exceeds the analytic selection-weighted expectation by +2.0 / +0.6 / +1.5 subjects. Re-weighting the *measured* signatures' population sizes by the *analytic* `sel_c` reproduces the analytic figure (71.97 vs 70.88 at n = 500; 72.82 vs 73.03 at n = 700; 71.73 vs 71.00 at the null) — the difference is in the weights, not in the rules' sizes. The handoff's §5b attribution of E|Ĥ| to population-versus-sample is **not supported**; the handoff itself flagged that mechanism as "less crisp" and to be checked, and the check goes the other way.

**Sensitivity is likewise predominantly between-rule at n = 500** (+0.0055 of the +0.0068), with a small within-rule component (+0.0013 ± 0.0008, 1.6 SE); at n = 700 the whole gap is 0.0012 and neither component is distinguishable from zero.

**PPV is the mixed case, and the one where the handoff's explanation holds:** at n = 500 the within-rule component (+0.0047 ± 0.0016, 2.9 SE) is three times the between-rule component (+0.0014) — realized sample purity at the selected rule *does* exceed its population purity, in the direction the prediction document records (L2283–2292). At n = 700 both components are below 0.001. Specificity shows the same small within-rule excess (+0.0021 ± 0.0006 at n = 500; +0.0010 ± 0.0005 at the null), consistent with the same mechanism (the selected rule's realized composition is slightly more favourable than its population composition), and it is what the handoff already placed in category (a).

So: **a mixed result, not one story.** Within-rule sample realization is real but small, and visible only on the composition metrics (PPV, specificity, marginally sensitivity) at n = 500; it does **not** explain the size gap, which is a between-rule effect — the analytic selection distribution puts less weight on larger rules than the search does. Magnitudes: between-rule +2.0 subjects at n = 500 (2.8% of E|Ĥ|), within-rule PPV +0.005 (1.2%).

**The naive bias (handoff §5c).** Realized aggregate vs analytic (`"resample"`): **74.27 ± 0.58 vs 75.44 (−1.17)** at n = 500; **77.09 ± 0.56 vs 75.53 (+1.55)** at n = 700; **74.85 ± 0.57 vs 76.22 (−1.37)** at the null. Within signatures the realized naive bias varies far more than the aggregate discrepancy — from 66 to 80 at n = 500 (SE ≈ 3 each), 70 to 79 at n = 700 — and it is strongly and negatively correlated with the signature's population size and purity (corr with `size_pop` −0.59 / −0.80 / −0.36; with `PQg_pop` −0.54 / −0.95 in the two alternative cells): small, impure winners carry the largest winner's curse, as the joint-normal law predicts. The per-candidate analytic naive bias is not stored (only the aggregate `Enaive_bias`), so a within-signature analytic-vs-realized comparison is not available from the `.rds`; what *is* established is that the aggregate is a mixture over signatures whose components differ by ±5–9 around it, so **a shift of a few percent in the signature mix moves the aggregate by ≈1 — the size of the discrepancy — and the between-rule weight differences documented above are of exactly that order.** At n = 700 the 13 one-factor winners alone carry a naive bias of 86.8 against 77.0 for two-factor winners. The sign flip across n is therefore *consistent with* a mix effect and *not* with a systematic definitional cause, but this is not yet a demonstration: it would need the analytic per-candidate naive bias, which the wrapper computes internally and does not return.

---

## 5. §5 — measured vs analytic selection distribution over signatures

Analytic `sel_c` at the driver's `(30, 10)`, both gates, aggregated by signature (labels mapped by the same clause parser: `var <= v` → `var<=`, `var > v` → `var>`, bare binary → `=1`, `var != 1` → `=0`; **no ambiguous labels** — every analytic candidate maps to exactly one signature). Analytic signatures with positive mass: 183 / 181 / 184; measured distinct: 158 / 153 / 160.

**Top analytic signatures beside measured frequency** (alt500; `"split"` is within 0.001 of `"resample"` on every row and is omitted here — the full table is in the log and `.rds`):

| signature | analytic | measured | n | | signature | analytic | measured | n |
|---|---:|---:|---:|---|---|---:|---:|---:|
| age> & preanti<= | 0.039 | 0.036 | 36 | | age> & wtkg<= | 0.020 | 0.024 | 24 |
| age> & cd40> | 0.032 | 0.023 | 23 | | age<= & preanti> | 0.019 | 0.019 | 19 |
| age> & cd80<= | 0.032 | 0.025 | 25 | | cd80<= & preanti> | 0.018 | 0.015 | 15 |
| age<= & age> | 0.031 | 0.036 | 36 | | cd40> & preanti> | 0.018 | 0.018 | 18 |
| age> & cd40<= | 0.029 | 0.033 | 33 | | cd40> & wtkg> | 0.015 | 0.008 | 8 |
| age> & cd80> | 0.023 | 0.028 | 28 | | cd80> & preanti> | 0.015 | 0.009 | 9 |
| age> & wtkg> | 0.023 | 0.019 | 19 | | age<= & wtkg> | 0.015 | 0.011 | 11 |
| preanti<= & preanti> | 0.021 | 0.012 | 12 | | | | | |

alt700: the top three coincide exactly (`age<= & age>` 0.074 vs 0.065; `preanti<= & preanti>` 0.057 vs 0.057; `age> & preanti<=` 0.033 vs 0.035). null500: `age<= & age>` 0.032 vs 0.027; `age<= & preanti>` 0.026 vs 0.027; `preanti<= & preanti>` is the one visible outlier (0.019 vs 0.008).

| cell | analytic mass on signatures **never selected** (resample / split) | measured frequency on signatures **absent from the analytic family** | of which `str2` | other absent signatures | TVD (resample / split) | TVD excl. `str2`, renormalised | multinomial noise floor of the TVD at this replicate count (mean, 95%) |
|---|---:|---:|---:|---|---:|---:|---:|
| alt500 | 0.0139 / 0.0144 | 0.0460 (21 sigs) | 0.0440 (19) | `drugs=1 & symptom=0`, `karnof<= & race=1` (1 each) | 0.189 / 0.191 | **0.172** | 0.136 (0.119–0.157) |
| alt700 | 0.0147 / 0.0141 | 0.0430 (21) | 0.0310 (15) | `age>` (3), `drugs=1 & homo=0` (3), `gender=0 & wtkg>` (2), `karnof<= & symptom=1` (2), 2 singles | 0.174 / 0.173 | **0.160** | 0.132 (0.115–0.149) |
| null500 | 0.0190 / 0.0194 | 0.0521 (23) | 0.0501 (21) | `drugs=1 & symptom=0`, `karnof<= & race=1` (1 each) | 0.198 / 0.198 | **0.178** | 0.138 (0.122–0.158) |

**Why the absent signatures are absent.** (i) `str2` — 3–5% of measured mass — because the wrapper's family was enumerated without `str2` while the driver included it (§1); a family difference, and a known one now. (ii) The handful of non-`str2` absentees (0.2–1.2%) are single-factor binaries or pairs that the population enumeration drops at a floor (`age>` alone at n = 700: a one-factor rule whose population prevalence sits at the `n.min/n` boundary; `drugs=1 & symptom=0` etc.: pairs below `Pg ≥ 0.12` on `df_super` but above 60 subjects in some sample). These are floor-boundary cases, not representability failures. Conversely the 1.4–1.9% analytic mass the search never selected is spread over many low-probability candidates (none above 0.003), consistent with 1000 replicates simply not sampling them.

**Overlap.** Total variation distance 0.17–0.20 over all signatures; 0.16–0.18 once the `str2` signatures are removed from the measured side and both renormalised. The multinomial sampling floor — the TVD one expects between 1000 draws from the analytic distribution and that distribution itself — is 0.13–0.14 (95% up to 0.16). The excess over the floor, **≈ 0.03–0.04**, is the genuine between-rule difference, and it is not spread uniformly: the search over-selects the age-band and age × cd40/cd80/wtkg signatures and under-selects `preanti<= & preanti>` and several `cd40`/`cd80` × `wtkg`/`preanti` pairs, at the alternative and the null alike.

---

## 6. What this settles, and what it leaves open

**Settled (facts).**
- The realized winners live on the same covariates and the same top signatures as the analytic family; the inventory's quoted winner is atypical; two thirds of winners use `age` or `preanti`.
- Every realized rule resolves on `df_super` with the package's own resolver; `betaHhat_H` is already a population functional of the realized rule (identity to 10⁻¹⁵).
- Handoff §5b: the **E|Ĥ| gap is between-rule** (+2.0 / +0.6 / +1.5 subjects), not within-rule — realized size does not exceed population size within signatures (it is 0.5 *below* at n = 500 and at the null). The population-versus-sample explanation does not carry E|Ĥ|.
- Handoff §5b: **sensitivity** is predominantly between-rule at n = 500; **PPV and specificity** carry a small genuine within-rule excess at n = 500 (+0.005, +0.002; 3 SE), which is the population-versus-sample effect the inventory's B5 row 4 predicted — real, but 1% in size.
- The measured and analytic selection distributions differ by a TVD of ≈ 0.03–0.04 above the sampling floor, with 3–5% of measured mass on `str2` rules the analytic family does not contain because the wrapper's confounder list omitted `str2`.

**Open (facts stated, no recommendation).**
- Whether the between-rule weight difference (search favouring larger, age-band rules) is a family effect (population vs sample quantile cuts; the missing `str2`; the floor-boundary absentees) or an argmax effect (statistics-keyed near-duplicate removal, pre-screen, or the `rmin` translation) is not separated by this task; both act on the same signatures.
- Handoff §5c: the naive-bias discrepancy is *consistent with* a signature-mix effect — within-signature naive bias spans ±5–9 around the aggregate and is strongly tied to the winner's population size/purity — but the per-candidate analytic naive bias is not stored, so the within-signature analytic comparison was not made.
- All of the above is one DGM, two sample sizes, one outcome type.

---

## 7. Close-out

`git status --porcelain` before staging: the four new files only. Staged by explicit path: `sgdef_selection_2026-08-29.R`, `.log`, `.rds`, this report. Commit hash in the follow-up line below. **No push. No install** (nothing in `R/` changed).

Commits: `69ccc84d` task doc; __COMMITS__.

---

## 8. Verdict (ten lines)

1. All 2,713 distinct realized rules across the three cells resolve on `df_super` with the package's own resolver; none failed.
2. The realized winners are two-factor rules on the analytic family's own top signatures; `age`/`preanti` appear in two thirds of them — the inventory's quoted winner is atypical.
3. **Size is between-rule, not within-rule:** within signatures realized `n_harm` matches population size to ±2 (pooled −0.5 ± 0.2), while the rules the search picks are +2.0 / +0.6 / +1.5 subjects larger in the population than the analytic argmax's selection-weighted expectation.
4. **PPV (and specificity) are within-rule at n = 500:** +0.005 ± 0.002 sample-over-population purity at the selected rule — the inventory's B5 mechanism, real but ~1%.
5. Sensitivity sits between: mostly between-rule at n = 500, indistinguishable from zero at n = 700.
6. Hence handoff §5b is **partly right and partly wrong**: right for PPV/specificity, wrong for E|Ĥ|; the evidence is mixed by quantity and it should be recorded that way.
7. The measured-vs-analytic selection distributions have TVD 0.17–0.20, of which 0.13–0.14 is 1000-replicate sampling noise and 3–5% is `str2` — the driver ran 13 confounders, the wrapper's family 12; the residual ≈ 0.03–0.04 is a genuine tilt toward larger, age-band rules.
8. Handoff §5c stays open: within-signature naive bias varies ±5–9 around the aggregate and tracks winner size/purity, so the ±1.2–1.6 aggregate discrepancy and its sign flip are consistent with a mix effect, unproven without the per-candidate analytic term.
9. `betaHhat_H` in the payloads is already the population effect of the realized rule; the "effect" row of §4 is an identity, not a comparison.
10. No `R/` change, no re-runs, one script, one `.rds`, one report; one DGM, two n, one outcome type.
