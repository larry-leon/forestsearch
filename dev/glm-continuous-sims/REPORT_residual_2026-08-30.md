# REPORT — the two candidate explanations for the remaining residual: `se_g`'s prevalence scaling (A) and the dedup asymmetry (B)

**Task:** `dev/tasks/cc_task_residual_2026-08-30.md` (commit `80544be1`)
**Repo / branch:** `~/Documents/GitHub/forestsearch` (pop-os), `feature/glm-extension`
**Follows:** `REPORT_oc_wrapper_sgdef_2026-08-29.md` (the between-rule finding), `REPORT_oc_wrapper_confs_2026-08-30.md` (the corrected 13-variable family). Both **stand**; their `.rds` are read, not redone.
**Category:** no `R/` change, no export, no edit to any package file, driver or document. Artefacts, all under `dev/glm-continuous-sims/`: `residual_2026-08-30.R` (parts `A`, `B <cell>`, `Bsum`), `residual_hypA_2026-08-30.{rds,log}`, `residual_hypB_{alt500,alt700,null500}_2026-08-30.log`, `residual_hypB_summary_2026-08-30.log`, `residual_hypB_2026-08-30.rds`, this report. No simulations, no replicates, no renders, no push, no install.

---

## 1. Provenance (§1 — GATE passed)

```
pop-os · /home/larryleon/Documents/GitHub/forestsearch · feature/glm-extension
HEAD at start 5ff1487a · git status --porcelain: empty
git log -4: 5ff1487a 202ef15b c946def6 e5b15126
packageVersion forestsearch 0.3.0 (DESCRIPTION 0.3.0); vi.grf.min default NULL
```

`~/Downloads/cc_task_residual_2026-08-30.md` was not present under that name; the stem arrives with hyphens stripped as **`cc_task_residual_20260830.md`** (one match) → `dev/tasks/cc_task_residual_2026-08-30.md`, committed **`80544be1`**.

DGMs rebuilt deterministically as the 08-29 / 08-30 scripts build them and gated against the payloads' `truth` (alt: `effect_Q`, `prevalence_Q`; null: `Q` empty, `effect_Qc`). The corrected families were **re-enumerated in every part and asserted `identical()`/`all.equal()` to the stored record** (`lab`, `Pg`, `PQg`, `beta_g`, `se_g`, `M` = 1696 / 1890 / 1696). Realized-rule memberships were re-resolved with `.fs_resolve_membership()` and asserted equal to the stored `Pg_pop` to 10⁻¹².

Compute: part A 4 min; part B 4 `fs_oc_predict()` runs per cell (2 analogues × 2 gates, seed 20260825, 2×10⁵ draws), cells run as three concurrent processes, longest 50 min; `Bsum` seconds. (An operator note: the first B launch died at the family-enumeration stage when the launching shell was torn down, and the second lost the alt cells' last gate to the same cause; the third launch, detached, ran to completion. Every run is deterministic from stored inputs and the same seed — the partial runs' lines, where they existed, agree with the final logs to every printed digit.)

---

## 2. Hypothesis A — is `se_g`'s prevalence scaling wrong for candidates unlike Q?

### 2.1 What the ratio is, algebraically

The wrapper (`R/fs_oc_family.R` L256, L404) sets `seQ1000 = sqrt(V_eff[Q] / (1000·πQ))` and `se_g = seQ1000 · sqrt(1000/n) · sqrt(πQ/Pg)` (the task's `sqrt(2)` is `sqrt(1000/n)` at n = 500). That collapses to

`se_scaled = sqrt(V_eff[Q] / (n·Pg))`, while `fs_dgm_scale(dgm, regions = g)` gives `se_direct = sqrt(V_eff[g] / (n·Pg))` (`R/fs_dgm_scale.R` L234),

so **`ratio = se_scaled / se_direct = sqrt(V_eff[Q] / V_eff[g])`** — n-free, and the same number at n = 500 and n = 700 (asserted in the script). The hypothesis is therefore exactly: *does the per-subject effective variance `V_eff[g] = 2·(V_arm0 + V_arm1)` with `V_arm_w = σ² + V_g[μ_w]` drift with the candidate?* On this fixture σ² = 16 256 and the whole-population `V[μ0]` is 729 (`fs_dgm_scale` S row), so the within-region term can move `V_eff` by at most ≈ 4 % and the ratio by at most ≈ 2 % in either direction; that bound is what the measurement fills.

Null cell: the reference row is S (`V_eff[S]` = 67 940, `πQ` = 1) and `ratio = sqrt(V_eff[S] / V_eff[g])`.

### 2.2 The distribution — every distinct realized rule (922 / 879 / 912), and the whole corrected family (1696 / 1890 / 1696)

| cell | population | min | 1 % | 5 % | 25 % | **median** | 75 % | 95 % | 99 % | max | selection-weighted mean |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| alt500 | realized rules | 0.9924 | 0.9936 | 0.9950 | 0.9976 | **1.0007** | 1.0042 | 1.0082 | 1.0098 | 1.0119 | 1.0009 (by replicate count) |
| alt500 | analytic family | 0.9923 | 0.9941 | 0.9954 | 0.9972 | **0.9987** | 1.0024 | 1.0068 | 1.0090 | 1.0106 | 1.0007 (`sel_c`) |
| alt700 | realized rules | 0.9919 | 0.9926 | 0.9945 | 0.9969 | **1.0002** | 1.0042 | 1.0088 | 1.0107 | 1.0127 | 1.0005 |
| alt700 | analytic family | 0.9922 | 0.9936 | 0.9950 | 0.9972 | **0.9988** | 1.0026 | 1.0069 | 1.0092 | 1.0120 | 1.0006 |
| null500 | realized rules | 0.9951 | 0.9962 | 0.9979 | 1.0008 | **1.0040** | 1.0074 | 1.0119 | 1.0132 | 1.0150 | 1.0040 |
| null500 | analytic family | 0.9957 | 0.9972 | 0.9984 | 1.0004 | **1.0022** | 1.0057 | 1.0099 | 1.0125 | 1.0142 | 1.0038 |

Everything lies in **[0.992, 1.015]**; the selection-weighted mean is within 0.1 % of 1 at both alternatives and 0.4 % at the null.

### 2.3 The trend — with prevalence `Pg` and with composition distance from Q (purity `PQg`; Jaccard with Q also computed)

| cell | population | cor(ratio, Pg) | cor(ratio, PQg) | cor(ratio, Jaccard) | lm slope on Pg (SE) | lm slope on PQg (SE) | mean ratio by Pg quartile Q1 → Q4 | **large & low-purity** (Pg ≥ Q3, PQg ≤ 0.4): n, mean, min |
|---|---|---:|---:|---:|---:|---:|---|---|
| alt500 | realized | **−0.003** | +0.124 | +0.092 | −0.0003 (0.0042) | +0.0021 (0.0005) | 1.0009 · 1.0008 · 1.0017 · 1.0011 | 116, **1.0006**, 0.9924 |
| alt500 | family | −0.232 | +0.149 | −0.000 | −0.0053 (0.0005) | +0.0029 (0.0004) | 1.0005 · 1.0004 · 1.0001 · 0.9984 | 296, **0.9984**, 0.9952 |
| alt700 | realized | **+0.009** | +0.070 | +0.044 | +0.0017 (0.0068) | +0.0012 (0.0006) | 1.0010 · 1.0006 · 1.0005 · 1.0010 | 126, **1.0009**, 0.9937 |
| alt700 | family | −0.202 | +0.140 | +0.001 | −0.0045 (0.0005) | +0.0027 (0.0004) | 1.0004 · 1.0006 · 0.9999 · 0.9987 | 333, **0.9988**, 0.9952 |
| null500 | realized | +0.037 | — | — | +0.0050 (0.0044) | — | 1.0036 · 1.0043 · 1.0048 · 1.0042 | — |
| null500 | family | −0.208 | — | — | −0.0047 (0.0005) | — | 1.0037 · 1.0036 · 1.0034 · 1.0018 | — |

By purity band (alt500 realized / family): ≤ 0.2: 0.9990 / 0.9983; 0.2–0.4: 1.0019 / 1.0003; 0.4–0.6: 1.0019 / 1.0002; 0.6–0.8: 1.0015 / 1.0001; 0.8–0.95: 0.9981 / 0.9984; > 0.95: 0.9997 / 1.0002 (alt700 the same to the third decimal).

Reading: over the rules the search actually selects there is **no prevalence trend** (|cor| ≤ 0.04, slopes indistinguishable from zero); over the whole family there is a statistically detectable but numerically negligible one — the top prevalence quartile (Pg 0.43–0.92, candidates the search never selects) sits 0.16 % low. The purity slope is +0.002–0.003 per unit of PQg: impure candidates have a ratio ≈ 0.1–0.2 % *below* the pure ones, i.e. the wrapper's `se_g` is *smaller* than the direct value for them by 0.1–0.2 % — the sign the hypothesis needs, at a magnitude three orders below anything that could move the argmax. The "large, low-purity" cell of the hypothesis is at 1.0006 / 1.0009 on the realized rules and 0.9984 / 0.9988 on the family.

### 2.4 Controls

| control | alt500 | alt700 |
|---|---|---|
| Q itself | 1.000000 | 1.000000 (by construction) |
| family candidates with PQg ≥ 0.95 | n = 15, mean 1.0002, range [0.9995, 1.0007] | n = 20, mean 0.9999, range [0.9980, 1.0007] |
| realized rules with PQg ≥ 0.95 | n = 29, mean 0.9997, range [0.9979, 1.0006] | n = 34, mean 0.9990, range [0.9975, 1.0005] |

The controls agree with the body of the distribution to the third decimal: the highest-purity candidates are not distinguishable from the impure ones on this ratio.

### 2.5 A side column, outside the task's definition

`fs_dgm_scale()` idealises the region size and arm split as fixed (its roxygen "Idealisations"); `fs_scale_se(jensen = TRUE)` inflates for the random counts. Against the Jensen-inflated direct sd the wrapper's `se_g` is **0.97–1.00** (median 0.987 / 0.987 / 0.990; cor with Pg **+0.55 / +0.53 / +0.55**) — the wrapper understates the unconditional sd by ~1.3 % on average and by ~2.5 % for the *smallest* candidates. That is a real, prevalence-trending omission, but (i) it is not the hypothesis as posed (the hypothesis is about `V_eff`, not the count randomness), (ii) its sign is *opposite* to what the size gap needs — small candidates are the ones under-dispersed, so correcting it would raise the noise on *small* rules and, if anything, tilt the analytic argmax toward *smaller* winners — and (iii) it is 1–2.5 %, on a scale where the gap is a few subjects out of 70. Recorded as a fact; no claim is made for it.

### 2.6 Verdict on A

**Refuted.** The ratio is `sqrt(V_eff[Q]/V_eff[g])`, sits in [0.992, 1.015] over every realized rule and every family member in all three cells, has no prevalence trend on the realized rules, and is 1.0006–1.0009 exactly where the hypothesis predicted it to be below 1. Per-subject effective variance is candidate-invariant to within 2 % on this fixture because σ² = 16 256 dominates the ≤ 750 of within-region mean variance; the handoff's "0.86 % spread" was the right reading. `se_g` does not understate the SEs of large, low-purity candidates.

---

## 3. Hypothesis B — does the dedup asymmetry account for the gap?

### 3.1 The key, from source

`remove_near_duplicate_subgroups()` (`R/subgroup_consistency_helpers.R` L1102–1131):

```r
cols_to_check <- 2:min(10, ncol(df))                          # L1110
for (i in cols_to_check) if (is.numeric(df_rounded[, i]))
  df_rounded[, i] <- round(df_rounded[, i] / tolerance) * tolerance   # L1114-1116, tolerance = 0.001
dup_key <- apply(key_cols, 1, function(x) paste(x, collapse = "_"))   # L1120
keep_rows <- !duplicated(dup_key)                             # L1122
```

Columns 2:10 of `found.hrs` are, by `format_search_results()` (`R/subgroup_search.R` L878–879): **`K, n, E, d1, m1, m0, HR, L(HR), U(HR)`**. On the GLM / continuous path (`evaluate_combination_with_status()` L640–647; `fit_glm_for_subgroup()` L815–820): `E = n0 + n1 = n`, `d1 = n1` (the treated count), **`m1 = m0 = NA`**, `HR` = the MD estimate, `L/U = estimate ∓ 1.96·se`. So on this path the key is **`(K, n_g, n1_g, β̂_g, ŝe_g)` to 0.001**, and `found.hrs` is ordered `(-HR, K)` (L883, stable) before the dedup, so "first" = highest estimate → fewest factors → enumeration order. It runs at `subgroup_consistency_main.R` L579 for every focus but `maxeff`, on the candidates that passed the effect screen (`HR >= hr.threshold`, L530–533), and has no off switch.

**What the key does in a trial.** Two candidates with *different* sample memberships have different `n` or, if the same `n`, different `n1` or a different `β̂` (a mean of a different set of outcomes; equality to 0.001 has probability ≈ 0). So the sample key collapses, in practice, exactly the candidates whose **sample memberships coincide within K** — near-twins in the population that happen to draw the same subjects — and nothing else. That is the mechanism the analytic family lacks: it collapses identical *population* memberships (189 / 208 / 189 of them), which never coincide in a sample unless they are identical.

### 3.2 The two analogues built

**v1 — literal.** Every key column's population counterpart at the same tolerance and rounding code: `(K, n·Pg, n·Pg, n·Pg/2, NA, NA, β_g, β_g − 1.96·se_g, β_g + 1.96·se_g)`. Because `Pg = |g|/5000`, `β_g = τQc + b·|g∩Q|/|g|` and `se_g` is a function of `Pg`, matching to 0.001 is matching **`(K, |g|, |g∩Q|)` exactly** (null: `(K, |g|)`) — asserted in the script. Representative: first after ordering `(−β_g, K, enumeration index)`, as the package orders.

**v2 — near-twin (labelled approximate).** The population stand-in for "sample memberships coincide": collapse, within K, pairs whose super-population symmetric difference `d` satisfies `(1 − d/5000)^n ≥ ½`, i.e. `d ≤ 6` at n = 500 and `d ≤ 4` at n = 700, keeping the higher `β_g` (greedy in the package's order). This replaces a per-replicate random event by a deterministic threshold; the pairs it removes coincide in roughly half of trials and survive as two candidates in the other half.

| cell | M | **v1** groups / removed → M | Jaccard(dropped, representative): median · mean · max · share ≥ 0.99 | **v2** `d_max` / pairs / removed → M | min symdiff between any two candidates | pairs with P(coincide) ≥ 0.05 |
|---|---:|---|---|---|---:|---:|
| alt500 | 1696 | 11 / **12** → 1684 | 0.25 · 0.23 · 0.44 · **0.00** | 6 / 10 / **10** → 1686 | 3 | 102 (d ≤ 29) |
| alt700 | 1890 | 17 / **19** → 1871 | 0.21 · 0.22 · 0.45 · **0.00** | 4 / 7 / **7** → 1883 | 3 | 78 (d ≤ 21) |
| null500 | 1696 | 348 / **458** → 1238 | 0.14 · 0.18 · 0.81 · **0.00** | 6 / 10 / **10** → 1686 | 3 | 102 |

**Fidelity.** The literal key does **not** reproduce the sample mechanism: the candidates it collapses are *unrelated* rules that happen to share `(K, |g|, |g∩Q|)` — `age <= 28 & wtkg > 67` with `age <= 30 & wtkg <= 67` (Jaccard 0.000), `age > 28 & homo != 1` with `wtkg > 75 & cd40 <= 349` (0.126) — which in a trial have different `β̂` and are never merged. At the null, where `β_g` is constant and the key reduces to `(K, |g|)`, it removes 27 % of the family (458), carrying **42 % of the stored selection mass**; that is an over-collapse with no counterpart in the search, and its results (below) are reported as an artefact of the key, not as a test. The near-twin version removes what the sample key actually removes: the `str2` / `preanti` twins (`age > 37 & str2` vs `age > 37 & preanti > 0`, 3 subjects apart — `str2` is a coarsening of `preanti`, confs report §3.6), 7–10 pairs carrying 1.1–1.4 % of the stored selection mass, at most half of which the search collapses in any given replicate. **Both are approximations; v2 is the closer one, and it is an upper bound on the mechanism's effect** (it removes every near-twin, the search only the half that coincide).

### 3.3 The re-runs — both gates, seed 20260825, 2×10⁵ draws, `(c1, c2) = (30, 10)`

MC SE of `EnH` from the `sel_c`-weighted spread of `n·Pg`: 0.028 (n = 500), 0.031 (n = 700). "between" = population size of the realized rules (72.90 / 73.65 / 72.55, sgdef §4) − analytic `EnH`; "re-weighted" = measured signatures' population sizes weighted by the analytic `sel_c` (the sgdef check).

**alt500** (measured: E|Ĥ| 72.34 ± 0.41, sens 0.1706, PPV 0.4055, naive bias 74.27 ± 0.58):

| family | gate | M | det | E\|Ĥ\| | **between** | re-weighted | E[sens] | E[PPV] | naive bias | gap to measured: E\|Ĥ\| · sens · PPV · naive |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| corrected (stored) | resample | 1696 | 0.99906 | 70.786 | **+2.111** | 71.88 | 0.1638 | 0.4000 | 75.53 | +1.55 · +0.0068 · +0.0054 · −1.26 |
| v1 literal | resample | 1684 | 0.99900 | 70.821 | +2.077 | 71.90 | 0.1639 | 0.4000 | 75.47 | +1.51 · +0.0067 · +0.0054 · −1.19 |
| **v2 near-twin** | resample | 1686 | 0.99894 | 70.841 | **+2.057** | 71.93 | 0.1635 | 0.3991 | 75.45 | +1.49 · +0.0070 · +0.0064 · −1.18 |
| corrected (stored) | split | 1696 | 1.00000 | 70.781 | +2.116 | 71.83 | 0.1635 | 0.3994 | 75.38 | +1.55 · +0.0070 · +0.0061 · −1.11 |
| v1 literal | split | 1684 | 1.00000 | 70.805 | +2.092 | 71.93 | 0.1638 | 0.4001 | 75.40 | +1.53 · +0.0068 · +0.0054 · −1.13 |
| v2 near-twin | split | 1686 | 0.99999 | 70.846 | +2.052 | 71.95 | 0.1638 | 0.3998 | 75.37 | +1.49 · +0.0067 · +0.0057 · −1.10 |

**alt700** (measured: 73.63 ± 0.42, 0.1229, 0.4016, 77.09 ± 0.56):

| family | gate | M | det | E\|Ĥ\| | **between** | re-weighted | E[sens] | E[PPV] | naive bias | gap: E\|Ĥ\| · sens · PPV · naive |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| corrected (stored) | resample | 1890 | 0.99976 | 73.043 | **+0.608** | 72.68 | 0.1219 | 0.4034 | 75.57 | +0.58 · +0.0009 · −0.0019 · +1.52 |
| v1 literal | resample | 1871 | 0.99982 | 73.218 | +0.433 | 72.69 | 0.1235 | 0.4084 | 75.44 | +0.41 · −0.0006 · −0.0068 · +1.65 |
| **v2 near-twin** | resample | 1883 | 0.99987 | 73.073 | **+0.578** | 72.70 | 0.1221 | 0.4036 | 75.58 | +0.55 · +0.0008 · −0.0020 · +1.50 |
| corrected (stored) | split | 1890 | 1.00000 | 73.083 | +0.568 | 72.73 | 0.1222 | 0.4038 | 75.49 | +0.54 · +0.0007 · −0.0023 · +1.59 |
| v1 literal | split | 1871 | 1.00000 | 73.274 | +0.377 | 72.74 | 0.1237 | 0.4087 | 75.40 | +0.35 · −0.0009 · −0.0071 · +1.69 |
| v2 near-twin | split | 1883 | 1.00000 | 73.073 | +0.578 | 72.67 | 0.1220 | 0.4033 | 75.57 | +0.55 · +0.0009 · −0.0017 · +1.52 |

**null500** (measured: 72.05 ± 0.40, naive bias 74.85 ± 0.57):

| family | gate | M | false decl. | E\|Ĥ\| | **between** | re-weighted | naive bias | gap: E\|Ĥ\| · naive |
|---|---|---:|---:|---:|---:|---:|---:|---|
| corrected (stored) | resample | 1696 | 0.99694 | 70.902 | **+1.649** | 71.60 | 76.32 | +1.15 · −1.47 |
| v1 literal (over-collapse) | resample | 1238 | 0.99494 | 73.133 | −0.583 | 72.10 | 72.73 | −1.09 · +2.12 |
| **v2 near-twin** | resample | 1686 | 0.99697 | 70.958 | **+1.592** | 71.68 | 76.26 | +1.09 · −1.41 |
| corrected (stored) | split | 1696 | 0.99999 | 70.829 | +1.721 | 71.60 | 76.07 | +1.22 · −1.22 |
| v1 literal (over-collapse) | split | 1238 | 0.99999 | 72.989 | −0.438 | 72.09 | 72.47 | −0.94 · +2.38 |
| v2 near-twin | split | 1686 | 0.99997 | 70.899 | +1.651 | 71.65 | 76.10 | +1.15 · −1.25 |

### 3.4 Reading

- **M before → after:** 1696 → 1686 / 1890 → 1883 / 1696 → 1686 (v2); 1684 / 1871 / 1238 (v1).
- **The between-rule gap (v2, the closer analogue):** +2.11 → **+2.06**, +0.61 → **+0.58**, +1.65 → **+1.59** — moves by 0.03–0.06 subject, ≈ 2 MC SE, i.e. **2–5 % of the gap**, in the right direction. Since v2 removes every near-twin while the search removes about half of them per replicate, the mechanism's expected contribution is about **half of that: ≈ 0.02–0.03 subject, 1–3 % of the gap**. The re-weighted check moves by ≤ 0.08 subject and stays where it was.
- **The literal key at the alternatives** (12 / 19 unrelated candidates): +2.11 → +2.08 and +0.61 → +0.43 (n = 700 moves by 0.18 because two of its 19 removed candidates carried 2.3 % of the mass; PPV moves +0.005 in the *wrong* direction there). At the null its −0.58 is the artefact of removing 27 % of the family — the analytic E|Ĥ| jumps 2.2 subjects *past* the measurement and the naive bias falls 3.6 units *past* it. That result quantifies what a real quarter-of-the-family collapse would do, not what the dedup does.
- **E|Ĥ|, PPV, sensitivity, naive bias against measurement:** the E|Ĥ| gap +1.55 / +0.58 / +1.15 → +1.49 / +0.55 / +1.09 (v2); sensitivity and PPV move by ≤ 0.001 at both alternatives (v2); naive bias by ≤ 0.08 (v2), the n = 500 / n = 700 sign flip untouched (−1.18 / +1.50 / −1.41). **Direction: toward measurement. Magnitude: a few percent of each residual, within two MC SE of no movement.**

### 3.5 Verdict on B

**It does not account for the gap.** The statistics-keyed dedup is, on this path, an identical-sample-membership collapse within K; the population candidates it can touch are the 7–10 `str2`/`preanti` near-twin pairs (1.1–1.4 % of selection mass), and removing all of them moves the between-rule gap by 0.03–0.06 subject of the +2.11 / +0.61 / +1.65 — **2–5 % as an upper bound, ≈ 1–3 % expected**. The literal population key is not a faithful analogue (it merges unrelated candidates and, at the null, over-collapses a quarter of the family); its alt-cell results are consistent with v2's and its null result is discarded as an artefact.

---

## 4. Verdict

1. **Hypothesis A — refuted.** `se_scaled/se_direct = sqrt(V_eff[Q]/V_eff[g])` ∈ [0.992, 1.015] over 2 713 realized rules and 5 282 family members; median 1.000–1.004; on the realized rules cor(ratio, Pg) = −0.003 / +0.009 / +0.037 and the large-low-purity cell is at 1.0006 / 1.0009; the controls (Q = 1.000000, purity ≥ 0.95 at 0.999–1.000) agree. The prevalence scaling is right to within 2 %, with no systematic drift where the hypothesis needed one.
2. **Hypothesis B — accounts for ≈ 0.03–0.06 of +2.11, ≈ 0.03 of +0.61, ≈ 0.06 of +1.65** as an upper bound (2–5 %), about half that in expectation. Direction correct, magnitude negligible.
3. **Neither explains the residual.** Both refuted at the magnitudes that matter; the between-rule size gap (+2.1 / +0.6 / +1.6 subjects), the E|Ĥ| residual (+1.5 / +0.6 / +1.1, 3.6 / 1.4 / 2.7 SE), and the naive-bias discrepancy with its n = 500 / n = 700 sign flip (−1.2 / +1.5 / −1.4) are **unexplained** by `se_g`'s scaling and by the dedup. What this task did *not* test, and which remains on the argmax side by elimination: sample-quantile cuts (the search enumerates on the trial's own quantile grid, the analytic on the population's), the `rmin` translation, the pre-screen ordering, and the identification `σ_D = se_g` in the consistency gate. No recommendation.

---

## 5. Close-out

`git status --porcelain` before staging: the new files only. Staged by explicit path: `residual_2026-08-30.R`, `residual_hypA_2026-08-30.rds`, `residual_hypA_2026-08-30.log`, `residual_hypB_alt500_2026-08-30.log`, `residual_hypB_alt700_2026-08-30.log`, `residual_hypB_null500_2026-08-30.log`, `residual_hypB_summary_2026-08-30.log`, `residual_hypB_2026-08-30.rds`, this report. The three per-cell part files live in the session scratch directory and are merged into `residual_hypB_2026-08-30.rds`. **No push. No install** (nothing in `R/` changed).

Commits: `80544be1` task doc; `243801e6` scripts, outputs and this report; the hash line in the follow-up commit.

---

## 6. Summary (ten lines)

1. A is arithmetic: `se_scaled/se_direct = sqrt(V_eff[Q]/V_eff[g])`, n-free; on this fixture σ² = 16 256 caps the candidate dependence of `V_eff` at ≈ 4 %.
2. Measured: the ratio is in [0.992, 1.015] everywhere, median 1.000–1.004, selection-weighted mean within 0.1 % of 1 at the alternatives.
3. No prevalence trend on the realized rules (|cor| ≤ 0.04); a 0.2 % purity slope with the hypothesised sign and none of the magnitude; large-low-purity candidates at 1.0006–1.0009; controls at 1.000.
4. **A refuted.** (Side fact: the wrapper omits the Jensen count-randomness inflation, 1–2.5 %, trending the *wrong* way for the gap.)
5. The dedup key on the GLM path is `(K, n, n1, β̂, ŝe)` to 0.001 with `m1 = m0 = NA`; in a trial it collapses identical sample memberships within K and nothing else.
6. The literal population key reduces to `(K, |g|, |g∩Q|)`; it removes 12 / 19 / 458 candidates that are *unrelated* rules (Jaccard ≈ 0.2, share ≥ 0.99 = 0) — not a faithful analogue; its null result (−0.58, past the measurement) is an over-collapse artefact.
7. The near-twin analogue (symdiff ≤ 6 / 4 subjects, the `str2`/`preanti` twins) removes 10 / 7 / 10 candidates carrying 1.1–1.4 % of the mass — an upper bound on the search's collapse.
8. Re-runs (both gates, same seed): between-rule +2.11 → +2.06, +0.61 → +0.58, +1.65 → +1.59; sens/PPV move ≤ 0.001; naive bias ≤ 0.08, sign flip intact.
9. **B accounts for 2–5 % of the gap at most, ≈ 1–3 % expected.**
10. **Neither explains the residual; the residual is unexplained.** No `R/` change, one script, two `.rds`, one report; one DGM, two n, one outcome type.
