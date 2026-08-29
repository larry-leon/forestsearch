# REPORT — the last residual lead: population-quantile cuts against replicate-quantile cuts

**Task:** `dev/tasks/cc_task_residual_quantiles_2026-08-30.md` (commit `d9ccc24c`)
**Repo / branch:** `~/Documents/GitHub/forestsearch` (pop-os), `feature/glm-extension`
**Follows:** `REPORT_residual_2026-08-30.md` §4 (the lead), `REPORT_oc_wrapper_confs_2026-08-30.md` (the corrected family), `REPORT_oc_wrapper_sgdef_2026-08-29.md` (the between-rule finding), `REPORT_fs_family_report_2026-08-30.md` (15 of 30 thresholds differ).
**Category:** no `R/` change, no export, no edit to any package file, driver or document. Artefacts, all under `dev/glm-continuous-sims/`: `residual_quantiles_2026-08-30.R` (parts `stage1 <cell>`, `gate`, `stage2 <cell> <variant>`, `summary`), `residual_quantiles_stage1_{alt500,alt700,null500}_2026-08-30.{rds,log}`, `residual_quantiles_gate_2026-08-30.{rds,log}`, `residual_quantiles_stage2_{cell}_{located,full}_2026-08-30.{rds,log}` (six each), `residual_quantiles_summary_2026-08-30.log`, this report. No renders, no push, no install.

---

## 0. The answer in one paragraph

**The mechanism is real, measurable, and points the wrong way.** Building the search's cut matrix on each replicate instead of on `df_super` shifts the selection-weighted prevalence of the candidates the analytic argmax picks by **−0.43 / −0.82 / −0.38 subjects** (n = 500 alternative / n = 700 / null) — 20 % / 135 % / 23 % of the +2.11 / +0.61 / +1.65 between-rule gap, **with the opposite sign**. Re-running `fs_oc_predict()` on a replicate-averaged family confirms it in the functionals: `E|Ĥ|` falls by 0.9–1.4 subjects, the between-rule gap **widens** to +3.0 / +2.0 / +2.6, sensitivity, PPV and the naive bias all move away from measurement, and the inverted `c1` rises by 0.4–0.7. Replicate-quantile cuts do not explain the residual; they make the analytic family select *smaller* rules still. **The residual is unexplained, and by Larry's decision it is closed as a documented limitation.**

---

## 1. Provenance (§1 — GATE passed)

```
pop-os · /home/larryleon/Documents/GitHub/forestsearch · feature/glm-extension
HEAD at start d9ccc24c (after the task-doc commit; the tree before it was 5196ea35, clean)
git status --porcelain: empty
git log -4: d9ccc24c 5196ea35 774d768c 1e45219f
packageVersion forestsearch 0.3.0
```

`~/Downloads/cc_task_residual_quantiles_2026-08-30.md` was absent under that name; the stem arrives hyphen-stripped as **`cc_task_residual_quantiles_20260830.md`** → `dev/tasks/cc_task_residual_quantiles_2026-08-30.md`, committed alone as **`d9ccc24c`**.

**Installed version: 0.3.0.** The dispatch task (`cc_task_parallel_dispatch_2026-08-30.md`) ran first in this session and **stopped at its §6 gate** — its edit is preserved as a patch, not installed — so 0.3.0 is the version in force. Neither 0.3.0 nor the proposed 0.3.1 touches anything used here (`get_FSdata()`, `fs_oc_family_enumerate()`, `fs_oc_predict()`, `simulate_from_glm_dgm()`); the dispatch change is to the consistency screen's execution path only.

DGMs rebuilt and gated against the payloads' `truth` as in every prior task; the corrected families were re-enumerated in every part and asserted `identical()`/`all.equal()` to the stored record (`lab`, `Pg`, `M` = 1696 / 1890 / 1696). The population cut list rebuilt by the script's `cut_matrix()` on `df_super` was asserted `identical()` to the stored family's `cuts` (37 expressions).

Compute: stage 1 ≈ 12 min per cell (three concurrent); stage 2 six `fs_oc_predict()` pairs + six inversions, ≈ 45 min (six concurrent processes); ~1 h 15 min in all, inside the estimate. Stage 1 was run twice — the second run adds per-replicate storage for stage 2 and reproduces the first run's every printed digit (same seeds, deterministic).

---

## 2. Stage 1 — does the cut placement move membership at all? (§2 — the gate did **not** close)

### 2.1 Design

- **Replicates:** R = **200** per cell, `simulate_from_glm_dgm(dgm, n, seed = 20260825 + r)`, r = 1…200 — the corrected run's `20260825` seed family (the measured payloads themselves used `8316951 + sim_id`; the task asked for the 20260825 family and that is what was used).
- **The search's cut matrix per replicate:** `get_FSdata()` with the drivers' arguments (13 confounders, `conf.cont_jcuts = list(age = 10, preanti = 10)`, `collapse_cuts = TRUE`, defaults otherwise), in **exactly** the call `fs_oc_family_enumerate()` makes on `df_super` (`R/fs_oc_family.R` L296–316), with the replicate frame in place of `df_super`; the two-direction indicator matrix via `dummy()` and the same `.fs_oc_column_labels()`. Cut construction depends only on the covariates and `n` (LASSO off), so this is what `forestsearch()` sees on that replicate.
- **Locating a population candidate in a replicate — by clause specification, never by threshold.** Each clause is keyed as `(variable, operator, rank of its threshold within that variable-and-operator's cut list, count of that list)` for continuous clauses and `(variable, level)` for binaries; a negated clause `age > 25` carries the rank of `25` in the `age <=` list. A candidate is located when every clause's key exists in the replicate; its replicate prevalence is the AND of the mapped columns. A candidate one of whose clauses has no key is **"disappears"** for that replicate; replicate keys matching no population clause are **"appears"**.
- Per replicate and candidate: `ΔP = replicate prevalence − Pg`. Per-replicate `P(g∩Q)`, `P(g∩Qᶜ)` and a pairwise-valid running sum of the M×M overlap are kept for stage 2.

### 2.2 Label correspondence — *not* nearly perfect, and for two identifiable reasons

Population continuous cut counts (`df_super`, n = 5000): `age` 10, `preanti` 7, `wtkg` 3, `karnof` 2, `cd40` 4, `cd80` 4. Per replicate:

| variable | population count | replicate mode | range | share of replicates equal to population (alt500 / alt700) | why it differs |
|---|---:|---:|---|---|---|
| `age` | 10 | 10 | 10–10 | 1.000 / 1.000 | the J = 10 grid is specified, never collapses |
| `preanti` | 7 | 7 | 6–8 | 0.876 / 0.909 | the four `preanti <= 0` grid points collapse to one everywhere; the fifth/sixth (55.8, 214 in the driver's comment) sometimes merge or split at the `sd/√n` band |
| **`wtkg`** | **3** | **4** | 3–4 | **0.419 / 0.361** | the default four cuts (mean, median, Q1, Q3): on `df_super` the rounded mean and median coincide at 75 and the list is 3; in a sample the mean rounds to 74 against a median of 75 and the two stay apart (|Δ| = 1 > band 0.58), so the list is 4 |
| **`karnof`** | 2 | 2 | 2–2 | 1.000 / 1.000 as a count — **but in ≈ 40 % of replicates `karnof` has no `<=` cuts at all** (it is absent from the key set in 81 / 200 replicates at n = 500): the sample carries fewer than `cont.cutoff = 4` distinct values, `is.continuous()` returns FALSE, and `get_FSdata()` treats it as a **categorical** factor with level indicators | a type change, not a threshold change |
| `cd40` | 4 | 4 | 3–4 | 0.943 / 0.977 | occasional collapse of two default cuts |
| `cd80` | 4 | 4 | 4–4 | 1.000 / 1.000 | — |

Consequences: candidates located in **every** replicate: **509 / 1696, 563 / 1890, 509 / 1696**; mean located per replicate 1389 / 1579 / 1389 (82–84 %); located in none: 0. Candidates that disappear in at least one replicate: 1187 / 1327 / 1187 — those carrying `preanti` (594 / 668 / 594), `cd40` (408 / 450 / 408), `wtkg` (307 / 347 / 307), `karnof` (116 / 128 / 116); they carry **≈ 71 %** of the analytic selection mass. The `wtkg × karnof` candidates are the worst case, missing in 141 / 200 replicates (the union of the two effects). Fully corresponding replicates (all M located): **45 / 48 / 45** of 200.

"Appears" keys: the always-present ones are clause keys with no population *candidate* — `age <= 25`, `age > 47`, `preanti > 1055`, `hemo` (each below the 0.12 size floor on `df_super` as a single clause), `karnof <= 90` / `> 90` (the rank-1 `karnof` cut, likewise floor-excluded; it "appears" in the 119 replicates where `karnof` is continuous) — and the count-mismatch keys (`wtkg|…|4`, `preanti|…|6` or `|8`, `cd40|…|3`). Nothing appears that is not explained by the floor or by the collapse/type mechanics above.

**What this means for the estimate.** `ΔP` is computed only where a candidate is located, so stage 1's shift is *conditional on correspondence*; in the other replicates the search's family is structurally different (an extra `wtkg` cut; `karnof` as a factor) — a family-composition difference this task's design cannot express, and one that `REPORT_fs_family_report_2026-08-30.md`'s "15 of 30 thresholds differ" did not distinguish from threshold placement. Restricting to the 509 / 563 / 509 always-located candidates gives the same sign and a slightly larger magnitude (−0.53 / −1.28 / −0.46 subjects), so the conclusion does not hinge on the conditioning.

### 2.3 The shift

Per-candidate mean `ΔP` over replicates:

| cell | mean | median | sd | range | unweighted, in subjects (× n) | by K = 1 / K = 2 | by `Pg` quartile Q1 → Q4 | cor(E[ΔP], `Pg`) | lm slope on `Pg` (SE) |
|---|---:|---:|---:|---|---:|---|---|---:|---|
| alt500 | +0.00031 | −0.00012 | 0.0067 | [−0.023, +0.036] | +0.16 | −0.00009 / +0.00033 | −0.00117 · −0.00093 · +0.00113 · +0.00222 | +0.19 | +0.0077 (0.0010) |
| alt700 | +0.00017 | −0.00009 | 0.0062 | [−0.024, +0.035] | +0.12 | −0.00002 / +0.00017 | −0.00119 · −0.00093 · +0.00064 · +0.00215 | +0.21 | +0.0075 (0.0008) |
| null500 | +0.00031 | −0.00012 | 0.0067 | [−0.023, +0.036] | +0.16 | −0.00009 / +0.00033 | −0.00117 · −0.00093 · +0.00113 · +0.00222 | +0.19 | +0.0077 (0.0010) |

(null500 equals alt500 to every digit: same covariate frame, same seeds — the DGM's outcome model does not enter the cut matrix.) There **is** a systematic dependence on `Pg`: small candidates (the two lower quartiles, `Pg` < 0.30 — where the winners live) are ≈ 0.1 % of the sample *smaller* under replicate cuts, large ones ≈ 0.2 % larger. That is the expected geometry: a candidate is the intersection of two clauses each fixed at a sample quantile, and near the size floor the intersection of two sample-quantile sets is on average slightly smaller than the population-quantile intersection evaluated on the same sample.

**The one that matters — `sel_c`-weighted, in subjects, against the gap (the gate's arithmetic):**

| cell | n | `sel_c`-weighted E[ΔP], resample / split | × n = **subjects** | between-rule gap | **share of the gap** |
|---|---:|---|---:|---:|---:|
| alt500 | 500 | −0.00086 / −0.00084 | **−0.430 / −0.418** | +2.11 | **−20.4 % / −19.8 %** |
| alt700 | 700 | −0.00118 / −0.00119 | **−0.824 / −0.830** | +0.61 | **−135 % / −136 %** |
| null500 | 500 | −0.00077 / −0.00074 | **−0.384 / −0.370** | +1.65 | **−23.3 % / −22.4 %** |

**Gate.** The rule was: skip §3 if the shift is under 20 % of the gap *and* shows no systematic dependence on K or `Pg`. It fails on both clauses in every cell — the shift is at or above 20 % (135 % at n = 700), it has a `Pg` gradient, and it is **wrong-signed**. §3 was entered, as the task directed for "a wrong-signed effect of the right magnitude".

---

## 3. Stage 2 — the replicate-averaged family, and what it predicts (§3)

### 3.0 The approximation, stated before the numbers

The analytic framework requires a **fixed** family: one candidate set, one `Pg`, one `β_g`, one `se_g`, one overlap matrix, and the argmax distribution follows from a single multivariate normal. What the search actually faces is a family that varies with the replicate — different thresholds, sometimes a different cut count, sometimes a different variable type. That object is outside the framework. Averaging the family's fields over replicates asks a narrower question: **does the *mean* sample-quantile family close the gap?** It is a first-order test of the mechanism — sufficient to see the sign and rough size of what the mechanism does to the argmax — and **not** a proof of what a replicate-varying family would do, **not** a constructor, and **not** proposed for adoption.

Where the averaging is least defensible, and what it could do: (i) **the overlap matrix** — the mean of replicate Gram matrices is not the Gram matrix of any family, and dividing it by `sqrt(P̄_g P̄_h)` does not give a correlation matrix of anything realisable; a covariance built from it can be non-PSD (it *is*, badly, in one variant below), and `fs_sym_root()` then silently clamps negative eigenvalues, distorting the joint law of the draws; (ii) **`se_g`** — kept at the corrected family's formula `sqrt(V_eff[Q] / (n·P̄_g))` at the replicate-mean prevalence, which the residual task showed is right to 2 % for any fixed membership but says nothing about a membership that itself varies; (iii) **`β_g`** — `τ_Qc + b·P̄Q_g` with `P̄Q_g` the replicate mean of `P(g∩Q)/P(g)`, the purity of the *average* candidate rather than the average purity of the realised ones (Jensen-level differences, second order). Either of (i) or (ii) could move `E|Ĥ|` by a few tenths of a subject in either direction; neither could plausibly flip a 1–1.4-subject movement's sign, which is the finding.

### 3.1 Two variants, one primary

| variant | fields | replicates per field | min eigenvalue of ρ | status |
|---|---|---|---|---|
| **`full`** (primary) | every field the mean over the **fully corresponding** replicates only (all M candidates located) | **45 / 48 / 45** | −0.0000 (PSD by construction: an average of Gram matrices under a diagonal congruence) | consistent; fewer replicates |
| `located` (sensitivity) | each candidate's `Pg`, `P(g∩Q)` over the replicates in which it was located (59–200); `ovl[i, j]` over the replicates in which *both* were located (pairwise) | 59–200 | **−0.47 / −0.37 / −0.47** | more data, but the pairwise-averaged overlap is far from PSD; `fs_sym_root()` clamps, so the joint law is distorted — reported, not relied on |

`M` unchanged (1696 / 1890 / 1696); every other field as in §3.0; `fs_oc_predict()` both gates and `fs_oc_invert()` (resample gate, the corrected run's targets) at seed 20260825, 2×10⁵ draws.

### 3.2 Results — corrected (stored) → replicate-averaged, against measurement

**`full` variant (primary).** "between" = population size of the realized rules (72.90 / 73.65 / 72.55) − analytic `E|Ĥ|`.

| cell | gate | E\|Ĥ\| | **between-rule** | E\|Ĥ\| gap to measured | E[sens] (meas) | E[PPV] (meas) | naive bias (meas) | declaration |
|---|---|---:|---:|---:|---|---|---|---:|
| alt500 | resample | 70.79 → **69.87** | +2.11 → **+3.03** | +1.55 → +2.46 | 0.1638 → 0.1603 (0.1706) | 0.4000 → 0.3968 (0.4055) | 75.53 → 76.09 (74.27) | 0.99906 → 0.99904 |
| alt500 | split | 70.78 → 69.91 | +2.12 → +2.99 | +1.55 → +2.43 | 0.1635 → 0.1602 | 0.3994 → 0.3964 | 75.38 → 75.92 | 1.00000 → 1.00000 |
| alt700 | resample | 73.04 → **71.66** | +0.61 → **+1.99** | +0.58 → +1.97 | 0.1219 → 0.1190 (0.1229) | 0.4034 → 0.4024 (0.4016) | 75.57 → 76.31 (77.09) | 0.99976 → 0.99979 |
| alt700 | split | 73.08 → 71.68 | +0.57 → +1.97 | +0.55 → +1.94 | 0.1222 → 0.1192 | 0.4038 → 0.4028 | 75.49 → 76.22 | 1.00000 → 1.00000 |
| null500 | resample | 70.90 → **70.00** | +1.65 → **+2.55** | +1.15 → +2.05 | — | — | 76.32 → 76.88 (74.85) | 0.99694 → 0.99706 |
| null500 | split | 70.83 → 69.95 | +1.72 → +2.60 | +1.22 → +2.09 | — | — | 76.07 → 76.62 | 0.99999 → 0.99999 |

Inverted `c1` (resample gate), corrected → replicate-averaged: alt500 targets 0.80 / 0.90 / 0.95: 92.04 → 92.48, 84.95 → 85.31, 79.14 → 79.51; alt700: 93.19 → 93.78, 86.64 → 87.18, 81.28 → 81.82; null 0.05 … 0.95: 133.24 → 134.00, 125.70 → 126.49, 117.11 → 117.75, 101.64 → 102.19, 87.29 → 87.75, 80.17 → 80.58, 74.43 → 74.84. Every threshold moves **up** by 0.4–0.6 — the direction of a noisier, smaller-winner family — and every declaration rate is within 0.0002 of the stored one.

**`located` variant (sensitivity, non-PSD overlap clamped):** E|Ĥ| 70.79 → 69.99 / 73.04 → 71.77 / 70.90 → 70.07 (resample); between-rule +2.91 / +1.88 / +2.48; sensitivity, PPV, naive bias and the inversions move the same way by 85–95 % of the `full` amounts (full tables in `residual_quantiles_summary_2026-08-30.log`). The two variants agree in sign on every quantity and within 0.15 subject on E|Ĥ|; the distortion of the clamped covariance is not what drives the result.

### 3.3 Reading

- **Direction: away from measurement, on every quantity, in every cell, under both gates and both variants.** E|Ĥ| falls 0.9 / 1.4 / 0.9 subjects; sensitivity and PPV fall (the measured values sit above the corrected prediction at n = 500, so the gap opens); the naive bias rises at n = 500 and the null (measured is *below* the prediction — gap opens) and rises at n = 700 (measured is above — gap narrows by 0.7, the one quantity on which the mechanism helps, and only because it is the cell whose naive-bias sign is the anomaly).
- **Magnitude: the between-rule gap widens from +2.11 / +0.61 / +1.65 to +3.03 / +1.99 / +2.55**, i.e. the mechanism, if it acts as its mean, accounts for **−43 % / −226 % / −55 %** of the gap. Stage 1's −20 % / −135 % / −23 % was the first-order membership shift alone; the argmax amplifies it because the shift is concentrated on the small candidates the argmax favours (§2.3's `Pg` gradient).
- The mean shift in `Pg` over the whole family is +0.0003 / +0.0002 / +0.0003 — the *family* barely moves; the *winners* do.

---

## 4. Verdict (§4)

1. **Stage 1's `sel_c`-weighted shift, in subjects, with its sign:** **−0.43** (alt500, gap +2.11, −20 %), **−0.82** (alt700, gap +0.61, −135 %), **−0.38** (null500, gap +1.65, −23 %). Negative in every cell: replicate-quantile cuts make the candidates the analytic argmax selects *smaller* in their own sample than the population cuts do, with a clear `Pg` gradient (slope +0.0077 ± 0.0010 per unit prevalence) and no K dependence.
2. **§3 ran. The replicate-averaged family accounts for none of the gap and widens it:** +2.11 → +3.03, +0.61 → +1.99, +1.65 → +2.55 subjects (primary variant; the sensitivity variant within 0.15). E|Ĥ|, sensitivity, PPV and the inverted `c1` all move away from measurement; the naive bias moves away at n = 500 and the null and toward at n = 700 by 0.7 of a 1.5 discrepancy.
3. **This does not explain the residual. The residual is unexplained.** Five mechanisms have now been measured — sample realization (within-rule, sgdef report), the confounder list, `se_g`'s prevalence scaling, the dedup asymmetry (≤ 5 %), and now cut placement (wrong sign, would widen the gap) — and the between-rule size gap of +2.1 / +0.6 / +1.6 subjects, the E|Ĥ| residual of +1.5 / +0.6 / +1.1, and the naive-bias discrepancy with its n = 500 / n = 700 sign flip stand as before. By Larry's decision of 2026-08-30 the question is **closed as a documented limitation** of the analytic prediction on this fixture: no further mechanism is proposed, and none is implied by anything measured here. One fact from §2.2 is recorded for the limitation's wording, not as a lead: on a fifth to two-fifths of replicates the search's family is not the population family relabelled — `wtkg` has four cuts instead of three, `karnof` is a categorical factor instead of a continuous cut — and a fixed family cannot represent that.

No recommendation.

---

## 5. Close-out

`git status --porcelain` before staging: the new files only. Staged by explicit path: `residual_quantiles_2026-08-30.R`; `residual_quantiles_stage1_{alt500,alt700,null500}_2026-08-30.rds` and `.log`; `residual_quantiles_gate_2026-08-30.rds` and `.log`; `residual_quantiles_stage2_{alt500,alt700,null500}_{located,full}_2026-08-30.rds` and `.log`; `residual_quantiles_summary_2026-08-30.log`; this report. **No push. No install** (nothing in `R/` changed; 0.3.0 stays installed).

Commits: `d9ccc24c` task doc; scripts, outputs and this report in the next commit; its hash in the follow-up line.

---

## 6. Summary (ten lines)

1. R = 200 replicates per cell, seeds 20260825 + r; the search's cut matrix rebuilt per replicate with `get_FSdata()` exactly as the population family builds it; candidates matched by clause specification (variable, operator, quantile rank), never by threshold.
2. Label correspondence is 82–84 % per replicate and full in only 45–48 of 200: `wtkg` keeps four default cuts in a sample where `df_super` collapses to three (rounded mean = median), and `karnof` becomes a **categorical** variable in ≈ 40 % of replicates (< 4 distinct values); the affected candidates carry 71 % of the selection mass.
3. Per-candidate prevalence shift: family mean +0.0003, but with a `Pg` gradient — the small candidates the argmax favours are ≈ 0.1 % smaller under replicate cuts, the large ones ≈ 0.2 % larger.
4. **`sel_c`-weighted shift: −0.43 / −0.82 / −0.38 subjects against gaps of +2.11 / +0.61 / +1.65 — 20 % / 135 % / 23 %, wrong sign.** Gate not closed; §3 run.
5. Replicate-averaged family (primary: the 45–48 fully corresponding replicates, PSD; sensitivity: pairwise means, non-PSD and clamped): both agree.
6. **E|Ĥ| 70.79 → 69.87, 73.04 → 71.66, 70.90 → 70.00; between-rule gap +2.11 → +3.03, +0.61 → +1.99, +1.65 → +2.55** — the mechanism widens the gap by 43–226 %.
7. Sensitivity, PPV, the naive bias (except at n = 700) and every inverted `c1` (+0.4–0.7) move away from measurement; declaration rates unchanged to 0.0002.
8. The approximation (averaging a replicate-varying family) is stated before the numbers; it can move E|Ĥ| by tenths, not flip a sign; not a constructor, not proposed.
9. **Cut placement does not explain the residual; it points the other way. The residual is unexplained and, by decision, closed as a documented limitation.** No lead 6.
10. No `R/` change; one script, nine `.rds`, eleven logs, one report; one DGM, two n, one outcome type; 0.3.0 installed throughout.
