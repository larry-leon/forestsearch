# REPORT — ACTG175 binary (OR) campaigns `orfs`, `orgrf`, `ordina`: Stage 1 / Gate 1

Date: 2026-09-17 (UTC). Machine: `pop-os` (64 physical cores, 251 GB; R 4.6.1). Branch `feature/glm-extension`.
Task: `dev/tasks/TASK_actg175_binary_campaign_2026-09-17.md` (committed first and alone, `60ca1b8c`).
Directory: `quarto/simulations/actg175/binary_020/`, written `<dir>`.

**Outcome: Gate 1 is GREEN on every gate, and the advance go DOES NOT APPLY.** The §1.7 projection
for Stage 2 at the chosen worker count is **19.48 h**, above the dispositions' 12 h threshold, so
this record stops at Gate 1 as the dispositions require. Nothing was reduced to fit under the
threshold — no cell, replicate count, gate or knob was changed.

---

## 1. Provenance — §1.1 GATE PASS

```
pop-os
feature/glm-extension
98475a05 IS in HEAD
[tracked modifications: 0]
[R / Rscript / quarto processes: none]
```

First commit, alone: `60ca1b8c dev/tasks: TASK_actg175_binary_campaign_2026-09-17.md as received`.
`<§1.1 HEAD>` = `60ca1b8c`.

**Finding (not a gate failure).** One orphaned `bash` process (PID 2086198, parent 1, ~11 d 22 h
elapsed) survives from an earlier session in a different working directory: it is the memory-sampler
loop of a `gbsg_020` run whose renders are long gone, still writing a running maximum to a dead
scratchpad path every 3 s. No `R`, `Rscript` or `quarto` process exists, so the §1.1 assertion holds
as written. It was left alone (out of this task's scope) and is recorded here so it is not
rediscovered as a surprise; it also matches `pgrep -f quarto` on its stale argument string, which is
why a naive process check reports a false positive.

## 2. Package check — §1.2 GATE PASS

```
git diff --quiet 064fce91..HEAD -- R/ DESCRIPTION NAMESPACE   ->  unchanged
packageDescription("forestsearch")$Built  ->  R 4.6.1; ; 2026-09-17 04:47:31 UTC; unix
two doFuture workers (multisession)       ->  both report the identical string
```

No `R/` change and no install anywhere in Stage 1. `forestsearch` 0.3.5.

## 3. The template — §1.3

`<dir>/sim_fs_mr_field_or_template.qmd`, committed `0cc4c015` (with `7484e62b` and `d5835deb`
below). Created by copying the MD field template and changing only the named blocks; the document's
own §sec-ledger quotes **both** sources for every block, and is the authoritative ledger. Sources:

- **[MD]** `../continuous/sim_fs_maxeffCons_mr_field_md_template.qmd` — the structure copied.
- **[S]** `<dir>/maxeffCons_mr_coverage_sweep_or075.qmd` (`1d42f6da`) — the study driver, this
  task's source for the data recipe, seeds, truths and detector settings.
- **[S10]/[S15]** `../binary/mr_coverage_sweep_or10.qmd` / `mr_coverage_sweep_or15.qmd` — consulted
  for `target_or_h` and calibration settings only.

| block | from | change |
|---|---|---|
| knobs | [MD] :93–374 | renamed `FS_OR_*`, plus `FS_OR_TARGET` (0.75 \| 1.0 \| 1.5) and `FS_OR_N`; `FS_MD_KNOISE`, `FS_MD_BUILD`, `FS_MD_QUICKRUN`, `FS_MD_FB*`, `FS_MD_UNIFORM`, `FS_MD_DMIN_GRF`, `FS_MD_JOIN_SKIP` dropped (pinned to the study's / the campaign's literal) |
| data prep, DGM, truth, eval frame | [S] :233–341 | verbatim except `target_effect = FS_OR_TARGET` ([S] :245; reachable range `k_inter_range = c(0.3, 1.5)` is [S] :222, unchanged in [S10] :222 / [S15] :222) |
| seeds | [S] :119–133, :447–456 | verbatim: the pre-generated table indexed by global `sim_id`, and the per-replicate `RNGkind("L'Ecuyer-CMRG")` / `seedit` |
| per-replicate data prep | [S] :457–484 | verbatim (id, factor→numeric before every fit, confounders, the noise offset) |
| oracle | [S] :360–380 | `.logit_or_ci()` verbatim, recorded under [MD]'s `or_*` names, not [S]'s `ora_*`, so `fs_sim_bias_coverage(target = "oracle")` reads it unrenamed |
| thresholds and search | [S] :175–207 | `sg_focus`, `effect_neighborhood`, `selection_rule` from the knobs (the dispositions: `effMaxSG`, 0.20, `neighborhood`); everything else the study's |
| GRF and DINA arguments | [S] :205–207 | verbatim (`dmin.grf = 0.0`, frontier / effect / depth 2; DINA `effect`, `dina_args = list()`) |
| MR | [MD] :343–348 | the field list replaces [S]'s `ci_method = "ij"` ([S] :516); `consistency_method = "resample"` is [S] :179, which MR requires |
| recorder | [MD] :516–636 | naming unchanged (`nv_*`, `mr_*`, `or_*`, `fld_*`), **plus** `seed` ([S] :512, a data-level Gate 2 column) and `C_dagger_H/Hc`, `C_ddagger_H/Hc` from the study's truth table, per replicate, beside `betaHhat_*` |
| `betaHhat` attachment | [S] :703–707 | verbatim |
| `meta` and poolability keys | [MD] :1135–1175, :1210–1225 | plus `target_or_h`, `design_tag`, the five truths, `sg_quantile`, `n_super`, `eval_seed`, `seed_scheme`, `outcome_type`, `effect_measure` and the thresholds ([S] :716–748) |
| stem and results dir | [MD] :190–216 | new root `mr_or_harm/`, design token `or075`/`or100`/`or150` and `n` in the stem, mandatory campaign tag; `.refuse_if_tracked()` live at both save sites, so no committed study path can be written |

**Three named drops**, each with its reason, all in the document's ledger:

1. `fs_dgm_scale()` and the SCALE CHECK — the function refuses ratio measures (`R/fs_dgm_scale.R:130`).
   `dgm_scale <- NULL`, so [MD]'s own `is.null()` guards skip the check and leave `scale` out of the
   payload. `fs_mr_oc_summary()` is **kept**: the truth list carries `effect_Q` / `effect_Qc` /
   `prevalence_Q` / `beta_inter` aliases beside the study's own names, which is what that function
   reads (`R/fs_mr_oc_summary.R:131–190`).
2. The `fs_identification_summary()` anchor / partner / proxy triple — that triple is the MD design's
   age / preanti / str2; the binary true region {wtkg > q70} ∩ {cd40 > q70} has no binary proxy and
   `proxy` has no default (`R/fs_identification_summary.R:222`). The classification table, the
   covariate-frequency figure and the realized-rule table are kept, re-pointed at wtkg / cd40.
3. The FB join — there is no committed binary FB bundle. `fb_mode` is pinned `"none"`; the `fb_*`
   recorder columns are kept all-NA for schema parity, and [MD]'s own `fb_in_tables` test drops the
   FB rows automatically.

**Two changes made after the first commit, both forced by the smoke and both named in the ledger:**

- `7484e62b` — **the band-width guard becomes a NOTE.** [MD] makes a non-band `sg_focus` with
  `FS_MD_NBHD` set a hard error, sound there because its default ε (0.10) was already
  `forestsearch()`'s own, so leaving it unset reproduced every committed non-band render. This
  template's default is the campaign's 0.20, so the same rule would make the **study's** rule
  (`maxeffCons` at ε 0.10) unreachable — and §1.5(a) must render exactly that. The study driver
  itself sets `effect_neighborhood = 0.10` alongside `sg_focus = "maxeffCons"` ([S] :151–154), so a
  NOTE that ε is forwarded but inert replaces the stop. ε is still recorded in `meta` and still tags
  the stem.
- `d5835deb` — **working vs effect scale on a ratio measure** (§5 below).

## 4. Scripts — §1.4

`<dir>/scripts_or/`, committed `f25162f7` (plus `d5835deb`). Copied from
`../continuous/scripts_mddina/` with the OR fields and `meta` items:

- `mem_sampler.sh` — unchanged.
- `smoke_identity.R` — four modes for §1.5's four renders: `recipe`, `fs`, `grf`, `dina`. §1.5(d)'s
  construction checks run in **every** mode.
- `gate2.R` — the per-cell Stage 2 gate, taking the campaign tag as its 4th argument: 2,000 rows with
  `sim_id` exactly 1–2000, batch/combined column identity, the `meta` set, the §1.5(d) checks, MR
  failures ≤ 40, and — for `orgrf` and `ordina` — the same-draws check against `orfs`'s bundle in the
  same cell **in both directions** (the per-replicate seed, `n_true`, the eight oracle columns and
  the truths, within 1e-8 relative).
- `run_or.sh` — the three campaigns in the task's order (`orfs`, `orgrf`, `ordina`), each over the six
  cells in the task's order (0.75/500, 0.75/2000, 1.5/500, 1.5/2000, 1.0/500, 1.0/2000); two batch
  renders plus a combine per cell, memory sampled per render, per-cell explicit-path commits,
  `LOG_or_progress.txt` / `HALT_or.md` / `logs_or/` and the Gate 2 record. **No hard-coded model or
  co-author trailer on any commit message** (the task's carried fix).

## 5. Smoke — §1.5 GATE PASS (all four renders)

Four renders of 20 replicates (`sim_id` 1–20) at `target_or_h = 0.75`, n = 500, tag `orsmoke`,
20 workers. Checker: `scripts_or/smoke_identity.R`.

### 5(a) The data recipe — GATE PASS

With the knobs at the study's rule, against the committed `fs_mr_n500_res.rds` for `sim_id` 1–20:

```
truth targets match the committed study within 1e-08 relative (max 1.41e-15)
  or_causal 0.6665537395 | marg_H 0.7499999955 | marg_Hc 0.6560116720 | cde_H 0.7321189340 | cde_Hc 0.6313905111
  prevalence(H) 0.096320 | k_inter 0.1480184
seed identical on every row (sim_id 1 seed 1530735852, committed 1530735852)
n_true identical on every row (sim_id 1: 47 vs 47)
oracle columns (or_* here vs ora_* committed, 8 columns) within 1e-08 relative
  on all 15 rows the committed bundle scores (max 6.46e-15)
```

**One named difference, reported not gated.** The committed study computes its oracle *after* its
detection return ([S] :628–631), so its `ora_*` columns are NA on its five NO-DETECTION rows
(`sim_id` 2, 3, 6, 7, 16). This template follows [MD] and computes the oracle *before* the search, so
it is filled on all 20 replicates. The check therefore scores the 15 rows the committed bundle has;
the extra rows are strictly more information on a rule-independent quantity, not a disagreement about
a value.

**Reported beside, not gated** (rule-dependent, and the MR construction differs — field vs IJ): the
same realized rule on **20 of 20** rows, status agreeing on 20, and the naive OR within 3.2e-15. So on
this path the field block is add-only in the draw stream as well.

### 5(b)–(c) The three identifiers under the campaign rule — GATE PASS

Every meta assertion passed for all three: the rule (`effMaxSG`, ε 0.20, `neighborhood`), the
thresholds (0.90 / 0.80 / 0.90 with `adverse_outcome = TRUE`), `target_or_h` 0.75 and `design_tag`
`or075` on the calibrated `alt` branch, the five truths plus prevalence, the MR construction set
(field, complement, `field_scale_complement = "selected"`, `two_term`, reselection, 5,000 draws),
`pkg_version` 0.3.5 and host `pop-os`, and each identifier's own arguments (GRF `dmin` 0, frontier,
depth 2, effect; DINA effect, `list()`).

**Zero factor-comparison warnings and zero NA-membership warnings on every one of the four renders —
in fact `warn_msg` is NA on all 80 replicates.** No CONFIG-ERROR row anywhere.

**DINA's floors, as applied.** `forestsearch()` derives the proposal floor as
`m_diff = log(hr.threshold)` for every non-Gaussian family (`R/forestsearch_helpers.R:1434–1437`) and
`.dina_collect_candidates()` drops any candidate with `mean_tau < m_diff`
(`R/dina_subgroup.R:748–749`). So the floor as applied **is** the OR-scale effect threshold 0.90 on
the harm side, expressed on the link scale: the smallest proposed τ̂ over the declared replicates is
−0.10535327 = log(0.90) + 7.2e-06, i.e. **OR 0.900007**, and
1 ≤ `admitted_n` ≤ `dina_proposed_n` on every declared replicate (searched 7,920–8,442; proposed
median 344; admitted median 187).

### 5(d) The constructions — GATE PASS on every declared replicate of (b) and (c)

Nine `fld_Hc_*_s` and nine `fld_joint_s_*` columns present and finite on every declared replicate
whose complement field block was filled (15 / 15, 20 / 20, 19 / 19, and 15 / 15 in (a)); zero
degenerate notes. `fld_Hc_lo1s_s ≤ fld_Hc_up1s_s` and `fld_Hc_lo2s_s ≤ fld_Hc_hi2s_s` throughout.
The Bonferroni harm bound is **identical** between `joint` and `joint_s` wherever the draw counts
agree (max |diff| 0.00e+00, all rows). Every bound is a positive OR (42 columns checked per render).
p̂(Ĥ) recorded on every declared replicate with a field block (mean 0.12–0.19; share below 0.5
0.89–1.00 — a pronounced tie regime on this design). γ and γ_s in [0.025, 0.026] throughout.

### 5 facts

| render | declared | `fit_mr_secs` mean / med / max | `fld_H_secs` mean / med | `fld_Hc_secs` mean / med | `sim_id` 1 selection (n_sel) | oriented log-OR | MR family | `fld_H_lo1s` | `fld_Hc_up1s_s` |
|---|---|---|---|---|---|---|---|---|---|
| (a) FS, study rule (`maxeffCons`, ε 0.10) | 15/20 | 27.0 / 33.2 / 37.7 | 11.7 / 11.6 | 2.28 / 2.28 | `{symptom} & {wtkg <= 79}` (61) | 0.8661 | 2238 | 0.2812 | 1.1186 |
| (b) FS, campaign rule | 15/20 | 33.0 / 40.6 / 48.0 | 18.8 / 18.7 | 2.70 / 2.69 | `{race} & !{karnof <= 95}` (73) | 0.8303 | 2238 | 0.2937 | 1.1332 |
| (c) GRF, campaign rule | 20/20 | 28.8 / 28.7 / 31.0 | 16.2 / 16.1 | 1.22 / 1.20 | `{age <= 28} & {preanti <= 777.2}` (89) | 0.6484 | 1051 | 0.2920 | 1.1375 |
| (c) DINA, campaign rule | 19/20 | 21.6 / 16.1 / 73.3 | 14.7 / 12.0 | 0.75 / 0.50 | `{age <= 28} & {preanti <= 777.2}` (89) | 0.6484 | 1157 | 0.3028 | 1.1313 |

**Beside S0 §7's fits.** Rows (b), (c) and (c) reproduce S0 §7's six fits at `sim_id` 1 **exactly** —
the same rule, the same size, the same oriented log-OR, the same MR family, and the same
`field$lower_1s` and `complement$upper_1s_s` to four decimals (S0 §7: FS `{race} & !{karnof <= 95}`
(73), 0.8303, 2238, 0.2937, 1.1332; GRF `{age <= 28} & {preanti <= 777.2}` (89), 0.6484, 1051, 0.2920,
1.1375; DINA the same selection, 1157, 0.3028, 1.1313). Row (a), under the study's own rule,
reproduces the committed payload's selection instead (`{symptom} & {wtkg <= 79}`, 61, naive OR
2.377622). GRF and DINA select the same region here; they diverge in their families (1,051 vs 1,157)
and in their bounds.

### 5 — three faults found by the smoke, all in transplanted identity-scale logic (`d5835deb`)

The data recipe was faithful on the first render; the failures were all in the checkers and in the
template's ratio lines, and all had one cause. `R/fs_mr_inference.R:480–488` states it: `lambda_mean`,
`lambda_sd`/`se_field` and the seven Λ* quantiles are on the **working** scale (the log-OR here),
while `est2` and every bound are on the **effect** scale (the OR).

1. **The identities.** [MD]'s additive forms are identity-scale specializations. Rewritten on the log
   scale they hold to 2e-16: `log(est2) + lambda_mean = log(β̃)` for both the field and field-s;
   `log(lo1s) = log(β̃) − q95`; `log(up1s) = log(β̃ᶜ) − q05`. The template's Λ-SD ratios now divide by
   SD(log est), not SD(est), and the Table-2 footnote records that every stored SE — Wald, IJ and the
   field's Λ-SD — shares the log-OR footing.
2. **Positivity.** `lambda_mean` is a log-scale correction and is routinely negative (13 of 15 rows), so
   it is not a bound. The checkers' positivity lists now cover the field-s **bounds** only. The
   template's own positivity guard was already correct.
3. **DINA's floor** — the log-scale reading recorded in §5(b)–(c) above.

## 6. Calibration — §1.6

The most expensive configuration, FS at `target_or_h = 1.5`, n = 2000, at `FS_OR_WORKERS` = 16, 32
and 63, with 3 × workers replicates, memory sampled every 5 s by `scripts_or/mem_sampler.sh`.
Threads pinned to 1 (`OMP`/`OPENBLAS`/`MKL`/`VECLIB`).

| W | replicates | render wall | reps/min | reps/min **per worker** | `fit_mr_secs` mean / median / p90 | id / field / complement (mean s) | peak summed RSS | k = wall / (waves × mean) |
|---|---|---|---|---|---|---|---|---|
| 16 | 48 | 487 s | 5.91 | 0.3696 | 138.8 / 140.7 / 148.5 | 10.5 / 34.2 / 10.0 | 28.0 GB | 1.169 |
| 32 | 96 | 599 s | 9.62 | 0.3005 | 158.9 / 159.9 / 177.2 | 13.3 / 39.1 / 13.1 | 55.4 GB | 1.257 |
| 63 | 189 | 1110 s | 10.22 | 0.1622 | 296.1 / 316.9 / 343.8 | 25.3 / 68.9 / 26.0 | 108.8 GB | 1.250 |

**A discarded and repeated measurement, recorded.** The first W = 32 render (tag `orcal32`) overlapped
two of my own single-core summary dry-runs inside its window, so it was discarded and re-measured
alone (tag `orcal32b`). The two agree: wall 598 s vs **599 s**, `fit_mr_secs` mean 159.0 vs **158.9 s**,
peak 56.1 vs 55.4 GB. The table above uses the clean repeat. W = 16 and W = 63 were never overlapped
(timestamps checked).

**A fourth render, to settle the projection.** The n = 500 side had no measurement at 63 workers, and
the projection's verdict turned on whether the heavy n = 2000 contention also applies there. One more
calibration render — FS, `target_or_h = 1.5`, **n = 500**, W = 63, 189 replicates (tag `orcal63n500`) —
settled it: **325 s wall, `fit_mr_secs` mean 46.4 s** (median 48.9, p90 59.5), peak 68.8 GB, 172 of 189
declared. Against S0 §7's one-worker FS n500 wall of 39.6 s that is an inflation of **1.172**, versus
**3.272** at n = 2000 (296.1 / 90.5). Cumulative Stage 1 compute ≈ 1.1 h, inside the kickoff's 3 h
ceiling.

### W = 63, and what limits it

**Chosen: W = 63.** It maximizes throughput — 10.22 reps/min against 9.62 at W = 32 and 5.91 at
W = 16 — and gives the shortest projection at every configuration (§7). Peak summed RSS is 108.8 GB
of 251 GB (43%), leaving headroom.

**What limits it is shared memory bandwidth / last-level-cache contention, not cores and not RAM
capacity.** The evidence:

- threads are pinned to 1, so each replicate's own work is single-threaded and unchanged by W;
- 63 of 64 physical cores are used, so the cores are not oversubscribed;
- peak RSS is 43% of RAM, so capacity is not the limit;
- yet the per-replicate wall inflates 138.8 → 158.9 → 296.1 s from 16 → 32 → 63 workers — **2.13×**
  from 16 to 63 — and the field block, which dominates the replicate, tracks it exactly
  (34.2 → 39.1 → 68.9 s);
- and the inflation grows with each replicate's working set: **1.172 at n = 500 against 3.272 at
  n = 2000**, relative to the same one-worker references. That n-dependence is the signature of
  shared-bandwidth and cache pressure, not of scheduling or of memory exhaustion.

Consequently throughput saturates between W = 32 and W = 63: doubling the workers buys **6.2%**
throughput and doubles the memory, and per-worker efficiency falls from 0.3005 to 0.1622 reps/min.

## 7. Projection, ceiling and timeout — §1.7

**Model.** Each configuration's per-replicate cost is S0 §7's one-worker wall scaled by the
**measured** inflation at the matching n and at W = 63 — 1.172 at n = 500, 3.272 at n = 2000 — both
measured on FS at `target_or_h = 1.5`, the most expensive design point, so the projection is
conservative for the other two. A render's wall is `waves × per-replicate + fixed`, with
`waves = ceil(1000/63) = 16` for a 1,000-replicate batch and `fixed` fitted on the two W = 63
calibrations (222 s from n = 2000, 186 s from n = 500; **220 s** used). A cell is two batch renders
plus one combine render (250 s, no replicates): `2 × (16 × per-rep + 220) + 250`.

| configuration | S0 §7 one-worker | inflation | projected per replicate | projected wall per cell | × 3 design points |
|---|---|---|---|---|---|
| FS n 500 | 39.6 s | 1.172 | 46.4 s (**measured**) | 2,175 s | 6,525 s |
| FS n 2000 | 90.5 s | 3.272 | 296.1 s (**measured**) | 10,165 s | 30,495 s |
| GRF n 500 | 26.6 s | 1.172 | 31.2 s | 1,687 s | 5,062 s |
| GRF n 2000 | 48.3 s | 3.272 | 158.0 s | 5,747 s | 17,241 s |
| DINA n 500 | 28.1 s | 1.172 | 32.9 s | 1,744 s | 5,232 s |
| DINA n 2000 | 11.1 s | 3.272 | 36.3 s | 1,852 s | 5,556 s |

- `orfs` 37,020 s = **10.28 h**; `orgrf` 22,303 s = **6.20 h**; `ordina` 10,787 s = **3.00 h**.
- **Stage 2 projection: 70,110 s = 19.48 h** at W = 63.
- **Ceiling (1.5 × the projection): 105,165 s = 29.21 h** cumulative render wall.
- **Longest projected batch render:** FS n 2000, 1,000 replicates — 4,958 s. **Per-render timeout
  (2 ×): 9,915 s = 2.75 h**, comfortably above the 20-minute floor.

**At the other two worker counts** (same model, the measured n = 2000 inflation at each W, and the
W = 63 n = 500 inflation as an upper bound at lower W): W = 32 → 22.33 h; W = 16 → 37.69 h. W = 63 is
the fastest of the three.

**Sensitivity.** The inflations were measured on `target_or_h = 1.5`, which declares on 99.5% of
replicates at n = 2000 and so runs the field block almost always. If the other two design points cost
20% less, the total is 16.88 h; 35% less, 14.93 h. **Both remain above 12 h**, so the verdict does not
turn on that assumption. The single dominant term is the three FS n = 2000 cells at 30,495 s = 8.47 h
between them.

### Advance-go decision — STOP AT GATE 1

Every Stage 1 gate is green, but the projection of **19.48 h** is **above** the dispositions' 12 h
threshold. The dispositions are explicit on both halves of this: *"if every Stage 1 gate is green and
the §1.7 projection for Stage 2 is under 12 hours, do not stop at Gate 1 … Otherwise stop at Gate 1
and report. Never reduce cells, replicates, gates or knobs to fit under the threshold."*

So Stages 2 and 3 were **not run**, and nothing was reduced: the six cells, 2,000 replicates per cell,
three identifiers, every gate and every knob stand exactly as the dispositions set them. Stage 2 is
ready to launch, unattended, with:

```
ORSG_WORKERS=63 ORSG_TIMEOUT=9915 ORSG_CEILING=105165 \
  setsid nohup bash quarto/simulations/actg175/binary_020/scripts_or/run_or.sh \
  > quarto/simulations/actg175/binary_020/logs_or/runner.log 2>&1 &
```

## 8. Stage 3 material prepared during Stage 1, committed but NOT run

Written while the calibration renders occupied the machine, and committed so the work is not lost.
**Status: dry-run verified only.** Neither has been run against production bundles, because Stage 2
did not run.

- `<dir>/summary_actg175_or.qmd` (§3.1) — renders end to end against the 20-replicate smoke bundles
  (`ORSG_SUMMARY_TAGS=orsmoke,orsmoke,orsmoke ORSG_SUMMARY_GLOB=res_1_20`) and writes all four
  artifacts. The dry run found and fixed two real faults: `c("name", numeric_vector)` coerced every
  identification and timing extract row to NA with n = 0, and the extract was built from a snapshot
  taken before the study-comparison rows were appended, so those rows never reached the CSV.
- `<dir>/scripts_or/current_status_regen.R`, `status_inventory.R`, `check_current_status.sh` and
  `<dir>/status_curated.md` (§3.4) — transplanted from `../continuous/scripts_mdsgnb20/`; syntax
  checked, **not yet executed**, so `current_status.md` does not exist.

## 9. Findings

1. **The projection exceeds the advance-go threshold** by a wide margin (19.48 h against 12 h), driven
   almost entirely by the three FS n = 2000 cells (8.47 h). The proximate cause is that this design's
   per-replicate cost at n = 2000 inflates 3.27× over the S0 one-worker reference once 63 workers
   compete for memory bandwidth, and the harm design point declares on 99.5% of replicates so the
   field block almost never short-circuits.
2. **Working vs effect scale is a real transplant hazard on a ratio measure.** Three separate checks
   inherited from the MD machinery were wrong on the OR path for one reason, and all three were silent
   failures in the sense that they would have passed as "findings" rather than as checker bugs. The
   record of which columns live on which scale is now in the template's summary-prep comment, in both
   checkers and in `COLUMNS_or.md`.
3. **[MD]'s band-width guard does not survive a changed default ε.** Its hard error was sound only
   because its default equalled the package default. Any future transplant that changes the ε default
   must revisit that guard (§3).
4. **The committed study records its oracle only on declared replicates**, because the refit sits after
   the detection return. Rule-independent though the oracle is, the study's bundles cannot be compared
   on it row-for-row without restricting to its declared rows (§5(a)).
5. **A pronounced tie regime on this design.** p̂(Ĥ) averages 0.12–0.19 and lies below 0.5 on 89–100%
   of declared replicates across all three identifiers at n = 500 — the selection is far from settled,
   which is exactly the regime the field construction is meant for, and worth carrying into the
   reading of Stage 2's bounds.
6. **An orphaned process from an earlier session** (§1), harmless but noted.
7. **`fs_dgm_scale()` has no ratio-measure path** and `fs_identification_summary()` requires a proxy
   covariate. Neither blocked anything here — both were dropped with reasons — but both are gaps if
   the binary path is to reach parity with the continuous one.

## 10. Commits

```
60ca1b8c  dev/tasks: the task document, as received, alone
0cc4c015  §1.3 the template
f25162f7  §1.4 scripts_or/ (sampler, smoke checker, gate2.R, runner)
7484e62b  §1.3 the band-width guard becomes a NOTE
d5835deb  §1.5 working vs effect scale: template ratios and both checkers
```

Stage 1's smoke and calibration payloads are untracked under `<dir>/mr_or_harm/` (600 KB) with raw
logs under `<dir>/logs_or/` (64 KB); §3.5's closeout, which deletes them, was not reached, and they
are the evidence behind §5 and §6.
