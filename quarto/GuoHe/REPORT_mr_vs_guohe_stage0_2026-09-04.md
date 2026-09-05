# REPORT — MR gate vs Guo & He (2021): Stage 0 (discovery)

Date: 2026-09-04. Task: `dev/tasks/TASK_mr_vs_guohe_2026-09-04.md` (committed `92e3b741`).
Repo state at discovery: branch `feature/glm-extension`, HEAD `92e3b741cacee00f29c130e9f2770cf26b905aa0`, forestsearch 0.3.5.
Decisions: Larry's message left D1–D5 as the placeholder (`[fill in, or "defaults"]`); Stage 0 needs no decision, and this report records what each default would mean against the discovered source. No compute was run; no file under `R/` was touched. `.rds` bundles were opened read-only for structure inspection.

## GATE 0: GREEN

All of 0a–0d resolved from source. None of the four stop conditions holds:
the replication source is unambiguous (two complementary, both-current documents covering disjoint scenario sets); the MR gate is callable on an external family with a pure argmax and no screen, arithmetic untouched; the continuous path is not needed (both replications are survival) and is present anyway; per-replicate reproduction is possible in every cell (per-replicate rows are stored outright, plus deterministic seed regeneration).

Two findings that are not stop conditions but need resolution before Stage 2 (§ Findings): the MR engine does not return per-candidate re-selection frequencies, and `M_eff` / the "(A6) diagnostic" named in Stage 2 have no in-repo definition.

---

## 0a — The replication source

`grep -ril "guo" quarto/ dev/ R/ vignettes/ inst/ tests/` resolves to `quarto/GuoHe/` (plus archived/legacy hits under `quarto/archive/` and the MR/dina cross-references). Two replication lineages exist there, **both current and complementary** — they cover disjoint sections of Guo & He (2021), JASA 116(535):1498–1506:

| Lineage | Document | Harness / driver | Bundles | Committed |
|---|---|---|---|---|
| §5.1, Tables 3–6 (disjoint families) | `quarto/GuoHe/guohe_reproduction.qmd` | `guohe_reproduction_sim.R`, `guohe_reproduction_run.R` | `guohe_repro_t35_beta2_00..05.rds` (6), `guohe_repro_t6_k02/06/10/12.rds` (4) | `8a98e05d` "GuoHe Section 5 reproduction: harness, 2000-rep results, and comparison qmd" |
| §5.2, Table 7 (nested post-hoc family) | `quarto/GuoHe/guohe_sec52.qmd` | `guohe_sec52_sim.R`, `guohe_sec52_run.R`, `guohe_sec52_truth.R` | `guohe_repro_t7_beta2_00..05.rds` (6) + truth caches `guohe_sec52_truth_beta2_00..05.rds` (6) + `sens_B500/guohe_repro_t7_beta2_01.rds` (B=500 sensitivity) + `guohe_repro_t7_beta2_00_pilot.rds` | `a9e099df` "GH reproduce Sec 5.2" |

`guohe_sec52.qmd:9` (header comment in `guohe_sec52_sim.R:8-10`) states the relationship explicitly: "The Section 5.1 reproduction (Tables 3-6, disjoint families) is complete and committed. Section 5.2 is the nested regime closest to the `forestsearch` maxeff family." There is no ambiguity about currency: neither supersedes the other; together they are "the replicated Guo–He simulations" (16 scenario cells). `test_guohe_algorithm3*.qmd`, `test_guohe_comparator.qmd` and the `quarto/archive/` hits are estimator-level validation / legacy, not the replication.

**Match record.** Both qmds embed the published targets in the sim files (`GH_TABLE_3/5/6` at `guohe_reproduction_sim.R:246-280`; `GH52_TABLE_7` at `guohe_sec52_sim.R:96-103`, "never retyped") and render reproduced-vs-published with combined Monte Carlo SE and per-cell flags. Verdicts, from source:

- §5.1 (`guohe_reproduction.qmd:264-277`): "every cell reproduces within combined Monte Carlo error", licensing "the estimator (`guohe_algorithm3` / `guohe_adaptive_r`) is faithfully implemented on disjoint, equal-probability families up to k = 12".
- §5.2 (`guohe_sec52.qmd:451-476`): "Twenty-six of thirty cells reproduce within two combined Monte Carlo standard errors; the four exceptions are the proposed columns of the β₂ = 0.1 row. The naive column ... agrees with published throughout"; proposed columns carry "a one-sided offset of order 0.010 ... with B excluded by direct experiment"; licenses faithful implementation "on nested, order-statistic-indexed families of ≈151 heavily-overlapping candidates at n = 400, with a small unexplained conservatism deficit in the bound".

Rendered `guohe_reproduction.html` / `guohe_sec52.html` exist locally but are gitignored (`quarto/.gitignore`: `GuoHe/*.html`), as are the run logs (`GuoHe/*.log`). The 16 result bundles and 6 truth caches **are** git-tracked (28 tracked `.rds` under `quarto/GuoHe/`).

## 0b — Replication design, from source

Common to both sections (quoted from `guohe_reproduction_sim.R` and `guohe_sec52_truth.R`/`_sim.R`):

- **Outcome type: survival** (Cox). Hazard λ₀(t)·exp(β·D), λ₀ unit exponential; D ~ Bernoulli(0.5); censoring log C ~ Unif(−1.25, 1.00), ~40–41%. `gh_sim_data()` at `guohe_reproduction_sim.R:63-77`; `gh52_sim_data()` in `guohe_sec52_truth.R`.
- **Within-subgroup estimator:** treatment-only Cox regression coefficient per candidate, `coxph(Surv(time, event) ~ treat)` on the subgroup rows (`gh_subgroup_fits()` :88-102; `gh52_subgroup_fits()` `guohe_sec52_sim.R:132-146`), model-based SE for the naive CI.
- **Selection rule: pure argmax**, no threshold. Naive: `which.max(...)` (`guohe_reproduction_sim.R:108`, `guohe_sec52_sim.R:151`). The engine's selection is r-invariant full-sample argmax, asserted per replicate (`guohe_sec52_sim.R:219-230`).
- **Guo–He bootstrap size:** B = 2000 in every committed cell (inspected in all 16 bundles). The paper is silent; 2000 is recorded in each bundle as an inference (`guohe_reproduction_run.R:25-26`). Level 0.05 one-sided; r grid {1/3, 1/12, 1/21, 1/30}; `min_events = 5`. Common random numbers across the r grid (`boot_seed = seed + 500000L`, `guohe_reproduction_sim.R:150-157`, `guohe_sec52_sim.R:199`).
- **Tuning-parameter rule:** fixed-r columns everywhere; the adaptive Algorithm 2 (`guohe_adaptive_r`, v = 5, adaptive B = 2000) runs in §5.1 only (`ad_*` columns in t35/t6 bundles); Table 7 has no Adaptive column by design.
- **Replicates:** 2000 per scenario, complete in all 16 cells (`n_rep_requested == n_rep_used == 2000` verified in every bundle — so row m ↔ seed `base + m` alignment holds with no dropped replicates).
- **Seed scheme:** `RNGkind("Mersenne-Twister", "Inversion", "Rejection")`; replicate seed = scenario `seed_base + m`, `set.seed(seed)` at the top of the one-rep function, data drawn first, then `boot_seed = seed + 500000L` for the engine. §5.1 base = `1000000L + as.integer(sum(utf8ToInt(id)) * 997L)` (`guohe_reproduction_run.R:100`); §5.2 base = `1000000L + as.integer(sum(utf8ToInt(id))) * 100003L` (`guohe_sec52_run.R:155`). Verified against stored `seed_base` (e.g. t7_beta2_00: 92,902,757 = 1e6 + 919·100003). Known, documented quirk: §5.1's 997 spacing makes adjacent t35 scenarios share 1003/2000 replicate seeds (coupled datasets across β₂; per-cell unbiased, joint reading correlated — `guohe_sec52_followup2_CC_NOTE.md` §1, "known, benign feature of the §5.1 results; do not retrofit or re-run"). §5.2 has zero sharing.

Section-specific:

- **§5.1**: K disjoint subgroups, membership equal-probability 1/k, columns `S1..Sk`. Tables 3–5: k = 2, n = 400, β₂ ∈ {0,…,0.5}. Table 6: k ∈ {2,6,10,12}, n = 200k, all β = 0. **Orientation: `orient = -1`** — inferential effect is the negated log-HR, best subgroup = argmax(−coef) (`guohe_reproduction_sim.R:34-39,163`). Truth θ_{k̂} = −β[selected], known exactly.
- **§5.2**: nested family S(c) = {W ≤ c}, c ∈ [30,60], data-determined order-statistic cutpoints {30} ∪ {Wᵢ ∈ (30,60]}, expected size 151, random per replicate, materialized as 0/1 columns `cand_001…` (`gh52_candidates()` `guohe_sec52_sim.R:116-123`). **Orientation: `orient = +1`** — argmax of the raw coefficient (resolution recorded at `guohe_sec52_sim.R:29-41`). Truth θ_{k̂} = β(ĉ) from the gate-cleared truth caches (n_big = 2e6, c_step = 0.25, seeds 20260721+…, all gates PASS per `guohe_sec52_production.log`); exact β(c) ≡ 0 at β₂ = 0 (`beta_exact` scoring basis), isotonic-smoothed simulated curve otherwise, with recorded residual anchor offsets +0.0030/+0.0017/−0.0008/−0.0012/−0.0010 for β₂ = 0.1…0.5.

**Per-replicate storage** (from `readRDS` of the bundles):

- t35/t6 rows: `k, n, beta2, cens_rate, naive_sel, naive_beta_s, naive_cover, naive_dist, naive_bias, r1..r4_{cover,dist,bias}, ad_r, ad_cover, ad_dist, ad_bias`. Naive point/lower are exact arithmetic on stored columns (point = `naive_bias + naive_beta_s`; lower = `naive_beta_s − naive_dist`). **The engine's selected index is NOT stored** (r-column cover/dist/bias reference the engine-selected truth `truth[sel_i]`, `guohe_reproduction_sim.R:173-179`); see 0d for recovery.
- t7 rows: `beta2, n, n_cand, cens_rate, c_hat_gh, c_hat_naive, sel_agree, n_sel, gamma_s, naive_point, naive_lower, naive_cover, naive_dist, naive_bias, gamma_s_naive, r1..r4_{cover,dist,bias}, rep_elapsed_sec`. G&H per-replicate estimate and bound are exact arithmetic on stored columns: debiased = `r*_bias + gamma_s`; lower = `gamma_s − r*_dist`; selected cutpoint `c_hat_gh` stored.
- Every bundle also stores `seed_base`, `B`, `r_grid`, `elapsed_sec`, full `sessionInfo`, and the cell summary; t7 additionally `orient`, `truth_identity`.

**Wall-clock** (production logs, 120 forked workers on this host; logs untracked but present):

| Cell group | wall/cell | per-replicate |
|---|---|---|
| t35 (k=2, incl. adaptive at B=2000) | 16.2–16.4 min | ≈ 58 core-s |
| t6 k=6 / k=10 / k=12 | 48.9 / 82.2 / 100.2 min | ≈ 176 / 296 / 361 core-s |
| t7 (≈151 candidates, 4 r's, B=2000) | 192.8–196.1 min | 682–693 s (stored `rep_elapsed_sec`) |

## 0c — The MR gate

**Entry point:** `fs_mr_inference()` at `R/fs_mr_inference.R:352` (internal, `@keywords internal`; last touched by `8694e3dc` "Expose mean_r / mean_r_c…"). Signature, quoted:

```r
fs_mr_inference <- function(df, candidates, spec, selected_members,
                           admission,
                           t_confirm = NULL, confirm_rule = c("point", "ci"),
                           reselection = c("maxcons", "maxeff", "maxSG",
                                           "minSG", "effMaxSG", "effMinSG"),
                           effect_neighborhood = 0.10,
                           selection_rule = c("neighborhood", "pareto", "both"),
                           draws = 2000L,
                           multiplier = c("poisson", "gaussian", "rademacher"),
                           include_complement = FALSE,
                           ci_method = c("ij", "wald"),
                           seed = NULL)
```

Component map, all in `R/fs_mr_inference.R` unless noted:

- Effect vector + influence matrix: `.fs_mr_assemble()` (:88) — per candidate, `.fs_mr_pieces()` (:61) dispatches to `.consistency_cox_pieces()` (`R/consistency_resample.R:55`) for survival: one `coxph(Surv ~ treat)` fit, `beta_hat` = raw treat log-HR, `dfbeta = residuals(fit, type="dfbeta")`, `sigma_D = sqrt(sum(dfbeta^2))` (the Lin–Wei robust SE). Candidates with < 6 members are dropped (:95) — irrelevant here (family sizes ≥ ~150).
- Multipliers: `.fs_mr_multipliers()` (:110); `"poisson"` = `rpois(n·draws, 1) − 1`, i.e. centred Poisson (the D3 default).
- Re-selection per draw: `.fs_mr_select()` (:132); `maxeff` = `passers[which.max(beta[passers])]` (:171) — pure argmax on the perturbed coefficient.
- Bias terms: `selection_bias = mean(sel_bias)` over winner draws `ok_H`, `fixed_bias = mean(P[sel, ok_H])` (:479-483); `beta_deb = beta_naive − selection_bias − fixed_bias` (:485).
- IJ variance with finite-B correction: `.fs_mr_ij_var()` (:195), `hat_V = tilde_V − (N/B_ok)·mean(rb²)` (:204); SE resolution with fallback chain ij → ij_raw → wald in `.fs_mr_se_from_ij()` (:219); intervals two-sided 95% + one-sided `lower_1s` (:502, :591-597). Return carries `se_ij, se_wald, var_ij, ij_source, ij_draws, selection_bias, fixed_bias, selection_rate, mean_r, timing_seconds`.

**External family: YES.** `candidates` is "a named list of integer row-index vectors, one per screened candidate subgroup" (:246-248) — no forestsearch object required. The docs warn direct calls bypass `.validate_mr_configuration()` and require the family to be ranked on the inferential coefficient (:330-345) — exactly satisfied here: both replications select the argmax of the same Cox coefficient MR perturbs.

**Pure argmax, no screen: YES.** `reselection = "maxeff"` plus `admission = list(effect_floor = NULL, consistency = NULL)` takes the unrestricted path — "every estimable candidate is admissible" (:424-429). No arithmetic change needed. (`"maxcons"` without a consistency floor is refused at :435-441; not our rule.)

**Orientation: no `orient` argument exists.** `maxeff` maximizes the raw coefficient. §5.2 (orient = +1, argmax raw coef) is directly compatible. §5.1 (orient = −1, argmax of the negated coef) requires an adapter-level data recode — swapping the treatment labels (`treat ↦ 1 − treat`) negates the treatment-only Cox coefficient and its dfbeta exactly, making `beta_hat` the oriented effect with the gate's arithmetic untouched. Stage 1 design point, recorded here.

**Outcome type: survival path present** (`spec$outcome_type = "survival"` → `.consistency_cox_pieces`, :62-65). The continuous/MD path is **not needed** (both replications are survival) and is present on this branch anyway (`R/consistency_resample.R:95` "GLM pieces (binary OR/RR/RD, continuous MD, count IRR)").

## 0d — Per-replicate reproduction

**Possible in every cell**, by two independent routes:

1. **Stored outright** (§5.2): the t7 rows carry naive point/lower and gamma_s, c_hat, and per-r cover/dist/bias from which the G&H debiased estimate and lower bound are exact arithmetic (0b). Nothing needs re-running for the §5.2 G&H and naive columns.
2. **Seed regeneration** (both sections): each replicate's dataset regenerates bit-for-bit from `set.seed(seed_base + m)` followed by the harness's draw sequence (RNGkind pinned; drivers state "a parallel run reproduces a serial run bit-for-bit", `guohe_sec52_run.R:27-30`; replicate counts are complete so row m ↔ seed alignment is intact). This is what MR needs to run on the identical replicates.

§5.1 gap and its closure: t35/t6 rows do not store the engine's selected index, so the per-replicate G&H debiased/lower cannot be read off the stored columns alone in t35 (in t6 all truths are 0, so debiased = `r*_bias`, lower = `−r*_dist` directly). Recovery without any bootstrap: regenerate the replicate's data from its seed, refit the k per-candidate Cox models, take the (r-invariant) argmax — then `truth[sel]` is known and debiased/lower follow from the stored bias/dist columns exactly. Cost ≈ k coxph fits per replicate (trivial; part of the paired Stage 2 pass anyway).

**D5 is therefore NOT triggered:** no re-run of `guohe_algorithm3()` is needed to obtain per-replicate G&H values or CIs in any cell.

## Decisions D1–D5 against the discovered source

- **D1 (default: all scenarios)** = 16 cells: t35 ×6, t6 ×4, t7 ×6. Note the §5.1 t35 seed-coupling quirk above affects joint readings across β₂ rows there; §5.2 is clean. §5.2 is the regime closest to the forestsearch maxeff family.
- **D2 (default)**: pure argmax on the replication's estimator = `reselection = "maxeff"`, unrestricted admission; orientation per section as in 0c (§5.1 via treat recode, §5.2 native). "Second MR row using the package's default gate": the package default `reselection = "maxcons"` is impossible on this family (no consistency floor → refused, `R/fs_mr_inference.R:435`), so the default of "no" is also the only coherent choice.
- **D3 (default)**: `draws = 5000L, multiplier = "poisson"` — both supported engine arguments (defaults are 2000/poisson).
- **D4 (default: no)**: the replication's existing bootstrap values are not stored per draw in the bundles; a full-bootstrap row is not free. Consistent with the default.
- **D5**: moot — per-replicate G&H values recoverable without re-running (0d).

## Findings requiring resolution before Stage 2 (not Gate 0 stops)

1. **Re-selection frequencies are not exposed.** Stage 2 asks per replicate for "re-selection frequencies p̂ for the top three candidates". `fs_mr_inference()` computes the per-draw `winner` vector (:449-457) but returns only `selection_rate`; the winner counts are discarded. Options for Stage 1, in the task's preference order: (a) the one permitted R/ addition (`mr_gate_fixed_family()`) calls the existing internals (`.fs_mr_assemble`, `.fs_mr_multipliers`, `.fs_mr_select`, `.fs_mr_ij_var`) directly — unchanged — and keeps the winner tabulation; or (b) an add-only optional argument on `fs_mr_inference()` returning the winner table (existing behaviour unchanged at the default). Larry's call at Gate 1 if (a) is not preferred.
2. **`M_eff` and the "(A6) diagnostic" have no in-repo definition.** `grep -rn "M_eff|(A6)" R/ dev/` returns nothing. These Stage 2 record fields come from the task/companion-paper notation, not from source; their definitions must be supplied (or pointed at) before Stage 2's record schema is final. Not needed for Stages 0–1.
3. §5.1 t35 G&H per-replicate values require the cheap selection-recovery pass described in 0d — to be folded into the Stage 2 paired run and verified against the stored `r*_cover` indicators as a bit-level cross-check.

## Stage 1 preview (no action taken)

The gate accepts the external family, so per Stage 1a the adapter belongs in the comparison `.qmd` with nothing under `R/` — except that Finding 1 may justify the one permitted R/ function to expose winner frequencies. Identity 1b's variance target: `sigma_D² = sum(dfbeta²)` is the Lin–Wei robust variance by construction (`R/consistency_resample.R:86`), which is what the sandwich-variance identity will check against `coxph(..., robust = TRUE)`.

— End of Stage 0. Stopped at Gate 0 per task protocol.
