# REPORT — applied OC evaluation, stage 0 (2026-08-31)

**Task:** `dev/tasks/cc_task_oc_applied_stage0_2026-08-31.md` (commit `d4681023`)
**Verdict up front (§3.5): the sign gate BINDS — at every rung tested.**

## Provenance

- Host `pop-os` · `/home/larryleon/Documents/GitHub/forestsearch` · branch
  `feature/glm-extension` · installed version **0.3.1** ✓
- First attempt (HEAD `a2a661fb`) stopped at the §1 gate: three untracked
  paths (commit `c7413bdb` records it). Larry then instructed (by chat,
  2026-08-31): park them in `~/fs_parked_2026-08-31/` and re-run. Parked;
  two files inside `_payloads_2026-08-31/` turned out to be **tracked**
  (`selected_subgroups_continuous.rds`, the template payload) and were
  restored into the working tree via `git restore` (nothing lost — both
  survive in git and in the parked copy). Re-run HEAD: `6c071c74`, tree
  clean (`git status --porcelain` empty).
- Task document committed alone as `d4681023`; `~/Downloads` stem arrived
  hyphen-stripped as `cc_task_oc_applied_stage0_20260831.md`, committed as
  `cc_task_oc_applied_stage0_2026-08-31.md`.
- Script: `dev/glm-continuous-sims/stage0_oc_applied_2026-08-31.R`; full log
  `stage0_oc_applied_2026-08-31.log`; results object
  `stage0_oc_applied_2026-08-31.rds`. R 4.6.1.

## §2 The analyst's spec and the anchor

### §2.1 The spec, quoted from `quarto/applications/actg175/analysis_actg175_continuous_compare_all.qmd`

Twelve confounders — six continuous, six binary; **no `str2`** (:87–89):

```r
cont_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
bin_vars  <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
confounders.name <- c(cont_vars, bin_vars)
```

Outcome and thresholds (:94, :107–108, :135–139, :116, :127–128, :142, :83):

```r
adverse_outcome    <- "y_decline"   # FS adverse-scale outcome        (:94)
fs_md_threshold    <- 10    # screening: MD >= 10 indicates harm      (:107)
fs_md_consistency  <- 5     # consistency: MD >= 5                    (:108)
fs_consistency_method <- "resample"                                   (:116)
fs_conf.cont_jcuts <- list(age = 10, preanti = 10, wtkg = 10,
                           karnof = 10, cd40 = 10, cd80 = 10)         (:127-128)
fs_pconsistency    <- 0.90                                            (:135)
fs_n_min           <- 60L                                             (:136)
fs_d0_min          <- 10L                                             (:137)
fs_d1_min          <- 10L                                             (:138)
fs_splits          <- 500L                                            (:139)
fs_maxk            <- 2L                                              (:142)
analysis_seed <- 8316951L                                             (:83)
```

with `adverse_outcome = TRUE`, `cut_type = "default"`, `cont.cutoff = 4L` at
the call site (:320, :331, :348). The MR-anchor section's configuration
(:415–417, :432–434) — **the fixed-family configuration the OC evaluation
attaches to**, for the reason the document gives at :380–391 (MR and the
lasso/GRF front ends are mutually exclusive under `subgroup_method =
"consistency"`):

```r
sg_focus         = "maxeffCons",
selection_rule   = "neighborhood",
...
use_lasso        = FALSE,
use_grf          = FALSE,
use_dina         = FALSE,
```

### §2.2 Data prep, replicated in scratch (qmd :225–251)

Arms 1/3; `treat = 1` for ZDV+ddI; drop missing `cd420`;
`y_decline = cd40 − cd420`; six binaries as factors. **N = 1083.**

### §2.3 The anchor

Payload check: `_payloads_2026-08-31/analysis_actg175_continuous_compare_all/`
holds `selected_subgroups_continuous.rds` (tracked) — its row 8
(`maxeffCons / neighborhood`, subgroup `!{cd40 <= 507} & {age <= 37}`,
N_H 66, effect 87.91667, Pcons 0.95) is the **comparison run with
`use_lasso = use_grf = TRUE`** (qmd :144–145 via `...` pass-through :330–334),
not the fixed-family configuration. `comparison_continuous.rds` (untracked;
now parked) holds the same eight front-ends-on runs and no `fs_anchor` — the
qmd saves only `selected_tidy` and `comparison` (:725–728). **Not
recoverable from the committed payloads**, so the single authorized anchor
call ran: the MR-anchor call minus `mr_inference`, with
`parallel_args = list(plan = "sequential")`.

- **Wall-clock: 8.1 s** (sequential; the rendered document's multisession
  run with MR draws was 61.9 s).
- **Ĥ = `{age <= 37} & !{cd40 <= 507}`**, n(Ĥ) = **66**,
  **T̂_obs = 87.916667** (fitted MD), **p.consistency = 0.95**.
- Cross-checks: identical to the rendered HTML's anchor block
  (`definition : {age <= 37} & !{cd40 <= 507}`, `n selected : 66`) and to
  the payload's row-8 subgroup — the fixed-family and front-ends-on runs
  happen to select the same rule on this data.
- *GATE (an anchor exists, non-empty): PASS.*

## §3 The orientation gate

1. **ITT MD on `y_decline`** (replicated data): **−27.590878**
   (SE 7.8889, p = 4.9e-4). Negative — on average the combination arm
   *reduces* CD4 decline.
2. **Clause mapping** (verified subject-for-subject against the anchor's
   membership on the analysis frame — 66 of 66, exact):
   - `{age <= 37}` → `age = 37` (fixed cut, `process_cutpoint()` numeric
     path: `var <= 37`)
   - `!{cd40 <= 507}` → `cd40 = list(type = "greater", value = 507)`
     (`var > 507`)

   Anchored DGM: `generate_glm_dgm(model = "alt", k_treat = 1, n_super =
   5000, seed = 8316951)`, `factor_vars = bin_vars`, `continuous_vars =
   cont_vars`. The linear relation on the identity link:
   `m_tau[Q] = beta_treat + k_inter` with baseline-GLM
   `beta_treat = −26.978725`, so `k_inter(q) = q + 26.978725`. Readback at
   q = 20: `fs_dgm_scale` gives `m_tau[Q] = 20.0000000000` ✓.
3. **Signed table** (one build per rung; `m_tau[Qc] = beta_treat` is
   constant by construction — `k_inter` shifts Q only):

   | q | m_tau[Q] | m_tau[Qc] | sign(Q) == sign(Qc) |
   |---|---|---|---|
   | 2.5 | 2.5 | −26.978725 | **no** |
   | 5 | 5 | −26.978725 | **no** |
   | 7.5 | 7.5 | −26.978725 | **no** |
   | 10 | 10 | −26.978725 | **no** |
   | 15 | 15 | −26.978725 | **no** |
   | 20 | 20 | −26.978725 | **no** |
   | 30 | 30 | −26.978725 | **no** |
   | 40 | 40 | −26.978725 | **no** |
   | 87.916667 (T̂_obs) | 87.916667 | −26.978725 | **no** |

4. **Enumeration attempt at q = 20** — stop message captured verbatim:

   > the region effects m_tau do not share a sign; the oriented (abs) reading of the scale table is not defined.

5. ***Verdict, stated plainly: the sign gate BINDS — at every rung tested
   (q ∈ {2.5, 5, 7.5, 10, 15, 20, 30, 40, 87.916667}).*** On this applied
   DGM the complement's effect is the fitted (beneficial, negative)
   treatment effect while every planted Q effect is positive, so the
   opposite-sign state is not an edge case here — it is the generic state
   of the applied ladder. The queued signed-orientation task's premise is
   met.

## §4 The family and Q's representability (null DGM, n = 1083, max_M = 10000)

`generate_glm_dgm(model = "null", ...)` (same data, covariates, seed;
`flag_harm` all zero), `fs_oc_family_enumerate()` at the analyst's spec:

- **Stage counts:** cut columns **108** · enumerated **5886** · dropped:
  empty 255, minp 0, rmin 215, size 707 · kept **4709** · duplicates 201 ·
  **M = 4508**. Floors as resolved: n.min = 60 (Pg ≥ 0.0554), minp = 0.025,
  rmin = 5. Prevalence range 0.0556–0.9196.
- **Ĥ on `df_super`:** 317 of 5000 (P = 0.0634), by the reconstructed flag
  (memberships depend only on covariates, so the null build carries them).
- **Nearest member:** `age <= 37 & cd40 > 506` (Pg = 0.0636) — **max purity
  0.9969 and max Jaccard 0.9969**, both at the same member. Ĥ's `cd40`
  threshold (507) lies just off the population grid (nearest grid cut 506),
  so **no member is identical to Ĥ** — the pure discovery grid expresses
  the found rule to within ~0.3% of its membership, reported here as a
  feature of the evaluation; Q is not forced in.

## §5 Cost of the evaluation — measured

At the realized M = 4508, n = 1083, analyst's gate (c1 = 10, c2 = 5,
resample, pcons = 0.90):

- **Enumeration:** 43.4 s.
- **One `fs_oc_predict()` at draws = 2e4:** **400.1 s** wall-clock, peak
  memory ≈ **3.8 GB** (R gc, Ncells+Vcells max used). Sanity of the null
  run itself: detection rate 1.000 with E[beta(Hhat)] = 26.98 — under this
  null the common oriented effect (|−26.98|) clears c1 = 10 by
  construction, so the family always declares; the structural-null OC
  design question is stage 2's, not stage 0's.
- **Linear extrapolation to the 2e5 production draw set: ≈ 4001 s ≈ 67 min
  per draw set.**
- **Stage-2 projection** (assumptions stated: 10 rungs + null = 11 ladder
  sets; two sweeps at ~10 grid points each = 20 sets; inversions ~2 × 10
  evaluations = 20 sets; all at 2e5 draws): **(11 + 20 + 20) × 67 min ≈
  56.7 hours single-threaded.** M×M memory: **155 MB per M×M matrix**
  (overlap, covariance; the symmetric root is O(M³) ≈ 9.2e10 flops at this
  M). This projection is the compute go/no-go input for the evaluation
  task; per the task, no part of the evaluation itself was run. Note the
  per-set cost is CPU-bound and the host has 128 cores — if stage 2
  parallelizes across draw sets, wall-clock divides accordingly (≈ 51 sets
  → under 2 h at ~30 concurrent sets, memory permitting at ~4 GB/set).

## §6 The scale table (anchored DGM at q = 20)

σ = 122.810629 (σ² = 15082.45; the baseline GLM's residual SD on
`y_decline`).

| region | n_g | P_g | m_tau | V_eff |
|---|---|---|---|---|
| Q | 317 | 0.0634 | 20.000000 | 66739.290 |
| Qc | 4683 | 0.9366 | −26.978725 | 66423.382 |
| S | 5000 | 1.0000 | −24.000273 | 69277.855 |

**V_eff[Q] / V_eff[S] = 0.963357** — the effective variances are within 4%
of each other, so the prevalence-scaled `se_g` approximation (the §8
condition) is unlikely to matter much at this DGM.

## Ten-line summary

1. Spec replicated from the compare-all document; N = 1083; ITT MD =
   −27.59 (beneficial on average).
2. The fixed-family maxeffCons anchor was not in the committed payloads;
   the one authorized call ran in 8.1 s.
3. Ĥ = `{age <= 37} & !{cd40 <= 507}`, n = 66, T̂_obs = 87.92,
   p.consistency = 0.95 — matching the rendered document.
4. The clause mapping into `generate_glm_dgm()` reproduces Ĥ's membership
   exactly, 66/66.
5. m_tau[Qc] = −26.98 at every rung while every planted m_tau[Q] > 0:
   **the sign gate binds everywhere on the applied ladder** — the queued
   signed-orientation task runs.
6. The enumeration stop message was captured verbatim at q = 20.
7. The family at the analyst's spec has M = 4508; Ĥ's nearest member
   agrees to purity/Jaccard 0.9969; no member is identical (507 vs grid
   cut 506).
8. Enumeration costs 43 s; a 2e4-draw predict costs 400 s and ~3.8 GB;
   2e5 draws ≈ 67 min per set.
9. Stage-2 total ≈ **57 CPU-hours** at ~51 draw sets (assumptions in §5);
   155 MB per M×M matrix; parallelizing across sets brings wall-clock
   under ~2 h on this host.
10. Larry decides next: (a) the **extension** — the sign gate binds, so
    the queued signed-orientation task proceeds; (b) the **evaluation's
    compute go/no-go** — the §5 number: ≈ 57 CPU-hours (≈ 67 min per
    2e5-draw set at M = 4508, ~4 GB per concurrent set).

## Deviations

- First attempt stopped at the §1 clean-tree gate (report commit
  `c7413bdb`, superseded by this file). Larry instructed the three
  untracked paths be parked in `~/fs_parked_2026-08-31/` rather than
  deleted; done. Two tracked files inside `_payloads_2026-08-31/` were
  restored into the working tree with `git restore` after the move so the
  tree stayed clean without discarding tracked content.
- The task-document commit (`d4681023`) was not repeated on the re-run,
  per Larry's instruction.
- None otherwise: no `R/` change, no renders, no push, no evaluation run
  beyond §5's single timing call.
