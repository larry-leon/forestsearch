# CC task — applied OC evaluation, stage 0: feasibility gates for anchoring the wrapper on the ACTG175 continuous analysis

**File:** `dev/tasks/cc_task_oc_applied_stage0_2026-08-31.md` · **Issued:** 2026-08-31 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1).
**Reads (in-repo):** `quarto/applications/actg175/analysis_actg175_continuous_compare_all.qmd` and its `_payloads/`; the continuous template document beside it; `R/fs_oc_family.R`, `R/fs_dgm_scale.R`, `R/generate_glm_dgm.R`.

**What this is.** Larry's decided use case: an analyst runs `forestsearch()` on the ACTG175 continuous analysis (adverse scale `y_decline`, thresholds `c1 = 10`, `c2 = 5`) and finds a subgroup Ĥ; the operating characteristics of that exact procedure are then evaluated with `fs_oc_predict()`/`fs_oc_grid()` under truths anchored to the analysis data — a ladder of planted effects for **Q = Ĥ** (decided: the found subgroup, its clauses over the raw covariates; **never forced into the enumerated family** — the family stays the pure discovery grid, and how well it can express Ĥ is reported, not adjusted), plus the structural null. This stage runs the feasibility gates and measures the costs; the evaluation itself is a later task with its own compute go/no-go built on this stage's numbers.

---

## ⚠ CATEGORY

**No `R/` change, no edit to any package file, driver or application document.** Writes: scratch scripts, their `.rds`/`.log`, and the report, under `dev/glm-continuous-sims/`; plus this task document. **Bounded compute authorized:** at most one fixed-family `forestsearch()` anchor call (§2, only if the committed payloads do not hold the anchor), one family enumeration, a handful of DGM builds, and one small timing draw (§5). **Estimate 30–45 minutes.** No simulation study, no renders.

**Unattended.** Gates stop the task; never ask, never work around.

---

## 1. Provenance and first commit — GATE

Standard block (`hostname; pwd; git rev-parse --abbrev-ref HEAD; git rev-parse --short HEAD; git status --porcelain; git log --oneline -4`) plus the installed version (expect **0.3.1**). *GATE:* branch `feature/glm-extension`, clean tree. Copy this document into `dev/tasks/` and commit it alone; the `~/Downloads` stem arrives with hyphens stripped, so report both names.

## 2. The analyst's spec and the anchor — from the repo, quoted

1. **The spec.** From the compare-all document, quote the parameter block verbatim: the twelve confounders (six continuous, six binary — note: no `str2`), outcome `y_decline` with `adverse_outcome = TRUE`, `effect.threshold = 10`, `consistency.threshold = 5`, `pconsistency.threshold`, `conf.cont_jcuts` (all six continuous at 10), `cut_type`, `maxk`, `n.min`, `d0.min`/`d1.min`, `fs.splits`, `consistency_method`, `seedit`. And the MR-anchor section's configuration (`sg_focus = "maxeffCons"`, `selection_rule = "neighborhood"`, `use_lasso = use_grf = use_dina = FALSE`) — **that fixed-family configuration is the one the OC evaluation attaches to**, for the same reason the document itself gives for MR.
2. **The data prep**, replicated exactly in scratch (arms 1/3, treat coding, drop missing `cd420`, `y_decline = cd40 − cd420`, binaries as factors). Report N.
3. **The anchor Ĥ.** Check the document's committed payload directory (`_payloads/analysis_actg175_continuous_compare_all/`: `selected_subgroups_continuous.rds`, `comparison_continuous.rds`, and anything else there) for the fixed-family maxeffCons result — the subgroup definition, its size, its fitted MD (call it `T̂_obs`), `p.consistency`. If recoverable, read it and say from which object. **If not recoverable, run the anchor once** — the MR-anchor call minus `mr_inference` (not needed here), `parallel_args = list(plan = "sequential")` — report its wall-clock, and record Ĥ, n(Ĥ), `T̂_obs`, `p.consistency`. *GATE:* an anchor exists by one route or the other; if the anchor declares nothing (no subgroup), stop and report — the design's Q needs a found subgroup.

## 3. The orientation gate — the numbers Larry asked for

1. The fitted ITT MD on `y_decline` (value and sign), from the replicated data.
2. Express Ĥ's clauses as `subgroup_vars` / `subgroup_cuts` for `generate_glm_dgm()` (the flexible cut format; quote the mapping and verify on the analysis frame that the reconstructed flag matches Ĥ's membership from the anchor exactly — subject-for-subject). Build the anchored DGM at `model = "alt"`, `k_treat = 1`, `n_super = 5000`, `seed = 8316951`, with `k_inter` set so the planted Q effect on the adverse scale is `q = 20` (the linear relation; verify via `fs_dgm_scale()` readback).
3. Report `m_tau[Q]` and `m_tau[Qc]` (signed) from `fs_dgm_scale(dgm)` at `q ∈ {2.5, 5, 7.5, 10, 15, 20, 30, 40}` and at `q` equal to `T̂_obs` — cheap, one build per rung — and state per rung whether `sign(m_tau[Q]) == sign(m_tau[Qc])`.
4. Attempt `fs_oc_family_enumerate()` once at a rung where the signs differ (if any) and capture the stop message verbatim.
5. *Verdict, stated plainly:* **the sign gate binds** (list the rungs where) **or does not.** This verdict is what the queued follow-on task keys on.

## 4. The family and Q's representability — via the null DGM, no `R/` change needed

The enumeration's memberships depend only on `df_super`'s covariates, not on the planted effect, so build the **null** DGM (`model = "null"`, same data, same covariates, same seed) and enumerate the family at the analyst's spec and the trial's N (the analysis N from §2 — state it), with `max_M` raised as needed (say 10000, deliberately). Report:

- the stage counts (cut columns, enumerated, each floor, kept, duplicates, **M**);
- Ĥ's membership evaluated directly on `df_super` (the reconstructed flag), and against every family member: the maximum purity `P(g ∩ Ĥ)/P(g)`, the maximum Jaccard with Ĥ, and the nearest member's label and prevalence. **Q is not forced in** — this is the measurement of how well the pure discovery grid expresses the found rule, reported as a feature of the evaluation;
- whether any member is *identical* to Ĥ (possible if Ĥ's thresholds happen to lie on the population grid — say which either way).

## 5. Cost of the evaluation — measured, not guessed

At the realized M: time the enumeration itself; then one `fs_oc_predict()` on the null family at the analyst's `(c1, c2)`, `"resample"`, with `draws = 2e4` and the default block, recording wall-clock and peak memory; extrapolate linearly in draws to the 2e5 production draw set. Project stage 2's total — roughly (8–10 rungs + the null) draw sets plus two sweeps and the inversions — as hours, with the M×M memory stated. **This projection is the compute go/no-go input for the evaluation task; do not run any part of the evaluation itself.**

## 6. The scale table

From the §3 anchored DGM at one rung: σ (the baseline GLM's residual SD on `y_decline`), the regions table (`Q`, `Qc`, `S`: `P_g`, `m_tau`, `V_eff`), and the ratio `V_eff[Q] / V_eff[S]` — the quick read on how much the prevalence-scaled `se_g` approximation is likely to matter at this DGM (the §8 condition travels here too).

## 7. Report

`dev/glm-continuous-sims/REPORT_oc_applied_stage0_2026-08-31.md`: provenance · §2's spec quotations and the anchor (source, Ĥ, n(Ĥ), `T̂_obs`, `p.consistency`) · §3's signed table and the verdict, with the stop message if captured · §4's family counts, M, and the purity/Jaccard of the nearest member · §5's measured costs and the stage-2 projection in hours · §6's scale table · ten-line summary ending with the two things Larry decides next: the extension (handled by the queued task) and the evaluation's compute go/no-go (the §5 number). Commit scripts, outputs and report by explicit path. **No push. No install. No `R/` change.**

## 8. Out of scope

- No edit to `R/`, to any application document, driver or payload. No forcing of Ĥ into any family.
- No evaluation run beyond §5's single timing call; no ladder, no sweeps, no calibration curve — those are the next task, after the go/no-go.
- No simulation replicates, no renders, no push.
