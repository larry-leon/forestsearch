# CC task — the applied OC evaluation: type-I, the effect ladder, and the calibration curve for the ACTG175 continuous analysis

**File:** `dev/tasks/cc_task_oc_applied_evaluation_2026-08-31.md` · **Issued:** 2026-08-31 by chat, commissioned by Larry (compute go given on ≈ 12 CPU-h / ≈ 2 h wall-clock)
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1).
**Reads (in-repo):** `REPORT_oc_applied_stage0_2026-08-31.md` + `stage0_oc_applied_2026-08-31.rds` (the anchor, the clause mapping, the costs); `REPORT_oc_signed_orientation_2026-08-31.md` (0.3.2); `quarto/applications/actg175/analysis_actg175_continuous_compare_all.qmd` (the analyst's spec and payload conventions).

**What this builds.** The deliverable of the applied use case: the analyst's exact procedure — `y_decline`, `c1 = 10`, `c2 = 5`, `pconsistency = 0.90`, `"resample"`, the twelve confounders, jcuts 10 on all six continuous, `n.min = 60`, `maxk = 2`, at n = 1083 — evaluated analytically under truths anchored to the analysis data with **Q = Ĥ = `{age <= 37} & !{cd40 <= 507}`** (n(Ĥ) = 66, `T̂_obs = 87.916667`, from stage 0; never forced into the family — the discovery grid stays pure, nearest member purity 0.9969). One effect ladder from essentially zero through the observed statistic; type-I from the zero-plus rung; power at fixed type-I; what gets declared per rung; the calibration curve for the observed statistic; the `se_g` band. Output: precomputed `.rds` per rung plus **one application document** that reads them, renders fast, and writes its payload under the applications conventions.

**The structural null is the zero-plus rung, by design.** Stage 0 demonstrated that the null *branch*'s oriented-absolute reading turns this population's −27 benefit into a +27 "harm" (false-declaration 1.000 — nonsense here, correct for the MD40 fixtures whose drivers hand the search a positive-adverse outcome). No further `R/` change: the operational structural null is the alternative branch at `q = 0.01`, correctly signed under 0.3.2 — the same truth as "no subgroup" to a hundredth of a CD4 unit. The document explains this in one paragraph, citing stage 0's number as the demonstration; **the null branch itself is not edited and not used here.**

---

## ⚠ CATEGORY

**No `R/` change, no edit to any existing package file, driver or application document.** Writes: the precompute script and its per-rung `.rds`/logs under `dev/glm-continuous-sims/`; **one new application document** `quarto/applications/actg175/analysis_actg175_continuous_oc_evaluation.qmd`, its rendered output (follow that directory's existing tracking convention for HTML) and its payload under `_payloads/` per the conventions; this task document.

**Compute — the go/no-go, given:** eleven draw sets at ≈ 67 min each (stage 0's measured 2e4-draw cost × 10) ≈ **12.3 CPU-hours; run rungs as concurrent background processes (≤ 12 at once, ≈ 4 GB each) for ≈ 1.5–2.5 h wall-clock**, plus enumeration (~43 s per rung), the diagnostic, and a fast render. *GATE:* after the first rung completes, recompute the projection from its measured time; if the total projects above **6 h wall-clock**, stop and report.

**Unattended.** Gates stop the task; never ask, never work around.

---

## 1. Provenance and first commit — GATE (dirt-tolerant form)

Standard block plus the installed version (expect **0.3.2**). *GATE:* branch `feature/glm-extension`; the stage-0 report commit (`209a8e85`) and the signed-orientation commits (`eb136a35`, `15ff5640`) in the log. **A dirty tree is not by itself a failure**: list any untracked or modified paths, leave them alone, and work around them — *stop only if* dirt touches this task's own paths (`dev/glm-continuous-sims/oc_applied_*`, the new qmd's path, its `_payloads/` subdirectory, `dev/tasks/`). Copy this document into `dev/tasks/` and commit it alone; the `~/Downloads` stem arrives with hyphens stripped, so report both names.

## 2. Rebuild the anchored context — GATE

From `stage0_oc_applied_2026-08-31.rds` and the report, not from memory: the data prep (N = 1083), the spec, the clause mapping (`age = 37`; `cd40 = list(type = "greater", value = 507)`), `beta_treat = −26.978725`, `k_inter(q) = q + 26.978725`, `T̂_obs = 87.916667`. Re-verify once, cheaply: the reconstructed flag matches the anchor membership 66/66 on the analysis frame. *GATE:* exact.

## 3. The rungs — build and gate

**q ∈ {0.01, 2.5, 5, 7.5, 10, 15, 20, 30, 40, 60, 87.916667}** — eleven rungs; `q = 0.01` is the zero-plus structural null (a planted effect of exactly zero is rejected by 0.3.2's orientation, deliberately), and the top rung is the observed statistic. Per rung: `generate_glm_dgm(model = "alt", k_treat = 1, n_super = 5000, seed = 8316951, k_inter = k_inter(q))` with the stage-0 covariate lists and subgroup cuts; then `fs_oc_family_enumerate()` at n = 1083, `max_M = 10000`.

*GATES per rung:* `abs(m_tau[Q]) − q` under 1e-9 and `m_tau[Qc] = −26.978725` unchanged (`fs_dgm_scale` readback); `orientation$s = +1`; M = 4508; `lab`, `Pg`, `PQg`, `se_g`, `sens_g`, `spec_g`, `ovl`, `memb` `identical()` to the first rung's; and `beta_g(q) − beta_g(0.01)` equal to `(q − 0.01)·PQg` to 1e-9. Any other movement means something depends on `k_inter` that source says does not — stop.

## 4. The evaluation — one draw set per rung

Per rung, `fs_oc_grid(family = <rung family>, n = 1083, c1 = sort(c(0:200, 87.916667)), c2 = 5, consistency_method = "resample", pconsistency = 0.90, draws = 2e5, block = 5e4, seed = 8316951)` — one draw set, every threshold swept against it. Launch rungs as concurrent detached processes (≤ 12), each logging; checkpoint one `.rds` per rung holding the grid `table` plus the full `fs_oc_predict` objects at the named points only — the analyst's `c1 = 10`, `c1 = 87.916667`, and (once known) `c1_05` and `c1_10`. **Do not save per-threshold `results` lists.** A rung that dies re-runs alone.

Then, from the stored tables:

1. **Type-I:** the zero-plus rung's `det_rate` against `c1`; `c1_05` and `c1_10` = the smallest integer `c1` with rate ≤ 0.05 / 0.10 (1-unit resolution is the document's resolution — state it); and the headline: the rate at the analyst's `(c1 = 10, c2 = 5)`.
2. **The ladder:** per rung, `det_rate` at the analyst's `c1`, power at `c1_05` and `c1_10`, and the smallest integer `c1` with rate ≥ 0.80 / 0.90 / 0.95 (read from the grid; no separate inversion calls — the one-block `fs_oc_invert` at this M costs ~7 GB and adds nothing at this resolution).
3. **What is declared:** per rung at the analyst's `c1`: `EnH`, `Eppv`, `Esens`, `Espec`, `EbetaH`, `Enaive_bias`, `mass_below` with MC SEs.
4. **The calibration curve:** `P(T ≥ T̂_obs | q)` = the grid's `det_rate` at `c1 = 87.916667`, per rung — no extra machinery. Report it as a table and mark where it crosses 0.05, 0.5, 0.95 (interpolating between rungs, stated as interpolation).
5. **Coherence checks, reported:** the ladder's power columns monotone in `q`; the zero-plus rung's declared rules (if any mass) carry `EbetaH ≈ −27` — the honest false-positive reading.

## 5. The `se_g` band — three rungs

At q ∈ {0.01, 20, 87.916667}: `fs_dgm_scale(dgm, regions = <the 4508 memb columns>)`, `se_direct = sqrt(V_eff[g]/(n·Pg))`, ratio to `se_g`; range, median, and the three purity bands (`PQg ≥ 0.95`, mid, `< 0.25`), per rung. Reported, not acted on; if the band at the top rung exceeds ~5%, add one sensitivity run there (the family with `se_g` replaced by `se_direct`, one draw set, the §4 outputs beside the primary) and label it exactly as the breadth ladder did: a sensitivity, not a constructor, not adopted.

## 6. The application document

`quarto/applications/actg175/analysis_actg175_continuous_oc_evaluation.qmd` — **reads the stored `.rds` only, computes nothing heavy, types no number** (inline R throughout). Sections: (1) the analysis and its anchor — the fixed-family configuration and why the OC evaluation attaches to it (the document's own MR-anchor rationale), Ĥ, `T̂_obs`, consistency; (2) the anchored DGM — the mixture stylization stated plainly, the clause mapping, and Q's representability (purity 0.9969, `cd40` 507 vs grid 506, not forced); (3) type-I under the structural null, with the one-paragraph orientation explanation and stage 0's 1.000 as the demonstration of why the null branch is not used; (4) the ladder — power at fixed type-I, the declaration curves figure (rungs + zero-plus, guides at `c1_05`, the analyst's `c1`, and 0.80); (5) what is declared per rung; (6) the calibration curve — figure with the observed `T̂_obs` line, and the model-based caveat verbatim: *a consonance curve for the planted effect under the anchored two-point mixture model; not a posterior, not adjusted for anything outside the model*; (7) the `se_g` band; (8) **Limitations** — carry, adapted to this document: the between-rule `E|Ĥ|` gap (~2 subjects at MD40/MD120, not scaling with the harm); Q-not-in-family with this analysis's numbers; the mixture stylization; the fixed-family conditioning (front ends off). Payload chunk per the conventions (`table labels meta extras est_scale built_at forestsearch_version`, `est_scale = "md"`, `results_dir`/`dirout` NULL resolution). Render (quarto path per the repo's convention); the render must be fast (< 2 min) since it reads stored objects — a chunk that recomputes a draw set is a defect.

## 7. Close-out

- Commit by explicit path: the precompute script, the per-rung `.rds` (if any single file exceeds 20 MB, leave it untracked and say so), the qmd, its rendered output per the directory's convention, the payload, and the report. **No push. No install. No `R/` change.**
- Report `REPORT_oc_applied_evaluation_2026-08-31.md`: provenance · §2–§3 gate arithmetic · the type-I table · the ladder table · the calibration table · the `se_g` bands (and the sensitivity if it ran) · per-rung wall-clocks and the concurrency used · the render check · commits · ten-line summary ending with the three numbers Larry will read first: type-I at the analyst's `(10, 5)`; the smallest `q` with 80% power at `c1_05`; and `P(T ≥ 87.92 | q)` at the rung nearest the crossing.

## 8. Out of scope

- No `R/` change — in particular no edit to the null branch; its orientation limitation is documented prose here and a recorded item for the record, not a task.
- No simulation replicates, no gbsg (blocked on a survival scale layer), no second population, no edit to the compare-all or template documents.
- No `fs_oc_invert()` runs (grid-resolution inversions only, stated); no fine-grid calibration beyond the eleven rungs.
- No push, no install.
