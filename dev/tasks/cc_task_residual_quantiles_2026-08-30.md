# CC task — the last residual lead: population-quantile cuts against replicate-quantile cuts

**File:** `dev/tasks/cc_task_residual_quantiles_2026-08-30.md` · **Issued:** 2026-08-30 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1).
**Reads:** `REPORT_oc_wrapper_sgdef_2026-08-29.md` + `sgdef_selection_2026-08-29.rds`; `REPORT_oc_wrapper_confs_2026-08-30.md` + `oc_wrapper_grid_corrected_2026-08-30.rds`; `REPORT_residual_2026-08-30.md`; `REPORT_fs_family_report_2026-08-30.md`.

**Where the question stands.** The analytic prediction's `E|Ĥ|` is low, and the `sgdef` tabulation established the gap is **between-rule**: within every covariate signature the realized values match their population values, so it is not sample realization. The search weights *larger* rules than the analytic `sel_c` does, by **+2.11 / +0.61 / +1.65 subjects**. Four mechanisms have now been eliminated — sample realization, the confounder list, `se_g`'s prevalence scaling, and dedup asymmetry (≤ 5%).

**This lead is CC's own**, from `REPORT_residual_2026-08-30.md` §4, and it is the one with a measured magnitude already behind it. `fs_oc_family_enumerate()` runs `get_FSdata()` on `df_super`, so cuts land at **population** quantiles — that is what makes the family fixed and deterministic, and it was the point of the design. The search runs `get_FSdata()` on **each replicate**, so cuts land at that replicate's quantiles. `REPORT_fs_family_report_2026-08-30.md` measured the size of that difference on real data: **15 of 30 continuous thresholds differ between the ACTG frame and `df_super`.**

**Chat does not predict the direction, and neither should this task.** A sample-quantile cut selects a fixed *fraction* of its own sample; a population-quantile cut selects a varying one. Which way that pushes the argmax through K-way intersections of correlated covariates is not obvious, and chat has been wrong on this residual twice. **Measure it.**

**Larry's decision, 2026-08-30: this is the last attempt.** If it does not account for the gap, the residual is **closed as a documented limitation** — leads 2–5 are not to be picked up, and no lead 6 is to be proposed.

---

## ⚠ CATEGORY

**No `R/` change of any kind.** No new export, no default change, no edit to any package file, driver or document. Writes: scratch scripts, their `.rds`/`.log`, and the report, under `dev/glm-continuous-sims/`. Plus this task document.

**Compute:** §2 is minutes. §3 is the expensive part — Monte Carlo over replicates, then `fs_oc_predict()` at three cells. **Estimate 2–4 hours.** §2's gate can end the task before §3 is entered. No simulation study, no renders.

**Unattended.** Gates stop the task; never ask, never work around.

---

## 1. Provenance and first commit — GATE

Standard block plus the installed version (**0.3.0**, or **0.3.1** if the dispatch task has landed first — state which, and note that neither changes results here). *GATE:* branch `feature/glm-extension`, clean tree. Copy this document into `dev/tasks/` and commit alone; the `~/Downloads` stem arrives with hyphens stripped, so report both names.

---

## 2. Stage 1 — does the cut placement move membership at all? — GATE

Cheap, and it decides whether §3 is worth running.

Draw **R replicates** from the DGM at each cell's `n` (R = 200 is enough; state what you used and the seed — **use the corrected run's 20260825 family of seeds**). For each replicate, build the search's cut matrix with `get_FSdata()` and the drivers' `forestsearch_args`, exactly as `forestsearch()` would.

Then, for each candidate in the **corrected analytic family** (`M = 1696 / 1890 / 1696`), identified by its clause specification — variable, direction, quantile index, not its threshold value — report:

| quantity | what it answers |
|---|---|
| label correspondence | does the replicate family contain the same *labels*? It should, nearly perfectly — the quantile grid is specified, not learned. Report any label that appears or disappears, and why. |
| `Δ prevalence` = replicate sample prevalence − population `Pg` | the membership shift, per candidate per replicate |
| `Δ prevalence` stratified by **K** and by `Pg` | whether the shift is concentrated in the large or small candidates |
| `Δ prevalence` weighted by **`sel_c`** | **the one that matters** — whether the shift falls on the candidates that actually win |
| `E[Δ prevalence] · n`, `sel_c`-weighted, summed | the shift expressed in **subjects**, directly comparable to +2.11 / +0.61 / +1.65 |

Report the sign and magnitude plainly, per cell.

*GATE:* if the `sel_c`-weighted shift in subjects is **small relative to the gap** (say under 20% of it) **and** shows no systematic dependence on K or `Pg`, the mechanism cannot produce the observed gap. **Report that, skip §3, and go straight to §4.** That is a complete and useful answer.

If the shift is comparable to the gap, or is large with the wrong sign, continue to §3 — a wrong-signed effect of the right magnitude is as informative as a right-signed one.

---

## 3. Stage 2 — the sample-quantile family, and what it predicts

Build a **replicate-averaged family**: the same nine fields, but with each quantity estimated by Monte Carlo over the R replicates' own cut placements rather than computed once on `df_super`.

- `Pg`, `PQg`, `sens_g`, `spec_g` — replicate means of the realized quantities;
- `ovl` — the replicate mean of the pairwise overlap matrix. Accumulate a running mean of one M×M matrix; do **not** store R of them;
- `beta_g`, `se_g` — from `fs_dgm_scale()`, consistent with how the corrected family computes them. State exactly what you did, because this is the one field where "average the family" is least well defined;
- `M` — unchanged unless label correspondence in §2 said otherwise.

Re-run `fs_oc_predict()` on this family at each cell — n ∈ {500, 700} alternative and the n = 500 null, **both gates, same seeds as the corrected run (20260825)** — and report against the corrected run and against measurement:

- the **between-rule size gap**, against +2.11 / +0.61 / +1.65;
- `E|Ĥ|`, PPV, sensitivity, naive bias;
- declaration rate and the inverted `c1`, so it is visible whether anything *else* moved;
- whether the direction and magnitude move toward measurement, away, or not at all.

**State the approximation, prominently, in these terms.** The analytic framework requires a **fixed** family; a replicate-varying family is outside it. Averaging the nine fields over replicates tests whether the *mean* sample-quantile family closes the gap — a first-order test of the mechanism, **not a proof and not a proposed constructor.** Say where the averaging is least defensible (the overlap matrix and `se_g` are the candidates) and what that could do to the answer. **Do not present the averaged family as a fix, and do not propose adopting it.**

---

## 4. Verdict

State, without hedging:

1. Stage 1's `sel_c`-weighted shift in subjects, per cell, against +2.11 / +0.61 / +1.65 — with its sign.
2. If §3 ran: how much of the gap the replicate-averaged family accounts for.
3. **Whether this explains the residual.** "It does not, and the residual is unexplained" is a valid and expected outcome. **Do not manufacture a cause, and do not propose a further mechanism** — Larry has closed the question either way, and a lead 6 is explicitly out of scope.

No recommendation.

---

## 5. Report

`REPORT_residual_quantiles_2026-08-30.md`: provenance · §2's tables, with the gate's arithmetic shown · §3 if it ran, with the approximation caveat stated before the numbers, not after · the three verdicts · ten-line summary. Commit scripts, outputs and report by explicit path. **No push.** No `devtools::install()`.

---

## 6. Out of scope

- No `R/` change, no new constructor, no edit to any driver, application or document. `oc_wrapper_verification.qmd` is not re-pointed here.
- **No further residual mechanism** — leads 2–5 from `REPORT_residual_2026-08-30.md` §4 are closed by decision, and no new one is to be proposed.
- No breadth testing on a second DGM, no binary/OR work.
- No re-derivation of the `sgdef` tabulation or the corrected enumeration — read the stored `.rds`.
- No simulations beyond §2's replicate draws, no renders, no push.
