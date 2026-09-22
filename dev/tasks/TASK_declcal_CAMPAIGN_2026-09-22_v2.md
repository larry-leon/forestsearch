# TASK v2 — Finite-sample evaluation of the calibrated declaration screen

**Filename note.** This file was previously issued as
`TASK_declaration_calibration_evaluation_2026-09-22_v2.md` and never reached the machine — its name was too
close to the implementation task's and the wrong file was downloaded. **The content is unchanged**; only the
filename and the §0 destination path differ. There is no third version of this document.

**Supersedes** `TASK_declaration_calibration_evaluation_2026-09-21.md` (v1) in full. Do not execute v1.
Three substantive changes, all forced by findings in `REPORT_declaration_calibration_2026-09-22.md`:

- **§2.2 is new and mandatory.** The executed screen compares `round(rate_closed, pconsistency.digits) >= p_star`
  with `pconsistency.digits = 2`, so its effective cutoff is `z_{(1+0.895)/2} = 1.6211`, not `1.6449`.
  v1's fidelity gate compared the search's indicator against `T >= 1.6449` and would have **stopped
  spuriously** on any replicate with a candidate whose rate falls in `[0.895, 0.900)`. Over 18,000
  replicates on families of order 1,800 that band will be hit. Fixed here.
- **§2.3 is new.** The field is computed on **every** replicate, not only on declaring ones, and
  `keep_field_matrix` is `FALSE`. v1 left both implicit and the pricing assumed `nullmr`'s
  declaring-replicates-only pattern, which does not support the §7.3 diagnostics.
- **§9 carries measured pricing.** The costing appendix supplied real numbers; the ceiling arithmetic is
  now stated from them rather than deferred entirely to the pilot.

Everything else is reproduced unchanged so this file is the single authoritative set.

**Repo:** `larry-leon/forestsearch`
**Branch:** `feature/glm-extension`
**Machine:** Pop!_OS / System76, 64 workers (match the `nullc125` / `nullmr` configuration so the measured
cost transfers)
**Prerequisite:** commits `8a318c10`, `c5e9f91e`, `b0507556` present, and the package installed with
`devtools::install()` — parallel workers see only the installed package, never `load_all()`.
**Kind:** simulation campaign. **This is compute.** Nothing runs until Larry's go, and the go carries the
wall-clock ceiling of §9.
**R/ call-out:** this task makes **no change under `R/`**. It calls the exported
`fs_declaration_calibration()` and the `fs_mr_inference()` capture formals added by the prerequisite task.
If any `R/` change turns out to be needed, that is a STOP and a report, not an edit.

---

## 0. First action

1. Copy this file verbatim to `dev/tasks/TASK_declcal_CAMPAIGN_2026-09-22_v2.md`
2. `git add` that explicit path and commit:
   `docs(tasks): add declaration-calibration finite-sample evaluation task v2 (2026-09-22)`
3. Record `git rev-parse HEAD` and the installed package version into the report.
4. Assert against the **installed** package:
   - `"fs_declaration_calibration" %in% getNamespaceExports("forestsearch")`
   - `keep_declaration_field` and `keep_field_matrix` are formals of `fs_mr_inference` with default `FALSE`
   If either fails, STOP: the prerequisite is not installed.

No `git fetch`, `git pull`, or `git push`. Every `git add` names explicit paths. Never stage a pre-existing
untracked file.

---

## 1. The question this campaign answers

Whether, in finite samples, the calibrated screen holds the null declaration rate **near alpha** where the
theory says it should (boundary null) and **at or below alpha** elsewhere (interior null) — and what it
costs in power.

"Declaration" is the family-wise event: **at least one candidate is admitted**. It is not the accuracy of
the declared subgroup, and classification metrics are not part of this campaign.

Context that raises the stakes: on the GBSG fit of the prerequisite report the `p* = 0.90` screen had
family-wise size **0.94** over 765 candidates, and the calibrated rule at `alpha = 0.05` admitted none.
That is one application at one threshold setting. This campaign is what decides whether the calibrated
screen is a usable operating rule or only a diagnostic — so it must not be reported as either until it runs.

---

## 2. The two screens, evaluated on the same run

Restated self-contained; do not consult the KB or the manuscript for these.

- `T(g) = { beta_hat(g) - c_cons } / sigma_D(g)`, with `c_cons` the consistency threshold (c2) and
  `sigma_D(g)^2 = sum_{i in g} db[g,i]^2` the robust variance the screen already uses.
- **Calibrated screen [Eq. 7]:** declare if `max_g T(g) >= kappa_hat(alpha)`, with
  `kappa_hat(alpha)` the empirical `(1 - alpha)` quantile of `Mstar[b] = max_g Zstar[b,g]`,
  `Zstar[b,g] = (sum_{i in g} xi[b,i] * db[g,i]) / sigma_D(g)`. Implied rate
  `pstar_implied = 2 * pnorm(kappa_hat(alpha)) - 1`. **No rounding** — this is a new rule with no digits
  convention.
- **Diagnostic [Eq. 8]:** `alpha_FW_hat = mean( Mstar[b] > cutoff )`, computed at **both** cutoffs of §2.2.

### 2.1 The cost design — one search per replicate (mandatory)

Both screens differ only in the constant on the right of an inequality the run already evaluates. So:

- Run the forest search **once** per replicate.
- Capture the pre-reduction family with per-candidate `beta_hat(g)`, `sigma_D(g)`, and `Mstar`, via the
  `fs_mr_inference()` capture and `fs_declaration_calibration()`.
- Evaluate **both** screens post hoc from the same `T(g)` vector and the same multiplier draws.

**Do not run a second search, a second resampling pass, or a second campaign arm per screen.** Any design
that re-runs the search per screen doubles the campaign for nothing and is a STOP.

### 2.2 The conventional arm must be evaluated the way the screen executes it

`R/subgroup_consistency_helpers.R:1468` compares `round(rate_closed, pconsistency.digits) >= p_star`, and
`pconsistency.digits` defaults to 2 and is not passed by `forestsearch()`. Consequences, all mandatory:

- **`declared_conv` is computed by the rounded-rate rule**, `round(max(0, 2*pnorm(T) - 1), digits) >= p_star`,
  with `digits` read from the fit — **not** by `T >= 1.6449`.
- **`declared_conv_exact`** is recorded alongside, using `T >= qnorm((1+p_star)/2) = 1.644854`, so the
  rounding's reach is measurable rather than assumed.
- **`n_band`** per replicate: the count of candidates whose closed-form rate lies in `[0.895, 0.900)` —
  the band where the two disagree. Report its distribution.
- **`alpha_FW_hat` is computed at both cutoffs:** `1.644854` (nominal) and `1.621` (the effective cutoff
  of the rounded rule, `qnorm((1+0.895)/2)`). The nominal one **understates** the executed screen's true
  family-wise size; reporting both quantifies by how much.
- **Fidelity gate (this is what the cost design rests on).** `declared_conv` — the *rounded* rule, on the
  **post-reduction** family — must equal the search's own declaration indicator, replicate by replicate, in
  every cell. Assert exact agreement. Disagreement is a STOP with the disagreeing replicate indices
  recorded. Do **not** relax the gate to `declared_conv_exact`; that is the comparison v1 got wrong.

### 2.3 The field runs on every replicate

- Call with `keep_declaration_field = TRUE` and **`keep_field_matrix = FALSE`**. Storing a `B x ~1800`
  matrix for 18,000 replicates is not feasible; `Mstar` plus the per-candidate summaries is sufficient for
  everything §5 records.
- The field is computed on **every** replicate, including those where no subgroup was declared. Two
  independent reasons:
  - §7.3 needs the implied-`p*` distribution and `mean(alpha_FW_hat)` **across replicates**, which do not
    exist on replicates that were skipped;
  - `kappa_hat` can fall below `1.6449` on a small or strongly correlated family, in which case the
    calibrated rule can declare where the conventional one did not — so restricting to declaring replicates
    would bias the calibrated rate downward by construction.
- `forestsearch()` runs MR only when a subgroup is identified. If the field therefore cannot be captured on
  a non-declaring replicate through `forestsearch()`, call `fs_mr_inference()` directly on that replicate's
  enumerated family. If **neither** route reaches a non-declaring replicate without an `R/` change, that is
  a **STOP** — report it; do not substitute a declaring-replicates-only design.

---

## 3. Family convention

- The maximum for `kappa_hat` ranges over the family **prior to any fitted-summary reduction**
  (`family = "prereduction"`, the settled ruling of handoff v2 §5). Primary convention for every reported rate.
- The post-reduction family is recorded in parallel as a secondary quantity so the gap is measurable. On the
  GBSG fit that gap was exactly zero (765 vs 764); do not assume it stays zero at other n.
- **Record `G_pre` and `G_post` every replicate.** The costing read puts `G_pre` around 1,711–1,830 on this
  grid — an order of magnitude above the GBSG fit's 765, so `kappa_hat` is expected higher here.
- The pre-reduction family deliberately excludes only outcome-dependent filters. The per-arm event minima
  (`d0.min`/`d1.min`) are **not** applied, because event counts are outcome data and a family filtered on
  them would not be covariate-measurable — which is exactly what the theorem forbids. The size minimum
  `n.min` **is** applied and is covariate-measurable. Record the floors in force in every cell's payload and
  assert they equal the `gbsg_020` template's values.

---

## 4. Cells

DGM: the `gbsg_020` template. Replicates: **2,000 per cell**. **Forest search only — one identifier, not
three.** Identification plus the standardized field only: **no MR post-selection correction, no intervals,
no bounds** beyond the field capture.

### Block A — boundary null (the "close to alpha" test) — 3 cells

| cell | DGM | c1 | c2 | n |
|---|---|---|---|---|
| A1 | complete null, HR 1 in every subgroup | 1.0 | 1.0 | 500 |
| A2 | complete null, HR 1 in every subgroup | 1.0 | 1.0 | 1000 |
| A3 | complete null, HR 1 in every subgroup | 1.0 | 1.0 | 1500 |

`c1 = c2` so the relevance floor does not bind and the null sits exactly at `c_cons`. Theorem 2 predicts a
calibrated rate near alpha. At 2,000 replicates the 95% Wilson band around 0.05 is about 0.040–0.060.
Convergence with n is expected; the rate may run **above** alpha at n = 500, where the smallest candidates
are least Gaussian — §5 records what is needed to diagnose that rather than just observe it.

### Block B — interior null (validity) — 6 cells

The six `nullid` designs, re-run at `c1 = c2 = 1.0` with the calibration on:

| cell | DGM | c1 | c2 | n |
|---|---|---|---|---|
| B1–B3 | uniform benefit HR 0.657, no planted region | 1.0 | 1.0 | 500 / 1000 / 1500 |
| B4–B6 | uniform benefit HR 0.721, no planted region | 1.0 | 1.0 | 500 / 1000 / 1500 |

Prediction: at or below alpha, falling with n.

### Block C — power cost (optional; switched on in the kickoff) — 4 cells

| cell | DGM | c1 | c2 | n |
|---|---|---|---|---|
| C1–C2 | planted harm region at HR 1.5 | 1.0 | 1.0 | 1000 / 1500 |
| C3–C4 | planted harm region at HR 2.0 | 1.0 | 1.0 | 1000 / 1500 |

**Design-feasibility firewall, pre-flight, before any Block C replicate runs.** Compute the planted
region's expected prevalence under the template and check it against the floor actually in force
(`min(60, 0.10 * n)` where that is the rule). A cell whose planted region falls below the floor in a
material share of replicates is **not run**: record the arithmetic and the exclusion. This exclusion is the
firewall working, not a failure — the rest of the campaign proceeds. n = 500 is deliberately absent; do not
add cells, the cell count is fixed by this document.

`c1 > c2` settings (for example `nullc125`'s 1.25 / 1.00) are conservative by construction and are **not** a
test of "close to alpha". Not in this campaign.

**GRF and DINA are excluded.** They are conditional on their proposed families, so a family-wise
declaration rate over an enumerated family is not the same quantity for them. Record in one line; not a
follow-up task.

---

## 5. Per-replicate record schema

One row per replicate. Columns, exactly:

```
cell_id, dgm, n, c1, c2, rep, seed,
G_pre, G_post,
max_T_pre, max_T_post,
kappa_hat_05, kappa_hat_10,
pstar_implied_05, pstar_implied_10,
alpha_FW_hat_1645, alpha_FW_hat_1621,
Mstar_q90, Mstar_q95, Mstar_q99,
declared_conv, declared_conv_exact, declared_cal05, declared_cal10,
n_admitted_conv, n_admitted_cal05, n_band,
sg_size_declared, sg_size_argmax,
B_cal, mult_law, pconsistency_digits, floors_id, wall_sec, status
```

Notes on four of these, because they are what make the payload reusable:

- `max_T_pre` / `max_T_post` mean any conventional `p*` can be re-evaluated later without re-running.
- `Mstar_q90/95/99` mean any alpha in that range can be re-evaluated later without re-running.
- `n_band` quantifies the §2.2 rounding band's reach.
- `sg_size_declared` and `sg_size_argmax` connect an over-rate at n = 500 to small, least-Gaussian
  candidates. Record them even where nothing was declared (NA).

`status` is `"ok"`, `"abort_time"`, or `"error"` with the condition message stored alongside.

---

## 6. Seeds and draw parity

- One documented seed stream per `(cell_id, rep)`; a single seed-offset convention recorded in the payload
  metadata. The earlier r-hat episode was partly a seed-grid artifact under a changed offset — do not
  introduce a second offset convention here.
- Both screens read the same replicate data and the same multiplier draws. Automatic under the §2.1
  post-hoc design; assert it by construction, with no re-draw anywhere.
- `B_cal` is fixed per campaign and recorded. §8 measures the sensitivity and selects it by rule.
- Multiplier law: the production law (centred Poisson), the same in every cell. The prerequisite report
  measured the Gaussian-vs-Poisson discrepancy at n = 400 as indistinguishable from Monte-Carlo error, so
  the production law is the right choice here and no Gaussian arm is needed.

---

## 7. Reporting

Write `REPORT_declaration_calibration_evaluation_2026-09-22.md` in `dev/reports/`, matching where the
prerequisite report was placed.

### 7.1 Primary table — cells × screen, proportions with Wilson 95% intervals

One table, this shape, no prose between rows:

| cell | DGM | n | conventional (as executed) | conventional (exact z) | calibrated alpha = 0.05 | calibrated alpha = 0.10 |
|---|---|---|---|---|---|---|
| A1 | complete null | 500 | p̂ [L, U] | p̂ [L, U] | p̂ [L, U] | p̂ [L, U] |
| … | | | | | | |

Wilson intervals, not Wald. Follow the table with a **short plain-language reading** — a few bullets, one
claim each — not a narrative.

### 7.2 Second table — the calibration's own quantities

| cell | G_pre (median, IQR) | G_post (median, IQR) | kappa_hat_05 (median, IQR) | implied p\* (median, IQR) | mean alpha_FW_hat at 1.6449 | mean alpha_FW_hat at 1.621 | median n_band |
|---|---|---|---|---|---|---|---|

### 7.3 The three free checks

- **Eq. 8 as an estimator.** At each Block A cell, `mean(alpha_FW_hat_1621)` against the **as-executed**
  conventional declaration rate on the same replicates, and `mean(alpha_FW_hat_1645)` against the exact-z
  rate. Report both pairings, their differences, and a paired Monte-Carlo standard error. A material gap is
  a **finding to record**, not a gate.
- **How strict the calibration is.** Per cell, the distribution of `pstar_implied_05`: min, 5%, 25%, 50%,
  75%, 95%, max, plus the fraction of replicates with `pstar_implied_05 > 0.90`.
- **Whether the calibrated rule can ever declare where the conventional one did not.** The count of
  replicates with `declared_cal05 == 1 & declared_conv == 0`, and the same at alpha = 0.10. Expected zero or
  near-zero on families of this size; a non-zero count is the case §2.3's second reason anticipates and must
  be reported, not discarded.

### 7.4 Also in the report

- The pilot's measured numbers (§8) and the ceiling arithmetic against §9's projection.
- Block C's pre-flight prevalence arithmetic and any cell excluded by the firewall.
- Every gate with its measured value, not just pass/fail.
- An **OPEN ITEMS** block for documentation gaps. Documentation gaps do not stop the run; untrustworthy
  results do.

### 7.5 Closeout (last action)

Regenerate `quarto/simulations/gbsg_020/current_status.md` to catalogue the new campaigns, their payload
locations, reading conventions and open work, stating the commit it describes — with the machine-checkable
post-condition that the stated pin equals HEAD at commit time. Every file carrying payload or summary
content is tracked; nothing another repo would need is left in a session scratchpad.

---

## 8. Stage 1 — pilot, and the rule that continues automatically

The pilot prices the campaign from a measured clock. It is not a question.

- Cell A2 (complete null, n = 1000), **200 replicates**, at `B_cal` in `{500, 1000, 2000}` on the same 200
  replicates (field assembled once at 2000 and sub-sampled — do not re-draw).
- Measure and record: per-replicate wall-clock split into search and field capture; worker count;
  `G_pre` / `G_post` distributions; the declaration rate at each `B_cal`; and the realized fraction of
  replicates that declared conventionally (the quantity §9's projection depends on).
- **`B_cal` selection rule (no chat round trip):** the smallest `B_cal` whose calibrated declaration rate is
  within 0.005 of the `B_cal = 2000` value; if none is, use 2000. Record the three rates and the choice.
- **Projection:** measured median per-replicate wall-clock × the replicate count of the approved blocks
  (9 cells × 2,000 = 18,000 core; 13 × 2,000 = 26,000 with Block C) at the pilot's worker count.

**Then, without stopping:** if the projection is at or below the kickoff's ceiling, continue straight into
Stage 2. If it exceeds the ceiling, **stop and report** with the projection and its measured inputs. This is
a stop-on-failure gate, not a stop-to-ask gate — a green pilot needs no round trip and the machine does not
idle.

---

## 9. Measured pricing, and Stage 2's abort discipline

From the prerequisite report's costing appendix — measured values only, all on pop-os at 64 workers:

- `nullmr` (this grid, MR on, 6 cells × 3 identifiers × 2,000 = 36,000 replicate-searches):
  **19,858 s**, i.e. **0.552 s per replicate-search**, with MR run on 21,014 of 36,000 declarations (58%).
- `nullc125` (same shape, MR off): 4,002 s; FS runs alone 1,184 s.

Projection for this campaign, stated with its assumptions so the pilot can correct it:

- Core (18,000 replicates): `18,000 × 0.552 = 9,936 s = 2.8 h`.
- The field must run on **every** replicate (§2.3), not 58% — scale the MR share by about `1/0.58 = 1.72`:
  **roughly 4.7 h core, 6.9 h with Block C.**
- Two things push this up and the pilot, not this document, settles them: `G_pre` here is ~1,711–1,830
  against the 765 the 0.552 s/replicate figure did not have to cover at every replicate, and the §7.3
  diagnostics add per-replicate work.

**A 12-hour ceiling leaves roughly 1.7× headroom on the core blocks and 1.7× with Block C.** That is the
recommended ceiling; the kickoff sets the number.

Abort discipline:

- **Per-replicate hard cap:** 10 × the pilot's median per-replicate wall-clock. On overrun, abort that
  replicate, record `status = "abort_time"`, continue.
- **Per-cell gate:** if more than 1% of a cell's replicates end `abort_time` or `error`, stop that cell,
  record it, continue with the remaining cells. Rates are never reported from a cell that tripped this.
- **Campaign hard cap:** 1.5 × the §8 projection. On overrun, stop, write the report from whatever cells
  completed, and mark the rest incomplete.
- Checkpoint each cell's payload to disk as it completes, so a stop never loses completed cells.
- If a live campaign belonging to another workstream is running on the machine, record it and do not compete
  for cores: report and stop.

---

## 10. Build posture

- **Transplant-first, author-never.** Copy the committed `nullid` campaign script in
  `quarto/simulations/gbsg_020/` and change named lines: the DGM effect, `c1`/`c2`, `n`, the identifier set
  (FS only), the capture formals, and the record schema of §5. Do not author a parallel campaign script from
  scratch. Record which committed file was transplanted and the lines changed.
- Campaign names: `declcal_bnull` (Block A), `declcal_inull` (Block B), `declcal_power` (Block C), under
  `quarto/simulations/gbsg_020/`, alongside `nullid`, `nullc125` and `nullmr`.
- tidyverse style. No new hard dependency. Installed package only.
- **No `R CMD check`, no `rcmdcheck`, no vignette build, no full testthat suite.**

---

## 11. Commit plan

Explicit named paths on every `git add`. No push — Larry pushes via GitHub Desktop.

1. `docs(tasks): add declaration-calibration finite-sample evaluation task v2 (2026-09-22)` — §0.
2. `feat(sims): add declcal boundary/interior null campaign scripts (transplanted from nullid)` — the
   campaign scripts and the pre-flight feasibility script, before any run.
3. `data(sims): add declcal pilot payload and measured costing` — the pilot payload and its numbers.
4. `data(sims): add declcal campaign payloads` — Stage 2 payloads, per block.
5. `docs(report): record declaration-calibration finite-sample evaluation` — the report.
6. `docs(sims): regenerate gbsg_020 current_status at <SHA>` — the closeout.

---

## 12. Explicitly out of scope

- No `R/` change of any kind.
- No MR post-selection correction, no intervals, no bounds, no field constructions beyond the capture.
- No classification or accuracy metrics; declaration rates only.
- No `c1 > c2` cells, no `nullc125` re-run, no GRF, no DINA.
- No re-run of the identifiers campaign, and no re-verification of anything already committed.
- No change to any threshold default; `c1 = c2 = 1.0` here is a campaign setting, not a package default.
- The `pconsistency.digits` rounding is handled as §2.2 specifies and is **not** fixed, changed, or
  proposed against here. It belongs to the threshold-specification workstream.
