# CC task — breadth, stage 2: one measured cell at the forecast harm, scored against the locked forecast

**File:** `dev/tasks/cc_task_oc_breadth_stage2_2026-08-31.md` · **Issued:** 2026-08-30 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1).
**Reads:** `REPORT_oc_breadth_ladder_2026-08-30.md` — **the locked forecast, commit `c9cb0ca2`** — and its `.rds` (the forecast rung's family, grid table and `fs_oc_predict` objects); the MD40 alternative n = 500 driver and its tracked payload; the MD40 null n = 500 payload; `REPORT_oc_wrapper_confs_2026-08-30.md` (how measured columns are computed from payloads) and `REPORT_oc_wrapper_sgdef_2026-08-29.md` (how realized rules are evaluated on `df_super`).

**Larry's decision, 2026-08-30 — go on Stage 2.** Stage 1 issued, before any replicate was drawn, a forecast at the rung `q = 120` (oriented Q effect 120 against MD40's 40): at `c1* = 135.741` the search declares with power **0.800 ± 0.001**, the null false-declares at **0.0387 ± 0.0004**, and the declared rule is a fragment of Q — `E|Ĥ|` 72.8, PPV 0.775, sensitivity 0.333, specificity 0.953, `E[β(Ĥ)]` 98.9 against a Q effect of 120, naive bias 62.3 — with the split gate agreeing to ≤ 0.001 everywhere. This task measures that cell and scores the forecast. **The forecast is not re-derived, re-read from memory, or adjusted; every predicted number is read from the locked artefacts.**

**The scoring is at thresholds the search was not run at**, which is possible because `c1` enters the wrapper only as `T ≥ c1` on the winner's statistic. §2 establishes from source whether the same is true of the search; §3 runs the direct case as well so the equivalence is a gate, not an assumption.

---

## ⚠ CATEGORY

**No `R/` change of any kind.** No new export, no default change, no edit to any package file, existing driver, application or document. Writes: **one new driver** (a copy of the MD40 alternative n = 500 driver with the DGM's interaction and the output stem changed, under the MD40 driver's own directory), its payload(s), scratch scoring scripts and their `.rds`, and the report under `dev/glm-continuous-sims/`. Plus this task document.

**Compute — the go/no-go, stated.** Two runs of the new driver at the MD40 run's replicate count (expect 1,000), possibly a third of the null driver (§4). The dispatch report measured ≈ 5 s per `forestsearch()` call at this fixture under the driver's `parallel_args`, so ≈ 1.4 h **sequential** per run; the MD40 driver's outer parallelism divides that by its worker count. **Before launching, state the projection from the MD40 run's own recorded wall-clock and the driver's worker setting. *GATE:* if the projection for all runs exceeds 6 hours, stop and report instead of launching.** No renders.

**Unattended.** Gates stop the task; never ask, never work around.

---

## 1. Provenance and first commit — GATE

Standard block (`hostname; pwd; git rev-parse --abbrev-ref HEAD; git rev-parse --short HEAD; git status --porcelain; git log --oneline -6`) plus the installed version (expect **0.3.1**). *GATE:* branch `feature/glm-extension`, clean tree, `8ff85f9c` in the log, and `git show --stat c9cb0ca2` names `REPORT_oc_breadth_ladder_2026-08-30.md`. Copy this document into `dev/tasks/` and commit it alone; the `~/Downloads` stem arrives with hyphens stripped, so report both names.

## 2. Source gates — before anything runs

**(a) The forecast, read not typed.** From the Stage 1 `.rds` files, read to full precision: `k_inter` at the forecast rung (the ladder's `k40 + s·(120 − 40)`), `c1*`, `c1_05`, and every predicted quantity of the forecast table under both gates, plus the forecast rung's grid `table` (`det_rate` against `c1`) and the null grid or inversions used for the null curve. Print them with 6 decimals in the report's §0 as the values being scored.

**(b) The MD40 driver.** Quote from the MD40 alternative n = 500 driver: how the DGM is built (`calibrate_glm_interaction()` — the same call Stage 1 §2 reproduced); how per-replicate seeds are assigned (the seed table and its base); the replicate count; the `forestsearch_args` (expect `effect.threshold = 30`, `consistency.threshold = 10`, the 13-confounder list, `parallel_args = list(plan = "sequential")`); the outer parallel plan and worker count; and **exactly what the recorder stores per replicate** — the winner's fitted effect, its definition (`sg_def`), its sample size (`n_harm`), its membership or enough to reconstruct it, `p.consistency`, and the OC block. *GATE:* the per-replicate **winner's fitted effect** and **rule definition** must be stored, since §4 scores from them. If either is missing, stop and report — the recorder would need extending, which is not authorised here.

**(c) Where `effect.threshold` enters the search.** Establish from `R/forestsearch_main.R` and the files it calls, quoting each site, every place `effect.threshold` (or the value derived from it) is used. The expectation, from the wrapper's own reduction, is that it acts **only as a filter on the candidate table before the consistency loop** — dropping candidates whose fitted effect is below it — and nowhere in the consistency evaluation, the near-duplicate key, the early-stopping rule, the size floors, or the GRF screen. Under that expectation, raising the threshold from 30 to `c1*` removes only candidates below the winner, so the winner at `c1*` is the winner at 30 whenever its fitted effect is ≥ `c1*`, and post-hoc scoring is exact. **State whether the expectation holds.** If it does not, say where it fails; §3's direct run then becomes the measured cell and §4's post-hoc scoring is reported as a cross-check only — carry on either way.

**(d) The null payload.** Confirm the MD40 null n = 500 payload stores the per-replicate winner's fitted effect. If it does, §4 scores the null post hoc and **no null run is made.** If it does not, a third run — the null driver at its own seeds, unchanged except the output stem — is **authorised by this task** as part of the compute go/no-go; say so in the report.

## 3. The measured cell — two runs

**The new driver.** Copy the MD40 alternative n = 500 driver to a sibling file whose stem replaces `md40` (or the equivalent tag) with `md120`, and change exactly three things: the DGM construction, the output stem/scenario tag, and — for run 2 only — `effect.threshold`. Everything else, including the seed table, the replicate count, `forestsearch_args`, the recorder and the outer parallel plan, is byte-identical; show the diff against the MD40 driver in the report.

**The DGM.** Build with `generate_glm_dgm()` by the exact linear route Stage 1 used for the rung — every argument identical to the MD40 build's direct form except `k_inter` at the value read in §2(a). *GATE before launch:* `fs_dgm_scale()` gives `abs(m_tau[Q]) = 120` to 1e-9 and `abs(m_tau[Qc])` equal to MD40's; the family enumerated from this DGM at n = 500 is `identical()` to the stored forecast-rung family on `lab`, `Pg`, `PQg`, `beta_g`, `se_g`, `sens_g`, `spec_g`, `M`.

**Run 1 — the comparability run.** `effect.threshold = 30`, the driver's own. This is the run scored post hoc at every `c1`, and the one comparable to the MD40 measured columns.

**Run 2 — the direct run.** `effect.threshold = c1*` at full precision, everything else identical to run 1, **same per-replicate seeds.**

Launch detached and monitored, sequentially or together as the box allows; record each run's wall-clock. If a run dies, report and stop.

**Paired replicates.** Because the seed table is the MD40 run's, each replicate here samples the same `df_super` rows and treatment assignment as the corresponding MD40 replicate; state whether the residual draws are also aligned (whether the RNG stream between sampling and outcome generation is consumed identically), since if they are, the two cells differ per replicate only by the shift of Q's treated outcomes.

## 4. Scoring — computed from the payloads, never typed

Reuse the corrected run's and the `sg_def` report's machinery for evaluating a realized rule on `df_super` (population size, purity, sensitivity, specificity, true effect). For **run 1**, per replicate: `declared_30`, the winner's fitted effect `T̂`, its definition, `n_harm`, and its population functionals.

1. **Power at `c1*`.** `declared_c1* = declared_30 & (T̂ ≥ c1*)`; the rate with its binomial SE — **against the forecast's 0.800 ± 0.001.**
2. **The measured declaration curve.** The rate of `declared_30 & (T̂ ≥ c1)` for `c1` from 0 to 300 by 1, beside the forecast rung's predicted `det_rate` at the same `c1` — the figure. Report the largest absolute gap and where it sits.
3. **The direct run as the gate on post-hoc scoring.** From run 2: its declaration rate with SE; and per replicate, whether run 2 declares iff `declared_c1*` from run 1, and whether the winner's definition and membership are identical. *GATE:* if §2(c) said the expectation holds, agreement must be complete — report any disagreement with the replicate ids and the two winners, and if there is any, treat run 2 as the measured cell. If §2(c) said it does not hold, run 2 is the measured cell and this comparison is the record of how far post-hoc scoring is off.
4. **Type-I at `c1*`.** From the null payload (or the authorised null run): `declared_30 & (T̂_null ≥ c1*)`, the rate with its binomial SE — **against 0.0387 ± 0.0004**; and the null's measured curve against `c1` beside the null grid's predicted curve.
5. **The declared rule, conditional on `declared_c1*`:** mean `n_harm` (sample) and mean population size of the realized rules; sample and population PPV, sensitivity, specificity; `E[β(Ĥ)]` as the population true effect of the realized rule; naive bias as `T̂ − β(Ĥ)`. Each with its SE, beside the forecast row under both gates. **Carry the limitation verbatim:** the population-versus-sample offset on PPV and sensitivity is expected (the driver's sample proportions at a selected rule are biased upward), and at MD40 `E|Ĥ|`, PPV, sensitivity and the naive bias inherit a between-rule gap of +2.11 subjects that five mechanisms failed to explain. Compute the between-rule gap **at this harm** the way the `sg_def` report did — measured-frequency-weighted population size against analytic `sel_c`-weighted population size — and report it beside +2.11. Whether it grew, shrank or held is a finding; **no mechanism is to be proposed.**
6. **Composition.** Measured selection mass on rules with `PQg ≥ 0.95` beside the forecast's 0.345; the top 15 realized rules by frequency beside the top 15 analytic by `sel_c`, with `Pg`, `PQg`, `β_g`; total variation distance beside its multinomial noise floor, as the `sg_def` report did. Note that Q itself is not a family member (Stage 1) and report the realized rules' closest analytic labels.
7. **At the driver's `c1 = 30`** — the saturated picture at this harm: measured declaration rate, `E|Ĥ|`, PPV, sensitivity, `E[β(Ĥ)]`, naive bias beside the forecast rung's row at `c1 = 30`, and beside MD40's measured columns from the corrected report.

## 5. Verdict

State, without hedging, for each scored prediction: **within noise** (measured − predicted inside 2 combined SEs), **population-versus-sample** (the known offset, for PPV and sensitivity), or **a gap** (with its size in the quantity's units):

1. Power at `c1*` — measured vs 0.800.
2. Type-I at `c1*` — measured vs 0.0387.
3. `E[β(Ĥ)]` and specificity at `c1*`.
4. PPV, sensitivity, `E|Ĥ|`, naive bias at `c1*`, and the between-rule gap at this harm against +2.11.
5. The declaration curves, alternative and null, over `c1`.

Then one paragraph: does the wrapper predict the crossover regime — a threshold with 80% power and small type-I — as well as it predicted the saturated one? Whatever the answer, **no recommendation, no mechanism, no next task.**

## 6. Report

`dev/glm-continuous-sims/REPORT_oc_breadth_stage2_2026-08-31.md`: §0 the forecast values as read (6 decimals) and the commit they were read from · §2's quotations and the three verdicts on source · §3's driver diff, DGM gates, timings and the pairing statement · §4's tables in order · §5 · ten-line summary.

Commit the driver, its payloads under the MD40 driver's tracking convention (if the MD40 payload is tracked, track these; if a single file exceeds 50 MB, leave it untracked and say so), the scoring scripts, their `.rds`, and the report, by explicit path. **No push. No `devtools::install()`.**

## 7. Out of scope

- No `R/` change, no recorder extension, no edit to the MD40 drivers, the OC-wrapper files, or any document. `oc_wrapper_verification.qmd` is not touched.
- No `n = 700`, no other rung, no further null run beyond §2(d)'s conditional authorisation.
- No change to `se_g` — Stage 1's direct-`V_eff` sensitivity stays a sensitivity and is not re-run here.
- No residual mechanism; the between-rule gap at this harm is recorded, not explained.
- No re-derivation of the forecast, no push, no install, no renders.
