# status_curated — the hand-maintained sections of `current_status.md`

This file is **included verbatim** by `scripts_or/current_status_regen.R`. Its four marked blocks
are placed as follows: `preamble` after the pin, `before-inventory` (§1–§2) before §3,
`after-inventory-table` after §3.5's table, and `after-inventory` (§4–§7) after §3. Keep the
markers, and keep every number that the generator can derive **out** of this file — the generator
reads the directory and git, never a report and never a chat record.

<!-- curated:preamble -->
- **What this directory is.** `quarto/simulations/actg175/binary_020/` holds the ACTG175 binary /
  odds-ratio design: the committed study that produced the supplement's Figures S9 and S10, its
  21 per-cell payloads and coverage grid, and the three campaigns that re-run that design under
  the current constructions — `orfs` (forest search), `orgrf` (GRF) and `ordina` (DINA).
- **Two things live here side by side, and they are not the same study.** The committed study ran
  `maxeffCons` at ε 0.10 with infinitesimal-jackknife intervals only, 1,000 replicates per cell,
  over a seven-point sample-size sweep at the single protective design point OR 0.75. The three
  campaigns run `effMaxSG` at ε 0.20 with the field constructions, 2,000 replicates per cell, at
  two sample sizes and three design points. **Neither supersedes the other**, and no table mixes
  them without saying so.
<!-- curated:before-inventory -->
## 1. The committed study

- **Driver:** `maxeffCons_mr_coverage_sweep_or075.qmd`, the producer of the supplement's Figures
  S9 and S10, with its figure fragment `_sim_mr_coverage_or075.qmd` and its render. **Not edited**
  by any campaign in this directory.
- **Design:** ACTG175 arms 1 (ZDV+ddI) vs 3 (ddI); the week-20 adverse outcome
  `y_neg = 1 - 1{cd420 > cd40}`, analysed directly, so **OR > 1 is harm**. The planted harm region
  is H = {wtkg > q70} ∩ {cd40 > q70}, calibrated to a marginal OR of 0.75 in H against a
  complement that inherits the fitted ACTG175 effect. Thresholds: effect 0.90, consistency 0.80,
  p-consistency 0.90.
- **Payloads:** `mr_sweep/maxeffCons_actg175_or075_seedtab_s1000/` — 21 per-cell bundles
  (3 identifiers × 7 sample sizes) and the coverage grid. Read by the campaigns' summary for the
  study-comparison table and by the Stage 1 data-recipe check; **never rewritten**. The campaigns'
  template refuses to write any git-tracked bundle path.
- **Its constructions:** the unadjusted plug-in, the per-replicate oracle, and the two-term
  de-biased multiplier-resampling estimate with its IJ interval. It carries **no field columns**.

## 2. The three campaigns

- **Rule, in every cell:** `sg_focus = "effMaxSG"`, `effect_neighborhood = 0.20`,
  `selection_rule = "neighborhood"`, with the study's thresholds unchanged and
  `adverse_outcome = TRUE`.
- **Design points,** each built by the study driver's calibration with only `target_effect`
  changed, so the prevalence, the interaction route and the complement references are the study's:
  **protective** (OR 0.75, the supplement's), **harm** (OR 1.5) and **borderline null** (OR 1.0 —
  the planted region *at* the null against a protective complement, as the committed OR 1.0 driver
  has it; the driver's homogeneous `dgm_model = "null"` branch is not used anywhere here).
- **Sizes:** n = 500 and n = 2000, the ends of the study's sweep. Six cells per campaign,
  eighteen in all.
- **Identifiers:** `orfs` runs first and is the reference the other two are checked against —
  the three campaigns share every draw cell for cell, and Gate 2 asserts that the data-level
  columns (the per-replicate seed, the true-region size and the oracle) are identical in both
  directions.
- **Constructions:** the field one-sided lower bound on Ĥ, the field-s one-sided upper bound on
  Ĥᶜ, and their Bonferroni pair, with the unadjusted, oracle and IJ two-term constructions as
  references. Excluded everywhere: the IJ winner-only and winner-floor variants, the κ / uniform
  calibration, covariate adjustment and tuned inflation factors. FB was never run and never
  joined (there is no committed binary FB bundle).
- **One template, one runner.** `sim_fs_mr_field_or_template.qmd` runs every cell of every
  campaign through `FS_OR_*` environment knobs; `scripts_or/run_or.sh` sequences the three
  campaigns over the six cells, two batches and a combine per cell, gating each cell with
  `scripts_or/gate2.R` before it commits.
<!-- curated:after-inventory-table -->
**Reading the table.** Each file is counted once, first matching rule wins, so the rows sum to the
total. `orsmoke` and `orcal*` rows are Stage 1 artefacts and are deleted at closeout, so they are
expected to be absent in a closed-out state. The `mr_sweep/` rows are the committed study's and are
tracked; everything the campaigns write lands under `mr_or_harm/` or in the directory root. Raw
render logs under `logs_or/` are deliberately untracked.
<!-- curated:after-inventory -->
## 4. Reading conventions on this path

- **Everything is an odds ratio, and OR > 1 is harm.** `adverse_outcome = TRUE`, so the analysis
  runs directly on the adverse outcome and no orientation flip is applied anywhere.
- **Two scales coexist and are never mixed** (`R/fs_mr_inference.R:480–488`). The **effect** scale
  is the OR and carries every estimate and every bound. The **working** scale is the log-OR and
  carries the stored Wald and IJ standard errors, the field's `lambda_mean` and `se_field`, and the
  Λ* quantiles. Consequences that matter when reading any table here: `lambda_mean` is a log-scale
  correction and is routinely **negative** — it is not a bound; every bias in SD units and every
  SE-to-SD ratio is computed against SD(log est); and the field identities are
  log(est2) + λ̄ = log(β̃), log(lo1s) = log(β̃) − q95 and log(up1s) = log(β̃ᶜ) − q05.
- **Bounds are read by location against a ladder, never as significance at OR = 1.** The ladder is
  τ = 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 2.0. On Ĥ a one-sided 95% lower bound at or above τ reads
  "harm of at least τ"; on Ĥᶜ a one-sided 95% upper bound at or below τ reads "harm of at most τ".
- **Targets are named on every coverage row.** The conditional rows are scored against β(Ĥ) /
  β(Ĥᶜ), the exact per-replicate population OR at the *realized* region. The two population
  references are θ† (the marginal OR in the block) and θ‡ (the controlled direct effect), recorded
  per replicate as `C_dagger_*` / `C_ddagger_*`; the oracle row is scored against θ†, and both are
  reported beside the harm block, as the committed study did.
- **DINA's proposal floor is on the link scale.** `forestsearch()` derives it as
  `log(effect.threshold)` for every non-Gaussian family, so the floor as applied is the OR-scale
  threshold 0.90 on the harm side, and `dina_tau_min` is a log-OR.

## 5. The conditional-family constraint

GRF's and DINA's candidate families are generated from fitted surfaces, so the fixed-family
condition does not hold for them: **every `orgrf` and `ordina` coverage figure is coverage of the
estimand conditional on the proposed family.** FS's family is the prespecified cut grid. Every
table, caption and extract row in this directory says which, and comparisons across the three
identifiers are **descriptive** — the identifier, the family construction and the set of detected
replicates all differ, so rows are read side by side and never ranked.

## 6. Scope

These are operating characteristics on one binary design at three effect sizes. They do not verify
condition (A3), no construction is promoted on this design's performance, and they supersede
nothing in the committed study, which ran a different rule and different constructions.

## 7. Open work

- **Only the two ends of the study's sweep are run.** The study swept n = 500 … 2000 by 250; these
  campaigns run n = 500 and n = 2000. The five interior sizes are open under the current rule and
  constructions.
- **One design family.** H is the study's {wtkg > q70} ∩ {cd40 > q70} at a 9.6% prevalence in all
  three design points; the study driver notes an alternative single-cut definition (`cd40` at q75,
  ~25% prevalence) that has never been run under the field constructions.
- **No genuine global null.** The borderline-null design point plants the region *at* the null
  against a protective complement; the driver's homogeneous `dgm_model = "null"` branch, under
  which the declaration rate *is* the false-positive rate, is not exercised here.
- **The other effect measures on the binary path** (RR, RD) are untouched by this directory.
- **The identification summary engine is not used here.** `fs_identification_summary()` takes an
  anchor / partner / proxy triple, which is the MD design's age / preanti / str2; the binary true
  region has no binary proxy, so this directory reports the classification table, the
  covariate-frequency figure and the realized-rule table instead.
