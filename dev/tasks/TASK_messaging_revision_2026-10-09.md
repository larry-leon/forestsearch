# TASK — Tone down user-facing interpretive text; cite the methods paper (forestsearch 0.4.0)

- Date: 2026-10-09
- Repo: `larry-leon/forestsearch`, branch `feature/glm-extension`, at or after the inventory report commit.
- Written against: `dev/reports/REPORT_messaging_inventory_2026-10-09.md` (block ids R1–R6, D1–D20, V0–V13 below refer to it). Re-read every cited line from current source before editing; if a site has moved, use the current one and record it.
- Citation used throughout:
  - Rd / vignettes: León, L. F. and Anderson, K. M. (2026). Inference for standard trial estimands after data-driven subgroup discovery. arXiv:2609.38361 (link `https://arxiv.org/abs/2609.38361`).
  - R strings (printed output) must stay ASCII: write `Leon & Anderson (2026), arXiv:2609.38361`.
- Gates are stop-on-failure.

## R/ callout

- Printed output (`R/forestsearch_methods.R`, `R/guohe_algorithm3.R`): **messaging only.** Header text, caveat paragraphs and one label change; the p-hat-conditional caveats are removed. No computed value, no element of any returned object, and no printed number changes (Gate 2).
- Roxygen in the R/ files listed below: **documentation only** (every diff line `#'`).
- Nothing changes behaviour or the method.

## Authorizations (kickoff)

- One local `R CMD build` + `R CMD check --as-cran` (Step 7) with the Stage 1 temp-library setup (CRAN `weightedsurv`, `latentcor`; binaries allowed for latentcor's Fortran dependencies; scratch library only).
- Running the full test suite once (Step 6).
- One upload: `devtools::check_win_devel(pkg = ".")` (Step 8). No CRAN submission, no push.
- Editing test expectations **only** where a test asserts text this task removes or rewrites (record each).

## First action

1. Copy this file to `dev/tasks/TASK_messaging_revision_2026-10-09.md` and commit it alone.

## Preconditions

- P1. `git status --porcelain` clean. Record HEAD.
- P2. Baseline for Gate 2: render the `print` and `summary` chunks of `vignettes/survival-post-selection.qmd` (same code, installed HEAD in a scratch library) and save the two outputs to the scratchpad as `baseline_print.txt`, `baseline_summary.txt`.

## Step 1 — Printed output (`R/forestsearch_methods.R`, `R/guohe_algorithm3.R`)

- **R1** header (`:196`): replace
  `"\nPost-selection inference (certified products):\n"`
  with
  `"\nPost-selection inference (selection-adjusted; Leon & Anderson 2026):\n"`.
- **R3, R4, R5** (`.fs_mr_caveats()`, `:146-172`): remove all three conditional caveats (the two-sided coverage caveat and both p-hat pole notes). Replace with one fixed note, printed once whenever the MR block prints (both `print()` and `summary()`):
  ```
    One-sided bounds are the primary products; the two-sided IJ interval is
    secondary. p-hat(H) is a descriptive diagnostic. Read each bound by its
    location relative to a clinically meaningful effect size. Methods and
    simulation operating characteristics: Leon & Anderson (2026),
    arXiv:2609.38361.
  ```
- **R2** (`summary()` paragraph, `:230-235`): delete it (the fixed note above replaces it; it must not print twice).
- **R6** (`print.guohe_a3()`, `:515`): change the label `[CERTIFIED]` to `[Guo-He bound]`. `[EXTENSION]` and `[HEURISTIC]` unchanged.
- Leave every other printed line (the bound lines, SE lines, re-selection lines) exactly as is.

## Step 2 — Roxygen (documentation only); regenerate `man/`

Edit the roxygen source; then `roxygen2::roxygenise()`.

- **D1** `print.forestsearch()` `@section Post-selection inference:` — replace the whole section body with:
  ```
  #' When the fit carries multiplier-resampling results (`mr_inference = TRUE`,
  #' attached as `x$mr_inference`) and the field block ran
  #' (`ci_method = "field"`, the default), a further block reports: the
  #' one-sided 95% lower bound for the identified subgroup (the field); the
  #' one-sided 95% upper bound for its complement (field-s, the studentized
  #' complement field); the Bonferroni pair for a joint two-subgroup
  #' statement; the infinitesimal-jackknife (IJ) two-term two-sided interval,
  #' as a secondary summary; and the re-selection frequency
  #' \eqn{\hat p(\widehat H)} as a descriptive diagnostic.  Bounds are read by
  #' their location relative to a clinically meaningful effect size, not as
  #' tests at the null.  With `ci_method = "ij"` only the two-sided interval is
  #' shown; without MR results the output is unchanged.  The constructions and
  #' their simulation operating characteristics are described in León and
  #' Anderson (2026), [arXiv:2609.38361](https://arxiv.org/abs/2609.38361).
  ```
- **D2** `summary.forestsearch()` section — replace body with:
  ```
  #' Carries the same block as [print.forestsearch()] (see there) and adds the
  #' log-scale standard errors (field, field-s with the naive complement SE
  #' beside it, and IJ two-term), the three largest re-selection frequencies
  #' with their labels, and `p_hat_sum` over the re-selection family.  Absent MR
  #' results, the output is unchanged.
  ```
- **D3** `@param return_reselection`: replace `(they carry the stable-pick frequency the coverage caveat is read by)` with `(\eqn{\hat p} is reported by print() and summary() as a descriptive diagnostic)`.
- **D4** `@param ci_method`: `whose one-sided products are the certified ones` → `whose one-sided bounds are the primary reported products`; replace `(method proposal, \`dev/tasks/TASK_mr_field_vs_guohe_2026-09-05.md\`)` with `(León and Anderson, 2026)`.
- **D5** `@param field_uniform`: delete `, method proposal \`dev/tasks/TASK_mr_field_uniform_2026-09-05.md\``; replace `Guarantee as documented there:` with `Experimental option. Intended properties:`.
- **D6** `@param field_complement`: delete `(method proposal \`dev/tasks/TASK_mr_field_complement_2026-09-06.md\`)`.
- **D7** `@param field_decompose`: delete `; PROPOSAL_complement_field_scale 2026-09-08 v2, Stage 1`.
- **D8** `@param field_scale_complement`: delete `; PROPOSAL_complement_ field_scale_2026-09-08_v2 s5, variant R1` (keep `(the percentile-t root)`).
- **D9** `@param ij_residual`: delete `method proposal \`dev/tasks/TASK_complement_refinements_2026-09-06.md\`,`; replace `which repairs the tie regime at zero cost` with `retained for comparison only and not recommended`; append after the `"winner"` description `, likewise retained for comparison only`.
- **D10** `@param field_recovery`: delete `(\`dev/tasks/TASK_field_recovery_2026-09-09.md\`, R-1)`.
- **D11** `forestsearch()` `mr_inference_args$ci_method`: `whose one-sided products are the certified ones` → `whose one-sided bounds are the primary reported products`.
- **D12** `fs_sim_bias_coverage()`: `Conventions follow the committed template's diagnostics section:` → `Conventions:`; delete `, TASK_mr_field_complement_2026-09-06`.
- **D13**: delete `; TASK_complement_refinements_2026-09-06`.
- **D14**: delete `(TASK_continuous_field_mac_2026-09-07)`.
- **D15** `fs_oc_predict()`: delete `(\code{dev/glm-continuous-sims/sigma_d_diagnostic_2026-08-29.R})` and the comma before it, keeping the sentence grammatical.
- **D16** `fs_mr_oc_summary()`: `used by the analytic documents` → `used in simulation summaries`; `exactly as the batch documents apply it` → `as the simulation harness applies it`; `This is what the batch documents report.` → `This is the detected-only summary.`
- **D17** `guohe_algorithm3()` `@return`: `the certified one-sided bound` → `the one-sided bound of Guo and He's Algorithm 3`; `print.guohe_a3()` description: `labelled by their inferential status (certified / extension / heuristic)` → `labelled by their status (Guo-He bound / extension / heuristic)`.
- **D18** `fs_betaHhat_one()`: `This is the estimand Theorem 2 certifies coverage for, evaluated at` → `This is the conditional estimand of León and Anderson (2026, Theorem 2), evaluated at`.
- **D19** `.consistency_glm_pieces()`: replace `is unresolved; see \`dev/identifier-alignment/f4_ps_effect_gap_findings.qmd\`.` with `is unresolved.`
- **D20**: no change.

## Step 3 — `vignettes/survival-post-selection.qmd`

Keep every code chunk exactly as it is. Replace prose as follows (headings shown with `##`).

- **V1** title and `%\VignetteIndexEntry`: `Post-selection inference for an identified survival subgroup`.
- **V2 + V3** (introduction, both paragraphs) → 
  ```
  `forestsearch()` identifies a subgroup by searching over candidate subgroups,
  so a treatment effect estimated within the identified subgroup is affected by
  the selection: a confidence interval computed as if the subgroup had been
  prespecified tends to under-cover. With `mr_inference = TRUE`,
  `forestsearch()` reports selection-adjusted estimates and bounds for the
  identified subgroup and its complement, based on multiplier resampling of the
  candidate effects. The methodology, its assumptions and its simulation
  operating characteristics are described in León and Anderson (2026); this
  vignette shows how to obtain and read the reported quantities on a public
  dataset.
  ```
- In "The data and the fit": replace the paragraph beginning `The fit below is the applications configuration cut to vignette size` with:
  ```
  The fit below uses the consistency search without a LASSO/GRF/DINA screen,
  `sg_focus = "effMaxSG"` (the largest subgroup within 20% of the strongest
  observed effect), 200 consistency splits and 500 multiplier draws.
  `mr_inference = TRUE` attaches the post-selection results to the fit as
  `fs$mr_inference`; `field_recovery = TRUE` adds the membership-agreement
  diagnostics used below. Bootstrap and cross-validation are not run here; see
  `vignette("forestsearch")`.
  ```
- In "The identified subgroup": replace the sentence `family_status records how that family was built -- "no-front-end" here, meaning no fitted model on the observed data shaped it, which is the weakest form of selection the package will report on.` with `` `family_status` records how the candidate family was built; "no-front-end" means the family was enumerated directly rather than proposed by a model fitted to the data. ``
- **V4** heading `## The certified products` → `## The reported quantities`. Keep the `print(fs)` lead-in sentence and chunk. Replace the three coverage paragraphs and the `ci_method` paragraph with:
  ```
  - **Identified subgroup:** the one-sided 95% lower bound on the hazard ratio
    (the *field* construction).
  - **Complement:** the one-sided 95% upper bound on the hazard ratio (*field-s*,
    the studentized complement field), relevant to a benefit claim.
  - **Both together:** the Bonferroni pair, each side at level 0.025.
  - **Two-sided interval:** the infinitesimal-jackknife (IJ) two-term interval,
    reported as a secondary summary.
  - **p-hat(H):** the re-selection frequency, a descriptive diagnostic (below).

  The field bounds require `ci_method = "field"` (the default); with
  `ci_method = "ij"` only the two-sided interval is reported. León and Anderson
  (2026) describe each construction and report its performance in simulation.
  ```
- **V5**: `` `summary()` adds the standard errors and the leading re-selection frequencies: ``
- **V6**: `-- 0.80 and 0.85 are the usual reading aids in these records --` → `-- for example 0.80 or 0.85 --`.
- **V7, V8, V9**: delete all three (the bias bullets, the band-edge paragraph and the "This fit sits at the low pole" paragraph). In their place, after the `p-hat` code chunk:
  ```
  A small p-hat indicates that the exact boundary of the identified subgroup is
  not stable under perturbation; by itself it is neither evidence for nor
  against the bounds. León and Anderson (2026) discuss how the bias correction
  behaves across values of p-hat.
  ```
- "What p-hat alone cannot say": remove the fit-specific numbers from prose. Replace `Over a family of 1,744 candidates, most of which differ from Hhat by one cutpoint, an exact-match rate of 0.006 is compatible with` with `Over a large family in which most candidates differ from Hhat by a single cutpoint, a small exact-match rate is compatible with`. Replace the paragraph beginning `So the two numbers say different things` with:
  ```
  Read together, the two kinds of number say different things about the same
  draws: exact re-selection can be rare while the typical re-selection still
  retains much of the identified subgroup. That pattern points to an unstable
  boundary rather than an unstable region.
  ```
  Keep the "One limit" paragraph, changing `One limit, stated plainly.` to `One limit.`
- **V10, V11, V12**: delete the whole section `## What is not certified`.
- **V13**: replace the section `## Where the evidence lives` and its list with:
  ```
  ## References

  León, L. F. and Anderson, K. M. (2026). Inference for standard trial
  estimands after data-driven subgroup discovery. arXiv:2609.38361.
  <https://arxiv.org/abs/2609.38361>

  León, L. F., Jemielita, T., Guo, Z., Marceau West, R. and Anderson, K. M.
  (2024). Exploratory subgroup identification in the heterogeneous Cox model:
  a relatively simple procedure. *Statistics in Medicine*.
  <https://doi.org/10.1002/sim.10163>
  ```
  (Before writing the 2024 entry, confirm its title and authors from the package's own citation of it in `DESCRIPTION`, `README.md` or `vignettes/forestsearch.Rmd`; if they differ, use the package's version and record the difference.)

## Step 4 — `vignettes/forestsearch.Rmd`

- `:115`, `:121`: replace "versus chemotherapy" / "vs. chemotherapy" with "versus no hormonal therapy" / "vs. no hormonal therapy". Before editing, confirm from `?survival::gbsg` that `hormon` is hormonal therapy yes/no; record the help text line.
- `:126`: `# Load GBSG data (from the survival package)`.
- `:864`: remove `exp()` so the line prints `fs$grp.consistency$out_sg$result$hr[1]` with `sprintf("%.2f", ...)`. Gate 4: in the rendered vignette the printed "Unadjusted HR" equals the identified subgroup's HR shown elsewhere in the same document (record both).
- `:792`: `95%95%` → `95%`.

## Step 5 — DESCRIPTION and NEWS

- `Description:` after the existing `<doi:10.1002/sim.10163>` sentence add: `Post-selection inference for the identified subgroup follows León and Anderson (2026) <arXiv:2609.38361>.`
- `NEWS.md`, under `# forestsearch 0.4.0` → `### New features`, add: `* Printed post-selection output and documentation simplified; methodology reference León and Anderson (2026) <arXiv:2609.38361>.`

## Step 6 — Gates

- **Gate 2 (printed numbers unchanged):** re-render the vignette's `print` and `summary` chunks as in P2. Every line containing a digit in the baseline outputs appears identically in the new outputs, except lines removed by Step 1 (R2–R5 text). Record the diff.
- **Gate 3 (references gone):** `grep -rnE "dev/|TASK_|PROPOSAL_|REPORT_|NOTE_" man/ vignettes/forestsearch.Rmd vignettes/survival-post-selection.qmd` → no hits. `grep -rniE "certif" man/ vignettes/forestsearch.Rmd vignettes/survival-post-selection.qmd` → no hits. In `R/`, any remaining hits for either pattern are in `#` comments or `@noRd` blocks only (list them).
- **Gate 3b:** in `R/`, every changed line is either `#'` roxygen or within `.fs_mr_caveats()` / `.fs_print_mr_products()` / `print.guohe_a3()`; `NAMESPACE` byte-identical.
- **Gate 3c (ASCII):** `tools::showNonASCIIfile()` on every changed `R/` file reports nothing outside roxygen.
- Full test suite once (`devtools::test()`, `NOT_CRAN=true`, installed HEAD so workers match): `FAIL` count equals the 6 pre-existing `fs-oc` reference-file failures only (the worker-invariance test now passes after the reinstall). Update test expectations only where a test asserts removed/rewritten text; list each.

## Step 7 — Local release check

Same commands and environment as Stage 1 Step 9. Gate 7: 0 ERRORs, 0 WARNINGs; record NOTEs verbatim. Copy `00check.log` to `dev/reports/cran_messaging_00check_2026-10-09.log`.

## Step 8 — Commit, then one win-builder R-devel upload

- Commits (explicit paths): printed output; roxygen + man; survival vignette; forestsearch vignette; DESCRIPTION + NEWS; test expectation updates (if any); check log.
- Then `devtools::check_win_devel(pkg = ".")` once; record console output. Do not wait for the email.
- `cran-comments.md`: add `win-builder R-devel round 3 (after documentation revision): pending`; commit.

## Deliverable

`dev/reports/REPORT_messaging_revision_2026-10-09.md`: per-block before/after (one line each), gate results, Gate 2 diff, test changes, check status, upload output. Commit, then copy to `~/Downloads`. Also copy the re-rendered `survival-post-selection.html` and `forestsearch.html` to `~/Downloads` (overwriting the earlier copies).

## Post-conditions

- PC1. Changed files: this task doc, the report, the check log, `DESCRIPTION`, `NEWS.md`, `cran-comments.md`, the R/ files named in Steps 1–2, regenerated `man/*.Rd`, the two vignettes, and any test files listed in Step 6.
- PC2. `NAMESPACE` byte-identical; `Version: 0.4.0`.
- PC3. `git status --porcelain` clean; no tarball, `.Rcheck` or scratch library in the repo.
