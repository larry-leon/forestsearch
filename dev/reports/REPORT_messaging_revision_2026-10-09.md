# REPORT — Messaging revision; methods-paper citation (forestsearch 0.4.0)

- Date: 2026-10-09
- Task: `dev/tasks/TASK_messaging_revision_2026-10-09.md` (commit `4e962745`)
- Inventory used: `dev/reports/REPORT_messaging_inventory_2026-10-09.md`
  - Block ids below refer to it.
  - Every site was re-read from current source before editing. None had moved
    since the inventory, apart from the line shifts the edits caused.
- Branch: `feature/glm-extension`

## Summary

- All steps were completed.
- Every gate passes:
  - Gates 2, 3, 3b, 3c, 4 and 7.
  - The test gate.
- Local `R CMD check --as-cran`: **`Status: 1 NOTE`**. That NOTE is the
  environment-only HTML Tidy NOTE that Stage 1 already recorded.
- Win-builder R-devel round 3 was uploaded once. Its result is pending by email.
- One deviation, approved by Larry mid-task: see "Deviation" below. The first
  check flagged the task's `<arXiv:2609.38361>` DESCRIPTION wording with a CRAN
  incoming NOTE. I switched it to the arXiv DOI form before the upload.

## Commits (in order)

| Commit | Content |
|---|---|
| `4e962745` | task doc |
| `a65f9764` | Step 1: printed output. Staged as hunks limited to `.fs_mr_caveats()`, `.fs_print_mr_products()` (including their `@noRd` headers) and `print.guohe_a3()` |
| `392253f0` | Step 2: roxygen in 9 `R/` files, plus 11 regenerated `man/*.Rd`. Every `R/` diff line is `#'` |
| `89c8dcd3` | Step 3: `vignettes/survival-post-selection.qmd` |
| `c4d75fa8` | Step 4: `vignettes/forestsearch.Rmd` |
| `c373fd83` | Step 5: `DESCRIPTION` and `NEWS.md`, plus the regenerated `man/forestsearch-package.Rd` |
| `876cbbbc` | First check log (`Status: 2 NOTEs`) |
| `2405f11f` | Deviation fix: arXiv cited by DOI in `DESCRIPTION`/`NEWS.md`; `forestsearch-package.Rd` regenerated |
| `94e58a4a` | Re-run check log (`Status: 1 NOTE`), overwriting the first |
| `fedf7a46` | `cran-comments.md`: round 3 pending |
| (this report) | report |

There was no test-expectation commit, because no test asserted the removed
text.

## Preconditions

- **P1:** `git status --porcelain` was clean. HEAD before the task-doc commit
  was `bc8bf525`.
- **P2:** I built a baseline from HEAD `4e962745`:
  - Installed into the scratch library `<scratch>/rlib_base` with
    `R CMD INSTALL`.
  - Ran the vignette's `setup`/`data`/`fit` code verbatim, with
    `R_LIBS=<scratch>/rlib_base`.
  - Saved `print(fs)` and `summary(fs)` as `baseline_print.txt` and
    `baseline_summary.txt`.
  - This fit has p-hat(H) = 0.006, so the baseline printed R1, R3 and R5. It
    did not print R4.

## Per-block before → after

### Printed output (Step 1)

- **R1** `forestsearch_methods.R` header: `Post-selection inference (certified products):` → `Post-selection inference (selection-adjusted; Leon & Anderson 2026):`
- **R2** `summary()` "Certified: … Source: dev/notes/NOTE_…" paragraph (6 `cat` lines) → deleted.
- **R3/R4/R5** the conditional caveats in `.fs_mr_caveats()`: the two-sided
  coverage caveat, the stable-pick note and the unstable-pick note → replaced
  by one fixed 5-line note.
  - It is printed whenever the block prints, in both `print()` and
    `summary()`.
  - It uses the same `.fs_cat_caveat()` path, so each line keeps the existing
    trailing space.
  - The argument `p` is kept for interface stability but is no longer read.
  - The function's `@noRd` header was rewritten to match. The stale
    NOTE-quoting `#` comments inside the function were removed with the code
    they described.
- **R6** `print.guohe_a3()`: `[CERTIFIED]` → `[Guo-He bound]`. `[EXTENSION]`
  and `[HEURISTIC]` are unchanged.
- `.fs_print_mr_products()` `@noRd` header: "certified-products block" →
  "post-selection block". The `long` description no longer mentions the
  deleted paragraph.

### Roxygen → man (Step 2)

- **D1** `print.forestsearch()` §Post-selection inference: certified-products / NOTE / pole wording (37 lines) → the task's 14-line body. It cites León and Anderson (2026) with an `\href` to arXiv.
- **D2** `summary.forestsearch()` §Post-selection inference: first paragraph ("certified-products block … sourced from dev/notes/…") → the task's 5-line body.
  - Interpretation: I replaced the D2 block as the inventory scoped it
    (`:396-406`). The following `sens_H` paragraph had no flagged content and
    is kept unchanged.
- **D3** `@param return_reselection`: "(they carry the stable-pick frequency the coverage caveat is read by)" → "(\eqn{\hat p} is reported by print() and summary() as a descriptive diagnostic)".
- **D4** `@param ci_method`: two edits.
  - "one-sided products are the certified ones" → "one-sided bounds are the primary reported products".
  - "(method proposal, `dev/tasks/TASK_mr_field_vs_guohe_2026-09-05.md`)" → "(León and Anderson, 2026)".
- **D5** `@param field_uniform`: two edits.
  - The method-proposal `dev/tasks/…` path was removed.
  - "Guarantee as documented there:" → "Experimental option. Intended properties:".
- **D6** `@param field_complement`: "(method proposal `dev/tasks/TASK_mr_field_complement_2026-09-06.md`)" removed.
- **D7** `@param field_decompose`: "; PROPOSAL_complement_field_scale 2026-09-08 v2, Stage 1" removed.
- **D8** `@param field_scale_complement`: "; PROPOSAL_complement_ field_scale_2026-09-08_v2 s5, variant R1" removed. "(the percentile-t root)" is kept.
- **D9** `@param ij_residual`: three edits.
  - The method-proposal path was removed.
  - "which repairs the tie regime at zero cost" → "retained for comparison only and not recommended".
  - ", likewise retained for comparison only" was appended after the `"winner"` description.
- **D10** `@param field_recovery`: "(`dev/tasks/TASK_field_recovery_2026-09-09.md`, R-1)" removed.
- **D11** `forestsearch()` `mr_inference_args$ci_method`: "one-sided products are the certified ones" → "one-sided bounds are the primary reported products".
- **D12** `fs_sim_bias_coverage()`: two edits.
  - "Conventions follow the committed template's diagnostics section:" → "Conventions:".
  - The `TASK_mr_field_complement_2026-09-06` reference was removed.
- **D13** `fs_sim_bias_coverage()`: "; TASK_complement_refinements_2026-09-06" removed.
- **D14** `fs_sim_bias_coverage()`: "(TASK_continuous_field_mac_2026-09-07)" removed.
  - "Add-only" now ends with a colon, so the sentence still reads.
- **D15** `fs_oc_predict()`: "(\code{dev/glm-continuous-sims/sigma_d_diagnostic_2026-08-29.R})" removed.
  - Now reads "within a few percent with no prevalence trend, so the wrapper identifies …".
- **D16** `fs_mr_oc_summary()`: three edits.
  - "used by the analytic documents" → "used in simulation summaries".
  - "exactly as the batch documents apply it" → "as the simulation harness applies it".
  - "This is what the batch documents report." → "This is the detected-only summary."
- **D17** `guohe_algorithm3()` `@return`: "the certified one-sided bound" → "the one-sided bound of Guo and He's Algorithm 3".
  - `print.guohe_a3()`: "(certified / extension / heuristic)" → "(Guo-He bound / extension / heuristic)".
- **D18** `fs_betaHhat_one()`: "This is the estimand Theorem 2 certifies coverage for" → "This is the conditional estimand of León and Anderson (2026, Theorem 2)".
- **D19** `.consistency_glm_pieces()`: "is unresolved; see `dev/identifier-alignment/f4_ps_effect_gap_findings.qmd`." → "is unresolved."
- **D20**: no change, per the task.

`roxygen2::roxygenise()` rewrote exactly these files and nothing else:

- `dot-consistency_glm_pieces.Rd`
- `forestsearch.Rd`
- `fs_betaHhat_one.Rd`
- `fs_mr_inference.Rd`
- `fs_mr_oc_summary.Rd`
- `fs_oc_predict.Rd`
- `fs_sim_bias_coverage.Rd`
- `guohe_algorithm3.Rd`
- `print.forestsearch.Rd`
- `print.guohe_a3.Rd`
- `summary.forestsearch.Rd`
- (Step 5) `forestsearch-package.Rd`

### `vignettes/survival-post-selection.qmd` (Step 3)

Code chunks are byte-identical to HEAD: I extracted every chunk before and
after and diffed them, with no difference.

- **V1** title and `%\VignetteIndexEntry`: "Survival post-selection inference: the certified products" → "Post-selection inference for an identified survival subgroup".
- **V2+V3** introduction: the "declines to certify" and "dev/notes/… source of truth" paragraphs → the task's single paragraph citing León and Anderson (2026).
- "The data and the fit": the "applications configuration cut to vignette size …" paragraph → the task's paragraph, which points to `vignette("forestsearch")`.
- "The identified subgroup": the "weakest form of selection …" sentence → the task's `family_status` sentence.
- **V4**: heading "The certified products" → "The reported quantities".
  - The three coverage paragraphs and the `ci_method` paragraph → the task's 5-bullet list plus the `ci_method` paragraph.
  - The `print(fs)` lead-in and chunk are kept.
- **V5**: → "`summary()` adds the standard errors and the leading re-selection frequencies:".
- **V6**: "-- 0.80 and 0.85 are the usual reading aids in these records --" → "-- for example 0.80 or 0.85 --".
- **V7/V8/V9**: bias bullets, band-edge paragraph and low-pole paragraph → deleted, replaced by the task's 4-line p-hat paragraph.
- "What p-hat alone cannot say": "1,744 candidates … 0.006" → "a large family … a small exact-match rate".
  - The "So the two numbers say different things …" paragraph → the task's paragraph.
  - "One limit, stated plainly." → "One limit."
- **V10/V11/V12**: section "What is not certified" → deleted.
- **V13**: section "Where the evidence lives" → "## References", with the two entries.
  - **Recorded difference in the 2024 entry:** I checked it against the
    package's own citations in `README.md` (BibTeX) and `forestsearch.Rmd`
    References. Authors, journal and DOI match the task text. The package
    capitalises the subtitle "**A** relatively simple procedure", where the
    task has "a relatively simple procedure". I used the package's version.

### `vignettes/forestsearch.Rmd` (Step 4)

- `?survival::gbsg` help text line: `‘hormon’ hormonal therapy, 0= no, 1= yes`.
  This confirms the comparison is hormonal therapy versus none.
- `:115` "versus chemotherapy" → "versus no hormonal therapy"
- `:121` "vs. chemotherapy" → "vs. no hormonal therapy"
- `:126` "# Load GBSG data (included in forestsearch package)" → "# Load GBSG data (from the survival package)"
- `:792` "95%95% CIs" → "95% CIs"
- `:864` `sprintf("%.2f", exp(fs$grp.consistency$out_sg$result$hr[1]))` → `sprintf("%.2f", fs$grp.consistency$out_sg$result$hr[1])`

### DESCRIPTION / NEWS (Step 5)

- `Description:` gains "Post-selection inference for the identified subgroup
  follows León and Anderson (2026) <doi:10.48550/arXiv.2609.38361>." It sits
  after the `<doi:10.1002/sim.10163>` sentence. The DOI form follows the
  deviation below.
- `NEWS.md`, `# forestsearch 0.4.0` → `### New features`, last bullet:
  "* Printed post-selection output and documentation simplified; methodology
  reference León and Anderson (2026) <doi:10.48550/arXiv.2609.38361>."
- `Version:` remains `0.4.0`.

## Deviation (approved)

- The first local check (log committed at `876cbbbc`) returned
  `Status: 2 NOTEs`. The new one came from the task-specified DESCRIPTION text:

  ```
  * checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Larry Leon <larry.leon.05@post.harvard.edu>'

  The Description field contains
    subgroup follows León and Anderson (2026) <arXiv:2609.38361>.
  Please refer to arXiv e-prints via their arXiv DOI <doi:10.48550/arXiv.YYMM.NNNNN>.
  ```

- I paused before the upload and asked. Larry chose "Fix, re-check, upload".
- `<arXiv:2609.38361>` → `<doi:10.48550/arXiv.2609.38361>` in `DESCRIPTION`
  and `NEWS.md` (commit `2405f11f`). `forestsearch-package.Rd` was regenerated
  and now uses `\doi{}`.
- Rebuild and re-check: incoming feasibility is OK and the status is
  `1 NOTE`.
- Unchanged by the fix:
  - Rd and vignette text still read "arXiv:2609.38361", with the abs link.
  - Printed output still reads `Leon & Anderson (2026), arXiv:2609.38361`.

## Gates

### Gate 2 — printed numbers unchanged: PASS

- **Method:**
  - The working tree was installed into `<scratch>/rlib_new`.
  - The same script was run, writing `new_print.txt` and `new_summary.txt`.
  - For every baseline line containing a digit, I checked for an exact
    whole-line match in the new output.
- **Unmatched digit lines:** only R2–R5 text.
  - `print()`: 6 lines, all from R3/R5.
  - `summary()`: 8 lines, from R2/R3/R5.
- Every bound, SE, p-hat, top-3, `p_hat_sum`, membership-recovery and timing
  line is identical.

Gate 2 diff, `print()`:

```
16c16
< Post-selection inference (certified products):
---
> Post-selection inference (selection-adjusted; Leon & Anderson 2026):
23,31c23,27
<   Two-sided intervals are not certified: the IJ two-term interval's
<   harm-block coverage falls to 0.913-0.917 at 12.4% prevalence with n >= 1000
<   (0.971-0.981 at 31%).  Read two-sided statements at low prevalence and
<   large n with that caveat. 
<   p-hat(H) < 0.20 is the unstable-pick pole, where the harm-block correction
<   runs the other way -- over-correction, -0.11 to -0.28 log units at 12.4%
<   prevalence with n = 1500 -- so the one-sided lower bound on H is
<   conservative here and the two-sided interval is the exposed one.  See
<   dev/notes/NOTE_survival_products_2026-09-09.md. 
---
>   One-sided bounds are the primary products; the two-sided IJ interval is 
>   secondary. p-hat(H) is a descriptive diagnostic. Read each bound by its 
>   location relative to a clinically meaningful effect size. Methods and 
>   simulation operating characteristics: Leon & Anderson (2026), 
>   arXiv:2609.38361. 
```

Gate 2 diff, `summary()`:

```
41c41
< Post-selection inference (certified products):
---
> Post-selection inference (selection-adjusted; Leon & Anderson 2026):
56,61c56,60
<   Certified: the one-sided lower bound on H (field) and the one-sided
<   upper bound on Hc (field-s), and the Bonferroni joint pair at gamma =
<   0.025 each side.  Not certified: any two-sided interval.  p-hat(H) is a
<   recorded diagnostic -- no construction reads it.  Read every bound by
<   location against a clinically meaningful effect size, never as
<   significance at the null.  Source: dev/notes/NOTE_survival_products_2026-09-09.md.
---
>   One-sided bounds are the primary products; the two-sided IJ interval is 
>   secondary. p-hat(H) is a descriptive diagnostic. Read each bound by its 
>   location relative to a clinically meaningful effect size. Methods and 
>   simulation operating characteristics: Leon & Anderson (2026), 
>   arXiv:2609.38361. 
63,72d61
<   Two-sided intervals are not certified: the IJ two-term interval's
<   ... (R3 and R5 text, 9 lines, as in print() above)
```

The fixed note prints once in `summary()`. The deleted R2 paragraph did not
double it.

### Gate 3 — references gone: PASS

- `grep -rnE "dev/|TASK_|PROPOSAL_|REPORT_|NOTE_" man/ vignettes/forestsearch.Rmd vignettes/survival-post-selection.qmd` → no hits.
- `grep -rniE "certif" man/ vignettes/forestsearch.Rmd vignettes/survival-post-selection.qmd` → no hits.
- `R/` still has 58 hit lines for either pattern. I classified each line by
  script:
  - **51 in `#` comments.** Examples:
    - `forestsearch_methods.R:117` (the `TASK_field_recovery` comment in
      `.fs_mr_products()`)
    - `forestsearch_methods.R:309` (`print.forestsearch()`'s `See dev/notes/…`
      comment)
    - `forestsearch_main.R:3799-3844` (14 MR-engine comments)
  - **7 in `@noRd` blocks:**
    - `consistency_resample.R:125` (`.dfbeta_glm`)
    - `forestsearch_methods.R:60`, `:66`, `:79` (`.fs_mr_products()`,
      "certified")
    - `fs_mr_inference.R:1146`, `:1234`, `:1425`
  - **0 in code, 0 in exported roxygen.**
  - The full classified list is in the scratchpad (`g3_classified.txt`).
    The counts are reproducible with the two greps above over `R/`.

### Gate 3b: PASS

- In `R/`, every changed line that is not `#'` lies inside
  `.fs_mr_caveats()`, `.fs_print_mr_products()` or `print.guohe_a3()`.
  - Commit `a65f9764` holds all of them.
  - Commit `392253f0` has 0 non-`#'` lines.
- `NAMESPACE` is byte-identical to the pre-task copy. Checked after each of
  the 3 `roxygenise()` runs.

### Gate 3c — ASCII: PASS

- `tools::showNonASCIIfile()` on the 9 changed `R/` files reports 3 lines, all
  roxygen. Each is "León":
  - `betaHhat_truth.R:212`
  - `forestsearch_methods.R:224`
  - `fs_mr_inference.R:315`
- No non-ASCII outside roxygen. Printed strings use "Leon".

### Gate 4 — forestsearch.Rmd "Unadjusted HR": PASS

- Source: the rendered `forestsearch.html` from the release tarball's
  `inst/doc`.
- Key Findings prints `Unadjusted HR: 1.95` for "Definition: {er <= 0}",
  "Size: 82 (12.0%)".
- The same subgroup's HR elsewhere in the document:
  - top-subgroups listing rank 1: `1.951  82  45  … {er <= 0}`
  - identified-subgroups table: `{er <= 0} … 82 45 16 1.951 0.960`
  - estimates table, Questionable row: `HR (95% CI) 1.95 (1.04, 3.67)`
  - bootstrap summary: `Observed HR 1.951`
- These are equal at 2 dp. Before the fix, the line printed `exp(1.951)`,
  roughly 7.04.
- The rendered HTML no longer contains "chemotherapy" or "95%95%".

### Test suite (Step 6): PASS against the gate

- **Command:** `NOT_CRAN=true R_LIBS=<scratch>/rlib_new Rscript -e 'devtools::test()'`.
  - `rlib_new` holds the working tree installed, so `multisession` workers
    load the same build: `find.package` resolved to the scratch library.
  - The user library was not touched.
- **Result:** `[ FAIL 6 | WARN 38 | SKIP 3 | PASS 5868 ]`, duration 201.6 s.
- **The 6 failures are exactly the pre-existing fs-oc reference-file failures**
  (F1/F2 of `REPORT_cran_release_wb_followup_2026-10-09.md`):
  - `test-fs-oc-grid.R`: 4 failures in "fs_oc_predict() is identical to the
    0.2.4 pre-refactor reference on both gates".
  - `test-fs-oc-null.R`: 2 failures in "the alternative path is untouched by
    the null branch".
- `search-reproducibility` (worker invariance): ✔ 6/6.
- **Test-expectation changes: none.** No test asserts the removed or rewritten
  text. I grepped `tests/testthat/` for "certified", "stable-pick", "0.913",
  "NOTE_survival", "caveat" and "Post-selection inference": no hits.

### Gate 7 — local release check: PASS

- Environment as in Stage 1, Step 9.
  - Temp library `<scratch>/rlib` with CRAN `weightedsurv` 0.1.0
    (`Repository: CRAN`; it resolves ahead of the GitHub copy) and `latentcor`
    2.0.2.
  - latentcor's recursive dependencies were installed as CRAN macOS
    **binaries** into the same temp library, as in Stage 1, because there is
    no gfortran:
    `pcaPP fMultivar mnormt heatmaply plotly geometry fBasics sn magic lpSolve linprog RcppProgress dendextend reshape2 seriation webshot assertthat egg crosstalk timeSeries gss stabledist ca gclus qap registry TSP`.
- Build: `R CMD build <repo>` with vignettes, in `<scratch>/chk`.
  - Took 51 s elapsed (first run).
  - Produced `forestsearch_0.4.0.tar.gz`, 3,343,864 bytes.
- Check: `R CMD check --as-cran` with
  `R_LIBS=<scratch>/rlib _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE _R_CHECK_LIMIT_CORES_=TRUE OMP_THREAD_LIMIT=2 _R_CHECK_EXAMPLE_TIMING_CPU_TO_ELAPSED_THRESHOLD_=2.5 _R_CHECK_TEST_TIMING_CPU_TO_ELAPSED_THRESHOLD_=2.5`.
  - First run (before the DOI fix): `Status: 2 NOTEs`. The first NOTE is quoted
    under "Deviation"; the second is the HTML Tidy NOTE below.
  - Final run, after the DOI fix:
    - **`Status: 1 NOTE`**
    - `checking CRAN incoming feasibility ... OK`
    - `checking tests ... [19s/17s] OK`
    - `checking PDF version of manual ... OK`, so the "León" Rd text builds in
      LaTeX.
    - The NOTE, verbatim:

      ```
      * checking HTML version of manual ... NOTE
      Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
      Please obtain a recent version of HTML Tidy by downloading a binary
      release or compiling the source code from <https://www.html-tidy.org/>.
      ```

      This is environment-only, identical to Stage 1, and already explained.
- Log: `dev/reports/cran_messaging_00check_2026-10-09.log` holds the final
  run.

## Step 8 — win-builder R-devel round 3

- `devtools::check_win_devel(pkg = ".")` was run once, after commit `94e58a4a`.
  Console output:

  ```
  Building windows version of forestsearch (0.4.0)
  ℹ Using R-devel with win-builder.r-project.org.
  ── R CMD build ─────────────────────────────────────────────────────────────────
  * checking for file '/Users/larryleon/Documents/GitHub/forestsearch/DESCRIPTION' ... OK
  * preparing 'forestsearch':
  * checking DESCRIPTION meta-information ... OK
  * installing the package to build vignettes
  * creating vignettes ... OK
  * checking for LF line-endings in source and make files and shell scripts
  * checking for empty or unneeded directories
  * building 'forestsearch_0.4.0.tar.gz'
  ────────────────────────────────────────────────────────────────────────────────
  ℹ Check <'larry.leon.05@post.harvard.edu'> for the results in 15-30 mins
    (~09:04 PM).
  ```

- I did not wait for the email.
- `cran-comments.md` now has the line
  `* win-builder R-devel round 3 (after documentation revision): pending`
  (commit `fedf7a46`).
- There was no CRAN submission and no push.

## Post-conditions

- **PC1:** `git diff --name-only 4e962745^ HEAD` lists:
  - the task doc and this report
  - `dev/reports/cran_messaging_00check_2026-10-09.log`
  - `DESCRIPTION`, `NEWS.md`, `cran-comments.md`
  - the 9 `R/` files:
    - `betaHhat_truth.R`
    - `consistency_resample.R`
    - `forestsearch_main.R`
    - `forestsearch_methods.R`
    - `fs_bias_coverage.R`
    - `fs_mr_inference.R`
    - `fs_mr_oc_summary.R`
    - `fs_oc_predict.R`
    - `guohe_algorithm3.R`
  - 12 regenerated `man/*.Rd` files
  - the two vignettes

  No test files.
- **PC2:** `NAMESPACE` is byte-identical. `Version: 0.4.0`.
- **PC3:** `git status --porcelain` is clean.
  - The tarball, the `.Rcheck` directory and all scratch libraries
    (`rlib`, `rlib_base`, `rlib_new`) are under the session scratchpad, not
    the repo.
  - They will be deleted at the end of the task.
- Copied to `~/Downloads`:
  - this report
  - `survival-post-selection.html` and `forestsearch.html`, from the final
    release tarball's `inst/doc`, overwriting the earlier copies

## Side observations (not acted on)

1. The `R/` `#` comments still cite `dev/notes/NOTE_survival_products_…` and
   `TASK_*` (51 lines). An example is `print.forestsearch()`'s block comment
   above `.fs_mr_products()`. They never reach output, and Gate 3 allows them.
   They could be swept in a later cleanup task if desired.
2. `.fs_mr_products()`'s `@noRd` header still says "certified post-selection
   products". It is internal-only and was outside Step 1's function list, so I
   left it.
3. `forestsearch.Rmd:858` still hard-codes the ITT result
   `"HR (95% CI): 0.69 (0.54, 0.89)"`. It matches the rendered ITT row
   `0.69 (0.54, 0.89)`, but it is not computed.
