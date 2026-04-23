# forestsearch Quarto documents — catalog

This file catalogs every `.qmd` document in `applications/` and
`simulations/`, flags which are the current canonical versions vs.
older drafts, and identifies what's glaringly problematic. It is
not exhaustive code review — just enough for a reader to understand
the landscape.

Last reviewed: v0.2.0-dev branch, April 2026.


## Quick summary

**22 documents total** across `applications/` (14) and `simulations/`
(8). Roughly four groupings:

1. Canonical analyses of real data (GBSG, COPD demo) — 4 files
2. Canonical simulation studies (ACTG175, León 2024 reproduction) — 4 files
3. Validation / methodology tests — 3 files
4. Older drafts / near-duplicates — 11 files (candidates for removal or archival)


## Applications directory

### Current / canonical

| File | Lines | Purpose | Status |
|---|---:|---|---|
| `gbsg_analysis.qmd` | 2513 | **Main GBSG vignette.** Full walkthrough including FPR calibration appendix. Uses centralized parameter block ("All tunable parameters centralised here"). | **Current canonical** |
| `gbsg_analysis_cox_poisson.qmd` | 1097 | GBSG analysis comparing Cox HR and Poisson IRR on the same data — demonstrates GLM pipeline against survival pipeline. | **Current canonical for GLM comparison** |
| `count_data_demo.qmd` | 1086 | COPD exacerbation count demo (IRR with time offset). Exercises the count outcome path on synthetic data. | **Current canonical for count demo** |
| `biomarker_comparison.qmd` | 2084 | Continuous-biomarker methodology comparison: STEPP, MFPI, spline Cox on GBSG + simulated data. Not a ForestSearch tutorial — a comparative methods paper. | **Current** |

### Validation / methods studies

| File | Lines | Purpose | Status |
|---|---:|---|---|
| `forestsearch_scenario_validation.qmd` | 1417 | End-to-end simulation tests covering all four outcome types (binary, continuous, count, survival) across search configurations. | **Current canonical validation** |
| `validation_hte_tests_crump.qmd` | 1510 | Pre-screening for treatment effect heterogeneity via Crump et al. (2008) tests: OLS, logistic, Poisson+offset. Standalone methodology paper. | **Current** |
| `validation_glm_simulation_study.qmd` | 958 | GLM operating characteristics with calibrated thresholds (binary, continuous, count). | **Current** |

### Duplicates / older drafts

| File | Lines | vs. Canonical | Recommendation |
|---|---:|---|---|
| `gbsg_analysis_working.qmd` | 875 | **Earlier draft** of `gbsg_analysis.qmd`. Missing centralized parameter block. Later timestamp (Apr 22) but older structure; apparently kept as a development scratchpad. | **Archive or delete** unless actively used |
| `gbsg_analysis_working2.qmd` | 1343 | Intermediate draft between `working` and the main. Duplicates structure of the 2sim variant. | **Archive or delete** |
| `gbsg_analysis_2sim-approaches.qmd` | 1575 | Near-identical to `working2` but with `code-fold: show` and `warning: true`. Appears to be a development variant. | **Archive or delete** |
| `gbsg_analysis_cox_poisson_1.qmd` | 1097 | **Older version** of `gbsg_analysis_cox_poisson.qmd`. Only substantive diff: `grf_tune = TRUE` (the `_1`) vs `FALSE` (current). Same content otherwise. | **Delete** — cleanly superseded |
| `forestsearch_scenario_validation_100.qmd` | 1128 | Same structure as `forestsearch_scenario_validation.qmd` but with reduced N=100 and fewer sims for a quick version. | Either **delete and add a quick-mode param** to the canonical file, or keep but rename `*_quick.qmd` so intent is clear |
| `gbsg_poisson_simulations.qmd` | 710 | Duplicate of the same-named file in `simulations/` (67-line diff, appears to be copy kept here by mistake). | **Delete from `applications/`** — belongs in `simulations/` only |


## Simulations directory

### Current / canonical

| File | Lines | Purpose | Status |
|---|---:|---|---|
| `actg175_binary_benefit_simulations.qmd` | 3180 | **Canonical ACTG175 binary-benefit simulation.** Full paper-grade study at nsims_alt=5000, nsims_null=20000. Exercises the canonical GLM pattern (mixed continuous + factor covariates, is.RCT=TRUE, vi.grf.min=-0.2, sg_focus="maxSG"). | **Current canonical** |
| `actg175_continuous_simulations.qmd` | 756 | ACTG175 continuous outcome (CD4 change, MD). | **Current canonical** |
| `actg175_survival_benefit_simulations.qmd` | 2015 | ACTG175 survival outcome benefit search. | **Current canonical** |
| `gbsg_poisson_simulations.qmd` | 698 | Poisson rate simulation on GBSG. | **Current canonical** |
| `reproduce_table1_m1.qmd` | 1005 | Reproduction of León et al. (2024) Table 1, M1 scenario. Uses v0.2.0 collect_results as exported (`::`). | **Current canonical pressure test** |

### Duplicates / older drafts

| File | Lines | vs. Canonical | Recommendation |
|---|---:|---|---|
| `actg175_binary_benefit_simulations_1.qmd` | 1735 | **Much older version** of the canonical (Apr 17 vs Apr 19). Missing ~1400 lines of structure including Priority D diagnostic capture and subgroup-distribution plotting. | **Delete** — clearly superseded |
| `actg175_binary_benefit_simulations_100.qmd` | 3180 | Same as canonical but with `nsims_alt=100, nsims_null=100` for fast render. | Either **delete and add a `quick_mode` param** to the canonical, or keep but rename `*_quick.qmd` |
| `reproduce_table1_m1_1.qmd` | 1002 | Older version. Uses `forestsearch:::collect_results` (now exported as `::`). Also contains `dmin.grf = 4, frac.tau = 0.6` hardcodes the canonical has dropped. | **Delete** — cleanly superseded |
| `reproduce_table1_m1_2stage.qmd` | 1003 | Near-identical to canonical but with `use_twostage = TRUE` (canonical is FALSE) and `nsims = 100` (quick-render). Meaningful variant. | Either **keep but rename** `reproduce_table1_m1_twostage.qmd` to make intent clear, or **parameterize** in the canonical file |


## What's glaringly problematic

Beyond the duplicates, a few things stood out during the scan that are worth flagging:

### 1. Naming convention inconsistency for "quick-render" variants

`forestsearch_scenario_validation_100.qmd` and
`actg175_binary_benefit_simulations_100.qmd` use `_100` to denote
"quick-render with 100 simulations." Meanwhile `reproduce_table1_m1`
uses `_2stage` for a semantic variant. And `_1` appears to mean both
"older version" (gbsg_analysis_cox_poisson_1) and "first in a series."

The existing `smoke_*.qmd` and `test_*.qmd` files use a `params:` block
with `quick_mode: true/false` for the same purpose. This is the better
pattern. Converting the `_100` files to a single file with `params` would
eliminate two files and make intent clear.

### 2. `applications/gbsg_poisson_simulations.qmd` is misplaced

A copy of `simulations/gbsg_poisson_simulations.qmd` (with a 67-line
diff) exists in `applications/`. The right home is `simulations/`.

### 3. Four near-parallel versions of the GBSG analysis

`gbsg_analysis.qmd`, `gbsg_analysis_working.qmd`,
`gbsg_analysis_working2.qmd`, and `gbsg_analysis_2sim-approaches.qmd`
all share the same title, subtitle, and general structure. They are
at different stages of development. The canonical `gbsg_analysis.qmd`
is the most complete; the other three are drafts. Keeping all four
adds noise for anyone navigating the directory.

### 4. No version-in-header or README indicating which is canonical

When I landed in the directory cold, I had to diff files and read
timestamps to figure out which was canonical. A one-line comment at
the top of each canonical file (`# Canonical vignette — last updated
YYYY-MM-DD`) or a directory-level README would save future readers
and collaborators that detective work.


## What to keep — proposed minimal set

**Applications (keep 7):**
- `gbsg_analysis.qmd`
- `gbsg_analysis_cox_poisson.qmd`
- `count_data_demo.qmd`
- `biomarker_comparison.qmd`
- `forestsearch_scenario_validation.qmd` (with added `params$quick_mode`)
- `validation_glm_simulation_study.qmd`
- `validation_hte_tests_crump.qmd`

**Simulations (keep 5):**
- `actg175_binary_benefit_simulations.qmd` (with added `params$quick_mode`)
- `actg175_continuous_simulations.qmd`
- `actg175_survival_benefit_simulations.qmd`
- `gbsg_poisson_simulations.qmd`
- `reproduce_table1_m1.qmd` (with added `params$use_twostage` flag)

**Delete or archive (10):**
`gbsg_analysis_working.qmd`, `gbsg_analysis_working2.qmd`,
`gbsg_analysis_2sim-approaches.qmd`, `gbsg_analysis_cox_poisson_1.qmd`,
`forestsearch_scenario_validation_100.qmd`,
`applications/gbsg_poisson_simulations.qmd`,
`actg175_binary_benefit_simulations_1.qmd`,
`actg175_binary_benefit_simulations_100.qmd`,
`reproduce_table1_m1_1.qmd`, `reproduce_table1_m1_2stage.qmd`.

If you're cautious about deletion, move them to `applications/archive/`
and `simulations/archive/` instead. That preserves history without
cluttering the active directory.


## Should any be in the package README?

Probably two or three at most, and only as links for the curious reader.
The README is first-page real estate; it should describe what the
package does and how to use it, not enumerate development history.

A reasonable README-level listing would be:

```markdown
## Examples and simulations

In-depth walkthroughs live in the [`applications/`](applications/) and
[`simulations/`](simulations/) directories:

- [`gbsg_analysis.qmd`](applications/gbsg_analysis.qmd) — main vignette,
  exploratory subgroup identification on the German Breast Cancer
  Study Group dataset (Cox, time-to-event).
- [`gbsg_analysis_cox_poisson.qmd`](applications/gbsg_analysis_cox_poisson.qmd)
  — side-by-side Cox and Poisson rate models on the same data.
- [`count_data_demo.qmd`](applications/count_data_demo.qmd) — subgroup
  identification on recurrent-event count data (COPD exacerbations).
- [`actg175_binary_benefit_simulations.qmd`](simulations/actg175_binary_benefit_simulations.qmd)
  — full simulation study based on the ACTG175 HIV trial, showing
  operating characteristics under H0 and H1.
- [`reproduce_table1_m1.qmd`](simulations/reproduce_table1_m1.qmd) —
  reproduces the M1 block of León et al. (2024) Table 1 as a
  v0.2.0 pressure test.

For in-package smoke tests covering all four outcome types, see the
`quarto/` directory (`smoke_*.qmd`, `test_*.qmd`).
```

That's five bullets covering: main real-data analysis, GLM
demonstration, count-outcome tutorial, full simulation study, and a
canonical reproduction. A new user needs those five to understand the
package; the rest is interior plumbing.

Do not list the drafts, `_1`/`_100`/`_2stage` variants, or the methods
comparison (`biomarker_comparison.qmd`, `validation_*`) in the README.
They're valuable but niche; link from a pkgdown Articles index, not
the first page.


## Suggested follow-ups (not blocking)

1. **Delete or archive the 10 superseded files** listed above.
2. **Add `params$quick_mode` to the two `_100` files and remove the duplicates.**
3. **Add a one-line header comment** to each canonical file identifying
   it as canonical, so future readers don't have to diff.
4. **Add a directory-level `README.md`** (one paragraph each for
   `applications/` and `simulations/`) describing intent.
5. **Consider consolidating the three `gbsg_analysis_*` drafts**
   into a single canonical file (likely `gbsg_analysis.qmd` is
   already the target).
