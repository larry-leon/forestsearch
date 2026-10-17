# TASK: Extreme-subgroups article — Stage 2 prose edits (D1–D6, X1)

- Date: 2026-10-17
- Repo: `larry-leon/forestsearch`, branch `master`. The starting HEAD must be `ef9c829a` or a descendant of it.
- File changed: **only** `vignettes/articles/extreme_subgroups.Rmd` (a pkgdown-only article).
- Nature: prose and display-block edits, plus one knob change (`n_sims_null` 144 → 150). **No `R/` changes. No other R chunk code changes.**
- Report: `~/Downloads/REPORT_extreme_subgroups_article_edits_2026-10-17.md`

## Standing rules

- Do not run `git fetch`, `git pull` or `git push`. Every `git add` names explicit paths. Never stage pre-existing untracked files.
- Do NOT run `R CMD check`, `rcmdcheck`, `devtools::check()`, `build_vignettes()`, any `pkgdown::build_*()`, `rmarkdown::render()`, `quarto render` or testthat. Do NOT evaluate any article chunk.
- Gates are stop-on-failure, not stop-to-ask. On failure, do not commit the article. Write the report with a `FAILED GATE` section and stop.
- Hard time cap: each step 5 minutes; the whole task 30 minutes.

## Step 0 — Record the task

1. Copy this file from `~/Downloads/` to `dev/tasks/TASK_extreme_subgroups_article_edits_2026-10-17.md`.
2. `git add` that path.
3. Commit with the message `dev/tasks: add extreme-subgroups article edits task (2026-10-17)`.

## Step 1 — Baseline snapshot (outside the repo)

- Copy the article to a temporary file outside the repo (for example, under `$TMPDIR`) as `before.Rmd`.
- With base R (`Rscript --vanilla`), extract every `{r ...}` chunk (label + body) from `before.Rmd` and save the result for the Step 4 comparison.

## Step 2 — Edit mechanics

Apply the edits with a short Python or R script kept **outside the repo**, in the order listed.

- **Anchors.** An anchor matches a line if the line, with trailing whitespace stripped, *contains* the anchor string.
- **Uniqueness gate.** At the moment each edit is applied, its start and end anchors must each match exactly one line, and the end must come at or after the start. Otherwise, FAILED GATE.
- **REPLACE.** Replace the lines from the start line through the end line, inclusive, with the given block.
- **INSERT AFTER.** Insert the given block after the anchor line.
- **INSERT BEFORE.** Insert the given block before the anchor line.
- Blocks are given between `~~~~` fences. Insert the content verbatim, without the fences, and follow each block with one blank line.

### E2 (D1, D2, D3, D5 sentence) — replace the front matter

REPLACE. Start: `# Overview`. End: the line **immediately before** the line `# Setup`. The `# Setup` line itself is kept.

~~~~
# Overview

A central concern in exploratory subgroup analysis is the potential for
**extreme results driven entirely by chance**, not by genuine treatment
effect heterogeneity. When many subgroups are examined in a single trial,
small subgroups in particular can produce hazard ratio estimates and
confidence bounds that look dramatic simply because they are noisy.

This article evaluates that risk directly. Synthetic trials are generated
from a data generating mechanism (DGM) with a **uniform treatment
effect** (the true hazard ratio is the same for every patient), and a Cox
model is fitted in each of a set of pre-defined subgroups. Any
subgroup-to-subgroup difference in the results is therefore due to
chance.

The article is a self-contained workflow using the GBSG breast cancer
trial (`survival::gbsg`):

1. **Building a DGM** (Part I): fitting a Weibull AFT model to the GBSG
   data with `generate_aft_dgm_flex()` and constructing a super population.
2. **Censoring diagnostics and calibration** (Part II): checking that
   simulated trials reproduce the censoring pattern of the source data
   with `check_censoring_dgm()`, and calibrating it with
   `calibrate_cens_adjust()`.
3. **Generating synthetic trials** (Part III): drawing trials from the
   DGM with `simulate_from_dgm()`.
4. **Extreme subgroup effects under a uniform effect** (Part IV): building
   a uniform-effect DGM (`model = "null"`), simulating `r n_sims_null`
   trials, and fitting a Cox model stratified by `grade` in each
   pre-defined subgroup. Random benchmark subgroups (`random15`,
   `random20`, `random40`, `random60`) with no clinical meaning are
   included alongside the clinical subgroups, so that subgroups of the
   same size can be compared directly.

Parts I–III illustrate the DGM tools on a DGM *with* an embedded subgroup
effect (ER-low and pre-menopausal). Part IV builds the uniform-effect DGM
used for the main results. For each subgroup, Part IV summarises the
distribution across simulated trials of the HR point estimate and of the
upper 95% confidence bound UB(HR).

The **forestsearch** methodology is described in {{SIM2024_INTEXT}};
inference for trial estimands after data-driven subgroup discovery is
developed in León and Anderson (2026). See the References section.

## Installation

**forestsearch** is available on CRAN:

```r
install.packages("forestsearch")
```

The development version can be installed from GitHub:

```r
# install.packages("pak")
pak::pak("larry-leon/forestsearch")
```

{{WEIGHTEDSURV_SENTENCE}}

```r
library(forestsearch)
```

## Documentation

Full documentation, including the function reference and all articles,
is available on the pkgdown site:
<https://larry-leon.github.io/forestsearch/>.

The Getting Started vignette is also installed with the package:

```r
vignette("forestsearch")
```

The other articles, including this one, are available on the pkgdown
site only.

---
~~~~

Fill in the placeholders before writing the file:

- `{{SIM2024_INTEXT}}`: look up the entry in `vignettes/articles/references.bib` whose `doi` is `10.1002/sim.10163` (case-insensitive).
  - 3 or more authors: `León et al. (2024)`.
  - 2 authors: `<Surname1> and <Surname2> (2024)`.
  - Entry not found: `León et al. (2024)`. Record that it was not found in the report.
- `{{WEIGHTEDSURV_SENTENCE}}`: read `DESCRIPTION`.
  - `weightedsurv` appears in `Depends` or `Imports`: `Dependencies, including [**weightedsurv**](https://CRAN.R-project.org/package=weightedsurv), are installed automatically.`
  - Otherwise: delete the placeholder line (and its surrounding blank line, so only one blank line remains).
  - Record the matched field, or "absent", in the report.

### E1 (D6) — setup chunk

REPLACE. Start: `# Number of null-DGM simulated trials`. End: `n_sims_null <- 144L`.

~~~~
# Number of null-DGM simulated trials. 150 keeps this article quick to build.
# For an actual application use 5,000-10,000, after timing a pilot run of
# about 100 trials on your own machine (see "Number of simulated trials").
n_sims_null <- 150L
~~~~

### E3 (D6) — note on the number of simulated trials

INSERT AFTER. Anchor: `## Simulation Loop`. Put one blank line before the block.

~~~~
**Number of simulated trials.** This article runs `r n_sims_null`
simulated trials so that it builds quickly. That is enough to show the
qualitative pattern, but the tail summaries (the 1st and 99th
percentiles, and the distribution of UB(HR) given UB ≥ 1) are imprecise
at this size. For an actual application, run 5,000–10,000 trials. Before
committing to a long run, time a pilot run of about 100 trials on your
own machine (for example by wrapping the loop in `system.time()`) and
scale up from there.
~~~~

### X1 — Part III prose matches the DGM it uses

REPLACE. Start: `With the null DGM built, diagnosed, and calibrated`. End: `looping over subgroups inside each replicate.`

~~~~
With the DGM built, diagnosed, and calibrated, each simulated trial is a
single call to `simulate_from_dgm()`. The function samples from the
super population with replacement, assigns treatment, generates latent
survival times from the AFT model, applies the calibrated censoring
model, and truncates at `analysis_time`.

The single-trial example below uses the Part I DGM, which has an
embedded subgroup effect, so `flag_harm` marks that subgroup (ER-low and
pre-menopausal). The `r n_sims_null` replications that drive Part IV
follow the same pattern, but are drawn from the uniform-effect DGM and
loop over subgroups inside each replicate.
~~~~

### E4a (D4) — Part IV "Core Question" bullet

REPLACE. Start: `- The **upper confidence bound UB(HR)** — the upper end of the Cox 95%`. End: `"failure to show benefit" looks in the worst-case trial).`

~~~~
- The **upper confidence bound UB(HR)**, the upper end of the Cox 95%
  CI, summarised by its distribution across trials, together with the
  proportion of trials with UB(HR) below 1.0 and the distribution of
  UB(HR) among trials where it is at or above 1.0.
~~~~

### E4b (D4) — "Summary Statistics" UB paragraph

REPLACE. Start: `**Upper confidence bound UB(HR).** The upper end of the Cox 95%`. End: `what a safety reviewer sees in the worst-case trial for that subgroup.`

~~~~
**Upper confidence bound UB(HR).** The upper end of the Cox 95%
confidence interval. Its distribution across trials is summarised in the
same way, together with two reference summaries:

- **`Pr(UB < 1.0)`**: the proportion of trials in which the upper bound
  lies below 1.0. HR = 1.0 is a reference line, not a decision rule:
  where an upper bound sits relative to clinically meaningful effect
  sizes matters more than which side of 1.0 it falls on.
- **Distribution of UB(HR) given UB ≥ 1**: in trials where the upper
  bound is at or above 1.0, how far above 1.0 it lies.
~~~~

### E4c (D4) — "Reading the table"

REPLACE. Start: `**Reading the table.**`. End: `` `random40` it commonly exceeds 3.0. ``

~~~~
**Reading the table.** Every row shares the same true HR, so the median
of the HR ECI should be close to the ITT value throughout. The *spread*
of the ECI is what changes: narrow for the ITT and the large clinical
subgroups, wide for `random20` and `random15`. The UB(HR) columns show
the same pattern for the upper bound: as subgroup size falls, upper
bounds far above the true HR become common. `P99 UB | UB>=1` shows how
large the upper bound is in the worst 1% of those trials where it is at
or above 1.0.
~~~~

### E4d (D4) — UB single-panel caption

Substring replacement inside the line containing `{r ub-forest-single`. The substring must occur exactly once in the file.

- Old: `Values right of HR = 1.0 (dotted line) represent trials failing to show significant benefit.`
- New: `Values right of the dotted line are upper bounds at or above HR = 1.0.`

### E4e (D4) — paragraph before the UB panels

REPLACE. Start: `The upper confidence bound UB(HR) is the quantity a safety reviewer`. End: ``benefit (`Med UB | UB>=1`).``

~~~~
The upper confidence bound UB(HR) is the quantity a reviewer examines
when asking whether a subgroup shows harm or lack of benefit. The plots
below show the distribution of UB(HR) across simulations in the same
ECI format, with HR = 1.0 marked as a reference (grey dotted line). The
right-hand columns show the proportion of trials with UB(HR) at or above
1.0 (`UB>=1`) and the median UB(HR) among those trials (`mUB | UB>=1`).
~~~~

### E4f (D4) — UB Panel B lead-in

REPLACE. Start: ``The `Med UB | UB>=1` column is the key safety quantity``. End: `to demonstrate benefit.`

~~~~
The `mUB | UB>=1` column shows how large the upper bound is in the
trials where it is at or above 1.0.
~~~~

### E5 (D5) — References section

INSERT BEFORE. Anchor: `# Session Information`.

~~~~
# References

{{SIM2024_FULL}}

León, L. F. and Anderson, K. M. (2026). Inference for standard trial
estimands after data-driven subgroup discovery. arXiv:2609.38361
[stat.ME]. <https://arxiv.org/abs/2609.38361>

---
~~~~

`{{SIM2024_FULL}}` comes from the same `references.bib` entry, formatted as plain text:

- Format: `<Authors> (<year>). <Title>. *<journal>*, <volume>(<number>), <pages>. <https://doi.org/10.1002/sim.10163>`
- Omit any field the entry lacks.
- Authors: `Surname, I.` separated by commas, with `and` before the last.
- Entry not found: `León, L. F. et al. (2024). *Statistics in Medicine*. <https://doi.org/10.1002/sim.10163>`

Do not add a `bibliography:` field and do not use `[@key]` syntax.

## Step 3 — Placeholder check

Gate: the edited article contains no `{{`.

## Step 4 — Post-conditions (machine-checkable; record each as PASS/FAIL with the counts)

Count lines in the edited article with fixed-string matching.

**Expected 0:**
- `get_dgm_stratified`
- `BRAF`
- `Six Analysis Variants`
- `remotes::install_github`
- `vignette("methodology")`
- `144`
- `Med UB | UB>=1`
- `significan` (case-insensitive). If this is nonzero, list the lines and FAIL.
- `failure to show benefit`
- `16 topic areas`

**Expected exactly 1:**
- `n_sims_null <- 150L`
- `install.packages("forestsearch")`
- `# Overview`
- `## Setting: Uniform Treatment Effect`
- `## Subgroups Under Study`
- `**Number of simulated trials.**`
- `# References`

**Expected ≥ 1:**
- `2609.38361`
- `10.1002/sim.10163`

**Chunk invariance.** Extract the `{r ...}` chunks from the edited file as in Step 1.
- The same labels must appear, in the same order.
- Every chunk body must be byte-identical to the baseline, except `setup`.
- In `setup`, the only differing lines may be the three comment lines and the `n_sims_null` line from E1.
- All chunks must `parse()` without error.

**YAML.** The front-matter block (between the first two `---` lines) is byte-identical to the baseline.

## Step 5 — Commit

If every gate and post-condition passes:

1. `git add vignettes/articles/extreme_subgroups.Rmd`
2. Commit with this message:
   ```
   articles/extreme_subgroups: update for CRAN 0.4.0 (front matter, install/docs, UB framing, references, n_sims 150)
   ```

Then check:
- `git show --stat HEAD` lists only that file.
- `git status --porcelain` shows no modified tracked files.
- There are exactly two new commits since the starting HEAD: the task doc and the article.

## Report (`~/Downloads/REPORT_extreme_subgroups_article_edits_2026-10-17.md`)

1. Starting HEAD, and the SHAs of the two new commits.
2. The placeholder resolutions: the bib entry found or not, the authors used, and the weightedsurv field.
3. The post-condition table (PASS/FAIL with counts).
4. The full `git diff <starting HEAD> HEAD -- vignettes/articles/extreme_subgroups.Rmd`, verbatim.
5. Any FAILED GATE detail.

In the terminal, print one line only: the report path and overall PASS/FAIL.
