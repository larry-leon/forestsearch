# CC task — vignette introduction (v4), post-selection vignette aligned to the manuscripts, the two methodology articles (2026-10-28)

**Repo:** `larry-leon/forestsearch` — run CC from this repo's root, on branch `feature/vignettes-0.4.1`.
**Kind:** documentation. Touches `vignettes/forestsearch.qmd` (prose: the introduction and references), `vignettes/survival-post-selection.qmd`, `vignettes/articles/methodology.Rmd`, `vignettes/articles/paper_simulations.Rmd`, `vignettes/articles/references.bib`, `NEWS.md`. **No edit to `R/`, `man/`, `DESCRIPTION`, `tests/`, `_pkgdown.yml`, `.Rbuildignore`, or any other article.** If a gate can only be passed by changing something under `R/`, STOP and report.
**Inputs from `~/Downloads`:** this file and `forestsearch_v4.qmd` (md5 `fe9925e92b598a0f913af53446032ca7`, 626 lines). Stage 0's record (`dev/reports/vignettes_stage0_2026-10-27/`) supplies the comparison values.
**Compute:** four renders — the primary vignette (74 s measured), the post-selection vignette (unmeasured at its new settings; its three MR products cost about what the primary vignette's FS fit costs), `methodology.Rmd` (diagrams and tables), `paper_simulations.Rmd` (runs 2 × 44 simulations under its own parallel plan; unmeasured here, known to complete on the pkgdown GitHub Action). Hard caps: 20 min each, 30 min for `paper_simulations`. No campaign.
**Unattended.** Gates stop the task; never ask, never work around. No `git fetch` / `pull` / `push`. Explicit paths on every `git add`; never stage a rendered `.html`. No `R CMD check`, no `devtools::check()`, no test suite, no `pkgdown::build_site()`, no `load_all()`.

Standing decisions (Larry, 2026-10-27/28): the vignette's introduction is v4 (approved 2026-10-28); `survival-post-selection` runs at the manuscript's settings (p\* 0.90, 5,000 draws); the two articles get the minimal treatment for 0.4.1 — a pointer note, citations, deprecated arguments removed, stale statements corrected — not a rewrite.

---

## 0. First action — GATE

1. `git status --short`: record pre-existing untracked paths (`vignettes/forestsearch.html` and `vignettes/articles/forestsearch-2024.html` are expected); never stage them. *GATE:* branch `feature/vignettes-0.4.1`; `git merge-base --is-ancestor 0de9a8e5 HEAD` succeeds. Record HEAD.
2. Copy this file to `dev/tasks/cc_task_vignettes_align_2026-10-28.md` (hyphens may arrive stripped; the committed name is this one); `git add` that path only; commit **`dev: task — vignettes align (2026-10-28)`**.
3. *GATE:* `~/Downloads/forestsearch_v4.qmd` exists with md5 `fe9925e92b598a0f913af53446032ca7`.
4. Provenance to `dev/reports/vignettes_align_2026-10-28/provenance.txt` (host, HEAD, R, installed forestsearch version and `Built`, quarto CLI, pkgdown). *GATE:* installed forestsearch `0.4.0.9000` with `NEWS.md` byte-identical to the repo's.

---

## 1. The primary vignette's introduction — commit 2

1. `cp ~/Downloads/forestsearch_v4.qmd vignettes/forestsearch.qmd`.
2. *GATE (prose only):* `knitr::purl()` the committed version (`git show HEAD:vignettes/forestsearch.qmd > /tmp/fs_head.qmd`) and the new file; the two purled scripts are byte-identical. `git diff` is confined to the block before `## The data` and to the `## References` section. No identity gate is re-run: the code has not changed.
3. Render under `/usr/bin/time -v`, cap 20 min: `quarto render vignettes/forestsearch.qmd`. *GATE:* exit 0; greps of `Error in`, `Warning`, `NaN`, `NA (` in the HTML all 0; the rendered Table 1 cells equal the committed render's (extract both tables' text with `pandoc -t plain` and `diff` the Table 1 block — identical). Record wall-clock and HTML size.
4. Commit `vignettes/forestsearch.qmd`: **`docs: vignette introduction — regulatory framing and the statistical pieces in bullets; references`**.

---

## 2. `vignettes/survival-post-selection.qmd` — commit 3

Apply these edits exactly (old → new; each old string occurs once in the file at HEAD). Nothing else in the file changes.

**A. Introduction** (lines 26–35):
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
→
```
`forestsearch()` identifies a subgroup by searching over candidate subgroups,
so a treatment effect estimated within the identified subgroup is affected by
the selection: a confidence interval computed as if the subgroup had been
prespecified tends to under-cover. With `mr_inference = TRUE`,
`forestsearch()` reports selection-adjusted estimates and bounds for the
identified subgroup and its complement, based on multiplier resampling of the
candidate effects. The constructions, their assumptions and their simulation
operating characteristics are described in León and Anderson (2026a); the
search that identifies the subgroup — the closed-form consistency screen, the
selection rule and the declaration calibration — in León and Anderson
(2026b). This vignette shows how to obtain and read the reported quantities on
a public dataset; `vignette("forestsearch")` runs the same analysis as the
manuscript's Table 1 and reads it alongside the two model-based identifiers.
```

**B. The settings sentence** (lines 57–59):
```
The fit below uses the consistency search without a LASSO/GRF/DINA screen,
`sg_focus = "effMaxSG"` (the largest subgroup within 20% of the strongest
observed effect), 200 consistency splits and 500 multiplier draws.
```
→
```
The fit below uses the consistency search without a LASSO/GRF/DINA screen,
`sg_focus = "effMaxSG"` (the largest subgroup within 20% of the strongest
observed effect), the closed-form consistency screen at p\* = 0.90 with every
screened candidate evaluated (`stop_threshold = NULL`), and 5,000 multiplier
draws — the settings of Table 1 of León and Anderson (2026b), so this fit is
the forest-search fit of `vignette("forestsearch")`.
```

**C. The call** (in chunk `fit`):
```
  mr_inference_args = list(draws = 500L, field_recovery = TRUE),
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
  hr.threshold = 1.0, hr.consistency = 1.0, pconsistency.threshold = 0.80,
  fs.splits = 200,
  details = FALSE, quiet = TRUE, plot.sg = FALSE,
```
→
```
  mr_inference_args = list(draws = 5000L, field_recovery = TRUE),
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
  hr.threshold = 1.0, hr.consistency = 1.0, pconsistency.threshold = 0.90,
  stop_threshold = NULL,
  details = FALSE, quiet = TRUE, plot.sg = FALSE,
```

**D.** (lines 116–118):
```
`ci_method = "ij"` only the two-sided interval is reported. León and Anderson
(2026) describe each construction and report its performance in simulation.
```
→
```
`ci_method = "ij"` only the two-sided interval is reported. León and Anderson
(2026a) describe each construction and report its performance in simulation.
```

**E.** (lines 131–134):
```
lower bound in the high 0.8s or above, against a naive subgroup estimate well
above 1, says the adjusted evidence does not establish clinically meaningful
harm -- the naive estimate was largely selection.
```
→
```
lower bound in the high 0.8s or above, against an unadjusted subgroup estimate
well above 1, says the adjusted evidence does not establish clinically
meaningful harm -- the unadjusted estimate was largely selection.
```

**F.** (chunk `locations`):
```
cat(sprintf("naive HR on Hhat            %.3f\n", g$naive$est))
cat(sprintf("de-biased HR on Hhat        %.3f\n", g$debiased$est))
```
→
```
cat(sprintf("unadjusted HR on Hhat       %.3f\n", g$naive$est))
cat(sprintf("corrected (two-term) HR     %.3f\n", g$debiased$est))
```

**G.** (lines 151–153):
```
The distance between the naive estimate and the de-biased one is what
selection cost; the spread across the constructions is what the adjustments
pay for it.
```
→
```
The distance between the unadjusted estimate and the corrected one is what
selection cost; the spread across the constructions is what the adjustments
pay for it.
```

**H.** (lines 174–175):
```
against the bounds. León and Anderson (2026) discuss how the bias correction
behaves across values of p-hat.
```
→
```
against the bounds. León and Anderson (2026a) discuss how the bias correction
behaves across values of p-hat.
```

**I. References** (lines 217–219):
```
León, L. F. and Anderson, K. M. (2026). Inference for standard trial
estimands after data-driven subgroup discovery. arXiv:2609.38361.
<https://arxiv.org/abs/2609.38361>
```
→
```
León, L. F. and Anderson, K. M. (2026a). Inference for standard trial
estimands after data-driven subgroup discovery. arXiv:2609.38361.
<https://doi.org/10.48550/arXiv.2609.38361>

León, L. F. and Anderson, K. M. (2026b). Interpretable subgroup identification
across survival and generalized-linear outcomes. arXiv:2610.28379.
<https://doi.org/10.48550/arXiv.2610.28379>
```

*GATE after the edits:* `grep -n 'Anderson (2026)' vignettes/survival-post-selection.qmd` returns nothing; `grep -n 'fs.splits\|0.80\|500L' vignettes/survival-post-selection.qmd` returns nothing.

**Render**, cap 20 min, under `/usr/bin/time -v`: `quarto render vignettes/survival-post-selection.qmd`. *GATE:* exit 0; HTML greps (`Error in`, `Warning`, `NaN`) 0. **Recorded, not gated:** the printed selected subgroup, the printed `unadjusted / corrected / field lower / field-s upper` lines and the `print(fs)` bounds, set beside stage 0's S3 values (`ladder_quantities.csv`): region `{er <= 0} & {pgr <= 26}`; 2.221839 / 1.276181 / 0.633561 / 0.803570; Bonferroni 0.508684 / 0.839594. The only difference between this call and the Table 1 call is `field_recovery = TRUE`; if any quantity differs from S3, say which and by how much — it tells us whether the recovery diagnostics consume draws. The p-hat and recovery values are recorded as rendered (they have no stage-0 counterpart).

Commit `vignettes/survival-post-selection.qmd`: **`docs: post-selection vignette at the manuscript's settings (p* 0.90, 5,000 draws, closed-form screen); cites 2026a/2026b`**.

---

## 3. `vignettes/articles/methodology.Rmd` and `references.bib` — commit 4

The article documents the 2024 methodology. Its two argument listings (fenced ```` ```r ```` blocks, display only — not evaluated) still show deprecated arguments and some wrong defaults.

1. **Pointer note.** Insert immediately before the line `## Introduction` (line 52), followed by a blank line:
   ```
   > **About this article (updated October 2026).** This article describes the methodology of León,
   > Jemielita, Guo, Marceau West and Anderson (2024), *Statistics in Medicine*,
   > <https://doi.org/10.1002/sim.10163>, as the package first implemented it: the split-sample consistency
   > screen evaluated by random half-splits, bootstrap bias correction, and GRF variable-importance screening
   > of the candidate factors. The current package evaluates the consistency screen in closed form
   > (`consistency_method = "resample"`, the default; random half-splitting is only the fallback, so
   > `fs.splits` has no effect otherwise), reports selection-adjusted estimates and one-sided bounds by
   > multiplier resampling (`mr_inference = TRUE`; León and Anderson 2026a, arXiv:2609.38361), and calibrates
   > the declaration of a subgroup at a protected level (León and Anderson 2026b, arXiv:2610.28379). The GRF
   > variable-importance screen is removed as of 0.4.1 (`vi.grf.min` is deprecated and has no effect). The
   > current methodology is documented in `vignette("forestsearch")` and `vignette("survival-post-selection")`;
   > the argument listings below are kept for the 2024 workflow and are not the complete current signature
   > of `forestsearch()` — see `?forestsearch`.
   ```
2. **Both listings** (the ```` ```r ```` blocks starting at lines 937 and 1048):
   - the header comment `# All forestsearch() parameters shown with their default values.` → `# forestsearch() arguments of the 2024 workflow with their defaults; not the complete current signature (see ?forestsearch).` — in both blocks (the second block's first comment line may differ in wording; replace whichever line claims that all parameters are shown).
   - delete the `max_n_confounders` line and the `vi.grf.min` line in each block (four deletions).
   - **reconcile the `# default …` comments** against `formals(forestsearch::forestsearch)`: for every line of the form `  <name> = <value>, # default <text>`, the text must deparse-equal the formal's default. Known mismatches to fix: `use_grf` (`# default TRUE` → `# default FALSE`, both blocks; the `<-- FOCAL` marker stays), `stop_threshold` (`# default 0.95` → `# default pconsistency.threshold`). Fix every other mismatch the check finds, and list each in the report (old comment → new). The argument *values* on the left are not changed, except that a line naming an argument that is no longer a formal of `forestsearch()` is deleted and reported (`workers` and `show_message` are elements of `parallel_args`, not formals — leave them).
   - *GATE:* a script in the record directory parses both blocks and asserts: every named argument (other than the two `parallel_args` elements) is a current formal; every `# default` comment equals the deparsed default; neither `vi.grf.min`, `max_n_confounders` nor `max.minutes` appears.
3. **Stale statements:**
   - line 1181: ``of consistency splits (`fs.splits`, default 400--1000), the number of`` → ``of consistency splits (`fs.splits`, default 1000; used only when the closed-form consistency evaluation is unavailable for a candidate), the number of``
   - lines 1172–1173: ``- **Minimum events:** 10--12 per treatment arm (default `d0.min = 12`,`` / ``  `d1.min = 12`)`` → ``- **Minimum events:** 10 per treatment arm (default `d0.min = 10`,`` / ``  `d1.min = 10`)``
4. **Citations.** `grep -n '2609.38361\|2610.28379' vignettes/articles/references.bib`. If absent, append:
   ```
   @misc{Leon_Anderson_2026a,
     title         = {Inference for standard trial estimands after data-driven subgroup discovery},
     author        = {Le{\'o}n, Larry F. and Anderson, Keaven M.},
     year          = {2026},
     eprint        = {2609.38361},
     archivePrefix = {arXiv},
     primaryClass  = {stat.ME},
     doi           = {10.48550/arXiv.2609.38361}
   }
   @misc{Leon_Anderson_2026b,
     title         = {Interpretable subgroup identification across survival and generalized-linear outcomes},
     author        = {Le{\'o}n, Larry F. and Anderson, Keaven M.},
     year          = {2026},
     eprint        = {2610.28379},
     archivePrefix = {arXiv},
     primaryClass  = {stat.ME},
     doi           = {10.48550/arXiv.2610.28379}
   }
   ```
   If present under other keys, use those keys. Then add the two keys to the YAML `nocite:` block after `@Wager_2014`.
5. **Render**, cap 20 min: `Rscript -e 'rmarkdown::render("vignettes/articles/methodology.Rmd")'` (set `RSTUDIO_PANDOC` as in task 1). *GATE:* exit 0; the rendered References section lists both 2026 entries (grep the HTML for `2609.38361` and `2610.28379`); the HTML stays untracked.
6. Commit `vignettes/articles/methodology.Rmd` and `vignettes/articles/references.bib`: **`docs: methodology article — 2024-workflow note, deprecated arguments removed, defaults reconciled, 2026a/2026b cited`**.

---

## 4. `vignettes/articles/paper_simulations.Rmd` — commit 5

1. **Pointer note.** Insert immediately before the line `# Introduction` (line 50), followed by a blank line:
   ```
   > **About this article (updated October 2026).** This article reproduces the simulation framework of
   > León, Jemielita, Guo, Marceau West and Anderson (2024), *Statistics in Medicine*,
   > <https://doi.org/10.1002/sim.10163>, and its published benchmarks, on the current package. The
   > simulation settings below are those of the 2024 workflow; two of them no longer act: the consistency
   > screen is evaluated in closed form, so `fs.splits` has no effect unless the closed form is unavailable
   > for a candidate, and the GRF variable-importance screen is removed as of 0.4.1 (`vi.grf.min` is
   > deprecated and has no effect). The operating characteristics of the current methodology — the
   > selection-adjusted bounds, the two-subgroup claim and the declaration calibration — are reported in León
   > and Anderson (2026a, arXiv:2609.38361) and (2026b, arXiv:2610.28379); the methodology is documented in
   > `vignette("forestsearch")` and `vignette("survival-post-selection")`.
   ```
2. **Deprecated argument in evaluated code.** In chunk `define-params`, delete the line `  vi.grf.min = -0.2,` (the only deprecated argument passed in evaluated code; `run_simulation_analysis()` forwards `fs_params` to `forestsearch()`, which would warn once per simulation). *GATE:* `grep -n 'vi.grf.min\|max_n_confounders\|max.minutes' vignettes/articles/paper_simulations.Rmd` returns only lines inside `eval=FALSE` chunks or prose — list them in the report; none in an evaluated chunk.
3. **Render**, cap 30 min, under `/usr/bin/time -v`: `Rscript -e 'rmarkdown::render("vignettes/articles/paper_simulations.Rmd")'`. The article sets its own parallel plan in its `setup-parallel` chunk — do not change it. *GATE:* exit 0. Record wall-clock, and the `n_sims` the article ran (44 + 44). The HTML stays untracked. If the cap is hit: record "exceeded cap", keep the edits (they do not depend on the render), and say so in the report — the GitHub Action is where this article is known to build.
4. Commit `vignettes/articles/paper_simulations.Rmd`: **`docs: simulation article — 2024-workflow note, deprecated vi.grf.min removed from the evaluated settings`**.

---

## 5. `NEWS.md` — commit 6

Under `### Documentation` in the development section, append:
```
* `vignette("survival-post-selection")` now runs at the manuscript's settings (consistency level 0.90,
  5,000 multiplier draws, the closed-form consistency screen) and cites both 2026 papers; the
  methodology and simulation articles on the package website carry a note placing them as the 2024
  workflow, with the deprecated arguments removed from their listings.
```
Commit `NEWS.md`: **`docs: NEWS — post-selection vignette settings and the article notes`**.

---

## 6. Report — `dev/reports/REPORT_vignettes_align_2026-10-28.md`, record `dev/reports/vignettes_align_2026-10-28/` — commit 7

Bullets, one item each: provenance; §1 purl-identity and render (wall-clock, size, Table 1 block identical); §2 the recorded post-selection quantities beside S3, with the `field_recovery` reading; §3 the reconciliation list (every `# default` comment changed, every line deleted), the bib outcome, render; §4 the grep listing and render wall-clock (or cap); deviations verbatim; OPEN ITEMS (findings only). Commit the record directory (scripts, logs, render logs, `provenance.txt`) and the report: **`dev: report — vignettes align (2026-10-28)`**. No push.

Copy to `~/Downloads/vignettes_align_2026-10-28/`: the report; the four rendered HTMLs (`forestsearch.html`, `survival-post-selection.html`, `methodology.html`, `paper_simulations.html`); `vignettes/survival-post-selection.qmd`, `vignettes/articles/methodology.Rmd`, `vignettes/articles/paper_simulations.Rmd` as committed; `git diff 0de9a8e5 -- vignettes/articles/references.bib NEWS.md > diffs.txt`. `ls -la` into the report's last section.

---

## 7. Out of scope

- No edit to `R/`, `man/`, `DESCRIPTION`, `tests/`, `_pkgdown.yml`, `.Rbuildignore`, `README.md`; no edit to `forestsearch-2024.Rmd`, `extreme_subgroups.qmd` or the three treatment-effect articles; no edit to the code chunks of `vignettes/forestsearch.qmd`.
- No `R CMD check` / `devtools::check()`, no test suite, no `pkgdown::build_site()`, no `devtools::install()`; no re-run of the stage-0 ladder or the task-1 identity gate; no push.
