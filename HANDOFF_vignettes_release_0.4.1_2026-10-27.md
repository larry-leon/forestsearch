# HANDOFF — vignette/pkgdown update and the 0.4.1 release (from the chat of 2026-10-27)

**For:** the next forestsearch chat(s): (1) the vignette and pkgdown-article update, (2) the CRAN 0.4.1 submission. **Written:** 2026-10-27, ~15:10 PT. **State verified against GitHub read-only at writing time.** Everything else is from the committed reports named below — re-verify from source before relying on any number or line reference; the project KB is stale for R code.

---

## 0. Read this first

- **Target: submit forestsearch 0.4.1 to CRAN on Monday 10 November 2026** (Larry's decision, 2026-10-27; reverses the earlier "December batch"). Four weeks after 0.4.0's acceptance (2026-10-10); the justification for a sub-two-month update is the user-visible bug fix for GitHub #7 (`xlog` never worked) and the `forestploter (>= 1.2.0)` requirement — say so in `cran-comments.md`. Version number: **0.4.1 proposed** (no new API: bug fix, two deprecations, performance, documentation); Larry has not yet confirmed 0.4.1 vs 0.5.0 — ask once.
- **The release should carry the updated vignettes.** The window before 10 Nov is for that work. Everything else for the release is a half-day task (§5).
- **What's already public:** `master` on GitHub and the pkgdown site (deployed automatically by `.github/workflows/pkgdown.yaml` on every push to `master`; verified live 2026-10-27 15:00 PT: `reference/forestsearch.html` carries both arXiv IDs, the three manuscript-aligned examples and the "Options not tested" section; NEWS shows the `max.minutes` and `fit_causal_forest_glm()` entries). CRAN is only for `install.packages()` users and the CRAN landing page.
- **Primary methodology references (verbatim, suffixes by posting date — Larry's decision):**
  - **2026a:** León, L. F. and Anderson, K. M. (2026a). *Inference for standard trial estimands after data-driven subgroup discovery.* arXiv:2609.38361 (29 Sep 2026). doi:10.48550/arXiv.2609.38361. — post-selection inference: MR/IJ, field bounds, complement, two-subgroup claim; survival (GBSG) application in its Sec. 5.
  - **2026b:** León, L. F. and Anderson, K. M. (2026b). *Interpretable subgroup identification across survival and generalized-linear outcomes.* arXiv:2610.28379 (20 Oct 2026). doi:10.48550/arXiv.2610.28379. — the search itself: forest search as closed-form stability selection on splitting consistency, family-wise false-declaration rate, two-subgroup claim; survival, binary and continuous applications (GBSG, ACTG175) and simulations built from them; supplement S3 has the ACTG175 illustrations in full, S5 the identifiers' settings.
  - Both PDFs are in the project KB as `LeonAnderson_2026_InferenceStandardTrialEstimands.pdf` and `LeonAnderson_2026_InterpretableSubgroupIdentificationSurvival.pdf` (the second file is the paper + supplement). The KB is not on CC's machine: anything CC must read travels `~/Downloads` → `dev/tasks/`.
- Protocol unchanged: three-layer (Larry decides → chat writes task docs → CC executes; CC never pushes/fetches/merges; Larry pushes via GitHub Desktop). Kickoffs short, in a fenced code block, "clear CC first, then paste". Task docs committed first and alone. Gates stop-on-failure. No `R CMD check`, vignette build or full suite inside a task (vignette *rendering* of the file being edited is fine as a gate). Never overwrite legacy analysis files under `quarto/applications/**`. Branch for anything that touches `DESCRIPTION`, `man/` or `vignettes/`; Larry's `devtools::check()` on the branch, then merge.

---

## 1. State at handoff (verified)

| item | state |
|---|---|
| forestsearch `master` | `30441802` ("dev: report — pre-check tidy"); pushed 2026-10-27 ~14:26 PT; fast-forward of `feature/docs-refs-maxminutes` |
| Remote branches | **`master` and `gh-pages` only** — Larry deleted 18 stale branches today; nothing unmerged was lost (checked per branch) |
| Local repo | `master` only; no worktrees (two leftover worktrees under `.claude/worktrees/` removed; see §6 archive) |
| Version | `0.4.0.9000`; CRAN carries 0.4.0 (accepted 2026-10-10) |
| `devtools::check()` on `30441802`, `--as-cran`, NOT_CRAN | **0 errors / 0 warnings / 0 notes**, 11m 53s; examples OK (incl. `--run-donttest`), tests 403 s, both vignettes rebuilt in 45 s |
| pkgdown site | live from `30441802` (`gh-pages` = `9c313e1c`) |
| Installed build on Pop!_OS | `30441802` tree via `devtools::install()` (CC, during the pre-check tidy); vignettes not installed |
| fs-glms-interpretable `main` | `a772a2fe` — untouched today; the roadmap placement (prior handoff §3.9) is still unconfirmed |
| Project KB | full; nothing filed today (the task docs and reports are in the repo) |

Commits on `master` from today, in order (all CC; merged/pushed by Larry):

| SHA | what |
|---|---|
| `a1b07bc9` | task doc `dev/tasks/cc_task_docs_refs_maxminutes_2026-10-27.md` |
| `8758718a` | refs: 2026a/2026b in `forestsearch()` roxygen (Description + `@references`), `DESCRIPTION`, README BibTeX |
| `f35399a5` | deprecate `max.minutes` (both formals); remove it from every package call site; NEWS |
| `a2d1c70a` | "not tested" sentence on `adjust_covariates`/`ps_method`/`ps_adjust_method`/`ps_hat` (16 man pages) + `@section` in `?forestsearch` |
| `9fd492d4` | `@examples` replaced by three manuscript-aligned examples |
| `2dfeea63` | `DESCRIPTION`: `Suggests: speff2trial` |
| `885e0be2` | man regenerated; `tests/testthat/test-max-minutes-deprecation.R` |
| `b36da9dc` | report `dev/reports/REPORT_docs_refs_maxminutes_2026-10-27.md` |
| `725d64b2`, `c5dfe0b1` | task doc `dev/tasks/cc_task_prechec_tidy_2026-10-27.md` (committed blank, then with Larry's A/B/C dispositions) |
| `8f94d117` | `Config/roxygen2/version: 8.1.1` accepted |
| `81d3dab3` | `fs_family_report()` `max.minutes` row note updated |
| `936b6734` | `fit_causal_forest_glm()` deleted (dead since the `vi.grf.min` deprecation); NEWS "Internal" |
| `de298322` | man regenerated (`fit_causal_forest_glm.Rd` deleted) |
| `30441802` | report: pre-check tidy section appended |

---

## 2. What the package now does differently vs CRAN 0.4.0 (cumulative; the 2026-10-26 handoff §2 covers the first three)

1. `vi.grf.min` / `max_n_confounders` deprecated (inert, warn on non-default). 2. Lean search-stage fits A (survival `coxph.fit`) and F (binary OR `glm.fit`), bit-identical, −44 % / −54 % per bootstrap replicate. 3. `xlog` works; `forestploter (>= 1.2.0)` API (#7).
4. **`max.minutes` deprecated** in `forestsearch()` and `subgroup.search()`. Mechanics (deviation from the task doc, approved mid-task; recorded in the report): `forestsearch()` warns on a **non-default value** (`!missing(max.minutes) && !identical(as.numeric(max.minutes), 3)`), because `args_call_all <- mget(names(formals()))` always carries it and bootstrap/CV/`fs_fdr_report`/MRCT/calibrate-null replay that list — a bare `!missing()` would have warned once per replicate. `forestsearch()` no longer forwards it (`search_args$max.minutes <- NULL`, `forestsearch_main.R` ≈L3477). `subgroup.search()` warns on `!missing()`. Removed from `vignettes/forestsearch.Rmd`, `vignettes/survival-post-selection.qmd`, `vignettes/articles/methodology.Rmd` (×2), `run_simulation_analysis()` defaults. The legacy `quarto/applications/gbsg/analysis_gbsg_survival_effMaxSG.qmd` (L580) still passes it — left by rule; it would warn (non-default? it passes 3 → no warning from `forestsearch()`, so silent).
5. **References**: `?forestsearch` Description now says inference follows 2026a and the search methodology is in 2026b; `@references` lists 2026a, 2026b first; `DESCRIPTION` cites 2024, 2026a, 2026b; `print.forestsearch` cites 2026a only (MR page). README has `@misc` BibTeX for both (with `doi` fields — kept, Larry's decision).
6. **"Not tested" note**, verbatim first sentence on the four `@param` blocks across 16 man pages: *This option is implemented but has not been tested against the manuscript analyses or the simulation studies of León and Anderson (2026a, 2026b), which use the unadjusted treatment-only fits.* Plus `@section Options not tested against the manuscript analyses:` in `?forestsearch`. No "Not validated" label (Larry). If the sentence reads as noise on helper pages (`get_split_hr_fast`, `build_cox_formula`, …), trim to a cross-reference in this workstream.
7. **Examples** (`\dontrun{}`): three calls mirroring the manuscript analyses, every argument checked against the lean-build re-renders' resolved-parameter tables: (1) GBSG survival, effMaxSG/neighborhood 0.10, HR thresholds 1.0/1.0, p* 0.90, `conf_force = c("er <= 0", "pgr <= 0")`, `fs.splits = 1000`, MR `draws = 5000L`, `seedit = 8316951`; (2) ACTG175 binary OR (`y_binary_neg`, adverse), thresholds 1.0/1.0/0.90, `conf.cont_jcuts = list(cd40 = 10, wtkg = 10)`, `fs.splits = 500`, MR draws 5000; (3) ACTG175 continuous MD (`y_decline`), thresholds 10/5/0.90, ε = 0.20, six `conf.cont_jcuts`, field bounds with complement, `field_scale_complement = "selected"`, `plan = "sequential"`. Examples 2–3 guarded by `requireNamespace("speff2trial")`.
8. `Config/roxygen2/version` is 8.1.1 — **the restore-after-`document()` step is gone**. `fit_causal_forest_glm()` removed (the GLM GRF path in `grf.subg.harm.glm()` has always fitted its forest inline, `R/grf_subg_harm_glm.R` ≈L460).

---

## 3. The vignette / pkgdown workstream (the next chat)

### 3.1 What exists (from source at `30441802`)

Two kinds of document share `vignettes/`, split by `.Rbuildignore` (`^vignettes/articles$`):

| kind | files | built by `R CMD check`? | published by |
|---|---|---|---|
| **CRAN vignettes (2)** | `vignettes/forestsearch.Rmd` (924 lines, "Exploratory Subgroup Identification with ForestSearch", GBSG; `VignetteBuilder: knitr`); `vignettes/survival-post-selection.qmd` (228 lines, GBSG, `quarto`) | yes (45 s today) | CRAN + pkgdown |
| **pkgdown articles (6)** | `vignettes/articles/{methodology, treatment_effect_definitions, biomarker_effects, causal_effects_brief_review, paper_simulations, extreme_subgroups}` | no | pkgdown only |

`_pkgdown.yml` groups them: Getting Started (`forestsearch`), Methodology (`articles/methodology`), Post-Selection Inference (`survival-post-selection`), Treatment Effect Definitions (three articles), Simulation Studies (`paper_simulations`, `extreme_subgroups`).

What the two vignettes currently are:
- `forestsearch.Rmd`: pre-manuscript "every argument with its default" walkthrough — GBSG prep, KM, a `grf.subg.harm.survival()` preliminary, then `forestsearch()` with all ~86 formals listed (≈L290–395), results tables, full bootstrap (`forestsearch_bootstrap_dofuture`), KM by subgroup, 10-fold CV. Cites only León et al. (2024). It is the document CRAN users read first and the one furthest from the manuscripts.
- `survival-post-selection.qmd`: already manuscript-aligned to 2026a (GBSG fit with `mr_inference = TRUE`, `draws = 500L`, `fs.splits = 200`, reading bounds by location, p̂ as a diagnostic). Cites 2026a. Needs 2026b added and a check that its wording matches the final arXiv text.

### 3.2 Scope Larry asked for (2026-10-27)

"Incorporate the two manuscripts and expand on the examples (the survival, continuous and binary analyses)." Concretely, as a proposal for the next chat to put to Larry before any task doc:

1. **`forestsearch.Rmd` (getting started):** re-centre on the manuscript GBSG analysis (Example 1 of `?forestsearch`): the search as 2026b describes it (cut grid, splitting consistency screen, the three thresholds, `sg_focus`/selection rule), then MR inference as 2026a describes it; keep bootstrap and CV sections but shorten, and drop the all-86-formals listing in favour of the manuscript call with the handful of non-default arguments explained. Cite 2024, 2026a, 2026b.
2. **A new GLM vignette** (binary OR and continuous MD on ACTG175, i.e. Examples 2–3): `speff2trial` is already in Suggests; `cont.cutoff` handles the 0/1 covariates numerically for binary (the fixed-family document's note); the continuous one uses factors. This is the "expand on the examples" item — nothing on CRAN currently shows a GLM outcome end-to-end.
3. **`survival-post-selection.qmd`:** add 2026b; align terminology with the final arXiv text (field bounds, complement, two-subgroup claim).
4. **Articles:** `methodology` and `paper_simulations` are the two that cite results — bring their citations to 2026a/2026b and check nothing contradicts the papers; `extreme_subgroups` was updated 2026-10-17/18 (see the pkgdown-articles-update memory). The other three (treatment-effect definitions) are reference material; leave unless wrong.
5. **Trim the 16-page "not tested" sentence** to a cross-reference on helper pages if Larry finds it noisy once seen on the site.

Constraints to design around:
- **CRAN build time.** Vignettes must build in minutes on CRAN's machines. `survival-post-selection.qmd` already runs reduced settings (`fs.splits = 200`, `draws = 500L`); a GLM vignette needs the same (the binary fixed-family run took 394 s *per render* at full settings on 102 cores — not CRAN-viable). Either reduced settings with a sentence saying the manuscript values, or precomputed results shipped under `inst/` and loaded (`R CMD check` then only renders text). Decide with Larry; reduced settings keep the vignette honest and self-contained, which is his stated preference for analysis documents.
- **Parallelism:** CRAN limits to 2 cores; `forestsearch.Rmd` already caps at 2 (L263).
- **Self-containment:** no reads from `dev/` scratch; wrapper calls visible and few; interpretation inside the vignette, not a separate memo.
- **Branch + check:** `feature/vignettes-0.4.1` (or similar); Larry's `devtools::check()` before merge; the pkgdown Action publishes on push.
- **Sources of truth for the analyses:** `quarto/applications/gbsg/analysis_gbsg_survival_effMaxSG.qmd`, `quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd`, `quarto/applications/actg175/analysis_actg175_continuous_intervals.qmd` (and their `rerun_2026-10-26/*_lean/` renders, which are the current numbers — Larry has the HTMLs). Numbers quoted in a vignette must come from the vignette's own run, not be pasted from those documents.

### 3.3 Reading list for the next chat

Both papers from the KB (abstract + application sections + supplement S3/S5 of 2026b); `vignettes/forestsearch.Rmd` and `survival-post-selection.qmd` from source; `?forestsearch` as rendered on the site (`reference/forestsearch.html`); `_pkgdown.yml`; `REPORT_docs_refs_maxminutes_2026-10-27.md` (for what moved today).

---

## 4. Standing decisions from today (binding)

- Release 0.4.1 targeted **10 Nov 2026**; version 0.4.1 proposed, Larry to confirm. The deprecated formals (`vi.grf.min`, `max_n_confounders`, `max.minutes`) stay through 0.4.1 and are removed in the release after.
- 2026a/2026b suffixes by posting date; README BibTeX keeps `doi`; no "Not validated" label anywhere.
- `max.minutes`: value trigger in `forestsearch()`, `!missing()` in `subgroup.search()`, never forwarded — the pattern for any future deprecation of a formal that `args_call_all` carries.
- Roxygen stamp 8.1.1 accepted: `document()` no longer needs a `DESCRIPTION` restore.
- Never overwrite legacy analysis files; `quarto/applications/**` is not edited by tasks, only recorded.
- Housekeeping precedent: before deleting a branch, check `git worktree list` (a branch checked out in a leftover `.claude/worktrees/*` worktree cannot be deleted, and `git worktree remove --force` destroys uncommitted files there — archive first). `gh-pages` is the site; never delete.

---

## 5. Release checklist for 0.4.1 (one CC task + Larry's checks; ~half a day)

1. Branch `release/0.4.1` from `master`. Task doc.
2. `DESCRIPTION`: `Version: 0.4.1`; `Date:` if present. `NEWS.md`: `# forestsearch (development version)` → `# forestsearch 0.4.1` with the date; re-read the Deprecated / Performance / Internal / Bug fix sections for a CRAN reader (they were written incrementally).
3. `cran-comments.md`: what changed since 0.4.0 (bug fix #7 + `forestploter >= 1.2.0`; two deprecations with warnings; performance; documentation incl. the new vignette); why 4 weeks after 0.4.0; test environments (local Pop!_OS R 4.6.1; win-builder devel/release; R-hub if used); `0 errors | 0 warnings | 0 notes`; reverse dependencies: none on CRAN (confirm with `tools::package_dependencies("forestsearch", reverse = TRUE)` or `revdepcheck` — forestploter is a *forward* dependency, not a revdep).
4. `devtools::document()`; `devtools::check(remote = TRUE, manual = TRUE)`; `devtools::check_win_devel()` + `check_win_release()` (results by e-mail, ~30 min); optionally `rhub::rhub_check()`.
5. `devtools::spell_check()`; `urlchecker::url_check()` (the two arXiv DOIs and the pkgdown URL).
6. Merge, push, `devtools::release()` or `devtools::submit_cran()` — Larry only. After acceptance: tag `v0.4.1` on GitHub (the pkgdown Action also runs on `release`), bump to `0.4.1.9000`, NEWS "(development version)" heading back.
7. Remove the deprecated formals in the *following* release, not this one.

---

## 6. Where everything lives

- Task docs: `dev/tasks/cc_task_docs_refs_maxminutes_2026-10-27.md`, `dev/tasks/cc_task_prechec_tidy_2026-10-27.md`. Report (both tasks): `dev/reports/REPORT_docs_refs_maxminutes_2026-10-27.md` with `dev/reports/docs_refs_2026-10-27/` (man MD5 listings). Test: `tests/testthat/test-max-minutes-deprecation.R`.
- Prior handoff: `HANDOFF_bootstrap_efficiency_followons_2026-10-26.md` (repo root, if Larry placed it; otherwise `~/Downloads`); its §5 queue is superseded by §7 below.
- **Archive outside the repo:** `/home/larryleon/forestsearch-archive/replication-check_2026-10-27/` — 48 files, 21 MB, the uncommitted output of a 0.2.0-era worktree (`wt/replication-check`, commit `5bad3df7`, the `mr_in_replicates` RNG-neutrality check; its `REPLICATION_FINDINGS.md` §9.2 records a `ps_hat`-in-bootstrap gap that is now documented on `?forestsearch` under `ps_hat` and is one reason for today's "not tested" note on the `ps_*` arguments). The commit itself is in `master`.
- Project KB (`claude/…`): full; nothing filed today. The two arXiv PDFs are project files (names in §0).
- Machines: Pop!_OS/System76 (128 cores; documents use `n_cores` = 102), R 4.6.1, roxygen2 8.1.1, forestploter 1.2.0, speff2trial installed; GitHub Desktop; CC from the repo root.

---

## 7. Open queue (nothing blocks)

**Before 10 Nov (the vignette chat):** §3.2 items 1–5, each its own task doc, on a branch, Larry's check before merge.
**Release task:** §5.
**Later / optional:** mirror the binary-OR fast path in `run_subgroup_sims.R` (behaviour change in `R/`, needs an identity gate; ≈3 s per binary replicate); three legacy SHF documents passing binding `vi.grf.min` (leave or annotate); two data-conditional skips in `test-cross-outcome-parity.R` L262/L358 (the full suite today reported 0 skips — check whether those two ran); consistency-stage lean items B/E/G/H only if more speed is wanted; fs-glms-interpretable roadmap placement (unconfirmed since 2026-10-26; `README.md` reproduce block stale by decision).
