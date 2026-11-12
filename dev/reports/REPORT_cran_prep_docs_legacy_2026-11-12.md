# REPORT — CRAN preparation, part 3: documentation corrections, legacy GBSG wrapper, GBSG inventory, 2026-11-12

Task: `dev/tasks/TASK_cran_prep_docs_legacy_2026-11-12.md` (committed as received at the start SHA). Run unattended from `~/Documents/GitHub/forestsearch` on branch `cran-prep`. Larry's governing decision held throughout: **no function's behaviour changes**; the only code change is the Part 1 deletion.

## 1. Provenance

| Item | Value |
|---|---|
| Branch | `cran-prep` (existed; checked out from `master` at the same commit `28e71ffa`; `git diff --quiet master cran-prep -- R man tests DESCRIPTION NAMESPACE NEWS.md vignettes _pkgdown.yml` exit 0) |
| Start SHA | `d93e6621` — `dev: add CRAN preparation task, part 3` (parent `28e71ffa`) |
| `DESCRIPTION` version at start | 0.4.1.9008 (as expected) → 0.4.1.9009 |
| fs-post-selection | `main` at `5d165acce104a336e11093f1881cfdb572348462`, `git status --porcelain` empty — unchanged at the end |
| fs-glms-interpretable | `main` at `dc94a17e0e93653d5d4a9fb83bacca93bb8f87ea`, `git status --porcelain` empty — unchanged at the end |
| R / tools | R 4.6.2; devtools 2.5.2, roxygen2 8.1.1, testthat 3.3.2; survival 3.8-12 |
| Scratch | `/tmp/cran_prep_docs_legacy_2026-11-12/` (reproducers, exports, check). Nothing installed into the user or site library. |

**Commits** (all by explicit path; nothing fetched, pulled, pushed or merged):

| SHA | Commit | Part |
|---|---|---|
| `d93e6621` | `dev: add CRAN preparation task, part 3` | Step 0 (start SHA) |
| `add6cabc` | `refactor: remove four uncalled legacy GBSG functions` | 1 |
| `dee93483` | `docs: legacy GBSG wrapper marked superseded; censoring step documented` | 2 |
| `386155f3` | `docs: synthetic-data help pages match the code; tests pin the documented behaviour` | 3 |
| `37de4ab1` | `docs: DINA count scale, dmin.grf gate, plan-dependent splitting, fs_fdr_report calibration size` | 4 |
| `aa667898` | `chore: version and NEWS for part 3` | 6 step 1 |
| (this commit) | `dev: CRAN preparation report, part 3` | 6 step 7 |

**How the evidence was gathered.** The main session made the edits and ran the gates. A read-only workflow of five research agents (p1 callers, p2 censoring, p3 synthetic-data wording, p4 audit sites, p5 inventory) read the source at the start SHA via `git show`/`git grep` (never the working tree), and four of them were re-derived by independent adversarial verifiers. Their findings, with every correction folded in, are the basis of the sentence-by-sentence checks in Parts 2-5. Optional reproducers (synthetic data or `survival::gbsg`, sequential, `pkgload::load_all()` of a `git archive` of the start SHA, `Rscript --vanilla` from `/tmp/cran_prep_docs_legacy_2026-11-12/`, `timeout 120`) took about 111 s of measured R wall clock in total (main session 3.3 s; p1 7.1 s; p2 15.0 s + verifier 3.4 s; p3 11.1 s + 5.9 s; p4 40.5 s + 18.1 s; p5 3.2 s + 3.1 s), against the 15-minute cap; the longest single run was 18.05 s, against the 2-minute cap. Before the commit, the assembled report was checked against the repository by two more read-only reviewers (quoted help text against the start and head; SHAs, counts, line numbers and file lists); they found no misquote, and their five factual corrections were applied.

## 2. Status per part

| Part | Status | One line |
|---|---|---|
| Step 0 | **PASS** | Clean tree; `cran-prep` checked out; branch and `master` hold the same package; task committed as `d93e6621`; version 0.4.1.9008. |
| 1 | **PASS** | All four functions cleared the caller gate and were removed (definitions, 3 Rd pages, 2 NAMESPACE exports, 3 `_pkgdown.yml` entries); every remaining reference fixed, none BLOCKED; grep, NAMESPACE and test gates pass. |
| 2 | **PASS** | (a), (b), (c) confirmed at source; Superseded sections on `setup_gbsg_dgm()` and `calibrate_k_inter()` (replacement named: `find_k_inter_for_target_hr()`), notes on `fs_fdr_report()` and `validate_k_inter_effect()`, source comment at the censoring line; five sentences adjusted to the code (§4.3); comment-only gate passes. |
| 3 | **PASS** | Help pages of the three synthetic-data functions corrected; seven sentences adjusted to the code (§5); four tests pin the help, 0 failures; comment-only gate passes. |
| 4 | **PASS** | A-024, A-062, INV-5, A-156 documented; six departures from the task's wording, each forced by the code (§6); comment-only gate passes; no label or run-time change. |
| 5 | **PASS** (read-only, nothing changed) | 7 GBSG-specific functions inventoried with callers in package, tests, articles, quarto and both sibling repositories; leftovers and failing examples listed (§7). |
| 6 | **PASS** | 0.4.1.9009 and NEWS; code-identity gate PASS; `man/` diff lists only edited/removed pages; full suite 0 failures (652 blocks, 8,478 expectations); `check()` 0 errors / 0 warnings / 0 notes; report written and copied. |

## 3. Part 1 — four uncalled legacy GBSG functions removed

Line numbers in this section are at the start SHA `d93e6621` unless marked "head".

### 3.1 Caller gate (Step 1) — every function CLEARED

Two independent searches were made at the start SHA: the main session's `git grep` (token-aware, so the kept `.create_gbsg_dgm_` is not counted as `create_gbsg_dgm`) and a separate read-only agent (p1) that also parsed every R/ file with `codetools::findGlobals`, scanned rendered `.html`/`.tex`, text-extracted PDFs with `pdftotext`, and unzipped `.docx`/`.zip` members. They agree.

| Function | Definition (start) | Exported | Callers found | Verdict |
|---|---|---|---|---|
| `create_gbsg_dgm()` | `R/sim_aft_gbsg.R:199-231` | yes (`NAMESPACE:116`) | none in `R/`, `tests/`, `vignettes/`, `inst/`, `quarto/` (current or archive), fs-post-selection or fs-glms-interpretable | **CLEARED** |
| `simulate_from_gbsg_dgm()` | `R/sim_aft_gbsg.R:789-819` | yes (`NAMESPACE:293`) | none; one archive quarto file mentions it in two comment lines | **CLEARED** |
| `.simulate_from_gbsg_dgm_()` | `R/sim_aft_gbsg.R:823-954` | no (no Rd) | only `simulate_from_gbsg_dgm()` (`:813`), removed with it | **CLEARED** |
| `get_dgm_with_output()` | `R/sim_aft_gbsg.R:1076-1137` | no (`@keywords internal`) | none | **CLEARED** |

- Search scope and counts (p1): `tests/` 103 tracked files including the two `.rds` fixtures (read): 0 hits. `inst/` 30 files: 0. `quarto/` 1,111 in-scope source files (922 current, 189 archive by the directory rule of the earlier report): 0 callers. Rendered quarto `.html`/`.tex`/`.pdf` (1,405 files): 0 mentions of the four names. No indirect call (`do.call`, `match.fun`, `get`, `getExportedValue`, a `paste0`-built name) anywhere.
- `dev/` is outside the gate's list (build-ignored). It holds old code copies that call their own local copies, e.g. `dev/glm/scratch/R_legacy/R/sim_aft_gbsg.R:813`, `dev/code-working/legacy_code/R/sim_aft_gbsg.R:949`. Not current callers.
- Starting inventory `dev/reports/REPORT_cran_prep_seed_legacy_2026-11-11.md` Part 3 (§5.3) said the same ("No caller anywhere searched").

### 3.2 Helpers (Step 3)

No internal helper is referenced only by the removed functions, so none was removed.
- `DEFAULT_N_SUPER`, `SEED_BASE` (defaults of `.create_gbsg_dgm_`, start `:242`, `:245`), `DAYS_PER_MONTH` (`:273`), `cut_numeric`, `cut_size` (`:281-284`, `:298`) are all used by the kept builder `.create_gbsg_dgm_`.
- `.simulate_from_gbsg_dgm_` reaches no package object.
- `R/globals.R`: by free-symbol use (what `globalVariables()` suppresses) no declared name is used only by the four functions. By textual token one name, `t.sim` (`R/globals.R:208`), occurs only inside `.simulate_from_gbsg_dgm_` (`:940`, `:945`, both as `df$col` or a string). Kept: "when in doubt, keep it", and the code-identity gate then shows no `globalVariables()` difference.

### 3.3 What was removed (Step 2)

| Item | Removed |
|---|---|
| `R/sim_aft_gbsg.R` | start lines 82-232 (`create_gbsg_dgm` roxygen and definition), 741-956 ("Simulation Function" section: `simulate_from_gbsg_dgm`, `.simulate_from_gbsg_dgm_`), 1048-1139 (`get_dgm_with_output`). The file went from 1,316 to 863 lines at `add6cabc`. |
| `man/` | `create_gbsg_dgm.Rd`, `simulate_from_gbsg_dgm.Rd`, `get_dgm_with_output.Rd` (deleted by `devtools::document()`) |
| `NAMESPACE` | `export(create_gbsg_dgm)`, `export(simulate_from_gbsg_dgm)` (the other two were never exported) |
| `_pkgdown.yml` | the three reference entries `create_gbsg_dgm`, `simulate_from_gbsg_dgm`, `get_dgm_with_output` (start `:333`, `:338`, `:339`) |
| `tests/` | nothing: no `test_that()` block tests only one of the four (0 references in `tests/`) |

`create_gbsg_dgm`'s roxygen carried `@importFrom survival survreg coxph Surv` and `@importFrom stats coef quantile rexp`; every one is declared elsewhere (`R/generate_aft_dgm_main.R:189-190`), so no import line left `NAMESPACE`.

### 3.4 References fixed (Step 5)

None of the remaining references was evaluated code (the only evaluated one, `:813`, was inside a removed function), so nothing was BLOCKED.

| Site (start) | Old | New | Note |
|---|---|---|---|
| `R/setup_gbsg_dgm.R:44` | `@inheritParams create_gbsg_dgm` | the ten `@param` entries (`model` … `verbose`) copied word for word from `create_gbsg_dgm`'s roxygen | Load-bearing: without it 10 of 12 arguments would be undocumented (a check WARNING). `man/setup_gbsg_dgm.Rd`'s `\arguments` came out byte-identical. |
| `R/setup_gbsg_dgm.R:25-31` (description, details) | "This is the replacement for `create_gbsg_dgm()`: it accepts exactly the same arguments and produces the same numeric output, but returns an object of class `"aft_dgm_flex"` instead of `"gbsg_dgm"`." / "Internally the function calls `create_gbsg_dgm()` and then:" | "It builds the mechanism with the package's internal GBSG builder and returns it as an object of class `c("aft_dgm_flex", "gbsg_dgm", "list")`." / "Internally the function calls that builder and then:" | Code: it calls `.create_gbsg_dgm_` (start `:106`), takes two more arguments than the removed function (`:99-100`) and its class keeps `"gbsg_dgm"` (`:182`); the old text said otherwise. |
| `R/setup_gbsg_dgm.R:60` (`@return`) | "with all fields from `create_gbsg_dgm()` plus:" | "with the builder's fields (`df_super_rand`; the true hazard ratios `hr_causal`, `hr_H_true`, `hr_Hc_true`, `AHR`, `AHR_H_true`, `AHR_Hc_true` and the list `hazard_ratios`; `model_params`, `cens_params`, `subgroup_info`, `effect_modifiers`, `analysis_vars`, `model_type`, `n_super` and `seed`) plus:" | Field list read from the builder's `result <- list(...)` (start `R/sim_aft_gbsg.R:675-733`). |
| `R/setup_gbsg_dgm.R:75` (`@seealso`) | `\link{create_gbsg_dgm}` | `\link{generate_aft_dgm_flex}` | a list: points to the general generator, as the task says |
| `R/setup_gbsg_dgm.R:6`, `:17` (`#` header) | "calls create_gbsg_dgm()" / "wraps create_gbsg_dgm()" | "calls the builder .create_gbsg_dgm_()" / "wraps .create_gbsg_dgm_()" | comments |
| `R/sim_aft_gbsg.R:9-13` (`#` header) | key-function list naming the three removed functions | `.create_gbsg_dgm_()` (internal builder used by `setup_gbsg_dgm()`, `calibrate_k_inter()`, `validate_k_inter_effect()`), `calibrate_k_inter()`, `validate_k_inter_effect()` | comments |
| `R/sim_aft_gbsg.R:233-234` (`#`) | "used by calibrate_k_inter(), get_dgm_with_output(), scan_k_inter_effect(), and setup_gbsg_dgm(). No deprecation warning." | "used by setup_gbsg_dgm(), calibrate_k_inter() and validate_k_inter_effect()." | `scan_k_inter_effect` is defined nowhere |
| `R/sim_aft_gbsg.R:976` (`calibrate_k_inter` `@param ...`) | "passed to `\link{create_gbsg_dgm}`" | "passed to the internal GBSG builder that `\link{setup_gbsg_dgm}` also calls: `k_z3`, `z1_quantile`, `n_super`, `use_rand_params` and `seed`, with the meanings given there." | **Differs from the task's "point to generate_aft_dgm_flex()"**: the dots go to `.create_gbsg_dgm_` (start `:1024-1027`), whose arguments are `setup_gbsg_dgm`'s; `eff_at` already passes `model`, `k_treat`, `k_inter`, `cens_type`, `verbose`. |
| `R/sim_aft_gbsg.R:1216` (`validate_k_inter_effect` `@param ...`) | "passed to create_gbsg_dgm" | "passed to the internal GBSG builder that `\link{setup_gbsg_dgm}` also calls: `k_treat`, `k_z3`, `z1_quantile`, `n_super`, `cens_type`, `use_rand_params` and `seed`, with the meanings given there." | same reason (`:1254` passes `model`, `k_inter`, `verbose`) |
| `\dontrun` examples: `R/sim_aft_gbsg.R:990`, `:995` (`calibrate_k_inter`), `:1153` (`print.gbsg_dgm`), `R/oc_analyses.R:123` (`compute_dgm_cde`), `:413` (`summarize_simulation_results`) | `create_gbsg_dgm(...)` | `setup_gbsg_dgm(...)` with the same arguments | **Differs from the task's "point to generate_aft_dgm_flex()"**: these examples read `hr_H_true`/`AHR_H_true`, call `print()` or `compute_dgm_cde()` on a `gbsg_dgm` object; `generate_aft_dgm_flex()` returns `c("aft_dgm_flex", "list")` without those fields and there is no `print.aft_dgm_flex`. `setup_gbsg_dgm()` forwards the same arguments to the same builder, so it is the code-accurate substitute. Not evaluated (`\dontrun`). |
| `R/oc_analyses.R:109` (`compute_dgm_cde` `@param dgm`) | "e.g., from `create_gbsg_dgm` or `generate_aft_dgm_flex`" | "e.g., from `generate_aft_dgm_flex` or `setup_gbsg_dgm`" | |
| `R/run_simulation_analysis.R:6` (`#`), `:198` (description) | "coupled to simulate_from_gbsg_dgm()" | "coupled to the legacy GBSG simulator (since removed)" | the next sentence already points to `simulate_from_dgm` and `generate_aft_dgm_flex` |
| `README.md:131` | `` `generate_aft_dgm_flex()`, `create_gbsg_dgm()`, `run_simulation_analysis()` `` | `` `generate_aft_dgm_flex()`, `simulate_from_dgm()`, `run_simulation_analysis()` `` | text only; `README.Rmd` is not tracked, nothing re-knit |

`CLAUDE.md`, `NEWS.md` (before Part 6), `DESCRIPTION` and `tests/` had no reference.

### 3.5 Mentions listed, not edited (Step 6)

- **vignettes/** (pkgdown only): `vignettes/articles/biomarker_effects.qmd:66` (prose: "`create_gbsg_dgm()` and `generate_aft_dgm_flex()`") and its rendered copy `vignettes/articles/biomarker_effects.html:2413`.
- **quarto/**: current 0; archive 1 file, `quarto/legacy calibration files/run_simulation_analysis.R:6`, `:76` (comments naming `simulate_from_gbsg_dgm()`; an archived copy of `R/run_simulation_analysis.R`).
- **fs-post-selection** (HEAD `5d165acc`), every one prose naming `create_gbsg_dgm()`; 0 mentions of the other three names:
  - the requested supplement sentence "The simulation study's data-generating mechanisms are `create_gbsg_dgm()` and `generate_aft_dgm_flex()`": `manuscript_rev1/postselection_subgroup_inference_supplement.qmd:1184`, its `.tex:3296` and `.pdf`; the same sentence in `manuscript/…supplement.qmd:1182` (`.tex:3280`, `.pdf`), `jrssb_submission/…supplement.pdf`, `jrssb_submission/tex_zipped.zip` (member `.tex:3280`);
  - `manuscript{,_rev1}/_survey_2026-08/report_C_fb_mr_inventory.md:833` (grep output in a fenced block);
  - submission PDFs: `jasa_submission_package/ACC_form.pdf`, `jasa_submission_package/postselection_subgroup_inference_main{,_blinded}.pdf`, `jasa_v5_8/postselection_subgroup_inference_main{,_blinded}.pdf`;
  - archive/legacy (75 rows): `manuscript/_archive/jasa-2026` and `pre-trim`, `manuscript_rev1/_archive/jasa-2026` and `pre-trim`, `manuscript_legacy`, `jasa_v5_8_legacy` (`.qmd`, `.tex`, `.pdf`, `.docx` and `.zip` members).
- **fs-glms-interpretable** (HEAD `dc94a17e`), all under `dev/`, all prose naming `create_gbsg_dgm()`: `dev/reference/post_selection/manuscript_jrssb_initialsubmit/postselection_subgroup_inference_supplement.qmd:1182` (`.tex:3280`, `.pdf`), `dev/arxiv/postselection_v1/build/postselection_subgroup_inference_supplement.tex:3280`, `dev/arxiv/postselection_v1/package/postselection_subgroup_inference_supplement_body.tex:3037`, PDFs `dev/arxiv/postselection_v1/postselection_arxiv_v1.pdf`, `dev/reference/articles_books/{LeonAnderson_2026,LeonAnderson_2609.38361v1,postselection_subgroup_inference_supplement}.pdf`, the survey report `:833`; archive (26 rows) under `dev/reference/post_selection/manuscript_jrssb_initialsubmit/_archive/`.

### 3.6 Gates (Step 7) and commit (Step 8) — **PASS**

- `git grep -nE 'create_gbsg_dgm|simulate_from_gbsg_dgm|get_dgm_with_output' -- R man tests NAMESPACE _pkgdown.yml 'README*'` at `add6cabc` returns only `.create_gbsg_dgm_` lines (`R/fs_fdr_report.R:486` comment; `R/setup_gbsg_dgm.R:6`, `:17`, `:132`; `R/sim_aft_gbsg.R:10`, `:85`, `:660`, `:801`); the token-aware form `(?<![.\w])(create_gbsg_dgm|simulate_from_gbsg_dgm|get_dgm_with_output)(?!\w)|\.simulate_from_gbsg_dgm_` returns nothing.
- `git diff d93e6621 -- NAMESPACE`: exactly `-export(create_gbsg_dgm)` and `-export(simulate_from_gbsg_dgm)`.
- `test-betaHhat-contract.R`, `test-knoise-population.R`, `test-plan-restored.R` (`NOT_CRAN=true`): 44 blocks, 236 expectations, 0 failures, 0 errors, 0 skips; 23.45 s.
- `devtools::document()` wrote 7 pages and deleted the 3; `DESCRIPTION` was not touched. Commit `add6cabc` (17 paths).

## 4. Part 2 — legacy GBSG wrapper: superseded, censoring step documented

No code changed. Line numbers are at the start SHA `d93e6621` unless marked.

### 4.1 Confirmations (Step 1) — (a), (b), (c) hold

Each was confirmed by the main session's reading, by the p2 agent (three runs, 14.98 s) and by an independent verifier (one run, 3.39 s) that re-derived every line and reproduced every measured value.

**(a) The wrapper's censoring covariate term is the fitted term × (event-model scale / censoring-model scale).** `R/sim_aft_gbsg.R`:
- `:375-379`, `:383` — `fit_aft <- survival::survreg(Surv(y, event) ~ z_true, data = dfa, dist = "weibull")`; `sigma <- fit_aft$scale` (the event-model scale).
- `:626-632` — `z_cens <- as.matrix(df_big[, covs_main])`; `fit_cens <- survival::survreg(Surv(y, 1 - event) ~ z_cens, data = df_big, dist = "weibull")` (`covs_main` = `treat, z1..z5`, `:332`).
- `:634`, `:636` — `sigma_cens <- fit_cens$scale`; `gamma_cens <- stats::coef(fit_cens)[-1]`.
- `:637-638` — `b0_cens <- -gamma_cens / sigma_cens`; `b_cens <- -b0_cens * sigma`, i.e. `b_cens = (sigma / sigma_cens) * gamma_cens`.
- `:646-647` — `linC1.conf`/`linC0.conf` = `z_cens_{1,0} %*% b_cens`.
- Consumed: `R/setup_gbsg_dgm.R:157-158` renames them `lin_pred_cens_1`/`_0`; `:175-179` sets `model_params$censoring <- list(type, mu = cens_params$mu, tau = cens_params$sigma)`; `R/simulate_from_dgm.R:320-322` picks the arm's predictor and `:326-328` draws `logC_sim <- params$censoring$mu + cens_adjust + params$censoring$tau * epsilon_cens + lin_pred_cens`.
- Measured (p2 run, `setup_gbsg_dgm(model = "alt", k_inter = 2)`): sigma 0.713167, sigma_cens 0.436289, ratio 1.634621; max |`lin_pred_cens_1` − ratio × fitted term| = 4.2e-17 (against 0.048 without the ratio).

**(b) `generate_aft_dgm_flex()` takes each subject's censoring predictor from the censoring fit itself, with no such factor.** `prepare_censoring_model()` is called at `R/generate_aft_dgm_main.R:358-368` (after `calculate_hazard_ratios()`, `:345`) on `df_work`, one row per source patient with the observed treatment (`R/generate_aft_dgm_helpers.R:459`). AIC route: `fit_cens <- get_best_survreg(comparison)` (`:326`), `lin_pred_cens_0 <- predict(fit_cens, newdata, type = "linear") - mu_cens` (`:361-362`, `_1` at `:369-370`); force-fit route `:237`, `:265-273`; analytical and uniform routes set 0. The roxygen states the contract: `lin_pred_cens_*` "store the covariate contribution only — i.e. gamma_c' X" (helpers `:9-13`). Measured: force-fit Weibull, max |`lin_pred_cens_1` − X·gamma_c| = 1.9e-15. So the comparison sentence is kept.

**(c) True hazard ratios and the `k_inter` calibration do not use the censoring predictor.** In `.create_gbsg_dgm_` every true value is computed before the censoring fit, from stacked potential outcomes with `event = 1L`: `set.seed(seed)` `:478`, `epsilon` `:498`, `logT_1`/`logT_0` `:501`, `:505`, `df_po` `:509-514`, `hr_causal` `:517-520`, `hr_H_true`/`hr_Hc_true` `:527-543`, `AHR*` from `loghr_po` `:556-563`; `hazard_ratios` (`:666-673`) uses no censoring quantity. `calibrate_k_inter`'s `eff_at` returns `AHR_H_true` or `hr_H_true` (`:1028`), and `.calibrate_by_root` uses only `effect_fun(k) - target` (`R/calibrate_helpers.R:58`, `:74`). `compute_dgm_cde()` uses only `theta_0`, `theta_1` and the flag (`R/oc_analyses.R:134-165`). Measured: Weibull-minus-uniform `setup_gbsg_dgm` builds differ by exactly 0 in all six `hazard_ratios` and both CDEs; `calibrate_k_inter(1.5)` is `identical()` (1.111548258) under both censoring types.

### 4.2 Which exported functions get what (Steps 2-4)

Call graph from code lines (comments, roxygen and `formals()` reads excluded), after Part 1:

| Exported function | Chain to `.create_gbsg_dgm_()` | Returns | Documentation added |
|---|---|---|---|
| `setup_gbsg_dgm()` | `R/setup_gbsg_dgm.R:106` | the DGM (truths, super-population, scaled censoring predictor) | `@section Superseded` (Step 2) |
| `calibrate_k_inter()` | `eff_at` → `R/sim_aft_gbsg.R:1024` | a number | `@section Superseded`, own role (Step 3) |
| `fs_fdr_report()` | `setup_gbsg_dgm()` `R/fs_fdr_report.R:257`, `:263`; `calibrate_k_inter()` `:259` | `fdr` rates from censored simulated trials (`simulate_from_dgm()` `:270-271`) and `hr_structure` truths | `@note` (Step 4) |
| `validate_k_inter_effect()` | `R/sim_aft_gbsg.R:1254` | a table of true HR/AHR values only | `@note` (Step 4). It builds its mechanisms through the wrapper, but the censoring step does not reach its output, which its note says. |

No other exported function reaches the builder. `simulate_from_dgm()`, `fs_build_eval_frame()`, `compute_dgm_cde()`, `run_simulation_analysis()`, `calibrate_cens_adjust()`, `check_censoring_dgm()` and `print.gbsg_dgm()` process a wrapper DGM only if a caller passes one; they were not given a note.

**Replacement for `calibrate_k_inter()` (Step 3).** `find_k_inter_for_target_hr()` is exported (`NAMESPACE:150`). Its help: title "Find k_inter Value to Achieve Target Harm Subgroup Hazard Ratio" (`R/find_k_inter_main.R:1`), "Uses numerical root-finding to determine the interaction parameter (k_inter) that achieves a specified target hazard ratio in the harm subgroup" (`:3-5`), `@return` "`dgm`: Object of class "aft_dgm_flex"" (`:41`), `@seealso` "`generate_aft_dgm_flex` for DGM generation" (`:56`); its code builds `generate_aft_dgm_flex()` (`:106`, `:139`). It is therefore named. It targets `$hazard_ratios$harm_subgroup` (`:120`), the Cox-based HR, and has no AHR option, which the section says.

### 4.3 Help text: the task's wording and what was written

**`setup_gbsg_dgm()` — `@section Superseded`.** Sentences changed from the task's wording because the code differs (evidence from §4.1 and the p2 runs):

| Task sentence | Written | Why |
|---|---|---|
| "…`simulate_from_dgm()`, which take any survival dataset." | "…which take any survival dataset with a 0/1 treatment column." | flex reads `data[[treatment_var]]` and uses it numerically in the design (`R/generate_aft_dgm_helpers.R:459`). |
| "On the GBSG data they give the same outcome model, true hazard ratios and `k_inter` calibration once the caller has built the GBSG covariates; there is no drop-in replacement call." | adds "at `k_z3 = 1` and `use_rand_params = FALSE` (`generate_aft_dgm_flex()` has neither argument)" | `k_z3` scales the z3 coefficient (`R/sim_aft_gbsg.R:416`), `use_rand_params` refits on `meno == 0` (`:400-409`); flex's formals (`R/generate_aft_dgm_main.R:218-241`) have neither. At those defaults the p2 run measured HR/AHR/CDE differences 0, mu diff 0, max |gamma diff| 5.6e-17, calibrated k within 2.2e-7 (tol 1e-6), as the earlier report §5.4 found. |
| "The censoring step differs: here the covariate term of the censoring model is the fitted term multiplied by the ratio of the event-model scale to the censoring-model scale." | "The censoring step differs: here (with the default `cens_type = "weibull"`) a Weibull censoring model is fitted to the resampled super-population, in which treatment is the arm assigned for simulation, not the observed treatment, and its covariate term is the fitted term multiplied by the ratio…" | The factor is not the only difference. The wrapper fits on `df_big` (`n_super` resampled rows, `:481-482`, `:601`, `:628-632`), whose `treat` is the assigned arm (`:590-599`; fitted treat coefficient 0.0028, unrelated to `hormon`); flex fits on the source rows with the observed treatment (force-fit treat coefficient 0.063). Removing the factor alone would not reproduce flex. The arm is assigned by position in the resample (the first `round(n_super/2)` rows are treated, `:586-599`), so "assigned for simulation", not "at random" (verifier). Scoped to the default `cens_type = "weibull"`: under `"uniform"` no censoring model is fitted (`:625`) and `simulate_from_dgm()` cannot simulate from the result (OPEN ITEMS). |
| "True values are not affected; simulated censoring times are." | "The true hazard ratios stored in the object and the `k_inter` calibration are not affected; simulated censoring times are, and with them the observed times and event indicators and anything computed from them, such as the beta(Hhat) targets scored on a frame from `fs_build_eval_frame()`." | `fs_build_eval_frame()` simulates the super-population with censoring (`R/betaHhat_truth.R:848-851`) and the `fs_betaHhat_*` targets default to `y_sim`/`event_sim`. Measured (p2 run 3, verifier run R4): the Cox HR at the true flag on the evaluation frame is 3.2612 with the factor and 3.3043 without (H), 0.6250 and 0.6254 (Hc), while `t_true` and `treat_sim` are identical — hence "observed times and event indicators" rather than "any quantity computed from simulated data" (verifier). |
| "The "Extreme Subgroup Effects" article shows the general machinery on the GBSG data." | "The "Extreme Subgroup Effects" article on the package website shows the general machinery on the GBSG data, with a uniform treatment effect." | The article (`vignettes/articles/extreme_subgroups.qmd`, navbar text "Extreme Subgroup Effects", `_pkgdown.yml:58`; YAML title "Extreme Subgroup Effects Under a Uniform Treatment Benefit") runs `generate_aft_dgm_flex()` on `survival::gbsg` (`:532-585`) with `model = "null"` (`:543`) and `subgroup_vars = NULL` (`:544`); it embeds no harm subgroup and calibrates no `k_inter`. It is build-ignored, so website only. |

Two of the task's seven sentences are as given ("This function is kept so that published simulations can be reproduced." and "This function will be deprecated in a future release."); the other five are the rows above. Plain text, no badge, no run-time message.

**`calibrate_k_inter()` — `@section Superseded`** (new; nothing replaced):
> This function is kept so that published simulations can be reproduced; it calibrates `k_inter` for `setup_gbsg_dgm()`. New work should use `find_k_inter_for_target_hr()`, which calibrates `k_inter` for `generate_aft_dgm_flex()` on any survival dataset with a 0/1 treatment column. It targets the Cox-based harm-subgroup hazard ratio, as `use_ahr = FALSE` does here; there is no packaged equivalent of `use_ahr = TRUE`. The calibrated quantity does not depend on the censoring model, so the censoring step of the GBSG wrapper (see `setup_gbsg_dgm`) does not affect the calibrated value. This function will be deprecated in a future release.

**`fs_fdr_report()` — `@note`** (new):
> This function builds its null data through the superseded GBSG wrapper (`setup_gbsg_dgm`, with `calibrate_k_inter` for the `"harm"` null), in which the censoring model is fitted to the resampled super-population and its covariate term is multiplied by the ratio of the event-model scale to the censoring-model scale; this changes the simulated censoring times, and with them the rates in `fdr`, but not the true hazard ratios in `hr_structure`.

**`validate_k_inter_effect()` — `@note`** (new):
> This function builds its mechanisms through the superseded GBSG wrapper (see `setup_gbsg_dgm`). That wrapper's censoring step, in which the covariate term of the censoring model is multiplied by the ratio of the event-model scale to the censoring-model scale, does not affect the hazard ratios this function reports.

### 4.4 Source comment (Step 5)

Above the unchanged line `b_cens <- -b0_cens * sigma` (head `R/sim_aft_gbsg.R`, in `.create_gbsg_dgm_`):
```r
    # b_cens is the fitted censoring term times sigma / sigma_cens. Kept as is so
    # that published simulations reproduce; true values (hazard ratios, k_inter
    # calibration) do not depend on it. See ?setup_gbsg_dgm, section Superseded.
```

### 4.5 Gate (Step 7) and commit (Step 8) — **PASS**

`git diff -U0 add6cabc -- <file> | grep -E '^[+-]' | grep -vE '^(\+\+\+|---)' | grep -vE "^[+-][[:space:]]*#"` is empty for `R/fs_fdr_report.R`, `R/setup_gbsg_dgm.R`, `R/sim_aft_gbsg.R`. `devtools::document()` rewrote `setup_gbsg_dgm.Rd`, `calibrate_k_inter.Rd`, `fs_fdr_report.Rd`, `validate_k_inter_effect.Rd`; each passes `tools::checkRd()`. Commit `dee93483`.

## 5. Part 3 — synthetic-data help pages match the code; tests pin them

No code changed (`R/synthetic_data_perturbation.R` and two example lines elsewhere). Every proposed sentence was checked against the code by a read-only agent (p3, three runs, 11.05 s) and re-derived by an independent verifier (two runs, 5.92 s). Line numbers are at the start SHA.

### 5.1 `generate_bootstrap_with_noise()`

Code facts the wording rests on: `set.seed(seed)` `:506`; continuous auto-detection `names(data)[vapply(data, is.numeric, logical(1))]` minus `id_var` (`:515-517`, never reads `cat_vars`); categorical auto-detection `is.factor(x) || is.logical(x) || (is.numeric(x) && length(unique(x)) <= 10L)` minus `c(continuous_vars, id_var)` (`:522-528`; character columns never selected); first draw `sample(nrow(data), n, replace = TRUE)` `:542`; noise sd `sd(x) * noise_level` `:549`, clamp `:554-556`, round if every non-missing value is whole `:559-560`; change probability `noise_level / 2` `:578`; `is.ordered(...) || var == "grade"` weighting only when `n_levels > 2` (`:589`, `:597`, `:604`).

| Piece | Task wording | Written | Changed because |
|---|---|---|---|
| Description | "Creates a bootstrap sample from a dataset and adds controlled noise. Numeric columns receive Gaussian noise. Factor and logical columns, and any column named in `cat_vars`, are changed with a small probability. Numeric-coded categories (for example 0/1 indicators) are treated as continuous unless they are named in `cat_vars`. This function is useful…" | "Creates a bootstrap sample from a dataset and adds controlled noise. Continuous columns receive Gaussian noise; categorical columns are changed with a small probability. By default every numeric column except `id_var` is continuous and every factor or logical column is categorical. Numeric-coded categories (for example 0/1 indicators) are therefore treated as continuous unless `continuous_vars` is given without them; naming them in `cat_vars` alone adds the categorical change without removing the noise. This function is useful…" (last sentence as given) | "Numeric columns receive noise" and "factor and logical columns are changed" hold only under auto-detection (noise goes to `continuous_vars` only, `:546`; factor/logical only when `cat_vars` is NULL, `:520-529`). "Treated as continuous unless they are named in `cat_vars`" is **false**: with `continuous_vars` NULL a column named in `cat_vars` stays continuous too (`:514-518`; measured: in both sets), which the task's own `cat_vars` text says. |
| `@param continuous_vars` | as given | as given | TRUE (`:515-517`) |
| `@param cat_vars` | "…If NULL (default), factor and logical columns are used, together with numeric columns that have 10 or fewer unique values and are not in `continuous_vars`. When `continuous_vars` is also NULL, every numeric column is already in the continuous set, so no numeric column is treated as categorical. When `cat_vars` is given…" | "…If NULL (default), factor and logical columns are used, together with numeric columns that have 10 or fewer unique values; columns in `continuous_vars` and `id_var` are left out, and character columns are not used. When `continuous_vars` is also NULL, every numeric column other than `id_var` is already in the continuous set, so no numeric column is treated as categorical. When `cat_vars` is given…" (rest as given) | `:528` drops `id_var` and `continuous_vars` from every type, not only numeric; character columns are never selected (`:523-524`; measured). |
| Details, continuous | add "Integer-valued columns are rounded after the noise is added, so a 0/1 column treated as continuous is usually unchanged at a small `noise_level`." | "Integer-valued columns (every non-missing value a whole number) are rounded after the noise is added, so a 0/1 column treated as continuous is usually unchanged at a small `noise_level`." | The test is `all(x == round(x), na.rm = TRUE)` (`:559`), not `is.integer`; the parenthesis says so. Measured: 0 of 2,000 rows changed at 0.1. |
| Details, categorical, first sentence | "Applies to the columns in `cat_vars`, given or auto-detected as described for that argument." | as given | TRUE (`:566`) |
| Details, ordered item | "Ordered factors, or a column named `grade`, with more than two levels: selection is weighted toward adjacent levels." (old: "Ordered factors: weights selection toward adjacent levels") | as given | TRUE (`:589`, `:597`, `:604`; weights `1/(d + 0.5)`) |
| `@note`, third item | "…A numeric `grade` column is processed as categorical only when it is named in `cat_vars`, or when `continuous_vars` is given without it." | "…A numeric `grade` column is processed as categorical only when it is named in `cat_vars`, or when `cat_vars` is NULL, `continuous_vars` is given without it and it has 10 or fewer unique values." (first sentence as given) | With `cat_vars` given without `grade`, `grade` is in neither set (measured: identical to the resample); auto-detection needs `<= 10L` (`:524`). |
| Examples | `data(gbsg, package = "survival")` → `gbsg <- survival::gbsg`; comment "all gbsg columns are numeric, so all are treated as continuous" | `gbsg <- survival::gbsg`; comment "(all gbsg columns are numeric, so all except pid, the id_var, are treated as continuous)" | `pid` is the default `id_var` (`:485`), excluded (`:517`) and reset to `1:n` (`:631-633`). |
| Source comment above auto-detection | "By design (L. Leon, 2026-11-12): under auto-detection every numeric column is continuous. Not a defect; see the help page." | "…every numeric column except id_var is continuous…" | `:517` |

Left unchanged on the page as the task directs ("Everything else on the page stays"), but inaccurate — see OPEN ITEMS: Details item "Preserves integer type for variables that appear to be integers" (an integer column comes back double, measured); `@note` item "Missing values (NA) … are not imputed" (a selected NA in a categorical column is replaced by a level; an ordered factor or `grade` with > 2 levels and an NA stops with "replacement has length zero", measured); `@return` "same structure as the input data" (integer → double).

### 5.2 `generate_bootstrap_synthetic()`

| Piece | Task wording | Written | Changed because |
|---|---|---|---|
| `@param ordinal_vars` | "Optional character vector naming which of the `cat_vars` are ordinal. A variable listed in both `cat_vars` and `ordinal_vars` is changed to a neighbouring level instead of a random other level. A variable listed only in `ordinal_vars` is not perturbed: it is kept as resampled and is named in the 'kept as-is' message. Known limitation: for a numeric-coded variable at its lowest or highest level, the new value is not always the adjacent level." | "Optional character vector naming which of the `cat_vars` are ordinal. A numeric or character variable listed in both `cat_vars` and `ordinal_vars` is changed to a neighbouring level (in sorted order) instead of a random other level. A factor listed in both is set to NA, with a warning, at each row selected for change, unless its levels are "1", "2", ... in that order. A variable listed only in `ordinal_vars` is not perturbed: it is kept as resampled and is named in the 'kept as-is' message. Known limitation: for a numeric-coded variable at its lowest or highest level, the new value is not always the adjacent level (with codes 1, ..., k, a row at 1 that is selected for change keeps its value half the time, and a row at k moves to any of 1, ..., k - 1 with equal probability)." | **The second task sentence is false for a factor.** `c()` at `:163-165` drops the factor class, so the integer code drawn at `:168` is assigned to the factor and gives NA with "invalid factor level, NA generated" (measured 291/291 and, by the verifier, 313/313 changed rows), unless the level labels are the integers in order. Neighbours come from `sort(unique())` (`:154`). The limitation's mechanism is `sample()`'s length-1 rule: one number `v >= 1` draws from `1:v` (measured on codes 1:3, 1:5; codes 0:2 always move to 1). |
| Both examples | add at the `ordinal_vars` line: "listed only in ordinal_vars, so it is resampled but not perturbed; list it in cat_vars as well to perturb it" | as given, as two comment lines above `ordinal_vars = c("grade"),` and `ordinal_vars = c("education"),` | TRUE (Ex1 `:23-24`, Ex2 `:45-46`) |
| Source comment at the categorical loop | "Behaviour kept as is by decision (2026-11-12): the loop runs over cat_vars only; see the help page for ordinal_vars." | as given | `:145` |
| Source comment at `sample(possible_new, 1)` | "Known limitation, documented in the help page; behaviour kept as is by decision (2026-11-12)." | as given, plus "(At an end level possible_new holds one number v, and sample(v, 1) draws from 1:v when v >= 1.)" | states the mechanism for the next reader |

### 5.3 `generate_gbsg_bootstrap_general()`

| Piece | Task wording | Written | Note |
|---|---|---|---|
| Description (the page had only a title, so the Rd `\description` repeated it) | add "`grade` is passed in `ordinal_vars` only, so it is resampled but not perturbed." | "Calls `generate_bootstrap_synthetic()` on `survival::gbsg` with fixed variable lists. `grade` is passed in `ordinal_vars` only, so it is resampled but not perturbed." | The first sentence is added so the new description paragraph says what the function does (`:366-378`). TRUE (`:372`; T1). |
| `@section Superseded` | as given | as given | TRUE: the wrapper's body is the single call `:368-378`; the first example `:20-29` has the same lists and `preserve_bounds = TRUE` is the default (`:61`); `identical()` TRUE (T2). |

### 5.4 The same example pattern elsewhere

- `R/pareto_frontier_table.R:124` `data(gbsg, package = "survival")` → `gbsg <- survival::gbsg` (fails as written in a fresh session: "data set 'gbsg' not found").
- `R/calibrate_k_treat.R:55` — the text is `data(gbsg)` after `library(survival)` (`:54`), not the package form; it works with the same warning. Replaced by `gbsg <- survival::gbsg` as directed.
- Listed, not edited: the task's pointer `vignettes/articles/extreme_subgroups.qmd:28` is a `stop()` version check; the `data(gbsg, package = "survival")` line is at `:282`, after `library(survival)` at `:272`, so it works with the warning.

### 5.5 Gates and commit (3.6) — **PASS**

Comment-only gate against `dee93483` empty for `R/synthetic_data_perturbation.R`, `R/pareto_frontier_table.R`, `R/calibrate_k_treat.R`; `devtools::document()` rewrote the five pages, each passing `tools::checkRd()`; the two test files: 0 failures (below). Commit `386155f3`.

### 5.6 Tests that pin the help (new file `tests/testthat/test-synthetic-data-documented-behaviour.R`)

Expected values come from the resample, not from the functions: after `set.seed(seed)` the first random call in `generate_bootstrap_synthetic()` (start `:98`, `:101`) and `generate_bootstrap_with_noise()` (`:506`, `:542`) is `sample(nrow(data), n, replace = TRUE)`, and nothing between draws, so each test rebuilds the row indices itself. Seeds and bounds were set before the first run and not changed; no help wording needed correcting after a test run.

| Test | Setting | Expectation | Measured | Result |
|---|---|---|---|---|
| T1 | `generate_gbsg_bootstrap_general(n = 200, seed = 11)` | `grade` `identical()` to `survival::gbsg$grade[idx]` | identical | PASS (1 expectation) |
| T2 | wrapper vs `generate_bootstrap_synthetic()` with the first example's lists, `n = 200, seed = 11, noise_level = 0.15` | `identical()` | identical | PASS (1) |
| T3 | `generate_bootstrap_with_noise()` under auto-detection; made-up 2,000 rows: `pid` (the default `id_var`), continuous `x`, numeric 0/1 `b` (`rbinom`, integer), 3-level factor `f`; `seed = 7`, `noise_level = 0.1` | `b` equals its resample (`expect_equal`: the continuous step returns the integer column as double, so `identical()` would fail on type); share of `f` rows changed in (0, 0.15) | `b` changed in 0 of 2,000 rows; `f` changed in 4.9% | PASS (3) |
| T4 | the same data, `continuous_vars = "x"`, `cat_vars = "b"` | share of `b` rows changed in (0, 0.15) | 4.65%; `f` identical to its resample | PASS (2) |

Run with `test-gbsg-bootstrap-general.R` (`NOT_CRAN=true`, `devtools::test(filter = …)`): 7 blocks, 29 expectations (new file 4 / 7; existing 3 / 22), 0 failures, 0 errors, 0 skips; 2.99 s. The p3 agent's draft of the same tests (seed 11) passed 7 of 7, and its verifier re-ran it.

## 6. Part 4 — four audit findings handled by documentation

No code changed, no printed label changed, no run-time warning added. Sites were located by content (the task's line numbers were taken at `e92d3c96`/`28e71ffa`; at the start SHA they had moved). Each row's documentation-only option comes from `dev/reports/REPORT_cran_prep_2026-11-11_part3_detail.md`. Every sentence was checked against the code by a read-only agent (p4, four runs, 40.5 s) and an independent verifier (one run, 18.05 s), whose corrections were applied before the commit. Line numbers are at the start SHA.

### 6.1 A-024 — DINA's Poisson quantity

Code: `dina_fit()` and `dina()` have no offset or exposure formal (`R/dina.R:157-167`, `R/dina_df.R:84-90`; `dina(..., offset = "fu")` fails with "unused argument"); the within-arm Poisson fit is of the raw count (`R/dina.R:541`), and the step-2 `glm.fit` carries only the nuisance offset (`:843-847`). So tau(x) = log E[Y|x,W=1] − log E[Y|x,W=0] = log IRR(x) + log(E[fu|x,W=1]/E[fu|x,W=0]). `forestsearch()` never passes `offset.name` to DINA (`R/forestsearch_helpers.R:1534-1537`; use_dina path `R/forestsearch_main.R:3216`). Measured (p4 run 1, true IRR 2): treated follow-up halved — DINA mean tau-hat 0.006, log mean-count ratio −0.006, offset log IRR 0.688; equal follow-up — 0.733, 0.739, 0.739.

| Site | Old | New |
|---|---|---|
| `R/dina_subgroup.R:138-142` (`dina_subgroup()` `@param m_diff`, "Scale") | "…`log(IRR)` for `poisson`, and the raw mean difference…" | "…the log mean-count ratio for `poisson` (DINA takes no exposure offset, so this equals `log(IRR)` only when mean follow-up does not depend on arm given the covariates), and the raw mean difference…" |
| `R/forestsearch_main.R:444-446` (`forestsearch()` `@param dina_args`) | "…`log(IRR)` for `poisson`, and the raw mean difference…" | "…the log mean-count ratio for `poisson` (DINA takes no exposure offset -- `offset.name` is not passed to its fit -- so this equals `log(IRR)` only when mean follow-up does not depend on arm given the covariates), and the raw mean difference…" |
| `R/dina_subgroup_refit.R:41-46` (description) | "…(log-HR for Cox, log-OR for binomial, log-IRR for Poisson, mean difference for Gaussian), directly comparable…" | "…(log-HR for Cox, log-OR for binomial, the log mean-count ratio for Poisson, mean difference for Gaussian), directly comparable … The Poisson fit takes no exposure offset, so its effect equals log(IRR) only when mean follow-up does not depend on arm within the subgroup (given any adjustment covariates); its `effect_scale` label remains `"log-IRR"`." |

Changed from the task's wording:
- "only when follow-up does not depend on arm" → "only when **mean** follow-up does not depend on arm": by the identity above the exact condition is equal mean follow-up; arm-independent follow-up is sufficient, not necessary.
- For the refit, "given the covariates" → "within the subgroup (given any adjustment covariates)": the refit is unadjusted by default (`confounders = "none"`, `R/dina_subgroup_refit.R:155`; `rhs_terms <- c(treatment, confounders)` `:334`; `glm` with no offset `:366-369`).

Not changed: the printed label `poisson = list(label = "log-IRR", ratio = TRUE)` (`R/dina_subgroup_refit.R:447`) and the `print()` that derives "IRR" from it (`:495-497`), as the task directs; the `@return` lines `:103`, `:107` that describe that literal label; and the `#` comment `R/dina_subgroup.R:402` ("log-IRR for poisson"), which is not a help page (listed in OPEN ITEMS). `R/dina.R:32` and `:953` already say "log mean ratio".

### 6.2 A-062 — `dmin.grf`

Code: `dmin.grf` thresholds a DR harm effect, `mean(ctrl[S]) - mean(trt[S])` of the doubly robust scores (`R/grf_subgroup_labels.R:274`, `:307`; gate `:361-362`); for survival the forest is `causal_survival_forest(..., horizon = tau.rmst)` with grf's default RMST target (`R/grf_helpers.R:74`, `:84`; grf 2.6.1), so it is an RMST difference. Under `grf_selection = "frontier"` with `grf_select_statistic = "effect"` (both defaults under `subgroup_method = "grf"`): if no DR candidate reaches it, `.grf_frontier_select()` returns NULL and the no-selection result has no `$candidates` (`R/grf_main.R:294-305`), so the effect re-selection is skipped (`R/forestsearch_helpers.R:1916`) and no subgroup is returned; otherwise the family attached at `R/grf_main.R:334` is filtered by `n.min` only (`:272-279`) and re-ranked against the admission effect floor (`R/forestsearch_helpers.R:1773`). Fallbacks: nothing scorable (`:1762`), or an error outside the ranking step (outer `tryCatch` `:1917-1927`), keeps the native DR winner, which did pass `dmin.grf`; an empty admitted set, including an error inside the ranking step (`:1777-1781`), returns no subgroup (`:1782-1792`).

| Site | Old | New |
|---|---|---|
| `R/forestsearch_main.R:555` `@param dmin.grf` | "Numeric. Minimum events for GRF. Default 0.0." | "Numeric. Harm threshold on GRF's doubly robust (DR) scale, not a number of events. It is compared with a subgroup's DR harm effect, the subgroup mean of the DR score for control minus the subgroup mean of the DR score for treatment (the two columns of `policytree::double_robust_scores()`, computed for every patient): for survival a difference in restricted mean survival time up to the GRF horizon (see `frac.tau`); for other outcomes a difference on the scale of the outcome given to GRF (see `grf_count_transform`). Under `grf_selection = "tree"` (and in `use_grf` screening) the policy-tree subgroup with the largest DR harm effect must reach it; under `grf_selection = "frontier"` with `grf_select_statistic = "dr"` it is applied to every frontier candidate. Under `grf_select_statistic = "effect"` it is applied first as an existence gate on the DR scale and is then not applied per candidate; see `grf_select_statistic`. Default 0.0." |
| `R/forestsearch_main.R:542-554` `@param grf_select_statistic` | (no statement about `dmin.grf`) | inserted before "Ignored for `grf_selection = "tree"`…": "Under `"effect"`, `dmin.grf` is applied first, on the DR scale, as an existence gate: if no DR candidate reaches it, no subgroup is returned. Otherwise it is not applied per candidate: the whole DR-candidate family (filtered by `n.min`, not by `dmin.grf`) is re-ranked, and admission uses the effect floor from `hr.threshold` (the resolved `effect.threshold` for GLM outcomes). If no candidate can be scored on the effect, the DR winner, which did reach `dmin.grf`, is kept; if candidates are scored but none clears the effect floor, no subgroup is returned." |

Changed from the task's wording: "is then not applied per candidate" carries the fallback qualifier in `grf_select_statistic` (`R/forestsearch_helpers.R:1762`, `:1782-1792`). The verifier corrected two of the first draft's sentences before the commit: the DR harm effect is the subgroup mean of the two per-action DR-score columns (computed for every patient), not a mean over control-arm and treated-arm patients (`R/grf_subgroup_labels.R:251-260`); and a failed re-ranking returns no subgroup rather than keeping the DR winner. The existence-gate wording is exact for `dmin.grf >= 0` (the default 0.0); with `dmin.grf < 0` and `sg_focus = "hrMaxSG"`/`"hrMinSG"` an empty band can make `.grf_sg_def_from_candidate()` fail (code reading, not run; OPEN ITEMS). Not changed: the `#` comment `R/forestsearch_main.R:2457-2458` ("RMST scale (event count threshold over a window)"), which is not a help page (OPEN ITEMS).

### 6.3 INV-5 — plan-dependent consistency under literal splitting

Code: the plain loop (sequential plan, a direct `subgroup.consistency()` call with `parallel_args = list()`, or missing `workers`) seeds the caller's RNG once (`R/subgroup_consistency_main.R:557-559`) and draws every candidate's splits from that stream (`:848-853`); the batched path gives candidate m its own L'Ecuyer-CMRG stream (`:1030-1034`, `future.seed = candidate_seeds[batch_indices]` `:1051-1056`). Literal splitting runs under `"split"`, under `"resample"` for IRD (`.glm_resample_supported()` admits `c("OR","RR","RD","MD","IRR")`, `R/subgroup_consistency_helpers.R:343-345`; `.consistency_glm_pieces` returns NULL for IRD, `R/consistency_resample.R:249`), and for a candidate whose closed form is NA (`helpers:1492-1496`, `:1520-1523`, two-stage `:1890`). A candidate with a closed form draws no RNG (`R/consistency_resample.R:456`). The `sg_focus = "maxeff"` fast path (`R/subgroup_consistency_main.R:818-846`) runs before the branch, so every plan runs the same loop. The detail report's reproducer (10.5 s) measured Pcons different on 9 of 13 candidates under `"split"` and identical under `"resample"`.

| Site | Old | New |
|---|---|---|
| `R/subgroup_consistency_main.R:322-329` (`subgroup.consistency()` `@param parallel_args`, `plan`) | "…Results under `consistency_method = "resample"` are identical on every path (the closed form consumes no RNG). Under `"split"` the plain loop draws from a single RNG stream seeded once, whereas the batched backend gives each candidate its own L'Ecuyer-CMRG stream indexed by global position, so a sequential plan's `p.consistency` values differ from a parallel plan's within Monte-Carlo noise at the stated `fs.splits`; invariance to worker count holds within any plan." | "…Under `consistency_method = "resample"` a candidate whose closed form is available consumes no RNG, so its result is identical on every path. Literal splitting -- `"split"`, an IRD effect under `"resample"`, or a candidate whose closed form is not available and falls back to splitting -- draws, in the plain loop, from a single RNG stream seeded once, whereas the batched backend gives each candidate its own L'Ecuyer-CMRG stream indexed by global position, so a sequential plan's `p.consistency` values can differ from a parallel plan's within Monte-Carlo noise at the stated `n.splits`, and admission and selection can differ too; invariance to worker count holds within any plan (with a non-NULL `seed`). Under `sg_focus = "maxeff"` every plan runs the same plain loop." |
| `R/subgroup_consistency_main.R:774-776` (`#` comment) | "Under consistency_method = "resample" no RNG is consumed and every path gives identical results." | "Under consistency_method = "resample" a candidate with an available closed form consumes no RNG and gives the same result on every path. Literal splitting ("split", IRD under "resample", or a candidate whose closed form is NA and falls back to splitting) draws from that single stream here but from per-candidate streams on the batched path, so Pcons -- and possibly admission and selection -- can differ between a sequential and a parallel plan. The maxeff fast path runs the same loop on every plan." |
| `forestsearch()` `@param parallel_args` (after `R/forestsearch_main.R:377`) and `@param consistency_method` (after `:1048`) | — | the same added sentence: "A sequential and a parallel plan can give different consistency proportions, and possibly a different selection, whenever literal splitting runs (`consistency_method = "split"`, an incidence-rate difference (`effect_measure = "IRD"`) under `"resample"`, or a candidate whose closed form is not available), except under `sg_focus = "maxeff"`, which runs the same loop on every plan." |

Changed from the task's wording:
- "a risk or rate difference under `"resample"`" → "an incidence-rate difference (`effect_measure = "IRD"`) under `"resample"`": RD has a closed form (identity-link binomial, `R/consistency_resample.R:285-291`) and splits only when that fit fails, which "a candidate whose closed form is not available" already covers.
- Added "except under `sg_focus = "maxeff"`": without it "whenever literal splitting runs" is false for maxeff (`:818-846`).
- `fs.splits` → `n.splits` in the rewritten `subgroup.consistency()` sentence: that function's argument is `n.splits`. "differ" → "can differ" (a candidate whose splits all agree is identical across plans; 4 of 13 in the reproducer), and "(with a non-NULL `seed`)": with `seed = NULL` the batched path uses `future.seed = TRUE` (`:1030-1034`, `:1054-1058`). Both from the verifier.
- The new sentence does not repeat the adjacent claim (`R/forestsearch_main.R:1045`) that "propensity-adjusted effects" fall back to split; the code does not do that (OPEN ITEMS).

### 6.4 A-156 — `fs_fdr_report()` calibration size

Code: `k_harm <- calibrate_k_inter(target_hr_harm = 1.0, model = "alt", use_ahr = FALSE)` (`R/fs_fdr_report.R:259-260`) passes no `n_super`, so `eff_at` builds on `DEFAULT_N_SUPER` = 5000 (`R/sim_aft_gbsg.R:35`, `:242`, `:1024-1027`); the harm DGM is then built on `n_super` (default 25000, `:188`, `:263`). Inside H both arms share one epsilon (`:498-506`) and T_1/T_0 = exp(gamma_treat + k·gamma_zh) for every patient, with coefficients fitted to `survival::gbsg` (`:270`, `:375-379`), so the root k* = −gamma_treat/gamma_zh = 0.5679168436 is structural. Measured (p4 run 2): k at 5000 = 0.5679169493, at 2000 = 0.5679168477, at 25000 = 0.5679164737; k(25000) − k(5000) = −4.76e-7, inside `uniroot`'s tol 1e-6. (The Cox `hr_H_true` at a fixed k does depend on size: 1.006568 at 5000, 1.001587 at 25000.)

| Site | Old | New |
|---|---|---|
| `R/fs_fdr_report.R:82` `@param n_super` | "Integer; super-population size for the DGM. Default 25000." | same, plus "Under `null = "harm"`, `k_inter` is calibrated by `calibrate_k_inter` on the GBSG builder's default super-population size (5000), not on `n_super`; the calibrated value does not depend on that size beyond the root-finding tolerance (1e-6), because the root is fixed by the AFT coefficients fitted to `survival::gbsg`." |

Changed from the task's wording: "does not depend on that size" → "does not depend on that size beyond the root-finding tolerance (1e-6)": the returned k moves by 4.8e-7 between 5000 and 25000.


### 6.5 Gate and commit — **PASS**

Comment-only gate against `386155f3` empty for `R/dina_subgroup.R`, `R/dina_subgroup_refit.R`, `R/forestsearch_main.R`, `R/fs_fdr_report.R`, `R/subgroup_consistency_main.R`; `devtools::document()` rewrote `dina_subgroup.Rd`, `dina_subgroup_refit.Rd`, `forestsearch.Rd`, `fs_fdr_report.Rd`, `subgroup.consistency.Rd`, each passing `tools::checkRd()`. Commit `37de4ab1`.

## 7. Part 5 — inventory of GBSG-specific functions (read-only)

Nothing was changed for this part. Rule: a function is GBSG-specific if it reads `survival::gbsg` (or builds or requires GBSG columns) and has no argument through which another dataset can be given. Read by a read-only agent (p5, one R run, 3.16 s) at the start SHA and re-checked by an independent verifier (one R run, 3.07 s, plus its own call scanner), which confirmed every row and count and added the items folded in below; Parts 1-4 change no function body, so the inventory holds at the branch head, and the **File:line column is at the branch head** (`aa667898`) — the definition lines were re-located there by parsing `R/`.

**Counting rule (quarto/ and sibling repositories).** A call is the name followed by `(` in code (`.R`/`.r` files; `{r}` chunks of `.qmd`/`.Rmd`/`.rmarkdown`), with comments and strings blanked, excluding chunks with `eval = FALSE` / `#| eval: false`; R-level `if` guards are counted. Archive = any directory component matching `(^|[^a-z])(archive[sd]?|legacy|broken|old|superseded|deprecated|obsolete)([^a-z]|$)`, ignoring case; current counts are given with the archive count in parentheses. The scanner reproduces the earlier verified counts exactly (setup 281/297, calibrate 277/282). Sibling repositories at HEAD (fs-post-selection `5d165acc`, fs-glms-interpretable `dc94a17e`).

| Function | File:line (head) | Exported | In CRAN 0.4.0 | What is hard-wired | Callers: package code | Callers: tests | Callers: vignettes and articles | Callers: current quarto (files / call lines; archive in parentheses) | Callers: fs-post-selection; fs-glms-interpretable (files / call lines) | General equivalent, and what a caller must add to match it |
|---|---|---|---|---|---|---|---|---|---|---|
| `.create_gbsg_dgm_` | `R/sim_aft_gbsg.R:85` | no (untitled `@keywords internal`, no Rd) | yes, internal (file blob identical at `v0.4.0`) | No data argument (formals `model, k_treat, k_inter, k_z3, z1_quantile, n_super, cens_type, use_rand_params, seed, verbose`). Reads `survival::gbsg`; `y = rfstime / 30.4375`, `event = status == 1`, `treat = hormon`; `z1 = er <= quantile(er, z1_quantile)`, `z2/z4/z5 = cut_numeric(age/pgr/nodes)` (quartile codes 1-4), `z3 = meno == 0`, `zh = treat·z1·z3`, harm `z1 & z3`; `v1..v7` (incl. `cut_size(size)`, `grade3`); outcome covariates `treat, z1..z5 (+ zh)`; `use_rand_params` refits on `meno == 0`; Weibull censoring fit with the `sigma / sigma_cens` factor; class `gbsg_dgm` | `setup_gbsg_dgm` (`R/setup_gbsg_dgm.R`), `eff_at` in `calibrate_k_inter`, `validate_k_inter_effect` | none directly; indirectly via `test-knoise-population.R`, `test-betaHhat-contract.R`, `test-plan-restored.R` | none directly; under the evaluated calls of `paper_simulations.qmd` | 0 / 0 (0 / 0); mentioned on 19 lines in 9 files | 0; 0 | `generate_aft_dgm_flex()`: the caller builds `y`, `event`, `treat`, `z1..z5` from `survival::gbsg` as the builder does (months scale for `mu` to match), passes `z2, z4, z5` as `continuous_vars` (binary `z1`, `z3` work either way), and gives the subgroup as `subgroup_vars = c("er","meno")`, `subgroup_cuts = list(er = list(type = "quantile", value = 0.25), meno = 0)`. No flex argument matches `k_z3 != 1` or `use_rand_params = TRUE`; the censoring step differs (Part 2). |
| `setup_gbsg_dgm` | `R/setup_gbsg_dgm.R:136` | yes (`NAMESPACE:285`) | yes (`v0.4.0:NAMESPACE:259`) | Calls `.create_gbsg_dgm_` and has no data argument; literal defaults `n_super = 5000L`, `seed = 8316951L`; renames `df_super_rand` columns into `df_super`; class `c("aft_dgm_flex", "gbsg_dgm", "list")` | `fs_fdr_report` (two calls); named in a `stop()` string in `R/fs_dgm_feasibility.R` (not a call) | `test-knoise-population.R` (helper, 5 blocks), `test-betaHhat-contract.R` (helper `.small_surv_dgm`), `test-plan-restored.R` via `fs_fdr_report` (`skip_on_cran`) | `paper_simulations.qmd` evaluated at `:507`, `:614`, `:645`, `:699` (not evaluated: `:375`, `:1756`, `:274`); `treatment_effect_definitions.qmd:1034` (`eval=FALSE`); pkgdown only | **281 / 297** (13 / 17) | fs-post-selection **30 / 42**: current 27 / 36 (`analysis/simulations/gbsg/resampling_comparisons` 18 / 18; `dgm_truth_targets{,_h10,_h15}.qmd` in `manuscript_rev1`, `manuscript`, `manuscript_v40`, 3 / 6 each), archive (`manuscript_legacy`) 3 / 6; fs-glms-interpretable **3 / 6** (`dev/reference/post_selection/manuscript_jrssb_initialsubmit/dgm_truth_targets*.qmd`) | `generate_aft_dgm_flex()` + `simulate_from_dgm()`: everything in the row above, plus the `k_random_noise` columns drawn by the caller onto `dgm$df_super`; read `dgm$hazard_ratios` instead of `hr_H_true` / `AHR_H_true`. |
| `calibrate_k_inter` | `R/sim_aft_gbsg.R:654` (nested `eff_at` `:673`) | yes (`NAMESPACE:92`) | yes (`v0.4.0:NAMESPACE:78`) | Every candidate DGM is `.create_gbsg_dgm_`; `...` goes there (no data argument); `model == "alt"` enforced | `fs_fdr_report` | `test-betaHhat-contract.R` (`.small_surv_dgm`) | `paper_simulations.qmd` evaluated at `:604`, `:635` (not evaluated `:1749`, `:273`, `:1753`); prose in `extreme_subgroups.qmd:202`, `:485` | **277 / 282** (10 / 10) | fs-post-selection **30 / 42** (same files as above); fs-glms-interpretable **3 / 6** | `find_k_inter_for_target_hr()` (exported): the caller adds the data and covariate build. It targets the Cox HR only (no `use_ahr`), returns a list, has no `...` (so a non-default seed or other flex argument cannot be forwarded; the default seeds coincide, 8316951), and has different defaults (range `c(-10, 10)`, `tol = 0.001`, `verbose = TRUE`). |
| `validate_k_inter_effect` | `R/sim_aft_gbsg.R:796` | yes (`NAMESPACE:310`) | yes (`v0.4.0:NAMESPACE:285`) | `.create_gbsg_dgm_(model = "alt", k_inter = k, ...)` per value; no data argument | none | none | `paper_simulations.qmd:683`, evaluated | 0 / 0 (0 / 0) | 0; 0 | `sensitivity_analysis_k_inter()` (exported): data and DGM arguments, `plot = FALSE`, an evenly spaced grid; for AHR columns build each DGM and read `hazard_ratios$AHR_harm` / `AHR_no_harm`. No PASS/CHECK print-out. |
| `print.gbsg_dgm` | `R/sim_aft_gbsg.R:715` | S3 method (`NAMESPACE:48`) | yes (`v0.4.0:NAMESPACE:43`) | Header "GBSG-Based AFT Data Generating Mechanism"; prints `effect_modifiers$k_z3`, `subgroup_info$definition`, the ER threshold; class `gbsg_dgm` is set only by the GBSG builder | dispatch only (no `print.aft_dgm_flex` exists) | none (a `gbsg_dgm` stub in `test-fs-dgm-feasibility.R` is not printed) | implicit `print(dgm_alt)` at `paper_simulations.qmd:523` | 0 explicit; 0 implicit found (0) | 0; 0 | none; print `dgm$hazard_ratios`. |
| `fs_fdr_report` | `R/fs_fdr_report.R:195` | yes (`NAMESPACE:184`) | yes (`v0.4.0:NAMESPACE:170`; the GBSG wiring is the same) | Both nulls fixed to the GBSG wrapper (`setup_gbsg_dgm` ×2, `calibrate_k_inter(target_hr_harm = 1.0)`); no DGM argument (`df_analysis` sets only `n`); fallback confounders `.FS_GBSG_CONFOUNDERS` (`R/fs_fdr_report.R:34`); survival only | none | `test-plan-restored.R` (`skip_on_cran`) | none | **1 / 2** (`quarto/calibration/fs_fdr_report_smoke.qmd:109`, `:168`) (0) | 0; 0 | No complete equivalent. Partial: `fpr_calibration()` with `null_fn = function(seed) simulate_from_dgm(<flex DGM>, ...)`; the caller builds both nulls (`generate_aft_dgm_flex(k_inter = 0)` and a `find_k_inter_for_target_hr(1.0)` DGM), passes `c1`/`c2` and `fs_params`, and computes B(g) from per-replicate MR. Equality of its rate with A not verified. |
| `generate_gbsg_bootstrap_general` | `R/synthetic_data_perturbation.R:390` | yes (`NAMESPACE:206`) | yes (`v0.4.0:NAMESPACE:191`), **broken there**: it loaded with `data("gbsg", package = "survival", envir = environment())`, which leaves `gbsg` NULL, so the call stopped ("'data' must be a data frame"); fixed after the release in `67ce5b41` (A-004) | `gbsg <- survival::gbsg`; fixed variable lists (continuous `age, size, nodes, pgr, er, rfstime`; categorical `meno, hormon, status`; ordinal `grade`; id `pid`); default `n = 686`; no data argument | none | `test-gbsg-bootstrap-general.R`; `test-synthetic-data-documented-behaviour.R` (new, T1-T2) | none | 0 / 0 (0 / 0) | 0; 0 | `generate_bootstrap_synthetic(data = survival::gbsg, ...)` with the wrapper's lists: the same data (`identical()`, T2). |


**Considered and not listed:** `cut_numeric` (a general quantile coder); `cut_size` (general, GBSG default; listed below); `.rename` in `setup_gbsg_dgm`; the `fs_fdr_report` helpers (`.fs_fdr_one_rep`, `.fs_fdr_run_null`, `.fs_fdr_hr_row`, `print.fs_fdr_report`), reached only from it (`.fs_fdr_hr_row` reads `gbsg_dgm` fields first; its fallbacks `hr$HR_H`/`hr$HR_Hc` name fields no builder sets, so on a flex DGM those columns would be NA); `generate_bootstrap_synthetic`, `detect_variable_types`, `find_k_inter_for_target_hr`, `sensitivity_analysis_k_inter`, `find_quantile_for_proportion` (take `data`; GBSG only in examples), `calibrate_k_treat` (data inside `base_args`) and `calibrate_cens_adjust` (takes a `dgm`); `fs_dgm_feasibility` (rejects any non-`glm_dgm`). The constants `SEED_BASE`, `DAYS_PER_MONTH`, `DEFAULT_N_SUPER` and `.FS_GBSG_CONFOUNDERS` are not functions.

### 7.1 General functions that carry a GBSG leftover

Line numbers at the branch head.
- **`generate_bootstrap_with_noise()`** (exported; identical at `v0.4.0`): default `id_var = "pid"` (`R/synthetic_data_perturbation.R:539`); special case `var == "grade"` for adjacent-level weighting (`:653`; documented in the Part 3 note).
- **`run_simulation_analysis()`** (exported): default `confounders_base = c("v1", ..., "v7")` (`R/run_simulation_analysis.R:398`), the GBSG wrapper's analysis factors; its `@param` states no default. A DGM without `v1..v7` needs `confounders_base` passed.
- **`fs_to_guohe()`** (exported): defaults `time = "time_months", event = "status", treatment = "hormon"` (`R/fs_to_guohe.R:339`; help `:280` "The defaults match the GBSG analyses"). It stops with a clear message when they are absent.
- **`plot.subgroup_sims_summary()`** (S3 method): for `metric = "hr"` the default x label is "Hazard Ratio (Cox, stratified by grade)" (`R/plot_subgroup_sims.R:296`), but `run_subgroup_sims()`'s default fit is unstratified. A label leftover only.
- **`cut_size()`** (internal, `man/cut_size.Rd`): default `breaks = c(20, 50)`, GBSG tumour size in mm (`R/sim_aft_gbsg.R:69`); its only caller is the GBSG builder.
- **Fallbacks for the `gbsg_dgm` object layout** (no GBSG data; they work on any DGM): `compute_dgm_cde()` tries `df_super_rand` first and `flag.harm` first (`R/oc_analyses.R:134`, `:145`); `interpret_estimation_table()` runs its CDE fallback only when `dgm$df_super_rand` exists (`R/simulation_tables.R`), so never for a flex DGM; `get_dgm_hr()`, `invert_dgm_hrs()`, `build_classification_table()`, `build_estimation_table()` fall back to `gbsg_dgm` fields.
- **`.extract_grf_estimates_gen()`** (internal, `R/run_simulation_analysis.R:1314`): falls back to the legacy GBSG simulator's output names `"y.sim"`/`"event.sim"` when present (`:1449`, `:1451`).
- **Non-function couplings:** `R/globals.R:159-167` (`z_age` … `z_hormon`, the flex column names for GBSG); `R/globals.R:170` block "GBSG dataset column names (used in examples/tests)" and `:197` block "GBSG DGM internal variables" (its header names a file `sim_aft_gbsg_refactored.R` that does not exist); `.FS_GBSG_CONFOUNDERS` (`R/fs_fdr_report.R:34`).


### 7.2 Help-page examples that load GBSG in a way that fails as written

Measured in a fresh `Rscript --vanilla` (R 4.6.2, survival 3.8-12): `data(gbsg, package = "survival")` warns "data set 'gbsg' not found" and loads nothing (in survival's `Rdata.rds` `gbsg` is filed under `cancer`); bare `gbsg` then fails with "object 'gbsg' not found"; `library(survival); data(gbsg)` warns but finds `gbsg` through the attached package. `survival` is in Imports, so `library(forestsearch)` does not attach it. All examples below are in `\dontrun`. Line numbers at the start SHA.

- **`data(gbsg, package = "survival")` then bare `gbsg`** — fixed in Part 3: `R/synthetic_data_perturbation.R:438` (`generate_bootstrap_with_noise`), `R/pareto_frontier_table.R:124` (`pareto_frontier_table`; `:125` is also a placeholder `forestsearch(gbsg, ...)`).
- **`library(survival); data(gbsg)`** — works with the warning; replaced in Part 3 as directed: `R/calibrate_k_treat.R:54-55`.
- **Bare `gbsg` with no loading line** — fail; not edited (outside Part 3's list; OPEN ITEMS). Each `forestsearch(gbsg, ...)` here also names `"grade3"`, which is not a `survival::gbsg` column:
  - `R/bootstrap_summaries_helpers.R:28` (`format_bootstrap_table`)
  - `R/forestsearch_cross_validation.R:2002` (`print.fs_kfold`), `:2031` (`print.fs_tenfold`)
  - `R/plot_sg_distribution.R:76` (`plot_sg_distribution`; `results_alt` also undefined)
  - `R/plot_sg_results.R:866` (`print.fs_sg_plot`), `:912` (`plot.fs_sg_plot`)
  - `R/render_forestplot.R:372`, `:375` (`print.fs_forestplot`), `:410`, `:413` (`plot.fs_forestplot`)
  - `R/run_subgroup_sims.R:660` (`run_subgroup_sims`; `dgm_uniform`, `subgroups`, `cal_uniform`, `k_treat_uniform` also undefined)
  - `R/subgroup_consistency_helpers.R:1031` (`sg_consistency_out`), `:1361` (`evaluate_subgroup_consistency`), `:1645` (`evaluate_consistency_twostage`); `:1223` has the pattern but its block is `@noRd`.
- **Calls to the functions removed in Part 1** — fixed in Part 1: `R/oc_analyses.R:123` (`compute_dgm_cde`), `:413` (`summarize_simulation_results`; its next line, `run_simulation_analysis(dgm, nsim = 10)`, fails on its own: there is no `nsim` argument), `R/sim_aft_gbsg.R:990`, `:995` (`calibrate_k_inter`), `:1153` (`print.gbsg_dgm`).
- **A GBSG frame the example never builds**: `R/fs_fdr_report.R:148`, `:159`, `:164` (`gbsg_df`).
- **Loads GBSG correctly but names a column it lacks**: `R/summary_utility_functions.R:742-745` (`sg_tables`): `library(survival)` then `confounders.name` with `"grade3"`; by code reading `forestsearch()` stops ("confounders.name references columns not in df.analysis"); not run.

## 8. Part 6 — version, NEWS, code identity, suite, check

### 8.1 Version and NEWS (Step 1)

`DESCRIPTION` 0.4.1.9008 → 0.4.1.9009; `NEWS.md` heading → `# forestsearch 0.4.1.9009 (development)`. New subsections in the development block: **Removed** (the four functions), **Superseded** (`setup_gbsg_dgm()`, `calibrate_k_inter()`, `generate_gbsg_bootstrap_general()`, and the censoring-step note), **Documentation** (every help-page correction of Parts 1-4), and an **Internal** entry for the new test file. Commit `aa667898`.

### 8.2 Code-identity gate (Step 2) — **PASS**

`dev/reports/REPORT_cran_prep_docs_legacy_2026-11-12_code_identity.R` exports `R/` at `d93e6621` and at `aa667898` with `git archive` into `/tmp/cran_prep_docs_legacy_2026-11-12/code_identity/`, parses every file with `parse(keep.source = FALSE)` and compares the deparsed top-level expressions by assigned name. Output: `…_code_identity.txt`.

| Gate | Result |
|---|---|
| Names only at the start = the functions Part 1 removed | `create_gbsg_dgm`, `simulate_from_gbsg_dgm`, `.simulate_from_gbsg_dgm_`, `get_dgm_with_output` — **PASS** (no helper was removed) |
| Names only at the head = none | 0 — **PASS** |
| Names in both with a different deparsed body = none | 0 of 798 — **PASS** |
| Other top-level differences only `globalVariables()` removals | 19 non-assignment expressions at each end, 0 differ — **PASS** (no globals name removed) |

125 `R/` files at both ends; 829 → 825 top-level expressions. Negative control: the same script on `a977f3d7` → `55e59792` (the RNG-restore commit of task part 2) reports 15 changed bodies and 3 new names and FAILS, so it detects body changes.

### 8.3 Help pages (Step 3) — **PASS**

`git diff --stat d93e6621 -- man/` lists 20 pages, each a page of a function this task removed or edited:

| Part | Pages |
|---|---|
| 1 removed | `create_gbsg_dgm.Rd`, `simulate_from_gbsg_dgm.Rd`, `get_dgm_with_output.Rd` |
| 1 edited (references) | `calibrate_k_inter.Rd`, `compute_dgm_cde.Rd`, `print.gbsg_dgm.Rd`, `run_simulation_analysis.Rd`, `setup_gbsg_dgm.Rd`, `summarize_simulation_results.Rd`, `validate_k_inter_effect.Rd` |
| 2 | `setup_gbsg_dgm.Rd`, `calibrate_k_inter.Rd`, `fs_fdr_report.Rd`, `validate_k_inter_effect.Rd` |
| 3 | `generate_bootstrap_with_noise.Rd`, `generate_bootstrap_synthetic.Rd`, `generate_gbsg_bootstrap_general.Rd`, `pareto_frontier_table.Rd`, `calibrate_k_treat.Rd` |
| 4 | `dina_subgroup.Rd`, `dina_subgroup_refit.Rd`, `forestsearch.Rd`, `subgroup.consistency.Rd`, `fs_fdr_report.Rd` |

20 files, +219 / −328 lines.

### 8.4 Full suite (Step 4)

| Run | Baseline (`REPORT_cran_prep_seed_legacy_2026-11-11.md` §3.5, not re-run) | This branch | Result |
|---|---|---|---|
| Full suite (`devtools::test()`, `NOT_CRAN=true`, in the repository) | 648 blocks, 88 files, 8,471 expectations, 0 failures, 3 skips, 38 warnings; 686 s | **652 blocks, 89 files, 8,478 expectations, 0 failures, 0 errors, 3 skips, 38 warnings**; wall **694 s** (10:57:30–11:09:04; cap 30 min) | **PASS** |
| Difference | — | +4 blocks, +1 file, +7 expectations: exactly `test-synthetic-data-documented-behaviour.R` (T1-T4: 1 + 1 + 3 + 2). Part 1 removed no test (none tested only a removed function). The skips are the same 3 (`test-cross-outcome-parity.R` ×2, `test-cv-no-subgroup-edges.R`); the 38 warnings are unchanged (`test-cross-outcome-parity.R` 12, `test-return-shape-contract.R` 9, `test-mrct_inject_dgm.R` 6, `test-downstream-consumers.R` 5, `test-get-fsdata-diagnostics.R` 4, `test-input-validation.R` 1, `test-subgroup-sims.R` 1). | |

- **Tree tested.** The suite started at 10:57:30 on HEAD `386155f3` with Part 4's edits (first draft) and the version/NEWS edits in the working tree. The R/ code then was the final code: Part 4's later corrections and Part 6 changed comments, `DESCRIPTION` and `NEWS.md` only, and the code-identity gate (§8.2) shows no expression differs between the start SHA and `aa667898` apart from the four removals. Tests load the package with `load_all()`, so comment text cannot affect them.
- **Contention.** `check()` ran at the same time (from 10:59), so both wall clocks include contention, as in the two earlier reports.


### 8.5 `devtools::check(document = FALSE, error_on = "never")` (Step 5)

| Run | Baseline (`REPORT_cran_prep_seed_legacy_2026-11-11.md` §3.5, not re-run) | This branch (`aa667898`) | Result |
|---|---|---|---|
| `devtools::check(document = FALSE, error_on = "never")` | 0 errors / 0 warnings / 0 notes; 1,121 s | **0 errors / 0 warnings / 0 notes**, Status OK, `forestsearch 0.4.1.9009`; wall **1,035 s** (10:59:32–11:16:47; R's Duration 15 m 26 s; examples 36 s, `--run-donttest` 63 s, tests 11 min, vignette rebuild 75 s; cap 60 min) | **PASS** |

- **Check route** (`REPORT_cran_prep_2026-11-11.md` §6, with the `-z` correction of `REPORT_cran_prep_seed_legacy_2026-11-11.md` §3.5): an export of the working tree at `aa667898`, `git ls-files -z -co --exclude-standard`, 9,682 files listed and 9,682 copied, into `/tmp/cran_prep_docs_legacy_2026-11-12/check/forestsearch`; `check()` ran there. Every file changed since the start SHA was `cmp`-identical to the repository's except one: `dev/reports/REPORT_cran_prep_docs_legacy_2026-11-12_code_identity.txt`, which the code-identity script was writing while the export ran; `dev/` is build-ignored (`.Rbuildignore` `^dev$`), so the check is unaffected. The local `legacy_vignette/` (git-excluded) was not exported, so the known NOTE did not arise.
- No repair was needed (the task allowed at most two for help-page or roxygen syntax errors).

## 9. Final gates (Step 8)

Checked on the tree of the report commit (in the final printout of this run):

- **`git status --porcelain` is empty** after the report commit.
- **The two sibling repositories are unchanged.** fs-post-selection `main` at `5d165acc`, fs-glms-interpretable `main` at `dc94a17e`, each with an empty `git status --porcelain` and no file newer than the start commit's task file — as at Step 0.
- **`git diff d93e6621 --stat` lists only files this report names, and nothing under `vignettes/` or `quarto/`**: the 38 package files of Parts 1-4 and 6 plus the report and its two code-identity files (Appendix).
- **Nothing pushed or merged**; no `fetch`/`pull`. The report is copied to `~/Downloads/REPORT_cran_prep_docs_legacy_2026-11-12.md`.

## OPEN ITEMS

1. **The supplement and other prose still name `create_gbsg_dgm()`.** fs-post-selection `manuscript_rev1/postselection_subgroup_inference_supplement.qmd:1184` (and its `.tex:3296`, `.pdf`) says the simulation DGMs are "`create_gbsg_dgm()` and `generate_aft_dgm_flex()`"; the function no longer exists at 0.4.1.9009. Same sentence in `manuscript/…supplement.qmd:1182`, the JRSSB submission and fs-glms-interpretable `dev/` copies, and in `vignettes/articles/biomarker_effects.qmd:66` (+ `.html:2413`). Read-only here; listed in §3.5. A rewording could say `setup_gbsg_dgm()`.
2. **Help text left as it was, though inaccurate** (the task said "everything else on the page stays", or the site was outside the listed ones):
   - `generate_bootstrap_with_noise()` Details "Preserves integer type for variables that appear to be integers" — an integer column treated as continuous comes back double (measured); `@return` "same structure as the input data" likewise; `@note` "Missing values (NA) … are not imputed" — a selected NA in a categorical column is replaced by a level, and an ordered factor (or `grade`) with more than two levels and an NA stops with "replacement has length zero" (measured).
   - `generate_bootstrap_with_noise()` has no `is.numeric` guard: a factor named in `continuous_vars` stops with "Calling var(x) on a factor x is defunct" (measured); `generate_bootstrap_synthetic()` has one (start `:113`).
   - (start lines) `R/dina_subgroup.R:402` `#` comment "log-IRR for poisson"; `R/forestsearch_main.R:2457-2458` `#` comment "(event count threshold over a window)" for `dmin.grf`.
   - `forestsearch()` `@param consistency_method` (start `R/forestsearch_main.R:1045`) and `R/subgroup_consistency_helpers.R:336-338` say "propensity-adjusted effects" fall back to `"split"`; `.glm_resample_supported()` checks only `effect_measure` (`:343-345`), so they do not.
   - `R/calibrate_k_treat.R:9-11` calls `calibrate_k_inter()` the flex "natural sibling" of `calibrate_k_treat()`; it builds only GBSG-wrapper DGMs. The flex calibrator is `find_k_inter_for_target_hr()`.
   - `find_k_inter_for_target_hr()`: `@return` lists a `method` element the returned list lacks (`R/find_k_inter_main.R:43` vs `:166-173`); `@param continuous_vars` says "standardized", but flex's default `standardize = FALSE` applies.
   - `vignettes/articles/paper_simulations.qmd:388-445` calls an `eval: false` `generate_aft_dgm_flex()` block "equivalent" to `setup_gbsg_dgm()` (it uses median splits; the wrapper uses quartile codes), and `:459-463` says `setup_gbsg_dgm` "remains the correct entry point", which now conflicts with its Superseded section. Vignettes are read-only here.
   - `compute_dgm_cde()` `@param dgm` still says the object "Must contain `df_super_rand`"; the code falls back to `df_super` (`R/oc_analyses.R:134`), so a `generate_aft_dgm_flex()` object (which Part 1's reference fix now names first) works. Part 1 changed only the example list in that sentence.
   - Twelve `\dontrun` examples (14 lines; §7.2) use bare `gbsg` with no loading line, plus the same pattern at `R/subgroup_consistency_helpers.R:1223`, which is on no help page, and `summarize_simulation_results()`'s example passes `nsim = 10`, which `run_simulation_analysis()` does not take.
3. **Behaviour documented, not changed (Larry's decision):** the factor-coded ordinal NA in `generate_bootstrap_synthetic()`, its end-level `sample()` draw, and its loop over `cat_vars` only. The factor-to-NA behaviour is stated on the help page but not pinned by a test (the task's T1-T4 do not cover it).
4. **`dmin.grf < 0` edge case (code reading, not run):** with `sg_focus = "hrMaxSG"`/`"hrMinSG"` an empty band gives a zero-row frontier selection (`R/grf_subgroup_labels.R:380-391`) and `.grf_sg_def_from_candidate()` then errors (`:405`). The new `dmin.grf` text is exact for `dmin.grf >= 0`, which includes the default.
5. **`setup_gbsg_dgm(cens_type = "uniform")`** stores `model_params$censoring = list(type = "uniform", mu = NULL, tau = NULL)` and `simulate_from_dgm()` then stops with "invalid arguments" (p2 verifier run R3). Pre-existing; not in scope.
6. **`create_gbsg_dgm`'s `@details`** (subgroup definition z1 = 1 & z3 = 1, the model specification and the k_inter algebra) went with its help page. `setup_gbsg_dgm()` now documents its own arguments and return fields, but not that model description; it can be carried over if wanted.
7. **`t.sim`** in `R/globals.R:208` is now referenced by no package code (it was used only inside the removed `.simulate_from_gbsg_dgm_`). Kept ("when in doubt, keep it"); harmless.
8. **Reverse dependencies** of the two removed exports were not checked (both were in CRAN 0.4.0).
9. **Not rebuilt:** the pkgdown site (the `_pkgdown.yml` reference index lost three entries; nothing else).

## Appendix — files changed since the start SHA

From `git diff --name-status d93e6621` on the tree of the report commit (A = added, M = modified, D = deleted). None is under `vignettes/` or `quarto/`. The task file itself is the start commit.

| Status | Path | Part |
|---|---|---|
| M | `DESCRIPTION` | 6 |
| M | `NAMESPACE` | 1 |
| M | `NEWS.md` | 6 |
| M | `R/calibrate_k_treat.R` | 3 |
| M | `R/dina_subgroup.R` | 4 |
| M | `R/dina_subgroup_refit.R` | 4 |
| M | `R/forestsearch_main.R` | 4 |
| M | `R/fs_fdr_report.R` | 2, 4 |
| M | `R/oc_analyses.R` | 1 |
| M | `R/pareto_frontier_table.R` | 3 |
| M | `R/run_simulation_analysis.R` | 1 |
| M | `R/setup_gbsg_dgm.R` | 1, 2 |
| M | `R/sim_aft_gbsg.R` | 1, 2 |
| M | `R/subgroup_consistency_main.R` | 4 |
| M | `R/synthetic_data_perturbation.R` | 3 |
| M | `README.md` | 1 |
| M | `_pkgdown.yml` | 1 |
| M | `man/calibrate_k_inter.Rd` | 1, 2 |
| M | `man/calibrate_k_treat.Rd` | 3 |
| M | `man/compute_dgm_cde.Rd` | 1 |
| D | `man/create_gbsg_dgm.Rd` | 1 |
| M | `man/dina_subgroup.Rd` | 4 |
| M | `man/dina_subgroup_refit.Rd` | 4 |
| M | `man/forestsearch.Rd` | 4 |
| M | `man/fs_fdr_report.Rd` | 2, 4 |
| M | `man/generate_bootstrap_synthetic.Rd` | 3 |
| M | `man/generate_bootstrap_with_noise.Rd` | 3 |
| M | `man/generate_gbsg_bootstrap_general.Rd` | 3 |
| D | `man/get_dgm_with_output.Rd` | 1 |
| M | `man/pareto_frontier_table.Rd` | 3 |
| M | `man/print.gbsg_dgm.Rd` | 1 |
| M | `man/run_simulation_analysis.Rd` | 1 |
| M | `man/setup_gbsg_dgm.Rd` | 1, 2 |
| D | `man/simulate_from_gbsg_dgm.Rd` | 1 |
| M | `man/subgroup.consistency.Rd` | 4 |
| M | `man/summarize_simulation_results.Rd` | 1 |
| M | `man/validate_k_inter_effect.Rd` | 1, 2 |
| A | `tests/testthat/test-synthetic-data-documented-behaviour.R` | 3 |
| A | `dev/reports/REPORT_cran_prep_docs_legacy_2026-11-12.md` | 6 |
| A | `dev/reports/REPORT_cran_prep_docs_legacy_2026-11-12_code_identity.R` | 6 |
| A | `dev/reports/REPORT_cran_prep_docs_legacy_2026-11-12_code_identity.txt` | 6 |
