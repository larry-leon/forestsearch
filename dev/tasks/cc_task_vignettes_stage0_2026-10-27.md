# CC task — vignettes stage 0: Table 1 provenance, package facts, and the timing ladder (2026-10-27)

**Repo:** `larry-leon/forestsearch` — run CC from this repo's root.
**Kind:** read-only inventory (forestsearch at HEAD; the sibling repo fs-glms-interpretable read-only) plus one small measurement (the timing ladder of §3). **No edit to `R/`, `vignettes/`, `man/`, `DESCRIPTION`, `NEWS.md`, `_pkgdown.yml`, `.Rbuildignore`, `tests/`.** Writes: this task doc under `dev/tasks/`; the record directory `dev/reports/vignettes_stage0_2026-10-27/`; the report `dev/reports/REPORT_vignettes_stage0_2026-10-27.md`.
**Compute:** the ladder of §3 — at most nine single-process `Rscript` runs, `plan = "sequential"`, each under a hard cap of 20 minutes, the ladder under a cumulative cap of 3 hours. Approved by Larry 2026-10-27 ("proceed as you recommend"). No campaign, nothing parallel.
**Unattended.** Gates stop the task; never ask, never work around. No `git fetch` / `pull` / `push`. Explicit paths on every `git add`; never stage pre-existing untracked files. No `R CMD check`, no vignette build, no test suite, no `load_all()`.

---

## Why this task exists (context for the record)

- The 0.4.1 release (target 10 Nov 2026) carries revised vignettes. The new primary vignette, `vignettes/forestsearch.qmd`, is to **run** the analysis behind Table 1 of León & Anderson (2026b), *Interpretable subgroup identification across survival and generalized-linear outcomes*, arXiv:2610.28379 — the three-identifier GBSG analysis (FS, DINA, GRF) with the selection-adjusted products — live in the render, for reproducibility. That vignette's task is written on the facts this task records: the exact calls that produced Table 1, the package's current formals and result-object layout, and the measured single-core cost at the manuscript's draw counts and at reduced ones.
- Standing decisions (Larry, 2026-10-27), recorded here, **none executed by this task**: D1 the current `vignettes/forestsearch.Rmd` moves as-is to `vignettes/articles/forestsearch-2024.Rmd` with a pointer note (no .Rmd → .qmd conversion). D2 the new primary vignette is `.qmd`. D3 all three identifiers if the ladder shows they fit the build budget; the FS column guaranteed. D4 the vignette runs the analyses — nothing precomputed is shipped; whether the CRAN build runs the manuscript's draw counts outright or reduced counts with a `params` switch is decided on this ladder. D5 `survival-post-selection.qmd` is kept and aligned (p\* 0.90, 2026b added, final arXiv wording). D6 the ACTG175 GLM document is a pkgdown article for 0.4.1. D7 version 0.4.1.
- Larry's clarification (2026-10-27), to be verified from source in §2.1 and stated in the vignette later: the analyses use the **closed-form calculation** of the splitting-consistency probability; random half-splitting (`fs.splits`) is only the fallback in the rare cases where the closed form is not evaluable. So the search is deterministic given the data and the only Monte Carlo in Table 1 is the multiplier draws and the field's outer/inner draws.

**Table 1 of 2026b (p. 19), the comparison target of §3** — hazard-ratio scale; unadjusted: two-sided 95% Wald; field and field-s: one-sided at nominal 0.95; Bonferroni pair: region / complement, joint 0.95. Whole trial: 686 patients, 299 events, HR 0.695 (0.544, 0.888).

| | FS | DINA | GRF |
|---|---|---|---|
| Region | ER 0 and PR ≤ 26 | grade 3 and PR ≤ 20 | ER 0 |
| N (events) | 75 (41) | 113 (61) | 82 (45) |
| Unadjusted, 95% CI | 2.222 (1.169, 4.221) | 1.315 (0.746, 2.317) | 1.951 (1.037, 3.673) |
| Corrected | 1.276 | 1.076 | 1.111 |
| Field, lower bound | 0.634 | 0.638 | 0.510 |
| Complement N (events) | 611 (258) | 573 (238) | 604 (254) |
| Unadjusted, 95% CI | 0.605 (0.465, 0.788) | 0.641 (0.489, 0.839) | 0.615 (0.472, 0.801) |
| Corrected | 0.645 | 0.658 | 0.663 |
| Field-s, upper bound | 0.804 | 0.842 | 0.846 |
| Bonferroni pair | 0.509 / 0.840 | 0.572 / 0.887 | 0.404 / 0.874 |

Specification per §6.1 / S3.1 / S5 of the paper: seven covariates (age, menopausal status, tumour size, grade, positive nodes, PR, ER); ten candidate cuts on each receptor with zero forced in as a cut of its own; candidates are conjunctions of at most two cuts, more than 60 patients, FS additionally ≥ 10 events per arm; screening and per-split thresholds both HR 1.0; consistency screen p\* = 0.90; rule max_𝒜 N(ε) at ε = 0.20; 5,000 multiplier draws; the field's Monte Carlo at 1,000 outer × 500 inner draws. Supplement Table S13, FS funnel: 2,211 enumerated, 1,410 fittable, 121 clearing the screening threshold, 9 passing the consistency screen, 3 on the frontier, 1,744 in the re-selection family.

---

## 0. First action — branch, task doc, provenance — GATE

1. `git status --short`: record every pre-existing untracked or modified path; never stage them. STOP if any touches `dev/tasks/`, `dev/reports/vignettes_stage0_2026-10-27/` or `dev/reports/REPORT_vignettes_stage0_2026-10-27.md`.
2. *GATE:* current branch is `master` and `git merge-base --is-ancestor 174778d09aea19311d402b1c52f11ad9b55cfb86 HEAD` succeeds (the handoff commit is in history).
3. `git switch -c feature/vignettes-0.4.1` — the workstream's branch; every vignette task that follows runs on it; Larry checks and merges at the end of the workstream.
4. Copy this file from `~/Downloads` to `dev/tasks/cc_task_vignettes_stage0_2026-10-27.md` (the downloaded name may arrive with its hyphens stripped; the committed name is this one — report both). `git add` that path only; commit `dev: task — vignettes stage 0 (2026-10-27)`.
5. Provenance to `dev/reports/vignettes_stage0_2026-10-27/provenance.txt`: hostname, `pwd`, HEAD, branch, `nproc`, `R.version.string`, `.libPaths()[1]`, `packageDescription("forestsearch")$Version` and `$Built`, `packageVersion()` of `survival`, `grf`, `future`, `doFuture`, `quarto` (R package, if installed), `knitr`, `rmarkdown`, and `quarto --version` from the shell.
   *GATE:* installed forestsearch is `0.4.0.9000` **and** the installed build is HEAD's tree: `cmp NEWS.md "$(Rscript -e 'cat(system.file("NEWS.md", package = "forestsearch"))')"` is silent (NEWS last changed at `936b6734`, before the install from `30441802`). If `NEWS.md` is not installed, compare the installed and repo `DESCRIPTION` after dropping the `Built:` / `Packaged:` / `Remote*:` lines. A mismatch is a STOP: the ladder would measure a build that is not HEAD.

---

## 1. Inventory A — where Table 1 came from (read-only, sibling repo)

Location: `~/Documents/GitHub/fs-glms-interpretable/` (found there 2026-10-26; `~/GitHub/fs-glms-interpretable/` does not exist on this machine). *GATE:* the directory exists and is a git repo. Record `git -C <dir> rev-parse HEAD` and `git -C <dir> status --short` (expected HEAD `a772a2fe`; a different HEAD is recorded, not a stop). **Nothing in that repo is edited, built or rendered.**

1. **Listing.** `quarto/gbsg/` and `quarto/actg175/`: every `.qmd`, its title, whether a rendered `.html` is committed, and its `_payloads/` entries (name, size, `built_at` and `forestsearch_version` where the payload carries them). One table per directory in the record.
2. **Locate Table 1's producers.** Grep the `.qmd`, rendered `.html` and payload directories under `quarto/gbsg/` for the published values — region: `1.276`, `0.634`, `1.076`, `0.638`, `1.111`, `0.510`; complement: `0.804`, `0.842`, `0.846`; Bonferroni: `0.509`, `0.840`, `0.572`, `0.887`, `0.404`, `0.874`. For `.rds` payloads, `readRDS` each and grep the `capture.output(str(x, list.len = 1e4))` text. Record which document and payload carries which column, the chunk names, the render dates, and whether a single document ("trio" or similar) assembles all three. *GATE:* the FS column is located (`analysis_gbsg_mr.qmd` is the candidate the 2026-09-22 declaration-calibration task names: c1 = c2 = 1.00, p\* 0.90, region `{er <= 0} & {pgr <= 26}`, N 75, HR 2.22); if it cannot be located, STOP and report. DINA and GRF columns not located is a finding, not a stop.
3. **Transplant verbatim** into `dev/reports/vignettes_stage0_2026-10-27/table1_calls/`, one file per column — `fs_call.R`, `dina_call.R`, `grf_call.R` — the data-preparation chunk(s) and the `forestsearch()` call that produced that column, with a header comment giving the source path, the chunk names and the sibling repo's HEAD. If a column was produced by a chain other than one `forestsearch()` call (a standalone identifier call followed by an MR function, say), transplant that chain. Do not edit the transplanted code beyond the header comment. If the document prepares the data from a file rather than `survival::gbsg`, note it in the header; the ladder uses `survival::gbsg` with the same preparation (GATE A of §3 covers the consequence).
4. **Settings table**, one row per column, values exactly as written in the document (where the document sets a value in a settings chunk and passes a name, quote both): `subgroup_method` / identifier; the consistency-method argument and its value (name per §2.1); `fs.splits`; `pconsistency.threshold`; `hr.threshold`; `hr.consistency`; `sg_focus`; `selection_rule`; `effect_neighborhood`; `n.min`; `d0.min`; `d1.min`; `maxk`; `conf.cont_jcuts`; `conf_force`; `collapse_cuts`; `seedit`; `mr_inference` and every element of `mr_inference_args` (draws, field outer/inner draws, complement scaling, Bonferroni, `declaration_c0`, anything else); `parallel_args`; every identifier-specific argument (`dina_args`, `dmin.grf`, `tune_grf`, `grf_*`, …); `max.minutes`, `vi.grf.min`, `max_n_confounders` if passed (all three are deprecated and inert or warn in the installed build — record, do not edit).
5. From the same documents or their reports: the render wall-clock if stated, and the `n_cores` / plan used, so the ladder's single-core numbers can be read against them.

---

## 2. Inventory B — package facts from source (forestsearch at HEAD; the installed build for run-time facts)

Record in `dev/reports/vignettes_stage0_2026-10-27/package_facts.md`, quoting file and line for every code fact.

1. **The consistency-evaluation method — GATE.** Find the formal of `forestsearch()` (and of the consistency routine it calls) that selects the closed-form calculation of p_cons versus random half-splitting: grep `R/` for the candidates `consistency_method`, `closed`, `analytic`, `fs.splits`. Record: the formal's name, allowed values and default in `forestsearch()`; the exact condition under which random splitting is used as a fallback (quote the lines); whether `fs.splits` has any effect when the closed form is evaluable; what the result object records about the method used per candidate (a flag, a count of fallbacks) and the element path. *GATE:* a closed-form path exists in the installed build and is reachable from `forestsearch()`. If not, STOP and report — the vignette's text depends on it.
2. **Formals.** `formals(forestsearch::forestsearch)` — every formal with its default (`deparse`d), to `formals_forestsearch.txt`; `names(formals())` of the MR routine (`forestsearch:::fs_mr_inference` or its current name) and the documented elements of `mr_inference_args` from the roxygen in `R/forestsearch_main.R`.
3. **MR products and where they live.** After rung S1 of §3 has run, `str()` its result and record the element path of every Table 1 quantity: unadjusted HR and 95% CI in Ĥ and in Ĥᶜ; N and events in each; the corrected (two-term) estimate in each; the field lower bound on β(Ĥ); the field-s upper bound on β(Ĥᶜ); the Bonferroni pair; the MR family size; the draw counts actually used (MR, field outer, field inner); the seed — and which of these `print.forestsearch`, `summary`, and any `fs_*` reporting helper prints. Record whether `declaration_c0` (κ̂, implied p\*, the diagnostic) is reachable through `forestsearch(mr_inference_args = …)`, and what it returns.
4. **Build facts**, copied verbatim into the record: `DESCRIPTION` (whole file); `.Rbuildignore`; `_pkgdown.yml`; the `vignettes/` listing (recursive) with each file's YAML header and `%\Vignette*` lines; whether `grf`, `quarto`, `knitr`, `rmarkdown`, `speff2trial`, `policytree` are in Imports or Suggests; any file under `tests/testthat/` that references `vignettes/` or `inst/doc` (none expected).

---

## 3. Compute — the timing ladder (approved 2026-10-27)

Purpose: on one core, measure (a) the FS search at Table 1's specification, (b) FS + MR at three draw counts, (c) DINA and GRF at the bottom and top rungs — so the vignette task can choose, on measured numbers, between the manuscript's settings outright and reduced draws with a `params` switch.

**Setup**
- `dev/reports/vignettes_stage0_2026-10-27/ladder.R`, one rung per invocation; driver `run_ladder.sh` runs the rungs in the order below as `timeout 1200 /usr/bin/time -v Rscript ladder.R <rung> > logs/<rung>.log 2>&1`, appends elapsed per rung to `ladder.csv`, and stops launching new rungs once the cumulative elapsed exceeds 3 h or a rung is killed by `timeout` (record "exceeded cap"; a cap exceedance is a finding, not a task failure — except at S0, GATE A). If `/usr/bin/time` is absent, record `system.time()` elapsed and `sum(gc()[, "max used"])` from inside R and say so.
- Every rung uses the **installed** package (`library(forestsearch)`), the transplanted calls of §1.3 **verbatim**, with exactly these overrides: `parallel_args = list(plan = "sequential")`, `quiet = TRUE`, `details = FALSE`, `plot.sg = FALSE` (and `plot.grf = FALSE` where such a formal exists); `mr_inference` and the draw arguments as each rung states, under the argument names §2.1–2.3 established. "Package default field" means the field's outer/inner draw arguments left unset. `saveRDS` the resolved argument list per rung (minus the data frame). Data: `survival::gbsg` prepared exactly as the transplanted chunk prepares it; no file reads.
- Warnings per rung are captured in the log and listed in the report (deprecation warnings are expected if a transplanted call passes `max.minutes` non-default, `vi.grf.min` or `max_n_confounders`); do not edit the transplanted call to silence them.

**Rungs, in order**

| rung | identifier | `mr_inference` | MR draws | field outer × inner |
|---|---|---|---|---|
| S0 | FS | FALSE | — | — |
| S1 | FS | TRUE | 500 | package default |
| S1r | FS | TRUE | 500 | package default — a repeat of S1, the determinism check |
| S2 | FS | TRUE | 2000 | package default |
| S3 | FS | TRUE | the transplanted call's own MR settings (expected 5000 and 1000 × 500; if different, record and keep the transplanted values — the comparison with Table 1 is then like for like) |
| D1 | DINA | TRUE | 500 | package default |
| G1 | GRF | TRUE | 500 | package default |
| D3 | DINA | TRUE | the transplanted call's own — only if S3 completed under cap |
| G3 | GRF | TRUE | the transplanted call's own — only if S3 completed under cap |

**Per rung, record** (to `ladder.csv` and a per-rung `.rds`): elapsed wall-clock and max RSS; the selected region (`sg.harm`); N and events in Ĥ and Ĥᶜ; the candidate-flow counters (`find.grps$filter_counts` and the `grp.consistency` counts, as the 2026-10-26 stage-0 scripts read them); p_cons for every candidate passing the screen, with the closed-form/fallback flag per candidate and the count of fallbacks (§2.1's element); the MR family size; the Table 1 quantities at §2.3's paths.

**GATE A (stop-on-failure, at S0):** the selected region is `er <= 0` & `pgr <= 26` with N = 75 and 41 events. Otherwise STOP and report, with the full `sg.harm`, the counters and the p_cons table: the transplant or the installed build does not reproduce the manuscript's search, and the chat must know before anything is written. Recorded beside it, not gated: the counters against Table S13 (2,211 / 1,410 / 121 / 9 / 3 / 1,744) — the definitions may differ; show them side by side with one line on any definitional difference found in the code.

**Findings, not gates**
- S1 vs S1r: `all.equal()` on the selected region and every Table 1 quantity — identical or not, differing elements listed.
- S3 vs Table 1, FS column: the nine quantities, published vs rung, with the difference on the log scale for the ratio quantities. Unadjusted HR and CI are expected exact; the corrected estimate and the bounds exact only if the seed and draw layout match the document's — say which, from §1.4.
- D1 / D3 vs the DINA column (grade 3 & `pgr <= 20`, N 113, 61 events) and G1 / G3 vs the GRF column (`er <= 0`, N 82, 45 events): region, N, events, quantities. A different region is recorded with the identifier's own seed and forest facts, not a stop.
- Budget arithmetic in the report: "full" = S3 + D3 + G3 elapsed (S3 includes the search); "reduced" = S1 + D1 + G1; the search alone = S0. Beside today's reference points: both current vignettes 45 s together; `devtools::check()` 11 m 53 s with tests 403 s.

---

## 4. Report — `dev/reports/REPORT_vignettes_stage0_2026-10-27.md`

Bullets, one item each: provenance; Inventory A (the directory tables; which document produced which column; the settings table; the transplant file list); Inventory B (§2.1 quoted; the formals file; the MR element paths; `declaration_c0` reachability; the build facts); the ladder table (rung × elapsed × max RSS × region × N/events × the nine quantities) with the comparisons of §3; the budget arithmetic; script deviations verbatim; OPEN ITEMS — anything that would need an `R/` change goes here as a finding, not fixed and not proposed as a task.

---

## 5. Close-out

1. Commit by explicit path: `dev/reports/vignettes_stage0_2026-10-27/` (scripts, logs, `ladder.csv`, `provenance.txt`, `package_facts.md`, `formals_forestsearch.txt`, `table1_calls/`; `.rds` files only if not gitignored and under 1 MB each — list what was and was not committed) and the report. Commit message `dev: report — vignettes stage 0 (2026-10-27)`. **No push.**
2. Copy to `~/Downloads/vignettes_stage0_2026-10-27/` (create it): the report; `ladder.csv`; `package_facts.md`; `formals_forestsearch.txt`; the three transplanted call files; `DESCRIPTION`, `.Rbuildignore`, `_pkgdown.yml`; `vignettes/forestsearch.Rmd`; `vignettes/survival-post-selection.qmd`; `vignettes/articles/methodology.*` and `vignettes/articles/paper_simulations.*`; and copies of the fs-glms-interpretable `.qmd` source(s) found in §1.2 (originals untouched). `ls -la` that directory into the report's last section. Larry attaches these to the chat.

---

## 6. Out of scope

- No edit to `R/`, `vignettes/`, `man/`, `DESCRIPTION`, `NEWS.md`, `_pkgdown.yml`, `.Rbuildignore`, `tests/`.
- No edit, build or render in fs-glms-interpretable; no re-render or re-run of any committed application document; no change to any committed number or payload.
- No `R CMD check`, vignette build or test suite; no parallel plan; no campaign; no push.
