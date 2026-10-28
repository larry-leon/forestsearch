# CC task — post-selection vignette stage 0: Table S2 provenance (fs-post-selection), package facts, MR and bootstrap timing (2026-10-28)

**Repo:** `larry-leon/forestsearch` — run CC from this repo's root, on branch `feature/vignettes-0.4.1`.
**Kind:** read-only inventory (fs-post-selection and fs-glms-interpretable read-only; forestsearch at HEAD) plus one small measurement. **No edit to `R/`, `vignettes/`, `man/`, `DESCRIPTION`, `NEWS.md`, `_pkgdown.yml`, `tests/`; nothing edited, built or rendered in any sibling repo.** Writes: this task doc under `dev/tasks/`; the record `dev/reports/postsel_stage0_2026-10-28/`; the report `dev/reports/REPORT_postsel_stage0_2026-10-28.md`.
**Compute:** the ladder of §3 — four MR fits (the three identifiers at the manuscript's settings, plus one FS variant), and 5 bootstrap replicates per identifier, all single-process `plan = "sequential"`. Unmeasured for these settings; stage 0 of 2026-10-27 measured 37 / 8 / 25 s for the three Table 1 fits. Hard cap 20 min per run, 2 h for the ladder. Approved by Larry 2026-10-28 (the post-selection vignette is to reproduce Table S2 of León & Anderson 2026a with the bootstrap illustrated at a small B).
**Unattended.** Gates stop the task — a failed gate is a stop even when the reason is clear; CC reports with a proposed fix and chat decides. No `git fetch` / `pull` / `push`; explicit `git add` paths; never stage pre-existing untracked files; no `R CMD check`, vignette build, test suite or `load_all()`.

**Table S2 of 2026a (Supplementary S6, p. 46), the comparison target.** De-biased subgroup hazard ratios for the GBSG trial, 95% confidence intervals in parentheses; Ĥ is `er <= 0 & size <= 35` for FS, `er <= 0` for GRF, `grade3 >= 1 & pgr <= 10` for DINA.

| Method | SG | n (%) | Naive | Full bootstrap | Multiplier |
|---|---|---|---|---|---|
| FS | Ĥ | 61 (9) | 2.54 (1.25, 5.17) | 1.94 (0.89, 4.21) | 1.44 (0.62, 3.36) |
| FS | Ĥᶜ | 625 (91) | 0.61 (0.47, 0.79) | 0.63 (0.43, 0.93) | 0.64 (0.40, 1.03) |
| DINA | Ĥ | 89 (13) | 1.60 (0.88, 2.91) | 0.89 (0.38, 2.07) | 1.11 (0.45, 2.75) |
| DINA | Ĥᶜ | 597 (87) | 0.63 (0.48, 0.82) | 0.66 (0.40, 1.11) | 0.66 (0.40, 1.09) |
| GRF | Ĥ | 82 (12) | 1.95 (1.04, 3.67) | 1.04 (0.49, 2.24) | 1.05 (0.47, 2.34) |
| GRF | Ĥᶜ | 604 (88) | 0.61 (0.47, 0.80) | 0.66 (0.39, 1.13) | 0.66 (0.41, 1.06) |

Settings per §6 of the paper: the same candidate family as Table 1 of 2026b (at most two cuts; ≥ 60 subjects with 10 events in each arm; `J = 10` cuts on `pgr` and `er`; `{er <= 0}` and `{pgr <= 0}` pre-specified), screening and consistency thresholds both log 1.0, p\* ≥ 0.9; the selection rule returns the largest-effect region whose size is maximal within a **10%** neighbourhood of the maximum; GRF at depth 2 with a DR-score floor of 0, DINA with a mean-effect floor of log 1.0; full bootstrap B = 1,000; multiplier resampling 5,000 draws; the multiplier column is the two-term estimate with the infinitesimal-jackknife two-sided interval.

---

## 0. First action — GATE

1. `git status --short`: record pre-existing untracked paths; never stage them. *GATE:* branch `feature/vignettes-0.4.1`; HEAD contains `99e8a64a` (`git merge-base --is-ancestor`). Record HEAD.
2. Copy this file to `dev/tasks/cc_task_postsel_stage0_2026-10-28.md`; `git add` that path only; commit **`dev: task — post-selection vignette stage 0 (2026-10-28)`**.
3. Provenance to `dev/reports/postsel_stage0_2026-10-28/provenance.txt` as at the last stage 0. *GATE:* installed forestsearch `0.4.0.9000` and `git diff 936b6734 HEAD -- R/ DESCRIPTION NAMESPACE inst/` is empty (the installed build is HEAD's `R/` tree; no NEWS identity check).

---

## 1. Inventory A — where Table S2 came from (read-only)

1. **The manuscript repo.** Find `fs-post-selection` under `~/Documents/GitHub/` or `~/GitHub/` (record which; *GATE:* found and a git repo; record HEAD and `git status --short`). Locate the source of Table S2: grep its `.qmd` / `.Rmd` / `.tex` / `.md` and any payload manifest for `De-biased subgroup hazard ratios`, `supp_gbsg`, and the published values `2.54`, `1.94`, `1.44`, `0.89, 4.21`. Record the fragment that renders the table, the payload or data file it reads, and the manifest entry (path, commit, date) that names where that payload was produced. Copy the fragment and the manifest excerpt into `dev/reports/postsel_stage0_2026-10-28/manuscript_source/`.
2. **The producing document.** Follow the manifest to the document and payload that computed the table — expected in forestsearch `quarto/applications/gbsg/` (`analysis_gbsg_survival_multimethod.qmd` and its `_payloads/`), possibly a dated re-run directory beside it; if the manifest points into fs-glms-interpretable or fs-post-selection instead, follow it there (read-only). *GATE:* the document and a payload carrying Table S2's FS row (2.54 / 1.94 / 1.44 at 2 dp, N 61) are located; record render date, `forestsearch_version` and the commit. If the bootstrap column came from a separate run or document (the B = 1,000 bootstrap may have been run apart from the MR fits), record that too.
3. **Transplants** into `dev/reports/postsel_stage0_2026-10-28/table_s2_calls/`, verbatim with a header giving source path, chunk names and commit: the data preparation; the three identifier `forestsearch()` calls (FS, DINA, GRF) as the document ran them for the MR column; the bootstrap call(s) — `forestsearch_bootstrap_dofuture()` or whatever the document used — with their arguments (B, seeds, plan, `mr_in_replicates`); and the code that assembled the naive / bootstrap / multiplier table (`mr_estimates_table()`, `summarize_bootstrap_results()` or document-local code).
4. **Settings table**, one row per identifier: every argument that differs from the Table 1 calls of stage 0 2026-10-27 (`table1_calls/`), with its value — expected at least `effect_neighborhood` (0.10 vs 0.20), `mr_inference_args$ci_method`, and whatever the DINA/GRF floors were — plus the bootstrap arguments. Quote both documents' values side by side.
5. Wall-clock and machine for the document's MR fits and its B = 1,000 bootstrap, if its report or payload records them.

## 2. Inventory B — package facts from source

Record in `package_facts_s2.md`, quoting file and line.

1. `forestsearch_bootstrap_dofuture()`: formals with defaults; what it needs from the fit object; whether it accepts a DINA or GRF fit (`subgroup_method` handling in `bootstrap_results()`), what each replicate re-runs (search + selection; MR stripped under `mr_in_replicates = FALSE`), how the de-biased estimate is formed from the replicates, and what object it returns.
2. `mr_estimates_table()` and `summarize_bootstrap_results()`: formals; whether `mr_estimates_table(fs, boot_results, est.scale = "hr")` yields the six Table S2 cells per identifier (naive, bootstrap-corrected, MR two-term with IJ interval, for Ĥ and Ĥᶜ) and in what format (gt / data frame); the element paths of the bootstrap-corrected estimate and its interval in the bootstrap object.
3. The IJ two-sided interval's paths in the MR object (`g$debiased$lower`, `$upper`, and for the complement), and whether they are produced under `ci_method = "ij"` and under `"field"` alike.
4. The print of the bootstrap summary helpers as rendered by `knitr::kable` or gt in a quarto vignette (gt is in Imports); note anything that would not render in a self-contained HTML vignette.

---

## 3. Compute — the ladder (approved 2026-10-28)

Single process, `plan = "sequential"`, installed package, the transplanted calls verbatim plus the overrides `parallel_args = list(plan = "sequential")`, `quiet = TRUE`, `details = FALSE`, `plot.sg = FALSE` (and `plot.grf = FALSE` where the formal exists). Driver as at the last stage 0 (`timeout 1200 /usr/bin/time -v Rscript ladder.R <rung>`), per-rung logs, `ladder.csv`, per-rung `.rds` of the args and the fit.

| rung | what | records |
|---|---|---|
| M-FS | the FS call as the document ran it (MR 5,000 draws) | region, N and events in Ĥ / Ĥᶜ; naive est and CI; two-term est and IJ CI, both regions; the field bounds if `ci_method = "field"`; `n_family`; wall-clock, peak RSS |
| M-DINA, M-GRF | the DINA and GRF calls as the document ran them | same |
| M-FS-ij | M-FS with `mr_inference_args$ci_method = "ij"` (if the document used `"field"`; the reverse if it used `"ij"`) | same, plus `all.equal()` of `g$debiased` (est, lower, upper, both regions) against M-FS, and the wall-clock difference |
| B-FS, B-DINA, B-GRF | the document's bootstrap call on each M fit with `nb_boots = 5` (seed as the document's) | per-replicate wall-clock (total / 5, and the first replicate separately), the bootstrap-corrected estimates the 5 replicates give, whether the call runs at all for DINA / GRF |

**GATE A (at M-FS):** the region is `{er <= 0} & {size <= 35}` with N 61; the naive HR and CI equal Table S2's FS row at 2 dp (2.54; 1.25, 5.17); the multiplier estimate and interval equal 1.44 (0.62, 3.36) at 2 dp. Otherwise STOP and report with the full `sg.harm`, counters and the settings in force. Recorded, not gated: M-DINA against `{grade3 >= 1} & {pgr <= 10}` N 89 and 1.11 (0.45, 2.75); M-GRF against `{er <= 0}` N 82 and 1.05 (0.47, 2.34); every quantity also against the located payload at full precision.

**Budget arithmetic in the report:** the three M fits' wall-clock; the FS-ij saving; per-replicate bootstrap cost per identifier, and from it the B each identifier can afford in 30 s, 60 s and 120 s of single-core time. No recommendation — the vignette task sets B on these numbers.

---

## 4. Report — `dev/reports/REPORT_postsel_stage0_2026-10-28.md`

Bullets: provenance; §1 (the fragment, the manifest, the producing document and payload, the settings table, transplant list); §2 facts; the ladder table with the Table S2 comparison and the budget arithmetic; deviations verbatim; OPEN ITEMS (findings only). Commit the record directory and the report: **`dev: report — post-selection vignette stage 0 (2026-10-28)`**. No push.

Copy to `~/Downloads/postsel_stage0_2026-10-28/`: the report, `package_facts_s2.md`, `ladder.csv`, `table_s2_calls/`, `manuscript_source/`, the producing document's `.qmd`, and `vignettes/survival-post-selection.qmd` as committed. `ls -la` into the report's last section.

---

## 5. Out of scope

No edit anywhere outside `dev/tasks/` and `dev/reports/`; no render; no bootstrap beyond `nb_boots = 5` per identifier; no re-run of any committed application document; no push.
