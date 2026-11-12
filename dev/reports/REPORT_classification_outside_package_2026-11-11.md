# REPORT — classification metrics outside `R/`, and the `paper_simulations` re-render, 2026-11-11

Task: `dev/tasks/TASK_classification_outside_package_2026-11-11_v2.md` (committed as received, `760a9bfc`), with **Amendment 1** to Part B appended verbatim (`7f69573f`; quoted in §6).
Outcome: **every gate passed, in every part.** Neither branch is pushed or merged.

Definition (León et al. 2024, Section 3): with no subgroup identified, the estimated complement is everyone. Under an alternative the replicate has sens = 0, ppv = 0, spec = 1, npv = |Hc|/N; under the null spec = 1, npv = 1, sens and ppv NA.

## 1. Provenance

| Item | forestsearch | fs-glms-interpretable |
|---|---|---|
| Host, R | `pop-os`, R 4.6.2 | same |
| Branch | `fix/npv-nofind-instances` (existed; checked out) | `main` at start; `fix/npv-nofind-instances` created (did not exist: **gate passed**) |
| `git status --porcelain` at start | empty (**gate passed**) | empty (**gate passed**) |
| Start SHA | **`d738e54d`** (HEAD before Step 0.4; after the empty-set task) | **`396144c`** |
| `DESCRIPTION` | 0.4.1.9006 | — |
| End HEAD | given in the closing message (the report commit cannot hold its own SHA) | **`dc94a17`** |
| fs-post-selection | not touched: `main` `5d165ac`, clean, before and after | |

**Step 0.3: the empty-set task was run.**
- The v2 task arrived at 16:02; `TASK_empty_set_nofind_2026-11-11.md` was not in `dev/tasks/` or `~/Downloads`, and nowhere on the machine. The session waited for it in the background, doing read-only preparation only. It appeared in `~/Downloads` at 17:24, and Larry's chat message confirmed it.
- It ran completely, with its own gates: `1f1af85b` (task), `d04e4cd8` (fix), `d738e54d` (report `dev/reports/REPORT_empty_set_nofind_2026-11-11.md`, copied to `~/Downloads`). **All its gates passed**: suite 613 blocks / 8,270 expectations / 0 failures (645 s); `check()` 0/0/0 (994 s).
- Step 0.4: `760a9bfc dev: add classification outside-package task v2` (`cmp`-identical to `~/Downloads`). Amendment 1, received in chat, appended verbatim in `7f69573f`.

**Commits on the forestsearch branch, from the start SHA:**

| SHA | Subject | Content |
|---|---|---|
| `760a9bfc` | dev: add classification outside-package task v2 | the task, as received |
| `7f69573f` | dev: classification outside-package task v2, Amendment 1 (Part B tolerance gate) | Amendment 1 appended |
| `bd2ed66f` | fix: classification no-find values in simulation document sources | Part A, 113 files |
| `7733a87d` | fix: stored unconditional NPV in MR payloads re-summarised | Part B: 62 payloads, script, CSV |
| (report) | dev: classification outside-package report and paper_simulations re-render | this report, Part R outputs |

## 2. Status

| Part | Result |
|---|---|
| Step 0 | **passed** (empty-set task run with all its gates; task committed; start SHA recorded) |
| R — re-render `paper_simulations` | **passed** (library gate; render finished in 84 s, cap 1 h; classification tables present) |
| A — document sources | **passed** (113 files, all in rows Q1, Q2, Q3, Q5, Q8, Q9) |
| B — stored summary in payloads | **passed** (62 written, 44 unchanged, 0 BLOCKED; changed set = CSV) |
| C — what other renders need | **passed** (read-only table, §7) |
| D — fs-glms-interpretable | **passed** (only `npv_uncond` and its MCSE changed in the CSVs; 3 checks run, 3 pass; nothing under `manuscript/`) |
| E — report and close | gates in §9 |

## 3. Summary for Larry

| Instance | Done / BLOCKED / not an instance / needs a render or re-run | What changed | What is still needed |
|---|---|---|---|
| **Q1** template `oc_agg_from_results()` NPV fill | **Done** (A1, `bd2ed66f`) | a no-find row stored with NPV missing now takes the row's `propHc_true`; 1 only when absent; spec fill stays 1 | nothing in source. Latent with the current package (rows already carry values). The tracked `.html` is Q6 |
| **Q2** template prose "npv = 1 for Ĥᶜ" | **Done** (A3) | one sentence: complement is everyone, sens(Ĥᶜ) = 1, ppv(Ĥᶜ) = \|Qᶜ\|/N | the tracked `actg175_binary_methods_template.html` still shows the old prose until re-rendered (Q6) |
| **Q3** `reproduce_table1_m1.qmd` ×2 NPV fill | **Done** (A2) | as Q1 | nothing (no tracked render; latent) |
| **Q4** batch templates storing `fs_mr_oc_summary()` | **Done** (package `71804ead`; payloads re-summarised in Part B, `7733a87d`) | 62 stored unconditional NPVs corrected in place (Δ −0.0736 to −0.0002) | nothing |
| **Q5** fs1-fs6 summary over pre-A-003 bundles | **Done** (A5) | after loading, no-find H1 rows re-scored spec = 1, npv = `propHc_true`; two-line comment | needs a render to show it; the render is not tracked and the 54 bundles are machine-local (`/media/larryleon/bigstorage/…`) |
| **Q6** 17 tracked actg175 renders | **needs a render or re-run** | none (no render allowed) | per render, §7 |
| **Q7** `paper_simulations` text dumps (2026-10-30) | **Done by Part R** | new record `dev/reports/paper_sims_rerun_2026-11-11/`; the 2026-10-30 records are kept as history (not edited, by decision) | the pkgdown site is not rebuilt |
| **Q8** prose "matches `fs.estimates.out()`" ×105 | **Done** (A3) | one sentence in each document's notation | their renders (16 `.html`, 1 `.rmarkdown`) keep the old prose until re-rendered |
| **Q9** four documents: NA then `na.rm` | **Done** (A4) | definition's values at each no-find return | 2 tracked renders re-simulate on every render (§7) |
| **Q10** per-replicate NA, no summary | not edited (by decision) | — | — |
| **Q11** reports quoting detected-only NPV | not edited (by decision) | — | — |
| **Q12** correct / not applicable | not an instance | — | — |
| **G1** `identifiers_metrics.R` | **Done** (D1/D2, `dc94a17`) | non-declaring npv = 1 − `n_true`/`n_sample`; CSV re-run: 47 of 63 `npv_uncond` change | — |
| **G2** `identifiers_summary.R` | **Done** | as G1; 224 of 322 change | — |
| **G3** `idsweep_metrics.R` | **Done** | as G1; 150 of 288 change | — |
| **G4** supplement Table `tbl-supp-classification` (display) | **needs a render — Larry's decision** (nothing under `manuscript/` edited) | the CSV it reads is corrected; 13 of 36 displayed Npv cells would change (§8.5) | decide, then re-render the supplement |
| **G5** checks pinning npv = 1 | **Done** (D3) | reconcile: npv removed from the identity (cannot be verified from its CSVs); `verify_oc_section.R:772` strings updated; `:488-489`, `:1015-1018` need no change; `conv.R`: `u_npv` removed | — |
| **G6** text stating the old convention | **Done for the D4 files**; history not edited | README ×2, writeup ×5 passages, `_0Y` draft ×2 | `drafted_sections.pdf` needs a re-render; manuscript wording (§8.5) is Larry's decision |

## 4. Part R — `paper_simulations` re-rendered with the fixed package

### 4.1 Precedent (2026-10-30)

- **Record:** `dev/reports/articles_finish_2026-10-30/results/paper_simulations.txt` and `_output.txt`, pandoc text of the **pkgdown** build `docs/articles/paper_simulations.html`, made by `pkgdown::init_site()` then `pkgdown::build_article("articles/paper_simulations", lazy = FALSE)` from the repo root (`dev/tasks/cc_task_articles_finish_2026-10-30_v2.md:85`, `logs/render_pkgdown_paper_simulations.log:1-3`), then `strip_and_gate.sh`. It is the later of the two 2026-10-30 records (the `articles_alignment` pair came from an earlier `rmarkdown::render()` of the `.Rmd`).
- **Source:** the `.qmd`, md5 `6d64628f…`, byte-identical to HEAD (`git diff 7afae118 HEAD -- vignettes/articles/paper_simulations.qmd` is empty).
- **Settings and seeds** (in the source): `nsims_alt = nsims_null = 44L` (`:32-33`); `seed_base = 8316951` (`:837, :845`; `simulate_from_dgm(seed = seed_base + sim_id)`); doFuture `seed = TRUE`; `plan(multisession, workers = n_workers)`, `n_workers = min(detectCores() - 1, 120)` = 120 on this host then and now; inner `plan = "sequential"`.
- **Versions then:** R 4.6.1, forestsearch **0.4.1** (R code of `c2967362`), pkgdown 2.2.1, quarto (R) 1.5.1, Quarto CLI 1.10.18, pandoc 3.1.3.
- **Recorded wall clocks:** pkgdown build 68.6 s (the record's producer); local Quarto render 53.3 s and 54.6 s; alignment `rmarkdown::render` 47 s. **Cap:** max(3 × 68.6 s, 1 h) = **1 h**.

### 4.2 Library gate, render, wall clock

- **Library:** `git archive HEAD` (at `7f69573f`; package code = `d738e54d`) to `/tmp/fs_src_2026-11-11/`; `R CMD INSTALL --no-test-load --library=/tmp/fs_lib_2026-11-11 /tmp/fs_src_2026-11-11` (15.7 s).
- **Gate** (`logs/library_gate.txt`), with `R_LIBS=/tmp/fs_lib_2026-11-11`: `packageVersion("forestsearch")` = **0.4.1.9006** = `DESCRIPTION`; `find.package()` = `/tmp/fs_lib_2026-11-11/forestsearch`; the same in a `plan(multisession, workers = 2)` worker (its own PID). **PASS.**
- **Scratch copy:** `/tmp/paper_sims_rerun_2026-11-11/pkg` (the archived package tree; the article md5 `6d64628f…` there, in the repo and in the 2026-10-30 render). Nothing was written under the repository's `vignettes/` or `docs/`.
- **Render** (the precedent's route): `cd /tmp/paper_sims_rerun_2026-11-11/pkg; R_LIBS=/tmp/fs_lib_2026-11-11 timeout 3600 Rscript -e 'pkgdown::init_site(); pkgdown::build_article("articles/paper_simulations", lazy = FALSE, quiet = FALSE)'` (`quiet = FALSE` only logs Quarto's output; the precedent ran `init_site()` as its own call). `R_LIBS` puts the temporary library first for the R session, Quarto's knitr process and the multisession workers.
- **Wall clock: 84 s measured** (18:23:17 to 18:24:41; `build_article` itself 70.8 s), against 68.6 s for the precedent's `build_article`. Exit 0; no error or warning in the log.
- **Gate:** the output holds the classification tables: the `format_oc_results()` tables under H1 and H0, the Type I table, and the `build_classification_table()` table. **PASS.** The page reports forestsearch 0.4.1.9006 (`sessionInfo()`), "Using 120 parallel workers", "Results: 88 rows".
- **Outputs** (`dev/reports/paper_sims_rerun_2026-11-11/`): `paper_simulations.html` and its two figures (`paper_simulations_files/figure-html/`, added with `git add -f` by explicit path because `.gitignore:4` ignores `*_files/`); `results/paper_simulations.txt`, `results/paper_simulations_output.txt`, `results/keyvalues.txt` (the precedent's paper_simulations gates 5, 6, 12a, 12b: all PASS) made by `strip.sh` (the precedent's method); `logs/` (render log, timing, library gate, install log, provenance); `accounting/` (scripts, `accounting_grouped.csv` with every number, `accounting_tables.md` with one table per block, `group_lines.txt`).

### 4.3 Accounting

**Settings and seeds are identical between the two renders.** Same source (md5 above), same `nsims`, `seed_base`, doFuture seeding and 120 workers, same pkgdown/Quarto/pandoc route and versions. What differs: R 4.6.1 → 4.6.2 and forestsearch 0.4.1 → 0.4.1.9006 (the other `sessionInfo()` packages have the same versions; the only `sessionInfo()` lines whose text changed are R's base packages and forestsearch).

**Method.** Every number visible on the page (prose, console output, tables; not code, not section numbers) of the 2026-10-30 build (`docs/articles/paper_simulations.html`, git-ignored, mtime 2026-10-30 14:55:38; its pandoc text is byte-identical to the committed record `paper_simulations.txt`, checked with `cmp`) is paired with the same number of the 2026-11-11 build. Lines are aligned by GNU diff with numbers masked, and table cells by their column labels. A number is a token not preceded by a letter, digit, `.` or `_`, so version strings (`R version 4.6.1`, `forestsearch_0.4.1`, `grid_4.6.1`) count at most their first two parts: the R and package version changes are reported in prose, not in the counts. The three `format_oc_results()` tables appear as `[TABLE]` in both text records, so their values come from the HTML tables. 849 lines on each side, all paired; 1,544 numbers.

- **Groups:** (i) `sens(Hc)`, `ppv(Hc)`, Spec, NPV; (ii) detection (Found / any(H) / Power / Type I error / empirical detection / "x of 44" and "x/44 estimable" counts) and PPV / ppv(H); (iii) everything else.
- "run" marks a run-dependent value: timing, the render date, `gc()` memory, the session date.
- Blocks are the article's sections. Below are the blocks with group (i) or (ii) numbers computed by the render, or with any change. The full per-block listing (every block, including Appendix B's 664 typed-in numbers, 255 of them in groups (i)/(ii), all unchanged) is `accounting/accounting_tables.md`; every number is in `accounting/accounting_grouped.csv`.

**Title block** (26 numbers; 2 changed)

| Block | Metric and cell | 2026-10-30 | 2026-11-11 | Difference |
|---|---|---|---|---|
| (iii) run | render date, month (`date: today`) | 10 | 11 | +1 |
| (iii) run | render date, day | 30 | 11 | −19 |

**5.3 Running Alternative Hypothesis Simulations** (94 numbers; 2 changed) and **5.4 Running Null Hypothesis Simulations** (91; 2 changed)

| Block | Metric and cell | 2026-10-30 | 2026-11-11 | Difference |
|---|---|---|---|---|
| 5.3 (iii) run | "Seconds and minutes forestsearch overall", seconds (verbose sim 1) | 2.5 | 2.159 | −0.341 |
| 5.3 (iii) run | same, minutes | 0.0417 | 0.036 | −0.0057 |
| 5.4 (iii) run | same, seconds | 1.632 | 1.684 | +0.052 |
| 5.4 (iii) run | same, minutes | 0.0272 | 0.0281 | +0.0009 |

**6.1 Operating Characteristics Summary Under HTEs** (72 numbers; 4 changed; Table 6 = `format_oc_results()`, H1)

| Group | Metric and cell | 2026-10-30 | 2026-11-11 | Difference |
|---|---|---|---|---|
| (i) | FS Spec | 0.98 | 0.99 | +0.01 |
| (i) | FS NPV | 0.98 | 0.97 | −0.01 |
| (i) | GRF Spec | 0.95 | 0.99 | +0.04 |
| (i) | GRF NPV | 1.00 | 0.95 | −0.05 |
| (ii) | FS Found / PPV | 0.89 / 0.80 | 0.89 / 0.80 | 0 / 0 |
| (ii) | GRF Found / PPV | 0.61 / 0.58 | 0.61 / 0.58 | 0 / 0 |
| (ii) | "FS: 39/44 (89%) estimable"; "Ĥ: 39 estimable" | 39, 89, 39 | 39, 89, 39 | 0 |

**6.1.1 Interpretation (H1)** (22 numbers; 0 changed): (ii) "39 of 44 simulations (88.6%)": unchanged.

**6.2 Operating Characteristics Under NULL** (67 numbers; 0 changed; Table 8)

| Group | Metric and cell | 2026-10-30 | 2026-11-11 | Difference |
|---|---|---|---|---|
| (i) | FS Spec / NPV | 0.99 / 1.00 | 0.99 / 1.00 | 0 / 0 |
| (i) | GRF Spec / NPV | 1.00 / 1.00 | 1.00 / 1.00 | 0 / 0 |
| (ii) | FS Found / PPV; GRF Found / PPV | 0.05 / 0.00; 0.02 / 0.00 | same | 0 |
| (ii) | "FS: 2/44 (5%) estimable"; "Ĥ: 2 estimable" | 2, 5, 2 | 2, 5, 2 | 0 |

**6.2.1 Interpretation (H0)** (20 numbers; 0 changed): (ii) "2 of 44 simulations (4.5%)", "2 false detections": unchanged.

**6.2.2 Classification Rate Details** (28 numbers; 4 changed; Table 10 = `build_classification_table()`)

| Group | Metric and cell | 2026-10-30 | 2026-11-11 | Difference |
|---|---|---|---|---|
| (i) | M Null sens(Hc) FS / GRF | 0.99 / 1.00 | 0.99 / 1.00 | 0 / 0 |
| (i) | M Null ppv(Hc) FS / GRF | 1.00 / 1.00 | 1.00 / 1.00 | 0 / 0 |
| (i) | **M Alt sens(Hc) FS** | 0.98 | 0.99 | +0.01 |
| (i) | **M Alt sens(Hc) GRF** | 0.95 | 0.99 | +0.04 |
| (i) | **M Alt ppv(Hc) FS** | 0.98 | 0.97 | −0.01 |
| (i) | **M Alt ppv(Hc) GRF** | 1.00 | 0.95 | −0.05 |
| (ii) | any(H) FS / GRF, null | 0.05 / 0.02 | same | 0 |
| (ii) | any(H) FS / GRF, alt | 0.89 / 0.61 | same | 0 |
| (ii) | ppv(H) FS / GRF, alt | 0.80 / 0.58 | same | 0 |

**7 Theoretical Subgroup Detection Rate Approximation** (52 numbers; 0 changed; includes Table 11, Type I)

| Group | Metric and cell | 2026-10-30 | 2026-11-11 | Difference |
|---|---|---|---|---|
| (i) | Type I table FS Spec / NPV; GRF Spec / NPV | 0.995 / 1.000; 0.998 / 1.000 | same | 0 |
| (ii) | Type I table FS Found / PPV; GRF Found / PPV | 0.045 / 0.000; 0.023 / 0.000 | same | 0 |
| (ii) | Empirical FS, GRF (H1) | 0.886, 0.614 | same | 0 |
| (ii) | Empirical FS, GRF (H0) | 0.045455, 0.022727 | same | 0 |

**7.1 Key Metrics** (10 numbers; 2 changed)

| Group | Metric and cell | 2026-10-30 | 2026-11-11 | Difference |
|---|---|---|---|---|
| (i) | "FS: … Spec =" | 0.961 | 0.991 | +0.030 |
| (i) | "GRF: … Spec =" | 0.961 | 0.991 | +0.030 |
| (ii) | FS / GRF Power; PPV | 0.750; 0.691 | same | 0 |
| (ii) | FS / GRF Type I Error | 0.0341 | same | 0 |

(The two lines are identical because the chunk's loop variable `analysis` is also a column name, so each line pools FS and GRF: an existing article defect, OPEN ITEM 6.)

**11 Computational Timing Summary** (35 numbers; 15 changed, all run): Time (sec) DGM H1 0.1 → 0.2, Cox calibration 3.9 → 3.8, AHR 1.4 → 1.5, validate 0.4 → 0.5, DGM H0 0.0 → 0.1, sims H1 5.9 → 5.8, sims H0 3.4 → 3.5, total 42.4 → 42.5; and 7 "% of Total" cells (DGM H1, Cox, AHR, validate, sims H1, sims H0, summarize H0; ±0.1 to ±0.4).

**17 Cleanup** (`gc()`; 12 numbers; 8 changed, all run): Ncells used 2,795,680 → 2,796,967 and max 4,652,291 → 4,652,385 (×2); Vcells used 5,848,236 → 5,850,028, gc trigger 12,303,779 → 12,285,035 (93.9 → 93.8 Mb), max used 10,186,483 → 10,170,863 (77.8 → 77.6 Mb).

**18 Session Info** (54 numbers; 2 changed, run): the R date 2026-06-24 → 2026-10-29 (R 4.6.1 → 4.6.2).

**Every other block: no number changed** (`accounting_tables.md`).

**Group lines:**

- **(i)** 166 numbers (22 computed by the render; 144 are the 2024 published values typed into Appendix B), **10 changed**, largest change **−0.05** (GRF NPV in Table 6 and GRF ppv(Hc) in Table 10, 1.00 → 0.95).
- **(ii)** 151 numbers (39 computed; 112 typed in Appendices), **0 changed**.
- **(iii)** 1,227 numbers (783 outside the appendix sections 15-16, 444 in them), **31 changed, all run-dependent** (render date 2, forestsearch-overall timing 4, timing table 15, `gc()` 8, session date 2). Largest change −18,744 (`gc()` Vcells gc trigger).

**Commits for each changed (iii) number** (between forestsearch 0.4.1, R code `c2967362`, and the rendered `d738e54d`), touching the functions the block calls. The prep's static call graph covers `R/` at `c2967362` and HEAD; listed, not judged:

| Block (changed numbers) | Functions the block calls | Commits touching them |
|---|---|---|
| Title block (render date) | none from forestsearch | none |
| 5.3, 5.4 ("forestsearch overall" seconds/minutes) | `run_simulation_analysis()` → `simulate_from_dgm()`, `forestsearch()`, `.extract_fs_estimates_gen()`, `.extract_grf_estimates_gen()` | `8c4ef5f6` (`simulate_from_dgm`, `.entry_delay_vector`), `808e284d` (`forestsearch`), `e2f9ec18` (both extractors), `d04e4cd8` (both extractors) |
| 11 Timing table | none from forestsearch (gt over recorded `Sys.time()`) | none |
| 17 Cleanup (`gc()`) | none | none |
| 18 Session Info | none | none |

For reference: group (i) changes come from blocks calling `format_oc_results()`, `build_classification_table()` and the key-metrics means over `run_simulation_analysis()` rows; commits touching those functions are `e2f9ec18`, `71804ead` and `d04e4cd8`.

## 5. Part A — document sources (`bd2ed66f`)

113 files, all `.qmd` under `quarto/`, `git diff --stat`: 589 insertions, 641 deletions. **Gate:** the changed set equals the union of the rows' file lists (Q1/Q2: 1, Q3: 2, Q5: 1, Q8: 105, Q9: 4 → 113), compared with `diff`. **PASS.** No `.html`, `.rmarkdown` or `.R` file and nothing under `dev/reports/` was touched.

| Row | Files | Change |
|---|---|---|
| A1 (Q1) | `quarto/simulations/actg175/binary_methods/actg175_binary_methods_template.qmd` (`oc_agg_from_results()`, `:356-376`; `hc` `:358-360`, fill `:362-365`) | `hc <- propHc_true` (or NA if absent); `npv_u <- ifelse(is.na(npv) & any.H == 0, ifelse(is.na(hc), 1, hc), npv)`; spec fill unchanged |
| A2 (Q3) | `quarto/simulations/gbsg/reproduce_table1_m1.qmd`, `quarto/simulations/reproduce_table1_m1.qmd` (`:641-643`; the two files were byte-identical and remain so) | the same fill with `ra$propHc_true` |
| A3 (Q2) | the template (`:119-128`) | 1 file; before/after below |
| A3 (Q8) | 105 files (13 `simulations/_legacy/`, 3 `actg175/_archive_pre_fsparams/`, 72 `actg175/binary_legacy/`, 17 `actg175/binary_methods/`) | one byte-identical 7-line passage → one 6-line sentence, exactly once per file; before/after below |
| A4 (Q9) | `quarto/dina/grf_dmin_calibration.qmd` (`:162-170`) | no-find return: sens 0 / ppv 0 / spec 1 / npv `mean(!truth)` (= \|Hc\|/N); under the null (`truth` all FALSE) sens, ppv NA and npv 1 |
| A4 (Q9) | `quarto/applications/validation_glm_simulation_study.qmd` (`:586-594`) | new `else if (!found)` branch: alternative (`result$true_H` not NA) sens 0, ppv 0, spec 1, npv `mean(flag_harm == 0)`; null spec 1, npv 1, sens/ppv NA. A found replicate without `df.est` keeps NA |
| A4 (Q9) | `quarto/archive/actg175_simulations_binary_working.qmd` (helper `:491-501`; returns `:503`, `:513`) | `fill_nofind()` at both no-find returns (`fs_res` NULL; no `sg.harm.id`): null (`dgm$model_type`) spec 1, npv 1; alternative sens 0, ppv 0, spec 1, npv `mean(flag_harm != 1)` |
| A4 (Q9) | `quarto/dina/dina_subgroup_signature_simulation.qmd` (`:237-248`) | alternative no-find: sensitivity 0, specificity 1, ppv 0 (NA if the true mask is empty). The document has no NPV. The null arm returns before scoring for every replicate (`:216`, "detection only"), and no null-arm classification is displayed, so it is unchanged |
| A5 (Q5) | `quarto/simulations/actg175/binary_legacy/_summaries/actg175_binary_m1_harm_fs1-fs6_summary.qmd` (end of `load-bundles`, `:173-182`) | re-score `results$alt` rows with `any.H == 0` and `propHc_true` present: spec = 1, npv = `propHc_true`; two-line comment |

**Before/after, cluster Q2** (template, 1 file):

> Before: "… is averaged **unconditionally** (a replication with no detected subgroup contributes sens $= 0$, ppv $= 0$ for $\hat H$ and spec $= 1$, npv $= 1$ for $\hat H^c$)."
>
> After: "… is averaged **unconditionally**.  When no subgroup is detected, the estimated complement $\hat H^c$ is everyone, so the replication contributes $\text{sens}(\hat H) = 0$ and $\text{ppv}(\hat H) = 0$ for the subgroup, and $\text{sens}(\hat H^c) = 1$ and $\text{ppv}(\hat H^c) = |Q^c|/N$ for the complement (León et al. 2024, Section 3)."

**Before/after, cluster Q8** (105 files, e.g. `quarto/simulations/actg175/binary_methods/actg175_binary_m1_harm_methods.qmd:1301-1308`):

> Before: "Classification metrics are averaged across **all** simulations (León et al. 2024, Table 1).  When no subgroup is detected ($\hat{Q} = \emptyset$), the simulation contributes $\text{sens}(\hat{Q}) = 0$ and $\text{ppv}(\hat{Q}) = 0$ (since no member of $Q$ is captured and no subjects are assigned to $\hat{Q}$).  This convention is implemented at the source in `run_simulation_analysis()`, matching the original `fs.estimates.out()` function used in León et al. (2024)."
>
> After: "Classification metrics are averaged across **all** simulations (León et al. 2024, Table 1).  When no subgroup is detected ($\hat{Q} = \emptyset$), the estimated complement $\hat{Q}^c$ is everyone, so the simulation contributes $\text{sens}(\hat{Q}) = 0$ and $\text{ppv}(\hat{Q}) = 0$ for the subgroup, and $\text{sens}(\hat{Q}^c) = 1$ and $\text{ppv}(\hat{Q}^c) = |Q^c|/N$ for the complement (León et al. 2024, Section 3)."

- Verification: `fs.estimates.out` now occurs in 0 tracked `quarto/**/*.qmd` (105 before); `npv $= 1$` in the template 0 (1 before); the new Q8 sentence in exactly 105 files.
- The 8 code-edited documents were `knitr::purl`ed and `parse()`d: all parse. The A1, A2, A4 and A5 code is the code the preparation tested on hand-built rows and on one real A5 bundle (FS Spec 0.9237 → 0.9341, NPV 0.9414 → 0.9309).
- **Not edited, by decision:** Q7, Q10, Q11, Q12.

## 6. Part B — stored summary in the payloads (`7733a87d`)

**Amendment 1 (verbatim, received in chat during the run):**

> TASK_empty_set_nofind_2026-11-11.md is now in ~/Downloads; proceed with Step 0.3. Amendment to Part B, authorized by Larry, record it verbatim in the report: (1) the per-file gate compares the recomputed $oc with the stored one using a numeric tolerance of 1e-8, not exact equality; a file is BLOCKED only if some value other than the unconditional NPV differs by more than that. (2) Write only the unconditional NPV cell into the stored payload; keep every other stored value bit for bit, and confirm after writing that the object is identical() to the original except that one cell. (3) Do not write a file whose unconditional NPV changes by 1e-8 or less.

- **Script:** `dev/reports/REPORT_classification_outside_package_2026-11-11_resummarise.R`.
  - `pkgload::load_all()` of this tree (0.4.1.9006), sequential.
  - It enumerates every tracked `.rds` (2,712) and keeps those with an `fs_mr_oc` object at `$oc`: **106**.
  - It recomputes `fs_mr_oc_summary(payload)` with all defaults, as the five batch templates call it (`sim_fs_mr_field_or_template.qmd:1220-1222, :1316-1317`; `sim_fs_maxeffCons_mr_field_md_template.qmd:1135, :1221`; `sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_1000.qmd:875, :968`, `…_md40_knoise0_n700_batch_1_1000.qmd:875, :968`, `…_md120_knoise0_n500_batch_1_1000.qmd:893, :986`). The function reads only `$results`, `$truth` and `$meta`.
  - Gate: Amendment 1 (1). Write: Amendment 1 (2)-(3), with the original's serialization version and compression.
- **Run:** 14.8 s (cap 15 min). Outcome: **62 written, 44 unchanged (\|change\| ≤ 1e-8), 0 BLOCKED.** The largest difference in any other value was 2.5e-14 (floating-point last bits on payloads built on another machine). A dry re-run after the writes finds all 106 unchanged.
- **Counts:** 40 in `quarto/simulations/actg175/binary_020/mr_or_harm/`, 6 in `…/binary_020/mr_or_harm_superseded_prev09632/`, 16 in `…/continuous/mr_md_harm/` — **62, as expected**, matching `REPORT_npv_nofind_instances_2026-11-11.md` §2.4. Change new − old: min −0.073628 (`fs_maxeffCons_mr_md120_knoise0_n500_c1star_res_1_1000.rds`, 0.7919 → 0.7183), mean −0.008990, max −0.0001535. The 44 unchanged are 34 payloads with no non-detected replicate and 10 null cells (|Hc|/N = 1).
- **Gate** (checked by a different route from the script: each written file against its `git show HEAD:` version): the changed `.rds` set equals the 62 files of `REPORT_classification_outside_package_2026-11-11_payloads.csv`. In each, restoring the old unconditional NPV gives an object `identical()` to HEAD's. Each keeps gzip header `1f8b0800000000000003`, XDR, serialization version 3. The CSV values are within 4.4e-16 of the stored doubles (15 significant digits). **PASS.**
- **Note:** a rewritten file's serialization header now records the writer R version 4.6.2. Values are kept bit for bit. Of the 62, 20 were first written with R 4.5.2 on another machine (32 of the 106 were); in 13 of those 20 a few internal flag bits of the stream are not written back by R 4.6.2 (143-180 decompressed bytes differ, against 7-8 in the others); `identical()` is unaffected.
- **Without Amendment 1:** under the original exact-equality gate the preparation's dry run would have BLOCKED 20 of the 62. Those 20 were built on that machine (R 4.5.2), whose recomputed `estimation` and `identification` values differ from this machine's in the last bits (≤ 2.5e-14).

## 7. Part C — what each other affected render needs (read-only; no render)

19 tracked renders: the 17 of row Q6 (16 in `binary_methods/`, including `actg175_binary_methods_template.html`, which is therefore not an extra row, and 1 in `binary_legacy/`), plus the tracked renders of 2 of the 4 Q9 documents. `validation_glm_simulation_study.qmd` and `archive/actg175_simulations_binary_working.qmd` have no render, tracked or on disk.

Common facts:
- (a) Every Q6 bundle is **untracked**, outside git on `/media/larryleon/bigstorage/forestsearch_sims/quarto/simulations/actg175/…` (abbreviated B), and present on this machine.
- (b) As committed, every Q6 source runs `sim_mode = "full"` or `"demo"`, so a plain render **re-simulates**, and a full render **overwrites its bundle**. `readRDS(BUNDLE_PATH)` runs only in `"import"` mode.
- (c) Every bundle stores the **pre-A-003 values** on H1 no-find rows (npv = 1, spec = `propHc_true`), with package 0.2.0 recorded. An import-mode render repeats the old printed numbers exactly unless the A5-style re-score is added after loading.
- (d) The installed user library is forestsearch **0.4.1.9002, pre-A-003**. A re-simulating render with the default library re-creates the old convention; only a library at `e2f9ec18` or later (0.4.1.9004+) stores |Hc|/N.
- (e) Recorded "total runtime" values are the sum of the timed components (bootstrap + LOO + simulation loops) printed in each render's Timing Summary, not full wall clocks.

| Render | Source `.qmd` | What it reads (path; tracked or not; present on this machine) | Does a render run simulations (cells × replicates) or only read stored results? | Do the stored rows hold pre-fix values, so a plain re-render would repeat the old numbers? | Recorded wall clock of a past render (source), or none |
|---|---|---|---|---|---|
| `binary_methods/actg175_binary_m1_harm_effMinSG-both_fs3_recommended.html` | same name (embedded source differs only in worker count `:336, :559`) | B`/binary/actg175_binary_m1_harm_fs3_effMinSG_both_n700_sims1000.rds`; untracked; present | runs: n 700, 1000 H1 + 1000 H0 × 2 analyses (FS, GRF), plus single-trial bootstrap and LOO; overwrites the bundle | yes: 97/97 H1 no-find rows pre-fix (FS s(Hc)/p(Hc) 0.9312/0.9395 → 0.9400/0.9306) | 55.60 min timed (html `:13389-13392`); bundle 777.2 s + 592.3 s |
| `binary_legacy/actg175_binary_m1_harm_effMinSG-both_fs3_recommended.html` (byte-identical to the row above) | `binary_legacy/…_fs3_recommended.qmd` (identical) | same bundle | same | same | same |
| `binary_methods/actg175_binary_m1_harm_eff_fs3.html` | same name | B`/binary/actg175_binary_m1_harm_fs3_eff_neighborhood_n700_sims1000.rds`; untracked; present | runs: n 700, 1000 + 1000 × 2 analyses; overwrites | yes: 95/95 | 31.94 min (`:13079-13082`) |
| `binary_methods/actg175_binary_m1_harm_methods.html` | same name, but rendered from an older version (sims 10, thresholds and `run_single_trial` since changed) | nothing stored (demo mode does not save) | runs: current source 100 + 100 × 10 arms at n 700, plus single trial | no stored rows; printed values from the pre-A-003 producer; a re-render gives different numbers (settings changed) | 3.53 min (`:9982-9985`) |
| `binary_methods/actg175_binary_m1_harm_methods_n1000_full.html` | same name | B`/binary_methods/actg175_binary_m1_harm_fs3_methods_both_n1000_sims1000.rds`; untracked; present | runs: n 1000, 1000 + 1000 × 10 arms, single trial; overwrites | yes: 1466/1466 (largest dina-eff ±0.0273) | 90.05 min (`:13008-13011`) |
| `binary_methods/actg175_binary_m1_harm_methods_n700_demo.html` | **no tracked `.qmd` of that name** (embedded source = `…_n700.qmd` with demo mode) | nothing behind this render; `…_n700.qmd` now points to an n700 full bundle (pre-fix, 959/959) | runs (demo 100 + 100 × 10 arms in the render; current `n700.qmd` is full mode) | no stored rows behind the html | 38.42 min (`:15139-15142`) |
| `binary_methods/actg175_binary_m1b_harm_methods_n1000_full.html` | same name; thresholds changed since the render | B`/binary_methods/actg175_binary_m1b_harm_fs3_methods_both_n1000_sims1000.rds`; present | runs: 1000 + 1000 × 10 arms; overwrites; new thresholds | yes: 2984/2984 (largest ±0.0599) | 75.24 min (`:12864-12867`) |
| `binary_methods/actg175_binary_m1c_harm_methods_n1000_full.html` | same name | B`/binary_methods/actg175_binary_m1c_harm_fs3b_methods_both_n1000_sims1000.rds`; present | runs: 1000 + 1000 × 10 arms; overwrites | yes: 1514/1514 | 57.31 min (`:12816-12819`) |
| `binary_methods/actg175_binary_m1c_harm_n1000_full.html` | same name | B`/binary_methods/actg175_binary_m1c_harm_fs3c_methods_both_n1000_sims1000.rds`; present | runs: 1000 + 1000 × 10 arms; overwrites | yes: 1554/1554 | 56.33 min (`:12821-12824`) |
| `binary_methods/actg175_binary_m1c_harm_n700_full.html` | same name | B`/binary_methods/actg175_binary_m1c_harm_fs3c_methods_both_n700_sims1000.rds`; present | runs: n 700, 1000 + 1000 × 10 arms; overwrites | yes: 1058/1058 | 56.56 min (`:15426-15429`) |
| `binary_methods/actg175_binary_m1d_harm_n1000_full.html` | same name | B`/binary_methods/actg175_binary_m1d_harm_fs3d_methods_both_n1000_sims1000.rds`; present | runs; overwrites | yes: 1611/1611 | 45.83 min (`:12812-12815`) |
| `binary_methods/actg175_binary_m1d_harm_n1000_full_adjusted.html` | same name | B`/binary_methods/actg175_binary_m1d_harm_adjusted_fs3d_methods_both_n1000_sims1000.rds`; present | runs; overwrites | yes: 1607/1607 | 47.17 min (`:12882-12885`) |
| `binary_methods/actg175_binary_m1e_harm_n1000_full_adjusted.html` | same name | B`/binary_methods/actg175_binary_m1e_harm_adjusted_fs3e_methods_both_n1000_sims1000.rds`; present; sources tracked `fdr_*.R` | runs: 1000 + 1000 × 10 arms + FDR calibration (20 null trials); overwrites | yes: 1495/1495 | 39.98 min (`:10244-10247`; FDR untimed) |
| `binary_methods/actg175_binary_m1e_harm_n1000_full_adjusted_split.html` | same name | B`/binary_methods/actg175_binary_m1e_harm_adjusted_split_fs3e_methods_both_n1000_sims1000.rds`; present | runs; overwrites | yes: 1673/1673 | 37.58 min (`:12858-12861`) |
| `binary_methods/actg175_binary_m1e_harm_n700_full_adjusted.html` | same name | B`/binary_methods/actg175_binary_m1e_harm_adjusted_fs3e_methods_both_n700_sims1000.rds`; present | runs (n 700) + FDR calibration; overwrites | yes: 1001/1001 | 34.52 min (`:10244-10247`) |
| `binary_methods/actg175_binary_m1f_harm_n1000_full_adjusted_resample.html` | same name | B`/binary_methods/actg175_binary_m1f_harm_adjusted_resample_fs3e_methods_both_n1000_sims1000.rds`; present | runs; overwrites | yes: 2682/2682 (fs-eff 0.9031/0.9733 → 0.9594/0.9170) | 33.43 min (`:10092-10095`) |
| `binary_methods/actg175_binary_methods_template.html` | same name (now also carries A1 and A3) | nothing stored (demo; `_bundles/` absent) | runs: 10 + 10 × 2 arms at n 700, plus calibration and single trial | no stored rows; printed with forestsearch 0.2.0 (fs-eff any(H) 0.50, ppv(Hc) 0.96: 5 of 10 are no-finds) | sims H1 107.17 s, H0 94.87 s (`:5135-5152`) |
| `quarto/dina/grf_dmin_calibration.html` | `quarto/dina/grf_dmin_calibration.qmd` (A4 edited) | nothing stored; simulates inline | runs on every render: 108 cells × 500 = 54,000 GRF fits | no stored rows; the old render's spec/npv were detected-only (document scorer). A re-render of the current source gives the definition's values | TOTAL 172.97 min (html `:3752-3777`) |
| `quarto/dina/dina_subgroup_signature_simulation.html` | `quarto/dina/dina_subgroup_signature_simulation.qmd` (A4 edited) | nothing stored; simulates inline | runs on every render: 12 cells, 3,600 replicates | no stored rows; as above (sens/spec now 0/1 on non-null no-finds) | none (`message()` suppressed; no log) |

- **Outside the task's list, also tracked:** `dev/temp_share/actg175_binary_m1_harm_methods_n700_full.html` (`3c953499`, 2026-06-04) is the full-mode render of the pre-fix n700 bundle (`B/binary_methods/actg175_binary_m1_harm_fs3_methods_both_n700_sims1000.rds`, 959/959 H1 no-find rows pre-fix): its H1/H0 loop times 31.44 and 26.69 min equal the bundle's `t_alt_sec` 1886.6 s and `t_null_sec` 1601.6 s, and it shows 10 "Classification Under H1" tables (total runtime 89.48 min, `:15545`). It is not a Q6 render (not under `quarto/`); it needs the same re-render as the bundle-backed rows.

## 8. Part D — fs-glms-interpretable (`dc94a17` on `fix/npv-nofind-instances`)

### 8.1 D1 (G1, G2, G3)

- `dev/payloads/R/identifiers_metrics.R:360-368`: `fill <- c(sens = 0, spec = 1, ppv = 0)`. A non-declaring replicate's npv is `1 - r$n_true / m$n_sample`. `emit()` with `mean(u)`, `mcse(u)`, and an `na_note` when `n_true` is missing (so an NA does not trip the script's NA stop).
- `dev/payloads/R/identifiers_summary.R:148-156`: the same through `u`, which also feeds `npv_uncond_mcse`. `na_notes` gains the NA case (`:213`). Header comment `:13-15` corrected.
- `dev/payloads/R/idsweep_metrics.R:149-155`: the same with `n <- m$n_sample`. `emit_mean()`'s `na_note` covers the NA case (`:117-118`). Header comment `:31-32` corrected.
- The sens, spec and ppv fills are unchanged.
- **No cell is NA:** `n_true` is never missing on a non-declaring replicate (13,306 in the 81 MR-grid bundles; 10,235 in the 288 idsweep bundles).

### 8.2 D2 re-run and gate

- Run the repository's way (`Usage` lines; README `:85-87, :186, :211`), base R only, from the repo root:
  - `Rscript dev/payloads/R/identifiers_metrics.R`: 49 s;
  - `Rscript dev/payloads/R/idsweep_metrics.R`: 17 s;
  - `Rscript dev/payloads/R/identifiers_summary.R <dir-of-46>`: 3 s. The 46 bundles were linked in scratch from the CSV's own `source_file` column: 37 in `bundles/`, 9 in `bundles/superseded/`, as README `:201-211` prescribes. A first attempt ran with an empty directory because of an argument-passing slip; the script stopped at its own count check and wrote nothing.
- Total 69 s (cap 30 min).
- The scripts do not call forestsearch. The installed version is 0.4.1.9002, recorded only.
- `COLUMNS.md` and `COLUMNS_idsweep.md` were regenerated byte-identical.
- **Gate** (old = `git show HEAD:`; compared by key, with row order unchanged):

| CSV | Rows × cols | Cells changed | Columns changed | Result |
|---|---|---|---|---|
| `identifiers_metrics.csv` | 46,998 × 23 | 94 in 47 rows, all `metric == npv_uncond` | `value`, `mcse` | **PASS** |
| `idsweep_metrics.csv` | 8,640 × 25 | 300 in 150 rows, all `npv_uncond` | `value`, `mcse` | **PASS** |
| `identifiers_summary.csv` | 322 × 104 | 448 in 224 rows | `npv_uncond`, `npv_uncond_mcse` | **PASS** |

### 8.3 D3 (G5) checks

| Check | Change | Run (cap 10 min) | Result |
|---|---|---|---|
| `dev/writeup/identifiers/R/reconcile_identifiers_summary.R:91-96` | npv removed from the identity `unconditional = conditional × declared + fill × (1 − declared)`. It cannot verify the new value from the CSVs it reads: the non-declaring replicates' \|Hc\|/N is in none of them. One-line comment; header `:13` notes "for sens, spec, ppv" | 0.2 s | **PASS** (every residual ≤ 1e-12) |
| `dev/writeup/identifiers/R/verify_oc_section.R:772` | the two literal strings it greps for in `idsweep_metrics.R` now match the new fill and the `npv` branch | 1.1 s | **PASS** ("holds 40; reproduces 131", GATE PASS). It reads `.qmd`/`.R` sources, not renders |
| `verify_oc_section.R:488-489`, `:1015-1018` | none needed: they read `npv_uncond` from the CSV and do not pin 1; same verdicts on the new values | (in the run above) | PASS |
| `dev/writeup/identifiers/R/conv.R` | console calculation over inline numbers. It cannot produce the new value, so `u_npv` is removed (`:38` → one-line comment); prints the conditional npv only; comments `:4-12` corrected | 0.2 s | **PASS** (runs) |

- All three were run again after the D4 text edits, with the same results.
- **Not in the G5 list, not run here:** `reconcile_complete_grid.R` and `reconcile_identifiers_metrics.R`. Both already fail at the start SHA for unrelated reasons, found in the read-only preparation: the 2026-09-24 s5rerun swap, and duplicate keys (README `:203`).

### 8.4 D4 (G6) text corrected

| File | Passage |
|---|---|
| `dev/payloads/identifiers/README.md` | `:248-251` (convention bullet) and `:276` (column table: "non-declaring scored 0/1/0/1") |
| `dev/writeup/identifiers/identifier_evaluation_writeup.qmd` | `:467` (table row npv = 1), `:470-474` (rationale "no assignment … matches León et al. (2024, Table 1)"), `:479-480` ("posts spec = 1 and npv = 1"), `:862-864` and `:871-874` ("no new extraction" now for sens/spec/ppv; NPV's unconditional form with the mean \|Hc\|/N of the non-detected replicates) |
| `dev/writeup/identifiers/draft/_0Y_recovery_and_location.qmd` | `:9-13` (convention paragraph), `:62-66` ("perfect benefit purity") |

- Not edited: decision records, requests, reports and status updates (history), e.g. `REQUEST_fs_classification_metrics_2026-09-12.md:131, 161, 218`, `dev/manuscript_status_update.md:965-967`, `REPORT_oc_comprehensive_2026-09-13.md:322`. See also OPEN ITEMS 8-9.

### 8.5 D5 — manuscript impact, for Larry's decision (nothing under `manuscript/` edited or rendered)

**Supplement Table `tbl-supp-classification`, Npv column**, via `manuscript/_R/oc_figures.R` `oc_tbl_classification()` on the old and new `identifiers_metrics.csv`. Displayed "value (MCSE)"; every other column of the table is identical.

| CSV line | Cell | Id | Old | New | Difference |
|---|---|---|---|---|---|
| 9704 | 12.4%, HR 1.50, n=500 | FS | 0.947 (0.001) | 0.937 (0.001) | −0.010535 |
| 10450 | 12.4%, HR 1.50, n=500 | DINA | 0.938 (0.001) | 0.923 (0.001) | −0.015177 |
| 11196 | 12.4%, HR 1.50, n=500 | GRF | 0.926 (0.001) | 0.925 (0.001) | −0.000202 |
| 12688 | 12.4%, HR 1.50, n=1000 | FS | 0.967 (0.001) | 0.964 (0.001) | −0.0030955 |
| 13434 | 12.4%, HR 1.50, n=1000 | DINA | 0.956 (0.001) | 0.937 (0.001) | −0.0197485 |
| 14180 | 12.4%, HR 1.50, n=1000 | GRF | 0.947 (0.001) | 0.946 (0.001) | −0.000302 |
| 15672 | 12.4%, HR 1.50, n=1500 | FS | 0.982 (0.001) | 0.981 (0.001) | −0.001432 |
| 16418 | 12.4%, HR 1.50, n=1500 | DINA | 0.968 (0.001) | 0.943 (0.001) | −0.024680 |
| 17164 | 12.4%, HR 1.50, n=1500 | GRF | 0.962 (0.001) | 0.962 (0.001) | −0.000288 |
| 18656 | 12.4%, HR 1.75, n=500 | FS | 0.954 (0.001) | 0.949 (0.001) | −0.005761 |
| 19402 | 12.4%, HR 1.75, n=500 | DINA | 0.945 (0.001) | 0.933 (0.001) | −0.011037 |
| 20148 | 12.4%, HR 1.75, n=500 | GRF | 0.935 (0.001) | 0.935 (0.001) | −0.000066 |
| 21640 | 12.4%, HR 1.75, n=1000 | FS | 0.976 (0.001) | 0.976 (0.001) | −0.000569 |
| 22386 | 12.4%, HR 1.75, n=1000 | DINA | 0.963 (0.001) | 0.952 (0.001) | −0.011275 |
| 23132 | 12.4%, HR 1.75, n=1000 | GRF | 0.955 (0.001) | 0.955 (0.001) | −0.000120 |
| 24624 | 12.4%, HR 1.75, n=1500 | FS | 0.990 (0.001) | 0.989 (0.001) | −0.000119 |
| 25370 | 12.4%, HR 1.75, n=1500 | DINA | 0.974 (0.001) | 0.962 (0.001) | −0.012418 |
| 26116 | 12.4%, HR 1.75, n=1500 | GRF | 0.967 (0.001) | 0.967 (0.001) | 0 |
| 33576 | 31%, HR 1.50, n=500 | FS | 0.837 (0.002) | 0.837 (0.002) | −0.000165 |
| 34322 | 31%, HR 1.50, n=500 | DINA | 0.808 (0.002) | 0.808 (0.002) | −0.000471 |
| 35068 | 31%, HR 1.50, n=500 | GRF | 0.816 (0.002) | 0.816 (0.002) | 0 |
| 35814 | 31%, HR 1.50, n=1000 | FS | 0.887 (0.002) | 0.887 (0.002) | 0 |
| 36560 | 31%, HR 1.50, n=1000 | DINA | 0.863 (0.002) | 0.863 (0.002) | 0 |
| 37306 | 31%, HR 1.50, n=1000 | GRF | 0.878 (0.002) | 0.878 (0.002) | 0 |
| 38052 | 31%, HR 1.50, n=1500 | FS | 0.937 (0.002) | 0.937 (0.002) | 0 |
| 38798 | 31%, HR 1.50, n=1500 | DINA | 0.918 (0.002) | 0.918 (0.002) | 0 |
| 39544 | 31%, HR 1.50, n=1500 | GRF | 0.928 (0.002) | 0.928 (0.002) | 0 |
| 40290 | 31%, HR 1.75, n=500 | FS | 0.866 (0.002) | 0.866 (0.002) | −0.000165 |
| 41036 | 31%, HR 1.75, n=500 | DINA | 0.824 (0.002) | 0.824 (0.002) | 0 |
| 41782 | 31%, HR 1.75, n=500 | GRF | 0.840 (0.002) | 0.840 (0.002) | 0 |
| 42528 | 31%, HR 1.75, n=1000 | FS | 0.911 (0.002) | 0.911 (0.002) | 0 |
| 43274 | 31%, HR 1.75, n=1000 | DINA | 0.882 (0.002) | 0.882 (0.002) | 0 |
| 44020 | 31%, HR 1.75, n=1000 | GRF | 0.899 (0.002) | 0.899 (0.002) | 0 |
| 44766 | 31%, HR 1.75, n=1500 | FS | 0.951 (0.002) | 0.951 (0.002) | 0 |
| 45512 | 31%, HR 1.75, n=1500 | DINA | 0.931 (0.002) | 0.931 (0.002) | 0 |
| 46258 | 31%, HR 1.75, n=1500 | GRF | 0.940 (0.002) | 0.940 (0.002) | 0 |

- **20 of 36 values change** (17 at 12.4%, 3 at 31%); **13 displayed values change**. Mean −0.00327, min −0.02468 (DINA 12.4% HR 1.50 n=1500). No displayed MCSE changes.
- Where declaration is below 1, the unconditional Npv now falls **below** Npv (c), e.g. DINA 12.4% HR 1.50 n=1500: 0.943 against 0.959. Before, it was above (0.968). `_S4:180-183` states no direction.

**The other two CSVs**, change in `npv_uncond` (new − old, all rows):

| CSV | Rows | Changed | Min | Mean | Max |
|---|---|---|---|---|---|
| `identifiers_summary.csv` | 322 | 224 | −0.081255 | −0.011943 | 0 |
| `idsweep_metrics.csv` | 288 | 150 | −0.082191 | −0.010133 | 0 |
| (and `identifiers_metrics.csv`, all 63 rows) | 63 | 47 | −0.081255 | −0.011520 | 0 |

The largest change everywhere is DINA at 12.4%, HR 1.00, n=1500, e.g. 0.970597 → 0.889343 in `identifiers_metrics.csv`.

**Every file:line under `manuscript/` that displays or states the old convention:**

| File:line | What |
|---|---|
| `manuscript/_R/oc_figures.R:377-400` (`oc_tbl_classification`; `:390-395` builds `Npv` from `npv_uncond`), fed by `oc_load` `:15-61` | displays (code) |
| `manuscript/_S4_operating_characteristics_supplement.qmd:16`, `:174-178` (chunk `tbl-supp-classification`), included at `manuscript/fs_glms_interpretable_supplement.qmd:98` | displays (source) |
| `manuscript/fs_glms_interpretable_supplement.tex:3005-3163` (longtable; caption `:3018-3022`; headers `Npv` `:3046`, `Npv (c)` `:3048`; rows `:3054-3161`) | displays (rendered) |
| PDFs, pp. 67-68: `manuscript/fs_glms_interpretable_supplement.pdf`, `manuscript/fs_glms_interpretable_supplement_blinded.pdf`, `manuscript/submissions/jasa_2026-10-13/{blinded,unblinded}_supplement.pdf` (the latter two are byte copies of the former two) | display |
| `manuscript/_05_operating_characteristics.qmd:25`; `manuscript/fs_glms_interpretable.tex:1125-1127`; main PDFs p. 25 (`fs_glms_interpretable{,_blinded}.pdf`, `submissions/jasa_2026-10-13/{blinded,unblinded}_main.pdf`) | wording: "a replicate that declares nothing being scored as making no assignment" — the old convention's rationale; no value stated |
| `manuscript/_S4_operating_characteristics_supplement.qmd:168-169`; `fs_glms_interpretable_supplement.tex:2995-2997`; supplement PDFs p. 66 | the same wording, directly above the table |
| `manuscript/fs_glms_interpretable_supplement.html:4023, :4057-4058, :4163…`; `manuscript/fs_glms_interpretable.html:505` | git-ignored, stale renders on disk (old values / wording) |

Consistent with the definition, so not listed as instances: `_02_identification.qmd:76`, `_S4…:60` (sens 0 and spec 1 for no region); `_S4…:41` and `supplement.tex:2763-2766` (a declaring replicate whose Ĥ is the whole trial); `_08_guohe_supplement.qmd:203, :221, :543` (NPV 1 by nesting).

## 9. Part E gates

| Gate | Result |
|---|---|
| `git diff d738e54d --stat -- R/ tests/ man/ vignettes/` | **empty**. The empty-set task's changes are before the start SHA: `git diff 96feee06 --stat -- R/ tests/ man/ vignettes/` lists only `R/forestsearch_cross_validation.R`, `R/run_simulation_analysis.R`, `tests/testthat/test-cv-mixed-and-fatal.R`, `tests/testthat/test-nofind-classification.R`, i.e. what the empty-set task changed. **PASS** |
| `git status --porcelain` empty in forestsearch and fs-glms-interpretable | checked after the report commit; given in the closing message |
| fs-post-selection unchanged | `main` `5d165ac`, clean, not touched. **PASS** |

## OPEN ITEMS

1. **G4, the supplement table (Larry's decision).** `identifiers_metrics.csv` is corrected; Table `tbl-supp-classification` (13 displayed Npv cells) and the "scored as making no assignment" wording (`_05:25`, `_S4:168-169`, and their `.tex`/PDF copies, including the 2026-10-13 JASA submission PDFs) are unchanged. The same table and wording are in the arXiv companion snapshot (`dev/arxiv/companion_v1/{build,package}/*.tex`, `companion_arxiv_v1.pdf`, `.tar.gz`) and the reference copies under `dev/reference/submission_items_preparation/`.
2. **Q6, the 17 actg175 renders** (§7), plus the tracked `dev/temp_share/actg175_binary_m1_harm_methods_n700_full.html` behind the pre-fix n700 bundle. The bundle-backed renders show the old values; they change only by (a) re-simulating with a library at 0.4.1.9004 or later (the user library is 0.4.1.9002), which overwrites each bundle (31.9-90.1 min timed each), or (b) switching to `"import"` and adding the A5 re-score after loading. The three demo renders (`…_methods.html`, `…_n700_demo.html`, the template) have nothing stored and must re-simulate; `…_n700_demo.html` has no tracked source of its name.
3. **Q9 tracked renders** re-simulate on every render (`grf_dmin_calibration`: 172.97 min recorded). Pre-existing points the A4 edits do not touch: in `validation_glm_simulation_study.qmd` and `archive/…_working.qmd`, a **found** null replicate is still scored against `flag_harm` (a would-be subgroup with no effect) rather than an empty H; `dina_subgroup_signature_simulation.qmd` keeps Jaccard NA on a no-find and its caption `:443-448` ("not-found reps contribute 0 via NA handling only where applicable") no longer describes spec; the template's `sens_u`/`ppv_u` fill 0 under H0 (no hypothesis check; not displayed).
4. **Q5** needs a render to show its correction; neither the render nor the 54 bundles are tracked. The comment at the summary's `:93-94` names a stale path.
5. **Part R record.** The pkgdown site is not rebuilt; the local, git-ignored `docs/articles/paper_simulations.html` still holds the 2026-10-30 build. Scratch left in place: `/tmp/fs_src_2026-11-11/`, `/tmp/fs_lib_2026-11-11/`, `/tmp/paper_sims_rerun_2026-11-11/`.
6. **`paper_simulations` key-metrics chunk** (`vignettes/articles/paper_simulations.qmd:1340-1355`): the loop variable `analysis` shadows the column, so the FS and GRF lines print the same pooled values (7.1 above). Article not edited.
7. **Part B, writer metadata.** The 62 rewritten payloads now carry R 4.6.2 in their serialization header (values bit for bit, `identical()` checked).
8. **fs-glms text not in the D4 list, left unchanged:** `dev/writeup/identifiers/oc_long_form.qmd:109-111, :290` ("scored as making no assignment"); `_0Y_recovery_and_location.qmd:74-75` ("specificity and benefit purity … should not be quoted alone": its reason for benefit purity is weaker now); the writeup's "*undefined, not missing*" at `:861-862`, which reports what the 2026-09-12 extraction states; the `_0Y` tables `:17-24`, `:41-45` hold conditional rates (sens, spec, PPV, NPV) under unconditional prose (pre-existing); `dev/writeup/identifiers/draft/drafted_sections.pdf` needs a re-render to show the D4 text.
9. **fs-glms generated column docs.** `idsweep_metrics.R:80` and `identifiers_metrics.R:487` describe `n_true` as feeding sizes and (idsweep) the sens/spec denominators; it now also feeds `npv_uncond`. Not changed, so `COLUMNS.md`/`COLUMNS_idsweep.md` stay byte-identical.
10. **Empty-set task open items** (its report): the GRF guard reading (fixed on the literal reading; BLOCKED on the other), the GRF frontier coding defect at `R/grf_main.R:226`, runner columns not reset, pre-existing stops (a non-empty FS set with NA members at `R/run_simulation_analysis.R:1133`; CV fold assignments {1, NA} or {0, NA} at `R/forestsearch_cross_validation.R:1363`/`:1412`), CV downstream medians, NEWS heading.
