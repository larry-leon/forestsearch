# REPORT — Inventory of `quarto/extreme_subgroups/` (read-only)

**Date:** 2026-10-19
**Task:** `dev/tasks/TASK_hta_move_and_es_inventory.md`, Part D
**Repo state when inventoried:** `master` after `6224bdfd` (Part C3)
**Method:** `find`, `git ls-files`, `git check-ignore`, `git log -1`, `grep`, `shasum -a 256`. Nothing in the folder was moved, renamed, deleted or rendered.

**Report location.** `REPORT_*` files are spread over several directories. The first `git ls-files | grep REPORT_ | head` hit is `dev/glm-continuous-sims/`, which holds reports for one workstream only. This report is in `dev/reports/` instead: that is the general reports directory (36 files) and also the task's fallback.

---

## D1. File table

There are 40 files. All are tracked: none is untracked or ignored. The folder has no `_files/`, `_freeze/` or `.quarto/` directories. `.gitignore:31` ignores `fixed_random/results/`, but that directory does not exist on disk.

| Path (relative to `quarto/extreme_subgroups/`) | Size | Status | Last commit | Category |
|---|---:|---|---|---|
| `build_deck.js` | 62.0 KB | tracked | 2026-08-12 — extreme subgroups slides | other (Node/pptxgenjs deck builder; writes `forestsearch_fixed_baseline_slides.pptx`) |
| `extreme_subgroups.html` | 5.82 MB | tracked | 2026-08-11 — extreme subgroups | rendered output |
| `extreme_subgroups.qmd` | 87.5 KB | tracked | 2026-08-11 — extreme subgroups | analysis .qmd |
| `extreme_subgroups_deck.pdf` | 252.5 KB | tracked | 2026-08-11 — create quarto extreme_subgroups | rendered output (PDF export of deck) |
| `extreme_subgroups_deck.pptx` | 430.1 KB | tracked | 2026-08-11 — create quarto extreme_subgroups | slides |
| `extreme_subgroups_design_comparison.html` | 2.12 MB | tracked | 2026-08-11 — fixed vs random comparison extreme subgroups | rendered output |
| `extreme_subgroups_design_comparison.qmd` | 20.4 KB | tracked | 2026-08-11 — fixed vs random comparison extreme subgroups | analysis .qmd |
| `extreme_subgroups_fixed.html` | 5.83 MB | tracked | 2026-09-02 — Housekeeping: parent-dir render refresh; archive legacy fixed/random versions | rendered output |
| `extreme_subgroups_fixed.qmd` | 92.4 KB | tracked | 2026-08-11 — Add extreme_subgroups_fixed.qmd: fixed-baseline companion vignette | analysis .qmd |
| `extreme_subgroups_section_6_6_4.pdf` | 628.2 KB | tracked | 2026-08-11 — create quarto extreme_subgroups | rendered output (PDF export of deck) |
| `extreme_subgroups_section_6_6_4.pptx` | 693.9 KB | tracked | 2026-08-11 — create quarto extreme_subgroups | slides |
| `extreme_subgroups_short.qmd` | 62.2 KB | tracked | 2026-08-11 — extreme subgroups | analysis .qmd |
| `extreme_subgroups_v1.qmd` | 61.0 KB | tracked | 2026-08-11 — create quarto extreme_subgroups | analysis .qmd |
| `extreme_subgroups_v2.qmd` | 60.0 KB | tracked | 2026-08-11 — create quarto extreme_subgroups | analysis .qmd |
| `FEASIBILITY_extreme_subgroups_HTA.docx` | 26.2 KB | tracked | 2026-09-15 — Create FEASIBILITY_extreme_subgroups_HTA.docx | rendered output (.docx; no source .qmd in folder) |
| `fixed_random.zip` | 6.89 MB | tracked | 2026-08-12 — extreme subgroup fixed vs random | other (archive; see D4) |
| `forestsearch_fixed_baseline_slides.pptx` | 757.3 KB | tracked | 2026-08-12 — extreme subgroups slides | slides |
| `fixed_random/actg175_extreme_subgroups_continuous_fixed.html` | 4.53 MB | tracked | 2026-09-03 — Arc C': retire kar_med completely (fifth site); Findings from the 10k | rendered output |
| `fixed_random/actg175_extreme_subgroups_continuous_fixed.qmd` | 12.8 KB | tracked | 2026-09-03 — Arc C': retire kar_med completely (fifth site); Findings from the 10k | analysis .qmd |
| `fixed_random/extreme_subgroups.html` | 5.72 MB | tracked | 2026-09-02 — Phase 2.1: lean formula environment in subgroup_cox(); render timing restored | rendered output |
| `fixed_random/extreme_subgroups.qmd` | 80.6 KB | tracked | 2026-09-02 — Phase 2: plot()/forest_height(); port extreme-subgroups vignettes to the wrappers | analysis .qmd |
| `fixed_random/extreme_subgroups_design_comparison.html` | 2.15 MB | tracked | 2026-09-02 — Phase 3: memo consolidation -- compare_subgroup_sims() replaces inline .col_stats() | rendered output |
| `fixed_random/extreme_subgroups_design_comparison.qmd` | 27.9 KB | tracked | 2026-09-02 — Phase 3: memo consolidation -- compare_subgroup_sims() replaces inline .col_stats() | analysis .qmd |
| `fixed_random/extreme_subgroups_fixed.html` | 8.39 MB | tracked | 2026-09-02 — Fixed vignette: analysis-vs-DGM note, payload documentation, workflow charts | rendered output |
| `fixed_random/extreme_subgroups_fixed.qmd` | 90.1 KB | tracked | 2026-09-02 — Fixed vignette: analysis-vs-DGM note, payload documentation, workflow charts | analysis .qmd |
| `fixed_random/extreme_subgroups_fixed_wrapper.html` | 5.73 MB | tracked | 2026-09-02 — update wrappers which is now  replaced | rendered output |
| `fixed_random/extreme_subgroups_fixed_wrapper.qmd` | 85.0 KB | tracked | 2026-09-02 — update wrappers which is now  replaced | analysis .qmd |
| `fixed_random/extreme_subgroups_slides.pdf` | 252.8 KB | tracked | 2026-09-02 — Relocate the Beamer deck into quarto/extreme_subgroups/fixed_random | slides |
| `fixed_random/extreme_subgroups_slides.qmd` | 35.1 KB | tracked | 2026-09-02 — Relocate the Beamer deck into quarto/extreme_subgroups/fixed_random | slides (beamer) |
| `fixed_random/extreme_subgroups_slides.tex` | 34.3 KB | tracked | 2026-09-02 — Relocate the Beamer deck into quarto/extreme_subgroups/fixed_random | slides (beamer `keep-tex` output) |
| `fixed_random/gbsg_extreme_subgroups_fixed.html` | 8.39 MB | tracked | 2026-09-03 — gbsg extreme subgroups rerun | rendered output |
| `fixed_random/gbsg_extreme_subgroups_fixed.qmd` | 90.1 KB | tracked | 2026-09-03 — gbsg extreme subgroups rerun | analysis .qmd |
| `fixed_random/legacy/extreme_subgroups.html` | 5.82 MB | tracked | 2026-09-02 — Housekeeping: parent-dir render refresh; archive legacy fixed/random versions | rendered output |
| `fixed_random/legacy/extreme_subgroups.qmd` | 87.5 KB | tracked | 2026-09-02 — Housekeeping: … archive legacy fixed/random versions | analysis .qmd |
| `fixed_random/legacy/extreme_subgroups_design_comparison.html` | 2.12 MB | tracked | 2026-09-02 — Housekeeping: … archive legacy fixed/random versions | rendered output |
| `fixed_random/legacy/extreme_subgroups_design_comparison.qmd` | 20.4 KB | tracked | 2026-09-02 — Housekeeping: … archive legacy fixed/random versions | analysis .qmd |
| `fixed_random/legacy/extreme_subgroups_fixed.html` | 5.83 MB | tracked | 2026-09-02 — Housekeeping: … archive legacy fixed/random versions | rendered output |
| `fixed_random/legacy/extreme_subgroups_fixed.qmd` | 92.4 KB | tracked | 2026-09-02 — Housekeeping: … archive legacy fixed/random versions | analysis .qmd |
| `fixed_random/sim_flow_part1.pdf` | 34.3 KB | tracked | 2026-09-02 — Relocate the Beamer deck into quarto/extreme_subgroups/fixed_random | rendered output (figure written by `extreme_subgroups_slides.qmd`) |
| `fixed_random/sim_flow_part2.pdf` | 30.7 KB | tracked | 2026-09-02 — Relocate the Beamer deck into quarto/extreme_subgroups/fixed_random | rendered output (figure written by `extreme_subgroups_slides.qmd`) |

Counts by category (40 files):

- 15 analysis .qmd
- 6 slides: 3 `.pptx`, plus the beamer `.qmd`, its `.tex` and its `.pdf`
- 0 R script
- 0 data
- 17 rendered output: 12 `.html`, 2 deck-export `.pdf`, 2 `sim_flow` `.pdf` and 1 `.docx`
- 0 freeze/cache
- 2 other: `build_deck.js` and `fixed_random.zip`

---

## D2. Per-.qmd detail

All GBSG documents use `data(gbsg, package = "survival")`, with `time_months = rfstime / 30.4375` and `treat = hormon`. None calls `source()`, `load()` or `read.*`. None reads from `dev/`. The only `dev/` mentions are comments that cite acceptance scripts: `dev/accept_phase1_bitident.R`, `dev/accept_phase3_memo_compare.R` and `dev/accept_slides_bitident.R`. `.param_or()` resolves Quarto `params`, with a default.

### Root level

| .qmd | `title` / `format` | Inputs | Outputs | Simulations |
|---|---|---|---|---|
| `extreme_subgroups.qmd` | "Extreme Subgroup Effects Under a Uniform Treatment Benefit" / `html` (self-contained) | `survival::gbsg` | none | Yes: an inline `foreach` loop over `simulate_from_dgm()`, `n_sims_uniform <- 10000L`. `baseline` is not passed, so the current default `"resample"` applies. `sim_config$n_per_trial = nrow(gbsg)` (super-population sampling). |
| `extreme_subgroups_fixed.qmd` | same title, subtitle "A Fixed-Baseline (Conditional-on-X) …" / `html` | `survival::gbsg` | none | Yes: inline loop, `n_sims_uniform <- 10000L`, `sim_config$baseline = "fixed"`. |
| `extreme_subgroups_short.qmd` | same title / `html` | `survival::gbsg` | none | Yes: inline `foreach`, `n_sims_null <- 10000L`. `baseline` is not passed (default `"resample"`). |
| `extreme_subgroups_v1.qmd` | same title / `html` | `survival::gbsg` | none | Yes: sequential `for` loop, `n_sims_null <- 200L`. `baseline` is not passed (default `"resample"`). |
| `extreme_subgroups_v2.qmd` | same title / `html` | `survival::gbsg` | none | Yes: `n_sims_null <- 200L`. `baseline` is not passed (default `"resample"`). |
| `extreme_subgroups_design_comparison.qmd` | "Fixed vs Random Baseline: A Design Comparison" / `html` | None. Per-subgroup results are a static data frame transcribed from the two rendered HTML vignettes, both at `n_sims_uniform = 10000`. | A `saveRDS(... "results/extreme_sims_%s_%d.rds")` appears only in an `eval: false` reproduction-notes chunk. | No. It compares `baseline = "resample"` with `"fixed"` descriptively. |

### `fixed_random/`

| .qmd | `title` / `format` | Inputs | Outputs | Simulations |
|---|---|---|---|---|
| `extreme_subgroups.qmd` | "Extreme Subgroup Effects Under a Uniform Treatment Benefit" / `html`; params `n_sims`, `save_payload`, `payload_dir`, `payload_file` | `survival::gbsg` | `saveRDS(sims_uniform, results/extreme_sims_resample_<n>_payload.rds)` when `save_payload` is TRUE (the default) | Yes: `run_subgroup_sims()`, `n_sims_uniform = .param_or("n_sims", 10000L)`, `baseline = "resample"` (`design_label = "resample"`). |
| `extreme_subgroups_fixed.qmd` | same title, fixed-baseline subtitle / `html`; same params | `survival::gbsg`; `readRDS(payload_path)` reads back the payload it just wrote | `saveRDS(sims_uniform, results/extreme_sims_fixed_<n>_payload.rds)` | Yes: `run_subgroup_sims()`, `n_sims` default 10000, `baseline = "fixed"`. |
| `gbsg_extreme_subgroups_fixed.qmd` | **byte-identical** to `extreme_subgroups_fixed.qmd` (D4) | same | same | same |
| `extreme_subgroups_fixed_wrapper.qmd` | same title, fixed-baseline subtitle / `html`; same params | `survival::gbsg` | `saveRDS(sims_uniform, results/extreme_sims_fixed_<n>_payload.rds)` | Yes: `run_subgroup_sims()`, `n_sims` default 10000, `baseline = "fixed"`. |
| `actg175_extreme_subgroups_continuous_fixed.qmd` | "Extreme Subgroup Effects Under a Uniform Treatment Effect: ACTG 175, Continuous Outcome" / `html` (`embed-resources: true`); params `n_sims`, `save_payload`, `payload_dir` | `speff2trial::ACTG175`, arms {1, 3}; outcome `cd40 - cd420` | `saveRDS(list(sim_hrs, sim_ubs, sim_ns, subgroups, effect, md_true, …), ./actg175_extreme_sims_continuous_fixed_<n>_payload.rds)` (`payload_dir` default `"."`) | Yes: `run_subgroup_sims()`, `n_sims` default 10000, `baseline = "fixed"` (N = 1083). |
| `extreme_subgroups_design_comparison.qmd` | "Fixed vs Random Baseline: A Design Comparison" / `html`; params `use_payloads`, `payload_dir`, `payload_random`, `payload_fixed` | `readRDS(results/extreme_sims_resample_10000_payload.rds)` and `readRDS(results/extreme_sims_fixed_10000_payload.rds)` (live mode, the default). Otherwise uses static numbers with `n_sims_r = n_sims_f = 10000L`. | none (an `str(readRDS(...))` display only) | No. It calls `compare_subgroup_sims()` on the two payloads. |
| `extreme_subgroups_slides.qmd` | same title / `beamer` (Madrid, 16:9, `keep-tex: true`) | `survival::gbsg`; `readRDS(results_rds)` when `results_rds` is non-NULL (default `NULL`). The `readRDS("results/extreme_sims_*.rds")` lines at L966/L1012–1013/L1041 are in `eval: false` chunks. | `sim_flow_part1.pdf/.png` and `sim_flow_part2.pdf/.png` (`rsvg_pdf`/`rsvg_png`); `saveRDS(sims_uniform, results_rds)` when set | Yes: `run_subgroup_sims()`, `n_sims_null <- 1000L` ("5000 for publication quality"), `baseline = "resample"`, `n = nrow(gbsg)`. |

### `fixed_random/legacy/`

These three are byte-identical to the root-level files of the same name (D4), so their details match the root rows: `extreme_subgroups.qmd` (10000, default resample), `extreme_subgroups_fixed.qmd` (10000, `"fixed"`) and `extreme_subgroups_design_comparison.qmd` (no sims).

---

## D3. External references

Command: `git grep -n "extreme_subgroups"`. This matches the folder path and also the unrelated article name `articles/extreme_subgroups`, so the results are split into two groups below.

### (a) Lines that point into `quarto/extreme_subgroups/` (`git grep -n "quarto/extreme_subgroups"`)

```
.gitignore:31:quarto/extreme_subgroups/fixed_random/results/
R/compare_subgroup_sims.R:4:# inline machinery (quarto/extreme_subgroups/fixed_random/
R/run_subgroup_sims.R:6:# vignettes (quarto/extreme_subgroups/fixed_random, commit 002f4f37):
dev/accept_phase2_payload_compare.R:21:#         default results_dir = quarto/extreme_subgroups/fixed_random/results
dev/accept_phase2_payload_compare.R:27:  "quarto/extreme_subgroups/fixed_random/results"
dev/accept_phase3_memo_compare.R:12:  "quarto/extreme_subgroups/fixed_random/results"
dev/diag_phase2_timing.R:26:res_dir <- "quarto/extreme_subgroups/fixed_random/results"
dev/glm-continuous-sims/REPORT_vi_grf_default_2026-08-30.md:75:  (directory count table mentioning quarto/extreme_subgroups/)
dev/tasks/TASK_extreme_subgroups_qmd_replace_2026-10-17.md:57:If `quarto/extreme_subgroups/extreme_subgroups_fixed.qmd` exists, record:
quarto/GuoHe/REPORT_guohe_record_payload_2026-09-11.md:37:Six pre-existing untracked files elsewhere (`quarto/extreme_subgroups/`, ...) left alone.
```

The two `R/` hits are comments only. No code reads the folder. The three `dev/*.R` scripts default to `fixed_random/results/`, which is absent on disk.

### (b) Checked locations from the task list

| Location | Result |
|---|---|
| `_pkgdown.yml` | L57 `href: articles/extreme_subgroups.html` (navbar) and L560 `- "articles/extreme_subgroups"` (articles index). Both name the **pkgdown article** `vignettes/articles/extreme_subgroups.qmd`, not this folder. |
| `vignettes/` | `vignettes/articles/extreme_subgroups.qmd` L226–227 (raw.githubusercontent URL of *itself*) and L242 (`quarto render extreme_subgroups.qmd`). Neither points into the folder. |
| `README*` | no match |
| `.Rbuildignore` | no `extreme_subgroups` match; `^quarto$` (L24) excludes the whole `quarto/` tree from the package build |
| `_quarto.yml` (`quarto/_quarto.yml`, the only one) | no match; project-level shared HTML format only, no render list |
| `.github/workflows/` (`pkgdown.yaml`) | no match |
| `docs/` | not tracked (`.gitignore:16 docs`) |
| `NEWS.md` | no match |

### (c) Other name-only matches (not paths into the folder)

These are listed for completeness. `dev/AI_review_updates/4-March_2026/{files,legacy}/_pkgdown.yml` hold old `_pkgdown.yml` copies with the same article entries. `dev/devtools_package.R:322,397` references `vignettes/articles/extreme_subgroups.Rmd` and `pkgdown::build_article("articles/extreme_subgroups")`. Five `dev/tasks/TASK_extreme_subgroups_*_2026-10-17.md` files concern the pkgdown article. `man/run_subgroup_sims.Rd:116` and `R/run_subgroup_sims.R:559` mention the "`extreme_subgroups` vignettes" in prose. `quarto/claude code notes/{CLAUDE.md,claude_code_guide.qmd,claude_md_how_to_use.qmd}` contain prose and examples. `quarto/archive/actg175_simulations_review.md:232` and `quarto/quarto_reorganize.sh:77` (`git mv extreme_subgroups.qmd guides/`) also match. `quarto/guides/build_deck.js:4,6,493,528` has a pkgdown URL at L493, and three binary decks in `quarto/guides/` also match.

---

## D4. Weight

- **Total size:** 77,023,303 bytes, about 73.5 MiB (`du -sh` gives 74M).
- **Files over 5 MB:** all are tracked. Together they are 62.1 MB, 81% of the folder.

  | File | Size | Tracked |
  |---|---:|---|
  | `fixed_random/gbsg_extreme_subgroups_fixed.html` | 8.39 MB | yes |
  | `fixed_random/extreme_subgroups_fixed.html` | 8.39 MB | yes |
  | `fixed_random.zip` | 6.89 MB | yes |
  | `extreme_subgroups_fixed.html` | 5.83 MB | yes |
  | `fixed_random/legacy/extreme_subgroups_fixed.html` | 5.83 MB | yes |
  | `extreme_subgroups.html` | 5.82 MB | yes |
  | `fixed_random/legacy/extreme_subgroups.html` | 5.82 MB | yes |
  | `fixed_random/extreme_subgroups_fixed_wrapper.html` | 5.73 MB | yes |
  | `fixed_random/extreme_subgroups.html` | 5.72 MB | yes |

- **`.qmd` files sharing a title or differing only by suffix:**
  - Title "Extreme Subgroup Effects Under a Uniform Treatment Benefit" (12 files):
    - root: `extreme_subgroups`, `_short`, `_v1`, `_v2`, `_fixed`
    - `fixed_random/`: `extreme_subgroups`, `_fixed`, `_fixed_wrapper`, `gbsg_extreme_subgroups_fixed`, `_slides`
    - `fixed_random/legacy/`: `extreme_subgroups`, `_fixed`

    The `_fixed` variants carry the subtitle "A Fixed-Baseline (Conditional-on-X) …". The others carry "A Simulation Study Based on Observed Trial Data".
  - Title "Fixed vs Random Baseline: A Design Comparison" (3 files): root, `fixed_random/` and `fixed_random/legacy/` `extreme_subgroups_design_comparison.qmd`.
  - Suffix-only families: `extreme_subgroups{,_short,_v1,_v2,_fixed}` at root; `extreme_subgroups_fixed{,_wrapper}` and `{,gbsg_}extreme_subgroups_fixed` in `fixed_random/`.
- **Byte-identical copies** (`shasum -a 256`):
  - `fixed_random/extreme_subgroups_fixed.qmd` = `fixed_random/gbsg_extreme_subgroups_fixed.qmd`. Their HTMLs differ by 28 bytes.
  - Root `extreme_subgroups.qmd`/`.html` and `extreme_subgroups_design_comparison.qmd`/`.html` each equal their `fixed_random/legacy/` copy. Root `extreme_subgroups_fixed.qmd` also equals its legacy copy. The two `extreme_subgroups_fixed.html` files are the same size but not identical.
  - `fixed_random.zip` (2026-08-12) holds the Aug-11 versions of `extreme_subgroups{,_fixed,_design_comparison}.{qmd,html}` plus `__MACOSX` resource forks: 14.6 MB uncompressed. Its sizes match the root and `legacy/` copies.
  - Outside the folder, `extreme_subgroups_deck.{pdf,pptx}` and `extreme_subgroups_section_6_6_4.{pdf,pptx}` are byte-identical to the copies in `quarto/guides/`. `build_deck.js` differs from `quarto/guides/build_deck.js`: they are different builders that write different `.pptx` outputs.

---

## D5. Vignette linkage

The forestsearch site publishes **no file from `quarto/extreme_subgroups/`**.

- **Published document:** `vignettes/articles/extreme_subgroups.qmd`, 2101 lines, last changed in `d6a6ba0d` on 2026-10-17.
- **Mechanism:** it is a **pkgdown article** (`articles/extreme_subgroups`), rendered by `pkgdown::build_article()`. It is listed in the `_pkgdown.yml` navbar (L57) and the "Simulation Studies" articles index (L560), and deployed through `.github/workflows/pkgdown.yaml`.
- **Relation to the folder:** the article is a separately maintained, adapted copy of the fixed-baseline document. `dev/tasks/TASK_extreme_subgroups_qmd_replace_2026-10-17.md` describes it as "Larry's fixed-baseline document (`extreme_subgroups_fixed.qmd`), adapted for pkgdown". Today it differs from root `extreme_subgroups_fixed.qmd` by 255 diff lines and from `fixed_random/extreme_subgroups_fixed.qmd` by 974.
- **What the site does not use:** no HTML is copied from the folder into `docs/` or `vignettes/`, and no site page links into the folder. The only outward link is the reverse direction: `quarto/guides/build_deck.js:493` points a deck at the published pkgdown URL.
- **Package build:** `.Rbuildignore` `^quarto$` keeps the whole folder out of the package tarball.
