# REPORT — ACTG 175 article, stage 0: provenance, reproduction, timing (2026-10-30)

Task: `dev/tasks/cc_task_articles_finish_2026-10-30.md`, Part 2. Branch `feature/articles-finish` (BASE_SHA `0bdc2036`).
Record: `dev/reports/actg175_article_stage0_2026-10-30/`.

**Outcome (v2): Part 2 complete; every gate passes.** Under v1, Gates A and B passed and Table S17 reproduced, all to
machine precision, and the part stopped at Gate C: 2026a's FS multiplier upper limit in Ĥ is **5.09 against Table S3's
5.08**. v2 (`dev/tasks/cc_task_articles_finish_2026-10-30_v2.md`) added the transplant check §2.2.4 (all 15 chunks
identical to `43b051b6`, so M-S3-FS was not re-run) and restated Gate C with a 0.01 tolerance on the interval limits:
**Gate C passes** (§11). M-S3-DINA and M-S3-GRF then reproduce Table S3 at 2 dp and the 0.2.0 payload to ~1e-14,
multiplier rows included. The three bootstrap rungs ran (3.1 / 1.3 / 4.0 s per replicate, single core). The v1 sections
below are kept as recorded; §11 is the v2 section.

## 1. Provenance

- forestsearch: installed **0.4.1** (built 2026-10-30 16:08:45 UTC, R 4.6.1); `git diff c2967362 HEAD -- R/ DESCRIPTION NAMESPACE inst/` empty.
  `speff2trial` 1.0.5 installed. All rungs `plan = "sequential"`, single process, on this Linux machine (128 logical cores).
- `fs-glms-interpretable`: HEAD `19c6ac41c1198dbc513f5c6978d0d6c5df1fc157` (2026-10-26), `git status --short` empty. Read only.
- `fs-post-selection`: HEAD `bc82bc32c8a6ef66a6fc252992278f102b44b819` (2026-09-15), `git status --short` empty. Read only.

## 2. Inventory A — 2026b, Tables S16 and S17 (§2.1)

- **Fragment:** `fs-glms-interpretable/manuscript/_11_applications_supplement.qmd`. Table S16 = chunk `tbl-supp-app-five`
  (lines 147–152), `knitr::kable(app_tbl_five(APP))`. **Table S17 is a hand-typed markdown table** (lines 162–174),
  not drawn from a payload at render. The prose with 3.584 / 1.773 / 0.686 / 71.1 / 22.1 is line 154.
- **Extract read:** `APP <- app_load(.app_dir)` over `dev/payloads/applications/`, by `manuscript/_R/app_tables.R`.
  ACTG 175 rows are in `applications_metrics.csv` (rows `actg175_binary` / `actg175_continuous` × `H` / `Hc`), with the
  whole-trial row in columns `whole_trial_*`. Extract built 2026-09-23 23:56:51 PDT (csv commit `664c907`).
- **Extract builder:** `dev/payloads/R/applications_metrics.R` (commit `3841897`). It reads the committed analysis payloads.
  The whole-trial row comes from each payload's `extras$reference_or` / `reference_md` (builder lines 84–86, 155), i.e.
  from the documents' own `reference-effects` chunk.
- **Producing documents** (§2.1.2): `quarto/actg175/analysis_actg175_binary_mr.qmd` → `analysis_actg175_binary_mr_payload.rds`
  (md5 `01c7f1b6…`, built 2026-09-23 23:54:32, forestsearch 0.3.5.9000, commit `d97ca1f`), and
  `analysis_actg175_continuous_mr.qmd` → `…_continuous_mr_payload.rds` (md5 `a3da7001…`, built 23:56:27, 0.3.5.9000,
  commit `935ae87`). Both payload md5s equal the extract's `source_md5`. **Table S17's declaration calibration is in
  the same two documents** (chunks `declcal` … `declcal-decision`), stored in payload element `declaration_calibration`.
  *GATE (§2.1.2): PASS.* The extract carries every Table S16 cell of both columns, including the whole-trial row, at
  printed precision.
- Copied to `manuscript_source_2026b/`: the fragment, `app_tables.R`, `applications_metrics.R`, `applications_metrics.csv`, `COLUMNS.md`.

## 3. Inventory B — 2026a, Table S3 (§2.2)

- **Fragment:** `fs-post-selection/manuscript/_actg175_supp.qmd` (included by the supplement, line 1128). It reads
  `_payloads/analysis_actg175_binary_multimethod_fixed_family_payload.rds` (md5 `04903a67e23a0a45c9d81dc5d9aaaa17`,
  1,926 bytes, committed `61e4ab7`). Caption string at line 84.
- **Manifest entry** (`manuscript/payload_manifest.qmd` lines 84–87): consumer `_actg175_supp (supp S7)`; "keep —
  post-audit (built 08-16); fd1ce74f: full-pool run diverges from published table"; generator forestsearch
  `quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd`.
- **Payload:** built 2026-08-16 09:50:40 PDT (16:50:40 UTC) on Darwin, R 4.5.2, **forestsearch 0.2.0**, 14 cores,
  B = 1000, MR draws 5000. Its `table` carries Table S3 at full precision (FS Ĥ: n 72, naive 3.584, FB 2.427, MR 1.794 (0.633, 5.082)).
  *GATE (§2.2.2): PASS.*
- **Commit as of the build:** the last commit touching the document before 09:50:40 is `fd1ce74f` (2026-08-15 21:33).
  **The render ran the text later committed at `43b051b6`** (2026-08-16 12:07). That commit's HTML, rendered by the
  same run, shows `fs_selection_rule <- "neighborhood"` and "2.43 (0.80", where `fd1ce74f` has `"both"`.
  `43b051b6` is identical to HEAD for this file. Transplants are taken from `43b051b6`.
- **Bootstrap in the same render: yes.** The payload's `meta$timings` records `boot_sec` per identifier from that run
  (FS 905 s, DINA 62.6 s, GRF 342 s at 14 workers), beside `mr_sec` (0.545 / 0.151 / 0.297 s).
- Copied to `manuscript_source_2026a/`: the fragment, the manifest and the payload.

## 4. Transplants (§2.1.3, §2.2.3)

Verbatim chunk bodies with a header (source path, chunk names, fence lines, sibling HEAD / commit). Script:
the scratchpad `transplant.py`, used as given; nothing in a transplant is edited.

| file | source chunks |
|---|---|
| `calls_2026b/binary_data_prep.R` | `label-map` 106–169, `setup` 177–212, `reference-effects` 219–255 |
| `calls_2026b/binary_fit.R` | `fit` 264–297 |
| `calls_2026b/binary_declcal.R` | `declcal`, `declcal-table`, `declcal-observed`, `declcal-decision` |
| `calls_2026b/binary_assembly.R` | `selection`, `extract`, `payload` (the Table S16 cell assembly) |
| `calls_2026b/continuous_*.R` | the same chunks of the continuous document |
| `calls_2026a/data_prep.R` | `setup`, `parallel-setup` (all the identification knobs), `load-packages`, `analysis-roles`, `data-treat`, `data-outcomes`, `data-covariates` |
| `calls_2026a/identifier_{fs,dina,grf}.R` | `forestsearch-main`; `dina-config` + `fs-dina-select`; `fs-grf-select` |
| `calls_2026a/bootstrap_calls.R` | `bootstrap`, `fs-dina-bootstrap`, `fs-grf-bootstrap` (`nb_boots = NB = 1000`, default seed) |
| `calls_2026a/table_assembly.R` | `export-section6-2-payload` |

The producing documents themselves are in `producing_docs/` (2026a's is the `43b051b6` text).

## 5. Settings and package facts (§2.3)

`settings.md` gives one table per analysis against the GBSG Table 1 calls and Examples 2–3. Points to note: 2026a
uses `effect_neighborhood = 0.10` and `use_twostage = TRUE` with `fs.splits = 500`, and its MR has only `draws = 5000`.
2026b binary matches the GBSG call except for the outcome-specific arguments and `max_subgroups_search = Inf`. 2026b
continuous adds `use_twostage = TRUE`, factor-coded binary covariates and a shorter MR argument list.

`package_facts.md` covers the element path of every Table S16 cell, what `print()`/`summary()` show and omit,
`fs_declaration_calibration()` with a negative `c0` on the MD path (accepted; D-cont ran it), and the
bootstrap/`mr_estimates_table()` interface (from code only, since the bootstrap rungs did not run).

## 6. The ladder (§2.4) — v1

Driver `ladder.R` (evaluates the transplants verbatim; overrides `quiet = TRUE`, `details = FALSE`, `plot.sg = FALSE`,
`plot.grf = FALSE`, `parallel_args = list(plan = "sequential")` by rewriting the call), run as
`timeout 1200 /usr/bin/time -v Rscript ladder.R <rung>`; comparisons by `compare.R`. Logs `logs/rung_<rung>.log`,
results `rungs/`, summary `ladder.csv`.

| rung | exit | wall (s) | fit / step (s) | peak RSS | outcome |
|---|---|---|---|---|---|
| M-bin | 0 | 65.9 | 64.2 (multiplier 11.4; field 17.1 + 4.7) | 1.62 GB | **GATE A PASS**: 27/27 Table S16 OR cells at 3 dp; region by membership; vs payload all 27 `all.equal` (1e-8), max \|Δ\| 3.6e-15 |
| M-cont | 0 | 101.2 | 99.4 (multiplier 19.6; field 23.4 + 7.3) | 2.23 GB | **GATE B PASS**: 24/24 MD cells at 1 dp; region by membership; vs payload 24/24, max \|Δ\| 2.6e-13 |
| D-bin | 0 | 1.7 | 0.2 | 0.27 GB | recorded: all 7 Table S17 OR cells match (T 2.494 / 2.494, κ̂ 2.889, implied 0.996, FŴ 0.838, claim κ̂ 3.611, FŴ 0.989); admits nothing; c0 table `all.equal` payload, \|Δ\| 0 |
| D-cont | 0 | 1.9 | 0.4 | 0.26 GB | recorded: all 7 MD cells match (1.701 / 2.043, 2.952, 0.997, 0.909, 3.625, 0.997); admits nothing; payload \|Δ\| 0 |
| M-S3-FS | 0 | 53.6 | 52.0 (MR 11.2) | 1.22 GB | **GATE C FAIL**: see below |
| M-S3-DINA, M-S3-GRF, B-S3-FS/-DINA/-GRF | — | — | — | — | not run under v1 (Part 2 stopped); run under v2, §11 |

Both 2026b regions print in the same label order as the extract (`!{cd40 <= 380} & !{wtkg <= 86}`;
`!{cd40 <= 507} & {gender}`). The documents' own anchor checks in `selection` passed (`selection chunk: ok`).

**GATE C detail (M-S3-FS against Table S3, 2 dp):**

| cell | this run (0.4.1) | Table S3 | payload (0.2.0) | match |
|---|---|---|---|---|
| region | `!{cd40 <= 380} & !{wtkg <= 86}` (membership = `wtkg > 86 & cd40 > 380`) | FS `wtkg > 86 & cd40 > 380` | `!{wtkg <= 86} & !{cd40 <= 380}` | yes (label order differs) |
| n | 72 / 1011 | 72 / 1011 | 72 / 1011 | yes |
| Ĥ naive | 3.5844 (1.3145, 9.7739) | 3.58 (1.31, 9.77) | identical (\|Δ\| ≤ 9e-15) | yes |
| Ĥ multiplier | 1.7945 (0.6321, **5.0947**) | 1.79 (0.63, **5.08**) | 1.7941 (0.6334, 5.0821) | **upper: 5.09 ≠ 5.08** |
| Ĥᶜ naive | 0.5868 (0.4539, 0.7585) | 0.59 (0.45, 0.76) | identical | yes |
| Ĥᶜ multiplier | 0.6164 (0.3699, 1.0271) | 0.62 (0.37, 1.03) | 0.6162 (0.3699, 1.0264) | yes |
| selection bias | 0.69221 | — | 0.69149 | (recorded) |
| re-selection family | 2794 | — | 2794 | (recorded) |

Printed-precision 13 / 14; against the payload at 1e-8: naive rows and counts equal, every MR quantity differs
(max \|Δ\| 0.0126, on the Ĥ upper limit). Stop record (§2.4: full `sg.harm`, counters, settings in force):
`logs/gateC_stop_record.txt`. In summary: `sg.harm` = `q7.0 "!{cd40 <= 380}"`, `q19.0 "!{wtkg <= 86}"`. Filter
counts: 3240 evaluated → 3090 variance → 3090 prevalence → 2975 redundancy → 2666 events → 2584 size → 2584 fittable →
160 clearing OR ≥ 1. Consistency: 137 evaluated, 5 passed. MR settings in force: `ci_method = "field"` (the 0.4.1
forest-search default), `confirm_rule = "point"`, `reselection = "effMaxSG"`, `selection_rule = "neighborhood"`,
`multiplier = "poisson"`, `draws = 5000`, `ij_draws = 4972`. `mr_inference_args` as received:
`list(draws = 5000L, family_native_neighborhood = NULL)`.

## 7. Budget arithmetic (§2.4) — v1 (superseded by §11.4)

Measured here, single core: 2026b binary fit 64.2 s and continuous 99.4 s (each including MR with the field, of
which the multiplier stage is 11.4 / 19.6 s); declaration calibrations 0.2 / 0.4 s; 2026a FS fit 52.0 s (MR 11.2 s).
**The per-replicate bootstrap cost was not measured** (B-S3 rungs not run), so the B each identifier can afford in
60 / 120 / 300 s of single-core time cannot be computed from this run. For reference only, and not single-core: the
2026a payload records B = 1000 at 14 workers on a Mac taking FS 905 s, DINA 62.6 s and GRF 342 s (0.905 / 0.063 /
0.342 s of wall time per replicate at 14 workers). The memory note on worker contention says not to scale such
figures to other (n, W) settings. No recommendation.

## 8. Deviations

1. **Stop at Gate C**, as above; §2.4's remaining five rungs and §2.5's commit and `~/Downloads` copy not done.
2. `ladder.R` wraps each document's `selection` chunk so that a failure of its own anchor `stop()` would be recorded
   rather than lose the cells. It did not fire: every `selection` chunk returned `ok`.
3. Rung M-S3-FS evaluated the 2026a `setup` chunk (`knitr::opts_chunk$set`, `timings <- list()`) in the
   transplant, because the identifier chunks write `timings$…`. The chunk is in `calls_2026a/data_prep.R`, verbatim.

## 9. OPEN ITEMS (findings only)

1. **Table S3's FS multiplier row under 0.4.1**: the Ĥ upper limit is 5.09 (payload 5.08 under 0.2.0); the estimate
   and complement agree at 2 dp. The payload predates the field and the `ci_method = "field"` default (0.4.1 records
   `ij_draws = 4972`, `mr_sec` 11.2 s against the payload's 0.545 s). Whether part B re-states the row under 0.4.1 or
   reproduces 0.2.0 is Larry's call. DINA and GRF against Table S3 are unmeasured.
2. **Table S17 is typed into `_11_applications_supplement.qmd`, not drawn from a payload.** Its values match the
   payloads' `declaration_calibration` exactly (D-bin, D-cont).
3. **2026a's render ran uncommitted text**: the payload (09:50) predates the commit of the text it ran (`43b051b6`,
   12:07). The manifest names `fd1ce74f`, which has `fs_selection_rule = "both"`.
4. **`print()`/`summary()` on GLM fits label the consistency table's `hr` column "HR" / "Hazard ratio"**
   (`R/forestsearch_methods.R:277`, `:499`). It prints 1.277 (the log OR, = log 3.584) on the binary fit and 71.096
   (an MD) on the continuous fit. Neither method prints the corrected point estimate, the unadjusted estimates, the
   complement's corrected estimate or the selection bias (`package_facts.md` §2).
5. The 2026a document's comment block cites "published supplement Table S11 … MR 1.80" from an earlier capped run.
   The current Table S3 (1.79) matches this payload, not that comment.

## 10. Record listing

v1: not copied to `~/Downloads` (§2.5 was not reached). v2: see §11.6.

## 11. v2

### 11.1 §2.2.4 — the transplants against `43b051b6`

`logs/s3_transplant_vs_43b051b6.diff`: `git show 43b051b6:quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd`
(md5 `ee724a298a4e7bdb5a22d96f6c80ba92`) was parsed by chunk label, independently of the fence lines recorded in the
headers, and each chunk body was compared with its transplant. **15 chunks compared, 0 differ** (`bootstrap`,
`fs-dina-bootstrap`, `fs-grf-bootstrap`; `setup`, `parallel-setup`, `load-packages`, `analysis-roles`, `data-treat`,
`data-outcomes`, `data-covariates`; `dina-config`, `fs-dina-select`; `forestsearch-main`; `fs-grf-select`;
`export-section6-2-payload`). The recorded fence lines agree too. Nothing was re-transplanted, and **M-S3-FS was not
re-run**; Gate C is evaluated on v1's M-S3-FS result (`rungs/M-S3-FS_result.rds`, `_compare.csv`).

### 11.2 Gate C as restated — PASS

`logs/gateC_v2.txt`. Region by membership identical to Table S3's (label order differs); n 72 / 1011; both naive rows
exact at 2 dp; both multiplier point estimates exact at 2 dp (1.79, 0.62). Interval limits, compared at 2 dp as integer
hundredths:

| cell | 0.4.1, full precision | 0.4.1, 2 dp | Table S3 | Δ (2-dp units) | payload (0.2.0) | 0.4.1 − payload |
|---|---|---|---|---|---|---|
| Ĥ MR lower | 0.6320850 | 0.63 | 0.63 | 0 | 0.6333541 | −0.00127 |
| Ĥ MR upper | 5.0947467 | 5.09 | 5.08 | **1** | 5.0821221 | **+0.0126** |
| Ĥᶜ MR lower | 0.3698888 | 0.37 | 0.37 | 0 | 0.3698846 | +4.2e-06 |
| Ĥᶜ MR upper | 1.0271492 | 1.03 | 1.03 | 0 | 1.0264276 | +0.00072 |

**How the tolerance was read (for Larry to confirm):** "within 0.01 in each interval limit" was applied to the
printed 2-dp values, as §2.4's own sentence states ("a one-unit difference in the second decimal of a limit is within
that drift"). At full precision the Ĥ upper limit is 0.0147 above the printed 5.08 and 0.0126 above the payload's
5.0821, both beyond 0.01. A full-precision reading of the tolerance would fail Gate C again on this one cell.

### 11.3 The remaining rungs

| rung | exit | wall (s) | fit / step (s) | peak RSS | outcome |
|---|---|---|---|---|---|
| M-S3-DINA | 0 | 4.5 | 2.97 (MR 1.04) | 0.49 GB | recorded: 14/14 Table S3 cells at 2 dp; region `{preanti >= 849.4} & {cd40 >= 338}` by membership, n 98 / 985; vs payload 16/16 `all.equal` (1e-8), max \|Δ\| 1.2e-14; family 147, selection bias 0.708 |
| M-S3-GRF | 0 | 21.0 | 18.03 (MR 5.43) | 0.82 GB | recorded: 14/14 at 2 dp; region `{wtkg > 84.3696} & {cd40 > 368.2}` by membership, n 95 / 988; vs payload 16/16, max \|Δ\| 1.0e-14; family 1289, selection bias 0.726 |
| B-S3-FS | 0 | 20.0 | 18.33 | 0.39 GB | 5 replicates, `boot_status ok`; per replicate (s) 2.87, 2.70, 3.91, 2.79, 3.05 |
| B-S3-DINA | 0 | 11.0 | 9.36 | 0.39 GB | runs; 0.94, 0.62, 3.06, 0.65, 1.02 |
| B-S3-GRF | 0 | 25.0 | 23.01 | 0.59 GB | runs; 4.17, 4.16, 3.98, 4.00, 3.85 |

DINA and GRF Ĥ / Ĥᶜ, 0.4.1 (= Table S3): DINA naive 1.94 (0.84, 4.47), MR 0.95 (0.29, 3.10); Ĥᶜ 0.60 (0.46, 0.78),
MR 0.63 (0.38, 1.06). GRF naive 2.25 (0.97, 5.20), MR 1.10 (0.43, 2.80); Ĥᶜ 0.58 (0.45, 0.76), MR 0.62 (0.37, 1.04).
Comparisons: `logs/compare_M-S3-DINA.txt`, `logs/compare_M-S3-GRF.txt`, `rungs/*_compare.csv`.

The bootstrap call as run (`logs/boot_timing_v2.txt`): the document's `fs.est`, `show_three = FALSE`, `details = FALSE`,
`digits = construction_digits`, with `nb_boots = 5L` and `parallel_args = list(plan = "sequential")` substituted; no
`seed` argument, so the default `seed = 8316951` (the document's own) and `mr_in_replicates = FALSE` (default) are in
force. The bias-corrected H₂ intervals at B = 5 (FS 2.20 (1.1e-06, 4.5e+06); DINA 0.49 (4.4e-04, 542); GRF 0.96 (0.30,
3.08)) are recorded only. B = 5 is a timing run and these are not estimates.

### 11.4 Budget arithmetic (single core, `plan = "sequential"`, this machine)

Fit wall-clock (identifier call with MR): FS 52.0 s (MR 11.2), DINA 2.97 s (MR 1.04), GRF 18.0 s (MR 5.43); 2026b
binary 64.2 s, continuous 99.4 s (§6).

| identifier | first replicate (s) | mean of 5 (s) | total / 5 (s) | B in 60 s | B in 120 s | B in 300 s |
|---|---|---|---|---|---|---|
| FS | 2.87 | 3.06 | 3.67 | 16 (19) | 32 (39) | 81 (97) |
| DINA | 0.94 | 1.26 | 1.87 | 32 (47) | 64 (95) | 160 (238) |
| GRF | 4.17 | 4.03 | 4.60 | 13 (14) | 26 (29) | 65 (74) |

B = ⌊budget / (total / 5)⌋, where total is the bootstrap call's own elapsed time including its setup and summary.
In parentheses: ⌊budget / mean per-replicate time⌋ (`results$tmins_iteration`). `logs/budget_v2.txt`. No recommendation.

### 11.5 Findings added in v2 (OPEN ITEMS)

6. **The 0.4.1-vs-0.2.0 multiplier drift is forest-search-only.** DINA's and GRF's MR rows reproduce the 0.2.0 payload
   to ~1e-14 under 0.4.1; FS's MR rows differ by up to 0.0126 (§11.2). The rows share the same `draws = 5000`. Not
   investigated (no `R/` reading was asked for); OPEN ITEM 1 stands with this added.
7. **A bootstrap replicate costs far less than the identifier fit**: FS 3.1 s per replicate against a 52 s fit, GRF 4.0
   against 18 s. The replicates run with `mr_in_replicates = FALSE`; where the rest of the difference goes was not
   measured.
8. The DINA replicate times are uneven (0.62–3.06 s); with B = 5 the mean is sensitive to the one slow replicate.

### 11.6 Commit and `~/Downloads` copy

Committed by path, alone: `dev: report — ACTG 175 article stage 0 (2026-10-30)`. Copied to
`~/Downloads/actg175_article_stage0_2026-10-30/`: this report, `ladder.csv`, `settings.md`, `package_facts.md`,
`calls_2026b/`, `calls_2026a/`, `manuscript_source_2026b/`, `manuscript_source_2026a/`, `producing_docs/` (the
producing documents' `.qmd`). The listing below was taken after the first copy; the report was then re-copied with this
listing appended, so its own size in the listing is the pre-append size.

```
$ ls -la ~/Downloads/actg175_article_stage0_2026-10-30/
total 88
drwxrwxr-x 7 larryleon larryleon  4096 Oct 30 17:36 .
drwxr-xr-x 8 larryleon larryleon 20480 Oct 30 17:36 ..
drwxrwxr-x 2 larryleon larryleon  4096 Oct 30 17:36 calls_2026a
drwxrwxr-x 2 larryleon larryleon  4096 Oct 30 17:36 calls_2026b
-rw-rw-r-- 1 larryleon larryleon   894 Oct 30 17:36 ladder.csv
drwxrwxr-x 2 larryleon larryleon  4096 Oct 30 17:36 manuscript_source_2026a
drwxrwxr-x 2 larryleon larryleon  4096 Oct 30 17:36 manuscript_source_2026b
-rw-rw-r-- 1 larryleon larryleon  5879 Oct 30 17:36 package_facts.md
drwxrwxr-x 2 larryleon larryleon  4096 Oct 30 17:36 producing_docs
-rw-rw-r-- 1 larryleon larryleon 19452 Oct 30 17:36 REPORT_actg175_article_stage0_2026-10-30.md
-rw-rw-r-- 1 larryleon larryleon  7151 Oct 30 17:36 settings.md

$ ls -la calls_2026b/
total 76
drwxrwxr-x 2 larryleon larryleon  4096 Oct 30 17:36 .
drwxrwxr-x 7 larryleon larryleon  4096 Oct 30 17:36 ..
-rw-rw-r-- 1 larryleon larryleon 10008 Oct 30 17:36 binary_assembly.R
-rw-rw-r-- 1 larryleon larryleon  7715 Oct 30 17:36 binary_data_prep.R
-rw-rw-r-- 1 larryleon larryleon  6627 Oct 30 17:36 binary_declcal.R
-rw-rw-r-- 1 larryleon larryleon  2312 Oct 30 17:36 binary_fit.R
-rw-rw-r-- 1 larryleon larryleon 11055 Oct 30 17:36 continuous_assembly.R
-rw-rw-r-- 1 larryleon larryleon  8350 Oct 30 17:36 continuous_data_prep.R
-rw-rw-r-- 1 larryleon larryleon  6638 Oct 30 17:36 continuous_declcal.R
-rw-rw-r-- 1 larryleon larryleon  2746 Oct 30 17:36 continuous_fit.R

$ ls -la calls_2026a/
total 56
drwxrwxr-x 2 larryleon larryleon  4096 Oct 30 17:36 .
drwxrwxr-x 7 larryleon larryleon  4096 Oct 30 17:36 ..
-rw-rw-r-- 1 larryleon larryleon  2371 Oct 30 17:36 bootstrap_calls.R
-rw-rw-r-- 1 larryleon larryleon 15619 Oct 30 17:36 data_prep.R
-rw-rw-r-- 1 larryleon larryleon  3235 Oct 30 17:36 identifier_dina.R
-rw-rw-r-- 1 larryleon larryleon  4126 Oct 30 17:36 identifier_fs.R
-rw-rw-r-- 1 larryleon larryleon  2852 Oct 30 17:36 identifier_grf.R
-rw-rw-r-- 1 larryleon larryleon 10131 Oct 30 17:36 table_assembly.R

$ ls -la manuscript_source_2026b/
total 96
drwxrwxr-x 2 larryleon larryleon  4096 Oct 30 17:36 .
drwxrwxr-x 7 larryleon larryleon  4096 Oct 30 17:36 ..
-rw-rw-r-- 1 larryleon larryleon 19447 Oct 30 17:36 _11_applications_supplement.qmd
-rw-rw-r-- 1 larryleon larryleon  9856 Oct 30 17:36 applications_metrics.csv
-rw-rw-r-- 1 larryleon larryleon 30156 Oct 30 17:36 applications_metrics.R
-rw-rw-r-- 1 larryleon larryleon 11799 Oct 30 17:36 app_tables.R
-rw-rw-r-- 1 larryleon larryleon 10626 Oct 30 17:36 COLUMNS.md

$ ls -la manuscript_source_2026a/
total 40
drwxrwxr-x 2 larryleon larryleon  4096 Oct 30 17:36 .
drwxrwxr-x 7 larryleon larryleon  4096 Oct 30 17:36 ..
-rw-r--r-- 1 larryleon larryleon 10778 Oct 30 17:36 _actg175_supp.qmd
-rw-r--r-- 1 larryleon larryleon  1926 Oct 30 17:36 analysis_actg175_binary_multimethod_fixed_family_payload.rds
-rw-r--r-- 1 larryleon larryleon 12724 Oct 30 17:36 payload_manifest.qmd

$ ls -la producing_docs/
total 216
drwxrwxr-x 2 larryleon larryleon  4096 Oct 30 17:36 .
drwxrwxr-x 7 larryleon larryleon  4096 Oct 30 17:36 ..
-rw-rw-r-- 1 larryleon larryleon 58863 Oct 30 17:36 analysis_actg175_binary_mr.qmd
-rw-rw-r-- 1 larryleon larryleon 89928 Oct 30 17:36 analysis_actg175_binary_multimethod_fixed_family_43b051b6.qmd
-rw-rw-r-- 1 larryleon larryleon 57696 Oct 30 17:36 analysis_actg175_continuous_mr.qmd
```
