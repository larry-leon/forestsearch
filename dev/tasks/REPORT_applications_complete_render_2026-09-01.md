# REPORT — complete applications at 0.3.5: LOO and K-fold exercised, main analysis guarded — 2026-09-01

**Task:** `dev/tasks/cc_task_applications_complete_render_2026-09-01.md` (committed on arrival, `55048e73`).
**Verdict up front: the complete documents are now canonical, and the headline analysis did not move.** The §6 guard is CLEAN — every shared COMPARED leaf `identical()` at full precision against the verified `_payloads_2026-09-01` reference; completeness added exactly the expected NEW-CONTENT (LOO/CV results) and nothing else changed. Projection gate passed with an order of magnitude to spare (realized CV/LOO cost: 33 s + 145 s wall on 102 workers, against the 3-hour ceiling).

## 1 Provenance and the §2 census

### 1.1 Provenance (§1 gate — PASS)

| item | value |
|---|---|
| branch | `feature/glm-extension`; tree clean at start; pre-task HEAD `164f7e07` |
| installed forestsearch | `0.3.5`; `vi.grf.min` default `NULL` |
| 09-01 comparison commits in log | `84d34a08`, `1ada23eb`, `7fef766b`, `86cbd5bb` — all present |
| park/manifest follow-on commit | **`164f7e07`** (parked comparison_continuous.rds + sha256 manifest + ignore rule) |
| commits this task | (1) `55048e73` task doc alone; (2) `f815135f` transplanted scripts; (3) `6ad5a724` flag flips; (4) the §6 final commit (section 6) |

### 1.2 Flag census (§2.1)

Quoted at pre-flip HEAD (`grep -n -E 'run_(loo|cv|cv_kfold|kfold|mr)\s*<-|^(K|Ksims)…'`):

- **FLAGGED — `gbsg/analysis_gbsg_survival_frozen_family.qmd`:**
  - `90:run_loo    <- FALSE     # leave-one-out (N-fold) honest out-of-sample subgroup estimate`
  - `98:run_cv_kfold <- FALSE    # repeated 10-fold cross-validation (stability diagnostics)` (with `99:Ksims       <- 50L`)
- **FLAGGED — `gbsg/analysis_gbsg_survival_multimethod.qmd`:**
  - `144:run_cv  <- FALSE    # K-Fold cross-validation (kfold-cv + diagnostics)` (with `148:Ksims    <- 50`, `149:Kfolds_cv <- 10`)
  - `145:run_loo <- FALSE    # Leave-one-out / N-fold OOB (oob-cv)`
  - (`180:run_mr <- TRUE` — exempt, already enabled, untouched)
- **UNFLAGGED — `actg175/template_actg175_continuous.qmd` and `actg175/analysis_actg175_continuous_compare_all.qmd`:** no LOO/CV machinery of any kind (no `run_*` assignments) — they render identically regardless and were not re-rendered.

### 1.3 Unit semantics and loop bounds (§2.2)

- **LOO** = `forestsearch_Kfold()` (`R/forestsearch_cross_validation.R:287`) at `Kfolds = nrow(df.analysis)` = **686** folds per document; each fold is a **full `forestsearch()` refit** on the n−1 training set, with the inner refit pinned `plan = "sequential", workers = 1` (the fold loop owns all parallelism; frozen_family qmd L91–95 documents exactly this). Both documents call it with `mr_in_replicates = FALSE` in effect (multimethod passes it explicitly; frozen_family gets the signature default).
- **K-fold CV** = `forestsearch_tenfold()` (`…:803`) at `sims = Ksims = 50`, `Kfolds = 10` = **500** training-fold refits per document (617-row training sets), same full-refit unit.
- frozen_family's LOO carries a **cache**: `loo_cache <- Sys.getenv("LOO_CACHE", unset = file.path(gh_dir, sprintf("cv_out_%s_%s_%s.rds", .cache_doc, .cache_focus, .cache_rule)))` with `gh_dir <- Sys.getenv("GH_DIR", unset = "../../../quarto/GuoHe")`; on a miss it computes and `saveRDS(cv_out, loo_cache)`. The resolved key is `cv_out_analysis_gbsg_survival_frozen_family_maxeff_neighborhood.rds` (focus from the fitted object = `"maxeff"`; rule = the `"neighborhood"` package default) — **no such file existed**, so LOO computed fresh; the only pre-existing cache (`cv_out_analysis_gbsg_survival_effMaxSG_hrMaxSG_neighborhood.rds`, tracked) belongs to a different document key and was not touched. Multimethod's LOO has no cache mechanism.

### 1.4 Seeding (§2.3) — reproducible, not a snapshot

Both mechanisms are **internally seeded at the entry points** (the documents pass no seeds, so the defaults are in force):

- `forestsearch_Kfold`: `seedit = 8316951L` (signature), `set.seed(seedit)` for fold assignment (`…:369`), and `.options.future = list(seed = TRUE)` (`…:466`) — parallel-safe L'Ecuyer streams per fold.
- `forestsearch_tenfold`: `seed = 8316951L` (signature), `set.seed(seed + 1000 * ksim)` per simulation (`…:958`), `.options.future = list(seed = TRUE)` (`…:952`).

The complete render is therefore a **reproducible artifact**. Corroborated three ways: the frozen_family qmd's own verification note ("**identical** per-subject predictions … at 5.2x the speed on 11 workers"); this render's LOO reproducing the 0.2.0-era baseline LOO verbatim (§5); and the probe fit reproducing the headline subgroup.

### 1.5 Output paths (§2.4)

Enumerated from source, all observed exactly: the two `.html` (tracked overwrites); the two document payloads at `gbsg/_payloads/<stem>/<stem>_payload.rds` (frozen_family's grown by `extras$cv`/`extras$loo`); and the **one new GuoHe file** `quarto/GuoHe/cv_out_analysis_gbsg_survival_frozen_family_maxeff_neighborhood.rds` (1,190 B; `CV_summary` length 686). Multimethod's CV/LOO write nothing to disk and — per its payload constructor — export **no CV/LOO fields** (extras = concordance/harm_rates only; `meta$timings` = fs/dina/grf blocks; `timings$kfold`/`timings$oob` stay document-local): its completeness lives in the HTML alone.

### 1.6 The 08-31 estimates (§2.5), quoted

Inventory report (`REPORT_applications_inventory_2026-08-31.md`): "LOO recompute cost (iii): **N = 686** for gbsg … So LOO(gbsg) = 686 × t_search, with fold-level parallelism over `loo_workers = n_cores` … Single-search time: **U — not stated in the 08-27 report and not measured in this task**"; per-document rows "≈ 4–5 min + 10-fold CV (U)" / "+ CV (U) + LOO". The probes below supply the U's.

### 1.7 Ordering and RNG (§2.6)

Both documents run their CV/LOO sections **after every main-analysis computation completes** (frozen_family: fit L256 → FB L353 → GH L384 → CV L419 → LOO L562 → summary table → payload; multimethod: three fit+bootstrap+MR blocks → CV L1974 → LOO L2137 → forest plot → payload). The CV/LOO entry points call `set.seed()` internally, which replaces the session RNG state — but nothing downstream of them draws random numbers (tables, a deterministic forest plot, the payload constructor), so enabling the flags **cannot reach back into the headline analysis**. The §6 guard confirms this empirically: every shared leaf identical.

## 2 Flag flips and transplant diffs

### 2.1 Flips (`6ad5a724` — exactly the four flag lines, `git diff --stat`: 2 files, 4 insertions, 4 deletions)

```
-run_loo    <- FALSE     # leave-one-out (N-fold) honest out-of-sample subgroup estimate
+run_loo    <- TRUE      # leave-one-out (N-fold) honest out-of-sample subgroup estimate
-run_cv_kfold <- FALSE    # repeated 10-fold cross-validation (stability diagnostics)
+run_cv_kfold <- TRUE     # repeated 10-fold cross-validation (stability diagnostics)
-run_cv  <- FALSE    # K-Fold cross-validation (kfold-cv + diagnostics)
+run_cv  <- TRUE     # K-Fold cross-validation (kfold-cv + diagnostics)
-run_loo <- FALSE    # Leave-one-out / N-fold OOB (oob-cv)
+run_loo <- TRUE     # Leave-one-out / N-fold OOB (oob-cv)
```

The symmetric reversal of `d1f0e88e` (08-31), alignment and comments preserved.

### 2.2 Transplants (`f815135f`) — named lines only (gate PASS)

- **Driver** (vs committed `render_driver_applications_0.3.5_2026-09-01.sh`): header comment; `SENTINEL` → `…_20260901_complete`; `LOG` default name; `DOCS` reduced to the two flagged gbsg documents; one gate-label string. Nothing else.
- **Comparator** (vs committed `compare_payloads_applications_0.3.5_2026-09-01.R` at `7fef766b`): header comment; output default name; pair set → the two gbsg pairs (ref `_payloads_2026-09-01`, new `_payloads_2026-09-01_complete`); side labels `A(08-31)/B(09-01)` → `A(reference)/B(complete)`; the **NEW-CONTENT** classification — three named rules: (a) a name present only on the complete side (previously STRUCTURAL "present in B only"), (b) a reference-`NULL` leaf that gained content (its mirror — complete-side content vanishing — stays a COMPARED difference), (c) a flag-off placeholder cell whose reference value literally encodes the disabled flag (`run_(loo|cv|cv_kfold) = FALSE`); the NEW-CONTENT output section; and the named-lines footer reduced to the gbsg documents. The stale comparison_continuous header note (not in this pair set) was dropped. All classes otherwise keep their `7fef766b` meanings.

## 3 Probe and projection (§4 gate — PASS)

Probe method (a truncated-bounds dry run — the fold unit is not directly callable): purl each document, evaluate top-down only until its main `fs` fit exists (stopping before bootstrap/GH/CV), then `forestsearch_tenfold(sims = 1, Kfolds = 10)` on a sequential 1-worker plan = ten real fold units. The multimethod probe fit reproduced the headline subgroup (`{er <= 0} & {size <= 35}`, Pcons 0.990). LOO units scale the fold unit by 685/617 ≈ 1.11 (linear-in-n).

| doc | fold unit (probe) | LOO 686 × unit × 1.11 | K-fold 500 × unit | seq CV/LOO | flags-FALSE wall | projected wall (conservative 10× parallel) |
|---|---|---|---|---|---|---|
| frozen_family | 0.719 s | 548 s | 360 s | 15.1 min | 3m10 | ≈ 6–7 min |
| multimethod | 3.437 s | 2618 s | 1719 s | 72.3 min | 10m44 | ≈ 20–22 min |
| **total** | | | | **87.4 min** | 13m54 | **≈ 27–29 min** |

Even the degenerate zero-parallelism bound (fully sequential CV/LOO + base walls ≈ **101 min**) is under the 3-hour ceiling → **proceed**. (Realized parallel efficiency beat the conservative assumption ~3×: §4 below.)

## 4 Render walls and staging

Both renders exit 0 (transplanted driver, sentinel `/tmp/fs_render_sentinel_20260901_complete`; every output `-newer` sentinel):

| doc | flags-FALSE (86cbd5bb render) | complete (this task) | Δ | realized CV/LOO cost inside it |
|---|---|---|---|---|
| frozen_family | 3:09.58 | **3:43.40** | +34 s | LOO **17.7 s** (686 folds, 102 workers, ≈31× vs seq projection) + K-fold **14.0 s** (payload `meta$timings`) |
| multimethod | 10:43.80 | **13:06.31** | +2m23 | K-fold CV **61.5 s** + LOO **83.2 s** (document timing table; ≈28–31× vs seq projections) |
| **total** | 13m53 | **16m50** | +2m57 | |

Staging census (`quarto/applications/gbsg/_payloads_2026-09-01_complete/`):

| file | bytes |
|---|---|
| `analysis_gbsg_survival_frozen_family/analysis_gbsg_survival_frozen_family_payload.rds` | 884,604 (was 3,116 flags-off — `extras$cv`/`extras$loo` added) |
| `analysis_gbsg_survival_frozen_family/cv_out_analysis_gbsg_survival_frozen_family_maxeff_neighborhood.rds` | 1,190 (the LOO cache, **moved** from `quarto/GuoHe/` — deviation D1) |
| `analysis_gbsg_survival_multimethod/analysis_gbsg_survival_multimethod_payload.rds` | 1,900 (no CV fields by construction, §1.5) |
| `html/analysis_gbsg_survival_frozen_family.html` | 2,590,112 |
| `html/analysis_gbsg_survival_multimethod.html` | 18,931,711 |

Readback: both payloads `forestsearch_version == 0.3.5`, `built_at` 15:25:42 / 15:38:39 (inside the render windows). **No file > 50 MB → nothing parked; no manifest entries needed.** Ignore-rule coverage check: `git check-ignore -v` on the hypothetical `_complete` comparison_continuous path matches the existing `.gitignore:30` rule (`_payloads*` glob) — **covered, no extension needed**. The `_complete` directory itself is not ignored (plain `git add`).

## 5 The §6 guard — CLEAN, and what completeness added

Comparator run against the committed `_payloads_2026-09-01` reference:

| pair | EXCLUDED-VOLATILE | STRUCTURAL | NONCOMP | **NEW-CONTENT** | **COMPARED differences** |
|---|---|---|---|---|---|
| frozen_family | 9 (six `meta/timings` leaves incl. `loo` NA→17.741 / `kfold` NA→13.967, `extras/mr/timing_seconds`, `built_at`, `forestsearch_version` 0.3.5→0.3.5) | 0 | 0 | **3** | **0 — identical** |
| multimethod | 12 (machine timestamp, nine timing leaves, `built_at`, version) | 0 | 0 | **0** | **0 — identical** |

**NEW-CONTENT inventory** (frozen_family — the complete rendering of what the flags gate):

1. `/table/Estimate[5]`: `— (run_loo = FALSE)` → **`2.22 (1.17, 4.22)`** — the LOO row of the corrections table, filled.
2. `/extras/cv`: NULL → the K-fold metrics gt table plus `identification_summary`, `pconsistency_distribution`, `fold_numeric_summary`, `original_agreement`.
3. `/extras/loo`: NULL → `row = {case "stable", hr "2.22 (1.17, 4.22)", n "75 (10.9%)"}`, `sens_metrics` (sens_H/sens_Hc/ppv_H/ppv_Hc all 1), `find_metrics` (Exact = 1), `cache` path.

Multimethod adds no payload fields (its CV/LOO content is HTML-only, §1.5). The named selection lines are identical on both sides: frozen `{pgr <= 32.5} & {er <= 0}`; multimethod `fs {er <= 0} & {size <= 35}` / `dina {grade3 >= 1} & {pgr <= 10}` / `grf {er <= 0}`. **Cross-vintage corroboration:** the fresh LOO reproduces the 0.2.0-era baseline verbatim — the 08-31 report §5.3 quotes the old document printing `LOO [stable]: out-of-sample harm HR 2.22 (1.17, 4.22) (n = 75 (10.9%))`; the 0.3.5 LOO prints the same line, and the doc's by-construction claim (LOO equals Naive 2.22 here) holds.

**Clean branch taken.** The §2.6 mechanism answer predicted exactly this: CV/LOO run after the main analysis and nothing downstream re-draws RNG.

## 6 Tracked-HTML end state, final commit, and the final gate

Per-document end state (all four now render-reproducible from committed source at 0.3.5):

| tracked document | vintage at HEAD after this task |
|---|---|
| `gbsg/analysis_gbsg_survival_frozen_family.html` | **0.3.5 complete** (this render, in place) |
| `gbsg/analysis_gbsg_survival_multimethod.html` | **0.3.5 complete** (this render, in place) |
| `actg175/template_actg175_continuous.html` | **0.3.5** (promoted by copy from committed `_payloads_2026-09-01/html/`; `cmp` byte-identical to the bundle) |
| `actg175/analysis_actg175_continuous_compare_all.html` | **0.3.5** (promoted likewise) |

Final commit contents: the `_complete` directory, the four tracked HTML, this report. Final-gate outputs: `git diff --stat HEAD -- '*_payloads_2026-08-31*' '*_payloads_2026-09-01/*'` **empty** (no pre-existing payload directory modified); `git status --porcelain quarto/GuoHe` **empty** (the new cache was moved out, the tracked old cache untouched); porcelain at close clean (quoted in the close-out); nothing deleted anywhere.

## 7 Ten-line verdict

1. **The complete documents are now canonical, and the headline analysis did not move** — the guard compared every shared payload leaf against the just-verified 0.3.5 reference and found all of them `identical()` at full precision.
2. Completeness added exactly three payload fields, all in frozen_family: the filled LOO table row, `extras$cv`, `extras$loo` — enumerated as NEW-CONTENT, nothing reclassified away.
3. The LOO out-of-sample harm HR is `2.22 (1.17, 4.22)`, n = 75 (10.9%), sens/ppv all 1 — reproducing the 0.2.0-era baseline LOO verbatim across two package generations.
4. Both mechanisms are internally seeded (8316951-derived), so the complete render is a reproducible artifact, not a snapshot.
5. The compute panic was unwarranted at 0.3.5 on this host: 686-fold LOO costs 17.7 s (frozen) / 83.2 s (multimethod) on 102 workers; the whole completeness increment was +2m57 of render wall against a 3-hour ceiling.
6. The four flag flips are committed configuration (`6ad5a724`), reversing `d1f0e88e`; complete-document source now reproduces its rendered output.
7. All four tracked application documents now sit at a coherent 0.3.5 vintage (two complete renders, two promoted verified renders).
8. frozen_family's LOO cache was moved into the `_complete` bundle to keep GuoHe clean; a re-render recomputes LOO in ~18 s, or the cache can be copied back (D1).
9. Multimethod's payload records no CV/LOO by construction — a payload-schema fact worth knowing before anyone expects the guard to see its CV results (they are HTML-only).
10. No `R/` edit, nothing deleted, no push; deviations below.

## 8 Deviations and notes

- **D1 — GuoHe cache moved into the bundle.** The frozen_family LOO wrote its keyed cache to `quarto/GuoHe/` (the document's designed behavior, enumerated in §1.5). The task's writes policy confines new outputs to the `_complete` directories, so after the render the cache was **moved** (never deleted) to `_complete/analysis_gbsg_survival_frozen_family/`, leaving GuoHe byte-identical to HEAD. Consequence: a future frozen_family render recomputes LOO (~18 s, seeded, identical numbers) instead of cache-loading; copying the committed file back to `quarto/GuoHe/` restores cache-hit behavior. Larry can choose to track it there instead.
- **D2 — one `_complete` directory, not two.** The unflagged actg175 documents were not re-rendered (per §5) and produced nothing, so no `actg175/_payloads_2026-09-01_complete/` exists; the task's "both `_complete` directories" collapses to the gbsg one.
- **D3 — transient `Rplots.pdf`** (44,341 B, `quarto/applications/gbsg/`, mtime 15:21 during the frozen render — a CV-plotting device byproduct): this task's own byproduct, removed at close; nothing pre-existing touched.
- **D4 — probe scratch uncommitted.** The probe driver and its output live in the session scratchpad (the task's commit list names only the task doc, transplanted scripts, flips, and the final commit); its numbers are §3.
- **D5 — placeholder-cell rule.** NEW-CONTENT rule (c) of §2.2 extends the task's "fields present only on the complete side" to the one shared cell whose reference value literally encodes the disabled flag (`/table/Estimate[5]`). Without it that cell would have been a COMPARED difference by the letter of the taxonomy while being precisely the content the flag gates; the rule is anchored to the literal `run_* = FALSE` text so it can never absorb a numeric drift.
- **D6 — `forestsearch_version` enumerated though equal** (0.3.5 both sides) — the volatile list enumerates by name, values printed; no information lost.
