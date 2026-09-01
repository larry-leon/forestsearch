# REPORT — applications re-render at 0.3.5, payload-to-payload vs the committed 0.3.1 reference — 2026-09-01

**Task:** `dev/tasks/cc_task_applications_render_0.3.5_2026-09-01.md` (committed on arrival, `84d34a08`).
**Verdict up front: the analyses did not change under the 0.3.1 → 0.3.5 revisions.** Comparison A is CLEAN — every COMPARED leaf `identical()` across all five payload pairs at full precision; Comparison B found 0 SUBSTANTIVE and 0 CLAUSE-ORDER hunks in any of the four documents. §6 **clean branch** taken: dated payload directories + this report committed. Full §6 detail in section 6.

## 1 Provenance, gates, and the §2 transplant

### 1.1 Provenance (§1 gate — PASS)

| item | value |
|---|---|
| branch | `feature/glm-extension`; `git status --porcelain` empty at start |
| pre-task HEAD | `af254e18` (bootstrap close-out report) |
| installed forestsearch | `0.3.5` (`Rscript -e 'packageVersion("forestsearch")'`) |
| efficiency-workstream close in log | `e31e612d` (assembly-skip verification report), `af254e18` (bootstrap close-out report) — both present |
| 08-31 render-task commits, identified from the log | `728ccf5a` (task doc), `d1f0e88e` (run_loo → FALSE flip), `8ad2c3f6` (payloads under `_payloads_2026-08-31`, **files >50 MB excluded from tracking: comparison_continuous.rds**), `2c2406a9` (the four 0.3.1 reference renders), `0e9d48fa` (render-and-compare report) |
| note on shas | the 08-31 report's own §8 names `713cd93e`/`e872763d` for the payloads/renders commits; those live on `backup/pre-size-fix-2026-08-31` — the branch history was amended into `8ad2c3f6`/`2c2406a9` when the >50 MB file was excluded. Content facts in that report (readback table, baseline commits) all verify against the amended commits. |
| `vi.grf.min` default in force | `NULL` (`deparse(formals(forestsearch::forestsearch)$vi.grf.min)` → `NULL`), as expected |
| toolchain | `/usr/lib/rstudio/resources/app/bin/quarto/bin/quarto` 1.9.38; R 4.6.1; host `pop-os` |
| R/ change set 0.3.1 → 0.3.5 | exactly `eb136a35` (0.3.2 signed orientation — the link this task newly pronounces on), `4dbf9f26` (0.3.3 fast paths), `3921ffdd` (0.3.4 medians-on-survivors), `f3975b99` (0.3.5 subset-vector fits), plus docs-only `ec116ffe` |

Commits this task: (1) `84d34a08` task doc alone; (2) `1ada23eb` scripts (`dev/tasks/render_driver_applications_0.3.5_2026-09-01.sh`, `dev/tasks/compare_payloads_applications_0.3.5_2026-09-01.R`), before any render; (2b) `7fef766b` comparator revision (deviation D2); (3) the final commit (section 6).

### 1.2 Document list and render commands (§2.1 — quoted from the committed 08-31 task)

From `dev/tasks/cc_task_applications_render_0.3.1_2026-08-31.md` (the transplant source; its committed report is `dev/tasks/REPORT_applications_render_0.3.1_2026-08-31.md`), "Four documents, cheap-first render order":

> 1. `quarto/applications/actg175/template_actg175_continuous.qmd`
> 2. `quarto/applications/actg175/analysis_actg175_continuous_compare_all.qmd`
> 3. `quarto/applications/gbsg/analysis_gbsg_survival_frozen_family.qmd`  *(after the flag edit below)*
> 4. `quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd`

Render invocation (08-31 task §3): `Use the RStudio-bundled binary: QUARTO=/usr/lib/rstudio/resources/app/bin/quarto/bin/quarto` … `/usr/bin/time -f "elapsed %E  maxRSS %MkB" $QUARTO render <path-to-qmd> --to html`, from the repo root. This task ran the same four, same order, same invocation, via the committed driver script; sentinel `/tmp/fs_render_sentinel_20260901`. The OC applied document has its own record and stayed out of scope.

### 1.3 Flag state at HEAD (§2.2 — PASS, nothing drifted)

`grep -n -E 'run_(loo|cv|cv_kfold|mr)\s*<-'` across the four `.qmd` at HEAD; every assignment hit, verbatim:

```
analysis_gbsg_survival_frozen_family.qmd:90:run_loo    <- FALSE     # leave-one-out (N-fold) honest out-of-sample subgroup estimate
analysis_gbsg_survival_frozen_family.qmd:98:run_cv_kfold <- FALSE    # repeated 10-fold cross-validation (stability diagnostics)
analysis_gbsg_survival_multimethod.qmd:144:run_cv  <- FALSE    # K-Fold cross-validation (kfold-cv + diagnostics)
analysis_gbsg_survival_multimethod.qmd:145:run_loo <- FALSE    # Leave-one-out / N-fold OOB (oob-cv)
analysis_gbsg_survival_multimethod.qmd:180:run_mr <- TRUE     # add mr_inference = TRUE to the forestsearch call
```

Template and compare_all: no such assignments (matching 08-31). frozen_family's `run_loo` is still the committed FALSE from `d1f0e88e`; `run_mr` TRUE is exempt (MR is core inference). No flag flipped in this task.

### 1.4 Restore list (§2.3) — the tracked paths renders overwrite

From `git ls-files` plus the documents' output declarations (08-31 task §3):

1. `quarto/applications/actg175/template_actg175_continuous.html`
2. `quarto/applications/actg175/analysis_actg175_continuous_compare_all.html`
3. `quarto/applications/gbsg/analysis_gbsg_survival_frozen_family.html`
4. `quarto/applications/gbsg/analysis_gbsg_survival_multimethod.html`
5. `quarto/applications/actg175/_payloads/template_actg175_continuous/template_actg175_continuous_payload.rds` (tracked 0.2.0 vintage)

All other render outputs land at untracked paths (`gbsg/_payloads/…`, `actg175/_payloads/analysis_actg175_continuous_compare_all/…` — neither existed before the renders).

### 1.5 Reference-set census (§2.4 — GATE PASS: both directories present, vintage 0.3.1)

| reference file (in-repo `_payloads_2026-08-31/`) | bytes | `forestsearch_version` | `built_at` |
|---|---|---|---|
| `gbsg/…/analysis_gbsg_survival_multimethod_payload.rds` | 1,901 | 0.3.1 | 2026-08-30 12:28:01 |
| `gbsg/…/analysis_gbsg_survival_frozen_family_payload.rds` | 3,113 | 0.3.1 | 2026-08-30 12:15:57 |
| `actg175/…/template_actg175_continuous_payload.rds` | 282,593 | 0.3.1 | 2026-08-30 12:06:05 |
| `actg175/…/selected_subgroups_continuous.rds` | 648 | (no version field, as documented 08-31) | — |

Every `built_at` matches the 08-31 report §7 readback table exactly. `comparison_continuous.rds` is **absent from the in-repo reference** (excluded by `8ad2c3f6` as >50 MB); the 08-31 original survives read-only at `~/fs_parked_2026-08-31/_payloads_2026-08-31/analysis_actg175_continuous_compare_all/comparison_continuous.rds` (110,436,997 bytes) and served as that one file's reference (deviation D1).

### 1.6 Triage definitions (§2.5 — quoted verbatim from the 08-31 task §5, reused in section 4)

> - **VOLATILE** — dates, wall-clock/timing tables, `built_at`, version stamps (0.2.0 → 0.3.1 is expected everywhere versions print), session/provenance echoes, worker counts, paths.
> - **CLAUSE-ORDER** — a subgroup definition string (`sg.harm` / `sg_harm` and kin) whose clause **set** is unchanged but whose order differs. Quote old and new verbatim.
> - **SUBSTANTIVE** — anything else: any estimate, CI, p-value, count, table value, membership size, or a definition whose clause set changed. Quote verbatim.

(Here the expected version stamp change is 0.3.1 → 0.3.5.)

## 2 Render walls beside the 08-31 figures

All four renders exit 0; every `.html` and payload `-newer` the sentinel; `git status --porcelain quarto/GuoHe` empty afterwards and `find . -name 'cv_out_*.rds' -newer <sentinel>` empty — no LOO/CV ran, GuoHe untouched.

| # | doc | 0.3.1 (08-31) | 0.3.5 (this task) | Δ wall | maxRSS |
|---|---|---|---|---|---|
| 1 | template_actg175_continuous | 0:11.74 | **0:10.50** | −11% | 782 MB |
| 2 | analysis_actg175_continuous_compare_all | 6:13.48 | **5:13.60** | −16% | 1862 MB |
| 3 | analysis_gbsg_survival_frozen_family | 3:09.55 | **3:09.58** | ±0% | 774 MB |
| 4 | analysis_gbsg_survival_multimethod | 12:01.57 | **10:43.80** | −11% | 918 MB |
| | **total** | **≈21m36s** | **≈19m17s** | **−10.7%** | |

**Dilution statement** — the first full-document measurement of the efficiency workstream's 0.3.2→0.3.5 changes. At the search core the speedup is large and visible in the documents' own echoes: compare_all's per-combo `time_search` dropped 0.124–0.133 min → 0.019–0.027 min (**≈5–6× faster**), and the template's `forestsearch overall` echo dropped 2.12 s → 0.662 s. But the search is a small slice of each document, so at document level the realized effect is **−10.7% total wall**:

- **compare_all** (−16%): eight combo fs walls 15.3–19.7 s → 8.8–15.5 s, but the MR anchor (61.9 → 52.2 s), GRF candidate generation, plotting and pandoc are untouched.
- **frozen_family** (±0%): bootstrap-dominated — FB 147.1 → 150.3 s, GH 30.4 → 28.2 s, while the sped-up fit itself is only 2.62 → 2.48 s of a 190 s document.
- **multimethod** (−11%): the gain is almost entirely the FS full bootstrap 300.3 → 217.9 s (the 0.3.3→0.3.5 replicate speedup at B=1000); the DINA (154.2 → 159.7 s) and GRF (193.1 → 198.0 s) bootstraps are unchanged, as expected — those methods' costs were not touched by the workstream.

## 3 Comparison A — payloads, field by field, full precision

Driver: `dev/tasks/compare_payloads_applications_0.3.5_2026-09-01.R` (committed `1ada23eb`, revised `7fef766b` — deviation D2). Method: recursive descent pruned by `identical()` (a subtree that passes `identical()` is wholly compared-equal; descent happens only where it fails, so every localized difference is exact); volatile fields enumerated by name rule with both values printed whether or not they differ; environments enumerated as non-comparable-by-type; closures/language objects compared by `deparse()` after `removeSource()` (captured environments are session artifacts).

### 3.1 Verdict per file — all five pairs CLEAN

| pair | EXCLUDED-VOLATILE | STRUCTURAL-VERSION | NONCOMP (env) | NORMALIZED-EQUAL | **COMPARED differences** |
|---|---|---|---|---|---|
| gbsg multimethod payload | 12 | 0 | 0 | 0 | **0 — identical** |
| gbsg frozen_family payload | 9 | 0 | 0 | 0 | **0 — identical** |
| actg175 template payload | 2 | 0 | 0 | 12 | **0 — identical** |
| compare_all selected_subgroups | 0 | 0 | 0 | 0 | **0 — bit-identical file-level** (`identical()` on the whole object) |
| compare_all comparison_continuous | 32 | 0 | 179 | 936 | **0 — identical** |

**STRUCTURAL-VERSION: none anywhere.** No field exists in one vintage only — the 0.3.1 and 0.3.5 payloads have identical shapes (no early-stop bookkeeping or other result-object evolution reaches these payloads; consistent with `args_call` keeping 85 formals both sides).

The NONCOMP and NORMALIZED-EQUAL entries are **confined entirely to plot objects** (`/plots/*`, `/plot_grid`, `/plot_combined_subsets/*` in comparison_continuous; the 12 template entries are gt formatter closures under `/extras/fs_tables/sg10_out/_formats`): ggplot2 4.x S7 objects were decomposed via attributes; their data, mappings (deparse-equal quosures) and every atomic leaf compared identical; the 179 environments are ggproto components (layers/scales/coords/facet), enumerated, not content-compared — the display layer is independently checked by Comparison B, where the rendered figures' diagram payloads are byte-identical.

### 3.2 The named selection lines (identical, both vintages, full precision)

- **gbsg frozen_family `sg_harm`:** `c(q6.1 = "{pgr <= 32.5}", q8.1 = "{er <= 0}")` — identical.
- **gbsg multimethod labels:** `fs = "{er <= 0} & {size <= 35}"`, `dina = "{grade3 >= 1} & {pgr <= 10}"`, `grf = "{er <= 0}"` — identical.
- **actg175 template `sg_harm`:** `c(q12.0 = "!{wtkg <= 82}", q18.0 = "!{cd40 <= 421}")` — identical.
- **actg175 compare_all, all eight combos** — identical including names, order, N_H, effect and Pcons at 17 significant digits:

```
effMaxSG / pareto        => !{cd40 <= 396} & !{wtkg <= 78}       (N_H=126, effect=51.465069356872597, Pcons=0.91)
effMaxSG / both          => !{cd40 <= 396} & !{wtkg <= 82.3284}  (N_H=87,  effect=82.782539682539706, Pcons=0.97)
effMinSG / pareto        => !{cd40 <= 507} & {age <= 37}         (N_H=66,  effect=87.9166666666667,   Pcons=0.95)
effMinSG / both          => !{cd40 <= 507} & {age <= 37}         (N_H=66,  effect=87.9166666666667,   Pcons=0.95)
eff / neighborhood       => !{cd40 <= 396} & !{wtkg <= 82.3284}  (N_H=87,  effect=82.782539682539706, Pcons=0.97)
maxSG / neighborhood     => !{cd40 <= 396} & !{wtkg <= 78}       (N_H=126, effect=51.465069356872597, Pcons=0.91)
minSG / neighborhood     => !{cd40 <= 507} & !{wtkg <= 73}       (N_H=61,  effect=82.268817204301101, Pcons=0.93)
maxeffCons / neighborhood=> !{cd40 <= 507} & {age <= 37}         (N_H=66,  effect=87.9166666666667,   Pcons=0.95)
```

Notably these are **payload-level `identical()`** results — no clause-order instability appears at the payload layer this time (both renders are 0.3.x; the 08-31 CLAUSE-ORDER hunks were 0.2.0-vs-0.3.1 artifacts).

### 3.3 The exclusion list, exhaustively (55 fields; every value printed)

**gbsg multimethod (12):** `/built_at` 2026-08-30 12:28:01 → 2026-09-01 14:37:09; `/forestsearch_version` 0.3.1 → 0.3.5; `/meta/machine/timestamp` 2026-08-30 19:28:00 UTC → 2026-09-01 21:37:09 UTC; and the nine wall-clock leaves under `/meta/timings/{fs,dina,grf}/`:

| leaf | A (0.3.1) | B (0.3.5) |
|---|---|---|
| fs/boot_sec | 300.284 | 217.921 |
| fs/mr_sec | 5.062 | 4.381 |
| fs/speedup | 59.3212169103122 | 49.7422962793883 |
| dina/boot_sec | 154.184 | 159.734 |
| dina/mr_sec | 0.357 | 0.323 |
| dina/speedup | 431.887955182039 | 494.53250773997 |
| grf/boot_sec | 193.102 | 197.966 |
| grf/mr_sec | 2.316 | 2.095 |
| grf/speedup | 83.3773747841115 | 94.4945107398555 |

The **result** fields inside `meta/timings` — `label`, `n_family`, `selection_bias` per method — were NOT excluded and compared identical (n_family 1744/84/858; selection_bias 0.556/0.357/0.609 both sides, inside the pruned-identical subtrees).

**gbsg frozen_family (9):** `/built_at` 2026-08-30 12:15:57 → 2026-09-01 14:26:36; `/forestsearch_version` 0.3.1 → 0.3.5; `/meta/timings/{loo,kfold}` NA → NA (both sides, enumerated for completeness); `/meta/timings/fit` 2.622 → 2.479; `/meta/timings/fb` 147.142 → 150.268; `/meta/timings/mr` 0.613 → 0.598; `/meta/timings/gh` 30.356 → 28.191; `/extras/mr/timing_seconds` 0.613 → 0.598.

**actg175 template (2):** `/built_at` 2026-08-30 12:06:05 → 2026-09-01 14:18:13; `/forestsearch_version` 0.3.1 → 0.3.5.

**compare_all selected_subgroups (0):** the whole file passed `identical()` — no exclusions needed (it carries no version/timing field, as documented 08-31).

**compare_all comparison_continuous (32):** per combo (×8): `/fs/<combo>/find.grps/time_search` (0.1233–0.1328 → 0.0190–0.0268 min — the raw search speedup) and `/fs/<combo>/minutes_all` (0.2552–0.3287 → 0.1462–0.2587); `/console/<combo>` (×8) — character vectors whose only differing elements are, in every combo, the two lines `Subgroup search completed in …` and `Seconds and minutes forestsearch overall = …`; `/diagnostics/<combo>/full` (×8) — text blobs whose only differing lines are the same two timing echoes (4 diff-lines each; `preview` and `summary` components identical). Every differing console/diagnostics line was volatile-pattern; the full line quotes are in the comparator output. **Nothing else in the 110 MB object differs**: `combos`, `ci_tab` (all eight tables), `fs` estimates/CIs/definitions/membership frames (`df.est`, `df.predict`, `df.test`), `grf_res`, `dina_res`, `mr_inference`, `args_call_all`, `threshold_config`, `errors` — all inside pruned-identical subtrees.

No exclusion is broader than a single named timing/version/timestamp leaf or a timing-line-only text diff; every excluded value is printed above or in the enumerated comparator output.

## 4 Comparison B — HTML against the tracked 0.3.1 baseline

Baselines: all four tracked `.html` at HEAD = the 0.3.1 reference renders (`2c2406a9`), snapshotted via `git show HEAD:<path>` before rendering. Method per 08-31 (report §5 preamble): strip all `data:` URIs (`sed -E 's/data:[a-zA-Z0-9/+-]+[^"]*//g'`), `diff -u`, classify every hunk with the §1.6 taxonomy.

| doc | hunks | VOLATILE | CLAUSE-ORDER | SUBSTANTIVE |
|---|---|---|---|---|
| template_actg175_continuous | 2 | 2 | 0 | **0** |
| analysis_actg175_continuous_compare_all | 19 | 19 | 0 | **0** |
| analysis_gbsg_survival_frozen_family | 4 | 4 | 0 | **0** |
| analysis_gbsg_survival_multimethod | 49 | 49 | 0 | **0** |

VOLATILE families (examples; every changed line audited):

- Timing echoes and timing tables: `Subgroup search completed in 0.03 → 0.01 minutes` (template), `Seconds and minutes forestsearch overall = 2.12 0.0353 → 0.662 0.011` (template), per-combo `18.2–19.7 → 13.6–15.5` ×8 plus `maxeffCons MR anchor elapsed: 61.9 → 52.2 s` and timing-table cells `275.4/4.6 → 227.6/3.8`, `362.4/6.0 → 303.0/5.0` (compare_all), `MR gate: 0.61 → 0.6 s`, timing-table `2.622→2.479 / 147.142→150.268 / 0.613→0.598 / 30.356→28.191` (frozen_family), `703.9 → 625.9 s` total and per-method bootstrap rows (multimethod).
- Version stamps read at render time: provenance caption `(forestsearch 0.3.1, read at render time) → (forestsearch 0.3.5, …)` (frozen_family ×1, multimethod gt subtitles ×4 — note "16 of 85 formals" / "7·10·11 of 36 knobs" counts unchanged, only the version string moved); sessionInfo `forestsearch_0.3.1 → forestsearch_0.3.5`.
- Render date `August 30 → September 1, 2026` + `dcterms.date` (multimethod); gt CSS id churn (`csbmbkduaa → jlzikwaowh` and kin); htmlwidget ids regenerated with **byte-identical digraph payloads** (both DiagrammeR trees: `age <= 50 / age <= 43 / er <= 0` and the single-split `er <= 0` — content unchanged).

**Every `Subgroup identified:` / `SELECTED` / estimate / CI / count line in all four documents is unchanged** (they appear only as context lines in the diffs). No SUBSTANTIVE hunk exists to quote; no CLAUSE-ORDER hunk exists — unlike 08-31's 12, consistent with both renders being 0.3.x. No cross-BLAS artifacts this time: same host, same BLAS, so the 08-31 report's pre-attributed print deltas (items 4–6) do not recur — those hunks are absent entirely.

## 5 Comparison outputs and scratch

Scratch kept outside the repository (deviation D4): baseline snapshots, the four diffs, the render log, and the full comparator output (`compare_A2.md`, 1,578 lines, containing the complete NONCOMP/NORMALIZED-EQUAL enumerations and every quoted volatile line) live in the session scratchpad, not the worktree. Everything deliverable from them is quoted in sections 2–4.

## 6 §6 branch taken and the final gate

**Clean branch** (no COMPARED differences, no SUBSTANTIVE hunks): both `_payloads_2026-09-01/` directories (payloads + `html/` with the four new renders) and this report committed as the final commit. The 0.3.5 set now sits beside the 08-31 set — nothing replaced, nothing deleted.

One tracking exception, matching the `8ad2c3f6` precedent exactly (deviation D3): the fresh `comparison_continuous.rds` (111,266,019 bytes — over GitHub's 100 MB hard limit) is **staged on disk in the dated directory but excluded from the commit**, the same exclusion the reference set itself carries.

Final-gate outputs (run at close; also verified immediately after the restore step):

- Restore verification: after `git checkout --` of the five §1.4 paths, `git diff --stat` (whole tree) was **empty** — every restored path byte-identical to HEAD; the default-path template payload read back `0.2.0 | 2026-08-17 09:17:37` (the tracked vintage, intact).
- `git diff --stat HEAD -- '*_payloads*'` — empty apart from the new dated directories: **no existing payload directory modified**. No `rm` of any `.rds` occurred anywhere; `_payloads_2026-08-31/`, the parked set, and the tracked template payload are untouched.
- `git status --porcelain` at close: only `?? …/_payloads_2026-09-01/analysis_actg175_continuous_compare_all/comparison_continuous.rds` (the deliberately-uncommitted >100 MB file inside this task's own declared directory) — quoted in the close-out below.
- Payload identity: every staged versioned payload reads back `forestsearch_version == 0.3.5` (built_at 14:18:13 / 14:26:36 / 14:37:09, matching the render minutes). `git check-ignore -v` on both dated dirs: no match (exit 1) — plain `git add` used.

## 7 Ten-line verdict

1. **The applications' analyses did not change under the 0.3.1 → 0.3.5 revisions** — the sentence this task was commissioned to produce, now measured rather than inferred.
2. Every COMPARED leaf across all five payload pairs is `identical()` at full precision: estimates, CIs, tables, subgroup definitions, memberships, meta, extras, MR outputs, consistency values.
3. That includes the 110 MB comparison_continuous object: all eight combos' fs results bit-identical; only 32 timing/version leaves and timing-echo lines differ.
4. The harm-subgroup selections are unchanged by name and value: gbsg `{pgr <= 32.5} & {er <= 0}` (frozen) and `{er <= 0} & {size <= 35}` / `{grade3 >= 1} & {pgr <= 10}` / `{er <= 0}` (multimethod FS/DINA/GRF); actg175 `!{wtkg <= 82} & !{cd40 <= 421}` (template) and all eight compare_all selections.
5. This closes the one previously unmeasured link: 0.3.1 → 0.3.2 (`eb136a35`, signed orientation) is now confirmed result-preserving on both applications, extending the guarded 0.3.2→0.3.5 chain back to 0.3.1.
6. No STRUCTURAL-VERSION drift: payload shapes are identical between vintages — no fields appeared or vanished.
7. The HTML layer agrees: 74 hunks across four documents, all VOLATILE (timings, version stamps, dates, id churn); zero SUBSTANTIVE, zero CLAUSE-ORDER.
8. Render effect of the efficiency workstream on full documents: −10.7% total (21m36s → 19m17s); the ≈5–6× search-core speedup is diluted by untouched DINA/GRF/MR/bootstrap costs, exactly as the task predicted.
9. The 0.3.5 payload set is committed beside the 08-31 reference as the current-codebase reference (comparison_continuous.rds on disk but untracked, per the reference set's own >50 MB precedent).
10. Deviations below; nothing was deleted, no `R/` file touched, no push.

## 8 Deviations and notes

- **D1 — comparison_continuous reference from the parked copy.** The in-repo `_payloads_2026-08-31/` excludes it (`8ad2c3f6`); the 08-31 original at `~/fs_parked_2026-08-31/…/comparison_continuous.rds` (110,436,997 B, parked 2026-08-31 per the standing record, read-only here) served as its reference. All other pairs used in-repo references.
- **D2 — comparator revised after its first pass** (`7fef766b`). First pass flagged 19 "differences"; inspection showed all three groups were comparator shortcomings, none a result change: (a) `timing_seconds` missing from the volatile name set (a wall-clock field, 0.613→0.598 s); (b) console character-vectors bypassed the volatile-line triage that blob leaves got (all 16 differing elements were the two timing-echo lines per combo); (c) ggplot2 4.x S7 objects (`typeof "object"`) fell to a type-mismatch fallback instead of being decomposed. The revision fixed all three; the re-run's output is the deliverable. Both script versions are in history.
- **D3 — the >100 MB payload staged but untracked** (section 6), mirroring the reference set's own exclusion. If Larry wants it tracked, LFS or an amended rule is his call.
- **D4 — scratch outside the repo.** The task's Writes list declares no scratch dirs, so baselines/diffs/comparator output went to the session scratchpad instead of `dev/tasks/_*` (unlike 08-31); everything deliverable is quoted herein.
- **D4b — transient `Rplots.pdf`.** The comparator's first pass printed ggplot objects while formatting values, which opened a null graphics device and dropped an untracked `Rplots.pdf` (52,857 B, mtime 14:38) at the repo root; it was this task's own byproduct and was removed at close. Nothing pre-existing was touched.
- **D5 — template payload grew 282,593 → 1,121,809 bytes with identical compared content.** The only non-identical content in that pair is the 12 gt formatter closures (deparse-equal, environments differ) plus built_at/version; the size delta lives in those closures' serialized environments — a serialization artifact, not a result change.
- **D6 — plot internals.** 179 ggproto environments (all inside plot objects) are enumerated as non-comparable-by-type rather than content-compared; their data/mappings and the rendered display layer (Comparison B, byte-identical diagram payloads) are verified. No taxonomy entry was silently dropped.
