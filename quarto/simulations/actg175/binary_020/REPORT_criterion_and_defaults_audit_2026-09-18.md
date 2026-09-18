# REPORT — The admission criterion and the defaults, at HEAD: read-only audit

**Task:** `dev/tasks/TASK_audit_criterion_and_defaults_2026-09-18.md` (committed alone at `bb3492d2`).
**Date:** 2026-09-18. **Branch:** `feature/glm-extension`. **Machine:** `pop-os`.

**HEAD at the read:** `bb3492d2f4aa778fed5db04ca68f1218e7c479c4` (the task-document commit; its parent
`82b421e9` is the partial pass's commit, and `R/` is unchanged between them — `R/` last moved at
`064fce91`, 2026-09-16).

**Tree state:** **dirty, untracked only.** `git status --porcelain` carries **no tracked modification**;
26 untracked entries, all campaign artefacts from the crashed OR campaign and the continuous logs
(`binary_020/*.html`, `binary_020/logs_or/`, ten `binary_020/mr_or_harm/*_d5000/` bundle directories,
`continuous/logs_{mddina,mdgrf,mdsgnb20}/`). Nothing untracked is under `R/`, so every quotation below
is from committed source.

**Category:** read-only. No `R/` edit, no fit, no compute, no install, no push. The only writes are the
task document and this report.

**Placement.** §5 says "beside the repository's existing `REPORT_*` files". There is no repository-root
`REPORT_*`; they live in campaign directories. This one is filed beside
`REPORT_binary_admission_check_2026-09-18.md`, the report it extends.

**Relation to the prior record.** `dev/tasks/TASK_binary_admission_check_2026-09-18.md` and
`quarto/simulations/actg175/binary_020/REPORT_binary_admission_check_2026-09-18.md` are read and cited.
Where that report already established a fact it is cited, not re-derived; every extension is marked.

---

## 1. The reference table

Every row is from current source. `FS` = `subgroup_method = "consistency"` (the ForestSearch path),
`GRF` = `subgroup_method = "grf"`, `DINA` = `subgroup_method = "dina"`. "Front door" is
`forestsearch()`; "helper" is the exported function that actually enforces the floor.

### 1a. Admission floors

| # | Argument | Aliases | Front-door default (`forestsearch_main.R`) | Helper default | Unit | Per arm? | Counts | Outcome types | Identifiers | Resolved at | Enforced at |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A1 | `n.min` | none | `60` (`:1246`) | `subgroup.search()` **`30`** (`subgroup_search.R:83`) | count of subjects | pooled | observations | all four | FS (search); GRF and DINA via their own `n_min`, see A2/A3 | `forestsearch_main.R:1368-1389` (adaptive branch only, when `n.min = NULL`) | `subgroup_search.R:613` (GLM), `:660` (survival) — `if (nx <= n.min) return(status 4L)` |
| A2 | `n.min` → GRF `n_min` | GRF spells it `n.min` in the config, `n_min` in the enumerator | inherits A1 | `grf_subgroup_labels.R:259`, `:285` default **`60L`** | count of subjects | pooled | observations | all four | GRF only | `grf_helpers.R:25-30` (`create_grf_config()`) | `grf_subgroup_labels.R:269`, `:302` — `if (nS < n_min || nS > n - 1L) next`; also `grf_helpers.R:178` `size >= n.min` on the policy-tree path |
| A3 | `n.min` → DINA `n_min` | `dina_args$n_min`; `dina_subgroup(n_min=)`; `n_min.frac` | inherits A1 (`forestsearch_helpers.R:1129`, `:1490`) | `dina_subgroup.R:320-321` **`n_min = 60L`, `n_min.frac = 0.10`** | count of subjects | pooled | observations | all four | DINA only | `dina_subgroup.R:416-421`; `.resolve_dina_args()` `forestsearch_helpers.R:1125-1132` | `dina_subgroup.R:746` — `if (n_S < n_min) next` |
| A4 | `n.min.frac` | none | `0.10` (`:1247`) | `dina_subgroup.R:321` `0.10` | proportion of N | — | observations | all four | FS (and DINA's own copy) | `forestsearch_main.R:1372-1377` — `n.min <- max(60L, ceiling(n.min.frac * N))`, **only when `n.min = NULL`** | not enforced directly; feeds A1 |
| A5 | `d0.min` | none | `10` (`:1266`) | `subgroup.search()` **`15`** (`subgroup_search.R:83`) | count | **yes — control arm** | **events** (binary `Y=1`; survival `dd`) | survival, binary. **Skipped for continuous and count** | **FS only** | not resolved; passed through `search_overrides` (`forestsearch_main.R:3002`) | `subgroup_search.R:598` (binary), `:654` via `meets_event_criteria()` (`:738-740`, survival) |
| A6 | `d1.min` | none | `10` (`:1267`) | `subgroup.search()` **`15`** (`subgroup_search.R:83`) | count | **yes — treatment arm** | as A5 | as A5 | **FS only** | as A5 (`:3003`) | as A5 |
| A7 | `d0.min`/`d1.min` fallback | — | as A5/A6 | — | count | yes | **subjects**, not events | "unknown GLM type" only (`outcome_type` neither `binary` nor `continuous`/`count`) | FS | — | `subgroup_search.R:603-607` |
| A8 | `minp` | none | `0.025` (`:1269`) | `subgroup.search()` **`0.05`** (`subgroup_search.R:85`) | proportion | pooled | column mean of each cut indicator, whole sample | all four | FS | forced to `0` under `sg_focus = "maxeff"` (`forestsearch_main.R:1549`) | `subgroup_search.R:572` via `meets_prevalence_threshold()` (`:719-722`) — `all(colMeans(x) >= minp)` |
| A9 | `rmin` | none | **not a `forestsearch()` formal** | `subgroup.search()` `5` (`subgroup_search.R:86`) | count of subjects | pooled | observations | all four | FS | set to `0` under `maxeff` (`forestsearch_main.R:3021`) | `subgroup_search.R:577` via `extract_idx_flagredundancy()` (`:500-513`) — each added factor must shrink membership by **more than** `rmin` |
| A10 | GLM fit floor | none | — | hard-coded `6L` | count of subjects | pooled | observations | binary, continuous, count | FS | — | `subgroup_search.R:872` (`fit_glm_for_subgroup`) and `:917` (`..._fast`) — `if (nrow(df_sg) < 6L) return(NULL)` → status 5 |
| A11 | split-half size floor | none | — | hard-coded `5` | count of subjects | pooled, **within each half** | observations | all four | FS consistency stage | — | `subgroup_consistency_helpers.R:271` (and the twostage copy at `:1696`) — `if (nrow(split1) < 5 || nrow(split2) < 5) return(NA_real_)` |
| A12 | split-half per-arm floor (GLM) | none | — | hard-coded `3L` | count | **yes, per arm, per half** | **subjects** | binary, continuous, count | FS consistency stage | — | `subgroup_consistency_helpers.R:282-283` (twostage copy `:1704-1705`) |
| A13 | split-half event floor (survival) | none | — | hard-coded `2` | count | pooled within half | **events** | survival | FS consistency stage | — | `subgroup_consistency_helpers.R:300` (twostage copy `:1718`) |
| A14 | `m1.threshold` | none | `Inf` (`:1256`) | `subgroup.consistency()` `Inf` (`subgroup_consistency_main.R:357`) | count (treatment-arm `m1`) | treatment arm | median/`m1` column | all four (`m1` is `NA` on GLM paths) | FS | — | `subgroup_consistency_main.R:530-531` (+ the NA drop at `:531`), pre-consistency |
| A15 | `max_subgroups_search` | `stop_Kgroups` (the `subgroup.consistency()` spelling) | `Inf` (`:1276`) | `subgroup.consistency()` **`stop_Kgroups = 200`** (`subgroup_consistency_main.R:368`) | count of candidates | — | candidates | all four | FS | `Inf` under `maxeff` (`forestsearch_main.R:1548`) | `subgroup_consistency_main.R` §truncation (pool ordered, then truncated) |
| A16 | GRF whole-sample arm floor | none | inherits A1 | `grf_helpers.R:671` | count | **yes, per arm, whole data** | observations | all four | GRF | — | `grf_helpers.R:688-691` — `if (n_treat < n.min || n_control < n.min) stop()` |
| A17 | `max.minutes` | none | `3` (`:1268`) | `subgroup.search()` `30` (`subgroup_search.R:84`) | minutes | — | — | all four | FS | — | **nowhere — inert.** `subgroup.search()` accepts and forwards it and never compares it; `fs_family_report.R:339-340` records this in source |

### 1b. Threshold and estimand defaults

| # | Argument | Aliases | Front-door default | Helper default | Scale | Outcome types / estimand | Identifiers | Resolved at | Compared at |
|---|---|---|---|---|---|---|---|---|---|
| B1 | `outcome_type` | none | `c("survival","binary","continuous","count")` → `"survival"` (`:1283-1284`) | — | — | — | all | `match.arg` | — |
| B2 | `effect_measure` | none | `NULL` (`:1285`) | — | — | **`survival`→`"HR"`; `binary`→`"RD"`; `continuous`→`"MD"`; `count`→`"IRR"`** | all | **two sites**: `forestsearch_main.R:1334-1341` and `:1731-1738` | — |
| B3 | `adverse_outcome` | none | `NULL` (`:1287`) | `TRUE` in the estimator closures | — | resolves `TRUE` for `binary` and `count`, `FALSE` for `continuous` and `survival` | all | **two sites**: `:1342-1344` and `:1717-1719` | `glm_effect_estimators.R:303-305` (flips `Y → 1-Y`), `:966` / `fit_glm_for_subgroup_fast` `subgroup_search.R:925` (flips `y → -y`) |
| B4 | `effect.threshold` (screening) | **`hr.threshold`** | `effect.threshold = NULL` (`:1248`), `hr.threshold = 1.25` (`:1250`) | `subgroup.search()` **`hr.threshold = 1.0`** (`subgroup_search.R:84`); `subgroup.consistency()` **`hr.threshold = 1.0`** (`subgroup_consistency_main.R:354`) | natural ratio for `HR`/`OR`/`RR`/`IRR`, converted to **log**; identity for `RD`/`IRD`/`MD` | survival `HR>=1.25`; OR/RR/IRR `log(1.25)`; **RD `0.05`**, **IRD `0.01`**, **MD `0.0`** | FS, DINA (as `m_diff`), GRF (via `admission$effect_floor`) | alias merge `:1352`; measure resolution `:1796-1885`; `screen_threshold` `subgroup_search.R:126-130` | `subgroup_search.R:634` (GLM) and `:675` (survival) — `if (!disable_effect_floor && hr <= hr.threshold) return(status 6L)`, i.e. **strict: a candidate must exceed it** |
| B5 | `consistency.threshold` (per split) | **`hr.consistency`** | `consistency.threshold = NULL` (`:1249`), `hr.consistency = 1.0` (`:1251`) | `subgroup.consistency()` `hr.consistency = 1.0` (`:355`), `consistency_threshold = NULL` (`:379`); `consistency_resample()` `hr.consistency = 1.0` (`consistency_resample.R:394`) | as B4 | survival `HR>1.0` (log `log(max(hr.consistency,0.001))`); OR/RR/IRR `log(1.0) = 0`; **RD/IRD/MD `0.0`** | FS only (GRF and DINA compute no Pcons —  `.fs_admission_applies()`, `forestsearch_helpers.R:2314-2315`) | alias merge `:1353`; remap `:1802-1804`, `:1827-1842`, `:1850-1856`; log conversion `:1884` | **GLM: `subgroup_consistency_helpers.R:294-295`, `>=` (non-strict). Survival: `:308`, `>` (strict).** |
| B6 | `pconsistency.threshold` (proportion of splits) | none | `0.90` (`:1257`) | `subgroup.consistency()` `0.9` (`:356`) | proportion in [0,1] | identical for every outcome type — it is a rate, not an effect | FS only | forced to `0` under `maxeff` (`forestsearch_main.R:1546`) | `subgroup_consistency_main.R` candidate qualification; two-stage screen derived at `:503-507` as `max(0.5, p* - 2.5*SE)` |
| B7 | `stop_threshold` | none | **promise** `= pconsistency.threshold` → `0.90` (`:1263`) | `subgroup.consistency()` `NULL` (`:369`) | proportion | as B6 | FS only | **reset to `NULL` for `sg_focus` in `hrMaxSG/hrMinSG/hr/maxSG/minSG`** (`:1636-1670`), silently when not user-set | `subgroup_consistency_main.R:832`, `:997` — `pcons >= stop_threshold` |
| B8 | `sg_focus` | GLM spellings `eff`/`effMaxSG`/`effMinSG` normalize to `hr`/`hrMaxSG`/`hrMinSG` | `"hr"` (`:1252`) | `subgroup.consistency()` `"hr"` (`:365`) | — | all | all three | `.normalize_sg_focus()`; `forestsearch_main.R:1536-1570` (maxeff overrides) | `.fs_admission_applies()` `forestsearch_helpers.R:2300-2318` |
| B9 | `dmin.grf` | none | `0.0` (`:1232`) | `grf_main.R:159` `0.0` | DR-score scale (natural) | all four | **GRF only** | `forestsearch_main.R:1974-1977` (`missing(dmin.grf)` → `0.0`; a no-op given the formal default is already `0.0`) | `grf_subgroup_labels.R:361` — `cand[cand$effect >= dmin, ]` (non-strict) |
| B10 | `fs.splits` | `n.splits` (the `subgroup.consistency()` spelling) | `1000` (`:1255`) | `subgroup.consistency()` **`n.splits = 100`** (`:358`) | count | all | FS | forwarded `forestsearch_main.R:3125` | — |
| B11 | `maxk` | none | `2` (`:1272`) | `subgroup.search()` `2` (`:88`); `subgroup.consistency()` **`7`** (`:362`) | count | all | all three | — | `generate_combination_indices()` |
| B12 | `use_twostage` | none | `TRUE` (`:1279`) | `subgroup.consistency()` **`FALSE`** (`:375`) | — | all | FS | `FALSE` under `maxeff` (`:1548`) | — |
| B13 | `selection_rule` / `effect_neighborhood` | — | `"neighborhood"` (`:1253`) / `0.10` (`:1254`) | same (`subgroup_consistency_main.R:366-367`) | proportion | all | all three | — | band logic in `sort_subgroups()` |

---

## 2. Part A — the admission criterion and its floors

### 2.1 The stages are separate, and only one of them has a per-arm floor

Four distinct stages gate a candidate. They are **not** the same floor applied three times.

| Stage | Where | Floors in force |
|---|---|---|
| **Candidate construction** (cut columns) | `get_fsdata()`, `forestsearch_main.R` §5 (`vi.grf.min`, `max_n_confounders`) | none of A1–A17; variable-importance filtering only, and only when `vi.grf.min` is non-`NULL` |
| **Screening / enumeration** | `evaluate_combination_with_status()`, `subgroup_search.R:553-700` | A8 `minp` → A9 `rmin` → **A5/A6 `d0.min`/`d1.min`** → **A1 `n.min`** → A10 (`< 6` rows) → B4 effect floor |
| **Pre-consistency pruning** | `subgroup_consistency_main.R:530-625` | A14 `m1.threshold`, near-duplicate removal, A15 `max_subgroups_search` |
| **Consistency** | `subgroup_consistency_helpers.R:254-311`, and a second private copy inside `evaluate_consistency_twostage()` at `:1685-1727` | A11 (`< 5` per half), A12 (`< 3` per arm per half, GLM), A13 (`< 2` events per half, survival) — **no `n.min`, no `d0.min`/`d1.min`** |

`subgroup.consistency()` has **no `n.min`, `d0.min` or `d1.min` formal at all** (`subgroup_consistency_main.R:352-385`), so those three floors are structurally absent from the consistency stage; a candidate that cleared them at screening is never re-checked.

### 2.2 Disposition of a candidate that fails — the distinction the task asks for

**Every screening-stage failure is silent.** `evaluate_combination_with_status()` returns
`list(status = <k>, result = NULL)` (`subgroup_search.R:568`, `:573`, `:579`, `:599`, `:606`, `:614`,
`:628`, `:635`, `:655`, `:661`, `:669`, `:676`, `:698`) — no condition is raised. The statuses are tallied into
`filter_counts` (`subgroup_search.R:275-300`) and **printed only when `details = TRUE`**
(`:188-219`, the filtering summary at `:206-218`). `filter_counts` is carried on the returned object via `format_search_results()`, so the
counts survive; the campaign runs with `details = FALSE`, so nothing is printed.

- A8, A9, A5/A6, A1, A10, B4 → **silently dropped**, counted.
- A11/A12/A13 (consistency stage) → the split returns `NA_real_`, which is **excluded from the Pcons
  denominator**, not counted as a failed split. A candidate whose splits are mostly `NA` therefore has
  its Pcons computed on the surviving splits.
- A14 `m1.threshold` with **all** candidates removed → `warning("All subgroups removed after filtering
  NA m1 values")` and an empty result (`subgroup_consistency_main.R:532-533`).
- Empty `hr.subgroups` → `warning("No valid hr.subgroups")` (`subgroup_consistency_main.R:514`).
- `subgroup.search()` erroring → `warning("Error in subgroup.search: ...")` (`forestsearch_main.R:3050`).

**No floor errors.** The only `stop()` on a floor is A16 (`grf_helpers.R:688-691`) and DINA's
`n_min > nrow(df)` guard (`dina_subgroup.R:420-421`), both of which are configuration errors, not
candidate rejections.

### 2.3 Per-arm behaviour

The only per-arm floors on a *candidate* are A5/A6, and they exist on the FS screening path only:

- **binary:** `d0_sg <- sum(yy[id.x == 1 & tt == 0])`, `d1_sg <- sum(yy[id.x == 1 & tt == 1])`;
  `if (d0_sg < d0.min || d1_sg < d1.min)` (`subgroup_search.R:595-600`). **Events per arm.**
- **survival:** `calculate_event_counts()` (`:727-733`) + `meets_event_criteria()` (`:738-740`),
  `d0 >= d0.min && d1 >= d1.min`. **Events per arm.**
- **continuous and count:** **skipped entirely.** `is_continuous <- outcome_type %in% c("continuous",
  "count")` (`:588-589`); the `else if (!is_continuous)` at `:603` is not taken, and `:609` records
  "Continuous: skip Status 3 entirely". **What admits a candidate on those paths instead is A1 `n.min`
  alone** (plus A8/A9/A10 and the effect floor) — there is no per-arm condition of any kind.
- **"unknown GLM type"** (A7): an `outcome_type` that is neither `binary` nor `continuous`/`count`,
  with an `estimator_fn` present, falls to a per-arm **subject** count against `d0.min`/`d1.min`
  (`:603-607`). Not reachable from `forestsearch()` (`match.arg` restricts `outcome_type` to the four),
  reachable from a direct `subgroup.search()` call.

Confirming the prior report (§3.1, §3.3): the binary path does implement at least `d0.min`/`d1.min`
events in each arm, at candidate admission, on the FS path.

### 2.4 Binary and GLM: the zero-cell question — **extension, not previously answered**

**There is no zero-cell guard anywhere on the fit path.** Neither `fit_glm_for_subgroup()`
(`subgroup_search.R:869-893`) nor `fit_glm_for_subgroup_fast()` (`:915-945`) inspects the 2×2 table;
neither does the binary closure (`glm_effect_estimators.R:297-337`). Reading the code path for a
candidate with an all-event arm:

- **`OR`** (`.estimate_or`, `glm_effect_estimators.R:463-527`): plain `stats::glm(..., binomial)`. Under
  separation `glm()` **does not error** — it returns a finite but divergent coefficient (of order ±25)
  with a correspondingly huge SE, and `fit$converged` may be either value. The returned
  `list(estimate, se, converged = fit$converged, ...)` is finite, so `fit_glm_for_subgroup()`'s only
  test, `is.na(res$estimate)` (`subgroup_search.R:879`, fast path `:929`), passes. **`converged` is discarded** —
  `fit_glm_for_subgroup()` never reads it (grep for `$converged` returns no hit in
  `subgroup_search.R`). The candidate reaches status 6 with a huge `hr`, so the effect floor
  **admits it**, and it enters the consistency stage. *This is the "retained with a degenerate
  estimate" disposition, and it is the code path behind the recorded Stage-2 halt.*
- **`RD`** (the **default** for binary, see §3.1) (`.estimate_rd`, `:346-455`): a three-tier fallback.
  Tier 1 identity-link binomial `stop()`s on non-convergence, tier 2 logistic G-computation likewise,
  and **tier 3 is unconditional raw proportions** (`:439-453`) returning `converged = FALSE` and a
  finite `p1 - p0`. So an RD candidate with a zero cell **never returns `NA`**; it returns a boundary
  estimate with `converged = FALSE`, which is again discarded.
- **`RR`** (`.estimate_rr`, `:534-600`) and the count estimators (`:615-800`) return `converged` the
  same way, with the same consumer.

`converged` **is** consulted elsewhere — `consistency_resample.R:279`, `:305`; `frontier_cis.R:241`,
`:249`, `:279`; `dina.R:848` — so the search path's omission is a local one, not a package-wide
convention.

This extends the prior report's §3.5 finding 3 ("an events-only floor cannot prevent separation") from
the floor side to the fit side: **nothing downstream of the floor catches it either.**

### 2.5 Identifier forwarding, per identifier

- **FS (`consistency`)**: receives A1, A5, A6, A8, A9, A17 via `search_overrides`
  (`forestsearch_main.R:2998-3016`) merged over `args_call_all` and filtered to
  `names(formals(subgroup.search))` (`:3042-3044`). All floors in force.
- **GRF**: `.forestsearch_grf_select()` is called with `n.min = n.min` and `dmin.grf = dmin.grf`
  (`forestsearch_main.R:2426-2427`) — **and nothing else from the floor set.** `d0.min` / `d1.min` do
  not appear anywhere in `grf_main.R`, `grf_helpers.R`, `grf_args.R` or `grf_subgroup_labels.R`. GRF's
  own floors are A2 (`nS < n_min`) and A16 (whole-sample per-arm size). **No per-arm event floor.**
- **DINA**: `.forestsearch_dina_select()` receives `n.min = n.min` (`forestsearch_main.R:2243`) and
  forwards it as `n_min` to both `dina_frontier()` and `dina_subgroup()`
  (`forestsearch_helpers.R:1462`, `:1490`). `d0.min` / `d1.min` do not appear in `dina*.R`. DINA's
  floors are A3 (`n_S < n_min`) and the effect floor `m_diff`. **No per-arm event floor.**

This confirms the prior report's §3.5 qualification 1, and extends it: the *size* floor **is**
forwarded to both (they inherit `forestsearch()`'s `n.min`, they do not fall back to their own `60L`
defaults), so what is missing on GRF and DINA is specifically the **per-arm** floor, not all floors.

### 2.6 Where a floor is bypassed

1. **Continuous and count bypass A5/A6 entirely** (`subgroup_search.R:588-609`). Size only.
2. **GRF and DINA bypass A5/A6 entirely** (§2.5). Size and effect only.
3. **`sg_focus = "maxeff"` disables the effect floor and `minp`**: `disable_effect_floor` is set from
   `!.admit_applies[["effect"]]` (`forestsearch_main.R:3015`), `rmin → 0` (`:3021`), `minp → 0`
   (`:1550`), `pconsistency.threshold → 0` (`:1546`), `max_subgroups_search → Inf` (`:1549`),
   `use_twostage → FALSE` (`:1548`). A warning is raised listing every override (`:1552-1562`), so this
   bypass is audible.
4. **The MR candidate family re-enumerates with a *different* size test and **no** per-arm floor**:
   `forestsearch_main.R:3372` uses `if (length(mem) >= n.min)` where the search uses
   `if (nx <= n.min) reject`, i.e. **MR admits a subgroup of exactly `n.min` that the search rejects**
   — an off-by-one between the two enumerations. The source comment at `:3353-3357` already records
   that `d0.min`/`d1.min` and `max_subgroups_search` are not replayed; the `>=` vs `>` discrepancy is
   **not** recorded there and is reported here as new.
5. **A17 `max.minutes` is inert on every path** (`fs_family_report.R:339-340`). Not a bypass of a
   floor so much as a floor that was never wired.

---

## 3. Part B — the defaults and their resolution

### 3.1 `effect_measure`: the resolution-site count, and the binary default

**Two sites, confirmed by search, not assumed.** `grep -n 'effect_measure' R/forestsearch_main.R`
returns exactly two assignment blocks:

- `R/forestsearch_main.R:1334-1341`:
  ```r
  if (outcome_type != "survival" && is.null(effect_measure)) {
    effect_measure <- switch(outcome_type,
      binary     = "RD",
      continuous = "MD",
      count      = "IRR"
    )
  }
  ```
- `R/forestsearch_main.R:1731-1738`: the identical `switch`, inside the
  `if (outcome_type != "survival")` block that builds the estimator closure.

The second is unreachable in practice — the first runs unconditionally earlier and leaves
`effect_measure` non-`NULL` — but both exist and would have to be changed together. No third site
exists anywhere in `R/` (the remaining hits are uses, `@param` text, or `if (outcome_type ==
"survival") "HR" else effect_measure` expressions).

**The answer the pending change rests on: for `outcome_type = "binary"` with `effect_measure` unset,
HEAD resolves `effect_measure = "RD"` — not `"OR"`.** `R/forestsearch_main.R:1336`.
Survival is handled separately and never carries an `effect_measure`; it is labelled `"HR"` at the
consumption sites (`:2003`, `:2039-2040`, `:2320`).

**`adverse_outcome` and `effect_measure` do *not* live in the same block.** They are adjacent at the
first site (`effect_measure` `:1334-1341`, then `adverse_outcome` `:1342-1344`) but **separated by
~14 lines and in the opposite order** at the second (`adverse_outcome` `:1717-1719`, then
`effect_measure` `:1731-1738`, with `.validate_outcome_threshold_config()` between them at
`:1700-1711`). So there are **four** duplicated resolution sites across two arguments, in two
different orders.

### 3.2 Duplicated defaults across layers (§1a of the task)

Every divergence found, front door vs. helper. The task's claim about `n.min` is confirmed; it is the
smallest of the set.

| Argument | `forestsearch()` | `subgroup.search()` | `subgroup.consistency()` | Reachable? |
|---|---|---|---|---|
| `n.min` | **60** `:1246` | **30** `subgroup_search.R:83` | (absent) | no — always forwarded |
| `d0.min` | **10** `:1266` | **15** `subgroup_search.R:83` | (absent) | no — always forwarded |
| `d1.min` | **10** `:1267` | **15** `subgroup_search.R:83` | (absent) | no — always forwarded |
| `minp` | **0.025** `:1269` | **0.05** `subgroup_search.R:85` | (absent) | no |
| `max.minutes` | **3** `:1268` | **30** `subgroup_search.R:84` | moot — inert | — |
| `hr.threshold` | **1.25** `:1250` | **1.0** `subgroup_search.R:84` | **1.0** `:354` | no |
| `maxk` | **2** `:1272` | 2 `subgroup_search.R:88` | **7** `:362` | no |
| `fs.splits` / `n.splits` | **1000** `:1255` | — | **100** `:358` | no |
| `max_subgroups_search` / `stop_Kgroups` | **Inf** `:1276` | — | **200** `:368` | no |
| `use_twostage` | **TRUE** `:1279` | — | **FALSE** `:375` | no |
| `stop_threshold` | **`pconsistency.threshold`** (0.90) `:1263` | — | **NULL** `:369` | no |
| `n.min` (DINA) | 60 | — | `dina_subgroup.R:320` **60L** | only on a direct `dina_subgroup()` call |
| `n.min` (GRF) | 60 | — | `grf_subgroup_labels.R:259`,`:285` **60L** | only on a direct call |
| `hr.consistency` | 1.0 `:1251` | — | 1.0 `:355` | — (agree) |
| `pconsistency.threshold` | 0.90 `:1257` | — | 0.9 `:356` | — (agree) |
| `m1.threshold` | Inf `:1256` | — | Inf `:357` | — (agree) |
| `sg_focus` | "hr" `:1252` | — | "hr" `:365` | — (agree) |

**The load-bearing question — which value a call actually receives.** `forestsearch()` builds
`search_overrides` with `n.min = n.min`, `d0.min = d0.min`, `d1.min = d1.min`, `minp = minp`,
`max.minutes = max.minutes`, `maxk = maxk` **explicitly** (`forestsearch_main.R:2998-3011`), then `modifyList(args_call_all, search_overrides)` (`:3042`) — and `args_call_all` is
the full formals capture (`:1401`), so even without the overrides every shared name is supplied. The
same holds for `subgroup.consistency()`: `filter_call_args(args_call_all, subgroup.consistency,
consistency_overrides)` (`:3180-3184`). **No divergent helper default is reachable through
`forestsearch()`.**

**Committed callers that reach a helper directly:** `dev/glm/glm_test_suite.qmd:430-436` and
`:465-470`, and their `_extended` twin — all four pass `n.min = 30, d0.min = 10, d1.min = 10`
explicitly. `tests/` contains no direct `subgroup.search()` call. **No committed caller reaches the
helper's default.** The divergence is latent, not active.

### 3.3 The three threshold arguments, per outcome type and estimand

Resolution is at `R/forestsearch_main.R:1796-1885`, one site, three branches.

| Outcome type | Resolved `effect_measure` | Screening default | Consistency default | `pconsistency` | Comparison scale |
|---|---|---|---|---|---|
| `survival` | `HR` | `1.25` → `log(1.25) = 0.2231` in `threshold_config` (`:2003`); **but the search compares on the natural scale**, see below | `1.0` → `log(max(1.0, 0.001)) = 0` | 0.90 | log for `threshold_config`, **natural** for the search |
| `binary`, `effect_measure = "RD"` (default) | `RD` | **`0.05`** (`:1800`) | **`0.0`** (`:1803`) | 0.90 | identity |
| `binary`, `effect_measure = "OR"` | `OR` | `log(1.25) = 0.2231` (`:1884`) | `log(1.0) = 0` (`:1885`) | 0.90 | log |
| `binary`, `effect_measure = "RR"` | `RR` | `log(1.25)` | `0` | 0.90 | log |
| `continuous` | `MD` | **`0.0`** (`:1851`) | **`0.0`** (`:1855`) | 0.90 | identity |
| `count` (default) | `IRR` | `log(1.25)` | `0` | 0.90 | log |
| `count`, `effect_measure = "IRD"` | `IRD` | **`0.01`** (`:1800`) | **`0.0`** (`:1803`) | 0.90 | identity |

**A scale asymmetry worth stating plainly.** On the survival path the screening comparison is on the
**natural HR scale** — `subgroup.search()` receives `hr.threshold = 1.25` (no `effect_threshold` is
set, `forestsearch_main.R:3006`), `screen_threshold` falls through to it (`subgroup_search.R:126-130`),
and `fit_cox_for_subgroup()` returns `exp(beta)` as `$hr`. On every GLM path the comparison is on the
**link scale** — `effect_threshold` is `log(...)` for ratio measures, and `fit_glm_for_subgroup()`
returns the raw coefficient. `threshold_config$screening` is `log(hr.threshold)` for survival, which
is the MR/admission scale, **not** the search's scale. The two are consistent because each consumer
reads the right one, but "the screening threshold" names two different numbers depending on which
object you read it from.

### 3.4 The identity-scale remapping, with the actual mapped values

`forestsearch_main.R:1796-1856`, verbatim values:

- `RD`: `1.25 → 0.05`, `1.0 → 0.0` (`:1800`, `:1803`, and the ratio-detection remap at `:1823-1824`).
- `IRD`: `1.25 → 0.01`, `1.0 → 0.0` (same lines).
- `MD`: `1.25 → 0.0`, `1.0 → 0.0` (`:1851`, `:1855`). **`MD` has no bounds check** — the
  `> 1.0` ratio-detection warning at `:1809` and `:1828` applies to `RD`/`IRD` only (`:1796`,
  and the comment at `:1807-1808` states why).
- Ratio measures (`OR`/`RR`/`IRR`): **no remap**; `effect_threshold <- log(effect_threshold)` and
  `consistency_threshold <- log(consistency_threshold)` (`:1884-1885`), preceded by an explicit
  `stop()` on a non-positive threshold (`:1865-1883`).

The task calls the survival-default remap "silent". **Partly:** the `!user_set_*` branches at
`:1799-1804` and `:1850-1856` raise nothing. The `user_set_threshold && effect_threshold > 1.0` branch
at `:1809-1825` **warns**, and remaps to the same default. `.validate_outcome_threshold_config()`
(`:79-184`) adds four further warnings — a binary-looking outcome under `survival`, a >10-level outcome
under `binary`, a ratio-scale threshold below 0.5, and an `RD`/`IRD` threshold above 1.0.

### 3.5 Where each threshold is compared, and whether the inequality is strict

| Threshold | Site | Inequality | Strict? |
|---|---|---|---|
| screening, GLM | `subgroup_search.R:634` | `if (!disable_effect_floor && glm_result$hr <= hr.threshold) reject` | **strict** — must exceed |
| screening, survival | `subgroup_search.R:675` | `if (!disable_effect_floor && cox_result$hr <= hr.threshold) reject` | **strict** |
| screening, re-filter in consistency | `subgroup_consistency_main.R:545`, `:548` — `hr.subgroups$HR >= hr.threshold` (threshold passed as `effect_threshold` on the GLM path, `forestsearch_main.R:3155`) | `>=` | **non-strict** |
| consistency stage entry | `forestsearch_main.R:3097` | `any(hr_values > check_threshold)` | **strict** |
| per-split, GLM | `subgroup_consistency_helpers.R:294-295` | `res1$estimate >= c && res2$estimate >= c` | **non-strict** |
| per-split, survival | `subgroup_consistency_helpers.R:308` | `hr.split1 > c && hr.split2 > c` | **strict** |
| Pcons vs `pconsistency.threshold` | `subgroup_consistency_main.R` qualification | `>=` | non-strict |
| Pcons vs `stop_threshold` | `subgroup_consistency_main.R:832`, `:997` | `pcons >= stop_threshold` | non-strict |
| GRF frontier | `grf_subgroup_labels.R:361` | `cand$effect >= dmin` | non-strict |
| DINA collector | `dina_subgroup.R:749` (and `:888` on the depth-2 path) | `if (mean_tau < m_diff) next` | non-strict |
| size floor | `subgroup_search.R:613`, `:660` | `if (nx <= n.min) reject` | **strict** — must exceed |
| MR family size |  `forestsearch_main.R:3372` | `if (length(mem) >= n.min)` | **non-strict** |

**Three strict/non-strict inconsistencies for the same conceptual floor:** screening (`<=` reject)
vs the consistency-stage re-filter (`>=` admit); per-split survival (`>`) vs per-split GLM (`>=`);
search size floor (`>`) vs MR family size floor (`>=`).

### 3.6 The consistency-stage entry condition, and what happens when it does not run

```r
  hr_values <- find.grps$out.found$hr.subgroups$HR
  check_threshold <- if (!is.null(consistency_threshold)) consistency_threshold else hr.consistency
  has_subgroups <- any(hr_values > check_threshold, na.rm = TRUE)
```
`forestsearch_main.R:3090-3097`, guarded by `find.grps` being non-`NULL`, not a `try-error`, and
carrying `out.found$hr.subgroups` (`:3083-3086`).

Note the threshold used here is the **consistency** threshold, not the screening one, applied to the
**screened** candidates' effects — so it is a second, weaker effect floor at the stage boundary. For
the GLM default (`RD`, consistency `0.0`) it reduces to "some candidate has `RD > 0`"; for survival
(`hr.consistency = 1.0`) to "some candidate has `HR > 1.0`" — while the screening floor that produced
those candidates was already `1.25`, so on survival with defaults the condition cannot fail if the
family is non-empty.

**When it does not run** (`:3104`, the whole block is `if (has_subgroups)`): `sg.harm`,
`grp.consistency`, `df.est_out` stay at their §7 initialisation of `NULL` (`:3065-3069`). **Downstream
consequences, not just the branch:**
- `forestsearch_bootstrap_dofuture()` and `forestsearch_Kfold()` require `sg.harm`; CV's required
  components are `c("df.est", "args_call_all", "sg.harm")` (`forestsearch_cross_validation.R:306`) and
  `.fs_cv_base_frame()` falls back to an ITT frame (`:82-99`, `:313`).
- MR is skipped with the explicit reason "no subgroup was identified, so there is nothing to de-bias"
  (`forestsearch_main.R:3330-3331`), and `mr_harm_confirmed` is `NA`.
- No condition is raised by the `has_subgroups == FALSE` path itself.

### 3.7 What `fpr_calibration()` and the OC functions require

- **`fpr_calibration()`** (`R/fpr_calibration.R`) takes the two thresholds as **`c1`** and **`c2`**
  (`:49-50`) and writes them into the `forestsearch()` call as the **legacy aliases**:
  `modifyList(fs_params, list(hr.threshold = c1, hr.consistency = c2, ...))` (`:249-255`), then
  `do.call(forestsearch, fs_call)` (`:297`). So it compares on **whatever scale `forestsearch()`
  resolves them to** — it passes natural-scale values and lets §3.3 do the conversion. Its own
  analytic layer, `compute_detection_probability_glm(theta, d_eff, c1, c2, effect_scale = "ratio")`
  (`:240-244`), is called with **`effect_scale = "ratio"` hard-coded**, regardless of `outcome_type` or
  the resolved `effect_measure`. It additionally requires `n_min` (default `60L`, `:188`) and
  documents that it "should match `fs_params$n.min`" (`:76-77`) — an unenforced coupling.
- **`fs_oc_family()`** (`R/fs_oc_family.R`) re-implements the size floor as `Pg >= n.min / n`
  (`:59-61`, `:385`) — **`>=`, where the search uses `>`** — and re-implements the `n.min = NULL`
  adaptive resolution itself (`:239-251`), duplicating `forestsearch_main.R:1368-1389`. It reads
  `rmin` from `subgroup.search()`'s default because it is not a `forestsearch()` formal (`:232`).
- **`fs_oc_predict()` / `fs_oc_grid()`** read `pconsistency` from
  `forestsearch_args$pconsistency.threshold`, falling back to
  `eval(formals(forestsearch)$pconsistency.threshold)` (`fs_oc_predict.R:173-174`,
  `fs_oc_grid.R:479-480`) — i.e. they read the front-door default reflectively rather than restating
  it. That is the one place in the package that cannot drift.

### 3.8 Identifier invariance: do `consistency`, `dina` and `grf` use the same estimand and threshold?

**No.** For a given `outcome_type` the three identifiers read the threshold from three different
places and, under the default `effect_measure`, compare three different quantities.

| Identifier | Reads the screening threshold from | Estimand compared | Under `binary` + default `RD` | Under `binary` + `OR` |
|---|---|---|---|---|
| `consistency` | `effect_threshold` (resolved, `forestsearch_main.R:3006`, `:3027`) | the resolved `effect_measure`, link scale | `RD >= 0.05`, identity | `log-OR > log(1.25)` |
| `dina` | **raw `hr.threshold`** (`forestsearch_main.R:2243`), converted at `forestsearch_helpers.R:1436` as `m_diff <- if (family == "gaussian") hr.threshold else log(hr.threshold)` | DINA's `tau_hat` on the **family link** (log-OR for `binomial`) | **`log-OR >= log(1.25)`** — *not* `RD`, and *not* `0.05` | `log-OR >= log(1.25)` — agrees |
| `grf` | `admission_resolved$effect_floor` = `threshold_config$screening` (`forestsearch_main.R:2026-2030`, `forestsearch_helpers.R:1664-1670`) under `grf_select_statistic = "effect"` (the default); **`dmin.grf` (0.0)** under `"dr"` | the resolved `effect_measure`, natural scale | `RD >= 0.05` — agrees | `OR >= 1.25` — agrees |

So under a **ratio** `effect_measure` the three agree. Under the **default binary estimand (`RD`)**,
DINA silently screens on a **log-odds-ratio** floor derived from `hr.threshold`, while FS and GRF
screen on a **risk-difference** floor. DINA's family comes from `.map_dina_family(outcome_type)`
(`forestsearch_helpers.R:1011`), which keys on `outcome_type` and **never consults `effect_measure`**.
Neither DINA nor GRF applies any consistency threshold — `.fs_admission_applies()` returns
`c(effect = TRUE, consistency = FALSE)` for both (`forestsearch_helpers.R:2314-2315`), and neither
computes a Pcons.

### 3.9 Is "not supplied" detectable? — **the most consequential finding in Part B**

The mechanism is `forestsearch_main.R:1350-1351`:

```r
  user_set_threshold   <- !is.null(effect.threshold)      || !missing(hr.threshold)
  user_set_consistency <- !is.null(consistency.threshold) || !missing(hr.consistency)
```

- For the **new** spellings (`effect.threshold`, `consistency.threshold`) detection is by `is.null()`
  on a `NULL` default — **robust**, survives any wrapper.
- For the **legacy** spellings (`hr.threshold`, `hr.consistency`) detection is by `missing()` —
  **lost through any wrapper that passes the argument explicitly**, because `missing()` reports how the
  call was written, not what value arrived.

**This is not hypothetical: the package's own bootstrap and CV are such wrappers.**
`args_call_all <- mget(names(formals()), ...)` (`:1400-1401`) captures **every** formal including
`hr.threshold = 1.25` and `hr.consistency = 1.0`; `bootstrap_analysis_dofuture.R:406` takes that list
as `args_FS_template`, `:558` copies it to `args_FS_boot`, and `:614` runs
`do.call(forestsearch, args_FS_boot)`. Inside the replicate both legacy arguments **are** supplied, so
`missing()` is `FALSE` and `user_set_* == TRUE`, even though the original analysis set neither.
`forestsearch_cross_validation.R:345` reads the same list.

Tracing the consequence for an analysis that left both thresholds at their defaults:

| `effect_measure` | Original call | Inside a bootstrap/CV replicate | Same? |
|---|---|---|---|
| `OR` / `RR` / `IRR` / survival | `log(1.25)`, `log(1.0)` | `log(1.25)`, `log(1.0)` — the ratio branch (`:1863-1885`) ignores `user_set_*` | **yes** |
| `RD` | screening `0.05`; consistency **`0.0`** | screening `0.05` (via the `user_set && > 1.0` warn-and-remap at `:1809-1825`, so **with a warning per replicate**); consistency **stays `1.0`** — `:1802` needs `!user_set_consistency`, and `:1827` needs `!user_set_consistency`, and `:1830`'s `> 1.0` is false at exactly `1.0` | **no** |
| `IRD` | screening `0.01`; consistency `0.0` | screening `0.01` (with warning); consistency **`1.0`** | **no** |
| `MD` | screening **`0.0`**; consistency **`0.0`** | screening **`1.25`** (`:1850` needs `!user_set_threshold`; **no** fallback branch for `MD`, since `:1809` is `RD`/`IRD`-only); consistency **`1.0`** | **no** |

On the identity-scale measures a replicate therefore runs with a consistency threshold of `1.0` — an
`RD` or `MD` of 1.0 per split — which no candidate can meet, and on `MD` a screening threshold of
`1.25` as well. `dmin.grf` avoids this because its resolved value is written back with
`.sync_args_call_all(..., "dmin.grf")` (`:1976-1977`); `effect_threshold` and `consistency_threshold`
are **locals, not formals**, so they cannot be synced, and `hr.threshold` / `hr.consistency` are never
overwritten with their resolved equivalents.

**Scope.** This does not touch the committed OR campaign, which passes `effect.threshold` and
`consistency.threshold` explicitly (§4), nor any survival or ratio-measure work. It bites an
identity-scale analysis (`RD`, `IRD`, `MD`) that leaves the thresholds at their defaults and then runs
bootstrap or CV. Stated as a finding; no fix is proposed here.

### 3.10 `git log -S` on the consistency default and on any `0.8` derivation

- `git log --all -S'hr.consistency = 1.0' -- R/forestsearch_main.R` → **one commit**, `f8ca1eae`
  ("CV update"). The default has been the literal `1.0` since.
- `git log --all -S'hr.consistency <-' -- R/` → `d5ef8c7c` (a local in `fs_family_report()`),
  `bc4690b3`, `0c36bb1e`. Both of the latter are the **alias merge**
  (`if (!is.null(consistency.threshold)) hr.consistency <- consistency.threshold`); `0c36bb1e`
  additionally wrote `args_call_all$hr.consistency <- hr.consistency`, which `bc4690b3` removed when it
  moved the alias merge ahead of the `args_call_all` capture.
- `git log --all -S'consistency_threshold <- 0' -- R/` → `af970d3a` (GLM phase 1), `4cff6762`,
  `0c36bb1e` — all three are the `→ 0.0` identity remap, never a derivation from another value.
- `git log --all -G'consistency.*0\.8|0\.8.*consistency' -- R/forestsearch_main.R` → `eb199177`,
  `b8df97dc`, `f8ca1eae`; inspecting each diff, **every `0.8` hit is the roxygen example value
  `pconsistency.threshold = 0.85`**. The only other `0.8` in current `forestsearch_main.R` are
  `frac.tau = 0.8` (`:1233`, a GRF quantile) and `floor(0.80 * n_phys)` (`:64`, the worker count).

**Answer: no derivation of the consistency threshold from any other quantity has ever existed in this
package's history, on any scale.** It has always been a literal default, remapped per estimand.

---

## 4. Part C — what the campaign template actually passes

**File:** `quarto/simulations/actg175/binary_020/sim_fs_mr_field_or_template.qmd` (committed at
`1f99dfe4`; matches HEAD). The sole `forestsearch()` call is `base_args` + `method_args` at `:650-684`,
invoked at `:693`.

Verbatim, with line numbers:

```
:277   or_threshold   <- 0.90   # effect.threshold      ([S] :192; OR scale, > 1 = harm)
:278   or_consistency <- 0.80   # consistency.threshold ([S] :193)
:279   pconsistency   <- 0.90   # ([S] :194)
:280   fs_splits <- 500L; maxk <- 2L; n_min <- 60L; d0_min <- 10L; d1_min <- 10L   # [S] :195
:290   dmin.grf             <- 0.0       # [S] :206, the study's literal (DR-score harm floor)
:291   dina_args            <- list()    # [S] :207
:306   outcome_type   <- "binary"    # [S] :174
:307   effect_measure <- "OR"        # [S] :175
```

```
:657     outcome_type    = outcome_type, effect_measure = effect_measure,
:661     effect.threshold       = or_threshold,
:662     consistency.threshold  = or_consistency,
:663     pconsistency.threshold = pconsistency, n.min = n_min,
:664     d0.min = d0_min, d1.min = d1_min, maxk = maxk,
```

| Item | Supplied or inherited | Value |
|---|---|---|
| `effect_measure` | **supplied explicitly** (`:657`, from `:307`) | `"OR"` |
| screening threshold | **supplied**, new spelling `effect.threshold` (`:661`, from `:277`) | `0.90` on the OR scale → `log(0.90) = -0.1054` |
| consistency threshold | **supplied**, new spelling `consistency.threshold` (`:662`, from `:278`) | `0.80` on the OR scale → `log(0.80) = -0.2231` |
| consistency proportion | **supplied** `pconsistency.threshold` (`:663`, from `:279`) | `0.90` |
| `n.min` | **supplied** (`:663`, from `:280`) | `60L` |
| `d0.min` | **supplied** (`:664`, from `:280`) | `10L` |
| `d1.min` | **supplied** (`:664`, from `:280`) | `10L` |
| `n.min.frac` | **inherited** — not passed | `0.10`, inert (`n.min` is supplied, so the adaptive branch at `:1368` is not taken) |
| `minp` | **inherited** | `0.025` |
| `rmin` | not a formal | `5` (`subgroup.search()` default) |
| `m1.threshold` | **inherited** | `Inf` — inert |
| `max.minutes` | **inherited** | `3` — inert everywhere (§2.6.5) |
| `max_subgroups_search` | **supplied** (`:668`, from `:283`) | `Inf` |
| `stop_threshold` | **supplied** (`:666`, from `:276`) | `NULL` — "pinned NULL = scan the full family" |
| `sg_focus` | **supplied** (`:659`, from `:179`) | `"effMaxSG"` (env-overridable, default `effMaxSG`) → normalizes to `hrMaxSG` |
| `adverse_outcome` | **supplied** (`:658`, from `:293`) | `TRUE` |
| `dmin.grf` | supplied **only on the `grf` branch** (`:683`, from `:290`) | `0.0` |
| `dina_args` | supplied **only on the `dina` branch** (`:681`, from `:291`) | `list()` — so DINA's frontier `n_min` inherits `n.min = 60` |

**Reported without judgment, as §4 asks:** the expectation on record — that the template passes `OR`
explicitly — **holds**; `:307` and `:657`. Three further facts about the file, stated because they are
what it says:

1. The screening threshold is **`0.90` on the OR scale**, i.e. **below the null**, and the consistency
   threshold is **`0.80`**, further below. Both are annotated "OR scale, > 1 = harm" (`:277`, `:354`).
   Under §3.5 the search admits any candidate with `log-OR > log(0.90)`.
2. `d0.min = 10` / `d1.min = 10` **are** passed on all three branches (they are in `base_args`), but
   per §2.5 only the `consistency` branch consults them. On the `grf` and `dina` branches they are
   accepted and ignored.
3. Because the template passes `effect.threshold` and `consistency.threshold` (the `NULL`-defaulted
   spellings), `user_set_*` is `TRUE` by `is.null()` and the §3.9 wrapper defect **cannot** affect
   this campaign.

---

## 5. Discrepancy list — current source vs. what the task document asserted

The task's §1b claims, each answered from source.

| # | Claim as written | Verdict at HEAD | Correction |
|---|---|---|---|
| 1 | per-arm floor for binary counts events per arm and rejects below `d0.min`/`d1.min` — `subgroup_search.R:596-598` | **confirmed**, line numbers off by a little | The block is `subgroup_search.R:591-600`; the test is at `:598`. |
| 2 | per-arm floor for survival does the same via `calculate_event_counts()`/`meets_event_criteria()` — `:653-654`, `:727`, `:738-739` | **confirmed** | Call at `:652-655`; `calculate_event_counts()` at `:727-733`; `meets_event_criteria()` at `:738-740`. |
| 3 | **the per-arm floor is skipped entirely for continuous and count** — `:609` | **confirmed** | `:609` is the comment; the governing predicate is `is_continuous <- outcome_type %in% c("continuous","count")` at `:588-589` and the `else if (!is_continuous)` at `:603`. |
| 4 | size floor rejects when `nx <= n.min` — `:613` (GLM), `:660` (survival) | **confirmed exactly** | — |
| 5 | **`n.min` defaults to 60 in `forestsearch()` but 30 in `subgroup_search()`** — `:1246`; `subgroup_search.R:83` | **confirmed exactly** | And **incomplete**: the same line, `subgroup_search.R:83`, also declares `d0.min = 15, d1.min = 15` against the front door's `10`/`10`, plus `minp` `0.05` vs `0.025` (`:85`), `max.minutes` `30` vs `3` (`:84`), `hr.threshold` `1.0` vs `1.25` (`:84`). `subgroup.consistency()` adds six more. **Sixteen divergences in total (§3.2), not one.** |
| 6 | `d0.min` and `d1.min` both default to 10 — `:1266-1267` | **confirmed exactly** | — |
| 7 | `sg_focus` defaults to `"hr"` — `:1252` | **confirmed exactly** | — |
| 8 | `effect.threshold` / `consistency.threshold` default `NULL`, resolved from `hr.threshold = 1.25` / `hr.consistency = 1.0` **at a single site** — `:1248-1251`, resolution `:1793-1794` | **confirmed with a correction to the site** | Defaults and line numbers are right. The **alias merge** is at `:1350-1353`; the **estimand resolution** is the block `:1796-1885`, not `:1793-1794`. It *is* a single site, but a 90-line one with three branches, and the resolved values are stored in **locals** (`effect_threshold`, `consistency_threshold`), not written back to the formals — which is what makes §3.9 possible. |
| 9 | `adverse_outcome` defaults `NULL` and resolves `TRUE` for binary and count only, at two sites — `:1287`, `:1342-1343` and `:1718-1719` | **confirmed**, off by one line | Default at `:1287`; sites at `:1342-1344` and `:1717-1719`. |
| 10 | ratio measures converted to log; RD, IRD and MD compared on identity, with a **silent** remap of the survival default — `:1781-1860` | **confirmed with two corrections** | The block runs `:1796-1885`, not to `:1860`. The remap is **silent only when the threshold was not user-set**; the `user_set && > 1.0` path **warns** (`:1809-1825`, `:1830-1842`). And `MD` has **no** bounds check at all (`:1807-1808` says why), so a ratio-scale value passed with `MD` is accepted unremapped. |

**Claims the task treated as open, now answered:**

| Question | Answer |
|---|---|
| Which `n.min` does a call actually receive at each entry point? | Always the front door's, on every path. `search_overrides` (`:3002-3007`) and the `args_call_all` merge (`:3041`) both supply it; `subgroup.consistency()` has no `n.min`; GRF and DINA receive it as `n_min` (`:2424`, `:2243`). |
| Does any committed caller reach the helper's default? | **No.** The only direct `subgroup.search()` callers in the repository are `dev/glm/glm_test_suite.qmd:430`/`:465` and the `_extended` twin, and all four pass `n.min`, `d0.min`, `d1.min` explicitly. `tests/` has none. |
| Exactly which outcome types bypass the per-arm floor, and what admits a candidate instead? | `continuous` and `count` (`:588-589`, `:603`, `:609`). What admits instead: `minp` (A8) → `rmin` (A9) → **`n.min` alone** (A1) → the `< 6` row floor (A10) → the effect floor. No per-arm condition of any kind. |
| How many `effect_measure` resolution sites? | **Two** (`:1334-1341`, `:1731-1738`), confirmed by search. The handoff's count of two is right. |
| What does the binary path resolve `effect_measure` to when unset? | **`"RD"`** (`forestsearch_main.R:1336`). Not `"OR"`. |
| Do the `effect_measure` and `adverse_outcome` resolutions live in the same duplicated blocks? | **No** — adjacent and in that order at the first site (`:1334-1341` then `:1342-1344`), separated and in the **reverse** order at the second (`:1717-1719` then `:1731-1738`, with a validator call between). |

---

## 6. Findings not asked about

Recorded as facts, with no task attached.

1. **`max.minutes` is inert.** `forestsearch()` defaults it to `3`, `subgroup.search()` to `30`, and no
   code path compares it. `fs_family_report.R:339-340` already states this in source, and
   `0409ab98` ("docs: flag max.minutes as inert, schedule v0.3.0 removal") records the intent. Any
   reading of a run's wall-clock behaviour that assumes a 3-minute search cap is wrong.
2. **`fs_family_report()` is an in-repo map of exactly this ground.** `R/fs_family_report.R:250-345`
   enumerates every candidate-family stage with its argument, its status
   (deterministic / disabled / inert / data-dependent) and a prose note citing file and line. It is
   generated per fitted object, so it reports the floors *that run actually used*. Both workstreams
   should be pointed at it before re-deriving anything from this table.
3. **`converged` is computed and discarded on the search path.** Every GLM estimator returns it
   (`glm_effect_estimators.R:476`, `:504`, `:516`, `:665`, `:707`, `:720`, `:753`, `:771`) and
   `fit_glm_for_subgroup()` tests only `is.na(res$estimate)` (`subgroup_search.R:879`, fast path `:929`). Six other
   consumers in the package do check it (§2.4). A non-convergent candidate is therefore indistinguishable
   from a convergent one at screening.
4. **RD never returns `NA`.** `.estimate_rd()`'s tier 3 (`glm_effect_estimators.R:439-453`) is an
   unconditional raw-proportions fallback with `converged = FALSE`. Combined with finding 3, an `RD`
   candidate with a degenerate table is admitted with a boundary estimate and no signal.
5. **`stop_threshold`'s documented default is inert under the default `sg_focus`.** The formal is the
   promise `stop_threshold = pconsistency.threshold` (`:1263`, `0.90`), and `:1636-1670` resets it to
   `NULL` for `sg_focus` in `hrMaxSG/hrMinSG/hr/maxSG/minSG` — which includes the **default** `"hr"` —
   silently when the user did not set it explicitly. It is meaningful for `maxeffCons` only.
6. **`user_explicit <- !missing(stop_threshold)`** (`:1640`) has the same wrapper fragility as §3.9,
   but is harmless: the resolved `NULL` is written back with `.sync_args_call_all()` (`:1669-1670`),
   so only the warning is suppressed in a replicate, not the behaviour.
7. **The `missing(dmin.grf)` block at `:1974-1977` is a no-op for the default.** Its comment describes
   a survival RMST-scale default that the formal no longer carries — `dmin.grf = 0.0` at `:1232`.
   `5cb003e3` ("Update forestsearch defaults: parallel_args, dmin.grf") is where it changed. The block
   still matters for a user-supplied value on a survival run.
8. **The MR family enumeration uses `>= n.min` where the search uses `> n.min`** (`:3372` vs
   `subgroup_search.R:613`/`:660`). The source comment at `:3353-3357` records the `d0.min`/`d1.min`
   and `max_subgroups_search` gaps but not this one. It makes MR's family a superset by at least the
   exactly-`n.min` subgroups.
9. **`fs_oc_family()` re-implements two things `forestsearch()` already resolves** — the adaptive
   `n.min` rule (`fs_oc_family.R:239-251` vs `forestsearch_main.R:1368-1389`) and the size floor, as
   `Pg >= n.min / n` (`:385`) where the search uses strict `>`. Its own comments cite the
   `forestsearch_main.R` section it mirrors, so the duplication is deliberate and documented.
10. **`fpr_calibration()` hard-codes `effect_scale = "ratio"`** in its analytic detection-probability
    call (`fpr_calibration.R:240-244`), independent of `outcome_type` and of the resolved
    `effect_measure`. Its simulation arm is unaffected (it delegates to `forestsearch()`), so this
    touches `P1`/`L_eff`/`fpr_corrected`, not `fpr_hat`.
11. **The consistency-stage entry condition is a second effect floor** (§3.6) and it uses the
    *consistency* threshold, not the screening one. On survival with defaults (screen `1.25`, entry
    `1.0`) it cannot fail on a non-empty family; on `RD` with defaults (screen `0.05`, entry `0.0`) it
    likewise cannot. It only binds when a user sets `consistency.threshold` above
    `effect.threshold`.
12. **The split-half floors have two implementations.** `run_single_consistency_split()`
    (`subgroup_consistency_helpers.R:254-311`) and the private
    `.run_single_consistency_split()` inside `evaluate_consistency_twostage()` (`:1685-1727`) carry
    byte-equivalent copies of A11 (`< 5` per half), A12 (`< 3` per arm per half) and A13 (`< 2` events
    per half), plus the same `>=` / `>` asymmetry between the GLM and survival comparisons
    (`:294-295` / `:308` vs `:1713` / `:1724`). The two-stage screen threshold
    `max(0.5, p* - 2.5*SE)` is likewise duplicated (`subgroup_consistency_main.R:503-507` and
    `subgroup_consistency_helpers.R:1735-1738`). Since `use_twostage = TRUE` is the front-door default
    (`forestsearch_main.R:1279`), **the copy that actually runs by default is the second one**, not
    the exported helper — which matters for anyone tracing a split-level floor from the exported
    surface.
13. **Prior report §6.7, still open and unrelated to this audit:** the Gate 2 record header gap
    (`[ -f ... ]` inside a `>> record` redirection) affects the committed `mdgrf` and `mddina` records
    too. Noted so it is not lost.

---

## 7. Scope

Read-only. `R/` unchanged. No fit, no simulation, no install, no push. The only commits are
`dev/tasks/TASK_audit_criterion_and_defaults_2026-09-18.md` (`bb3492d2`) and this report. No fix is
proposed and no follow-up task is attached; §3.9, §2.4 and §3.8 are reported as findings and stop
there.
