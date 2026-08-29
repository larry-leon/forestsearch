# REPORT — `vi.grf.min = -0.2` vs `NULL`: equivalence of results, and whether the ordering buys any efficiency

**Task:** `dev/tasks/cc_task_vi_grf_smoke_2026-08-30.md` (commit `01969fa2`)
**Repo / branch:** `~/Documents/GitHub/forestsearch` (pop-os), `feature/glm-extension`
**Follows:** `REPORT_vi_grf_default_2026-08-30.md` (stopped at its §2 gate, outcome (B), no edit). This is the evidence that gate asked for.
**Category:** no `R/` change, no default change, no edit to any driver, application or document; nothing already committed re-run or re-verified; synthetic data only. Artefacts under `dev/glm-continuous-sims/`: `vi_grf_smoke_2026-08-30.R` (the harness), `vi_grf_smoke_analysis_2026-08-30.R` (the membership second pass), their `.log` and `.rds`, this report. No renders, no push, no install.

---

## 1. Provenance (§1 — GATE passed)

```
pop-os · /home/larryleon/Documents/GitHub/forestsearch · feature/glm-extension
HEAD at start 695ecd39 · git status --porcelain: empty
git log -6: 695ecd39 81c9568a d5ef8c7c 424d27ba 8d6cc64c 06728783
packageVersion forestsearch 0.2.7 (installed 0.2.7)
vi.grf.min default in force: R/forestsearch_main.R:1216 `vi.grf.min = -0.2`; installed formal -0.2  (unchanged — the default task made no edit)
```

`ls ~/Downloads/cc_task_vi_grf_smoke_2026*`: one match, `cc_task_vi_grf_smoke_20260830.md` (hyphens stripped) → `dev/tasks/cc_task_vi_grf_smoke_2026-08-30.md`, committed **`01969fa2`**.

---

## 2. The comparison primitives used

Per paired run (same data, same `seedit`, `vi.grf.min = -0.2` vs `NULL`):

| what | how (from the fitted object) |
|---|---|
| declaration | `!is.null(fit$sg.harm)` |
| selected subgroup | `sort(which(fit$grp.consistency$sg.harm.id == 1))` — the subject index set, `identical()` |
| definition string | `fit$sg.harm` collapsed — used only to split *identical* from *clause-order only* |
| effect | `grp.consistency$out_sg$result$hr[1]`, `Pcons[1]` |
| candidate family | `find.grps$out.found$hr.subgroups` as a **set of rows keyed by the selected cut columns** — the `q<i>.<0/1>` column names are stable identities across the two settings (the same cut label gets the same `q` index in both runs; only the column *order* moves), so `paste(sort(qc[row == 1]))` is an exact membership key; compared with `identical()` after sorting |
| evaluated cut labels | `fit$confounders.evaluated`, compared as a set |
| consistency stage | `grp.consistency$n_candidates_total`, `n_candidates_evaluated`, `n_passed`, `early_stop_triggered`, `early_stop_candidate` |
| search-side fits | `find.grps$filter_counts$n_passed_sample_size` (combinations reaching a model fit) and `n_passed_hr` (clearing the effect screen) |
| wall-clock | `proc.time()` around the `forestsearch()` call, both settings; plus a **direct** timing of the same grf forest (`fit_causal_forest_glm()` / `causal_survival_forest()`) on the same cut matrix rebuilt from `confounders.evaluated` |

Classification: **identical** (membership and string equal) · **clause-order only** (membership equal, string differs) · **differs substantively** (membership or declaration differs).

---

## 3. Fixtures (n = 500, 20 seeds each; generators in the script)

| # | fixture | generator | targets | expectation |
|---|---|---|---|---|
| F1 | continuous, 8 covariates (6 continuous + 2 binary), real subgroup `age <= 50 & bm_hi == 1`, MD 1.5 | `gen_F1` | baseline | identical / clause-order |
| F2 | F1 + binary `rare` with 496/500 ones (prevalence 0.992; adding `{rare}` first removes 4 ≤ `rmin = 5` subjects) and a 0.4 interaction on it | `gen_F2` | `extract_idx_flagredundancy()` order dependence | **may differ — sensitive case** |
| F3 | F1 + `x_dup` (= `bm_hi`) + `x_near` (`bm_hi` with one subject flipped) | `gen_F3` | tie-breaking, statistics-keyed dedup | may differ |
| F4 | F1 + 2 noise covariates (10 total), **`max_n_confounders = 5`** | `gen_F4` | cap applied only inside the VI block | **expected to differ — positive control** |
| F5 | F1's covariates, outcome pure noise (`0.2·treat` + N(0,1)) | `gen_F5` | `vi_max > 0` / all-equal importances | identical |
| F6 | `.make_survival_data()` with `event = 0` in the treatment arm | `gen_F6` | survival zero-events guard skips VI | identical |
| F7 | `.make_survival_data()` (HR 2 in `age <= 50 & stage == 3`) | `gen_F7` | `causal_survival_forest` path | identical / clause-order |
| F8 | `.make_binary_data()` (OR 2.5) | `gen_F8` | the binary VI path | identical / clause-order |

Common arguments: `n.min = 30`, `maxk = 2`, `fs.splits = 30`, `pconsistency.threshold = 0.5` (main arm), `sg_focus = "maxeffCons"`, `use_lasso = use_grf = use_dina = FALSE`, `is.RCT = TRUE`, sequential plan, `stop_threshold` at its default (= `pconsistency.threshold`). Continuous: `effect.threshold = 0.5`, `consistency.threshold = 0`; survival: `hr.threshold = 1.25`, `hr.consistency = 1.0`; binary (OR): `hr.threshold = 1.10`, `hr.consistency = 1.0`. Extra arms: F1 under `sg_focus = "eff"` (early stopping cannot fire — reset to `NULL`, `forestsearch_main.R` L1575–1609); and a **strict early-stop arm** — F1, F2, F4, F7 at `pconsistency.threshold = 0.9` — because at 0.5 the first candidate clears the gate almost always and `n_candidates_evaluated` is 1 on both sides, which cannot discriminate.

---

## 4. The efficiency question

### 4.1 Chat's re-sort claim, verified from source

The VI order enters the search only as the **column order of `Z`**:

```r
# R/forestsearch_main.R
2820:        conf.screen <- conf.screen[seq_len(min(length(conf.screen), max_n_confounders))]
2844:  df.confounders <- df[, conf.screen, drop = FALSE]      # VI-ordered (or the FSdata order when vi.grf.min = NULL)
2845:  df.confounders <- dummy(df.confounders)
2851:  colnames(Z) <- names(df.confounders)
```

so the combination enumeration in `subgroup.search()` walks columns in VI order, which fixes (i) the order the per-candidate checks see the factors (`extract_idx_flagredundancy()` walks `seq_len(ncol(x))`, `R/subgroup_search.R` L500–514) and (ii) the row order of the raw results. Then **before** anything downstream reads them:

```r
# R/subgroup_search.R — format_search_results()
883:  hr.out <- data.table::setorder(hr.out, -HR, K)
```

and again at the entry of the consistency stage, before the evaluation loop:

```r
# R/subgroup_consistency_main.R
592:  found.hrs <- sort_subgroups_preview(found.hrs, sg_focus,
593:                                      selection_rule      = selection_rule,
594:                                      effect_neighborhood = effect_neighborhood,
595:                                      effect_log_scale    = effect_log_scale)
...
719:    for (m in seq_len(n_candidates)) {          # sequential path: walks found.hrs in ITS order
# R/subgroup_consistency_helpers.R — sort_subgroups_preview()
656:  # maxeffCons enumerates effect-descending, matching its selection key's
657:  # primary term -- which is what makes first-passer early stopping valid for
658:  # it (see subgroup.consistency()'s batch guard).
659:  if (sg_focus == "maxeffCons") {
660:    data.table::setorder(result_new, -HR, K)
661:    return(result_new)
662:  }
```

**The claim holds.** The scan order at the early-stop stage is `(-HR, K)` — fitted effect, then cut count — not VI order. VI order can survive only as the tie-break among candidates with *exactly* equal `HR` and `K` (`setorder` is stable), i.e. duplicates that the dedup step (L579) collapses anyway. The ordering therefore **cannot** reach the consistency loop except through the search-side redundancy walk (which changes *which* candidates exist, not their scan order) and the `max_n_confounders` truncation (L2820, which changes the family).

### 4.2 Measured

260 paired runs (8 fixtures × 20 seeds at `pconsistency = 0.5` under `maxeffCons`; F1 × 20 under `eff`; F1, F2, F4, F7 × 20 at `pconsistency = 0.9`), 458 s. Per-pair rows in `vi_grf_smoke_2026-08-30.rds$res`; the membership second pass in `vi_grf_smoke_analysis_2026-08-30.rds`.

**Classification (selected subgroup by subject membership):**

| arm | identical | clause-order only | **differs substantively** |
|---|---:|---:|---:|
| F1 maxeffCons 0.5 | 11 | 9 | 0 |
| F1 eff 0.5 | 10 | 10 | 0 |
| F1 maxeffCons 0.9 | 11 | 9 | 0 |
| **F2** maxeffCons 0.5 (the sensitive case) | 12 | 8 | **0** |
| F2 maxeffCons 0.9 | 12 | 8 | 0 |
| F3 maxeffCons 0.5 | 9 | 11 | 0 |
| **F4** maxeffCons 0.5 (positive control) | 8 | 6 | **6** |
| F4 maxeffCons 0.9 | 8 | 6 | **6** |
| F5 maxeffCons 0.5 | 5 | 15 | 0 |
| F6 maxeffCons 0.5 | 20 | 0 | 0 |
| F7 maxeffCons 0.5 | 10 | 10 | 0 |
| F7 maxeffCons 0.9 | 11 | 9 | 0 |
| F8 maxeffCons 0.5 | 5 | 15 | 0 |
| **all** | **132** | **116** | **12 (all F4)** |

*GATE (F4 must differ substantively): passed — 6 of 20 seeds in both F4 arms.* Wherever membership was identical, the fitted effect and `Pcons` were identical to machine precision (260/260 such pairs), and the declaration never differed. F6 (zero events in the treatment arm) produced no candidate on either side in all 20 seeds — identical trivially, and the survival guard's skip is confirmed by the timing (Δ = 0.0003 s, forest not fitted). F2 — the one demonstrated order dependence in the redundancy walk — **never changed the selected subgroup** in 40 pairs; its effect is visible one level down, in the family (next).

**The candidate family, two ways.** Keyed by *cut-set* (which cut columns a row selects), the families differ in most pairs: 0.7–10 rows per pair out of 23–625 (F2 the most, 10.2 of 624). **Every one** of those differing rows — 100 % across all 220 non-F4 pairs — is a two-cut candidate on **nested cuts of the same variable** (e.g. `{age <= 54} & {age <= 63}`), i.e. a redundant *encoding* of a membership that the family also carries as a single cut; `extract_idx_flagredundancy()` walks the pair in column order and flags it only when the wider cut comes second, so the VI order decides which encoding survives. Keyed by **subject membership**, the families are **identical in 220 of 220 non-F4 pairs** (137.8 = 137.8 memberships at F1, 600.5 = 600.5 at F2, …). The consistency-stage `n_candidates_total` therefore differs by those few duplicate rows in 30–75 % of pairs; `n_candidates_evaluated`, `early_stop_candidate` and `n_passed` are **identical in every pair** (see below); the search-side fits (`n_passed_sample_size`) are identical in every non-F4 pair.

**F4 (`max_n_confounders = 5`, 10 covariates → 25 cut columns).** At `-0.2` the cap keeps the **5 highest-importance cut columns** (seed 1: `age <= 55 | 54 | 46 | 63 | bm_hi` — four cuts of `age` and the true partner) and the family has 13 memberships against 197 at `NULL` (5–19 vs 52–379 candidates); the search fits 47 models instead of 1,171. In 12 of 40 pairs the winner changes (n\|Ĥ\| 63–125 vs 31–70; effect 0.77–1.28 vs 1.11–1.61 — the capped family loses the small, high-effect candidates), in the other 28 the true `age × bm_hi` rule survives the cap and the winner is the same. Note `fit$confounders.evaluated` is the **pre-cap** label vector (`confs_labels`) on both sides; the cap shows only in the `q` columns of `hr.subgroups` (5 vs 24–27).

**Efficiency — the counters.** `n_candidates_evaluated` and `early_stop_candidate`: paired difference **exactly 0 in all 240 pairs where early stopping applies** (every seed, every fixture, both `pconsistency` arms). At `pconsistency = 0.5` early stopping fires at candidate 1 on both sides in every seed (the `-HR` head clears 0.5), which cannot discriminate — hence the strict arm: at 0.9, F7 runs 1–12 candidates deep and stops at the **same** candidate on both sides in all 19 seeds where it stops (F1/F2/F4 still stop at 1). Under `eff` early stopping is reset and every candidate is evaluated (`n_evaluated = n_total`, 33–283), so the count differs only by the duplicate rows above. **The ordering never reaches the scan** — as §4.1 says it cannot.

**Efficiency — wall-clock** (`-0.2` minus `NULL`, mean per call, SE over 20 seeds; direct forest timing on the same cut matrix):

| arm | secs at `-0.2` | secs at `NULL` | **Δ total** (SE) | forest fit, timed directly | Δ − forest |
|---|---:|---:|---:|---:|---:|
| F1 maxeffCons 0.5 | 1.094 | 1.031 | **+0.063** (0.058) | 0.080 | −0.017 |
| F1 eff 0.5 | | | **+0.114** (0.018) | 0.072 | +0.042 |
| F1 maxeffCons 0.9 | 0.996 | 0.892 | **+0.104** (0.018) | 0.073 | +0.031 |
| F2 maxeffCons 0.5 | 1.095 | 0.956 | **+0.139** (0.023) | 0.075 | +0.064 |
| F2 maxeffCons 0.9 | 1.030 | 0.931 | **+0.099** (0.017) | 0.076 | +0.023 |
| F3 maxeffCons 0.5 | 1.217 | 1.070 | **+0.147** (0.026) | 0.075 | +0.072 |
| F4 maxeffCons 0.5 | 0.286 | 1.300 | **−1.014** (0.034) | 0.077 | −1.091 (the cap, not the ordering) |
| F4 maxeffCons 0.9 | 0.279 | 1.313 | **−1.034** (0.028) | 0.077 | −1.111 |
| F5 maxeffCons 0.5 | 1.017 | 0.875 | **+0.142** (0.017) | 0.077 | +0.065 |
| F6 maxeffCons 0.5 | | | **+0.0003** (0.002) | 0 (guard) | 0 |
| F7 maxeffCons 0.5 | 0.815 | 0.620 | **+0.195** (0.020) | 0.171 | +0.024 |
| F7 maxeffCons 0.9 | 0.827 | 0.619 | **+0.208** (0.023) | 0.163 | +0.045 |
| F8 maxeffCons 0.5 | 0.509 | 0.426 | **+0.084** (0.009) | 0.069 | +0.015 |

Outside F4 and the F6 guard, `-0.2` costs **+0.06 to +0.21 s per call** (mean +0.13 s), 3–11 SE above zero in every arm; the directly timed forest accounts for 0.07–0.17 s of it and the rest (0.01–0.07 s) is Section 5's own overhead (`apply(..., as.numeric)`, `variable_importance()`, the ordering). **There is no search-side saving to set against it**: the fit counts, the evaluated counts and the early-stop candidate are identical. Net: **positive** (the ordering does not pay for itself) in every fixture where it runs; the one negative net, F4, is the `max_n_confounders` truncation shrinking the family 8–25×, which is a different behaviour, not an ordering effect.

---

## 5. Verdict

1. **Equivalence.** 260 pairs: 132 identical, 116 clause-order only, **12 differ substantively — all 12 in F4**, none elsewhere. **F2**, the sensitive redundancy-walk case, never changed the selected subgroup (40/40 pairs membership-identical); its order dependence is real but confined to *which encoding* of a nested same-variable pair survives — 100 % of the cut-set differences across all non-F4 pairs are such pairs, and the families are identical as sets of subject memberships in 220/220. **F4** differs in 12/40 pairs by winner, size and effect, and in 40/40 by family (13 vs 197 memberships). The clause-order-only rate (36–75 % by fixture) is the prior report's observation that the VI order reaches the winner's clause order; it is cosmetic.
2. **Efficiency.** **No.** `n_candidates_evaluated` and `early_stop_candidate` differ by exactly 0 in every one of the 240 early-stopping pairs, including the strict arm where the scan runs 1–12 deep; the search-side fit count is identical in every non-F4 pair. The measured cost of the ordering is **+0.13 s per `forestsearch()` call on these fixtures (range +0.06 to +0.21 s, SE 0.01–0.06 s)**, of which the causal forest is 0.07–0.17 s. The source says why: `format_search_results()` L883 and `sort_subgroups_preview()` L660 re-sort by `(-HR, K)` before the consistency loop, so VI order cannot reach the scan.
3. **The `max_n_confounders` coupling.** From F4: with the cap set below the number of cut columns, `-0.2` keeps the top-importance cut columns (5 of 25 here) and searches a family 8–25× smaller, ~4× faster (0.28 vs 1.31 s), and in 30 % of seeds selects a different subgroup — larger and with a smaller effect — than the full family. Changing the default to `NULL` would make such a caller's cap **silently inert**: same argument, full family, different (and slower) result, with no warning, and `confounders.evaluated` would not reveal it either way (it is pre-cap). No caller in the repository sets `max_n_confounders` (prior report §2.4). That is the one behaviour that genuinely changes.

No recommendation.

---

## 6. Close-out

`git status --porcelain` before staging: the two scripts, their `.log` and `.rds`, and this report. Staged by explicit path; committed; hash in the follow-up line. **No push. No install.**

Commits: `01969fa2` task doc; the script/outputs/report commit and its hash-recording follow-up below.

---

## 7. Ten-line summary

1. Default in force confirmed `-0.2` (source L1216, installed 0.2.7); task doc `01969fa2`; 260 paired `forestsearch()` runs, synthetic data, 458 s.
2. Re-sort claim verified from source: `setorder(hr.out, -HR, K)` (`subgroup_search.R` L883) and `sort_subgroups_preview()` `setorder(-HR, K)` (`subgroup_consistency_helpers.R` L660) run before the loop at `subgroup_consistency_main.R` L719; VI order only sets `Z`'s column order (L2844–2851).
3. Selected subgroup by membership: 248/260 pairs identical or clause-order only; the 12 substantive differences are all F4 (`max_n_confounders = 5`).
4. F4 gate passed (6/20 per arm differ); F2 (redundancy walk) 0/40; F6 (survival guard) 20/20 identical, forest skipped.
5. Candidate families: identical as subject-membership sets in 220/220 non-F4 pairs; differ as cut-sets by 0.7–10 nested same-variable encodings — the redundancy walk's order dependence, made visible but membership-neutral.
6. `n_candidates_evaluated` and `early_stop_candidate`: Δ = 0 in all 240 early-stopping pairs, including the `pconsistency = 0.9` arm where F7 scans 1–12 deep.
7. Wall-clock: `-0.2` costs +0.06 to +0.21 s per call (mean +0.13 s, 3–11 SE), the forest 0.07–0.17 s of it; no search-side saving anywhere.
8. F4: cap keeps 5 of 25 cut columns, family 13 vs 197 memberships, 47 vs 1,171 fits, 0.28 vs 1.31 s; winner changes in 12/40; `confounders.evaluated` is pre-cap on both sides.
9. Under `eff` early stopping is reset and every candidate is evaluated; counts differ only by the duplicate encodings.
10. No `R/` change, no default change, nothing historical touched; no recommendation.
