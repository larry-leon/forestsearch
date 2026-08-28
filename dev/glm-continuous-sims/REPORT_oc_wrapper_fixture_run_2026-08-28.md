# REPORT — `fs_oc_predict()` on the MD40 fixture: `rmin` semantics, floor counts, and the three-column comparison

**Follows:** `dev/glm-continuous-sims/REPORT_oc_wrapper_build_2026-08-28.md` (commit `9e538ec2` / `438a2e81`)
**Script and output:** `dev/glm-continuous-sims/fixture_run_fs_oc_predict_2026-08-28.R`, `.log` (committed beside this report)
**Category:** read-only. No `R/` edit, no render, no push. The fixture is the deterministic MD40 rebuild (gate against the tracked n = 500 payload's committed truth: all four `|diff|` = 0; `fs_dgm_scale` regions identical to the payload's `$scale$regions`), and `forestsearch_args` are the driver's (`confounders.name` = the 12 analysis variables, `conf.cont_jcuts = list(age = 10, preanti = 10)`, `n.min = 60`, `maxk = 2`, `sg_focus = "maxeffCons"`, `effect.threshold = 30`, `consistency.threshold = 10`), exactly as for the M = 1601 enumeration in the build report.

---

## 1. From source

### (a) What `rmin` counts, and against what — `extract_idx_flagredundancy()` (`R/subgroup_search.R` L500–514)

```r
extract_idx_flagredundancy <- function(x, rmin) {
  n <- nrow(x)
  id.x <- rep(1, n)
  flag.redundant <- FALSE
  nx.last <- n
  for (m in seq_len(ncol(x))) {
    if (!flag.redundant) {
      id.x <- id.x * x[, m]
      nx.this <- sum(id.x)
      if (nx.last - nx.this <= rmin) flag.redundant <- TRUE
      nx.last <- nx.this
    }
  }
  list(id.x = id.x, flag.redundant = flag.redundant)
}
```

`x` is the analysis-sample indicator matrix restricted to the combination's selected columns (`evaluate_combination_with_status()` L562: `x <- zz[, selected_cols, drop = FALSE]`), so `n = nrow(x)` is the **full analysis sample size** and `nx.this` is the **count of subjects** in the running intersection after factor `m` is added. The test is `nx.last − nx.this <= rmin`: the number of subjects the *m*-th factor **removes from the running intersection** must exceed `rmin`, otherwise the combination is flagged redundant (status 2) and dropped. Two consequences the name does not advertise:

- The running count starts at `nx.last = n`, so the **first** factor is tested against the whole sample: a single factor covering all but `<= rmin` subjects is itself "redundant" (the same rule that drops near-full-sample cuts).
- The rule is **order-dependent** in principle (columns are visited in `selected_cols` order, i.e. increasing column index), and it counts **subjects**, in the sample the search is run on — not a proportion, and not a quantity defined on the DGM.

The wrapper's population reading (`.fs_oc_redundant()` in `R/fs_oc_family.R`) keeps the same walk, starting from proportion 1 and testing `p.last − p.this <= rmin / n`: "removes more than `rmin` subjects *of a trial of size n*". This is not the only possible reading; see (c) for the raw-count alternative.

### (b) Does `forestsearch()` pass an `rmin` through to `subgroup.search()`?

**No — the search always receives `subgroup.search()`'s own default `rmin = 5`, except under `sg_focus = "maxeff"`, where `forestsearch()` sets it to 0.** From source:

- `forestsearch()` has **no `rmin` formal and no `...`** (`R/forestsearch_main.R` L1140–1226; `grep "\.\.\."` over the formals region: none).
- `args_call_all <- mget(names(formals()), envir = environment())` (L1339–1340) — the call bundle is exactly the formals, so it cannot carry `rmin`.
- `search_args <- modifyList(args_call_all, search_overrides); valid_args <- names(formals(subgroup.search)); search_args <- search_args[names(search_args) %in% valid_args]; do.call(subgroup.search, search_args)` (L2946–2951) — `rmin` reaches the search only if `search_overrides` contains it.
- The only writer is L2925: `if (!.admit_applies[["effect"]]) search_overrides$rmin <- 0`, with the comment (L2922–2924) "maxeff: relax redundancy pruning to exact-duplicate only. rmin is not a forestsearch formal, so it is set here ... The default (subgroup.search rmin = 5) is preserved for every other focus."
- `subgroup.search()` formal: `rmin = 5` (`R/subgroup_search.R` L85).

So for the driver's `sg_focus = "maxeffCons"` the search ran with **`rmin = 5` subjects at n = 500**, which the wrapper translates to `5/500 = 0.01` in population proportion. (Note for the record: the wrapper honours an `rmin` entry in `forestsearch_args` if one is given; `forestsearch()` itself would ignore such an entry, since it is not a formal. No driver sets one.)

### (c) What `rmin` alone and `n.min` alone remove, over all 2628 combinations

Each floor was re-evaluated independently on every combination (`fixture_run_fs_oc_predict_2026-08-28.R` §3), using the same `Z` (72 columns from 36 cut expressions) the enumeration used:

| floor, evaluated alone on all 2628 | removed | singles / pairs | also empty-intersection | also below `n.min/n` |
|---|---|---|---|---|
| `rmin` — wrapper's rule, shrink `<= rmin/n = 0.01` in proportion | **99** | 0 / 99 | 0 | 22 |
| `rmin` — raw subject count on `df_super` (`rmin = 5` of 5000) | 85 | — | — | — |
| `n.min` — `Pg < 60/500 = 0.12` | **772** | 4 / 768 | 119 | — |
| empty intersection (status 0) | 119 | 0 / 119 | — | 119 (all) |
| `minp = 0.025` per factor | 0 | — | — | — |

Sequentially, in the search's order (empty → minp → rmin → size), the wrapper's counts reproduce: 119 / 0 / 99 / 631, kept 1779 (asserted equal in the script). The 631 in the build report is `n.min`'s marginal contribution after the 119 empty and the 22 rmin-and-size combinations have already gone; `n.min` alone is 772, of which 119 are the empty intersections (prevalence 0) — i.e. **653 genuinely small non-empty candidates**, plus 119 that are empty either way.

The raw-count reading (85) is smaller than the proportional one (99) because on the 5000-row population frame "more than 5 subjects" is `> 0.1%`, a much weaker requirement than "more than 5 of 500" (`> 1%`). The wrapper's choice is the one that means at `n` what the search's rule means at `n`; it is stated in the roxygen and the build report, and the difference is 14 candidates out of 2628.

---

## 2. The prediction on the fixture

`fs_oc_predict(family = <M = 1601 enumeration>, n = 500, c1 = 30, c2 = 10, consistency_method = "split", draws = 2e5, seed = 20260825)` — the seed is the document's `worked-predictions` seed, chosen so the only differences from the document column are the family and the RNG stream length (M = 1601 versus 16 columns of draws), not the seed.

| quantity | wrapper: M = 1601 population family, `"split"` gate | document `worked-predictions`: M = 16, `"split"` | measured n = 500 record (`m500`), driver `"resample"` screen |
|---|---:|---:|---:|
| det_rate | 1.000 (MC SE 0.0000: no non-declaring draw in 2×10⁵) | 0.964 | 1.00 |
| EnH | 70.8 | 132.8 | 72 |
| Esens | 0.164 | 0.329 | 0.17 |
| Espec | 0.870 | 0.767 | 0.87 |
| Eppv | 0.400 | 0.359 | 0.41 |
| Enpv | 0.664 | 0.695 | 0.66 |
| EbetaH | 31.75 | 31.19 | 31.76 (= `oracle − oracle_bias`, the document's own back-out) |
| Enaive_bias | 75.4 | 30.9 | 74.27 |
| M | 1601 | 16 | — (sample-realized, pre-screened family; not recorded as a count) |

Settings: `n = 500`, `c1 = 30`, `c2 = 10`, `draws = 2e5`, `seed = 20260825`; runtime 634 s for the prediction step on M = 1601 (the 2×10⁵ × 1601 draw matrices are ~2.6 GB each). Log: `fixture_run_fs_oc_predict_2026-08-28.log`.

**Read the comparison with this in mind.** The wrapper's gate is the document's **single-split gate** (`"split"`: `(W1 + W2 >= 2·c1) & (W1 >= c2) & (W2 >= c2)`, one half-sample pair standing in for the consistency rate), and the document column uses the same gate on its M = 16 family. The **measured column came from the driver's closed-form `"resample"` screen** (`consistency_method = "resample"`, `pconsistency.threshold = 0.90`, `fs.splits = 400` unused on that path) applied to sample-fitted candidates after the GRF pre-screen, the effect screen and the statistics-keyed near-duplicate removal, with cuts at *sample* quantiles. So the wrapper-vs-measured comparison carries (i) the gate difference (split vs resample), (ii) the family difference (population enumeration of the full cut space vs the sample-realized, pre-screened family), and (iii) the two comparison columns' own differences in what "E[beta(Hhat)]" is (the document's `Enpv`/`EbetaH` are population functionals; the measured `EbetaH` is the document's own back-out `oracle − oracle_bias` from the `m500` block, as `appendix-n500` defines it). The build report's §1 records why (i) is not closable in this task.

Both comparison columns were obtained by **evaluating the document's own chunks** (`anchor`, `worked-scenario`, `worked-predictions`, `appendix-n500`) in the script and reading `det_rate, EnH, ..., Enaive_bias`, `M`, `m500` and `bH500` from that environment — not transcribed from any report or from the task.

**What the table shows, stated without interpretation beyond the numbers.** With the family replaced by the population enumeration of the driver's actual cut space (M = 1601 instead of the document's stylized 16), the wrapper's column sits next to the measured column on every row — declaration 1.000 vs 1.00, E|Ĥ| 70.8 vs 72, sensitivity 0.164 vs 0.17, specificity 0.870 vs 0.87, PPV 0.400 vs 0.41, NPV 0.664 vs 0.66, E[β(Ĥ)] 31.75 vs 31.76, naive bias 75.4 vs 74.3 — where the document's M = 16 column was off by a factor of ~2 on E|Ĥ|, sensitivity and the naive bias. The two columns still differ in gate (split vs resample) and in family (population enumeration vs sample-realized, pre-screened), so this is not a validation of either gate; it is the statement that, at this fixture and these floors, the family size is what moved the document's Step-4 figures, and the gate difference did not move them visibly. The wrapper's top selections are all at the size floor (`Pg` 0.120–0.126) with purity 0.4–1.0 and true means 32–40; selection mass on rules with true mean below `c1` is 0.265; 333 of the 1601 candidates have single-split declaration probability above 0.5 (maximum 0.618).

---

## 3. Provenance

- Branch `feature/glm-extension`; HEAD before this work `438a2e81`; tree clean at start.
- Package as installed and as loaded by the script: `0.2.3` (`devtools::load_all`).
- Commit: 1334a6dd (this paragraph recorded in the follow-up commit) (script, log, this report — explicit paths). No push.
