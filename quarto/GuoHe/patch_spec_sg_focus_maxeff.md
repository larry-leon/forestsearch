# Patch spec — `sg_focus = "maxeff"`: strict effect-maximiser selection

**Target:** `forestsearch` (v0.2.0)
**Revision 2.** Supersedes both earlier drafts. Revision 1 listed five gates;
execution tracing found a sixth (candidate-pool truncation), which is
unconditional and which `stop_threshold = NULL` does **not** disable.

---

## 1. Requirement

Guo & He (2021) correct the **unconditional maximiser**

$$\hat H \;=\; \arg\max_{g\in\mathcal F}\hat\beta(g)$$

over a supplied family $\mathcal F$. Their estimand $\beta_s$, their offsets
$d_i=(1-n^{r-1/2})(\hat\beta_{\max}-\hat\beta_i)$ anchored to
$\hat\beta_{\max}$, and their limiting law $\max_{i\in H}T_i$ all reference that
functional.

`sg_focus = "maxeff"` must therefore select **the effect maximiser subject to no
other condition**, *by construction* — so that a dataset on which the auxiliary
conditions happen not to bind is indistinguishable, in code, from one on which
they do.

**Excluded** (each changes which candidate wins, or prevents the family being
fully evaluated): effect threshold, consistency threshold, consistency-based
early stopping, two-stage screening, candidate-pool truncation.

**Retained:** `n.min`, `d0.min`, `d1.min`. These are **estimability**
constraints — they determine which candidates are well-posed subgroups at all,
constituting $\mathcal F$ rather than screening within it. Guo & He carry the
same notion; a candidate whose within-subgroup model cannot be fit is not a
competitor.

---

## 2. The six gates

All verified in source. Gate 6 was found by execution, not by reading — see §5.

| # | gate | location | mechanism |
|---|---|---|---|
| **G1** | effect threshold | `subgroup_consistency_main.R:508,511` | `hr.subgroups[hr.subgroups$HR >= hr.threshold, ]`, applied *before* consistency |
| **G2** | consistency drop, single-stage | `subgroup_consistency_helpers.R:1439` | `if (isTRUE(p.consistency < pconsistency.threshold)) return(NULL)` |
| **G3** | consistency drop, two-stage | `subgroup_consistency_helpers.R:1779` | same test, `return(NULL)` |
| **G4** | early stopping | `stop_threshold` | halts the search once a candidate clears; for `sg_focus %in% c("hr","minSG")` `batch_size_parallel <- 1L`, so it stops at the **first passer in preview order** (`subgroup_consistency_main.R:750`) |
| **G5** | two-stage screening | `use_twostage`, `ts_params$screen.threshold` | Stage 1 prunes on consistency before effects are compared |
| **G6** | **candidate-pool truncation** | `subgroup_consistency_main.R:565-567` | `maxsgs <- min(nrow(found.hrs), stop_Kgroups); found.hrs <- found.hrs[seq_len(maxsgs), ]` — **unconditional**, with `stop_Kgroups = max_subgroups_search`, applied to the `sg_focus`-ordered preview list |

G4 and G6 do not filter a completed table: G4 stops the enumeration early, G6
never admits most of it. An argmax over a partially evaluated family is not an
argmax.

G4 also explains, independently of the lexicographic sort, why `"hr"` is not
Guo & He's functional: with `batch_size = 1` it returns the **first candidate in
`-HR` order clearing `stop_threshold`**, not the maximiser over the evaluated
set.

---

## 3. Edits

### Edit 1 — `R/forestsearch_main.R`, whitelist (~line 1053)

```r
valid_sg_focus <- c("hr", "hrMaxSG", "maxSG", "hrMinSG", "minSG", "maxeff")
```

`valid_sg_focus_user <- c(valid_sg_focus, "eff", "effMaxSG", "effMinSG")` picks
it up automatically.

### Edit 2 — `R/forestsearch_main.R`, regime enforcement

Immediately after `sg_focus` is validated, before any search argument is used:

```r
  # sg_focus = "maxeff" implements Guo & He's argmax primitive: the effect
  # maximiser over the enumerated family, subject to NO auxiliary selection
  # condition.  Consistency screening, effect thresholding, early stopping,
  # two-stage pruning and candidate-pool truncation would each change which
  # candidate wins (or stop the family being fully evaluated), so all are
  # disabled BY CONSTRUCTION rather than left to the caller.  Size/event minima
  # (n.min, d0.min, d1.min) are retained: they define which candidates are
  # estimable, not which one wins.
  if (identical(sg_focus, "maxeff")) {
    .ov <- character(0)
    if (!isTRUE(all.equal(pconsistency.threshold, 0)))
      .ov <- c(.ov, "pconsistency.threshold -> 0")
    if (!is.null(stop_threshold))       .ov <- c(.ov, "stop_threshold -> NULL")
    if (isTRUE(use_twostage))           .ov <- c(.ov, "use_twostage -> FALSE")
    if (is.finite(max_subgroups_search))
      .ov <- c(.ov, sprintf("max_subgroups_search %g -> Inf", max_subgroups_search))
    pconsistency.threshold <- 0
    stop_threshold         <- NULL
    use_twostage           <- FALSE
    max_subgroups_search   <- Inf
    if (length(.ov)) {
      warning("sg_focus = 'maxeff' selects the effect maximiser with no ",
              "auxiliary conditions; overriding: ",
              paste(.ov, collapse = ", "), ".", call. = FALSE)
    }
    if (identical(consistency_method, "split")) {
      warning("sg_focus = 'maxeff' evaluates consistency for EVERY candidate ",
              "(no screening, no early stopping, no truncation). With ",
              "consistency_method = 'split' that is fs.splits x 2 refits per ",
              "candidate and is usually impractical; 'resample' is a single fit ",
              "plus a closed form.", call. = FALSE)
    }
  }
```

Why each neutraliser is safe:

- `pconsistency.threshold <- 0` disables **G2/G3**: the tests are
  `p.consistency < 0`, false for any valid probability, and `isTRUE(NA < 0)` is
  `FALSE`, so inestimable-rate candidates are retained rather than dropped.
- `stop_threshold <- NULL` is the documented disable for **G4**.
- `use_twostage <- FALSE` disables **G5**.
- `max_subgroups_search <- Inf` makes **G6** inert:
  `min(nrow(found.hrs), Inf) == nrow(found.hrs)`.

**Do not** substitute a structural argument for the `max_subgroups_search`
override. Under `maxeff` the preview order is `-HR`, so top-*K*-by-`HR` would
contain argmax-by-`HR` *provided* preview `HR` and final `hr` are the same
quantity — which is unverified. `Inf` removes the dependency.

### Edit 3 — `R/subgroup_consistency_main.R`, whitelist (~line 971)

A **second, independent** whitelist. Omitting it means the fit proceeds and then
stops inside the consistency step.

```r
valid_sg_focus <- c("hr", "hrMaxSG", "maxSG", "hrMinSG", "minSG", "maxeff")
```

### Edit 4 — `R/subgroup_consistency_main.R`, G1 bypass (~lines 505-517)

`hr.threshold = 0` is **not** a usable neutraliser: `log(hr.threshold)` becomes
`-Inf` and may propagate (`forestsearch_helpers.R:1288`). Guard explicitly:

```r
  if (identical(sg_focus, "maxeff")) {
    # Guo & He's argmax ranges over the whole family; an effect floor would make
    # it a constrained maximiser.
    found.hrs <- hr.subgroups
  } else {
    # ... existing hr.threshold filtering, both branches, unchanged ...
  }
```

### Edit 5 — `R/subgroup_consistency_helpers.R`, `sort_subgroups()` (~line 547)

Insert **before** the `"hr"` branch:

```r
  # Effect maximiser, no auxiliary condition.  Distinct from "hr", which sorts
  # (-Pcons, -hr, K) and lets consistency dominate.  K breaks exact ties toward
  # the simpler rule; it cannot override the effect ordering.
  if (sg_focus == "maxeff") {
    data.table::setorder(result_new, -hr, K)
    return(result_new)
  }
```

No `effect_log_scale` handling: ordering is invariant under the monotone log
transform. (`hrMaxSG` needs `exp()` only because its *band* is multiplicative.)

### Edit 6 — `R/subgroup_consistency_helpers.R`, `sort_subgroups_preview()`

```r
  if (sg_focus == "maxeff") {
    data.table::setorder(result_new, -HR, K)
    return(result_new)
  }
```

Pre-consistency there is no `Pcons`, so this coincides with the existing `"hr"`
preview branch — intentionally. The rules diverge only once consistency exists.

### Docs

- `@param sg_focus` (`forestsearch_main.R` ~line 414): add `"maxeff"` — *the
  effect maximiser over the enumerated family, with no consistency threshold,
  effect threshold, early stopping, two-stage screening or candidate-pool
  truncation. Implements the argmax primitive of Guo and He (2021). Size/event
  minima still apply. Contrast `"hr"`, which ranks by consistency first, uses
  effect only as a tiebreaker, and under early stopping returns the first passer
  in effect order.*
- `@param pconsistency.threshold`, `stop_threshold`, `use_twostage`,
  `max_subgroups_search`: note each is overridden under `sg_focus = "maxeff"`.
- `subgroup_consistency_main.R` roxygen (~lines 97-99):
  `\item{\code{"maxeff"}}{Sort by \eqn{(-hr, K)}; pick top.  No consistency filter.}`

### Do NOT change

`.normalize_sg_focus()` (unknown values pass through; `"maxeff"` is already the
canonical gate token), `.validate_selection_rule()` (constrains only
non-`"neighborhood"` rules), `.fs_dg_reselection_from_focus()`
(`fs_debias_gate_methods.R:142` already maps `maxeff -> "maxeff"`).

---

## 4. Tests

Configuration as in `quarto/GuoHe/fit_sg_focus_baseline.R`, at
**`max_subgroups_search = 2000`** so T4 matches the recorded baseline.

**T1 — accepted.** `forestsearch(..., sg_focus = "maxeff")` completes, returns
non-`NULL` `sg.harm`, and emits the override warning.

**T2 — it is the global argmax.** Independently enumerate all $\le$ `maxk`
combinations over the dummy-expanded cut matrix; keep those meeting **`n.min`
AND `d0.min`/`d1.min`** (`subgroup_search.R:581,588` — the gate/bridge family
filters on `n.min` only, so an `n.min`-only reference is not a valid comparison);
fit each within-subgroup effect; take the argmax. `maxeff`'s `sg.harm` **must
equal it**. Assert on the subgroup label, not the effect value, so ties are
caught.

*This is the only test that verifies the requirement itself.* It must pass
regardless of whether any gate binds on this dataset.

**T3 — hostile settings do not move the selection.** Refit `maxeff` with
`pconsistency.threshold = 0.99`, `hr.threshold = 5`, `stop_threshold = 0.5`,
`use_twostage = TRUE`, `max_subgroups_search = 5`. The selection **must be
unchanged**, and the warning must name every override. If it moves, a gate was
missed.

**T4 — regression.** `"hr"`, `"effMaxSG"`, `"maxSG"` must reproduce
`quarto/GuoHe/sg_focus_baseline.rds` exactly (`{er <= 0} & {size <= 35}` / 61 /
2.5369; same; `{er <= 0}` / 82 / 1.9514) at the same `max_subgroups_search =
2000`.

**T5 — `Pcons` retained but not filtering.** Under `maxeff` the candidate table
must contain candidates with `Pcons` below the *user-supplied* threshold —
consistency still computed and reported, just not used to select.

**T6 — gate mapping.** `sg_focus = "maxeff"`, `debias_gate = TRUE` must report
`Re-selection rule = maxeff`.

**T7 — truncation inert.** Fit `maxeff` at `max_subgroups_search = 5`; selection
unchanged and the override named. (Contrast: at depth 50, `maxSG` returned
`NULL` on this dataset purely from truncation.)

---

## 5. Note on method

G6 was not found by reading source; it was found by tracing execution after
`maxSG` returned `NULL`. It is unconditional, sits three functions from where
selection appears to happen, and is invisible to a grep for the parameters that
seem to govern selection. T3 and T7 exist precisely so that a **seventh** gate,
if one exists, announces itself rather than hiding behind a dataset where it
does not bind. Treat a T3/T7 failure as a discovery, not as a test to adjust.

---

## 6. Downstream, once this lands

**Bridge** (`quarto/GuoHe/guohe_from_forestsearch.R`): tighten the scope gate to
`in_scope <- identical(focus, "maxeff")`, demoting `"hr"`/`"eff"` to force-only.
Add a **precondition assertion** that the bridge's own argmax over the
reconstructed family equals `fs$sg.harm`; if not, Guo & He is debiasing an object
the analysis does not report, and the result should be refused rather than
returned quietly.

**Do not re-run the r-sweep.** The bridge performs its own argmax and never reads
`sg_focus`, so all five arms are numerically identical to the committed results.
One confirmatory run at `r = 0.03` suffices to check the assertion.

**Residual, to document not fix.** With G1-G6 disabled, `maxeff` and Guo & He's
primitive coincide by construction. The remaining gap is elsewhere: the
enumerated family derives from **data-dependent cut points** (sample quantiles),
whereas Guo & He's Section 3 assumes a fixed index set $D$. This is the same
conditioning MR makes, so the two stay comparable — but it should be stated
rather than left implicit.
