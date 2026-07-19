# Patch spec — add `sg_focus = "maxeff"` (pure max-effect selection)

**Target:** `forestsearch` (v0.2.0)
**Scope:** 4 code edits + 2 doc updates + 5 tests. No new dependencies, no
change to any existing selection path.

---

## 1. Why

`forestsearch()` currently offers no selection rule that ranks candidates by
**effect alone**. The closest, `sg_focus = "hr"` (user alias `"eff"`), sorts
lexicographically by `(-Pcons, -hr, K)` — consistency is the *primary* key and
the effect only breaks ties. The `hrMaxSG`/`hrMinSG` family ranks by size within
an effect band. So the argmax-of-effect primitive is absent.

That gap matters for two reasons:

1. **Guo & He (2021) comparability.** Their correction is built on
   $\hat H = \arg\max_g \hat\beta(g)$ — the estimand $\beta_s$, the offsets
   $d_i=(1-n^{r-1/2})(\hat\beta_{\max}-\hat\beta_i)$, and the limiting law
   $\max_{i\in H}T_i$ all reference that functional. A lexicographic
   `(-Pcons, -hr)` order is not it. On the GBSG `sg_focus = "eff"` run the two
   happened to coincide, but only because `Pcons` saturates near 1 for many of
   the 1744 candidates, leaving `hr` to decide. That is a numerical accident of
   the dataset, not a property of the rule.
2. **The gate already expects it.** `.fs_dg_reselection_from_focus()`
   (`R/fs_debias_gate_methods.R:142`) contains `maxeff = "maxeff"`, and
   `.fs_dg_select()` (`R/fs_debias_gate.R:146`) implements
   `maxeff = passers[which.max(beta[passers])]`. The Tier-2 half of this option
   is already built; only the search side is missing.

**Semantics.** `sort_subgroups()` receives `result_new`, which is built from
`results_list_valid` — candidates that already cleared
`pconsistency.threshold` (`R/subgroup_consistency_main.R:938`, and the message
at line 941 reads "subgroups passed consistency threshold"). So `"maxeff"` is
**argmax of effect among candidates with `Pcons >= pconsistency.threshold`** — a
*constrained* maximiser, exactly matching the gate's
`passers[which.max(beta[passers])]`.

---

## 2. Edits

### Edit 1 — `R/forestsearch_main.R` (~line 1053)

```r
# before
valid_sg_focus <- c("hr", "hrMaxSG", "maxSG", "hrMinSG", "minSG")
# after
valid_sg_focus <- c("hr", "hrMaxSG", "maxSG", "hrMinSG", "minSG", "maxeff")
```

`valid_sg_focus_user` is built as `c(valid_sg_focus, "eff", "effMaxSG",
"effMinSG")`, so `"maxeff"` propagates into the error message automatically. Do
not add it separately.

### Edit 2 — `R/subgroup_consistency_main.R` (~line 971)

There is a **second, independent** whitelist here. Missing it means the fit
succeeds and then stops inside the consistency step.

```r
# before
valid_sg_focus <- c("hr", "hrMaxSG", "maxSG", "hrMinSG", "minSG")
# after
valid_sg_focus <- c("hr", "hrMaxSG", "maxSG", "hrMinSG", "minSG", "maxeff")
```

### Edit 3 — `R/subgroup_consistency_helpers.R`, `sort_subgroups()` (~line 547)

Insert **before** the existing `if (sg_focus == "hr")` branch:

```r
  # Pure max-effect: rank by effect alone among candidates that already
  # cleared pconsistency.threshold.  Distinct from "hr", which sorts
  # (-Pcons, -hr, K) and therefore lets consistency dominate.
  if (sg_focus == "maxeff") {
    data.table::setorder(result_new, -hr, K)
    return(result_new)
  }
```

No `effect_log_scale` handling is required: ordering is invariant under the
monotone log transform, so `-hr` is correct whether `hr` holds a ratio or a
log-effect. (`hrMaxSG` needs `exp()` only because its *band* is multiplicative.)

### Edit 4 — `R/subgroup_consistency_helpers.R`, `sort_subgroups_preview()`

Insert **before** that function's `if (sg_focus == "hr")` branch:

```r
  if (sg_focus == "maxeff") {
    data.table::setorder(result_new, -HR, K)
    return(result_new)
  }
```

Pre-consistency there is no `Pcons` column, so this is identical to the existing
`"hr"` preview branch — intentionally. The two rules differ only *after*
consistency is evaluated, which is precisely the distinction being added.

### Docs

- `R/forestsearch_main.R` `@param sg_focus` (~line 414): add `"maxeff"` —
  *ranks by effect alone among candidates passing `pconsistency.threshold`;
  contrast `"hr"`, which ranks by consistency first and uses effect only as a
  tiebreaker.*
- `R/subgroup_consistency_main.R` roxygen (~lines 97–99, the sort-order list):
  add `\item{\code{"maxeff"}}{Sort by \eqn{(-hr, K)}; pick top.}`

### Do NOT change

- `.normalize_sg_focus()` — unknown values pass through, and `"maxeff"` is
  already the canonical token on the gate side.
- `.validate_selection_rule()` — it only constrains non-`"neighborhood"` rules,
  so `maxeff` + default `neighborhood` passes.
- `.fs_dg_reselection_from_focus()` — line 142 already maps `maxeff -> "maxeff"`.

---

## 3. Tests

Run on GBSG, matching the existing analysis configuration
(`use_lasso = FALSE`, `use_grf = FALSE`, `maxk = 2`, `n.min = 60`,
`pconsistency.threshold = 0.90`, `conf_force = c("er <= 0", "pgr <= 0")`,
`conf.cont_jcuts = list(er = 10, pgr = 10)`).

**T1 — accepted.** `forestsearch(..., sg_focus = "maxeff")` completes without
error and returns a non-NULL `sg.harm`.

**T2 — semantics (the decisive test).** From the fitted object's candidate
table (locate it on `fs$grp.consistency`; it carries `Pcons`, `hr`, `N`, `K`),
verify the selected subgroup is the row attaining
`max(hr)` among rows with `Pcons >= 0.90`. Assert equality of the subgroup
label, not just of `hr`, so ties are caught.

**T3 — regression, existing focuses untouched.** Fit with each of `"hr"`,
`"effMaxSG"`, `"maxSG"` at a fixed `seedit` and confirm the selected subgroup and
its naive effect are **identical to the pre-patch values**. Capture those before
applying any edit.

**T4 — `hr` vs `maxeff` on GBSG.** Report both selected subgroups side by side.
If they agree, that confirms `Pcons` saturation is why the earlier
`sg_focus = "eff"` run coincided with an unconstrained argmax. If they differ,
the distinction is live on this dataset. Either outcome is a result — report it,
do not tune toward one.

**T5 — gate mapping.** Fit with `sg_focus = "maxeff"`, `debias_gate = TRUE` and
confirm the gate reports **`Re-selection rule = maxeff`**. This checks the
already-present mapping at `fs_debias_gate_methods.R:142` fires.

---

## 4. Follow-on (not part of this patch)

Once `"maxeff"` exists, `quarto/GuoHe/guohe_from_forestsearch.R` should tighten
its scope gate to

```r
in_scope <- identical(focus, "maxeff")
```

with `"hr"`/`"eff"` demoted to force-only, since `(-Pcons, -hr)` is lexicographic
rather than an effect argmax. The driver
`quarto/GuoHe/run_guohe_gbsg_maxeff.R` then switches to
`sg_focus = "maxeff"`, and its filename finally matches what it does.

**One caveat to document, not to fix.** Even `"maxeff"` is an argmax over a
**consistency-screened** set. Guo & He's offsets correct the argmax; they do not
correct the screen's own selection effect. That residual is real and is
something MR captures (its `pass <- which(bs >= t_g)` is recomputed per draw)
and Guo & He does not. Worth one sentence in the docs rather than silence.
