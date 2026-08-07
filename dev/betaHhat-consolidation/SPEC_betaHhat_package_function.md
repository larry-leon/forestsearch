# SPEC: `betaHhat` as package functions

**Target** `R/betaHhat_truth.R` (new), plus a contract-test document.

**Scope** Membership resolution, the record contract, and `n_eff` reporting.
**No estimand changes.** Each outcome family keeps exactly the target it
computes today.

Proposed home for this spec and its verification: `dev/betaHhat-consolidation/`.

---

## Why this moves into the package

Four near-identical copies exist and none is in `R CMD check`:

| file | state |
|---|---|
| `quarto/simulations/gbsg_redux/betaHhat_truth.R` | survival |
| `dev/identifier-alignment/rerun/betaHhat_truth.R` | **byte-identical** to the above |
| `quarto/simulations/actg175/binary/betaHhat_truth_glm.R` | binary/OR |
| `quarto/simulations/actg175/continuous/betaHhat_truth_md.R` | continuous/MD |

Zero occurrences of `betaHhat`, `beta_region` or `attach_betaHhat` in `R/`.
Zero in the audit's four `verify_*.R` scripts. No roxygen, no man page, no
tests.

Every defect found to date lives in logic that is identical across all four.

---

## The four defects being fixed

All measured, none inferred.

**D1 — `" & "` split before disjunction dispatch.** `get_dfpred()`
(`R/forestsearch_helpers.R`) tests `length(labels) == 1L && grepl("|", labels,
fixed = TRUE)` **first**, on the unsplit string, and its disjunction branch does
its own per-conjunct `" & "` split, paren strip and negation. Splitting first
makes `length > 1`, skips that branch, and shreds the rule. Corrected in the
binary and continuous modules at `e6f6024`; **still live in both survival
copies**.

**D2 — partial-`NA` membership is consumed rather than rejected.**
`evaluate_comparison()` warns and yields `NA` for a missing column.
`get_dfpred()` then propagates it differently by branch:

- conjunction: `in_harm <- in_harm & member` — `FALSE & NA = FALSE`,
  `TRUE & NA = NA`
- disjunction: `in_any <- in_any | in_cj` — `FALSE | NA = NA`,
  `TRUE | NA = TRUE`

`ifelse(in_harm, 0L, 1L)` then carries `NA` into `treat.recommend`. Measured
consequences differ by module and both are wrong: the survival module returns a
**finite hazard ratio for the wrong region** (13,777 rows; verified equal to the
rule with the noise clause deleted, to 7 dp), while
`betaHhat_one_md()` **raises an uncaught error** (`any(idx)` on an all-`NA`
vector is `NA`). Neither returns the all-`NA` record it should.

**D3 — the partition invariant is asserted nowhere.** $\widehat H$ and
$\widehat H^{c}$ must partition the frame, so `nH_eval + nHc_eval` must equal
`nrow(frame)` unconditionally; if no subgroup is identified, $\widehat H^{c}$ is
the ITT population. Under D2 the `NA` rows fall into neither side and the sum
silently drops below `N`. Only `acceptance_betaHhat_md.qmd` checks this, and
only incidentally.

**D4 — `n_eff` is not reported with the target.** `.coverage_meta()` computes
`ok <- is.finite(target) & is.finite(lo) & is.finite(hi)` and
`n_eff <- sum(ok)`, so a non-finite target is dropped with no error and no
warning. Measured effect: `C_betaHhat` computed on 52.2% of detections at
$k = 3$ and 33.9% at $k = 6$, printed beside `C_dagger` computed on ~100%. The
surviving subset is **not random** — it inflates coverage on the consistency arm
(0.9778 against a subsample band of [0.9493, 0.9715]) and deflates it on dina
and grf (0.6525 against [0.7125, 0.7775]).

---

## What is NOT in scope

- **No estimand changes.** Survival keeps the population `coxph` value on a
  fixed evaluation frame. Binary keeps its `glm` fit. Continuous keeps
  `compute_aor()` on `df_super`. Whether the potential-outcome form should
  replace the fitted value for MD/OR is a separate question and is deliberately
  excluded.
- **The population-noise `df_super` scheme** is separate, already designed, and
  fixes a different problem (a rule naming a column the frame lacks). This spec
  makes that failure *loud and correctly recorded*; it does not prevent it.
- **`get_dfpred()` itself is not modified.** It has callers well beyond these
  modules. The resolver wraps it.

---

## Design

Five layers. All defects live in layer 1.

### Layer 1 — `.fs_resolve_membership()` (internal)

The single shared piece, and the only place rule strings are interpreted.

```r
.fs_resolve_membership(df, rule, sg_def_struct = NULL)
```

Returns a list, never a bare vector:

```r
list(in_region = <logical(nrow(df))> or NULL,
     status    = "ok" | "unresolved" | "empty",
     missing   = <character(0) or the unresolvable column names>)
```

Contract:

1. **Structured first.** If `sg_def_struct` is supplied, evaluate membership via
   `.grf_evaluate_subgroup()` (`R/forestsearch_helpers.R:1778, 1781, 1784`) and
   do not touch the label string. This is the package's own answer for GRF
   membership, chosen because `get_dfpred()` cannot evaluate the disjunction
   string — see the comment at `:1772`. String parsing is the fallback, not the
   primary path.
2. **String path mirrors `get_dfpred()`'s dispatch order exactly.** Test
   `length(rule) == 1L && grepl("|", rule, fixed = TRUE)` **before** any
   `" & "` split, and pass the disjunction through intact. Fixes D1.
3. **Reject partial `NA`.** After resolution, if `any(is.na(in_region))`,
   return `status = "unresolved"`, `in_region = NULL`, and `missing` populated
   by testing each referenced column against `names(df)`. **Never return a
   vector containing `NA`.** Fixes D2. The resolver must not attempt repair —
   a rule the frame cannot express has no membership, and inventing one is
   worse than reporting none.
4. **Empty is distinct from unresolved.** A rule that resolves to zero members
   returns `status = "empty"` with a real all-`FALSE` vector. That is a valid
   region; the target is `NA` but the partition still holds
   (`nH = 0`, `nHc = N`).
5. `missing` is populated on `status = "unresolved"` so the caller can
   distinguish an unresolvable rule from a degenerate one without re-deriving
   it.

### Layer 2 — `.fs_region_effect()` (internal)

Per-family target on a region. Dispatch on `outcome_type` / `effect_measure`.
**Each branch reproduces exactly what its module computes today.**

| family | branch |
|---|---|
| survival | `survival::coxph(Surv(time, event) ~ treat)` on the evaluation frame, guarding `<5` events and single-arm with `NA` |
| binary | `stats::glm(y ~ treat, binomial)`, guarding separation and minimum cell counts with `NA` |
| continuous / count | `compute_aor(df, subset_indicator, effect_measure)` |

No branch may restate arithmetic available from `compute_aor()` or
`compute_ahr()`. The package already carries two sites for the GLM marginal
effect (`compute_aor()` and the `.effect_one()` closure inside
`.compute_glm_effects()`); do not add a third, and do not refactor those two
under this spec.

### Layer 3 — `fs_betaHhat_one()`

The H / Hᶜ pair, and the only place the invariant is enforced.

```r
fs_betaHhat_one(rule, frame, focus, outcome_type, effect_measure, ...)
```

Returns one record with a fixed schema, scale-agnostic so downstream
`paste0("betaHhat_", suffix)` code is shared verbatim:

```
betaHhat_H, betaHhat_Hc, nH_eval, nHc_eval, status, missing_cols
```

**HARD INVARIANT, asserted on every call:**

```
nH_eval + nHc_eval == nrow(frame)
```

Three cases, exhaustive:

| case | `nH_eval` | `nHc_eval` | targets | status |
|---|---|---|---|---|
| resolved | `sum(in_region)` | `N - nH_eval` | computed | `"ok"` |
| no subgroup identified | `0` | `N` | `betaHhat_Hc` = ITT effect | `"ok"` |
| unresolvable rule | `NA` | `NA` | both `NA` | `"unresolved"` |

The invariant is checked on the first two and skipped on the third, where both
counts are `NA` by construction. It must be an assertion that stops, not a
warning: a partition failure means membership is incoherent and every number
downstream of it is meaningless.

`focus` stays a required argument with no default, erroring when missing, until
the harm/benefit convention is settled. `focus = "harm"` reads
`treat.recommend == 0L` as in-region; `"benefit"` inverts.

### Layer 4 — `fs_betaHhat_table()` and `fs_attach_betaHhat()`

Deduplication is **by distinct realized rule**, not by replicate.
`fs_attach_betaHhat()` remains "the one call each engine adds to `run_cell()`".

Both return, as attributes and as columns in the bundle `meta`:

```
n_rules_total, n_rules_resolved, n_rules_unresolved,
n_reps_total,  n_reps_resolved,  n_reps_unresolved
```

Fixes D4 at the source. Under the population-noise `df_super` scheme
`n_reps_unresolved` should be zero, which makes it a regression check rather
than dead weight.

### Layer 5 — `fs_betaHhat_neff_parity()`

Generalise `betaHhat_neff_parity_or()` (`e6f6024`) to all three pathways and
wire it into the survival consumers, which never received it.

The coverage assembly **must refuse to print** a target whose `n_eff` differs
from its siblings in the same cell unless the caller passes an explicit
acknowledgement. Silence is what produced D4's measured consequence. `strict =
TRUE` stops; `strict = FALSE` warns, for a run whose cause is known.

Keep the diagnostic companion (`betaHhat_neff_report_*()`): a parity break can
mean an unresolvable rule *or* a legitimately degenerate region, and the two
need different responses.

---

## Export policy

Follow `compute_aor()`'s pattern: `@keywords internal` **and** `@export`, so the
functions are usable by anyone building simulations on `generate_glm_dgm()` or
`setup_gbsg_dgm()` without cluttering the reference index. `.fs_*` internals
stay unexported.

`Roxygen: list(markdown = TRUE)` is set, so write literal `% < > &` in roxygen
text and never Rd-escape. Keep `@section` titles plain.

---

## Contract tests

Static properties, which the working conventions admit as the standing
exception to "no `testthat` scaffolding". Everything else stays integration-style
Quarto.

**Every test must pair the assertion with a plausible-but-wrong form asserted to
fail.** Otherwise it cannot distinguish a correct implementation from one that
happens to be close. And every test region must have a value distinguishable
from zero, from the degenerate values, and from every other test region — a
collision makes the assertion vacuous.

| id | property | negative control |
|---|---|---|
| T1 | `nH_eval + nHc_eval == N` for every rule shape × every family | a deliberately partial-`NA` membership must **fail** the sum |
| T2 | a rule naming a missing column returns `status = "unresolved"` and an all-`NA` record | the pre-fix `get_dfpred()` output must be shown to contain `NA`, i.e. the guard is doing work |
| T3 | a GRF disjunction and its structured `sg_def` give **identical** membership | the split-first form must **fail** to resolve the same string |
| T4 | negation: `!{x <= c}` and `{x > c}` give identical membership; a rule and its negation partition | — |
| T5 | no subgroup ⇒ `nHc_eval == N` and `betaHhat_Hc` equals the ITT effect | a run that returns a subgroup must **not** give `nHc_eval == N` |
| T6 | `n_reps_unresolved` matches the count of unresolvable rules in a synthetic bundle | — |
| T7 | the parity guard **fires** on an injected dropped target and on a missing row; `strict = FALSE` warns without stopping | it must **pass** on a clean bundle |
| T8 | the same rule on the same frame gives identical membership regardless of `outcome_type` — membership is family-independent | — |
| T9 | repeated evaluation of a resolved target is bit-identical (`max - min == 0`, not a tolerance) | — |
| T10 | every family's target reproduces its **current** module's value on a fixed fixture, to `0` for exact families and to Monte Carlo error for survival | — |

T10 is the migration gate. Nothing may move until each family reproduces what
it produces today, because this spec changes no estimand and any movement is a
defect in the consolidation.

T3 needs a *realized* GRF disjunction, not a hand-written one. If none exists in
the committed bundles — 0 of 250,025 non-empty `sg_def` rows contain `"|"` — one
must be generated by running `forestsearch()` at settings that produce a
multi-leaf policy tree (`maxdepth >= 2`). A hand-written string tests the parser
and not the path.

---

## Migration

1. Land the package functions and the contract tests. Nothing else changes.
2. Run T10 against all three families on fixed fixtures. **Gate.**
3. Convert the four simulation modules to thin shims that call the package
   functions, keeping their current names so no sweep `.qmd` changes yet.
4. Re-point the sweeps at the package functions directly; delete the shims and
   the duplicate `dev/identifier-alignment/rerun/betaHhat_truth.R`.

Step 3 is what makes this safe to do incrementally: the shims let T10 run
against real sweeps before any `.qmd` is touched.

The duplicate survival copy should not survive step 4. Decide whether
`dev/identifier-alignment/rerun/` should retain a copy at all now that the audit
it served is complete.

`R CMD check --as-cran` clean is the bar for every commit that touches `R/`.
