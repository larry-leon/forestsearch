# CC BRIEF — execute the "gate"/"tier" rename (IMPLEMENTATION)

```
git checkout -b feature/mr-terminology
claude "Read dev/terminology-work/CC_BRIEF_terminology_rename.md and execute it."
```

**Prerequisite:** `dev/terminology-work/TERMINOLOGY_SURVEY.md` must exist. This
brief executes the map the maintainer approved from it. Where this brief and the
survey disagree, **this brief governs** --- the maintainer has since decided to
purge `gate` entirely, which the survey had proposed retaining in three places.

---

## 1. The approved map

### 1.1 API --- no deprecation shim required

Confirmed in the survey (F7, two independent evidence chains): none of these
shipped in CRAN 0.1.0. `gate` occurs zero times in the 0.1.0 reference manual,
and the family first appears 2026-06-08, eleven weeks after publication. **Rename
outright. Do not add `deprecated()` forwards.**

| current | new |
|---------|-----|
| `fs_debias_gate()` | `fs_mr_inference()` |
| `debias_gate` (formal of `forestsearch()`) | `mr_inference` |
| `debias_gate_args` | `mr_inference_args` |
| `out$debias_gate` | `out$mr_inference` |
| `harm_flag_debiased` | `mr_harm_confirmed` |
| `gate_estimates_table()` | `mr_estimates_table()` |

### 1.2 The harm-confirmation names --- also purged

The survey recommended keeping these on the grounds that they name the
harm-confirmation decision accurately. **The maintainer has decided otherwise:**
retaining `gate` for one sense preserves exactly the hazard that motivated the
purge --- a reader must still learn that the word means one specific thing here.

| current | new | what it is |
|---------|-----|-----------|
| `gate` (argument of `fs_debias_gate()`) | `confirm_rule` | selects the harm-confirmation rule: `"point"` or `"ci"` |
| `t_gate` | `t_confirm` | the threshold the rule is applied against, effect scale |
| `c_gate` (formal of `fs_fdr_report()`) | `c_confirm` | vector of such thresholds swept by the FDR report |
| `rep$fdr$c_gate` (result column) | `c_confirm` | as above |
| `.fs_dg_gate_null()` | `.fs_mr_confirm_null()` | near-null default for `t_confirm` |

The permitted values of `confirm_rule` (`"point"`, `"ci"`) are **unchanged** ---
only the argument name moves.

### 1.3 Internal prefix

`.fs_dg_*` → `.fs_mr_*` throughout. Enumerate them from the survey's §4 table and
rename all; report the count.

### 1.4 "tier" → FB / MR

"tier-1" becomes **full bootstrap (FB)**, "tier-2" becomes **multiplier
resampling (MR)**. Inside the package this is documentation-only plus one
user-visible string (`gate_estimates_table.R:106`, the subtitle). Confirm that
against the survey before assuming it.

---

## 2. What must NOT be touched

The survey found these; a name-based substitution would corrupt them.

1. **Calibrator tolerance gates** --- `calibrate_helpers.R`, `calibrate_k_treat.R`
   and four others use "gate" for a tolerance check that `stop()`s. A genuine
   pass/fail, unrelated to MR. **Leave every occurrence.**
2. **`glm_effect_estimators.R:342-441`** uses `tier1`/`tier2` for the risk-
   difference fallback ladder. These are the **only** `tier1`/`tier2` identifiers
   in `R/`, and they have nothing to do with FB/MR. **Leave them.**
3. **`man/*.Rd`** --- generated. Never hand-edit; regenerate with
   `devtools::document()`.
4. **`NEWS.md` entries for released versions.** The 0.1.0 section is the record
   of what shipped. Add a 0.2.0 entry; rewrite nothing above it.
5. **Control-flow guards** the survey classified as sense 6. Leave them.

**Do not use `sed` or any bulk substitution, even scoped.** Rename identifier by
identifier, checking each site. The survey's own inventory shows how many
substring false positives a pattern match produces.

---

## 3. Documentation --- the maintainer's priority

The maintainer's stated reason for the purge is that the current names require
looking up. **A rename that leaves the reader guessing has failed even if every
identifier changes.**

### 3.1 A definitions block in `?forestsearch`

Add an explicit `@section` to `forestsearch()`'s roxygen defining every term,
plain-text title (no markup --- roxygen rejects it). It must state, in the
reader's terms rather than the implementer's:

* **Full bootstrap (FB)** --- what it resamples, what it corrects, which function
  implements it (`forestsearch_bootstrap_dofuture()`).
* **Multiplier resampling (MR)** --- what it resamples, what it corrects, which
  function implements it (`fs_mr_inference()`), and that it is post-selection
  inference on a completed analysis, **not** part of subgroup identification.
* **Harm confirmation** --- that after a subgroup is identified, MR yields a
  de-biased estimate, and `confirm_rule` decides whether that estimate still
  indicates harm against `t_confirm`. Spell out both rules:
  `"point"` compares the de-biased point estimate; `"ci"` compares the one-sided
  95% selection-adjusted lower bound.
* **`t_confirm`** --- the threshold, its scale, and what `NULL` resolves to.
* **`mr_harm_confirmed`** --- what `TRUE`, `FALSE` and `NA` each mean.
  **`NA` must be documented explicitly**: it arises when MR did not run or could
  not be computed, and is *not* evidence against harm. The survey's F1 and F4 are
  both cases where a user could read `NA` as "harm not confirmed".

Also distinguish, in one sentence each, the four senses the survey separated:
admissibility criteria, selection rule, post-selection inference, harm
confirmation. This is the vocabulary that prevents the confusion being renamed
away from.

### 3.2 `@param` blocks

Every renamed argument gets a `@param` that says what it does, not what it is
called. `mr_inference` must state that the **default is `FALSE`** and that the
default analysis is therefore pure identification --- no MR, no harm
confirmation.

### 3.3 Cross-references

`fs_mr_inference()`, `mr_estimates_table()`, `fs_fdr_report()` and
`forestsearch_bootstrap_dofuture()` should `@seealso` one another, so a reader
landing on any one finds the others.

---

## 4. Blast radius outside `R/`

The survey measured ~218 files, **entirely the maintainer's own `quarto/`
notebooks** --- `tests/` and `vignettes/` are clean of these names. Since nothing
shipped, this is mechanical cost rather than a compatibility constraint.

* Update `quarto/` sources (`.qmd`) that set or read the renamed names.
* **Do not** touch rendered `.html` or `_freeze/` caches --- the survey found
  those inflate any count by ~30,000 lines and they regenerate.
* The six simulation cells in `quarto/simulations/gbsg_redux/` set
  `debias_gate`/`gate_draws` and read `gate_ok`. Rename the settings; where a
  local variable like `gate_ok` is the notebook's own, rename for consistency and
  say so.

---

## 5. Test coverage --- do this BEFORE the rename

The survey's standing recommendation, and the reason F1 survived: **`tests/` has
zero coverage of `debias_gate`.** Renaming an untested surface across 218 files
is the riskiest part of this plan.

Add `tests/testthat/test-mr-inference.R` covering, at minimum:

1. Default is off --- `formals(forestsearch)$mr_inference` is `FALSE`, and a
   default run returns `NULL` for `out$mr_inference`.
2. Enabled on the Cox/consistency path, it returns a populated object with the
   documented fields.
3. No subgroup identified → clean skip, `mr_harm_confirmed` is `NA`, no error.
4. `confirm_rule = "point"` and `"ci"` both run and can differ.
5. The re-selection rule derived from `sg_focus` matches the survey's mapping
   table for all eleven accepted foci.

Commit the tests **first**, on the pre-rename names, so they prove the rename is
behaviour-preserving. Then rename tests and code together.

Do **not** attempt to fix F1 (silent no-op on the default GLM path), F2
(recomputation inside bootstrap/CV), or F3 (admissibility replay) here. Test
around current behaviour; a test that documents F1 as the present behaviour is
correct and useful. Those are separate work.

---

## 6. Verification

* `devtools::document()` --- `man/` and `NAMESPACE` regenerate cleanly.
* `devtools::test()` --- no regressions; report before/after counts.
* `devtools::check(args = c("--as-cran", "--no-build-vignettes"))` --- must stay
  0/0/0.
* `git grep -n "debias_gate\|harm_flag_debiased\|t_gate\|c_gate\|\.fs_dg_"` over
  `R/`, `tests/`, `vignettes/` and `quarto/**/*.qmd` returns **nothing**.
* `git grep -n "gate"` over `R/` returns only the calibrator tolerance gates and
  sense-6 guards --- enumerate what remains and justify each.
* Render one simulation cell and one `quarto/methodology` document that used the
  old names, confirming they still execute.

## 7. Sequencing

1. Tests on pre-rename names (§5). Commit.
2. Internal-only renames --- `.fs_dg_*` → `.fs_mr_*`, `.fs_dg_gate_null()`.
   Commit. User-invisible.
3. API renames (§1.1, §1.2) plus their `quarto/` call sites. **One commit** ---
   the package and the notebooks must not be out of step.
4. Documentation (§3) and `NEWS.md`. Commit.

Stop and report if: the calibrator or `tier1`/`tier2` false positives are
affected; any check regresses; or a rename would collide with an existing name.
