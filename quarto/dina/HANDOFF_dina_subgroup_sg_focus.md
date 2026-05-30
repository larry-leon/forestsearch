# Handoff Note — Add `sg_focus` / `effect_neighborhood` selection modes to `dina_subgroup()`

**Purpose of this note:** seed a new chat to implement an enhancement to
`forestsearch::dina_subgroup()`. Paste this into the first message (and
attach the source files listed under "Files to attach"). It captures the
motivation, the design, the implementation plan, and the validation
path, so the new chat can build on the real code rather than
reconstructing context.

---

## 1. One-line task

Add a `sg_focus` argument (`"maxSG"`, `"minSG"`, `"effMaxSG"`,
`"effMinSG"`, `"eff"`) and an `effect_neighborhood` argument to
`dina_subgroup()`, mirroring `forestsearch()`'s selection semantics,
**defaulting to the current behavior** (additive, backward-compatible).

---

## 2. Why — the finding that motivates this

A simulation study comparing `forestsearch(sg_focus = "effMaxSG")` against
`dina_subgroup()` on a single-factor threshold DGM (Cox; true harm
subgroup `S = {x1 >= 0.30}`, matched harm thresholds
`effect.threshold = 1.25` for forestsearch and `m_diff = log(1.25) ~
0.223` for dina) showed **DINA systematically identifies much larger
subgroups than forestsearch**, with near-perfect sensitivity but very low
specificity (e.g. discovered `|S-hat|/n` approaching 1, specificity ~0.0–0.1).

This is **correct behavior, not a bug** — and the reason exposes the gap:

- **`dina_subgroup()` uses a fixed hard threshold.** It keeps every
  candidate whose subgroup-mean `tau_hat >= m_diff` and returns the
  largest. With `m_diff` well below the in-subgroup effect `tau_in`, the
  qualifying set balloons: the mean over `{x1 >= q}` equals
  `tau_in * P(x1 >= 0.30 | x1 >= q)`, which clears a low `m_diff` even at
  `q = -1` (the whole sample). So "largest qualifying" drifts to nearly
  everyone.
- **`forestsearch`'s `effMaxSG` uses a relative neighborhood.** It first
  restricts to candidates whose effect is within `effect_neighborhood`
  (default 0.10) of the **maximum observed effect** —
  `effect >= (1 - effect_neighborhood) * max(effect)` — then takes the
  largest of *those*. Because the band is anchored to the data's
  strongest signal rather than a fixed cutoff, it does not drift to the
  whole sample. **This neighborhood mechanism is most of why forestsearch
  returns the tighter subgroup, and it is exactly what DINA lacks.**

Adding the `sg_focus` family (with `effect_neighborhood`) to
`dina_subgroup()` closes this gap and makes the two methods genuinely
apples-to-apples (same selection rule, same band), isolating the only
remaining difference as the search space: DINA is univariate; forestsearch
searches up to `maxk` factors.

---

## 3. Design (proposed; confirm in the new chat)

`dina_subgroup()` already enumerates every `(covariate, direction,
threshold)` candidate with its subgroup-mean `tau_hat` and size `n_S`
(only candidates with `n_S >= n_min` and `mean tau_hat >= m_diff` are
currently eligible). The new `sg_focus` changes only the **selection step
over that enumerated, eligible set** — no change to the search itself.

Proposed modes (natural-parameter / harm scale; larger `tau_hat` = more
harm):

- **`"maxSG"`** — largest sample size among eligible candidates. **This is
  the current behavior; make it the default** so nothing changes unless
  the user opts in. (Today's tie-break: larger `tau_hat`.)
- **`"minSG"`** — smallest eligible subgroup (most concentrated).
- **`"effMaxSG"`** — among candidates within `effect_neighborhood` of the
  maximum `tau_hat`, return the **largest**. (Direct analog of
  forestsearch; the over-extension fix.)
- **`"effMinSG"`** — same band, return the **smallest**.
- **`"eff"`** — the single candidate with the largest `tau_hat` (subject
  to `n_min` / `m_diff`); tightest.

New companion arg: **`effect_neighborhood`** (numeric in `[0, 1)`, default
`0.10`), used only by `effMaxSG`/`effMinSG`. Inclusion test on the
natural scale: a candidate is in-band iff
`tau_hat >= (1 - effect_neighborhood) * max(tau_hat over eligible)`.
**Mind the scale:** forestsearch defines the band on the **natural
(ratio)** scale for ratio measures, not the log scale. DINA's `tau_hat`
is on the log/natural-parameter scale (log-HR for Cox), so decide and
document explicitly whether the band is applied to `tau_hat` (log) or
`exp(tau_hat)` (ratio) — to match forestsearch for Cox it should be the
**ratio** scale, i.e. band on `exp(tau_hat)`. This is the one subtle
correctness point; get it right and document it.

Backward-compatibility / validation rules:
- `sg_focus = "maxSG"` (default) must reproduce current results exactly.
- Accept the forestsearch GLM-vocabulary aliases too if cheap
  (`effMaxSG`/`effMinSG`/`eff` already read naturally; `maxSG`/`minSG`
  unchanged). See `.normalize_sg_focus()` in `forestsearch_helpers.R`.
- Loud errors on invalid `sg_focus` / out-of-range `effect_neighborhood`,
  per the package convention (`stop()` not silent).

---

## 4. Implementation pointers (real code to build on)

- **`R/dina_subgroup.R`** — the function. The search loop builds the
  candidate set; the qualifying test is literally `if (mean_tau < m_diff)
  next` and the current "largest, tie-break larger tau" selection is the
  `(n_S > best$n_subgroup) || (n_S == best$n_subgroup && mean_tau >
  best$mean_tau)` comparison. The refactor: **collect all eligible
  candidates into a table first, then apply an `sg_focus`-driven selection
  step** (rather than the running-best single pass), so the neighborhood
  band can be computed over the full eligible set.
- **`R/subgroup_consistency_helpers.R`** — `sort_subgroups()` (line ~468)
  and `sort_subgroups_preview()` (line ~542) implement forestsearch's
  selection/sort keys and the `effect_neighborhood` band + Pareto options.
  **Reuse or closely mirror this logic** so the two functions stay
  consistent; do not reinvent the band/tiebreak rules.
- **`R/forestsearch_helpers.R`** — `.normalize_sg_focus()` (the
  `effMaxSG -> hrMaxSG` alias map) and surrounding selection helpers.
- **`R/subgroup_consistency_main.R`** — `effect_neighborhood` default
  `0.10` and its inclusion semantics
  (`(1 - effect_neighborhood) * max(effect)`); the doc text there is the
  authoritative description to mirror in `dina_subgroup()`'s roxygen.

Also update: roxygen (`@param sg_focus`, `@param effect_neighborhood`,
`@return` if a selection-diagnostic field is added), `NEWS.md`, and the
`@seealso` cross-links. The package has **roxygen markdown ON** — write
literal `%` `<` `>` `&` (never `\%`), and keep `@section` titles plain (no
`\code{}`/backticks in the title); hand-escaping breaks Rd generation
under markdown.

---

## 5. Validation path

The comparison document **`dina_vs_forestsearch_signature_comparison.qmd`**
(built in the originating chat) is the natural validation artifact:

1. After adding the modes, re-run it with **matched `effMaxSG` on both
   sides** and the **same `effect_neighborhood`**. Expectation: DINA's
   discovered `|S-hat|/n` should drop from ~1 toward the true ~0.35, and
   specificity should rise substantially — converging toward
   forestsearch when `fs_maxk = 1` (both single-factor).
2. Confirm `sg_focus = "maxSG"` reproduces the pre-change DINA numbers
   (backward-compat check).
3. The standalone **`dina_subgroup_signature_simulation.qmd`** (DINA-only,
   with the `m_diff` sweep) is the other relevant artifact: the new
   `sg_focus` should interact with `m_diff` in the documented way (a low
   `m_diff` + `effMaxSG` should behave very differently from a low
   `m_diff` + `maxSG`).

Smoke-test convention for this codebase: living Quarto docs exercising the
full pipeline, `.tc()` accumulator + gt pass/fail table; no testthat
scaffolding. Add `sg_focus` cases to the existing
`dina_subgroup_refit_smoke_test.qmd` pattern (or a small dedicated
selection smoke test).

---

## 6. Files to attach to the new chat

Attach the **entire R codebase as a single zip** (e.g. `R.zip` containing
the `R/` directory) rather than individual files — it's less tedious and
guarantees the new chat has the complete, current, mutually-consistent
source. Always prefer the attached code over memory as the source of
truth.

Within that zip, the files most relevant to this task (where the new chat
should look first) are:

- `R/dina_subgroup.R` — the function being changed
- `R/subgroup_consistency_helpers.R` — `sort_subgroups()` /
  `sort_subgroups_preview()` (selection/sort + `effect_neighborhood` band
  to reuse or mirror)
- `R/forestsearch_helpers.R` — `.normalize_sg_focus()` + selection helpers
- `R/subgroup_consistency_main.R` — the authoritative `effect_neighborhood`
  semantics to mirror in the roxygen

Optionally also attach the validation Quarto docs (not in the codebase
zip) if you want them on hand:
`dina_vs_forestsearch_signature_comparison.qmd` and
`dina_subgroup_signature_simulation.qmd`.

---

## 7. Caveats / things not to over-read

- The neighborhood band is anchored to the **maximum observed** effect,
  which is noisy in finite samples, so `effMaxSG` inherits a mild
  winner's-curse flavor in where it anchors — fine for discovery, not to
  be read as an exact cut.
- This is a **discovery** enhancement (which subgroup is selected). It does
  **not** touch post-selection inference — the conditional-on-signature
  CIs in `dina_subgroup_bootstrap()` and the deferred selection-adjusted /
  bias-corrected bootstrap (Leon et al. 2024, Section 3.2) are separate.
- Keep scope tight: this is a selection-step addition to `dina_subgroup()`.
  Do not expand DINA's univariate search to multi-factor — the
  univariate-vs-`maxk` difference is intentional and is what the
  comparison study isolates.

---

## 8. Suggested first message for the new chat

> Add `sg_focus` (maxSG/minSG/effMaxSG/effMinSG/eff) and
> `effect_neighborhood` to `dina_subgroup()`, mirroring forestsearch's
> selection semantics, defaulting to current behavior (`maxSG`).
> Motivation: in the DINA-vs-forestsearch comparison, DINA's fixed-`m_diff`
> "largest qualifying" rule over-extends the subgroup (specificity ~0)
> whereas forestsearch's `effMaxSG` uses an `effect_neighborhood` band
> anchored to the max effect, giving a tighter subgroup. I want the same
> band option for `dina_subgroup()`. I've attached the full R codebase as a
> zip; the key files are `dina_subgroup.R`, `subgroup_consistency_helpers.R`
> (`sort_subgroups()`), `forestsearch_helpers.R` (`.normalize_sg_focus()`),
> and `subgroup_consistency_main.R` (`effect_neighborhood` semantics).
> Please unzip and audit the selection code first and propose the design
> before implementing, reusing `sort_subgroups()` logic where possible.
> Build on the attached codebase as the source of truth.
