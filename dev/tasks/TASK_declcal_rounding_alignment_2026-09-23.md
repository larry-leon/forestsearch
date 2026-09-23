# TASK — Align κ̂ / FŴ with the screen's rounded admission rule

**Repo:** `larry-leon/forestsearch`, branch `feature/glm-extension` — run CC from the forestsearch clone.

**Depends on** `TASK_pconsistency_digits_argument_2026-09-23.md` having landed: the calibration must be able to
read `pconsistency.digits` off a fit's `args_call_all`. Do not start until that is committed.

**Purpose.** The search screen admits on `round(Pcons, d) >= p*`, but `fs_declaration_calibration.R` computes
FŴ at the exact threshold `qnorm((1 + p*)/2)`. The screen's real bar is lower, so FŴ at the exact threshold is
not a valid ceiling for the procedure as implemented. This task puts FŴ and the implied-p\* output on the
screen's own scale.

**`R/` CALLOUT — this CHANGES BEHAVIOUR.** Reported FŴ values change (they rise), and the implied-p\* output
changes shape. It does **not** change the method: κ̂ is a quantile of the max statistic at level α and is
unchanged. What changes is the mapping between the analytic diagnostic and the screen it describes.

**Decision on 2026-09-23:** the effective-threshold version replaces the exact-threshold version outright. Do
**not** retain, report or store the exact-threshold FŴ alongside it.

**Kind:** package change. **Compute: negligible** — one GBSG fit plus the calibration. Hard abort at 30 min.

---

## 0. First action

1. Copy this file to `dev/tasks/TASK_declcal_rounding_alignment_2026-09-23.md`; `git add` that path; commit
   `docs(tasks): align kappa-hat/FW-hat with the rounded admission rule (2026-09-23)`.
2. Record HEAD, installed forestsearch version and build date, R version, platform.
3. `git status --short`: record pre-existing untracked files; never stage them.

No `git fetch` / `pull` / `push`. Explicit paths on every `git add`.

---

## 1. Baseline capture (before any edit)

On the **unmodified** package, run the GBSG application as `quarto/simulations/gbsg_app_null/run_gbsg_app_null.R`
(`e6477a22`) sets it up, at its own `p* = 0.90`, then its declaration calibration at `c0 = 0.75`, `alpha = 0.10`.

Record κ̂, FŴ, the implied p\*, and the admitted set. Save to a scratch `.rds` under `~/Downloads`, not the
repo, not committed. This is the comparator for Gates 1 and 2.

---

## 2. The effective threshold — one shared internal helper

Read the current source before editing.

Add a single internal helper that, given `p_star` and `digits`, returns the threshold on the **`Pcons` scale**
that `round(Pcons, digits) >= p_star` actually implements:

```
g          <- ceiling(p_star * 10^digits) / 10^digits   # p* rounded UP to the grid
pcons_eff  <- g - 0.5 * 10^(-digits)
```

Every site that needs this must call the helper — do not repeat the expression. The drift this prevents is real:
`quarto/simulations/gbsg_020/scripts_dinamr/declcalc0_run.R:414` already reproduces the rounded screen with its
own inline `round(rate, digits) >= p_star`.

Under `consistency_method = "resample"`, `Pcons = 2*pnorm(T) - 1`, so the z-scale threshold is
`qnorm((1 + pcons_eff)/2)`. That conversion belongs only on the resample path.

**Verify from source how the calibration behaves under `consistency_method = "split"`** and handle it
explicitly: on that path `Pcons` is `k / n_valid`, already discrete, and the closed-form `2*pnorm(T) - 1`
correspondence does not hold. If the calibration is resample-only, say so in the roxygen rather than computing
something undefined.

---

## 3. The change in `R/fs_declaration_calibration.R`

- **FŴ** — replace the exact threshold `z_pstar <- qnorm((1 + p_star)/2)` (currently around `:433`) with the
  effective one from §2. This is the only substantive change to what FŴ measures.
- **κ̂ — unchanged.** It is a quantile of the max statistic at level α; applying it compares it to T directly,
  with no rounding involved. Any change to κ̂ is a Gate 1 failure.
- **Implied p\* — becomes a settable pair.** Replace the bare `pstar_implied = 2*pnorm(kap) - 1` (currently
  around `:563`) with the information an analyst needs to actually run κ̂ as a p\* screen, at the fit's own
  `digits`:
  - the smallest settable `p*` whose effective threshold is at or above κ̂;
  - the effective threshold that setting achieves, on both the `Pcons` and z scales;
  - the gap to κ̂ in z units (positive = conservative);
  - the smallest `digits` at which that gap falls below 0.01 in z units, with the corresponding `p*`.
  - Where no `p*` at the fit's `digits` achieves κ̂ from below, report that plainly rather than returning the
    nearest value silently.
- **Reading `digits`.** Take it from the fit's `args_call_all`. Fall back to the `subgroup.consistency()`
  default of 2 when absent, matching the existing precedent at `declcalc0_run.R:316-320`, and record which
  route was used.

---

## 4. Roxygen and NEWS

- Document that FŴ is the family-wise size of the screen **as implemented**, at the rounded admission rule, and
  that it therefore depends on `pconsistency.digits`.
- Document the implied-p\* output as guidance on what to set, not as a mathematical identity.
- NEWS.md, under the current development version: FŴ and the implied-p\* output are now computed at the
  screen's rounded admission threshold; FŴ values rise relative to previous releases; κ̂ is unchanged.

---

## 5. Gates (stop on failure)

- **Gate 1 — κ̂ is untouched.** On the modified package, the GBSG calibration returns κ̂ numerically identical
  to the Step 1 baseline. Report to at least four decimals.
- **Gate 2 — FŴ moves in the right direction.** The new FŴ is strictly greater than the baseline value at the
  same `c0` and `alpha`. Report both. A new value at or below the baseline is a gate failure — the effective
  threshold is lower than the exact one, so the size must rise.
- **Gate 3 — the helper matches `round()` in fact, not in theory.** For `digits` in 2, 3, 4 and `p*` in
  {0.90, 0.95, 0.99, 0.9936}, take a dense grid of `Pcons` values straddling `pcons_eff` and assert that
  `round(Pcons, digits) >= p_star` agrees with `Pcons >= pcons_eff` at every point, or report every
  disagreement. R rounds halves to even, so exact boundary values may differ; report any such point rather
  than special-casing it.
- **Gate 4 — the settable pair round-trips.** For the GBSG fit, take the reported `(p*, digits)` pair, run
  `forestsearch()` at exactly those settings, and confirm the admitted set matches what κ̂ admits — or report
  the difference and its cause. This is the check that the guidance is usable rather than merely arithmetic.
- **Gate 5 — test suite and `R CMD check`.** `devtools::test()` passes; `R CMD check` clean with no new NOTEs,
  WARNINGs or ERRORs against a pre-change baseline. Any test asserting the old FŴ or implied-p\* values must be
  updated to the new definition, and each such update reported individually.

---

## 6. Commits

Explicit paths, in order: task doc; `R/` changes plus the new helper and roxygen; regenerated `man/`; NEWS.md;
the verification record beside the existing `REPORT_*` files.

Do not commit the Step 1 baseline `.rds`.

---

## POST-CONDITIONS (machine-checkable)

1. Gate 1: κ̂ identical to baseline to at least four decimals.
2. Gate 2: new FŴ strictly greater than baseline; both reported.
3. Gate 3: helper agrees with `round()` across the grid, or every disagreement is enumerated.
4. Gate 4: the reported `(p*, digits)` pair reproduces κ̂'s admitted set, or the difference is reported.
5. Gate 5: tests pass; `R CMD check` clean; every updated test listed.
6. No exact-threshold FŴ is retained, reported or stored anywhere in the output.
7. The effective threshold is computed in exactly one place; no site repeats the expression.
8. Files modified are confined to `R/fs_declaration_calibration.R`, the file holding the new helper, `man/`,
   `NEWS.md`, `dev/tasks/` and the verification record.

---

## OUT OF SCOPE

No change to the rounding design itself — the rounded value continues to decide admission, by decision on
2026-09-23. **No change to MR admission** (`R/fs_mr_inference.R:660-661`), which uses the same exact-threshold
bound and carries the same misalignment; that is a separate decision and is recorded, not fixed here. No
re-run of the GBSG p\* grid campaign. No change to `pconsistency.digits` defaults. Nothing written to
`fs-glms-interpretable`: re-running the applications under the new FŴ belongs to the applications chat.
