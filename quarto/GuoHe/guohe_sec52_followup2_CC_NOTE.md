---
bibliography: []
---

# Follow-up 2 to the §5.2 reproduction — two changes before the pilot

Both items must land **before** any coverage run. No coverage output exists
yet, so both are free now and expensive the moment a scenario `.rds` is
written. Truth step 0 is settled and its verdict stands; nothing below
disputes it.

Scope: `guohe_sec52_run.R` (item 1 and the item 2 cache guard),
`guohe_sec52.qmd` (item 3). `guohe_sec52_truth.R` is **supplied updated with
this note** — replace the file wholesale; do not hand-edit it.

Acknowledgement of the anchor report first: the arithmetic reproduces exactly
(combined $-0.550$, mean deviation $-0.00067$), and dividing each deviation by
its $z$ recovers $\mathrm{SE}(30)$ as $0.00302, 0.00297, 0.00291, 0.00287,
0.00280, 0.00292$ against the smoke value scaled by $1/\sqrt{10}$ of
$0.00291$. The variance fell exactly as $n_{big}^{-1/2}$ across all six
scenarios, which is independent corroboration that the truth machinery
behaves as intended. The fresh-noise adjudication is correct and no gate
tolerance was or should be widened.

---

## 1. Scenario seed bases are closer together than the replicate count

In `guohe_sec52_run.R`:

```r
base <- 1000000L + as.integer(sum(utf8ToInt(id)) * 997L)
```

Adjacent scenario ids differ by one character, so adjacent bases differ by
**997**, while replicates span `base + 1 … base + 2000`. Adjacent scenarios
therefore share **1003 of 2000** replicate seeds, and scenarios two apart
share 6:

| pair | 0.0–0.1 | 0.1–0.2 | 0.2–0.3 | 0.3–0.4 | 0.4–0.5 | two apart |
|---|---|---|---|---|---|---|
| shared seeds | 1003 | 1003 | 1003 | 1003 | 1003 | 6 |

Consequence: adjacent $\beta_2$ scenarios run on *coupled* datasets — identical
$W$, identical $D$, and identical uniforms driving the exponential, differing
only in rate. No individual cell is biased; each scenario's coverage remains
unbiased. What it affects is the joint reading — the six rows are correlated,
so the "cells beyond 2 combined SE" tally is not a count over independent
trials, and the combined-SE construction tacitly assumes cell independence.

**Replacement:**

```r
  # Scenario-specific seed base. The multiplier must EXCEED n_rep: adjacent
  # ids differ by one character, so the spacing between bases is the
  # multiplier itself, and a multiplier below n_rep makes adjacent scenarios
  # share replicate seeds (at 997 with 2000 replicates, 1003 of 2000 were
  # shared, coupling the datasets across beta2 and correlating the rows).
  base <- 1000000L + as.integer(sum(utf8ToInt(id))) * 100003L
```

Verified: spacing 100003, zero shared seeds at every pair, maximum
`base + n_rep` = 93,404,772, well inside `.Machine$integer.max`. Note the
parenthesis moves — `as.integer()` wraps only the character sum, so the
multiplication is integer-by-integer.

**This is inherited from the §5.1 driver**, which uses the identical
construction with the same 997 spacing and the same 2000 replicates. §5.1 is
already run and committed; do not retrofit or re-run it. Record the property
in the companion-paper hand-off as a known, benign feature of the §5.1
results, and note in the §5.2 report that §5.2 does not share it.

---

## 2. Exact truth curve at $\beta_2 = 0$

**Why.** At $\beta_2 = 0$ the DGM sets $b(w) = 0$ for every $w$, so within
*every* $S(c)$ the treatment-only Cox model is correctly specified with true
coefficient exactly zero: $\beta(c) \equiv 0$ on the whole grid. The
production cache instead scores against a simulated estimate whose anchor
deviation is $-0.0057$ — the **largest of the six in $\beta$ units**, in the
one scenario where the error is entirely avoidable. Because a single curve
scores all 2000 replicates, that offset is systematic and does not average
out. Propagating it through
$\mathrm{d}\,\text{cover}/\mathrm{d}\delta = \varphi(1.645)/\sigma_D$:

| $\sigma_D$ | coverage shift | vs MCSE (0.0049) | vs 2-combined-SE band (0.0138) |
|---|---|---|---|
| 0.15 | $-0.0039$ | 0.80× | 28% |
| 0.20 | $-0.0029$ | 0.60× | 21% |
| 0.30 | $-0.0020$ | 0.40× | 14% |

The sign matters: scoring $\gamma_s$ below truth *understates* coverage, in
the row where the published proposed columns sit at 0.947–0.962 and where the
conservatism check is the designated diagnostic. A systematic $-0.003$ there
could turn a genuinely conservative result into a near-nominal one — which the
brief designates as the red flag for a mis-specified family. That is a
false-alarm risk on the exact test the reproduction exists to run.

**What changed in `guohe_sec52_truth.R`** (supplied; replace wholesale):

- New field `beta_exact`: `rep(0, m)` when $\beta_2 = 0$, `NULL` otherwise.
  No closed form exists for $\beta_2 > 0$ — the misspecified Lin–Wei estimand
  is a weighted average of $b(\cdot)$ over $[0, c]$ with censoring-dependent
  weights.
- New fields `scoring_basis` (`"exact"` / `"smooth"`), `offset_anchor`
  ($\hat\beta(30) - \beta_2$, observable for every scenario), and
  `offset_mean` (curve-wide, where an exact curve exists).
- `gh52_truth_at(truth, c_hat, use = c("auto", "exact", "smooth", "raw"))`.
  `"auto"` is the new default and is what the harness should use: exact where
  available, isotonic otherwise. `"exact"` errors where no closed form exists.
  `"smooth"`/`"raw"` force the simulated curve, retained for diagnostics.
- **The gates are unchanged and still run on the simulated curve.** They
  remain a check on the DGM and the prefix sweep; only the *scoring* basis
  changes.
- Backward compatible: objects written before `beta_exact` existed return
  `NULL` for it, so `"auto"` falls back to `"smooth"` and old caches stay
  readable.

Verified by execution: at $\beta_2 = 0$, `"auto"` returns exactly zero at
every $\hat c$ while `"smooth"` returns $-0.014$; at $\beta_2 = 0.3$,
`"auto"` and `"smooth"` agree identically and `"exact"` errors; a
`beta_exact`-stripped object falls back correctly; all gates still pass.

**Actions:**

1. Replace `guohe_sec52_truth.R` with the supplied file.
2. Recompute the $\beta_2 = 0$ cache only:
   `Rscript quarto/GuoHe/guohe_sec52_truth.R --force --beta2=0`. The seed is
   unchanged, so the simulated curve regenerates bit-for-bit and only the new
   fields are added (~5.4 min). **Do not recompute the other five** — they are
   gate-cleared and unaffected.
3. Add a cache guard in `guohe_sec52_run.R`, beside the existing
   `all(tr$gates$pass)` check, so a stale cache fails loudly rather than
   silently reverting to offset scoring:

```r
  if (isTRUE(tr$beta2 == 0) && is.null(tr$beta_exact)) {
    stop("Truth cache for beta2 = 0 predates the exact-curve field ",
         "(beta_exact). Recompute it: Rscript quarto/GuoHe/guohe_sec52_truth.R ",
         "--force --beta2=0", call. = FALSE)
  }
```

4. No change is needed in `guohe_sec52_sim.R`: it calls `gh52_truth_at()` at
   its default, which is now `"auto"`. Confirm by reading that the calls in
   `gh52_naive()` and `gh52_one_rep()` pass no `use =` argument.

---

## 3. Ledger and report additions

In `guohe_sec52.qmd`:

- Amend inferred item 3 so the isotonic projection is described as the basis
  **where no closed form exists**, and add that $\beta_2 = 0$ is scored
  against the exact $\beta(c) \equiv 0$ and therefore carries no truth-curve
  Monte Carlo error.
- Add a new inferred item recording the **residual truth-curve offsets** for
  the five $\beta_2 > 0$ scenarios — read `offset_anchor` from each cache
  (production values $+0.0030, +0.0017, -0.0008, -0.0012, -0.0010$) — with the
  note that these are fixed per scenario, do not average out over replicates,
  and imply a systematic coverage shift of order 0.001–0.002, roughly a tenth
  of the comparison band.
- In the truth-curve section, state which basis each scenario uses (the new
  `scoring_basis` field) so the asymmetry is visible rather than implicit.

---

## 4. Considered and NOT implemented — raise if you disagree

$\beta(30) = \beta_2$ holds exactly for *every* scenario, not only
$\beta_2 = 0$, so the anchor deviation $\varepsilon(30) = \hat\beta(30) -
\beta_2$ is observable throughout and the whole curve could be shifted by it:
$\tilde\beta(c) = \hat\beta(c) - \varepsilon(30)$. Because the grid points are
nested prefixes of one sample, $\varepsilon(c)$ is strongly correlated with
$\varepsilon(30)$, so this removes much of the offset. A rough calculation
puts the break-even correlation at $c = 60$ near 0.71 against an
information-weighting heuristic that gives almost exactly that, so the
adjustment is variance-neutral at the far end and variance-reducing near
$c = 30$, where selection concentrates as $\beta_2$ grows.

It is **not implemented**: it is a further design change beyond what was
approved, it would require recomputing or patching all five remaining caches,
and the residual offsets it targets are already an order of magnitude below
the comparison band. Recorded here so the option is on file rather than
rediscovered later.

---

## 5. Sequence after these changes

1. Replace the truth file; recompute the $\beta_2 = 0$ cache only.
2. Apply the seed-base replacement and the cache guard in the driver.
3. Apply the ledger and report additions.
4. `Rscript quarto/GuoHe/guohe_sec52_run.R --pilot --B=2000 --cores=120`, and
   report the projection before production, as already agreed.

Conventional Commits on `feature/glm-extension`; the owner commits and pushes.
Report failures loudly; never fill a missing number with a plausible one.
