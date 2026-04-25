# Appendix Patch Instructions

How to apply the orientation-guide integrations + pkgdown URL hyperlinks
to your locally edited copy of `actg175_binary_benefit_simulations.qmd`,
without overwriting your edits.

These instructions assume your file currently has the Appendix structure
from the prior round (A.1 prose, A.2 Reference table, A.3 Worked
example).  If your file is in a different state, the anchor strings in
each step should still locate the right insertion / replacement points.

There are **5 changes** to apply, in order.  Each change identifies a
search anchor in your file and gives you the exact text to paste.

---

## Change 1 of 5: Add italics intro to the Appendix

**Where:** Right after the top-level Appendix heading, before A.1.

**Find this anchor in your file** (the Appendix heading):

```
# Appendix: Effect Measures and Outcome Direction Conventions {#sec-appendix-conventions}
```

**Insert the following block** as a new paragraph between that heading
and `## A.1`:

```
*This appendix is a quick-lookup reference for this vignette's
directional configuration.  For comprehensive coverage across all
endpoint types -- including the eight-scenario decision table,
worked code examples for binary, continuous, count, and survival
outcomes, common pitfalls, and verification tests -- see the
standalone* [*ForestSearch Outcome Orientation Reference*](https://larry-leon.github.io/forestsearch/articles/forestsearch_outcome_orientation_guide.html) *vignette.*
```

(One blank line above and below.)

---

## Change 2 of 5: Rewrite A.1 prose

**Where:** Replace the entire body of `## A.1 The unified
harm-search convention` (everything from that heading up to but not
including `## A.2 Reference table`).

**Replace with:**

```
## A.1 The unified harm-search convention

Every subgroup search in `forestsearch` uses a **single comparison
direction** throughout the pipeline.  In each candidate subgroup
$S$, the package checks whether the estimated treatment effect
satisfies

$$
\hat{\theta}(S) \;\geq\; c_1
$$

(`effect.threshold`), with the same direction reused for the
bootstrap consistency check, $\hat{\theta}(S_b) \geq c_2$.  This
comparison never flips.  The entire question of how to set up a
search reduces to: **how do I orient my data so the subgroup I
want to find produces $\hat{\theta}(S) > 0$ on the chosen effect
scale?**

Two transformations handle the orientations that don't satisfy
this naturally.  First, **benefit search** -- finding a subgroup
in which the original experimental treatment is *beneficial* --
is implemented by switching the treatment labels so the worse
arm becomes `treat = 1`.  A benefit subgroup of the original
experimental then appears as a harm subgroup of the
switched-treatment; magnitudes are inverted at reporting time
via `subgroup_notation = "benefit"`.  Second, when the user's
binary outcome is encoded on the positive scale (higher $Y$ =
better), the binary effect estimator can substitute $1 - Y$
before computing the OR -- this is what the `adverse_outcome`
flag controls in the FS estimator path.

The flag's *name* suggests an interpretive label ("the outcome
is adverse"), but mechanically what it does in the binary FS
estimator is **flip $Y$ to $1 - Y$ when `adverse_outcome = FALSE`,
and leave $Y$ as-is when `TRUE`** (the default for binary).  The
two framings -- pre-flip $Y$ at data prep and pass `TRUE`, or
keep $Y$ positive and pass `FALSE` -- are mechanically
equivalent for binary outcomes.  For **continuous** and **count**
outcomes the FS effect estimator does **not** auto-flip on the
flag (only the GRF candidate-cut step does); for those outcome
types you must encode $Y$ on the harm direction at data prep
and pass `adverse_outcome = TRUE`.  The simplest safe rule
across all GLM outcome types is therefore: **encode $Y$ on the
adverse scale at data prep and pass `adverse_outcome = TRUE`
everywhere**.  This is the convention used by the present
vignette; see @sec-stat-goals for the specific row that applies
and Section 1.1 for the data-prep step.

::: {.callout-tip title="DGM golden rule"}
Once the data is correctly oriented, build the DGM with **positive
`k_inter`**, a **positive calibration target**, and a **positive
`k_inter_range` grid**.  The DGM operates on the final, oriented
data; outcome and treatment manipulations (negation, switching)
happen at the data-prep step, *before* the DGM is constructed.
:::
```

---

## Change 3 of 5: Insert new A.2 (Effect scales) and renumber

**Where:** Between the (new) end of A.1 and the existing
`## A.2 Reference table` heading.

**Step 3a -- insert the new section** between the end of A.1's
DGM-golden-rule callout (`:::` close fence) and the line
`## A.2 Reference table`:

```
## A.2 Effect scales and threshold conventions

`forestsearch` works on the **log scale** internally for ratio
measures and on the **identity scale** for additive measures.
This affects how `effect.threshold` is interpreted:

| Measure type | Examples | What `effect.threshold = c₁` becomes |
|--------------|----------|---------------------------------------|
| Ratio | HR, OR, IRR, RR | $\log(c_1)$ -- e.g., 1.25 → 0.223 |
| Additive (bounded) | RD, IRD | $c_1$ directly; bounded in $[-1, 1]$ |
| Additive (unbounded) | MD | $c_1$ directly; unbounded |

For RD and IRD, the package warns and remaps if the threshold
is set outside the bounded range.  For MD (mean difference) the
threshold is genuinely unbounded -- values like
30 cells/$\mu$L are legitimate, and the package no longer
remaps them as it once did when MD was incorrectly treated as a
bounded measure.  In all cases, the comparison the package
performs internally is $\hat{\theta} \geq$ (transformed
threshold), with $\hat{\theta}$ on the same scale (log for
ratios, identity for additives).
```

**Step 3b -- renumber existing A.2 → A.3:**

Change the line:

```
## A.2 Reference table
```

to:

```
## A.3 Reference table
```

**Step 3c -- renumber existing A.3 → A.4:**

Change the line:

```
## A.3 Worked example: this vignette's configuration
```

to:

```
## A.4 Worked example: this vignette's configuration
```

---

## Change 4 of 5: Add A.5 Common pitfalls

**Where:** After the end of A.4 (Worked example) and before the
top-level heading `# Session Info {#sec-session}`.

**Find this anchor:**

```
# Session Info {#sec-session}
```

**Insert immediately above** that line:

```
## A.5 Common pitfalls

The standalone [*ForestSearch Outcome Orientation Reference*](https://larry-leon.github.io/forestsearch/articles/forestsearch_outcome_orientation_guide.html)
is the authoritative source for setup mechanics across all endpoint
types.  Two pitfalls show up often enough in practice to flag
here:

::: {.callout-warning title="Pitfall 1: calibration target on the wrong scale"}
Setting (e.g.) `cal_target_or = 0.5` together with a negative
`k_inter_range` produces a $Q$ that is only marginally more
extreme than the rest of the sample in the same direction --
no real contrast.  *Symptom:* near-100% detection paired with
very high false-positive rate (~97%).  *Fix:* use a target
above 1 with a positive `k_inter` grid for ratio measures, or
a positive target with a positive grid for MD.
:::

::: {.callout-warning title="Pitfall 2: `adverse_outcome` confused with outcome negation"}
Setting `adverse_outcome = FALSE` for a benefit search with
switched treatment and a *positive* binary outcome causes the
binary FS estimator to flip $Y \to 1 - Y$, reversing the OR
direction and breaking the $\hat{\theta} \geq c_1$ comparison.
*Symptom:* zero detection.  *Fix:* leave `adverse_outcome` at
its default `TRUE` for binary, *or* pre-flip $Y$ at data prep
(the convention used by this vignette -- see @sec-stat-goals).
:::

For the full list of pitfalls -- including MD-threshold remapping
behaviour, survival with positive events, and DGM construction
gotchas -- the eight-scenario decision table, the four
verification tests, and worked code examples for binary,
continuous, count, and survival outcomes, see the standalone
[*ForestSearch Outcome Orientation Reference*](https://larry-leon.github.io/forestsearch/articles/forestsearch_outcome_orientation_guide.html)
vignette.

```

(Leave one blank line between the closing of A.5 and the
`# Session Info` heading.)

---

## Change 5 of 5: Update setup-chunk cross-reference

**Where:** Inside the setup chunk, near the calibration-target
comment block.  Look for the line that reads:

```
# magnitude is the same -- both give 2.0; see Appendix A.2.)
```

**Change `A.2` to `A.3`:**

```
# magnitude is the same -- both give 2.0; see Appendix A.3.)
```

(This is the only narrative cross-reference that drifts when A.2
becomes Reference table → A.3.  The Quarto `@sec-` anchors
elsewhere in the file are unaffected by the A.x renumbering
because they reference top-level anchors, not subsection
numbers.)

---

## Verification after applying all 5 changes

After applying, your Appendix's heading sequence should be:

```
# Appendix: Effect Measures and Outcome Direction Conventions {#sec-appendix-conventions}
*<italics intro paragraph>*
## A.1 The unified harm-search convention
## A.2 Effect scales and threshold conventions
## A.3 Reference table
## A.4 Worked example: this vignette's configuration
## A.5 Common pitfalls
# Session Info {#sec-session}
```

Two quick sanity checks before rendering:

1. **Triple-backtick balance.**  Run
   `grep -c '^```' your_file.qmd` -- the count must be even (the
   number of fenced code chunks doubled).  The new content adds
   no R chunks, so this should match what you had before plus
   zero.

2. **A.x cross-references.**  Run
   `grep -n "Appendix A\\.[0-9]"` -- you should see exactly one
   match, the updated `A.3` reference at the top of the setup
   chunk.

If you'd prefer, you can diff your current file against the staged
copy in the artifact panel from my prior turn -- the staged copy
has all 5 changes applied -- and pull just the Appendix portion
across.
