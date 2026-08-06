# Closed-form verification scripts

## What these are, and what they are not

These are the **reference side of a fixture that has not been built**. They
verify that the closed forms are *correct* — that the algebra in the derivations
gives the numbers it claims. They do **not** test the package pipeline.

Nothing here will fail if `forestsearch()` regresses. A fixture comparing the
pipeline's output against these reference values is the missing piece; until it
exists, these scripts establish only that the target values are right.

The one exception is narrow and deliberate: `verify_residual_centering.R` calls
the package's shipped `.fs_mr_ij_var()` rather than a copy of it (see
[Provenance](#provenance)), so the IJ arithmetic it exercises is the real one.
Everything feeding that function is still synthetic.

## The scripts

### `verify_closedform_fixture.R` — V1 to V6

| block | what it checks |
|---|---|
| V1 + V2 | the exact multivariate-normal law of the perturbations, and the two-candidate selection bias |
| V3 | the `E[max]` constants for M iid standard normals (`1/sqrt(pi)`, `3/(2 sqrt(pi))`, 1.02938) |
| V4 | the M = 1 degenerate IJ limit, `V_tilde / sigma2_D -> 4` |
| V5 | the closed-form discriminator between the divide-by-`B` and divide-by-winners conventions for `bias_sel` under an active admission floor |
| V6 | the same convention question in the other direction: `V_tilde` is invariant to it, `V_hat` is not |

V5 and V6 are the audit's **F13** material — the denominator convention for
`bias_sel`. That question is undecided, and these scripts do not decide it; they
establish what each convention implies, so the choice can be made on the
numbers.

### `verify_eq9_alignment.R` — V8 and `verify_conditional_convention.R` — V9

The two scripts behind the F13 change, in order.

**V8** compares three arms: as-coded (`bias_sel` over winner draws, `bias_fix`
over `B`, IJ on winners); Eq. 9 literal (both over `B` with `D := 0` on a
no-winner draw, IJ over all `B`); and half-aligned. It establishes that the
as-coded arm is not centered (`mean(r)` = −0.169) and that Eq. 9 literal is
(+2.3e−17), and that the half-aligned arm is **worse than doing nothing**
(−0.367) — so the denominator change alone is a regression, not a partial
improvement.

**V9** adds the option V8 did not consider: **conditional**, both terms over the
winner draws, IJ on the same. It is equally centered (`mean(r)` = 9.7e−18),
leaves `tilde_V` unchanged (+0.0%) and leaves `selection_bias` unchanged, where
Eq. 9 literal multiplies `selection_bias` by the selection rate.

**The package implements the conditional convention.** Both repairs are
coherent and both center; they differ in estimand, and centering does not
choose between them. See `../NOTE_f13_open_questions.md`.

### `verify_residual_centering.R` — V7

Checks the residual-centering condition underlying the finite-B correction:
**when the correction's uncentered second moment `B^-1 sum_b r_b^2` equals
Wager–Hastie–Efron's centered `B^-1 sum_b (t_b - tbar)^2`.** They coincide iff
`mean(r_b) = 0` over the draws the IJ actually uses, which holds exactly at
selection rate 1 and fails as the admission floor starts to bind.

## Running them

Base R only, seeded, no arguments:

```sh
Rscript verify_closedform_fixture.R
Rscript verify_residual_centering.R
Rscript verify_eq9_alignment.R
Rscript verify_conditional_convention.R
```

Verified under R 4.3.3. All were re-run under R 4.6.1 when filed here; all
blocks passed.

**`verify_residual_centering.R`, `verify_eq9_alignment.R` and
`verify_conditional_convention.R` additionally require `forestsearch` to be
installed**, because of the provenance change below.
`verify_closedform_fixture.R` has no package dependency.

## Provenance {#provenance}

`verify_residual_centering.R`, `verify_eq9_alignment.R` and
`verify_conditional_convention.R` all originally `source()`d
`ij_var_verbatim.R`, a
byte-for-byte snapshot of `.fs_mr_ij_var()` from `R/fs_mr_inference.R`. That
snapshot was **deliberately not filed**; each now binds the shipped function
directly:

```r
.fs_mr_ij_var <- forestsearch:::.fs_mr_ij_var   # the SHIPPED function, not a copy
```

The reasoning is in that script's own header and is not repeated here.

**Edits made when filing**, in full:

| file | edit |
|---|---|
| `verify_residual_centering.R` | the `source()` line above |
| `verify_residual_centering.R` | header corrected — it still described the function as extracted verbatim, base R only, with no package dependencies, all three of which the line above makes false |
| `verify_residual_centering.R` | closing note qualified — it cites V6 without saying V6 is in the sibling script |
| `verify_eq9_alignment.R` | same `source()` replacement, and the matching header correction |
| `verify_conditional_convention.R` | same `source()` replacement, and the matching header correction |
| `verify_closedform_fixture.R` | none; as received |

Line-range references in `verify_residual_centering.R`'s header were checked
against the source rather than assumed: `.fs_mr_ij_var()` is at
`R/fs_mr_inference.R:195-206`, and the five reproduced lines match `434-457`
verbatim. Both were, and remain, accurate.
