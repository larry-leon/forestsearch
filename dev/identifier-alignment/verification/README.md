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
```

Verified under R 4.3.3. Both were re-run under R 4.6.1 when filed here; all
blocks passed on both.

**`verify_residual_centering.R` additionally requires `forestsearch` to be
installed**, because of the provenance change below. `verify_closedform_fixture.R`
has no package dependency.

## Provenance {#provenance}

`verify_residual_centering.R` originally `source()`d `ij_var_verbatim.R`, a
byte-for-byte snapshot of `.fs_mr_ij_var()` from `R/fs_mr_inference.R`. That
snapshot was **deliberately not filed**, and the `source()` line now reads:

```r
.fs_mr_ij_var <- forestsearch:::.fs_mr_ij_var   # the SHIPPED function, not a copy
```

A second copy of a package function that can drift from the original is exactly
the failure mode this audit exists to catch — a quantity defined once and
reconstructed at a second site. Binding the shipped function directly means V7
cannot silently verify a stale computation.

This is the only edit made to either script when filing; their headers and
caveats are as received.
