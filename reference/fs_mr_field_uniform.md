# Uniform (kappa) calibration of the field interval

Computes, at analysis time and from the trial's own influence structure,
the smallest widening factor `kappa` on `kappa_grid` such that the
two-sided field interval `c + kappa * (q - c)` attains `1 - alpha`
coverage uniformly over a winner-profile protection family: hypothetical
true fields equal to the gate's shrunk field `w` with the winner's entry
set to `max_{g != winner} w_g + delta * sigma_sel`, for `delta` on
`delta_grid` (in SE units of the winner). Nothing is tuned to any
simulation design.

## Usage

``` r
fs_mr_field_uniform(
  B = NULL,
  Sigma = NULL,
  w,
  sel,
  sigma_sel,
  p_hat,
  t_g = NULL,
  reselection = "maxeff",
  sz = NULL,
  effect_neighborhood = 0.1,
  selection_rule = "neighborhood",
  log_scale = TRUE,
  sdv = NULL,
  zcons_c = NULL,
  delta_grid = seq(0, 4, by = 0.5),
  mass = 0.99,
  M_cap = 40L,
  R_rep = 300L,
  R_out = 300L,
  R_in = 150L,
  alpha = 0.05,
  kappa_grid = seq(1, 3, by = 0.01),
  seed = NULL
)
```

## Arguments

- B:

  Optional n x K influence matrix (the gate's `asm$B`); multiplier
  perturbations are `crossprod(B, rnorm(n))`, covariance `t(B) %*% B`.

- Sigma:

  Optional K x K covariance, used when `B` is `NULL` (the reference
  checks); perturbations via its Cholesky root.

- w:

  Numeric K: the gate's shrunk field (`beta_hat` with the winner's entry
  replaced by the de-biased estimate), working scale.

- sel:

  Integer: the winner's index in `w`.

- sigma_sel:

  The winner's SE on the working scale (`sdv[sel]`).

- p_hat:

  Numeric K: re-selection frequencies from the gate's multiplier pass
  (`tabulate(winner)/draws`).

- t_g:

  Numeric K admission floors, or `NULL` for no admission filter (the
  gate's resolved `t_g`, restricted with the candidate set).

- reselection, sz, effect_neighborhood, selection_rule, log_scale:

  The gate's selection-map configuration; `zcons_c`, when non-`NULL`, is
  the consistency centre `c_cons` and `sdv` the per-candidate SEs so the
  `maxcons` score can be reconstructed on the restricted family.

- sdv:

  Numeric K per-candidate SEs (needed only for `maxcons`).

- zcons_c:

  Consistency centre `c_cons`, or `NULL`.

- delta_grid:

  Winner-separation grid in SE units (H2 default 0-4 by 0.5).

- mass, M_cap:

  Mass-carrying reduction controls (H3 as adjusted at the Gate 1
  adjudication, 2026-09-05: cap 40, so the \>= 0.99 mass target is
  reachable on enumerated forestsearch families; original default was
  12).

- R_rep, R_out, R_in:

  Monte Carlo sizes per delta (H4 defaults 300/300/150).

- alpha:

  Two-sided miscoverage target (0.05).

- kappa_grid:

  Candidate widening factors (default 1 to 3 by 0.01; the grid was
  extended past 2 at the Gate 1 adjudication after ceiling hits).

- seed:

  Optional integer; drawn once at entry.

## Value

List: `kappa` (kappa\*), `kappa_mcse` (Monte Carlo SE of kappa\* from
the binding profile's binomial error and the local slope of min-delta C2
in kappa), `M`, `mass_covered`, `keep` (indices), `minC1`, `C1` (profile
over `delta_grid`), `C2_k1` (profile at kappa = 1), `C2_kstar` (profile
at kappa\*), `kappa_grid`, `C2_min` (min-over-delta profile over
`kappa_grid`), `n_kept` (trials with a winner per delta), `delta_grid`,
`timing_seconds`.

## Details

The computation restricts to a mass-carrying candidate set (the smallest
`M <= M_cap` candidates by re-selection frequency whose cumulative share
of the re-selection mass is at least `mass`, the winner always
included), then per `delta` replays `R_rep` hypothetical trials with
Gaussian multipliers of covariance exactly `Sigma_hat` (via
`crossprod(B, xi)`, or a Cholesky root of a supplied `Sigma`), applies
the gate's own selection map (thresholded argmax under
`reselection = "maxeff"`, vectorized; any other rule via the per-draw
`.fs_mr_select` path), drops trials with no winner (conditioning on
detection, as in the real analysis), and runs the field procedure on
each hypothetical trial to obtain the coverage profiles C1(delta)
(one-sided) and C2(delta; kappa) (two-sided, widened).

The guarantee, as documented for the gate: the one-sided bound is
uniformly valid; the two-sided interval widened by the returned `kappa`
is uniformly valid over the winner-profile family; the plain quantile
interval is approximate.
