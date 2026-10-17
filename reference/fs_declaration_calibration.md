# Calibrated declaration threshold and family-wise size of the p-star screen

Post-hoc, opt-in, reported diagnostics of the consistency screen,
computed from multiplier draws the package already produces. The screen
admits a candidate `g` when its consistency rate, rounded to
`pconsistency.digits`, is at least `p_star`; on the resample path that
is exactly the standardized statistic
`T(g) = (beta_hat(g) - c_cons) / sigma_D(g)` clearing the effective
cutoff `z_pstar = qnorm((1 + pcons_eff) / 2)` (see the Rounded admission
rule section). This function reports how often the maximum of that
statistic over the candidate family would clear the screen's cutoff
under the null perturbation law (`fw_size`), and the cutoff that would
hold the family-wise declaration rate at `alpha` (`kappa_hat`).

## Usage

``` r
fs_declaration_calibration(
  fit,
  alpha = 0.05,
  family = c("prereduction", "reduced"),
  ...,
  c0 = NULL
)

# S3 method for class 'fs_declaration_calibration'
print(x, ...)
```

## Arguments

- fit:

  A `forestsearch` object fitted with the declaration field retained, or
  the list returned by
  [`fs_mr_inference()`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)
  with `keep_declaration_field = TRUE`.

- alpha:

  Target family-wise declaration rate, in (0, 1). Default `0.05`.

- family:

  `"prereduction"` (default) or `"reduced"`; see the Family section.

- ...:

  Only `p_star` and `c_cons` are read, and only when the fit's resolved
  admission set carries no consistency floor (for example
  `sg_focus = "maxeff"`); they then name the hypothetical screen whose
  size is evaluated. `c_cons` is on the comparison scale (log for ratio
  measures).

- c0:

  `NULL` (default) or a numeric vector of clinically specified,
  pre-specified protected null levels; see the Protected null level
  section. Named only (it follows `...`).

- x:

  An `fs_declaration_calibration` object.

## Value

An object of class `fs_declaration_calibration`: a list with
`kappa_hat`, `fw_size`, `alpha`, `p_star`, `digits`, `digits_source`,
`pcons_eff`, `z_pstar` (the effective z cutoff of the rounded screen),
`consistency_method`, `consistency_method_source`, `c_cons`, `c_screen`,
`B`, `multiplier_law`, `quantile_type`, `family_source`, `family_label`,
`n_family_prereduction`, `n_family_reduced`, `admitted_current` (the
candidates the executed screen admitted; `NULL` when the fit is a bare
[`fs_mr_inference()`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)
result), `admitted_calibrated`, `admitted_pstar` (the relabelled current
rule over the same family), `screened` (the candidates the consistency
screen evaluated, when known), `Mstar`, `beta_hat`, `sigma_D`, `T_hat`,
`column_sd`, `zstar_mean`, `field_cor` (families of at most 8
candidates, else `NULL`), and a `reduction` list recording the replay of
the near-duplicate reduction. When `c0` is given it also carries `c0`, a
list with `table` (one row per `c0`: `c0`, `c0_cmp`, `kappa_hat`,
`fw_size`, `pstar_settable`, `pstar_achievable`, `pcons_eff_settable`,
`z_eff_settable`, `z_gap`, `digits_fine`, `pstar_fine`, `z_gap_fine`,
`n_admitted_calibrated`, the 0.90 / 0.95 / 0.99 quantiles of `Mstar_c0`,
and `is_c2`), `admitted_calibrated` (a list indexed by `c0`), `Mstar_c0`
and `source` (`"capture"` or `"field_matrix"`). Every other element is
unchanged by `c0`.

## Details

It reads a fitted object and returns a new one. It does not modify its
input, does not re-run the search, does not call the consistency engine,
and does not change what the search admitted: `admitted_calibrated` is
the set the calibrated rule *would* admit, reported beside
`admitted_current`.

## Definitions

With shared multipliers `xi[b, i]` (one vector per draw `b`, reused
across every candidate) and the dfbeta influence `db[g, i]`:

- `Zstar[b, g] = sum_i xi[b, i] * db[g, i] / sigma_D(g)`, with
  `sigma_D(g)^2 = sum_i db[g, i]^2` – the robust scale the screen itself
  uses, never a model-based standard error.

- `Mstar[b] = max_g Zstar[b, g]` over the family: one-sided, in the harm
  direction.

- `kappa_hat(alpha)` is the empirical (`type = 1`) `1 - alpha` quantile
  of `Mstar`.

- `fw_size = mean(Mstar > z_pstar)`,
  `z_pstar = qnorm((1 + pcons_eff) / 2)`. Its threshold is the fit's own
  `p_star` at the fit's own `pconsistency.digits`, not `alpha`: it is
  the family-wise size of the screen as implemented, not of the
  calibrated rule.

- The calibrated admission rule is
  `beta_hat(g) >= max(c_screen, c_cons + kappa_hat * sigma_D(g))`, the
  current rule with `z_pstar` replaced by `kappa_hat`.

## Rounded admission rule

The screen admits on `round(Pcons, digits) >= p_star`, with `digits` the
fit's `pconsistency.digits`. With `p_star` rounded up to the
`10^-digits` grid (`g`), that is
`Pcons >= pcons_eff = g - 0.5 * 10^-digits`, a bar below `p_star` itself
(0.895 for `p_star = 0.90`, `digits = 2`). On the resample path
`Pcons = 2 * pnorm(T) - 1`, so the screen is
`T >= z_pstar = qnorm((1 + pcons_eff) / 2)`. `fw_size`, `admitted_pstar`
and the settable `p*` are computed at that effective threshold;
`fw_size` is therefore the family-wise size of the screen as implemented
and depends on `pconsistency.digits`. `kappa_hat` does not: it is a
quantile of the maximum statistic and is compared with `T` directly,
with no rounding.

`digits` is read from the fit's `args_call_all$pconsistency.digits`;
when absent (a bare
[`fs_mr_inference()`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)
result, or an older fit) the
[`subgroup.consistency()`](https://larry-leon.github.io/forestsearch/reference/subgroup.consistency.md)
default of 2 is used. `digits_source` records which.

The correspondence with `T` holds only under
`consistency_method = "resample"`. Under `"split"`, `Pcons` is a split
proportion `k / n_valid` and admission is not a threshold on `T`, so
`fw_size`, `z_pstar`, `admitted_pstar` and the settable-`p*` columns are
`NA` / `NULL` rather than computed; `kappa_hat` and the calibrated rule
are still returned. The method is read from
`args_call_all$consistency_method`; a bare
[`fs_mr_inference()`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)
result is treated as resample (the field is the closed-form statistic),
recorded in `consistency_method_source`.

## Family

`family = "prereduction"` (the default, and the only value meant for the
calibrated rule) takes the maximum over the multiplier-resampling
family: every enumerated candidate of at most `maxk` factors meeting the
size minimum, before any reduction keyed on fitted quantities. The
near-duplicate reduction keys on sample-fitted summaries, so a family
reduced by it is outcome-dependent; the pre-reduction family is
covariate-measurable, and is also the conservative choice (a larger
family raises the maximum). The per-arm event minima are not replayed in
that family, which makes it a superset of the one the search chose
among.

`family = "reduced"` removes the candidates the near-duplicate reduction
removed and is a diagnostic only, so the gap is measurable. Its result
is conditional on the realized family and is labelled so; it is never
meant for admission. It needs a `forestsearch` fit (the reduction lives
in the identifier) and the full field matrix
(`keep_field_matrix = TRUE`).

## Protected null level

`c1` and `c2` state the claim; `c0` states what the claim is protected
against: a clinically specified, pre-specified benefit level (for
example HR 0.75). The family-wise declaration rate is then controlled at
`alpha` whenever every candidate's true effect is at least as good as
`c0`, rather than only when every candidate sits at `c2`. With
`delta_g = (c_cons - c0_cmp) / sigma_D(g)`:

- `Mstar_c0[b] = max_g { Zstar[b, g] - delta_g }`;

- `kappa_hat(c0)` is its empirical (`type = 1`) `1 - alpha` quantile;

- `fw_size(c0) = mean(Mstar_c0 > z_pstar)`, at the rounded admission
  rule;

- admission is unchanged in form, `T(g) >= kappa_hat(c0)`;

- guidance on what to set to run `kappa_hat(c0)` as a `p*` screen at the
  fit's `digits`: `pstar_settable`, the smallest `p*` on the
  `10^-digits` grid whose effective threshold is at or above
  `kappa_hat(c0)`; that threshold as `pcons_eff_settable` (Pcons scale)
  and `z_eff_settable` (z scale);
  `z_gap = z_eff_settable - kappa_hat(c0)` (positive = conservative);
  and `digits_fine` / `pstar_fine` / `z_gap_fine`, the smallest `digits`
  (searched over 1 to 12) at which the gap falls below 0.01, with its
  `p*`. When no `p* <= 1` reaches `kappa_hat(c0)` at the fit's `digits`,
  `pstar_achievable` is `FALSE` and the settable columns are `NA`; the
  print method says so. This is guidance on what to set, not an
  identity: the settable screen is `kappa_hat(c0)` rounded up to what
  the grid can express, so it admits a subset of what `kappa_hat(c0)`
  admits.

If every candidate's true effect is `c0_cmp`, `T(g)` is centred at
`-delta_g`, so the null law of `max_g T(g)` is that of the shifted
maximum. At `c0 = c2`, `delta_g = 0` and every quantity is the unshifted
one.

`c0` is on the natural scale of the consistency threshold `c2` (the HR
on the survival path; the ratio for OR / RR / IRR; the difference for RD
/ MD) and is mapped exactly as `c2` is:
[`log()`](https://rdrr.io/r/base/Log.html) for ratio measures, identity
otherwise. `c0 <= c2` is required. A non-positive `c0` on a ratio path
(a log supplied by mistake) errors; on an identity path a wrong-scale
`c0` cannot be detected.

The shifted maxima are read from the capture when the fit carries them
for every requested `c0` (`declaration_c0` at fit time); otherwise they
are computed from the stored field matrix (`keep_field_matrix = TRUE`);
otherwise the call errors. It never falls back to the unshifted maximum.
`family = "reduced"` always computes from the field matrix.

## Inputs

The field is retained only on request. Fit with `mr_inference = TRUE`
and `mr_inference_args = list(keep_declaration_field = TRUE)` (add
`keep_field_matrix = TRUE` for the reduced diagnostic), or call
`fs_mr_inference(..., keep_declaration_field = TRUE)` and pass its
result. The multiplier law and `B` are those of the MR draws and are
recorded in the result.

## See also

[`fs_mr_inference()`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md),
[`fs_family_report()`](https://larry-leon.github.io/forestsearch/reference/fs_family_report.md),
[`fs_fdr_report()`](https://larry-leon.github.io/forestsearch/reference/fs_fdr_report.md).

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- forestsearch(df, ..., mr_inference = TRUE,
                    mr_inference_args = list(keep_declaration_field = TRUE))
fs_declaration_calibration(fit, alpha = 0.05)
} # }
```
