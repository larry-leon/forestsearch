# =============================================================================
# fs_declaration_calibration.R
#
# Calibrated declaration threshold (kappa-hat-alpha) and family-wise size
# diagnostic (FW-hat-alpha) for the p-star consistency screen.
# TASK_declaration_calibration_2026-09-22_v2.
#
# The closed-form consistency screen admits candidate g when
#   max(0, 2 * pnorm((beta_hat - c_cons) / sigma_D) - 1) >= p_star,
# the rate rounded to pconsistency.digits = d first, which is exactly
# T(g) = (beta_hat(g) - c_cons) / sigma_D(g) >= qnorm((1 + pcons_eff) / 2),
# pcons_eff = .fs_pcons_eff(p_star, d) (TASK_declcal_rounding_alignment_2026-09-23).
# Under multiplier resampling the standardized perturbation field is
#   Zstar[b, g] = sum_i xi[b, i] * db[g, i] / sigma_D(g),
# with one shared multiplier vector per draw.  Its family maximum Mstar[b]
# calibrates the declaration: kappa_hat(alpha) is the (1 - alpha) quantile of
# Mstar, and fw_size = mean(Mstar > qnorm((1 + pcons_eff) / 2)) is the
# family-wise size the p-star screen, as implemented, has on the realized family.
#
# Post-hoc, opt-in and reported: nothing here re-runs the search, mutates the
# fit, or feeds back into admission.
# =============================================================================


#' Assemble and column-standardize the declaration field
#'
#' The one place the standardization lives.  Callable on a supplied influence
#' matrix and multiplier matrix directly, so it can be exercised without a fit.
#'
#' @param db `N x G` influence matrix: column `g` holds candidate `g`'s
#'   treatment `dfbeta` in its members' rows and zero elsewhere.
#' @param xi `N x B` multiplier matrix, one column per draw, shared across all
#'   candidates (mean 0, variance 1).
#' @param keep_matrix Logical; return the `B x G` matrix `Zstar`.
#' @param cor_max Largest family for which the empirical correlation matrix of
#'   `Zstar` is returned.
#' @param shift `NULL` (default: no shift), a length-`G` numeric vector, or a
#'   `G x K` matrix of per-candidate shifts `delta[g, k]`.  When supplied, the
#'   shifted maxima `max_g { Zstar[b, g] - delta[g, k] }` are returned as
#'   `Mstar_shift` (`B x K`, column names those of `shift`).  They are taken
#'   from the full field before `keep_matrix` is consulted, so they are
#'   available without storing the matrix.
#' @return List with `Mstar` (length `B`), `Zstar` (or `NULL`), `sigma_D`,
#'   `column_sd`, `column_mean`, `zstar_mean`, `field_cor` (or `NULL`),
#'   `shared_multipliers`, `B`, `G`; and `Mstar_shift` only when `shift` is
#'   supplied.
#' @keywords internal
#' @noRd
.fs_decl_field <- function(db, xi, keep_matrix = TRUE, cor_max = 8L,
                           shift = NULL) {
  db <- as.matrix(db)
  xi <- as.matrix(xi)
  if (nrow(db) != nrow(xi)) {
    stop("db and xi must have the same number of rows (subjects).", call. = FALSE)
  }
  sigma_D <- sqrt(colSums(db^2))
  if (any(!is.finite(sigma_D) | sigma_D <= 0)) {
    stop("every candidate needs a positive, finite sigma_D.", call. = FALSE)
  }
  # One crossprod: draw b's single multiplier vector xi[, b] enters every
  # candidate's perturbation, so the multipliers are shared by construction.
  zs <- crossprod(xi, db)                         # B x G : D[b, g]
  zs <- zs / rep(sigma_D, each = nrow(zs))        # column g / sigma_D(g)
  n_g <- ncol(zs)
  m_star <- zs[cbind(seq_len(nrow(zs)), max.col(zs, ties.method = "first"))]
  m_shift <- if (is.null(shift)) NULL else .fs_decl_shifted_max(zs, shift)
  out <- list(
    Mstar = m_star,
    Zstar = if (isTRUE(keep_matrix)) zs else NULL,
    sigma_D = sigma_D,
    column_sd = apply(zs, 2L, stats::sd),
    column_mean = colMeans(zs),
    zstar_mean = mean(zs),
    field_cor = if (n_g <= cor_max) stats::cor(zs) else NULL,
    shared_multipliers = TRUE,
    B = nrow(zs), G = n_g
  )
  if (!is.null(m_shift)) out$Mstar_shift <- m_shift
  out
}


#' Shifted family maxima of a standardized field
#'
#' `M[b, k] = max_g { zs[b, g] - shift[g, k] }`.  With a zero shift column the
#' result is bit-identical to the unshifted maximum.
#'
#' @param zs `B x G` standardized field.
#' @param shift Length-`G` vector or `G x K` matrix.
#' @return `B x K` matrix, column names from `shift`.
#' @keywords internal
#' @noRd
.fs_decl_shifted_max <- function(zs, shift) {
  shift <- as.matrix(shift)
  if (nrow(shift) != ncol(zs) || !all(is.finite(shift))) {
    stop("shift must be finite with one row per candidate (", ncol(zs),
         ").", call. = FALSE)
  }
  n_b <- nrow(zs)
  rows <- seq_len(n_b)
  m <- vapply(seq_len(ncol(shift)), function(k) {
    zk <- zs - rep(shift[, k], each = n_b)
    zk[cbind(rows, max.col(zk, ties.method = "first"))]
  }, numeric(n_b))
  m <- matrix(m, n_b, ncol(shift))
  colnames(m) <- colnames(shift)
  m
}


#' Map a protected null level c0 to the comparison scale
#'
#' `c0` is supplied on the natural scale of the consistency threshold `c2`
#' (the HR itself on the survival path) and transformed exactly as `c2` is:
#' `log()` for ratio measures, the identity otherwise.  `c0 <= c2` is required
#' on the natural scale, i.e. `c0_cmp <= c_cons` on the comparison scale.
#'
#' @param c0 Numeric vector, natural scale.
#' @param c_cons `c2` on the comparison scale.
#' @param log_scale Logical; `TRUE` for ratio measures.
#' @return Numeric vector `c0_cmp`, named by `as.character(c0)`.
#' @keywords internal
#' @noRd
.fs_decl_c0_cmp <- function(c0, c_cons, log_scale) {
  if (!is.numeric(c0) || !length(c0) || any(!is.finite(c0))) {
    stop("c0 must be a finite numeric vector.", call. = FALSE)
  }
  if (anyDuplicated(c0)) stop("c0 values must be distinct.", call. = FALSE)
  if (is.null(c_cons)) {
    stop("c0 needs a consistency threshold c2 to shift from, and the ",
         "resolved admission set has none.", call. = FALSE)
  }
  if (is.null(log_scale)) {
    stop("the field records no log_scale, so c0 cannot be mapped to the ",
         "comparison scale.", call. = FALSE)
  }
  if (isTRUE(log_scale)) {
    if (any(c0 <= 0)) {
      stop("c0 = ", paste(c0[c0 <= 0], collapse = ", "), " is not on the ",
           "natural ratio scale: on this path c0 is a ratio (e.g. 0.75 for ",
           "HR 0.75), not its log.", call. = FALSE)
    }
    c0_cmp <- log(c0)
  } else {
    c0_cmp <- c0
  }
  bad <- c0_cmp > c_cons
  if (any(bad)) {
    c2_nat <- if (isTRUE(log_scale)) exp(c_cons) else c_cons
    stop("c0 must satisfy c0 <= c2 on the natural scale: c0 = ",
         paste(c0[bad], collapse = ", "), " exceeds c2 = ",
         format(c2_nat, digits = 6),
         ". A protected level worse than the harm criterion is not a null.",
         call. = FALSE)
  }
  stats::setNames(c0_cmp, as.character(c0))
}


#' Per-candidate shift matrix for protected levels c0
#' @param c0_cmp Named comparison-scale levels (from `.fs_decl_c0_cmp()`).
#' @param c_cons `c2` on the comparison scale.
#' @param sigma_D Length-`G` robust scales.
#' @return `G x K` matrix `(c_cons - c0_cmp[k]) / sigma_D[g]`.
#' @keywords internal
#' @noRd
.fs_decl_c0_shift <- function(c0_cmp, c_cons, sigma_D) {
  n_g <- length(sigma_D)
  d <- matrix((c_cons - rep(unname(c0_cmp), each = n_g)) / sigma_D,
              n_g, length(c0_cmp))
  colnames(d) <- names(c0_cmp)
  d
}


#' Effective Pcons threshold of the rounded admission rule
#'
#' The screen admits on `round(Pcons, digits) >= p_star`.  With `g` the
#' smallest point of the `10^-digits` grid at or above `p_star`, that is
#' `Pcons >= g - 0.5 * 10^-digits`, which this returns.  The one place the
#' expression lives; every site that needs the screen's threshold calls it.
#' The scaled `p_star` is snapped to 6 decimals before `ceiling()` so that a
#' grid value carrying representation error (`0.07 * 100` is
#' `7.000000000000001`) is not pushed up a grid step.
#'
#' @param p_star Numeric, the consistency threshold `p*`.
#' @param digits Integer, `pconsistency.digits`.
#' @return Numeric, the threshold on the `Pcons` scale.
#' @keywords internal
#' @noRd
.fs_pcons_eff <- function(p_star, digits) {
  g <- ceiling(round(p_star * 10^digits, 6)) / 10^digits
  g - 0.5 * 10^(-digits)
}


#' Settable p-star for a calibrated cutoff
#'
#' The smallest `p*` on the `10^-digits` grid whose effective threshold
#' (`.fs_pcons_eff()`, mapped to z by `qnorm((1 + .) / 2)`) is at or above
#' `kap`.  `NA` when no `p* <= 1` reaches `kap` at these `digits`.
#' @return List with `p_star`, `pcons_eff`, `z_eff`, `z_gap` (`z_eff - kap`).
#' @keywords internal
#' @noRd
.fs_decl_settable <- function(kap, digits) {
  none <- list(p_star = NA_real_, pcons_eff = NA_real_, z_eff = NA_real_,
               z_gap = NA_real_)
  if (!is.finite(kap)) return(none)
  step <- 10^(-digits)
  p_k <- 2 * stats::pnorm(kap) - 1
  # start at p_k on the grid; the effective threshold sits within a step, so
  # the searches below move at most a step or two, each via .fs_pcons_eff()
  g <- max(step, ceiling(round(p_k / step, 6)) * step)
  z_at <- function(p) stats::qnorm((1 + .fs_pcons_eff(p, digits)) / 2)
  while (g - step >= step && z_at(g - step) >= kap) g <- g - step
  while (g <= 1 + step / 2 && z_at(g) < kap) g <- g + step
  if (g > 1 + step / 2) return(none)
  g <- round(g, digits)
  z <- z_at(g)
  list(p_star = g, pcons_eff = .fs_pcons_eff(g, digits), z_eff = z,
       z_gap = z - kap)
}


#' Settable-pair columns for a vector of calibrated cutoffs
#'
#' At the fit's `digits`: the smallest settable `p*` reaching each `kap`, its
#' effective threshold on the `Pcons` and z scales, and the gap to `kap` in z
#' units (positive = conservative).  Then the smallest `digits` in
#' `1:digits_max` at which that gap falls below `gap_tol`, with its `p*`.
#' @keywords internal
#' @noRd
.fs_decl_settable_table <- function(kap, digits, gap_tol = 0.01,
                                    digits_max = 12L) {
  rows <- lapply(kap, function(k) {
    s <- .fs_decl_settable(k, digits)
    d_fine <- NA_integer_; p_fine <- NA_real_; gap_fine <- NA_real_
    for (d in seq_len(digits_max)) {
      sd <- .fs_decl_settable(k, d)
      if (is.finite(sd$z_gap) && sd$z_gap < gap_tol) {
        d_fine <- d; p_fine <- sd$p_star; gap_fine <- sd$z_gap
        break
      }
    }
    data.frame(pstar_settable = s$p_star,
               pstar_achievable = is.finite(s$p_star),
               pcons_eff_settable = s$pcons_eff,
               z_eff_settable = s$z_eff,
               z_gap = s$z_gap,
               digits_fine = d_fine, pstar_fine = p_fine,
               z_gap_fine = gap_fine)
  })
  do.call(rbind, rows)
}


#' Normalize a candidate label to an order-free key
#' @keywords internal
#' @noRd
.fs_decl_key <- function(labels) {
  vapply(strsplit(labels, " & ", fixed = TRUE),
         function(p) paste(sort(p), collapse = " & "), character(1))
}


#' Replay the near-duplicate reduction on a forestsearch fit
#'
#' Reads the fit's retained candidate table and applies the same filter and the
#' same `remove_near_duplicate_subgroups()` call the consistency stage applied,
#' to identify which members of the pre-reduction family the reduction removed.
#' Deterministic; nothing is refit.
#'
#' @param fs A `forestsearch` object.
#' @param fam_pre Character vector, the pre-reduction family labels.
#' @return `NULL` when no candidate table is available, otherwise a list with
#'   `applicable`, `removed`, `screened`, `admitted_current`, `n_unmatched`,
#'   `replay_check` and `note`.
#' @keywords internal
#' @noRd
.fs_decl_reduction <- function(fs, fam_pre) {
  hs <- fs$find.grps$out.found$hr.subgroups
  if (is.null(hs) || !nrow(hs)) return(NULL)
  hs <- as.data.frame(hs)
  if (identical(fs$sg_focus, "maxeff")) {
    return(list(applicable = FALSE, removed = character(0), screened = NULL,
                admitted_current = NULL, n_unmatched = NA_integer_,
                replay_check = NA,
                note = paste("sg_focus = \"maxeff\" uses exact membership",
                             "de-duplication, not the near-duplicate",
                             "reduction; the reduced-family diagnostic does",
                             "not apply.")))
  }
  names_z <- setdiff(names(hs), c("K", "n", "E", "d0", "d1", "m0", "m1", "HR",
                                  "L(HR)", "U(HR)", "grp"))
  row_key <- function(d) {
    ind <- as.matrix(d[, names_z, drop = FALSE]) == 1
    .fs_decl_key(apply(ind, 1L, function(r) paste(names_z[r], collapse = " & ")))
  }

  # The consistency stage's filter, on the comparison scale of the admission
  # floor: the survival HR column is natural-scale, the GLM column is already
  # on the comparison scale.
  floor_cmp <- fs$admission$effect_floor
  found <- hs
  m1_thr <- fs$args_call_all$m1.threshold
  if (!is.null(m1_thr) && is.finite(m1_thr)) {
    found <- found[!is.na(found$m1), , drop = FALSE]
    found <- found[found$m1 <= m1_thr, , drop = FALSE]
  }
  if (!is.null(floor_cmp)) {
    eff <- if (identical(fs$outcome_type, "survival")) log(found$HR) else found$HR
    found <- found[eff >= floor_cmp, , drop = FALSE]
  }
  kept <- if (nrow(found) > 1L) {
    as.data.frame(remove_near_duplicate_subgroups(found))
  } else {
    found
  }
  removed_rows <- found[!found$grp %in% kept$grp, , drop = FALSE]

  fam_key <- .fs_decl_key(fam_pre)
  to_family <- function(keys) fam_pre[fam_key %in% keys]
  kept_key <- row_key(kept)
  removed <- to_family(row_key(removed_rows))
  screened <- to_family(kept_key)

  res <- fs$grp.consistency$out_sg$result
  admitted_current <- if (is.null(res) || !nrow(res)) character(0) else {
    to_family(row_key(hs[match(as.character(res$g), as.character(hs$grp)), ,
                         drop = FALSE]))
  }
  n_total <- fs$grp.consistency$n_candidates_total
  list(applicable = TRUE, removed = removed, screened = screened,
       admitted_current = admitted_current,
       n_unmatched = sum(!kept_key %in% fam_key),
       replay_check = if (is.null(n_total)) NA else nrow(kept) == n_total,
       note = NULL)
}


#' Calibrated declaration threshold and family-wise size of the p-star screen
#'
#' Post-hoc, opt-in, reported diagnostics of the consistency screen, computed
#' from multiplier draws the package already produces.  The screen admits a
#' candidate `g` when its consistency rate, rounded to `pconsistency.digits`,
#' is at least `p_star`; on the resample path that is exactly the standardized
#' statistic `T(g) = (beta_hat(g) - c_cons) / sigma_D(g)` clearing the
#' effective cutoff `z_pstar = qnorm((1 + pcons_eff) / 2)` (see the Rounded
#' admission rule section).  This function reports how often the maximum of
#' that statistic over the candidate family would clear the screen's cutoff
#' under the null perturbation law (`fw_size`), and the cutoff that would hold
#' the family-wise declaration rate at `alpha` (`kappa_hat`).
#'
#' It reads a fitted object and returns a new one.  It does not modify its
#' input, does not re-run the search, does not call the consistency engine,
#' and does not change what the search admitted: `admitted_calibrated` is the
#' set the calibrated rule *would* admit, reported beside `admitted_current`.
#'
#' @section Definitions:
#' With shared multipliers `xi[b, i]` (one vector per draw `b`, reused across
#' every candidate) and the dfbeta influence `db[g, i]`:
#'
#' * `Zstar[b, g] = sum_i xi[b, i] * db[g, i] / sigma_D(g)`, with
#'   `sigma_D(g)^2 = sum_i db[g, i]^2` -- the robust scale the screen itself
#'   uses, never a model-based standard error.
#' * `Mstar[b] = max_g Zstar[b, g]` over the family: one-sided, in the harm
#'   direction.
#' * `kappa_hat(alpha)` is the empirical (`type = 1`) `1 - alpha` quantile of
#'   `Mstar`.
#' * `fw_size = mean(Mstar > z_pstar)`, `z_pstar = qnorm((1 + pcons_eff) / 2)`.
#'   Its threshold is the fit's own `p_star` at the fit's own
#'   `pconsistency.digits`, not `alpha`: it is the family-wise size of the
#'   screen as implemented, not of the calibrated rule.
#' * The calibrated admission rule is
#'   `beta_hat(g) >= max(c_screen, c_cons + kappa_hat * sigma_D(g))`, the
#'   current rule with `z_pstar` replaced by `kappa_hat`.
#'
#' @section Rounded admission rule:
#' The screen admits on `round(Pcons, digits) >= p_star`, with `digits` the
#' fit's `pconsistency.digits`.  With `p_star` rounded up to the
#' `10^-digits` grid (`g`), that is `Pcons >= pcons_eff = g - 0.5 * 10^-digits`,
#' a bar below `p_star` itself (0.895 for `p_star = 0.90`, `digits = 2`).  On
#' the resample path `Pcons = 2 * pnorm(T) - 1`, so the screen is
#' `T >= z_pstar = qnorm((1 + pcons_eff) / 2)`.  `fw_size`, `admitted_pstar`
#' and the settable `p*` are computed at that effective threshold; `fw_size`
#' is therefore the family-wise size of the screen as implemented and depends
#' on `pconsistency.digits`.  `kappa_hat` does not: it is a quantile of the
#' maximum statistic and is compared with `T` directly, with no rounding.
#'
#' `digits` is read from the fit's `args_call_all$pconsistency.digits`; when
#' absent (a bare `fs_mr_inference()` result, or an older fit) the
#' `subgroup.consistency()` default of 2 is used.  `digits_source` records
#' which.
#'
#' The correspondence with `T` holds only under
#' `consistency_method = "resample"`.  Under `"split"`, `Pcons` is a split
#' proportion `k / n_valid` and admission is not a threshold on `T`, so
#' `fw_size`, `z_pstar`, `admitted_pstar` and the settable-`p*` columns are
#' `NA` / `NULL` rather than computed; `kappa_hat` and the calibrated rule are
#' still returned.  The method is read from `args_call_all$consistency_method`;
#' a bare `fs_mr_inference()` result is treated as resample (the field is the
#' closed-form statistic), recorded in `consistency_method_source`.
#'
#' @section Family:
#' `family = "prereduction"` (the default, and the only value meant for the
#' calibrated rule) takes the maximum over the multiplier-resampling family:
#' every enumerated candidate of at most `maxk` factors meeting the size
#' minimum, before any reduction keyed on fitted quantities.  The
#' near-duplicate reduction keys on sample-fitted summaries, so a family
#' reduced by it is outcome-dependent; the pre-reduction family is
#' covariate-measurable, and is also the conservative choice (a larger family
#' raises the maximum).  The per-arm event minima are not replayed in that
#' family, which makes it a superset of the one the search chose among.
#'
#' `family = "reduced"` removes the candidates the near-duplicate reduction
#' removed and is a diagnostic only, so the gap is measurable.  Its result is
#' conditional on the realized family and is labelled so; it is never meant
#' for admission.  It needs a `forestsearch` fit (the reduction lives in the
#' identifier) and the full field matrix (`keep_field_matrix = TRUE`).
#'
#' @section Protected null level:
#' `c1` and `c2` state the claim; `c0` states what the claim is protected
#' against: a clinically specified, pre-specified benefit level (for example
#' HR 0.75).  The family-wise declaration rate is then controlled at `alpha`
#' whenever every candidate's true effect is at least as good as `c0`, rather
#' than only when every candidate sits at `c2`.  With
#' `delta_g = (c_cons - c0_cmp) / sigma_D(g)`:
#'
#' * `Mstar_c0[b] = max_g { Zstar[b, g] - delta_g }`;
#' * `kappa_hat(c0)` is its empirical (`type = 1`) `1 - alpha` quantile;
#' * `fw_size(c0) = mean(Mstar_c0 > z_pstar)`, at the rounded admission rule;
#' * admission is unchanged in form, `T(g) >= kappa_hat(c0)`;
#' * guidance on what to set to run `kappa_hat(c0)` as a `p*` screen at the
#'   fit's `digits`: `pstar_settable`, the smallest `p*` on the `10^-digits`
#'   grid whose effective threshold is at or above `kappa_hat(c0)`; that
#'   threshold as `pcons_eff_settable` (Pcons scale) and `z_eff_settable`
#'   (z scale); `z_gap = z_eff_settable - kappa_hat(c0)` (positive =
#'   conservative); and `digits_fine` / `pstar_fine` / `z_gap_fine`, the
#'   smallest `digits` (searched over 1 to 12) at which the gap falls below
#'   0.01, with its `p*`.  When no `p* <= 1` reaches `kappa_hat(c0)` at the
#'   fit's `digits`, `pstar_achievable` is `FALSE` and the settable columns
#'   are `NA`; the print method says so.  This is guidance on what to set, not
#'   an identity: the settable screen is `kappa_hat(c0)` rounded up to what the
#'   grid can express, so it admits a subset of what `kappa_hat(c0)` admits.
#'
#' If every candidate's true effect is `c0_cmp`, `T(g)` is centred at
#' `-delta_g`, so the null law of `max_g T(g)` is that of the shifted maximum.
#' At `c0 = c2`, `delta_g = 0` and every quantity is the unshifted one.
#'
#' `c0` is on the natural scale of the consistency threshold `c2` (the HR on
#' the survival path; the ratio for OR / RR / IRR; the difference for RD /
#' MD) and is mapped exactly as `c2` is: `log()` for ratio measures, identity
#' otherwise.  `c0 <= c2` is required.  A non-positive `c0` on a ratio path
#' (a log supplied by mistake) errors; on an identity path a wrong-scale
#' `c0` cannot be detected.
#'
#' The shifted maxima are read from the capture when the fit carries them
#' for every requested `c0` (`declaration_c0` at fit time); otherwise they are
#' computed from the stored field matrix (`keep_field_matrix = TRUE`);
#' otherwise the call errors.  It never falls back to the unshifted maximum.
#' `family = "reduced"` always computes from the field matrix.
#'
#' @section Inputs:
#' The field is retained only on request.  Fit with
#' `mr_inference = TRUE` and
#' `mr_inference_args = list(keep_declaration_field = TRUE)` (add
#' `keep_field_matrix = TRUE` for the reduced diagnostic), or call
#' `fs_mr_inference(..., keep_declaration_field = TRUE)` and pass its result.
#' The multiplier law and `B` are those of the MR draws and are recorded in
#' the result.
#'
#' @param fit A `forestsearch` object fitted with the declaration field
#'   retained, or the list returned by `fs_mr_inference()` with
#'   `keep_declaration_field = TRUE`.
#' @param alpha Target family-wise declaration rate, in (0, 1). Default `0.05`.
#' @param family `"prereduction"` (default) or `"reduced"`; see the Family
#'   section.
#' @param ... Only `p_star` and `c_cons` are read, and only when the fit's
#'   resolved admission set carries no consistency floor (for example
#'   `sg_focus = "maxeff"`); they then name the hypothetical screen whose size
#'   is evaluated.  `c_cons` is on the comparison scale (log for ratio
#'   measures).
#' @param c0 `NULL` (default) or a numeric vector of clinically specified,
#'   pre-specified protected null levels; see the Protected null level
#'   section.  Named only (it follows `...`).
#' @return An object of class `fs_declaration_calibration`: a list with
#'   `kappa_hat`, `fw_size`, `alpha`, `p_star`, `digits`, `digits_source`,
#'   `pcons_eff`, `z_pstar` (the effective z cutoff of the rounded screen),
#'   `consistency_method`, `consistency_method_source`, `c_cons`,
#'   `c_screen`, `B`, `multiplier_law`, `quantile_type`, `family_source`,
#'   `family_label`, `n_family_prereduction`, `n_family_reduced`,
#'   `admitted_current` (the candidates the executed screen admitted; `NULL`
#'   when the fit is a bare `fs_mr_inference()` result), `admitted_calibrated`,
#'   `admitted_pstar` (the relabelled current rule over the same family),
#'   `screened` (the candidates the consistency screen evaluated, when known),
#'   `Mstar`, `beta_hat`, `sigma_D`, `T_hat`, `column_sd`, `zstar_mean`,
#'   `field_cor` (families of at most 8 candidates, else `NULL`), and a
#'   `reduction` list recording the replay of the near-duplicate reduction.
#'   When `c0` is given it also carries `c0`, a list with `table` (one row per
#'   `c0`: `c0`, `c0_cmp`, `kappa_hat`, `fw_size`, `pstar_settable`,
#'   `pstar_achievable`, `pcons_eff_settable`, `z_eff_settable`, `z_gap`,
#'   `digits_fine`, `pstar_fine`, `z_gap_fine`, `n_admitted_calibrated`, the 0.90 / 0.95 / 0.99 quantiles of `Mstar_c0`,
#'   and `is_c2`), `admitted_calibrated` (a list indexed by `c0`), `Mstar_c0`
#'   and `source` (`"capture"` or `"field_matrix"`).  Every other element is
#'   unchanged by `c0`.
#' @seealso [fs_mr_inference()], [fs_family_report()], [fs_fdr_report()].
#' @examples
#' \dontrun{
#' fit <- forestsearch(df, ..., mr_inference = TRUE,
#'                     mr_inference_args = list(keep_declaration_field = TRUE))
#' fs_declaration_calibration(fit, alpha = 0.05)
#' }
#' @export
fs_declaration_calibration <- function(fit,
                                       alpha = 0.05,
                                       family = c("prereduction", "reduced"),
                                       ...,
                                       c0 = NULL) {
  family <- match.arg(family)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number in (0, 1).", call. = FALSE)
  }
  dots <- list(...)

  is_fs <- inherits(fit, "forestsearch")
  mr <- if (is_fs) fit$mr_inference else fit
  if (is.null(mr) || !is.list(mr)) {
    stop("fit carries no multiplier-resampling result. Fit with ",
         "mr_inference = TRUE and ",
         "mr_inference_args = list(keep_declaration_field = TRUE).",
         call. = FALSE)
  }
  fld <- mr$declaration_field
  if (is.null(fld)) {
    stop("fit has no declaration field: it is retained only when ",
         "keep_declaration_field = TRUE is set (forestsearch(): ",
         "mr_inference_args = list(keep_declaration_field = TRUE); ",
         "fs_mr_inference(): keep_declaration_field = TRUE). ",
         "Nothing is computed from partial inputs.", call. = FALSE)
  }

  meta <- fld$meta
  p_star <- meta$p_star
  c_cons <- meta$c_cons
  if (is.null(p_star) || is.null(c_cons)) {
    if (is.null(dots$p_star) || is.null(dots$c_cons)) {
      stop("the fit's resolved admission set has no consistency floor ",
           "(e.g. sg_focus = \"maxeff\"), so there is no p-star screen to ",
           "size. Supply p_star = and c_cons = through ... to evaluate a ",
           "hypothetical screen.", call. = FALSE)
    }
    p_star <- dots$p_star
    c_cons <- dots$c_cons
  }
  c_screen <- meta$c_screen

  # The screen admits on round(Pcons, digits) >= p_star; digits is the fit's
  # own pconsistency.digits, else the subgroup.consistency() default.
  aca <- if (is_fs) fit$args_call_all else NULL
  digits <- aca$pconsistency.digits
  digits_source <- "fit$args_call_all$pconsistency.digits"
  if (is.null(digits)) {
    digits <- eval(formals(subgroup.consistency)$pconsistency.digits)
    digits_source <- "subgroup.consistency() default (not recorded on the fit)"
  }
  digits <- as.integer(digits)
  cons_method <- aca$consistency_method
  cons_method_source <- "fit$args_call_all$consistency_method"
  if (is.null(cons_method)) {
    cons_method <- "resample"
    cons_method_source <- paste("assumed: the field is the closed-form",
                                "(resample) statistic; not recorded on the fit")
  }
  # Pcons = 2 * pnorm(T) - 1 holds on the resample path only; under "split"
  # Pcons is k / n_valid and the screen is not a threshold on T.
  on_T <- identical(cons_method, "resample")
  pcons_eff <- .fs_pcons_eff(p_star, digits)
  z_pstar <- if (on_T) stats::qnorm((1 + pcons_eff) / 2) else NA_real_

  fam_pre <- fld$family_id
  red <- if (is_fs) .fs_decl_reduction(fit, fam_pre) else NULL
  red_ok <- !is.null(red) && isTRUE(red$applicable)
  n_pre <- length(fam_pre)
  n_red <- if (red_ok) n_pre - length(red$removed) else NA_integer_

  if (identical(family, "prereduction")) {
    cols <- seq_len(n_pre)
    m_star <- fld$Mstar
    label <- "pre-reduction family (covariate-measurable)"
  } else {
    if (!red_ok) {
      stop("family = \"reduced\" needs a forestsearch fit whose ",
           "near-duplicate reduction can be replayed",
           if (!is.null(red$note)) paste0(": ", red$note) else ".",
           call. = FALSE)
    }
    if (is.null(fld$Zstar)) {
      stop("family = \"reduced\" needs the full field matrix: refit with ",
           "keep_field_matrix = TRUE alongside keep_declaration_field = TRUE.",
           call. = FALSE)
    }
    cols <- which(!fam_pre %in% red$removed)
    zs <- fld$Zstar[, cols, drop = FALSE]
    m_star <- zs[cbind(seq_len(nrow(zs)), max.col(zs, ties.method = "first"))]
    label <- paste("reduced family -- CONDITIONAL ON THE REALIZED FAMILY;",
                   "diagnostic only, not for admission")
  }

  kappa_hat <- stats::quantile(m_star, 1 - alpha, type = 1, names = FALSE)
  fw_size <- if (on_T) mean(m_star > z_pstar) else NA_real_

  bh <- fld$beta_hat[cols]
  sdv <- fld$sigma_D[cols]
  floor_at <- function(k) {
    fl <- c_cons + k * sdv
    if (!is.null(c_screen)) fl <- pmax(c_screen, fl)
    fl
  }
  admitted_calibrated <- names(bh)[bh >= floor_at(kappa_hat)]
  admitted_pstar <- if (on_T) names(bh)[bh >= floor_at(z_pstar)] else NULL

  c0_out <- if (is.null(c0)) NULL else
    .fs_decl_c0_block(fld, c0, c_cons, cols, family, alpha, z_pstar, bh, sdv,
                      floor_at, digits, on_T)

  field_cor <- meta$field_cor
  if (!is.null(field_cor) && length(cols) != n_pre) {
    field_cor <- field_cor[cols, cols, drop = FALSE]
  }

  out <- list(
    kappa_hat = kappa_hat,
    fw_size = fw_size,
    alpha = alpha,
    p_star = p_star,
    digits = digits,
    digits_source = digits_source,
    pcons_eff = pcons_eff,
    z_pstar = z_pstar,
    consistency_method = cons_method,
    consistency_method_source = cons_method_source,
    c_cons = c_cons,
    c_screen = c_screen,
    B = length(m_star),
    multiplier_law = meta$multiplier,
    quantile_type = 1L,
    family_source = family,
    family_label = label,
    n_family_prereduction = n_pre,
    n_family_reduced = n_red,
    admitted_current = if (red_ok) red$admitted_current else NULL,
    admitted_calibrated = admitted_calibrated,
    admitted_pstar = admitted_pstar,
    screened = if (red_ok) red$screened else NULL,
    Mstar = m_star,
    beta_hat = bh,
    sigma_D = sdv,
    T_hat = (bh - c_cons) / sdv,
    column_sd = meta$column_sd[cols],
    zstar_mean = if (length(cols) == n_pre) meta$zstar_mean
                 else mean(fld$Zstar[, cols]),
    field_cor = field_cor,
    shared_multipliers = meta$shared_multipliers,
    reduction = red
  )
  if (!is.null(c0_out)) out$c0 <- c0_out
  class(out) <- c("fs_declaration_calibration", "list")
  out
}


#' Per-c0 calibration block of fs_declaration_calibration()
#'
#' Reads the captured shifted maxima when they cover every requested `c0`
#' (pre-reduction family only), otherwise computes them from the stored field
#' matrix, otherwise errors.  Never substitutes the unshifted maximum.
#' @keywords internal
#' @noRd
.fs_decl_c0_block <- function(fld, c0, c_cons, cols, family, alpha, z_pstar,
                              bh, sdv, floor_at, digits, on_T = TRUE) {
  meta <- fld$meta
  c0_cmp <- .fs_decl_c0_cmp(c0, c_cons, meta$log_scale)
  keys <- names(c0_cmp)
  cap <- fld$Mstar_c0
  use_cap <- identical(family, "prereduction") && !is.null(cap) &&
    all(keys %in% colnames(cap)) && identical(meta$c_cons, c_cons)
  if (use_cap) {
    m_c0 <- cap[, keys, drop = FALSE]
    src <- "capture"
  } else if (!is.null(fld$Zstar)) {
    m_c0 <- .fs_decl_shifted_max(fld$Zstar[, cols, drop = FALSE],
                                 .fs_decl_c0_shift(c0_cmp, c_cons, sdv))
    src <- "field_matrix"
  } else {
    stop("c0 = ", paste(c0, collapse = ", "), " needs the shifted maxima: ",
         "refit with declaration_c0 = c(", paste(c0, collapse = ", "),
         ") (forestsearch(): mr_inference_args = list(keep_declaration_field ",
         "= TRUE, declaration_c0 = ...)), or with keep_field_matrix = TRUE",
         if (identical(family, "reduced")) " (required for family = \"reduced\")",
         ". The unshifted maximum is never substituted.", call. = FALSE)
  }
  kap <- apply(m_c0, 2L, stats::quantile, probs = 1 - alpha, type = 1,
               names = FALSE)
  fw <- if (on_T) apply(m_c0 > z_pstar, 2L, mean)   # as mean(Mstar > z)
        else rep(NA_real_, length(keys))
  adm <- lapply(stats::setNames(kap, keys),
                function(k) names(bh)[bh >= floor_at(k)])
  qs <- apply(m_c0, 2L, stats::quantile, probs = c(0.90, 0.95, 0.99),
              type = 1, names = FALSE)
  qs <- matrix(qs, nrow = 3L)
  tab <- data.frame(
    c0 = as.numeric(c0), c0_cmp = unname(c0_cmp),
    kappa_hat = unname(kap), fw_size = unname(fw),
    n_admitted_calibrated = lengths(adm, use.names = FALSE),
    Mstar_c0_q90 = qs[1L, ], Mstar_c0_q95 = qs[2L, ], Mstar_c0_q99 = qs[3L, ],
    is_c2 = unname(c0_cmp) == c_cons,
    row.names = NULL
  )
  # What to set to run kappa_hat as a p-star screen: defined only where the
  # screen is a threshold on T (resample path).
  st <- .fs_decl_settable_table(unname(kap), digits)
  if (!on_T) st[] <- lapply(st, function(v) rep(v[NA_integer_], length(v)))
  tab <- cbind(tab[, 1:4], st, tab[, -(1:4)])
  list(table = tab, admitted_calibrated = adm, Mstar_c0 = m_c0, source = src)
}


#' @rdname fs_declaration_calibration
#' @param x An `fs_declaration_calibration` object.
#' @export
print.fs_declaration_calibration <- function(x, ...) {
  f3 <- function(v) if (is.null(v)) "none" else formatC(v, digits = 4, format = "fg")
  f3v <- function(v) formatC(v, digits = 4, format = "fg")
  cat("Declaration calibration (post-hoc, reported; admission unchanged)\n")
  cat("  family source:       ", x$family_source, "\n")
  cat("  family:              ", x$family_label, "\n")
  if (identical(x$family_source, "reduced")) {
    cat("  CAVEAT: conditional on the realized (reduced) family; the reduction\n",
        "         keys on sample-fitted summaries. Diagnostic only.\n", sep = "")
  }
  cat("  n family (pre / red):", x$n_family_prereduction, "/",
      if (is.na(x$n_family_reduced)) "NA" else x$n_family_reduced, "\n")
  cat("  multiplier law, B:   ", x$multiplier_law, ",", x$B, "\n")
  cat("  alpha:               ", f3(x$alpha), "\n")
  cat("  kappa_hat(alpha):    ", f3(x$kappa_hat), "\n")
  cat("  p_star, digits:      ", f3(x$p_star), ",", x$digits,
      paste0("(", x$digits_source, ")"), "\n")
  cat("  screen threshold:     Pcons >=", f3(x$pcons_eff), "; T >=",
      f3(x$z_pstar), "(rounded admission rule)\n")
  cat("  fw_size (as screened):", f3(x$fw_size), "\n")
  if (!identical(x$consistency_method, "resample")) {
    cat("  NOTE: consistency_method = \"", x$consistency_method, "\": Pcons is ",
        "a split proportion, not 2 * pnorm(T) - 1, so fw_size and the settable ",
        "p-star are not defined.\n", sep = "")
  }
  cat("  c_cons, c_screen:    ", f3(x$c_cons), ",", f3(x$c_screen), "\n")
  cat("  admitted (current):  ",
      if (is.null(x$admitted_current)) "not recorded"
      else length(x$admitted_current), "\n")
  cat("  admitted (p-star rule, this family):   ", length(x$admitted_pstar), "\n")
  cat("  admitted (calibrated rule, this family):", length(x$admitted_calibrated), "\n")
  if (!is.null(x$c0)) {
    tb <- x$c0$table
    cat("  Protected null level c0 (source: ", x$c0$source, "):\n", sep = "")
    lab <- ifelse(tb$is_c2, paste0(f3v(tb$c0), " (= c2, unshifted)"),
                  f3v(tb$c0))
    shown <- data.frame(c0 = lab, kappa_hat = f3v(tb$kappa_hat),
                        fw_size = f3v(tb$fw_size),
                        pstar_set = ifelse(tb$pstar_achievable,
                                           f3v(tb$pstar_settable), "none"),
                        z_eff = f3v(tb$z_eff_settable),
                        z_gap = f3v(tb$z_gap),
                        digits_fine = tb$digits_fine,
                        pstar_fine = f3v(tb$pstar_fine),
                        n_admitted = tb$n_admitted_calibrated,
                        q95 = f3v(tb$Mstar_c0_q95))
    print(shown, row.names = FALSE, right = TRUE)
    cat("  pstar_set: smallest p* settable at digits =", x$digits,
        "whose rounded screen is at least kappa_hat (z_eff; z_gap = z_eff -",
        "kappa_hat,\n  positive = conservative); digits_fine / pstar_fine:",
        "smallest digits with z_gap < 0.01.\n")
    if (any(!tb$pstar_achievable & !is.na(tb$kappa_hat)) &&
        identical(x$consistency_method, "resample")) {
      cat("  No p* <= 1 at digits =", x$digits, "reaches kappa_hat for c0 =",
          paste(f3v(tb$c0[!tb$pstar_achievable]), collapse = ", "),
          "; use more digits (digits_fine).\n")
    }
  }
  invisible(x)
}
