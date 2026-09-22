# =============================================================================
# fs_declaration_calibration.R
#
# Calibrated declaration threshold (kappa-hat-alpha) and family-wise size
# diagnostic (FW-hat-alpha) for the p-star consistency screen.
# TASK_declaration_calibration_2026-09-22_v2.
#
# The closed-form consistency screen admits candidate g when
#   max(0, 2 * pnorm((beta_hat - c_cons) / sigma_D) - 1) >= p_star,
# which is exactly T(g) = (beta_hat(g) - c_cons) / sigma_D(g) >= z_{(1+p_star)/2}.
# Under multiplier resampling the standardized perturbation field is
#   Zstar[b, g] = sum_i xi[b, i] * db[g, i] / sigma_D(g),
# with one shared multiplier vector per draw.  Its family maximum Mstar[b]
# calibrates the declaration: kappa_hat(alpha) is the (1 - alpha) quantile of
# Mstar, and fw_size = mean(Mstar > z_{(1+p_star)/2}) is the family-wise size
# the conventional p-star screen actually has on the realized family.
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
#' @return List with `Mstar` (length `B`), `Zstar` (or `NULL`), `sigma_D`,
#'   `column_sd`, `column_mean`, `zstar_mean`, `field_cor` (or `NULL`),
#'   `shared_multipliers`, `B`, `G`.
#' @keywords internal
#' @noRd
.fs_decl_field <- function(db, xi, keep_matrix = TRUE, cor_max = 8L) {
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
  list(
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
#' candidate `g` when its standardized statistic
#' `T(g) = (beta_hat(g) - c_cons) / sigma_D(g)` clears `qnorm((1 + p_star) / 2)`;
#' that admission is exactly the closed-form consistency rule, relabelled.
#' This function reports how often the maximum of that statistic over the
#' candidate family would clear the conventional cutoff under the null
#' perturbation law (`fw_size`), and the cutoff that would hold the
#' family-wise declaration rate at `alpha` (`kappa_hat`).
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
#' * `fw_size = mean(Mstar > qnorm((1 + p_star) / 2))`.  Its threshold is the
#'   fit's own `p_star`, not `alpha`: it is the family-wise size of the screen
#'   as run, not of the calibrated rule.
#' * The calibrated admission rule is
#'   `beta_hat(g) >= max(c_screen, c_cons + kappa_hat * sigma_D(g))`, the
#'   current rule with `qnorm((1 + p_star) / 2)` replaced by `kappa_hat`.
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
#' @return An object of class `fs_declaration_calibration`: a list with
#'   `kappa_hat`, `fw_size`, `alpha`, `p_star`, `z_pstar`, `c_cons`,
#'   `c_screen`, `B`, `multiplier_law`, `quantile_type`, `family_source`,
#'   `family_label`, `n_family_prereduction`, `n_family_reduced`,
#'   `admitted_current` (the candidates the executed screen admitted; `NULL`
#'   when the fit is a bare `fs_mr_inference()` result), `admitted_calibrated`,
#'   `admitted_pstar` (the relabelled current rule over the same family),
#'   `screened` (the candidates the consistency screen evaluated, when known),
#'   `Mstar`, `beta_hat`, `sigma_D`, `T_hat`, `column_sd`, `zstar_mean`,
#'   `field_cor` (families of at most 8 candidates, else `NULL`), and a
#'   `reduction` list recording the replay of the near-duplicate reduction.
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
                                       ...) {
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
  z_pstar <- stats::qnorm((1 + p_star) / 2)

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
  fw_size <- mean(m_star > z_pstar)

  bh <- fld$beta_hat[cols]
  sdv <- fld$sigma_D[cols]
  floor_at <- function(k) {
    fl <- c_cons + k * sdv
    if (!is.null(c_screen)) fl <- pmax(c_screen, fl)
    fl
  }
  admitted_calibrated <- names(bh)[bh >= floor_at(kappa_hat)]
  admitted_pstar <- names(bh)[bh >= floor_at(z_pstar)]

  field_cor <- meta$field_cor
  if (!is.null(field_cor) && length(cols) != n_pre) {
    field_cor <- field_cor[cols, cols, drop = FALSE]
  }

  out <- list(
    kappa_hat = kappa_hat,
    fw_size = fw_size,
    alpha = alpha,
    p_star = p_star,
    z_pstar = z_pstar,
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
  class(out) <- c("fs_declaration_calibration", "list")
  out
}


#' @rdname fs_declaration_calibration
#' @param x An `fs_declaration_calibration` object.
#' @export
print.fs_declaration_calibration <- function(x, ...) {
  f3 <- function(v) if (is.null(v)) "none" else formatC(v, digits = 4, format = "fg")
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
  cat("  p_star, z_pstar:     ", f3(x$p_star), ",", f3(x$z_pstar), "\n")
  cat("  fw_size (at p_star): ", f3(x$fw_size), "\n")
  cat("  c_cons, c_screen:    ", f3(x$c_cons), ",", f3(x$c_screen), "\n")
  cat("  admitted (current):  ",
      if (is.null(x$admitted_current)) "not recorded"
      else length(x$admitted_current), "\n")
  cat("  admitted (p-star rule, this family):   ", length(x$admitted_pstar), "\n")
  cat("  admitted (calibrated rule, this family):", length(x$admitted_calibrated), "\n")
  invisible(x)
}
