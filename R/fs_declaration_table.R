#' Declaration at a protected level, as a table
#'
#' Reads the two [fs_declaration_calibration()] results a report needs from one fit (at
#' the protected levels `c0`, and at the claim threshold) with the fit's declaration
#' field, into one tidy table. Nothing is redrawn: the calibration reads the multiplier
#' draws the fit already made.
#'
#' `fw_size` is the family-wise size of the screen as executed, were every candidate to
#' carry a uniform effect of `c0`; at the claim threshold it is a boundary size.
#'
#' @param fit A [forestsearch()] fit (or a bare [fs_mr_inference()] result) with the
#'   declaration field retained:
#'   `mr_inference_args = list(keep_declaration_field = TRUE, declaration_c0 = ...)`.
#'   `declaration_c0` stores the shifted maxima for those `c0` levels; `keep_field_matrix
#'   = TRUE` lets them be computed for any.
#' @param alpha Target family-wise declaration rate. Default 0.10.
#' @param c0 `NULL` (the default, the claim-threshold row only) or a numeric vector of
#'   protected levels, on the natural scale of the consistency threshold (see
#'   [fs_declaration_calibration()]).
#' @param c0_primary Optional. The level, among `c0`, to mark as primary.
#'
#' @return A plain data frame with one row per `c0`, then the claim-threshold row:
#'   - `row`: `"c0"` or `"claim_threshold"`.
#'   - `c0`: the level (`NA` on the claim-threshold row).
#'   - `primary`: logical: the row is `c0_primary`.
#'   - `alpha`, as given.
#'   - `n_region`: the declared subgroup's size.
#'   - `T_region`, `T_max`: the statistic of the declared subgroup and the largest in the
#'     family.
#'   - `kappa_hat`: the calibrated cutoff.
#'   - `implied_level`: \eqn{2\Phi(\hat\kappa) - 1}.
#'   - `fw_size`: see above.
#'   - `n_admitted`: the candidates the calibrated rule admits at this row.
#'   - `admits_none`: logical: the calibrated rule admits nothing at any row.
#'   - `source`: where the shifted maxima came from, `"capture"` or `"field_matrix"` (`NA`
#'     on the claim-threshold row).
#'
#' @seealso [fs_declaration_calibration()].
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 500
#' d <- data.frame(id = seq_len(n), treat = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n))
#' d$time <- rexp(n, 0.1 * ifelse(d$treat == 1, ifelse(d$x1 < -0.6, 3, 0.45), 1))
#' d$event <- as.integer(d$time < 30)
#' d$time <- pmin(d$time, 30)
#' fit <- forestsearch(
#'   df.analysis = d, outcome.name = "time", event.name = "event", treat.name = "treat",
#'   id.name = "id", confounders.name = c("x1", "x2"), use_grf = FALSE, use_lasso = FALSE,
#'   maxk = 1, fs.splits = 200, sg_focus = "hrMaxSG", mr_inference = TRUE,
#'   mr_inference_args = list(ci_method = "ij", draws = 500, keep_declaration_field = TRUE,
#'                            declaration_c0 = c(0.6, 0.7)),
#'   parallel_args = list(plan = "sequential"), quiet = TRUE, details = FALSE,
#'   plot.sg = FALSE, show_candidate_summary = FALSE)
#' fs_declaration_table(fit, c0 = c(0.6, 0.7), c0_primary = 0.7)
#' }
#'
#' @export
fs_declaration_table <- function(fit, alpha = 0.10, c0 = NULL, c0_primary = NULL) {
  g <- if (inherits(fit, "forestsearch")) fit[["mr_inference"]] else fit
  fld <- if (is.list(g)) g[["declaration_field"]] else NULL
  if (is.null(fld))
    stop("fit carries no declaration field; fit with mr_inference = TRUE and ",
         "mr_inference_args = list(keep_declaration_field = TRUE).", call. = FALSE)
  if (!is.null(c0_primary) && (length(c0_primary) != 1L || !c0_primary %in% c0))
    stop("c0_primary must be one of the levels in c0.", call. = FALSE)

  # The reads of s17_column() in the ACTG 175 article.
  cal_c2 <- fs_declaration_calibration(fit, alpha = alpha)      # c0 = NULL: the claim threshold
  cal <- if (is.null(c0)) cal_c2 else fs_declaration_calibration(fit, alpha = alpha, c0 = c0)
  T_H <- cal$T_hat[[fld$family_id[fld$meta$selected_index]]]
  T_max <- max(cal$T_hat, na.rm = TRUE)
  tb <- if (is.null(c0)) NULL else cal$c0$table
  n_admitted <- c(if (!is.null(c0)) unname(lengths(cal$c0$admitted_calibrated)),
                  length(cal_c2$admitted_calibrated))
  kappa <- c(tb$kappa_hat, cal_c2$kappa_hat)
  data.frame(
    row = c(rep("c0", length(tb$c0)), "claim_threshold"),
    c0 = c(tb$c0, NA_real_),
    primary = c(if (is.null(c0_primary)) rep(FALSE, length(tb$c0)) else tb$c0 == c0_primary, FALSE),
    alpha = alpha,
    n_region = g[["n_selected"]],
    T_region = T_H,
    T_max = T_max,
    kappa_hat = kappa,
    implied_level = 2 * stats::pnorm(kappa) - 1,
    fw_size = c(tb$fw_size, cal_c2$fw_size),
    n_admitted = n_admitted,
    admits_none = all(n_admitted == 0L),
    source = c(rep(if (is.null(c0)) NA_character_ else cal$c0$source, length(tb$c0)), NA_character_),
    stringsAsFactors = FALSE)
}
