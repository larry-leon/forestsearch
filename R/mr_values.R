#' Multiplier-resampling (MR) inference quantities as a tidy data frame
#'
#' Extracts the MR quantities of a [forestsearch()] fit (its `mr_inference`
#' element) or of a bare [fs_mr_inference()] result into a plain data frame,
#' one row per region: `"H"` (the selected subgroup) and, when the complement
#' block exists, `"Hc"` (its complement). This is pure extraction: every value
#' already exists in `mr_inference` and nothing is recomputed, apart from the
#' optional `log()` of `scale = "log"`.
#'
#' Only the evaluated construction set is exposed: the naive estimate, the
#' two-term infinitesimal-jackknife (IJ) de-biased estimate, the field one-sided
#' bounds (lower on the selected subgroup, upper on its complement, plain and
#' studentized) and the Bonferroni pair, with the fit-level diagnostics.
#'
#' @param x A [forestsearch()] fit with a non-`NULL` `mr_inference` element, or
#'   the list returned by [fs_mr_inference()] (including its short variant,
#'   which yields an `"H"` row that is `NA` apart from `n_family`).
#' @param scale `"effect"` (the default) returns every column exactly as stored.
#'   `"log"` applies `log()` to the columns listed as *effect scale* below, for
#'   ratio measures (`HR`, `OR`, `RR`, `IRR`); it is an error for identity
#'   measures (`MD`, `RD`). The other columns are returned as stored under both
#'   settings.
#'
#' @return A plain `data.frame` with a fixed set of columns, whatever the fit
#'   computed: a quantity the fit did not compute is `NA`, and no column is ever
#'   dropped. The `"Hc"` row is present whenever `complement` is non-`NULL`; if
#'   the complement could not be fit, that row is all `NA`. Below, *working
#'   scale* means log for ratio measures and the effect itself for identity
#'   measures. Paths are relative to `mr_inference`.
#'
#'   - `row`: `"H"` or `"Hc"`.
#'   - `n`: Patients in the region: `n_selected` on `"H"`, `complement$n` on
#'     `"Hc"`. Count. `NA` in the short return, or when the complement could not
#'     be fit.
#'   - `naive_est`, `naive_lower`, `naive_upper`: The unadjusted estimate
#'     with its two-sided 95% robust interval: `naive$est`, `naive$lower`,
#'     `naive$upper`,
#'     and `complement$naive` on `"Hc"`. Effect scale. `NA` in the short return,
#'     or when the complement could not be fit.
#'   - `debiased_est`: The two-term de-biased estimate: `debiased$est`, and
#'     `complement$debiased$est` on `"Hc"`. Effect scale. `NA` as `naive_est`.
#'   - `ij_lower`, `ij_upper`: The two-sided 95% two-term IJ interval:
#'     `debiased$lower`, `debiased$upper`, and `complement$debiased` on `"Hc"`. Effect
#'     scale. Also `NA` under `ci_method = "wald"` or
#'     `ij_residual != "two_term"`, when those two stored limits are not this
#'     interval.
#'   - `ij_se`: The two-term IJ standard error: `debiased$se_ij_two_term`,
#'     and `complement$debiased$se_ij_two_term` on `"Hc"`. Working scale, under
#'     both settings of `scale`. It is the subgroup robust SE when the IJ
#'     variance is unusable (`debiased$ij_source` is `"wald_fallback"`). `NA` as
#'     `naive_est`.
#'   - `field_est2`: The field second-order point estimate: `field$est2`,
#'     and `field$complement$est2` on `"Hc"`. Effect scale. `NA` unless
#'     `ci_method = "field"`, and on `"Hc"` unless the complement field ran
#'     (`field_complement = TRUE` with `include_complement = TRUE`), or when
#'     that block returned a note instead of numbers.
#'   - `field_est2_s`: The studentized field-s counterpart on `"Hc"`:
#'     `field$complement$est2_s`. Effect scale. `NA` on `"H"`, and on `"Hc"` as
#'     `field_est2` or under `field_scale_complement = "none"`.
#'   - `field_lower_1s`: The field one-sided 95% lower bound on `"H"`:
#'     `field$lower_1s`. Effect scale. `NA` on `"Hc"`, and unless
#'     `ci_method = "field"`.
#'   - `field_upper_1s`: The field one-sided 95% upper bound on `"Hc"`:
#'     `field$complement$upper_1s`. Effect scale. `NA` on `"H"`, and on `"Hc"`
#'     as `field_est2`.
#'   - `field_upper_1s_s`: The field-s one-sided 95% upper bound on `"Hc"`:
#'     `field$complement$upper_1s_s`. Effect scale. `NA` as `field_est2_s`.
#'   - `bonf_lower_H`, `bonf_upper_Hc`: The Bonferroni pair (each bound at
#'     97.5%, jointly 95%), from the field's joint block: `field$joint$bonf_lower_H`
#'     on `"H"`, `field$joint$bonf_upper_Hc` on `"Hc"`. Effect scale. `NA` on the
#'     other row, and when `field$joint` is absent, which it is unless the
#'     complement field ran, or has no bounds.
#'   - `bonf_lower_H_s`, `bonf_upper_Hc_s`: The same pair from the field-s
#'     joint block, `field$joint_s`. Effect scale. `NA` on the other row, and
#'     when `field$joint_s` is absent, as for `bonf_lower_H`, or under
#'     `field_scale_complement = "none"`.
#'   - `selection_bias`: The selection-bias term: `selection_bias`. Working
#'     scale, under both settings of `scale`. Fit-level: `"H"` row only. `NA` in
#'     the short return.
#'   - `p_hat`: The selected candidate's re-selection frequency:
#'     `reselection$p_hat[[selected_index]]`. A probability. `"H"` row only. `NA`
#'     under `return_reselection = FALSE`, or in the short return.
#'   - `n_family`: The candidate family's size: `n_family`. Count. `"H"` row
#'     only.
#'   - `n_selected`: Patients in the selected subgroup: `n_selected`. Count.
#'     `"H"` row only. `NA` in the short return.
#'   - `timing_seconds`: Seconds in the main multiplier stage:
#'     `timing_seconds`. `"H"` row only. `NA` in the short return.
#'   - `field_timing_seconds`: Seconds in the field block: `field$timing_seconds`.
#'     `"H"` row only. `NA` unless `ci_method = "field"`.
#'   - `field_complement_timing_seconds`: Seconds in the complement field
#'     block: `field$complement$timing_seconds`. `"H"` row only. `NA` as
#'     `field_est2` on `"Hc"`.
#'
#' @seealso [forestsearch()] for the `mr_inference` switch;
#'   [fs_mr_inference()] for the computation; [mr_estimates_table()] to render
#'   the estimates as a table.
#'
#' @examples
#' set.seed(11)
#' n <- 300
#' d <- data.frame(time = rexp(n, 0.2), event = rbinom(n, 1, 0.8),
#'                 treat = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n))
#' family <- list(g1 = which(d$x1 > 0), g2 = which(d$x2 > 0),
#'                g3 = which(d$x1 + d$x2 > 0.3))
#' spec <- list(outcome_type = "survival", effect_measure = "HR",
#'              treat.name = "treat", outcome.name = "time",
#'              event.name = "event", offset.name = NULL,
#'              adjust_covariates = NULL, adverse_outcome = TRUE)
#' fit <- fs_mr_inference(
#'   df = d, candidates = family, spec = spec,
#'   selected_members = family$g1,
#'   admission = list(effect_floor = NULL, consistency = NULL),
#'   reselection = "maxeff", draws = 300L, seed = 4242L,
#'   ci_method = "field", include_complement = TRUE,
#'   field_R_out = 200L, field_R_in = 50L)
#' mr_values(fit)
#' mr_values(fit, scale = "log")[, c("row", "naive_est", "ij_se", "selection_bias")]
#'
#' @export
mr_values <- function(x, scale = c("effect", "log")) {
  scale <- match.arg(scale)
  g <- .mr_values_input(x)
  if (scale == "log" && isFALSE(.mr_get(g, "log_scale"))) {
    stop("scale = \"log\" applies to ratio measures only; this fit's measure (",
         .mr_get(g, "measure"), ") is on an identity scale.", call. = FALSE)
  }

  fb <- .mr_get(g, "field")
  cc <- .mr_get(g, "complement")
  fc <- .mr_get(fb, "complement")
  jt <- .mr_get(fb, "joint")
  js <- .mr_get(fb, "joint_s")

  # debiased$lower / $upper hold the two-term IJ interval only on this condition;
  # otherwise they are a robust-SE or a winner-variant interval, which is not
  # exposed.  fs_mr_inference.R: se (the harm side) and sec_used (the complement).
  residual <- .mr_get(g, "ij_residual")
  ci_method <- .mr_get(g, "ci_method")
  two_term_ci <- length(ci_method) == 1L && ci_method %in% c("field", "ij") &&
    (is.null(residual) || identical(residual, "two_term"))
  ij_limit <- function(block, limit) {
    if (two_term_ci) .mr_num(.mr_get(block, "debiased", limit)) else NA_real_
  }

  p_hat <- .mr_get(g, "reselection", "p_hat")
  sel <- .mr_get(g, "selected_index")
  p_hat_sel <- if (length(p_hat) && length(sel) == 1L && !is.na(sel)) {
    unname(p_hat[[sel]])
  } else {
    NA_real_
  }

  # The columns both rows read; each row's own columns are filled in below.
  region_row <- function(row, block, n, field_block) {
    data.frame(
      row = row,
      n = .mr_int(n),
      naive_est = .mr_num(.mr_get(block, "naive", "est")),
      naive_lower = .mr_num(.mr_get(block, "naive", "lower")),
      naive_upper = .mr_num(.mr_get(block, "naive", "upper")),
      debiased_est = .mr_num(.mr_get(block, "debiased", "est")),
      ij_lower = ij_limit(block, "lower"),
      ij_upper = ij_limit(block, "upper"),
      ij_se = .mr_num(.mr_get(block, "debiased", "se_ij_two_term")),
      field_est2 = .mr_num(.mr_get(field_block, "est2")),
      field_est2_s = NA_real_,
      field_lower_1s = NA_real_,
      field_upper_1s = NA_real_,
      field_upper_1s_s = NA_real_,
      bonf_lower_H = NA_real_,
      bonf_upper_Hc = NA_real_,
      bonf_lower_H_s = NA_real_,
      bonf_upper_Hc_s = NA_real_,
      selection_bias = NA_real_,
      p_hat = NA_real_,
      n_family = NA_integer_,
      n_selected = NA_integer_,
      timing_seconds = NA_real_,
      field_timing_seconds = NA_real_,
      field_complement_timing_seconds = NA_real_,
      stringsAsFactors = FALSE)
  }

  out <- region_row("H", g, .mr_get(g, "n_selected"), fb)
  out$field_lower_1s <- .mr_num(.mr_get(fb, "lower_1s"))
  out$bonf_lower_H <- .mr_num(.mr_get(jt, "bonf_lower_H"))
  out$bonf_lower_H_s <- .mr_num(.mr_get(js, "bonf_lower_H"))
  out$selection_bias <- .mr_num(.mr_get(g, "selection_bias"))
  out$p_hat <- p_hat_sel
  out$n_family <- .mr_int(.mr_get(g, "n_family"))
  out$n_selected <- .mr_int(.mr_get(g, "n_selected"))
  out$timing_seconds <- .mr_num(.mr_get(g, "timing_seconds"))
  out$field_timing_seconds <- .mr_num(.mr_get(fb, "timing_seconds"))
  out$field_complement_timing_seconds <- .mr_num(.mr_get(fc, "timing_seconds"))

  if (!is.null(cc)) {
    hc <- region_row("Hc", cc, .mr_get(cc, "n"), fc)
    hc$field_est2_s <- .mr_num(.mr_get(fc, "est2_s"))
    hc$field_upper_1s <- .mr_num(.mr_get(fc, "upper_1s"))
    hc$field_upper_1s_s <- .mr_num(.mr_get(fc, "upper_1s_s"))
    hc$bonf_upper_Hc <- .mr_num(.mr_get(jt, "bonf_upper_Hc"))
    hc$bonf_upper_Hc_s <- .mr_num(.mr_get(js, "bonf_upper_Hc"))
    out <- rbind(out, hc)
  }

  if (scale == "log" && isTRUE(.mr_get(g, "log_scale"))) {
    effect_cols <- c("naive_est", "naive_lower", "naive_upper", "debiased_est",
                     "ij_lower", "ij_upper", "field_est2", "field_est2_s",
                     "field_lower_1s", "field_upper_1s", "field_upper_1s_s",
                     "bonf_lower_H", "bonf_upper_Hc", "bonf_lower_H_s",
                     "bonf_upper_Hc_s")
    out[effect_cols] <- lapply(out[effect_cols], log)
  }
  out
}


#' The mr_inference list of a forestsearch fit or of a bare fs_mr_inference() result
#' @noRd
.mr_values_input <- function(x) {
  if (inherits(x, "forestsearch")) {
    g <- x[["mr_inference"]]
    if (is.null(g)) {
      stop("x is a forestsearch fit with no mr_inference results; run ",
           "forestsearch(mr_inference = TRUE). They are also absent when no ",
           "subgroup was identified or the MR step failed.", call. = FALSE)
    }
    return(g)
  }
  if (is.list(x) && all(c("selected_index", "settings") %in% names(x))) return(x)
  stop("x must be a forestsearch fit or a fs_mr_inference() result.", call. = FALSE)
}


#' Walk a path of list names with exact matching
#'
#' `$` matches partially, so reading `n` from a note-only block such as
#' `list(note = "...")` would return the note. A missing level gives `NULL`.
#' @noRd
.mr_get <- function(x, ...) {
  for (k in c(...)) {
    if (!is.list(x)) return(NULL)
    x <- x[[k]]
  }
  x
}


#' A stored scalar as numeric; NA_real_ when absent
#' @noRd
.mr_num <- function(x) if (length(x) == 1L) as.numeric(x) else NA_real_


#' A stored count as integer; NA_integer_ when absent
#' @noRd
.mr_int <- function(x) if (length(x) == 1L) as.integer(x) else NA_integer_
