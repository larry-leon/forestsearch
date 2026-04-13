# =============================================================================
# compare_fs_subgroups.R
#
# Compare subgroup membership between two forestsearch analyses.
# Produces a cross-tabulation of H / Hc membership with agreement
# statistics, formatted as a publication-ready gt table.
# =============================================================================

# Suppress R CMD check NOTE for gt column references used in cells_body()
utils::globalVariables("group")


#' Compare Subgroup Membership Across Two Analyses
#'
#' Cross-tabulates subgroup assignments from two analyses and produces
#' a formatted \code{gt} summary table with concordance statistics.
#' Accepts \code{forestsearch} result objects, GRF result lists, or raw
#' indicator vectors.  Non-detection (no subgroup found) is handled
#' automatically: all subjects are classified to the complement.
#'
#' @param fs1 Either a \code{forestsearch} result object, a
#'   \code{grf.subg.harm} result (with \code{data$treat.recommend}),
#'   or an integer/logical vector of subgroup indicators
#'   (0 = in H/Q, 1 = in Hc/Qc).
#' @param fs2 Same types as \code{fs1} for the second analysis.
#' @param label1 Character. Display label for the first analysis.
#' @param label2 Character. Display label for the second analysis.
#' @param id.name Character. Subject ID column for merge alignment.
#'   Only used when both inputs are result objects.
#'   If \code{NULL}, row-position alignment is assumed.
#' @param sg0_label Character or \code{NULL}. Label for H/Q
#'   (\code{treat.recommend == 0}).  If \code{NULL}, auto-generated
#'   from \code{subgroup_notation}.
#' @param sg1_label Character or \code{NULL}. Label for Hc/Qc
#'   (\code{treat.recommend == 1}).  If \code{NULL}, auto-generated.
#' @param subgroup_notation Character. \code{"harm"} uses H/Hc labels;
#'   \code{"benefit"} uses Q/Qc labels.  Only used when
#'   \code{sg0_label}/\code{sg1_label} are \code{NULL}.
#' @param title Character. Table title.
#' @param subtitle Character or \code{NULL}. Auto-generated if NULL.
#' @param font_size Numeric. Font size for gt table.
#'
#' @return A list with:
#'   \describe{
#'     \item{table}{A \code{gt} cross-tabulation table.}
#'     \item{crosstab}{The raw cross-tabulation matrix.}
#'     \item{concordance}{Agreement, kappa, counts.}
#'     \item{membership}{Per-subject assignments from both analyses.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Harm search: Cox vs Poisson (both forestsearch objects)
#' result <- compare_fs_subgroups(fs_cox, fs_pois,
#'   label1 = "Cox PH", label2 = "Poisson",
#'   subgroup_notation = "harm")
#'
#' # Benefit search: FS vs GRF (handles non-detection automatically)
#' result <- compare_fs_subgroups(fs_result, grf_result,
#'   label1 = "ForestSearch", label2 = "GRF",
#'   subgroup_notation = "benefit")
#'
#' # Raw vectors
#' result <- compare_fs_subgroups(
#'   c(rep(0L, 30), rep(1L, 70)),
#'   c(rep(0L, 25), rep(1L, 75)),
#'   label1 = "Method A", label2 = "Method B")
#' }
#'
#' @importFrom gt gt tab_header tab_spanner cols_label fmt_number
#'   tab_source_note tab_footnote tab_style cell_text cells_body
#'   cells_column_labels
#' @export
compare_fs_subgroups <- function(
    fs1,
    fs2,
    label1    = "Analysis 1",
    label2    = "Analysis 2",
    id.name   = NULL,
    sg0_label = NULL,
    sg1_label = NULL,
    subgroup_notation = c("harm", "benefit"),
    title     = "Subgroup Membership Comparison",
    subtitle  = NULL,
    font_size = 12
) {

  subgroup_notation <- match.arg(subgroup_notation)

  # ─── Default labels from notation ──────────────────────────────────────
  if (is.null(sg0_label)) {
    sg0_label <- if (subgroup_notation == "harm")
      "H (Questionable)" else "Q (Benefit)"
  }
  if (is.null(sg1_label)) {
    sg1_label <- if (subgroup_notation == "harm")
      paste0("H", "\u1d9c", " (Recommend)") else
      paste0("Q", "\u1d9c", " (No benefit)")
  }

  # ─── Extract indicator vectors from various input types ────────────────
  .extract_sg <- function(obj, obj_label) {
    # Raw vector
    if (is.numeric(obj) || is.logical(obj)) {
      return(list(sg = as.integer(obj), ids = NULL, sg_def = NULL,
                  detected = any(as.integer(obj) == 0L)))
    }
    # forestsearch object
    if (inherits(obj, "forestsearch")) {
      sg_def <- if (!is.null(obj$sg.harm) && length(obj$sg.harm) > 0)
        paste(obj$sg.harm, collapse = " & ") else "(none identified)"

      if (!is.null(obj$df.est) &&
          "treat.recommend" %in% names(obj$df.est)) {
        ids <- if (!is.null(id.name) && id.name %in% names(obj$df.est))
          obj$df.est[[id.name]] else NULL
        return(list(sg = obj$df.est$treat.recommend, ids = ids,
                    sg_def = sg_def, detected = any(obj$df.est$treat.recommend == 0L)))
      }
      # df.est is NULL → no subgroup detected → all Hc
      # N will be inferred from the other input in alignment
      return(list(sg = NULL, ids = NULL, sg_def = sg_def, detected = FALSE))
    }
    # GRF result (list with data$treat.recommend)
    if (is.list(obj) && !is.null(obj$data) &&
        "treat.recommend" %in% names(obj$data)) {
      ids <- if (!is.null(id.name) && id.name %in% names(obj$data))
        obj$data[[id.name]] else NULL
      sg_def <- if (!is.null(obj$sg.harm.id) && length(obj$sg.harm.id) > 0)
        paste(obj$sg.harm.id, collapse = " & ") else "(none identified)"
      return(list(sg = obj$data$treat.recommend, ids = ids,
                  sg_def = sg_def,
                  detected = any(obj$data$treat.recommend == 0L)))
    }
    stop(sprintf(
      paste0("%s must be a forestsearch object, a GRF result with ",
             "data$treat.recommend, or a numeric/logical vector."),
      obj_label
    ), call. = FALSE)
  }

  ex1 <- .extract_sg(fs1, "fs1")
  ex2 <- .extract_sg(fs2, "fs2")

  # ─── Handle non-detection: fill NULL sg with all-complement ────────────
  if (is.null(ex1$sg) && is.null(ex2$sg))
    stop("Both inputs have no subgroup assignments.", call. = FALSE)

  if (is.null(ex1$sg)) {
    n <- length(ex2$sg)
    ex1$sg <- rep(1L, n)
  }
  if (is.null(ex2$sg)) {
    n <- length(ex1$sg)
    ex2$sg <- rep(1L, n)
  }

  sg_def1 <- if (!is.null(ex1$sg_def)) ex1$sg_def else "(vector input)"
  sg_def2 <- if (!is.null(ex2$sg_def)) ex2$sg_def else "(vector input)"

  # ─── Align subjects ─────────────────────────────────────────────────────
  if (!is.null(id.name) && !is.null(ex1$ids) && !is.null(ex2$ids)) {
    merged <- merge(
      data.frame(id = ex1$ids, sg1 = ex1$sg),
      data.frame(id = ex2$ids, sg2 = ex2$sg),
      by = "id"
    )
    if (nrow(merged) == 0L)
      stop("No matching IDs between the two analyses.", call. = FALSE)
    if (nrow(merged) < length(ex1$sg) || nrow(merged) < length(ex2$sg))
      warning(sprintf("Merged %d of %d / %d subjects.",
                      nrow(merged), length(ex1$sg), length(ex2$sg)),
              call. = FALSE)
  } else {
    if (length(ex1$sg) != length(ex2$sg))
      stop(sprintf("Lengths differ (%d vs %d). Provide id.name or equal-length vectors.",
                   length(ex1$sg), length(ex2$sg)), call. = FALSE)
    merged <- data.frame(
      id  = seq_along(ex1$sg),
      sg1 = ex1$sg,
      sg2 = ex2$sg
    )
  }

  # ─── Subgroup labels ────────────────────────────────────────────────────
  sg_levels <- c(0L, 1L)
  sg_labels <- c(sg0_label, sg1_label)

  merged$sg1_lab <- factor(merged$sg1, levels = sg_levels, labels = sg_labels)
  merged$sg2_lab <- factor(merged$sg2, levels = sg_levels, labels = sg_labels)

  # ─── Cross-tabulation ───────────────────────────────────────────────────
  ct <- table(merged$sg1_lab, merged$sg2_lab,
              dnn = c(label1, label2))
  ct_mat <- as.matrix(ct)

  n_total    <- sum(ct_mat)
  n_agree    <- sum(diag(ct_mat))
  n_disagree <- n_total - n_agree
  agreement  <- n_agree / n_total

  # Cohen's kappa
  p_obs <- agreement
  row_marg <- rowSums(ct_mat) / n_total
  col_marg <- colSums(ct_mat) / n_total
  p_exp <- sum(row_marg * col_marg)
  kappa <- if (abs(1 - p_exp) < 1e-10) 1.0 else (p_obs - p_exp) / (1 - p_exp)

  concordance <- list(
    agreement  = agreement,
    kappa      = kappa,
    n_agree    = n_agree,
    n_disagree = n_disagree,
    n_total    = n_total
  )

  # ─── Build gt display table ─────────────────────────────────────────────

  # Cross-tab as data frame for gt
  ct_df <- as.data.frame.matrix(ct_mat)
  ct_df$Total <- rowSums(ct_mat)
  ct_df <- rbind(ct_df, c(colSums(ct_mat), n_total))
  rownames(ct_df)[nrow(ct_df)] <- "Total"

  # Add row label column
  gt_df <- data.frame(
    group = rownames(ct_df),
    ct_df,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  rownames(gt_df) <- NULL

  # Auto subtitle
  if (is.null(subtitle)) {
    det1 <- if (ex1$detected) sg_def1 else "(no detection)"
    det2 <- if (ex2$detected) sg_def2 else "(no detection)"
    subtitle <- sprintf("%s: %s vs. %s: %s", label1, det1, label2, det2)
  }

  # Build gt table
  tbl <- gt::gt(gt_df) |>
    gt::tab_header(
      title    = gt::md(paste0("**", title, "**")),
      subtitle = subtitle
    ) |>
    gt::cols_label(group = label1) |>
    gt::tab_spanner(
      label   = label2,
      columns = sg_labels
    ) |>
    gt::tab_style(
      style     = gt::cell_text(weight = "bold"),
      locations = gt::cells_column_labels()
    ) |>
    # Bold the Total row
    gt::tab_style(
      style     = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(
        rows = group == "Total"
      )
    ) |>
    # Bold the Total column
    gt::tab_style(
      style     = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(columns = "Total")
    ) |>
    # Bold the row label column
    gt::tab_style(
      style     = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(columns = "group")
    ) |>
    # Highlight concordant cells (diagonal) with light green
    gt::tab_style(
      style     = list(
        gt::cell_fill(color = "#d4edda"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = sg0_label,
        rows    = group == sg0_label
      )
    ) |>
    gt::tab_style(
      style     = list(
        gt::cell_fill(color = "#d4edda"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = sg1_label,
        rows    = group == sg1_label
      )
    ) |>
    # Highlight discordant cells (off-diagonal) with light red if > 0
    gt::tab_style(
      style     = list(
        gt::cell_fill(color = "#f8d7da"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = sg1_label,
        rows    = group == sg0_label
      )
    ) |>
    gt::tab_style(
      style     = list(
        gt::cell_fill(color = "#f8d7da"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = sg0_label,
        rows    = group == sg1_label
      )
    ) |>
    gt::tab_source_note(
      source_note = gt::md(sprintf(
        paste0(
          "**Agreement:** %d / %d (%.1f%%)
",
          "**Cohen's \u03ba:** %.3f
",
          "**%s H:** %s
",
          "**%s H:** %s"
        ),
        n_agree, n_total, 100 * agreement,
        kappa,
        label1, sg_def1,
        label2, sg_def2
      ))
    ) |>
    gt::tab_options(
      table.font.size   = font_size,
      heading.align      = "left",
      source_notes.font.size = gt::px(font_size - 1)
    )

  # ─── Return ─────────────────────────────────────────────────────────────
  list(
    table       = tbl,
    crosstab    = ct,
    concordance = concordance,
    membership  = merged
  )
}
