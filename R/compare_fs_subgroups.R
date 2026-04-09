# =============================================================================
# compare_fs_subgroups.R
#
# Compare subgroup membership between two forestsearch analyses.
# Produces a cross-tabulation of H / Hc membership with agreement
# statistics, formatted as a publication-ready gt table.
# =============================================================================

# Suppress R CMD check NOTE for gt column references used in cells_body()
utils::globalVariables("group")


#' Compare Subgroup Membership Across Two ForestSearch Analyses
#'
#' Cross-tabulates the subgroup assignments (\code{treat.recommend}) from
#' two \code{forestsearch} result objects and produces a formatted \code{gt}
#' summary table with concordance statistics.  Useful for assessing whether
#' different modelling approaches (e.g., Cox PH vs. Poisson rate model)
#' identify the same patients as belonging to the harm subgroup.
#'
#' @param fs1 A \code{forestsearch} result object (must contain
#'   \code{df.est$treat.recommend}).
#' @param fs2 A second \code{forestsearch} result object.
#' @param label1 Character. Display label for the first analysis.
#'   Default: \code{"Analysis 1"}.
#' @param label2 Character. Display label for the second analysis.
#'   Default: \code{"Analysis 2"}.
#' @param id.name Character. Name of the subject ID column in
#'   \code{df.est}.
#'   If \code{NULL} (default), row-position alignment is assumed (both
#'   \code{df.est} data frames must have identical row order and length).
#'   When specified, a merge on \code{id.name} is performed so that
#'   analyses with different row orders are handled correctly.
#' @param sg0_label Character. Label for H (harm / questionable benefit)
#'   subgroup (\code{treat.recommend == 0}).
#'   Default: \code{"H (Questionable)"}.
#' @param sg1_label Character. Label for Hc (recommend treatment)
#'   subgroup (\code{treat.recommend == 1}).
#'   Default: \code{"H\u1d9c (Recommend)"}.
#' @param title Character or \code{NULL}. Table title.
#'   Default: \code{"Subgroup Membership Comparison"}.
#' @param subtitle Character or \code{NULL}. Table subtitle.
#'   Default: auto-generated from \code{label1} and \code{label2}.
#' @param font_size Numeric. Base font size for the gt table.
#'   Default: 12.
#'
#' @return A list with components:
#'   \describe{
#'     \item{table}{A \code{gt} object ready for display or export.}
#'     \item{crosstab}{The underlying cross-tabulation matrix.}
#'     \item{concordance}{A list of agreement statistics: \code{agreement}
#'       (proportion), \code{kappa} (Cohen's kappa), \code{n_agree},
#'       \code{n_disagree}, \code{n_total}.}
#'     \item{membership}{A data frame with per-subject assignments
#'       from both analyses.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Compare Cox and Poisson ForestSearch analyses
#' result <- compare_fs_subgroups(
#'   fs_cox, fs_pois,
#'   label1 = "Cox PH (HR)",
#'   label2 = "Poisson (IRR)",
#'   id.name = "id"
#' )
#' result$table
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
    sg0_label = "H (Questionable)",
    sg1_label = "H\u1d9c (Recommend)",
    title     = "Subgroup Membership Comparison",
    subtitle  = NULL,
    font_size = 12
) {

  # ─── Input validation ────────────────────────────────────────────────────
  stopifnot(
    inherits(fs1, "forestsearch"),
    inherits(fs2, "forestsearch"),
    !is.null(fs1$df.est),
    !is.null(fs2$df.est),
    "treat.recommend" %in% names(fs1$df.est),
    "treat.recommend" %in% names(fs2$df.est)
  )

  df1 <- fs1$df.est
  df2 <- fs2$df.est

  # ─── Align subjects ─────────────────────────────────────────────────────
  if (!is.null(id.name)) {
    if (!id.name %in% names(df1) || !id.name %in% names(df2)) {
      stop("id.name '", id.name, "' not found in one or both df.est objects.",
           call. = FALSE)
    }
    merged <- merge(
      data.frame(id = df1[[id.name]], sg1 = df1$treat.recommend),
      data.frame(id = df2[[id.name]], sg2 = df2$treat.recommend),
      by = "id"
    )
    if (nrow(merged) == 0L) {
      stop("No matching IDs between the two analyses.", call. = FALSE)
    }
    if (nrow(merged) < nrow(df1) || nrow(merged) < nrow(df2)) {
      warning(
        sprintf("Merged %d of %d (analysis 1) and %d (analysis 2) subjects.",
                nrow(merged), nrow(df1), nrow(df2)),
        call. = FALSE
      )
    }
  } else {
    if (nrow(df1) != nrow(df2)) {
      stop(
        sprintf(
          paste0("Row counts differ (%d vs %d) and id.name is NULL.\n",
                 "  Provide id.name for merge-based alignment."),
          nrow(df1), nrow(df2)),
        call. = FALSE
      )
    }
    merged <- data.frame(
      id  = seq_len(nrow(df1)),
      sg1 = df1$treat.recommend,
      sg2 = df2$treat.recommend
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

  # Subgroup definitions
  sg_def1 <- if (!is.null(fs1$sg.harm)) {
    paste(fs1$sg.harm, collapse = " & ")
  } else {
    "(none identified)"
  }
  sg_def2 <- if (!is.null(fs2$sg.harm)) {
    paste(fs2$sg.harm, collapse = " & ")
  } else {
    "(none identified)"
  }

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
    subtitle <- sprintf("%s vs. %s", label1, label2)
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
