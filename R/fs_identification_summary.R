# Identification-structure summaries for simulation bundles ------------------
#
# Factors the aggregation, figures and table that the manuscript supplement
# builds inline (fs_identification_figures*.qmd) into reusable package code.
# The three fragments differ only in payload paths and captions; the engine --
# an anchor/partner/proxy classification of each realized rule -- is identical,
# and is implemented once here.
#
# Two deliberate departures from the inline fragments:
#
#   1. Covariates are recovered from `sg_def` via .fs_rule_columns(), falling
#      back to the `covs` column only when `sg_def` is absent.  The fragments
#      read `covs` directly, which is empty by design in the GLM/continuous
#      harness -- strsplit(NA, "+") yields NA, every `has()` test is FALSE, and
#      every replicate silently classifies as "anchor missed".  Keying on
#      `sg_def` makes the summary work on any bundle that records a rule.
#
#   2. No scales/kableExtra dependency.  Percent labels are formatted directly
#      and the table is returned as a data frame (optionally rendered with gt,
#      already an Import).  A LaTeX document can still style the returned frame
#      with kableExtra without the package depending on it.

# --- internal helpers -------------------------------------------------------

#' Percent labels without a scales dependency
#' @keywords internal
#' @noRd
.fs_id_pct <- function(x, digits = 1) {
  ifelse(is.na(x), NA_character_,
         paste0(formatC(100 * x, format = "f", digits = digits), "%"))
}

#' Coerce one sweep element to a per-replicate results frame
#'
#' Accepts a data frame, a bundle list carrying `$results`, or a path to an
#' `.rds` holding either.
#'
#' @keywords internal
#' @noRd
.fs_id_read <- function(x) {
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) {
      stop("identification summary: file not found: ", x, call. = FALSE)
    }
    x <- readRDS(x)
  }
  if (is.data.frame(x)) {
    return(x)
  }
  if (is.list(x) && !is.null(x$results)) {
    return(x$results)
  }
  stop("identification summary: each element must be a data frame, a bundle ",
       "with a `results` element, or a path to an .rds holding either.",
       call. = FALSE)
}

#' Covariate sets for each replicate's realized rule
#'
#' @keywords internal
#' @noRd
.fs_id_covsets <- function(results, sg_def_col = "sg_def",
                           covs_col = "covs") {
  n <- nrow(results)
  sg <- if (!is.null(results[[sg_def_col]])) {
    as.character(results[[sg_def_col]])
  } else {
    rep(NA_character_, n)
  }
  cv <- if (!is.null(results[[covs_col]])) {
    as.character(results[[covs_col]])
  } else {
    rep(NA_character_, n)
  }
  if (all(is.na(sg)) && all(is.na(cv))) {
    stop("identification summary: neither `", sg_def_col, "` nor `", covs_col,
         "` carries a realized rule; nothing to classify.", call. = FALSE)
  }
  lapply(seq_len(n), function(i) {
    if (!is.na(sg[i]) && nzchar(sg[i])) {
      .fs_rule_columns(sg[i])
    } else if (!is.na(cv[i]) && nzchar(cv[i])) {
      trimws(strsplit(cv[i], "+", fixed = TRUE)[[1L]])
    } else {
      character(0)
    }
  })
}

#' Structure-label keys for one anchor/partner/proxy triple
#'
#' @keywords internal
#' @noRd
.fs_id_keys <- function(anchor, partner, proxy) {
  c(true  = sprintf("%s + %s (true)", anchor, partner),
    proxy = sprintf("%s + %s (proxy)", anchor, proxy),
    other = sprintf("%s + other", anchor),
    only  = sprintf("%s only", anchor),
    miss  = sprintf("anchor missed (no %s)", anchor))
}


# --- single-frame aggregation ------------------------------------------------

#' Identification structure for one set of replicates
#'
#' Classifies every detected replicate's realized rule against a planted
#' `anchor` / `partner` / `proxy` triple, and returns covariate involvement,
#' the structure composition, and classification accuracy.
#'
#' @section The classification:
#'
#' Each detected replicate falls in exactly one of five mutually exclusive
#' categories, tested in this order: the anchor is absent (`anchor missed`);
#' the rule names the anchor and nothing else (`anchor only`); the rule pairs
#' the anchor with the true `partner`; the rule pairs it with the `proxy`; or
#' the rule pairs it with something else. The order matters -- a rule naming
#' both the partner and the proxy counts as the true pairing.
#'
#' @param results A per-replicate results data frame, a bundle carrying
#'   `results`, or a path to an `.rds` holding either.
#' @param anchor Character. The covariate the search is expected to recover.
#' @param partner Character. The true second covariate of the planted region.
#' @param proxy Character. A correlated covariate that may substitute for
#'   `partner`.
#' @param confounders Character vector. The analysis covariate pool, for the
#'   involvement vector. Defaults to the three named roles.
#' @param sg_def_col,covs_col Character. Columns holding the realized rule and
#'   its projected covariate names.
#'
#' @return A list with `det_rate`, `n_detected`, `n_total`, `involvement`,
#'   `structure` (a named share vector), and the four accuracy means.
#'
#' @keywords internal
#' @export
fs_identification_structure <- function(results, anchor, partner, proxy,
                                        confounders = NULL,
                                        sg_def_col = "sg_def",
                                        covs_col = "covs") {
  stopifnot(is.character(anchor), length(anchor) == 1L,
            is.character(partner), length(partner) == 1L,
            is.character(proxy), length(proxy) == 1L)
  r <- .fs_id_read(results)
  if (is.null(confounders)) confounders <- unique(c(anchor, partner, proxy))

  if (is.null(r$detected)) {
    stop("identification summary: `results` has no `detected` column.",
         call. = FALSE)
  }
  det <- r[r$detected %in% 1L, , drop = FALSE]
  keys <- .fs_id_keys(anchor, partner, proxy)
  lv <- unname(keys)

  if (!nrow(det)) {
    zero <- stats::setNames(rep(0, length(lv)), lv)
    return(list(
      det_rate = 0, n_detected = 0L, n_total = nrow(r),
      involvement = stats::setNames(rep(NA_real_, length(confounders)),
                                    confounders),
      structure = zero,
      sens = NA_real_, spec = NA_real_, ppv = NA_real_, npv = NA_real_))
  }

  cl <- .fs_id_covsets(det, sg_def_col = sg_def_col, covs_col = covs_col)
  has <- function(cv) vapply(cl, function(z) cv %in% z, logical(1))
  only_anchor <- vapply(cl, function(z) length(setdiff(z, anchor)) == 0L,
                        logical(1))

  structure <- ifelse(
    !has(anchor), keys[["miss"]],
    ifelse(only_anchor, keys[["only"]],
           ifelse(has(partner), keys[["true"]],
                  ifelse(has(proxy), keys[["proxy"]], keys[["other"]]))))

  .m <- function(v) if (is.null(v)) NA_real_ else mean(v, na.rm = TRUE)
  list(
    det_rate = mean(r$detected %in% 1L),
    n_detected = nrow(det),
    n_total = nrow(r),
    involvement = vapply(confounders, function(cv) mean(has(cv)), numeric(1)),
    structure = prop.table(table(factor(structure, levels = lv))),
    sens = .m(det$sens), spec = .m(det$spec),
    ppv = .m(det$ppv), npv = .m(det$npv))
}


# --- sweep across trial sizes ------------------------------------------------

#' Identification structure across a sweep
#'
#' Applies [fs_identification_structure()] across a named collection of
#' bundles -- typically one per trial size -- and returns a single object the
#' plot and table helpers consume. This is the package-level replacement for
#' the inline `fs_identification_figures*.qmd` fragments.
#'
#' @param bundles A named list. Each element is a results data frame, a bundle
#'   carrying `results`, or a path to an `.rds`. Names label the sweep points
#'   (for example `"500"` and `"1000"`) and are used on the x axis in the
#'   order supplied.
#' @param anchor,partner,proxy Character. The planted-structure triple.
#' @param confounders Character vector. The analysis covariate pool.
#' @param sweep_label Character. Axis label for the sweep dimension.
#' @param sg_def_col,covs_col Character. Passed through.
#'
#' @return An object of class `fs_identification`: a list with `per_point`
#'   (the raw aggregates), tidy `involvement` and `structure` data frames, an
#'   `accuracy` data frame, and the triple it was built with.
#'
#' @examples
#' \dontrun{
#' id <- fs_identification_summary(
#'   list(`500` = "results/fs_..._n500_combined_1_500.rds",
#'        `1000` = "results/fs_..._n1000_combined_1_500.rds"),
#'   anchor = "er", partner = "meno", proxy = "age",
#'   confounders = c("er", "age", "meno", "pgr", "nodes", "size", "grade"))
#' plot_fs_identification_involvement(id)
#' plot_fs_identification_structure(id)
#' fs_identification_table(id)
#' }
#'
#' @export
fs_identification_summary <- function(bundles, anchor, partner, proxy,
                                      confounders = NULL,
                                      sweep_label = "Trial size (n)",
                                      sg_def_col = "sg_def",
                                      covs_col = "covs") {
  if (!is.list(bundles) || !length(bundles)) {
    stop("identification summary: `bundles` must be a non-empty list.",
         call. = FALSE)
  }
  if (is.null(names(bundles)) || any(!nzchar(names(bundles)))) {
    stop("identification summary: `bundles` must be NAMED; the names label ",
         "the sweep points (for example \"500\", \"1000\").", call. = FALSE)
  }
  if (is.null(confounders)) confounders <- unique(c(anchor, partner, proxy))

  per_point <- lapply(bundles, fs_identification_structure,
                      anchor = anchor, partner = partner, proxy = proxy,
                      confounders = confounders,
                      sg_def_col = sg_def_col, covs_col = covs_col)
  names(per_point) <- names(bundles)
  pts <- names(per_point)
  keys <- .fs_id_keys(anchor, partner, proxy)

  involvement <- do.call(rbind, lapply(pts, function(p) {
    a <- per_point[[p]]
    data.frame(point = p, covariate = names(a$involvement),
               involvement = unname(a$involvement),
               stringsAsFactors = FALSE)
  }))
  involvement$point <- factor(involvement$point, levels = pts)

  structure <- do.call(rbind, lapply(pts, function(p) {
    a <- per_point[[p]]
    data.frame(point = p, category = names(a$structure),
               share = as.numeric(a$structure), stringsAsFactors = FALSE)
  }))
  structure$point <- factor(structure$point, levels = pts)
  structure$category <- factor(structure$category, levels = unname(keys))

  accuracy <- do.call(rbind, lapply(pts, function(p) {
    a <- per_point[[p]]
    data.frame(point = p, n_total = a$n_total, n_detected = a$n_detected,
               det_rate = a$det_rate, sens = a$sens, spec = a$spec,
               ppv = a$ppv, npv = a$npv, stringsAsFactors = FALSE)
  }))

  out <- list(per_point = per_point, involvement = involvement,
              structure = structure, accuracy = accuracy,
              anchor = anchor, partner = partner, proxy = proxy,
              confounders = confounders, keys = keys,
              sweep_label = sweep_label, points = pts)
  class(out) <- c("fs_identification", "list")
  out
}

#' @export
print.fs_identification <- function(x, ...) {
  cat("<fs_identification>\n")
  cat(sprintf("  triple : anchor = %s | partner = %s | proxy = %s\n",
              x$anchor, x$partner, x$proxy))
  cat(sprintf("  sweep  : %s\n", paste(x$points, collapse = ", ")))
  for (p in x$points) {
    a <- x$per_point[[p]]
    cat(sprintf("  %-8s detection %s of %d; anchor recovered %s\n", p,
                .fs_id_pct(a$det_rate), a$n_total,
                .fs_id_pct(a$involvement[[x$anchor]])))
  }
  invisible(x)
}


# --- figures -----------------------------------------------------------------

#' Covariate involvement across the sweep
#'
#' Reproduces the supplement's involvement figure: the share of detections
#' whose rule names the anchor, the true partner, and the proxy, against the
#' sweep dimension.
#'
#' @param x An `fs_identification` object.
#' @param roles Character vector. Which covariates to draw; defaults to the
#'   anchor/partner/proxy triple.
#' @param palette Named character vector of colours, or `NULL` for the
#'   colour-blind-safe default used in the supplement.
#'
#' @return A `ggplot` object.
#' @importFrom rlang .data
#' @export
plot_fs_identification_involvement <- function(x, roles = NULL,
                                               palette = NULL) {
  stopifnot(inherits(x, "fs_identification"))
  if (is.null(roles)) roles <- c(x$anchor, x$partner, x$proxy)
  labs_role <- c(sprintf("%s (anchor)", x$anchor),
                 sprintf("%s (true partner)", x$partner),
                 sprintf("%s (proxy)", x$proxy))
  if (is.null(palette)) {
    palette <- stats::setNames(c("#000000", "#0072B2", "#D55E00"), labs_role)
  }
  d <- x$involvement[x$involvement$covariate %in% roles, , drop = FALSE]
  d$role <- factor(d$covariate, levels = roles,
                   labels = labs_role[seq_along(roles)])

  ggplot2::ggplot(d, ggplot2::aes(x = .data$point, y = .data$involvement,
                                  colour = .data$role, group = .data$role)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::scale_colour_manual(values = palette) +
    ggplot2::scale_y_continuous(labels = function(v) .fs_id_pct(v, 0),
                                limits = c(0, NA)) +
    ggplot2::labs(x = x$sweep_label, y = "Involvement (% of detections)",
                  colour = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "top")
}

#' Composition of the identified subgroup by structure
#'
#' Reproduces the supplement's stacked-composition figure.
#'
#' @param x An `fs_identification` object.
#' @param min_label Numeric. Shares below this are drawn without a label.
#' @param palette Named character vector of fills, or `NULL` for the default.
#'
#' @return A `ggplot` object.
#' @export
plot_fs_identification_structure <- function(x, min_label = 0.03,
                                             palette = NULL) {
  stopifnot(inherits(x, "fs_identification"))
  if (is.null(palette)) {
    palette <- stats::setNames(
      c("#2c7fb8", "#d95f02", "#7fcdbb", "#c7e9b4", "#bdbdbd"),
      unname(x$keys))
  }
  d <- x$structure
  ggplot2::ggplot(d, ggplot2::aes(x = .data$point, y = .data$share,
                                  fill = .data$category)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(.data$share >= min_label,
                                  .fs_id_pct(.data$share, 0), "")),
      position = ggplot2::position_stack(vjust = 0.5), size = 3,
      colour = "grey15") +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_y_continuous(labels = function(v) .fs_id_pct(v, 0)) +
    ggplot2::labs(x = x$sweep_label, y = "Share of detections", fill = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "right")
}


# --- table -------------------------------------------------------------------

#' Identification-structure table
#'
#' The supplement's structure table: detection rate, the five structure shares,
#' anchor recovery, and the conditional partner-versus-proxy split, one column
#' per sweep point. Returned as a data frame by default so a LaTeX document can
#' style it however it likes; `as_gt = TRUE` renders it with gt.
#'
#' @param x An `fs_identification` object.
#' @param as_gt Logical. Render with gt instead of returning the data frame.
#' @param digits Integer. Digits on the percent labels.
#'
#' @return A data frame, or a `gt_tbl`. The accuracy sentence the supplement
#'   carries as a table footnote is attached as the `"accuracy_note"`
#'   attribute of the data frame.
#' @export
fs_identification_table <- function(x, as_gt = FALSE, digits = 1) {
  stopifnot(inherits(x, "fs_identification"))
  k <- x$keys
  pts <- x$points

  split_ratio <- function(a) {
    num <- a$structure[[k[["true"]]]]
    den <- num + a$structure[[k[["proxy"]]]]
    if (!is.finite(den) || den <= 0) return("-")
    sprintf("%.0f : %.0f", 100 * num / den, 100 * (1 - num / den))
  }
  col_for <- function(p) {
    a <- x$per_point[[p]]
    c(.fs_id_pct(a$det_rate, digits),
      .fs_id_pct(a$structure[[k[["true"]]]], digits),
      .fs_id_pct(a$structure[[k[["proxy"]]]], digits),
      .fs_id_pct(a$structure[[k[["other"]]]], digits),
      .fs_id_pct(a$structure[[k[["only"]]]], digits),
      .fs_id_pct(a$structure[[k[["miss"]]]], digits),
      .fs_id_pct(a$involvement[[x$anchor]], digits),
      split_ratio(a))
  }

  tab <- data.frame(
    Quantity = c("Detection rate (all replicates)",
                 sprintf("%s + %s (true partner)", x$anchor, x$partner),
                 sprintf("%s + %s (proxy partner)", x$anchor, x$proxy),
                 sprintf("%s + other factor", x$anchor),
                 sprintf("%s only", x$anchor),
                 sprintf("anchor missed (no %s)", x$anchor),
                 sprintf("Anchor recovered (%s, any pairing)", x$anchor),
                 sprintf("Conditional partner split, %s : %s",
                         x$partner, x$proxy)),
    stringsAsFactors = FALSE, check.names = FALSE)
  for (p in pts) tab[[p]] <- col_for(p)

  acc <- x$accuracy
  fmt <- function(v) paste(sprintf("%.2f", v), collapse = " / ")
  attr(tab, "accuracy_note") <- sprintf(
    paste0("Classification of the identified subgroup against the planted ",
           "membership, averaged over detections (%s): sensitivity %s, ",
           "specificity %s, PPV %s, NPV %s."),
    paste(pts, collapse = " / "), fmt(acc$sens), fmt(acc$spec),
    fmt(acc$ppv), fmt(acc$npv))

  if (!isTRUE(as_gt)) return(tab)

  gt::tab_footnote(
    gt::tab_header(
      gt::gt(tab),
      title = "Identification structure across the sweep",
      subtitle = sprintf("anchor %s | true partner %s | proxy %s",
                         x$anchor, x$partner, x$proxy)),
    footnote = attr(tab, "accuracy_note"))
}


# --- realized-rule structure, collapsed over cut values ----------------------

#' Covariate signatures of the realized rules
#'
#' Collapses realized rules to the **set of covariates** they name, discarding
#' the cut values. A rules-by-frequency table lists every distinct `sg_def`
#' string and so is dominated by cut noise -- a 500-replicate survival run
#' produced 583 rows, in which `{age <= 46} & {nodes <= 5}` and
#' `{nodes <= 3} & {age <= 46}` are separate entries despite being the same
#' pairing. This reports how often each covariate *combination* is selected,
#' which is the structural question.
#'
#' @section What a signature is:
#'
#' The covariate names in the rule, de-duplicated and sorted, joined by `" & "`.
#' Sorting is what makes the two orderings above collapse together. Two
#' consequences worth knowing:
#'
#' * A rule naming one covariate twice --- a band such as
#'   `!{nodes <= 3} & {nodes <= 5}` --- has signature `"nodes"`, of size 1. The
#'   `k` column gives the number of distinct covariates, and `n_rules` the
#'   number of distinct rule strings that collapsed into the row, so bands are
#'   visible rather than hidden.
#' * Logically equivalent spellings of the same condition (`!{meno}` versus
#'   `{meno == 0}`) collapse together, which a raw rule table also splits.
#'
#' @param x A results data frame carrying `sg_def`, a bundle with `results`, a
#'   path to an `.rds`, or a bare character vector of rules.
#' @param true_covariates Character vector. The true-region covariates, used to
#'   classify each signature as `"exact"` (the signature is exactly this set),
#'   `"partial"` (overlaps it) or `"none"`. `NULL` skips the classification.
#' @param detected_only Logical. Restrict to rows with `detected == 1` when a
#'   results frame is supplied.
#'
#' @return A data frame ordered by frequency: `covariates`, `k`, `n`, `share`,
#'   `n_rules`, and `match` when `true_covariates` is given.
#'
#' @export
fs_rule_covariate_pairs <- function(x, true_covariates = NULL,
                                    detected_only = TRUE) {
  rules <- if (is.character(x)) {
    x
  } else {
    r <- .fs_id_read(x)
    if (isTRUE(detected_only) && !is.null(r$detected))
      r <- r[r$detected %in% 1L, , drop = FALSE]
    r$sg_def
  }
  rules <- rules[!is.na(rules) & nzchar(rules)]
  if (!length(rules)) {
    return(data.frame(covariates = character(0), k = integer(0),
                      n = integer(0), share = numeric(0),
                      n_rules = integer(0), stringsAsFactors = FALSE))
  }

  sig <- vapply(rules, function(rl) {
    cv <- .fs_rule_columns(rl)
    if (!length(cv)) return(NA_character_)
    paste(sort(unique(cv)), collapse = " & ")
  }, character(1), USE.NAMES = FALSE)

  keep <- !is.na(sig)
  sig <- sig[keep]; rules <- rules[keep]
  tab <- table(sig)
  out <- data.frame(
    covariates = names(tab),
    n = as.integer(tab),
    stringsAsFactors = FALSE)
  out$share <- out$n / length(sig)
  out$k <- vapply(strsplit(out$covariates, " & ", fixed = TRUE), length,
                  integer(1))
  out$n_rules <- vapply(out$covariates, function(s)
    length(unique(rules[sig == s])), integer(1), USE.NAMES = FALSE)

  if (!is.null(true_covariates)) {
    tset <- sort(unique(true_covariates))
    out$match <- vapply(out$covariates, function(s) {
      cs <- strsplit(s, " & ", fixed = TRUE)[[1L]]
      if (setequal(cs, tset)) "exact"
      else if (length(intersect(cs, tset))) "partial"
      else "none"
    }, character(1), USE.NAMES = FALSE)
  }
  out <- out[order(-out$n, out$covariates), c("covariates", "k", "n", "share",
                                              "n_rules",
                                              if (!is.null(true_covariates)) "match")]
  rownames(out) <- NULL
  out
}

#' Selected-covariate-pair frequencies
#'
#' Bar-chart companion to [fs_rule_covariate_pairs()], in the same idiom as the
#' selected-covariate figure: horizontal bars, ordered by share, coloured by
#' how the signature relates to the true region.
#'
#' @param x An `fs_identification` object, or anything
#'   [fs_rule_covariate_pairs()] accepts.
#' @param true_covariates Character vector; taken from `x` when it is an
#'   `fs_identification` object.
#' @param top_n Integer. Show only the most frequent `top_n` signatures; the
#'   remainder are pooled into a single "other" bar so the shares still sum to
#'   1 and nothing is silently dropped.
#' @param detected_only Passed through.
#'
#' @return A `ggplot` object.
#' @importFrom rlang .data
#' @export
plot_fs_rule_pairs <- function(x, true_covariates = NULL, top_n = 12L,
                               detected_only = TRUE) {
  if (inherits(x, "fs_identification")) {
    if (is.null(true_covariates))
      true_covariates <- c(x$anchor, x$partner)
    stop("plot_fs_rule_pairs(): pass the results frame or bundle, not the ",
         "fs_identification object -- signatures are built from `sg_def`, ",
         "which that object does not retain.", call. = FALSE)
  }
  d <- fs_rule_covariate_pairs(x, true_covariates = true_covariates,
                               detected_only = detected_only)
  if (!nrow(d)) return(ggplot2::ggplot() + ggplot2::theme_void())

  if (nrow(d) > top_n) {
    head_d <- d[seq_len(top_n), , drop = FALSE]
    tail_d <- d[-seq_len(top_n), , drop = FALSE]
    pooled <- data.frame(
      covariates = sprintf("other (%d signatures)", nrow(tail_d)),
      k = NA_integer_, n = sum(tail_d$n), share = sum(tail_d$share),
      n_rules = sum(tail_d$n_rules), stringsAsFactors = FALSE)
    if (!is.null(true_covariates)) pooled$match <- "none"
    d <- rbind(head_d, pooled[, names(head_d), drop = FALSE])
  }

  d$covariates <- factor(d$covariates, levels = rev(d$covariates))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$covariates, y = .data$share))
  if (!is.null(true_covariates)) {
    p <- p + ggplot2::geom_col(ggplot2::aes(fill = .data$match)) +
      ggplot2::scale_fill_manual(
        values = c(exact = "#b3261e", partial = "#e0a030", none = "grey60"),
        breaks = c("exact", "partial", "none"),
        labels = c(exact = "exactly the true pair", partial = "overlaps it",
                   none = "neither"),
        name = NULL)
  } else {
    p <- p + ggplot2::geom_col(fill = "grey50")
  }
  p + ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(labels = function(v) .fs_id_pct(v, 0)) +
    ggplot2::labs(x = NULL, y = "Share of detected replicates",
                  title = "Selected-covariate-pair frequencies",
                  subtitle = "Realized rules collapsed over cut values") +
    ggplot2::theme_minimal(base_size = 12)
}
