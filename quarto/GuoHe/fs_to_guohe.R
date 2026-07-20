# fs_to_guohe.R
#
# Bridge from a `forestsearch(sg_focus = "maxeff")` fit to the Guo & He (2021)
# Algorithm 3 de-biasing estimator (`guohe_algorithm3.R`). This is the
# "post-hoc identified subgroups" application of Guo & He: forestsearch
# enumerates a large candidate family and selects the argmax; Guo & He supplies
# a de-biased effect and one-sided bound for the SELECTED subgroup, correcting
# the winner's-curse optimism of that selection.
#
# WHY A PURE POST-PROCESSOR. The full deduplicated maxeff candidate family
# survives on the returned fit at `fit$grp.consistency$out_sg$result`, one row
# per candidate, untruncated under `sg_focus = "maxeff"`. No re-enumeration and
# no package change are required.
#
# THE REAL FAMILY SCHEMA (verified by execution against a GBSG maxeff fit; the
# adapter is written to THIS, not to an assumed indicator-matrix layout):
#
#   out_sg$result columns: Pcons, hr, N, E, g, m, K, M.1, M.2, ...
#     * K            number of cuts (1 or 2 for maxk = 2)
#     * M.1 ... M.K  cut LABELS as character strings, e.g. "{er <= 0}",
#                    "{size <= 35}", or a negated cut "!{size <= 20}".
#     * M.j is "" (empty) for j > K (e.g. M.2 == "" when K == 1).
#
# The labels are NOT 0/1 indicator columns. Membership must be reconstructed:
# strip the outer braces and an optional leading "!" from each label, match the
# inner string to `fit$confounders.evaluated` (the brace-free label vector),
# read the aligned dummy column name from `fit$confounders.candidate` (q1..qJ
# on `fit$df.est`), and take that dummy (or its complement 1 - dummy when the
# label was negated). A candidate's membership is the intersection (AND) over
# its K cuts.
#
# MEMBERSHIP MUST MATCH FORESTSEARCH EXACTLY. forestsearch already computed the
# selected subgroup's per-subject membership as
# `fit$grp.consistency$sg.harm.id` (length nrow(df.est); NOT `fit$sg.harm.id`,
# which is empty). The adapter's materialized column for the selected candidate
# must equal it row for row. `fs_assert_membership()` enforces this as a hard
# error; a mismatch means the reconstruction is wrong and the run stops.
#
# Depends on {survival} (via guohe_algorithm3.R) and base R. Not part of the
# package surface. Source alongside guohe_algorithm3.R.

# ---------------------------------------------------------------------------
# label parsing
# ---------------------------------------------------------------------------

#' Parse one forestsearch cut label into a dummy-column reference
#'
#' @param lbl A cut label, e.g. "{er <= 0}", "!{size <= 20}", or "" (absent cut).
#' @param cand `fit$confounders.candidate` (dummy column names, q1..qJ).
#' @param evaluated `fit$confounders.evaluated` (brace-free labels, aligned to
#'   `cand`).
#' @return `NULL` for an empty label; otherwise a list with `dummy` (the column
#'   name) and `negate` (logical).
.fs_gh_parse_label <- function(lbl, cand, evaluated) {
  lbl <- trimws(lbl)
  if (!nzchar(lbl)) return(NULL)
  negate <- startsWith(lbl, "!")
  if (negate) lbl <- sub("^!", "", lbl)
  inner <- sub("^\\{", "", sub("\\}$", "", lbl))
  inner <- trimws(inner)
  j <- match(inner, evaluated)
  if (is.na(j)) {
    stop("cut label '", inner, "' not found in fit$confounders.evaluated; ",
         "cannot map to a dummy column.")
  }
  list(dummy = cand[j], negate = negate)
}

#' Membership (0/1 vector) for one candidate row of the family table
.fs_gh_candidate_membership <- function(row, df, cand, evaluated, mcols) {
  labels <- vapply(mcols, function(cn) as.character(row[[cn]]), character(1))
  parsed <- Filter(Negate(is.null),
                   lapply(labels, .fs_gh_parse_label, cand = cand,
                          evaluated = evaluated))
  if (!length(parsed)) return(integer(nrow(df)))
  cols <- lapply(parsed, function(p) {
    v <- as.integer(as.character(df[[p$dummy]])) == 1L
    if (isTRUE(p$negate)) !v else v
  })
  as.integer(Reduce(`&`, cols))
}

# ---------------------------------------------------------------------------
# family export
# ---------------------------------------------------------------------------

#' Materialize the maxeff candidate family as 0/1 membership columns
fs_export_maxeff_family <- function(fit, prefix = "sg_") {
  if (!inherits(fit, "forestsearch")) stop("fit must be a 'forestsearch' object.")
  if (!identical(fit$sg_focus, "maxeff")) {
    stop("fs_export_maxeff_family() requires sg_focus = 'maxeff'; got '",
         fit$sg_focus, "'.")
  }
  gc <- fit$grp.consistency
  if (is.null(gc) || is.null(gc$out_sg) || is.null(gc$out_sg$result)) {
    stop("No maxeff candidate family (grp.consistency$out_sg$result is NULL).")
  }
  fam <- as.data.frame(gc$out_sg$result)
  if (!nrow(fam)) stop("The maxeff candidate family is empty.")

  cand      <- fit$confounders.candidate
  evaluated <- fit$confounders.evaluated
  if (is.null(cand) || is.null(evaluated) || length(cand) != length(evaluated)) {
    stop("fit$confounders.candidate / $confounders.evaluated missing or ",
         "misaligned; cannot map labels to dummy columns.")
  }

  df <- as.data.frame(fit$df.est)
  miss <- setdiff(cand, names(df))
  if (length(miss)) {
    stop("Dummy columns absent from fit$df.est: ",
         paste(utils::head(miss, 10L), collapse = ", "),
         if (length(miss) > 10L) ", ..." else "", ".")
  }

  mcols <- grep("^M\\.[0-9]+$", names(fam), value = TRUE)
  if (!length(mcols)) stop("No cut-label columns (M.1, M.2, ...) in the family.")

  n <- nrow(df)
  cand_names <- sprintf("%s%04d", prefix, seq_len(nrow(fam)))
  mem <- matrix(0L, nrow = n, ncol = nrow(fam),
                dimnames = list(NULL, cand_names))
  for (i in seq_len(nrow(fam))) {
    mem[, i] <- .fs_gh_candidate_membership(fam[i, , drop = FALSE], df,
                                            cand, evaluated, mcols)
  }

  out_df <- cbind(df, as.data.frame(mem))
  selected <- cand_names[1L]  # out_sg$result is sorted selected-first
  list(data = out_df, candidates = cand_names, selected = selected, family = fam)
}

# ---------------------------------------------------------------------------
# membership integrity check
# ---------------------------------------------------------------------------

#' Assert the exported selected-subgroup membership matches forestsearch's own
fs_assert_membership <- function(exp, fit) {
  own <- fit$grp.consistency$sg.harm.id
  if (is.null(own)) own <- fit$grp.consistency$out_sg$sg.harm.id
  if (is.null(own)) stop("grp.consistency$sg.harm.id is NULL; cannot verify.")
  own <- as.integer(own)
  mine <- as.integer(exp$data[[exp$selected]])
  if (length(own) != length(mine)) {
    stop("Membership length mismatch: sg.harm.id has ", length(own),
         " rows, exported column has ", length(mine), ".")
  }
  d <- sum(own != mine)
  if (d != 0L) {
    stop("Exported selected-subgroup membership disagrees with ",
         "grp.consistency$sg.harm.id in ", d, " of ", length(own), " rows. ",
         "The reconstructed family does not match forestsearch's selection; ",
         "refusing to proceed.")
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# top-level bridge
# ---------------------------------------------------------------------------

#' Run Guo & He Algorithm 3 on a forestsearch maxeff fit
fs_to_guohe <- function(fit,
                        time = "time_months", event = "status", treatment = "hormon",
                        B = 2000L, r = 0.03, level = 0.05, seed = NULL,
                        min_events = 5L, parallel = FALSE, verify = TRUE) {
  if (!exists("guohe_algorithm3", mode = "function")) {
    stop("guohe_algorithm3() not found; source guohe_algorithm3.R first.")
  }
  exp <- fs_export_maxeff_family(fit)
  if (isTRUE(verify)) fs_assert_membership(exp, fit)

  need <- c(time, event, treatment)
  miss <- setdiff(need, names(exp$data))
  if (length(miss)) {
    stop("Columns absent from fit$df.est: ", paste(miss, collapse = ", "),
         ". Set time/event/treatment to match the forestsearch frame.")
  }

  # orient = +1: forestsearch maxeff selects the MAX-HR (harm) subgroup, which is
  # argmax(orient * logHR) under orient = +1 -- so Guo & He's argmax coincides
  # with the forestsearch selection and corrects that (harm) subgroup. (orient =
  # -1 would select argmax(-logHR), i.e. the most PROTECTIVE subgroup, the wrong
  # target.) Matches guohe_from_forestsearch.R's harm convention.
  gh <- guohe_algorithm3(
    data = exp$data, outcome = "survival", treatment = treatment,
    candidates = exp$candidates, time = time, event = event,
    orient = +1, B = B, r = r, level = level, seed = seed,
    min_events = min_events, parallel = parallel, diagnostics = TRUE
  )

  if (!identical(gh$selected, exp$selected)) {
    warning("Guo & He argmax (", gh$selected, ") differs from the ",
            "forestsearch selection (", exp$selected, "). The returned Guo & He ",
            "objects refer to ", gh$selected, ".")
  }

  list(gh = gh, export = exp, selected = gh$selected,
       naive_hr = unname(gh$naive), debiased_hr = unname(gh$debiased),
       bound_hr = unname(gh$bound_one_sided))
}
