#' Per-Subject Membership Matrix for Pareto Frontier Members
#'
#' For each candidate subgroup on
#' \code{fs$grp.consistency$out_sg$pareto_frontier}, returns a binary
#' indicator (one row per subject in \code{fs$df.est}, one column per
#' frontier member) marking which subjects belong to that subgroup.
#'
#' This lets the user compute downstream quantities the frontier table
#' does not surface directly: per-frontier-member sample sizes by
#' covariate stratum, set comparisons between frontier members ("which
#' subjects are in the selected subgroup but not on the runner-up?"),
#' and survival summaries by frontier member.
#'
#' @param fs A \code{forestsearch} result object.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{\code{flags}}{Integer matrix of dimension
#'       \code{n_subjects} x \code{n_frontier_members}, with 1 = subject
#'       in subgroup, 0 = not.}
#'     \item{\code{member_labels}}{Character vector of length
#'       \code{n_frontier_members}, each entry the
#'       "\code{cut1 & cut2 & ...}" definition of that frontier member.}
#'     \item{\code{is_selected}}{Logical vector of length
#'       \code{n_frontier_members}, marking which member is the selected
#'       subgroup.}
#'   }
#'
#' @examples
#' \dontrun{
#' fmf <- frontier_member_flags(fs)
#' # How many subjects in each frontier member?
#' colSums(fmf$flags)
#' # Which subjects are unique to the selected subgroup?
#' sel <- which(fmf$is_selected)
#' rowSums(fmf$flags) == 1L & fmf$flags[, sel] == 1L
#' }
#'
#' @seealso \code{\link{pareto_frontier_table}},
#'   \code{\link{compute_frontier_cis}}, \code{\link{plot_pareto_frontier}}.
#'
#' @importFrom data.table is.data.table
#' @export
frontier_member_flags <- function(fs) {

  out_sg   <- tryCatch(fs$grp.consistency$out_sg, error = function(e) NULL)
  frontier <- tryCatch(out_sg$pareto_frontier,    error = function(e) NULL)
  if (is.null(frontier) || !data.table::is.data.table(frontier) ||
      nrow(frontier) == 0L) {
    message("No Pareto frontier available on this forestsearch object.")
    return(invisible(NULL))
  }

  df.est <- tryCatch(fs$df.est, error = function(e) NULL)
  if (is.null(df.est) || !nrow(df.est)) {
    stop("frontier_member_flags(): fs$df.est is missing or empty.",
         call. = FALSE)
  }

  m_cols    <- grep("^M\\.", names(frontier), value = TRUE)
  K         <- nrow(frontier)
  n_subj    <- nrow(df.est)
  flags_mat <- matrix(0L, nrow = n_subj, ncol = K)
  labels    <- character(K)

  for (k in seq_len(K)) {
    row_k <- frontier[k, ]
    cuts  <- unlist(row_k[, m_cols, with = FALSE], use.names = FALSE)
    cuts  <- cuts[!is.na(cuts) & nzchar(cuts)]
    labels[k] <- if (length(cuts) == 0L) "<empty>" else
                 paste(cuts, collapse = " & ")

    keep <- rep(TRUE, n_subj)
    missing_any <- FALSE
    for (c_name in cuts) {
      if (!c_name %in% names(df.est)) { missing_any <- TRUE; break }
      col <- df.est[[c_name]]
      keep <- keep & !is.na(col) & col == 1L
    }
    if (!missing_any) flags_mat[keep, k] <- 1L
  }

  selected_m <- tryCatch(
    as.integer(out_sg$result[1, ]$m),
    error = function(e) NA_integer_
  )
  ft_m <- as.integer(frontier[["m"]])
  is_selected_vec <- !is.na(selected_m) & !is.na(ft_m) & ft_m == selected_m

  colnames(flags_mat) <- paste0("m", ft_m)

  list(
    flags         = flags_mat,
    member_labels = labels,
    is_selected   = is_selected_vec
  )
}
