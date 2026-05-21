# R/grf_args.R
# Argument-list builders for the GRF subgroup-identification calls invoked
# inside forestsearch().  These exist for two reasons:
#
#   1. As the single source of truth for how forestsearch() invokes
#      grf.subg.harm.glm() / grf.subg.harm.survival().  The call sites in
#      forestsearch_main.R should use these builders rather than inlining
#      the argument lists, so the contract cannot drift between locations.
#
#   2. To support a public extractor `extract_grf_args(fs)` that returns
#      the exact list forestsearch() used at GRF time, enabling
#      reproducibility audits, sensitivity analyses, and independent
#      standalone GRF runs.
#
# All three functions return a list suitable for
#   do.call(grf.subg.harm.glm,      args)        # glm path
#   do.call(grf.subg.harm.survival, args)        # survival path


# ---- Internal: GLM-path builder -------------------------------------------
# Mirrors the construction at forestsearch_main.R (pre-refactor lines
# 1322–1338) one-for-one.  Used both by the internal call site and by
# extract_grf_args() for the binary/continuous/count outcome paths.

.build_grf_glm_args <- function(data,
                                confounders.name,
                                outcome.name,
                                treat.name,
                                id.name,
                                outcome_type,
                                n.min,
                                dmin.grf,
                                is.RCT,
                                grf_depth,
                                seedit,
                                return_selected_cuts_only,
                                adverse_outcome,
                                offset.name         = NULL,
                                overdispersion      = "none",
                                grf_count_transform = "log") {

  # Map forestsearch outcome_type to the grf.subg.harm.glm vocabulary
  ot_grf <- if (outcome_type == "count")       "count"
            else if (outcome_type == "binary") "binary"
            else                               "continuous"

  out <- list(
    data                      = data,
    confounders.name          = confounders.name,
    outcome.name              = outcome.name,
    treat.name                = treat.name,
    id.name                   = id.name,
    outcome_type              = ot_grf,
    n.min                     = n.min,
    dmin.grf                  = dmin.grf,
    RCT                       = is.RCT,
    details                   = FALSE,
    maxdepth                  = grf_depth,
    seedit                    = seedit,
    return_selected_cuts_only = return_selected_cuts_only,
    adverse_outcome           = adverse_outcome
  )

  # Count-specific parameters
  if (ot_grf == "count") {
    out$offset.name         <- offset.name
    out$overdispersion      <- overdispersion
    out$grf_count_transform <- grf_count_transform
  }
  out
}


# ---- Internal: Survival-path builder --------------------------------------
# Mirrors the construction at forestsearch_main.R (pre-refactor lines
# 1303–1318) one-for-one.

.build_grf_survival_args <- function(data,
                                     confounders.name,
                                     outcome.name,
                                     event.name,
                                     id.name,
                                     treat.name,
                                     frac.tau,
                                     n.min,
                                     dmin.grf,
                                     is.RCT,
                                     grf_depth,
                                     seedit,
                                     return_selected_cuts_only) {
  list(
    data                      = data,
    confounders.name          = confounders.name,
    outcome.name              = outcome.name,
    event.name                = event.name,
    id.name                   = id.name,
    treat.name                = treat.name,
    frac.tau                  = frac.tau,
    n.min                     = n.min,
    dmin.grf                  = dmin.grf,
    RCT                       = is.RCT,
    details                   = FALSE,
    maxdepth                  = grf_depth,
    seedit                    = seedit,
    return_selected_cuts_only = return_selected_cuts_only
  )
}


# ---- Public: extract_grf_args ---------------------------------------------

#' Recover the GRF arguments used inside a forestsearch run
#'
#' Returns the exact list of arguments that \code{forestsearch()} passed
#' to its inner GRF subgroup-identification call
#' (\code{grf.subg.harm.survival()} for survival outcomes,
#' \code{grf.subg.harm.glm()} for binary, continuous, and count
#' outcomes).  The result is suitable for
#' \code{do.call(grf.subg.harm.glm, args)} or
#' \code{do.call(grf.subg.harm.survival, args)} to reproduce the
#' standalone GRF run with identical settings.
#'
#' Typical uses:
#' \itemize{
#'   \item Sensitivity analyses: modify one or two arguments and re-fit
#'         GRF without re-running the full forestsearch pipeline, e.g.
#'         \code{do.call(grf.subg.harm.glm,
#'           modifyList(extract_grf_args(fs), list(maxdepth = 3L)))}.
#'   \item Reproducibility audits: confirm what the inner GRF saw.
#'   \item Independent variable-importance or policy-tree extraction
#'         from a forest fit with the same configuration.
#' }
#'
#' This function reads from \code{fs$args_call_all} (populated by
#' \code{forestsearch()} via \code{mget()} of its formals) and from
#' \code{fs$outcome_type}.  It is the user-facing complement to the
#' internal builders that the \code{forestsearch()} call sites use.
#'
#' @param fs A fitted \code{forestsearch} object.
#'
#' @return A named list of arguments.  For survival outcomes the list
#'   targets \code{\link{grf.subg.harm.survival}}; for binary,
#'   continuous, or count outcomes it targets
#'   \code{\link{grf.subg.harm.glm}}.  The list carries an attribute
#'   \code{"grf_function"} holding the target function name as a
#'   character string.
#'
#' @examples
#' \dontrun{
#' # Reproduce the inner GRF run
#' args <- extract_grf_args(fs)
#' fn   <- match.fun(attr(args, "grf_function"))
#' grf_standalone <- do.call(fn, args)
#'
#' # Sensitivity: deeper trees, everything else unchanged
#' args_deep <- modifyList(args, list(maxdepth = 3L))
#' grf_deep  <- do.call(fn, args_deep)
#' }
#'
#' @seealso \code{\link{forestsearch}},
#'   \code{\link{grf.subg.harm.glm}},
#'   \code{\link{grf.subg.harm.survival}}.
#' @export
extract_grf_args <- function(fs) {
  if (is.null(fs$args_call_all)) {
    stop("extract_grf_args(): fs$args_call_all is missing. ",
         "Is `fs` a forestsearch() result?", call. = FALSE)
  }
  a  <- fs$args_call_all
  ot <- fs$outcome_type %||% a$outcome_type %||% "continuous"

  if (ot == "survival") {
    args <- .build_grf_survival_args(
      data                      = a$df.analysis,
      confounders.name          = a$confounders.name,
      outcome.name              = a$outcome.name,
      event.name                = a$event.name,
      id.name                   = a$id.name,
      treat.name                = a$treat.name,
      frac.tau                  = a$frac.tau,
      n.min                     = a$n.min,
      dmin.grf                  = a$dmin.grf,
      is.RCT                    = a$is.RCT,
      grf_depth                 = a$grf_depth,
      seedit                    = a$seedit,
      return_selected_cuts_only = a$return_selected_cuts_only
    )
    attr(args, "grf_function") <- "grf.subg.harm.survival"
  } else {
    args <- .build_grf_glm_args(
      data                      = a$df.analysis,
      confounders.name          = a$confounders.name,
      outcome.name              = a$outcome.name,
      treat.name                = a$treat.name,
      id.name                   = a$id.name,
      outcome_type              = ot,
      n.min                     = a$n.min,
      dmin.grf                  = a$dmin.grf,
      is.RCT                    = a$is.RCT,
      grf_depth                 = a$grf_depth,
      seedit                    = a$seedit,
      return_selected_cuts_only = a$return_selected_cuts_only,
      adverse_outcome           = a$adverse_outcome,
      offset.name               = a$offset.name,
      overdispersion            = a$overdispersion         %||% "none",
      grf_count_transform       = a$grf_count_transform    %||% "log"
    )
    attr(args, "grf_function") <- "grf.subg.harm.glm"
  }
  args
}


# Local null-coalesce helper (no rlang dependency)
`%||%` <- function(a, b) if (is.null(a)) b else a
