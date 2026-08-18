# =============================================================================
# fs_focus_tag(): one source of truth for sg_focus stem tags
# =============================================================================
# Result-file stems encode the focus that produced them.  Every simulation
# driver used to carry its own pasted switch() for that, and they drifted: the
# maxeffCons drivers fell through to a default arm and stem-tagged their output
# "maxeffCons" on the DINA path, where the run had actually ranked by the plain
# effect argmax.  Files named for a rule that did not run are not recoverable
# after the fact, so the mapping lives here and the drivers call it.
# =============================================================================

#' Stem tag for a (subgroup_method, sg_focus) pair
#'
#' Returns the short label used in result-file stems and plot annotations for
#' the rule that a given \code{(subgroup_method, sg_focus)} combination
#' actually runs -- not the spelling the caller passed.  The two differ
#' whenever an alias or a behavioural synonym is in play, which is precisely
#' when a hand-written tag goes wrong.
#'
#' The collapse is engine-specific:
#'
#' \describe{
#'   \item{\code{"consistency"}}{\code{eff}, \code{hr} and \code{maxcons} all
#'     tag as \code{"maxcons"} -- they name the consistency argmax.
#'     \code{maxeff} and \code{maxeffCons} stay distinct, because the
#'     consistency floor genuinely separates them here.}
#'   \item{\code{"dina"}, \code{"grf"}}{\code{eff}, \code{hr}, \code{maxcons},
#'     \code{maxeff} and \code{maxeffCons} all tag as \code{"eff"}.  Neither
#'     engine computes a Pcons, so the consistency qualifier has nothing to
#'     bind to and all five rank by \code{order(-eff)}.}
#' }
#'
#' The band foci collapse the same way on every engine: \code{hrMaxSG} tags as
#' \code{"effMaxSG"} and \code{hrMinSG} as \code{"effMinSG"}.  \code{maxSG} and
#' \code{minSG} are canonical foci in their own right and pass through.
#'
#' @param subgroup_method Character scalar. One of \code{"consistency"},
#'   \code{"dina"}, \code{"grf"}.  Any value other than \code{"consistency"}
#'   takes the effect-argmax collapse.
#' @param sg_focus Character scalar. Any accepted \code{sg_focus} spelling,
#'   alias or canonical.
#'
#' @return Character scalar; the stem tag.
#'
#' @examples
#' fs_focus_tag("consistency", "eff")   # "maxcons"
#' fs_focus_tag("dina",        "eff")   # "eff"
#' fs_focus_tag("dina",        "maxeffCons")  # "eff" -- not "maxeffCons"
#' fs_focus_tag("grf",         "hrMaxSG")     # "effMaxSG"
#'
#' @seealso \code{\link{forestsearch}} for the per-engine semantics of each
#'   focus.
#' @export
fs_focus_tag <- function(subgroup_method, sg_focus) {
  if (!is.character(subgroup_method) || length(subgroup_method) != 1L ||
      is.na(subgroup_method)) {
    stop("`subgroup_method` must be a single non-NA character string.",
         call. = FALSE)
  }
  if (!is.character(sg_focus) || length(sg_focus) != 1L || is.na(sg_focus)) {
    stop("`sg_focus` must be a single non-NA character string.",
         call. = FALSE)
  }

  if (identical(subgroup_method, "consistency")) {
    switch(sg_focus,
           eff      = ,
           hr       = ,
           maxcons  = "maxcons",
           effMaxSG = ,
           hrMaxSG  = "effMaxSG",
           effMinSG = ,
           hrMinSG  = "effMinSG",
           sg_focus)
  } else {
    switch(sg_focus,
           eff        = ,
           hr         = ,
           maxcons    = ,
           maxeff     = ,
           maxeffCons = "eff",
           effMaxSG   = ,
           hrMaxSG    = "effMaxSG",
           effMinSG   = ,
           hrMinSG    = "effMinSG",
           sg_focus)
  }
}
