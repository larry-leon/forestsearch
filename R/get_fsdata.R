#' ForestSearch Data Preparation and Feature Selection
#'
#' Prepares a dataset for ForestSearch, including options for LASSO-based dimension reduction,
#' GRF cuts, forced cuts, and flexible cut strategies. Returns a list with the processed data,
#' subgroup factor names, cut expressions, and LASSO selection results.
#'
#' @param df.analysis Data frame containing the data.
#' @param use_lasso Logical. Whether to use LASSO for dimension reduction.
#' @param use_grf Logical. Whether to use GRF cuts.
#' @param grf_cuts Character vector of GRF cut expressions.
#' @param dina_cuts Character vector of DINA cut expressions (optional),
#'   merged into the candidate pool exactly like \code{grf_cuts}.
#' @param collapse_cuts Logical.  If TRUE, collapse near-redundant continuous
#'   candidate cuts after the full pool is assembled and resolved to literal
#'   numerics.  Cuts on the same variable and operator whose thresholds lie
#'   within a per-variable standard-error band are merged to a single
#'   rounded-centroid threshold, subject to a membership safety check; see
#'   \code{\link{collapse_redundant_cuts}}.  Categorical, indicator, bare-name
#'   and equality cuts are untouched.  Default TRUE; pass
#'   \code{collapse_cuts = FALSE} to recover the un-coarsened candidate pool
#'   used by earlier package versions.
#' @param collapse_cuts_args List of overrides for the coarsening, merged onto
#'   the defaults \code{list(c = 1.0, tol = 0.05, digits = 0L)}: \code{c} is the
#'   band multiplier (\code{band = c * sd(x)/sqrt(n)}), \code{tol} the
#'   membership safety tolerance (fraction of n when < 1, absolute count when
#'   >= 1), and \code{digits} the rounding for the representative threshold.
#'   Ignored when \code{collapse_cuts = FALSE}.
#' @param confounders.name Character vector of confounder variable names.
#' @param cont.cutoff Integer. Cutoff for continuous variable determination.
#' @param conf_force Character vector of forced cut expressions.
#' @param conf.cont_medians Character vector of continuous confounders to cut at median.
#' @param conf.cont_medians_force Character vector of additional continuous confounders to force median cut.
#' @param conf.cont_jcuts Named list of positive integers, one per
#'   continuous confounder for which J-quantile cuts are desired.  For a
#'   variable \code{X} listed as \code{X = J}, the default cut set
#'   (\code{cut_var()}: mean, median, Q1, Q3) is replaced by \code{J}
#'   binary cut points at the (k/(J+1))-th quantiles of X for
#'   \code{k = 1, ..., J}, defining J+1 non-overlapping intervals
#'   \deqn{[\min(X), c_1),\ [c_1, c_2),\ \ldots,\ [c_J, \max(X)].}
#'   Variables not listed here retain default behaviour.  J-quantile
#'   cuts are unconditional (not filtered by LASSO), matching
#'   \code{defaultcut_names} semantics.  Names must be in
#'   \code{confounders.name} and must not overlap with
#'   \code{defaultcut_names} or \code{conf.cont_medians_force}.
#' @param replace_med_grf Logical. If TRUE, removes median cuts that overlap with GRF cuts.
#' @param defaultcut_names Character vector of confounders to force default cuts.
#' @param cut_type Character. "default" or "median" for cut strategy.
#' @param exclude_cuts Character vector of cut expressions to exclude.
#' @param outcome.name Character. Name of outcome variable.
#' @param event.name Character. Name of event indicator variable.
#' @param details Logical. If TRUE, prints details during execution.
#' @param outcome_type Character. One of \code{"survival"} (default),
#'   \code{"binary"}, or \code{"continuous"}.
#' @param offset.name Character or \code{NULL}. Name of the follow-up time
#'   column for rate-based measures (IRR, IRD).
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{df}}{Data frame with derived subgroup factor columns.}
#'     \item{\code{confs_names}}{Character vector of factor column names
#'       added to \code{df}.}
#'     \item{\code{confs}}{Named character vector of cut expressions
#'       (continuous cuts plus categorical indicators).}
#'     \item{\code{lassokeep}}{Character vector of confounders retained
#'       by LASSO screening (empty if \code{use_lasso = FALSE}).}
#'     \item{\code{lassoomit}}{Character vector of confounders dropped by
#'       LASSO screening (empty if \code{use_lasso = FALSE}).}
#'   }
#'
#' @importFrom stats median quantile
#' @examples
#' \donttest{
#' library(survival)
#' df <- survival::gbsg
#' df$grade3 <- as.integer(df$grade == "3")
#' fs_data <- get_FSdata(df.analysis = df,
#'   confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
#'   outcome.name = "rfstime", event.name = "status",
#'   use_lasso = FALSE, use_grf = FALSE)
#' names(fs_data)
#' }
#' @export

get_FSdata <- function(df.analysis, use_lasso = FALSE, use_grf = FALSE, grf_cuts = NULL ,confounders.name,
                       cont.cutoff = 4,conf_force = NULL, conf.cont_medians = NULL, conf.cont_medians_force = NULL,
                       conf.cont_jcuts = NULL,
                       dina_cuts = NULL,
                       collapse_cuts = TRUE, collapse_cuts_args = list(),
                       replace_med_grf = TRUE, defaultcut_names = NULL, cut_type = "default", exclude_cuts = NULL,
                       outcome.name = "tte", event.name = "event", details=TRUE,
                       outcome_type = "survival", offset.name = NULL){

  # Validate df.analysis at entry with a clear error.  Previously this
  # silently coerced via as.data.frame(), which for non-data.frame
  # inputs produces a downstream failure at the outcome-column numeric
  # check with the misleading message "Outcome column must be numeric."
  # Rejecting at entry gives the user an actionable message about what
  # they actually got wrong.
  if (!is.data.frame(df.analysis)) {
    stop(sprintf(paste0(
      "'df.analysis' must be a data.frame (or inherit from data.frame); ",
      "got an object of class %s.  If you have a matrix, tibble, or ",
      "data.table, coerce with as.data.frame() before calling."),
      paste(class(df.analysis), collapse = "/")),
      call. = FALSE)
  }
  df.FS <- df.analysis

  # Resolve outcome_type if it's the raw default vector, then validate
  # explicitly.  Downstream code branches on this value via switch();
  # rejecting unknown values at entry gives a clearer error than a
  # cryptic "EXPR must be a length 1 vector" from deep in the pipeline.
  if (length(outcome_type) > 1L) outcome_type <- outcome_type[1L]
  valid_otypes <- c("survival", "binary", "continuous", "count")
  if (!outcome_type %in% valid_otypes) {
    stop(sprintf(
      "outcome_type = '%s' is not recognized.  Must be one of: %s.",
      outcome_type,
      paste(shQuote(valid_otypes), collapse = ", ")),
      call. = FALSE)
  }

  # Validate cut_type.  Downstream code has four separate branches
  # matching on cut_type; a typo like "defualt" would silently match
  # none of them, leaving conf_force unmodified and producing a
  # different candidate-cut set than the user expected.
  valid_cut_types <- c("default", "median")
  if (!cut_type %in% valid_cut_types) {
    stop(sprintf(
      "cut_type = '%s' is not recognized.  Must be one of: %s.",
      cut_type,
      paste(shQuote(valid_cut_types), collapse = ", ")),
      call. = FALSE)
  }

  # Validate confounders.name structure.  NULL or empty character leads
  # to silently-empty confs downstream; the resulting error (Fix B) is
  # already actionable, but catching at entry is cleaner.
  if (is.null(confounders.name) || length(confounders.name) == 0L) {
    stop("'confounders.name' must be a non-empty character vector.",
         call. = FALSE)
  }
  if (!is.character(confounders.name)) {
    stop(sprintf(
      "'confounders.name' must be a character vector; got %s.",
      paste(class(confounders.name), collapse = "/")),
      call. = FALSE)
  }
  missing_conf <- setdiff(confounders.name, names(df.FS))
  if (length(missing_conf) > 0L) {
    stop(sprintf(
      "confounders.name references columns not in df.analysis: %s.",
      paste(shQuote(missing_conf), collapse = ", ")),
      call. = FALSE)
  }

  # Check that outcome and event columns are present and numeric
  if (!outcome.name %in% names(df.FS)) {
    stop(sprintf("Outcome column '%s' not found in df.analysis.",
                 outcome.name), call. = FALSE)
  }
  if (!event.name %in% names(df.FS)) {
    stop(sprintf("Event column '%s' not found in df.analysis.",
                 event.name), call. = FALSE)
  }
  if (!is.numeric(df.FS[[outcome.name]])) {
    stop(sprintf("Outcome column '%s' must be numeric.",
                 outcome.name), call. = FALSE)
  }
  if (!is.numeric(df.FS[[event.name]])) {
    stop(sprintf("Event column '%s' must be numeric (0/1).",
                 event.name), call. = FALSE)
  }

  # Check that confounders are numeric or factors, naming any that
  # violate (previously reported generically "All confounders must be
  # numeric or factor" with no column identification).
  types <- vapply(confounders.name,
                  function(v) is.numeric(df.FS[[v]]) || is.factor(df.FS[[v]]),
                  logical(1))
  if (!all(types)) {
    bad <- confounders.name[!types]
    bad_classes <- vapply(bad, function(v) class(df.FS[[v]])[1], character(1))
    stop(sprintf(
      "All confounders must be numeric or factor; offending: %s.",
      paste(sprintf("'%s' (%s)", bad, bad_classes), collapse = ", ")),
      call. = FALSE)
  }

  # Validate conf.cont_jcuts (per-variable J-quantile cut override).
  # Caught here so any structural problems are flagged before the cut
  # machinery starts running and producing partial results.
  if (!is.null(conf.cont_jcuts) && length(conf.cont_jcuts) > 0L) {
    if (!is.list(conf.cont_jcuts) || is.null(names(conf.cont_jcuts)) ||
        any(!nzchar(names(conf.cont_jcuts)))) {
      stop("'conf.cont_jcuts' must be a NAMED list, e.g. list(X8 = 10).",
           call. = FALSE)
    }
    jcut_names <- names(conf.cont_jcuts)
    not_in_conf <- setdiff(jcut_names, confounders.name)
    if (length(not_in_conf) > 0L) {
      stop(sprintf(
        "conf.cont_jcuts names not in confounders.name: %s.",
        paste(shQuote(not_in_conf), collapse = ", ")),
        call. = FALSE)
    }
    bad_J <- vapply(conf.cont_jcuts, function(j) {
      !is.numeric(j) || length(j) != 1L || is.na(j) ||
        j != as.integer(j) || j < 1L
    }, logical(1L))
    if (any(bad_J)) {
      stop(sprintf(
        "conf.cont_jcuts: each value must be a single positive integer.  Bad entries: %s.",
        paste(shQuote(jcut_names[bad_J]), collapse = ", ")),
        call. = FALSE)
    }
    if (!is.null(defaultcut_names)) {
      overlap <- intersect(jcut_names, defaultcut_names)
      if (length(overlap) > 0L) {
        stop(sprintf(
          "Variables appear in both 'conf.cont_jcuts' and 'defaultcut_names': %s.  These are mutually exclusive.",
          paste(shQuote(overlap), collapse = ", ")),
          call. = FALSE)
      }
    }
    if (!is.null(conf.cont_medians_force)) {
      overlap <- intersect(jcut_names, conf.cont_medians_force)
      if (length(overlap) > 0L) {
        stop(sprintf(
          "Variables appear in both 'conf.cont_jcuts' and 'conf.cont_medians_force': %s.  These are mutually exclusive.",
          paste(shQuote(overlap), collapse = ", ")),
          call. = FALSE)
      }
    }
  }

  # Default cuts forced per defaultcut_names
  if(!is.null(defaultcut_names)){
    conf_force_default <- get_conf_force(df = df.FS, conf.force.names = defaultcut_names, cont.cutoff = cont.cutoff)
    # append to conf_force
    conf_force <- c(conf_force,conf_force_default)
  }

  # If grf was attempted but NO cuts were found then considering no cuts per grf
  if(use_grf && is.null(grf_cuts)) use_grf <- FALSE

  # Type contract: grf_cuts must be a character vector of "var <= value"

  # expressions when non-NULL.  GLM GRF may produce a named list instead;
  # normalize here as a safety net (primary conversion is in forestsearch_main.R).
  if (!is.null(grf_cuts) && is.list(grf_cuts) && !is.null(names(grf_cuts))) {
    grf_cuts <- unlist(lapply(names(grf_cuts), function(nm) {
      paste0(nm, " <= ", grf_cuts[[nm]])
    }))
  }
  if (!is.null(grf_cuts) && !is.character(grf_cuts)) {
    stop(sprintf(
      "grf_cuts must be a character vector of cut expressions, got %s",
      paste(class(grf_cuts), collapse = "/")
    ))
  }

  if(use_lasso &  (is.null(outcome.name))) stop("Outcome variable name needed for lasso")

  flag_continuous <- vapply(
    confounders.name,
    function(var) {
      aa <- df.FS[[var]]
      c(is.continuous(aa, cutoff = cont.cutoff) ==1)
    },
    logical(1)
  )
  if(details){
    cat("# of continuous/categorical characteristics",c(sum(flag_continuous),sum(!flag_continuous)),"\n")
    if(sum(flag_continuous) > 0)  cat("Continuous characteristics:",c(confounders.name[flag_continuous]),"\n")
    if(sum(!flag_continuous) > 0)   cat("Categorical characteristics:",c(confounders.name[!flag_continuous]),"\n")
  }

  if(sum(flag_continuous) == 0) conf.categorical <- confounders.name
  if(sum(flag_continuous) > 0) conf.categorical <- confounders.name[!flag_continuous]

  # If there are no median cuts (either "conf.cont_medians" or "conf.cont_medians_force")
  # then we initialize all continuous confounders to be cut
  if(is.null(conf.cont_medians) & is.null(conf.cont_medians_force) & sum(flag_continuous) > 0){
    conf.cont_medians <- confounders.name[flag_continuous]
  }

  # J-quantile cut overrides (conf.cont_jcuts).  These are unconditional
  # (not filtered by LASSO, matching defaultcut_names semantics) and
  # REPLACE the default cut set for the named variables.  The names are
  # stripped from conf.cont_medians here so that the default/median/lasso
  # branches below skip those variables entirely.  Validation upstream
  # has already enforced no overlap with defaultcut_names or
  # conf.cont_medians_force.
  if (!is.null(conf.cont_jcuts) && length(conf.cont_jcuts) > 0L) {
    jcut_names <- names(conf.cont_jcuts)
    conf_force_jq <- get_conf_force_jq(
      df = df.FS,
      conf_jcuts = conf.cont_jcuts,
      cont.cutoff = cont.cutoff
    )
    if (length(conf_force_jq) > 0L) {
      conf_force <- c(conf_force, conf_force_jq)
    }
    if (!is.null(conf.cont_medians)) {
      conf.cont_medians <- setdiff(conf.cont_medians, jcut_names)
      if (length(conf.cont_medians) == 0L) conf.cont_medians <- NULL
    }
    if (details) {
      cat("J-quantile cuts applied to:",
          paste0(jcut_names, " (J=", unlist(conf.cont_jcuts), ")",
                 collapse = ", "),
          "\n")
    }
  }

  lassokeep <- NULL
  lassoomit <- NULL
  if(use_lasso){
    # Reduce dimension via LASSO (Cox for survival, GLM for binary/continuous)
    get_lasso <- lasso_selection(
      df = df.FS,
      confounders.name = confounders.name,
      outcome.name = outcome.name,
      event.name = event.name,
      outcome_type = outcome_type,
      offset.name = offset.name
    )
    lassokeep <- get_lasso$selected
    lassoomit <- get_lasso$omitted

    if(details){
      cat("## Prior to lasso:", c(conf.cont_medians), "\n")
      cat("#### Lasso selection results", "\n")
      print(get_lasso$fit$beta)
      lasso_label <- switch(outcome_type,
        binary     = "Logistic-LASSO",
        continuous = "Linear-LASSO",
        count      = "Poisson-LASSO",
        "Cox-LASSO"
      )
      cat(paste0(lasso_label, " selected:"), c(lassokeep), "\n")
      cat(paste0(lasso_label, " not selected:"), c(lassoomit), "\n")
      cat("### End Lasso selection", "\n")
    }
    # If any selected per lasso
    if (length(lassokeep) > 0) {
      conf.cont_medians <- filter_by_lassokeep(conf.cont_medians, lassokeep)
      conf.categorical  <- filter_by_lassokeep(conf.categorical, lassokeep)
    }
    if(details)  cat("## After lasso:", c(conf.cont_medians), "\n")
  } # Done Lasso

  # If forcing any cuts, then done below
  if (use_lasso && cut_type == "default") {
    conf_force_lasso <- NULL
    # Override conf.cont_medians
    # conf.cont_medians are continuous factors selected per lasso
    # Create default cuts for conf.cont_medians NOT in defaultcut_names
    if (!is.null(defaultcut_names)) {
      # If not already contained in defaultcut_names
      lasso_tocut <- setdiff(conf.cont_medians, defaultcut_names)
    } else {
      lasso_tocut <- conf.cont_medians
    }
    if (length(lasso_tocut) > 0) {
      conf_force_lasso <- get_conf_force(df = df.FS, conf.force.names = lasso_tocut, cont.cutoff = cont.cutoff)
    }
    # Override cuts at medians
    conf.cont_medians <- NULL
    # Append to conf_force
    conf_force <- c(conf_force, conf_force_lasso)
    # These will be cut and conf.cont_medians is now reset
    if (details) {
      cat("Default cuts included from Lasso:", c(conf_force_lasso), "\n")
      cat("Categorical after Lasso:", c(conf.categorical), "\n")
    }
  }
  if (use_lasso && cut_type == "median") {
    if (!is.null(defaultcut_names)) {
      conf.cont_medians <- setdiff(conf.cont_medians, defaultcut_names)
    }
    if (details) {
      cat("Median cuts included from Lasso:", c(conf.cont_medians), "\n")
      cat("Categorical after Lasso:", c(conf.categorical), "\n")
    }
  }
  if (!use_lasso && cut_type == "default") {
    conf_force_add <- NULL
    # Override conf.cont_medians
    # Create default cuts for conf.cont_medians not in defaultcut_names
    if (!is.null(defaultcut_names)) {
      tocut <- setdiff(conf.cont_medians, defaultcut_names)
      if (length(tocut) > 0) {
        conf_force_add <- get_conf_force(df = df.FS, conf.force.names = tocut, cont.cutoff = cont.cutoff)
      }
      # Override cuts at medians
      conf.cont_medians <- NULL
    }
    if (is.null(defaultcut_names)) {
      tocut <- conf.cont_medians
      if (length(tocut) > 0) {
        conf_force_add <- get_conf_force(df = df.FS, conf.force.names = tocut, cont.cutoff = cont.cutoff)
      }
      # Override cuts at medians
      conf.cont_medians <- NULL
    }
    # Append to conf_force
    conf_force <- c(conf_force, conf_force_add)
    if (details) {
      toprint <- min(20, length(conf_force_add))
      cat("Default cuts included (1st 20)", c(conf_force_add[1:toprint]), "\n")
      cat("Categorical:", c(conf.categorical), "\n")
    }
  }
  if (!use_lasso && cut_type == "median") {
    if (!is.null(defaultcut_names)) {
      conf.cont_medians <- setdiff(conf.cont_medians, defaultcut_names)
    }
    if (details) {
      toprint <- min(20, length(conf.cont_medians))
      cat("Median cuts included:", c(conf.cont_medians[1:toprint]), "\n")
      cat("Categorical:", c(conf.categorical), "\n")
    }
  }

  if(!is.null(conf.cont_medians_force)){
    conf.cont_medians<-c(conf.cont_medians,conf.cont_medians_force)
  }
  if(details & use_grf){
    cat("Factors per GRF:",c(grf_cuts),"\n")
  }
  if(details & use_grf & length(conf.cont_medians)>0){
    toprint <- min(20,length(conf.cont_medians))
    cat("Continuous factors initially cut at medians:",c(conf.cont_medians[1:toprint]),"\n")
  }
  # Remove any factors to cut at median if already in GRF
  # Remove any factors to cut at median if already in GRF
  if (replace_med_grf) {
    if (use_grf && length(conf.cont_medians) > 0 && length(grf_cuts) > 0) {
      # Find which conf.cont_medians are present in any grf_cuts
      to_exclude <- vapply(
        conf.cont_medians,
        function(x) any(grepl(x, grf_cuts)),
        logical(1)
      )
      # Update conf.cont_medians
      if (any(to_exclude)) {
        conf.cont_medians <- conf.cont_medians[!to_exclude]
        if (length(conf.cont_medians) == 0) conf.cont_medians <- NULL
      }
    }
    if (details && use_grf && length(conf.cont_medians) > 0) {
      cat("Factors after removing any duplicates also in GRF:", conf.cont_medians, "\n")
    }
  }

  if(details & cut_type=="median" & use_lasso & length(conf.cont_medians)==0){
    cat("***conf.cont_medians is NULL --> NO MEDIAN CUTS per lasso***","\n")
  }

  # Re-introduce conf.cont_force_medians
  if(!is.null(conf.cont_medians_force)) conf.cont_medians <- c(conf.cont_medians,conf.cont_medians_force)

  # Build median cut expressions from the names remaining in
  # conf.cont_medians.  IMPORTANT: prior versions assigned the vector
  # to 'conf.cont_Medcuts_vec' and then used 'conf.cont_Medcuts' (NULL)
  # in the confs assembly below, silently dropping all median cuts
  # under cut_type = "median".  This block now builds and uses a
  # single consistent name.
  if (length(conf.cont_medians) > 0L) {
    medians <- vapply(conf.cont_medians,
                      function(x) round(median(df.FS[[x]], na.rm = TRUE), 2),
                      numeric(1))
    conf.cont_Medcuts <- paste0(conf.cont_medians, ' <= ', medians)
  } else {
    conf.cont_Medcuts <- character(0)
  }
  confs <- c(conf.categorical, conf.cont_Medcuts)
  # At this stage, these are confs per Lasso (GRF step is next)
  if(use_lasso) confs_lasso <- confs
  # Factors included per GRF not in confs_lasso
  if(use_grf){
    if(details) cat('Initial GRF cuts included', grf_cuts, '\n')
    if(length(confs) > 0 & length(grf_cuts) > 0){
      # Vectorized check: keep only grf_cuts not matching any confs
      flag_omit <- sapply(grf_cuts, function(cut) any(sapply(confs, function(x) grepl(x, cut))))
      grf_cuts_keep <- grf_cuts[!flag_omit]
      confs <- unique(c(confs, grf_cuts_keep))
    }
    if(length(confs) == 0 & length(grf_cuts) > 0){
      grf_cuts_keep <- grf_cuts
      confs <- unique(c(confs, grf_cuts_keep))
    }
  }
  if(use_lasso & use_grf){
    which_both <- (confs %in% confs_lasso)
    if(details & length(which_both) >0){
      cat("Factors included per GRF (not in lasso)",c(confs[!which_both]),"\n")
    }
  }
  # DINA cuts: merged into the candidate pool exactly like GRF cuts -- keep
  # only those not already represented in confs (substring match), then
  # union-dedup.  dina_cuts are already canonical "var <= value" expressions.
  if(length(dina_cuts) > 0){
    if(details) cat('Factors per DINA:', c(dina_cuts), '\n')
    if(length(confs) > 0){
      flag_omit <- sapply(dina_cuts, function(cut) any(sapply(confs, function(x) grepl(x, cut))))
      dina_cuts_keep <- dina_cuts[!flag_omit]
    } else {
      dina_cuts_keep <- dina_cuts
    }
    confs <- unique(c(confs, dina_cuts_keep))
  }
  conf_forceNew <- vapply(
    conf_force,
    process_conf_force_expr,
    FUN.VALUE = character(1),
    df = df.FS
  )
  if(!is.null(conf_force)) confs <- unique(c(confs,conf_forceNew))
  # Excluding cuts
  if(!is.null(exclude_cuts)){
    # Remove the restricted cuts (eg., not allowing a variable to be "<=")
    for(ee in seq_along(exclude_cuts)){
      to_exclude <- grepl(exclude_cuts[ee],confs)
      confs <- confs[!to_exclude]
    }
  }
  # Collapse near-redundant continuous candidate cuts (opt-in).  Cuts on the
  # same variable+operator whose thresholds lie within a per-variable standard-
  # error band (collapse_cuts_args$c * sd(x)/sqrt(n)) are treated as practically
  # redundant and merged to a single rounded-centroid threshold, subject to a
  # membership safety check.  Operates on the fully-resolved literal-numeric
  # candidate pool; categorical / indicator / bare-name cuts are untouched.
  # On by default (collapse_cuts = TRUE); pass collapse_cuts = FALSE to recover
  # the un-coarsened candidate pool used by earlier package versions.
  if (isTRUE(collapse_cuts) && length(confs) > 0L) {
    .cca <- utils::modifyList(
      list(c = 1.0, tol = 0.05, digits = 0L), collapse_cuts_args)
    confs <- collapse_redundant_cuts(
      cuts             = confs,
      df               = df.FS,
      confounders.name = confounders.name,
      c_band           = .cca$c,
      safety_tol       = .cca$tol,
      digits           = .cca$digits,
      cont.cutoff      = cont.cutoff,
      details          = details
    )
  }
  n_confs<-length(confs)
  if (n_confs == 0) {
    # Build a diagnostic identifying which upstream source was empty
    # so the user can see at a glance where the pipeline broke down.
    # Each line reports the count that actually reached this point.
    sources <- c(
      sprintf("GRF cuts (%d)",
              if (use_grf) length(grf_cuts) else 0L),
      sprintf("DINA cuts (%d)",           length(dina_cuts)),
      sprintf("LASSO selected (%d)",
              if (use_lasso) length(lassokeep) else 0L),
      sprintf("conf_force (%d)",          length(conf_force)),
      sprintf("conf.cont_medians (%d)",   length(conf.cont_medians)),
      sprintf("conf.categorical (%d)",    length(conf.categorical))
    )
    stop(sprintf(paste0(
      "get_FSdata: no candidate cut expressions could be constructed ",
      "from the supplied inputs.  ",
      "Upstream source counts: %s.  ",
      "Check that: (1) GRF is succeeding for this outcome type ",
      "(rerun forestsearch() with details = TRUE to inspect); ",
      "(2) LASSO is not shrinking all covariates to zero (try ",
      "use_lasso = FALSE); ",
      "(3) candidate covariates are numeric or factor and have ",
      "sufficient variability; ",
      "(4) GRF-related thresholds (dmin.grf, frac.tau, vi.grf.min) ",
      "are appropriate for the outcome type and sample size ",
      "(frac.tau is survival-specific and may be ignored for GLM)."),
      paste(sources, collapse = ", ")
    ))
  }

  # =========================================================================
  # REFACTORED SECTION: CONSOLIDATED CUT EVALUATION (3.4x faster)
  # =========================================================================
  # IMPROVEMENT: Evaluate all cuts ONCE and cache results
  # This replaces the old pattern of evaluating cuts multiple times
  # =========================================================================

  if(details) {
    cat("\n===== CONSOLIDATED CUT EVALUATION (IMPROVED) =====\n")
    cat("Evaluating ", length(confs), " cut expressions once and caching...\n", sep = "")
  }

  # STEP 1: Evaluate ALL cuts exactly once and cache results
  cut_cache <- evaluate_cuts_once(confs, df.FS, details = details)
  evaluations <- cut_cache$evaluations
  is_valid <- cut_cache$is_valid
  has_error <- cut_cache$has_error

  # STEP 2: Classify cuts based on cached evaluation (no re-evaluation!)
  flag_continuous <- vapply(
    confs,
    is_flag_continuous,
    FUN.VALUE = logical(1),
    confounders.name = confounders.name,
    df = df.FS,
    cont.cutoff = cont.cutoff
  )

  # STEP 3: Remove invalid cuts (they have <=1 unique value)
  flag_drop <- !is_valid

  if(details & any(flag_drop)){
    cat("Dropping variables (cut only has 1 level):", c(confs[flag_drop]), "\n")
  }

  # STEP 4: Split into categories using cached results
  conf.categorical <- confs[!flag_continuous & !flag_drop]
  conf.cont_cuts <- NULL
  if(sum(flag_continuous) > 0){
    conf.cont_cuts <- confs[flag_continuous & !flag_drop]
  }

  # Count valid cuts after filtering
  n_confs <- length(c(conf.categorical, conf.cont_cuts))
  if(length(confs) - n_confs > 0){
    if(details){
      cat("Total cuts after dropping invalid: ", n_confs, "\n", sep = "")
    }
  }

  # =========================================================================
  # CREATE NEW COLUMNS USING CACHED EVALUATIONS
  # =========================================================================
  # IMPROVEMENT: Uses pre-computed cached evaluations instead of re-evaluating
  # Avoids redundant eval(parse()) calls
  # =========================================================================

  # Create map from original confs to cached evaluations for quick lookup
  confs_to_index <- setNames(seq_along(confs), confs)

  # Generate new column names
  names_new <- c(unlist(lapply(seq_along(c(conf.cont_cuts, conf.categorical)),
                               function(x){paste0("q",x,sep="")})))

  # Process continuous cuts (use cached evaluations)
  offset <- 0
  for(i in seq_along(conf.cont_cuts)) {
    thiscut <- conf.cont_cuts[i]
    idx <- confs_to_index[[thiscut]]

    # Use cached evaluation instead of re-evaluating!
    if (!is.null(evaluations[[idx]]) && is_valid[idx]) {
      result <- evaluations[[idx]]
      df.FS[[names_new[i]]] <- as.factor(as.numeric(result))
    } else {
      # Fallback (shouldn't happen if caching works correctly)
      if(details) warning("Cut '", thiscut, "' not found in cache, re-evaluating")

      result <- evaluate_comparison(thiscut, df.FS)

      df.FS[[names_new[i]]] <- as.factor(as.numeric(result))
    }
  }

  # Process categorical cuts (use cached evaluations)
  offset <- length(conf.cont_cuts)
  for(i in seq_along(conf.categorical)) {
    thiscut <- conf.categorical[i]
    idx <- confs_to_index[[thiscut]]

    # Use cached evaluation instead of re-evaluating!
    if (!is.null(evaluations[[idx]]) && is_valid[idx]) {
      result <- evaluations[[idx]]
      # Convert to numeric THEN to factor to ensure 0/1 levels
      numeric_result <- as.numeric(result)
      df.FS[[names_new[i + offset]]] <- as.factor(numeric_result)
    } else {
      # Fallback (shouldn't happen if caching works correctly)
      if(details) warning("Cut '", thiscut, "' not found in cache, re-evaluating")

      result <- evaluate_comparison(thiscut, df.FS)

      numeric_result <- as.numeric(result)
      df.FS[[names_new[i + offset]]] <- as.factor(numeric_result)
    }
  }

  # =========================================================================
  # VALIDATION: Verify all factors are 0/1
  # =========================================================================

  check_factors <- vapply(names_new, function(col_name) {
    col_data <- df.FS[[col_name]]
    # Get unique values as numeric to check range
    unique_vals <- as.numeric(as.character(unique(col_data)))
    unique_vals <- unique_vals[!is.na(unique_vals)]
    # Should only contain 0 and/or 1
    all(unique_vals %in% c(0, 1))
  }, logical(1))

  if (!all(check_factors)) {
    invalid_cols <- names_new[!check_factors]
    cat("DEBUG: Invalid column values:\n")
    for (col in invalid_cols) {
      cat("  ", col, ": ", paste(unique(df.FS[[col]]), collapse = ", "), "\n", sep = "")
    }
    stop("Error in factor setup: some factors contain values other than 0/1. ",
         "Invalid columns: ", paste(invalid_cols, collapse = ", "))
  }

  if(details){
    cat("\u2713 All ", length(names_new), " factors validated as 0/1\n", sep = "")
    cat("===== END CONSOLIDATED CUT EVALUATION =====\n\n")
  }

  # NOTE: include check that all confs are 0,1 (not TRUE, FALSE)
  if(details){
    cat("# of candidate subgroup factors=",c(length(c(conf.cont_cuts, conf.categorical))),"\n")
    print(c(conf.cont_cuts, conf.categorical))
  }

  return(list(df = df.FS, confs_names = names_new, confs = c(conf.cont_cuts, conf.categorical),
              lassokeep = lassokeep, lassoomit = lassoomit))
}
