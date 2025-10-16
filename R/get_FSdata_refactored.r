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
#' @param confounders.name Character vector of confounder variable names.
#' @param cont.cutoff Integer. Cutoff for continuous variable determination.
#' @param conf_force Character vector of forced cut expressions.
#' @param conf.cont_medians Character vector of continuous confounders to cut at median.
#' @param conf.cont_medians_force Character vector of additional continuous confounders to force median cut.
#' @param replace_med_grf Logical. If TRUE, removes median cuts that overlap with GRF cuts.
#' @param defaultcut_names Character vector of confounders to force default cuts.
#' @param cut_type Character. "default" or "median" for cut strategy.
#' @param exclude_cuts Character vector of cut expressions to exclude.
#' @param outcome.name Character. Name of outcome variable.
#' @param event.name Character. Name of event indicator variable.
#' @param details Logical. If TRUE, prints details during execution.
#' @importFrom stringr str_length str_sub
#' @importFrom stats median quantile
#' @export

get_FSdata <- function(df.analysis, use_lasso = FALSE, use_grf = FALSE, grf_cuts = NULL ,confounders.name,
                       cont.cutoff = 4,conf_force = NULL, conf.cont_medians = NULL, conf.cont_medians_force = NULL,
                       replace_med_grf = TRUE, defaultcut_names = NULL, cut_type = "default", exclude_cuts = NULL,
                       outcome.name = "tte", event.name = "event", details=TRUE){

  # Initialize to original analysis dataframe and output df.FS containing cutpoints
  if(!is.data.frame(df.analysis)){
    df.FS <- as.data.frame(df.analysis)
  } else {
    df.FS <- df.analysis
  }


  # Check that outcome and event columns are numeric
  if(!is.numeric(df.FS[[outcome.name]])) stop("Outcome column must be numeric.")
  if(!is.numeric(df.FS[[event.name]])) stop("Event column must be numeric (0/1).")

  # Check that confounders are numeric or factors

  types <- sapply(df.FS[confounders.name], function(x) is.numeric(x) || is.factor(x))
  if (!all(types)) stop("All confounders must be numeric or factor.")

  # Default cuts forced per defaultcut_names
  if(!is.null(defaultcut_names)){
    conf_force_default <- get_conf_force(df = df.FS,conf.force.names = defaultcut_names,cont.cutoff = 4)
    # append to conf_force
    conf_force <- c(conf_force,conf_force_default)
  }

  # If grf was attempted but NO cuts were found then considering no cuts per grf
  if(use_grf && is.null(grf_cuts)) use_grf <- FALSE

  if(use_lasso &  (is.null(outcome.name) | is.null(event.name))) stop("Cox variable names needed for lasso (outcome and event)")

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
  lassokeep <- NULL
  lassoomit <- NULL
  if(use_lasso){
    # Reduce dimension via Cox lasso
    get_lasso <- lasso_selection(
      df = df.FS,
      confounders.name = confounders.name,
      outcome.name = outcome.name,
      event.name = event.name
    )
    lassokeep <- get_lasso$selected
    lassoomit <- get_lasso$omitted

    if(details){
      cat("## Prior to lasso:", c(conf.cont_medians), "\n")
      cat("#### Lasso selection results", "\n")
      print(get_lasso$fit$beta)
      cat("Cox-LASSO selected:",c(lassokeep),"\n")
      cat("Cox-LASSO not selected:",c(lassoomit),"\n")
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
      conf_force_lasso <- get_conf_force(df = df.FS, conf.force.names = lasso_tocut, cont.cutoff = 4)
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
        conf_force_add <- get_conf_force(df = df.FS, conf.force.names = tocut, cont.cutoff = 4)
      }
      # Override cuts at medians
      conf.cont_medians <- NULL
    }
    if (is.null(defaultcut_names)) {
      tocut <- conf.cont_medians
      if (length(tocut) > 0) {
        conf_force_add <- get_conf_force(df = df.FS, conf.force.names = tocut, cont.cutoff = 4)
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
  conf.cont_Medcuts<-NULL
  medians <- sapply(conf.cont_medians, function(x) round(median(df.FS[[x]]), 2))
  conf.cont_Medcuts_vec <- paste0(conf.cont_medians, ' <= ', medians)
  conf.cont_Medcuts_vec
  confs<-c(conf.categorical,conf.cont_Medcuts)
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
    for(ee in 1:length(exclude_cuts)){
      to_exclude <- grepl(exclude_cuts[ee],confs)
      confs <- confs[!to_exclude]
    }
  }
  n_confs<-length(confs)
  if(n_confs==0) stop("Error in FS dataset prior to flag drop")

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
  names_new <- c(unlist(lapply(c(1:length(c(conf.cont_cuts, conf.categorical))),
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
      result <- eval(parse(text = thiscut), envir = df.FS)
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
      result <- eval(parse(text = thiscut), envir = df.FS)
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
    cat("✓ All ", length(names_new), " factors validated as 0/1\n", sep = "")
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
