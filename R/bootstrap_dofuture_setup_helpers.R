# Required packages for parallel bootstrap workers
BOOTSTRAP_REQUIRED_PACKAGES <- c(
  "data.table",    # Data manipulation in parallel workers
  "foreach",       # Foreach iteration framework
  "doFuture",      # Future-based parallel backend
  "survival"       # Cox model fitting for bias correction
)

#' Resolve parallel processing arguments for bootstrap
#'
#' If parallel_args not provided, falls back to forestsearch call's
#' parallel configuration. Always reports configuration to user.
#'
#' @param parallel_args List or empty list
#' @param forestsearch_call_args List from original forestsearch call
#' @return List with plan, workers, show_message
#' @keywords internal

resolve_bootstrap_parallel_args <- function(parallel_args, forestsearch_call_args) {
  # Use provided args if non-empty, otherwise inherit from forestsearch call
  if (length(parallel_args) == 0) {
    resolved_args <- as.list(forestsearch_call_args$parallel_args)
    message("Bootstrap parallel config: Using 'observed' data analysis forestsearch settings")
  } else {
    resolved_args <- parallel_args
    message("Bootstrap parallel config: Using user-provided settings")
  }

  # Report configuration to user
  max_cores <- parallel::detectCores()
  message("System max cores available: ", max_cores)
  message("Bootstrap will use: ", resolved_args$workers, " workers with '",
          resolved_args$plan, "' plan")

  resolved_args
}



#' Get all exported functions from ForestSearch namespace
#' @keywords internal
get_bootstrap_exports <- function() {
  # Automatically discover all exported functions
  ns <- asNamespace("forestsearch")
  ls(ns, all.names = FALSE)
}


#' Check that required packages are installed
#'
#' @param pkgs Character vector of package names.
#' @return Invisible TRUE if all packages are available.
#' @keywords internal
#' @noRd
ensure_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Required package(s) not installed: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
