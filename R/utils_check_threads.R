# =============================================================================
# Thread cap for grf under R CMD check
# =============================================================================
#
# grf trains forests with its own C++ thread pool (std::thread), which ignores
# future plans and OMP_THREAD_LIMIT. With num.threads = NULL grf uses every
# available core, so an example run under CRAN's check can exceed the 2-core
# CPU-to-elapsed limit. Under _R_CHECK_LIMIT_CORES_ the cap is 2 threads;
# otherwise NULL is returned and grf's own default applies unchanged.
#
# =============================================================================


#' Number of threads to pass to grf forest training calls
#'
#' @param user Optional thread count supplied by the caller. When non-NULL
#'   it is returned unchanged.
#' @return `user` when supplied; otherwise `2L` when the environment variable
#'   `_R_CHECK_LIMIT_CORES_` is set to anything other than `""`, `"false"` or
#'   `"FALSE"`; otherwise `NULL` (grf's default).
#' @noRd
.fs_grf_num_threads <- function(user = NULL) {
  if (!is.null(user)) return(user)
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (!chk %in% c("", "false", "FALSE")) return(2L)
  NULL
}
