# =============================================================================
# Read-only guard + environment capture
# =============================================================================
# Enforces the harness's central claim: the package source is not modified.
# Hashes are taken before and after the run and compared.

#' Hash every R source file the evaluation could plausibly touch
#'
#' Covers the installed package directory and, when the harness is run from
#' inside a clone, the git worktree's `R/` as well.
#'
#' @return Named character vector of SHA-256 digests, sorted by path.
fs_hash_sources <- function() {
  paths <- character(0)

  inst <- tryCatch(system.file(package = "forestsearch"), error = function(e) "")
  if (nzchar(inst) && dir.exists(inst)) {
    paths <- c(paths, list.files(inst, recursive = TRUE, full.names = TRUE))
  }

  # a sibling / parent worktree, if the eval lives inside the repo
  for (cand in c("R", "../R", "../../R")) {
    if (dir.exists(cand)) {
      paths <- c(paths, list.files(cand, pattern = "[.][Rr]$", full.names = TRUE))
    }
  }

  paths <- unique(normalizePath(paths[file.exists(paths)], mustWork = FALSE))
  if (!length(paths)) return(character(0))

  h <- vapply(paths, function(p) {
    tryCatch(digest::digest(file = p, algo = "sha256"),
             error = function(e) NA_character_)
  }, character(1))
  h[order(names(h))]
}

#' Compare two hash snapshots and stop on any difference
fs_guard_verify <- function(before, after) {
  if (!length(before) && !length(after)) {
    return(list(ok = NA, note = "no source files located to hash"))
  }
  added   <- setdiff(names(after), names(before))
  removed <- setdiff(names(before), names(after))
  common  <- intersect(names(before), names(after))
  changed <- common[!identical(before[common], after[common]) &
                      before[common] != after[common]]

  if (length(added) || length(removed) || length(changed)) {
    stop("READ-ONLY CONTRACT VIOLATED. ",
         "changed: ", length(changed), " added: ", length(added),
         " removed: ", length(removed), "\n",
         paste(utils::head(c(changed, added, removed), 10), collapse = "\n"))
  }
  list(ok = TRUE, n_files = length(common),
       note = sprintf("%d source files verified unchanged", length(common)))
}

#' Capture the execution environment so results are interpretable later
fs_capture_env <- function() {
  bl <- tryCatch(sessionInfo()$BLAS, error = function(e) NA_character_)
  data.frame(
    stringsAsFactors = FALSE,
    timestamp     = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    r_version     = paste(R.version$major, R.version$minor, sep = "."),
    platform      = R.version$platform,
    blas          = if (is.null(bl)) NA_character_ else bl,
    detect_cores  = parallel::detectCores(),
    avail_cores   = tryCatch(unname(parallelly::availableCores()),
                             error = function(e) NA_integer_),
    forestsearch  = tryCatch(as.character(utils::packageVersion("forestsearch")),
                             error = function(e) NA_character_),
    survival      = as.character(utils::packageVersion("survival")),
    data_table    = as.character(utils::packageVersion("data.table")),
    future        = tryCatch(as.character(utils::packageVersion("future")),
                             error = function(e) NA_character_)
  )
}
