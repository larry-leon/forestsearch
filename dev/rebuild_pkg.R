# rebuild_pkg.R
# Minimal, reliable full reinstall that avoids the
# "lazy-load database '.../pkg.rdb' is corrupt" error on macOS.
#
# The corruption comes from R CMD INSTALL overwriting a package's .rdb while
# that package is still loaded in a session (the live session holds the old
# .rdx index). The fix is simply: unload it first, then install. That is all
# this does -- no R restart required.
#
#   source("rebuild_pkg.R")
#   rebuild_pkg("/path/to/pkg")                # normal rebuild
#   rebuild_pkg("/path/to/pkg", fresh = TRUE)  # wipe an already-corrupt copy first

rebuild_pkg <- function(path = ".", document = TRUE, fresh = FALSE) {
  path <- normalizePath(path, mustWork = TRUE)
  pkg  <- unname(read.dcf(file.path(path, "DESCRIPTION"), fields = "Package")[1, 1])

  # unload so we never replace a .rdb this session still has open
  if (isNamespaceLoaded(pkg)) try(pkgload::unload(pkg), silent = TRUE)

  # clear a stale 00LOCK; with fresh = TRUE also delete an already-corrupt
  # installed copy so nothing bad survives into the reinstall
  for (lib in .libPaths()) {
    unlink(file.path(lib, paste0("00LOCK-", pkg)), recursive = TRUE, force = TRUE)
    if (fresh) unlink(file.path(lib, pkg), recursive = TRUE, force = TRUE)
  }

  if (document) devtools::document(path)
  devtools::install(path, args = "--no-byte-compile",
                    upgrade = "never", quick = FALSE)
}
