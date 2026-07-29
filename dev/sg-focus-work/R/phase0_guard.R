# =============================================================================
# Phase 0 -- read-only guard
# =============================================================================
# Reuses fs_hash_sources() / fs_guard_verify() / fs_capture_env() from the
# earlier efficiency evaluation rather than rewriting them.  Must be sourced
# with the repository root as the working directory (fs_hash_sources() probes
# "R", "../R", "../../R" relative to getwd()).

source("dev/efficiency-eval/R/00_guard.R")

PHASE0_GUARD_RDS <- "dev/sg-focus-work/out/guard_before.rds"

phase0_guard_snapshot <- function(path = PHASE0_GUARD_RDS) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  h <- fs_hash_sources()
  saveRDS(list(hashes = h, env = fs_capture_env(),
               timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), path)
  cat(sprintf("guard: snapshot of %d source files -> %s\n", length(h), path))
  invisible(h)
}

phase0_guard_check <- function(path = PHASE0_GUARD_RDS) {
  before <- readRDS(path)$hashes
  after  <- fs_hash_sources()
  v <- fs_guard_verify(before, after)
  cat(sprintf("guard: %s\n", v$note))
  invisible(v)
}
