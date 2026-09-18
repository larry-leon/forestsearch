# =============================================================================
# probe_threshold_sync_2026-09-18.R -- dev-side entry point
#
# The probe itself lives with the acceptance tests, at
# tests/testthat/helper-threshold-sync.R, so there is exactly one copy and the
# tests and the dev-side runs cannot drift.  This file only sources it.
#
#   source("dev/tasks/probe_threshold_sync_2026-09-18.R")
#   probe_threshold_sync(sync = TRUE)
# =============================================================================
source(file.path("tests", "testthat", "helper-threshold-sync.R"))
