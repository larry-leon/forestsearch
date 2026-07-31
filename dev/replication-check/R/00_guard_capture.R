# Guard snapshot + environment capture (read-only contract)
source("dev/efficiency-eval/R/00_guard.R")
which <- commandArgs(trailingOnly = TRUE)[1]
h <- fs_hash_sources()
saveRDS(h, file.path("dev/replication-check/out", paste0("guard_", which, ".rds")))
cat(which, ": hashed", length(h), "files\n")
if (identical(which, "before")) {
  env <- fs_capture_env()
  saveRDS(env, "dev/replication-check/out/env_capture.rds")
  print(t(env))
}
if (identical(which, "after")) {
  b <- readRDS("dev/replication-check/out/guard_before.rds")
  v <- fs_guard_verify(b, h)
  saveRDS(v, "dev/replication-check/out/guard_verdict.rds")
  str(v)
}
