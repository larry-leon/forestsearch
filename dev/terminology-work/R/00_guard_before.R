#!/usr/bin/env Rscript
# Guard: take the BEFORE snapshot of package source hashes.
# Reuses fs_hash_sources() from dev/efficiency-eval/R/00_guard.R per the brief.
# Writes only under dev/terminology-work/.

source("dev/efficiency-eval/R/00_guard.R")
dir.create("dev/terminology-work/out", showWarnings = FALSE, recursive = TRUE)

before <- fs_hash_sources()
env    <- fs_capture_env()

saveRDS(before, "dev/terminology-work/out/hash_before.rds")
write.csv(env, "dev/terminology-work/out/env_capture.csv", row.names = FALSE)

cat(sprintf("BEFORE snapshot: %d files hashed\n", length(before)))
print(env)
