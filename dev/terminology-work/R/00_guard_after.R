#!/usr/bin/env Rscript
# Guard: AFTER snapshot + verdict.  Fails loudly if any package source changed.

source("dev/efficiency-eval/R/00_guard.R")

before <- readRDS("dev/terminology-work/out/hash_before.rds")
after  <- fs_hash_sources()
saveRDS(after, "dev/terminology-work/out/hash_after.rds")

v <- fs_guard_verify(before, after)
cat("GUARD VERDICT\n")
cat("  ok      :", format(v$ok), "\n")
cat("  n_files :", v$n_files, "\n")
cat("  note    :", v$note, "\n")
