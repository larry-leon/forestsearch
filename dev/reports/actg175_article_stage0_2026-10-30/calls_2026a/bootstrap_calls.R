# ---- transplant (verbatim chunk bodies; fences removed) ----
# source: forestsearch/quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd
# forestsearch commit 43b051b6 (2026-08-16 12:07; the render that wrote the payload at 2026-08-16 16:50:40 UTC ran this text — its HTML is committed in the same commit; identical to HEAD; differs from fd1ce74f, the last commit before the build, in fs_selection_rule both->neighborhood)
#   chunk 'bootstrap': fence lines 770-791
#   chunk 'fs-dina-bootstrap': fence lines 1277-1300
#   chunk 'fs-grf-bootstrap': fence lines 1501-1523
# ----------------------------------------------------------

## ==== chunk 'bootstrap' ====
t0 <- proc.time()

fs_bc <- forestsearch_bootstrap_dofuture(
  fs.est = fs,
  nb_boots = NB,
  show_three = FALSE,
  details = FALSE,
  digits = construction_digits,
  parallel_args = list(
    plan = "multisession",
    workers = n_cores,
    show_message = TRUE
  )
)

plan("sequential"); gc(verbose = FALSE)
timings$bootstrap <- (proc.time() - t0)["elapsed"]

cat("\nBootstrap completed in",
    round(timings$bootstrap / 60, 1), "minutes\n")

## ==== chunk 'fs-dina-bootstrap' ====
# Runs only when DINA identified a subgroup.  The bootstrap bias-corrects an
# identified subgroup's effect, so there is nothing to do otherwise.
t0 <- proc.time()

fs_dina_bc <- forestsearch_bootstrap_dofuture(
  fs.est = fs_dina,
  nb_boots = NB,
  show_three = FALSE,
  details = FALSE,
  digits = construction_digits,
  parallel_args = list(
    plan = "multisession",
    workers = n_cores,
    show_message = TRUE
  )
)

plan("sequential"); gc(verbose = FALSE)
timings$fs_dina_bootstrap <- (proc.time() - t0)["elapsed"]

cat("\nDINA-mode bootstrap completed in",
    round(timings$fs_dina_bootstrap / 60, 1), "minutes\n")

## ==== chunk 'fs-grf-bootstrap' ====
# Runs only when GRF identified a subgroup.
t0 <- proc.time()

fs_grf_bc <- forestsearch_bootstrap_dofuture(
  fs.est = fs_grf,
  nb_boots = NB,
  show_three = FALSE,
  details = FALSE,
  digits = construction_digits,
  parallel_args = list(
    plan = "multisession",
    workers = n_cores,
    show_message = TRUE
  )
)

plan("sequential"); gc(verbose = FALSE)
timings$fs_grf_bootstrap <- (proc.time() - t0)["elapsed"]

cat("\nGRF-mode bootstrap completed in",
    round(timings$fs_grf_bootstrap / 60, 1), "minutes\n")

