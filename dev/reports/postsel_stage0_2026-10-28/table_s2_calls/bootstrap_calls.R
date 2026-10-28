# ---------------------------------------------------------------------------
# Table S2 of León & Anderson (2026a), Full-bootstrap column: forestsearch_bootstrap_dofuture() for FS, DINA, GRF -- verbatim transplant.
# Source: forestsearch quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd
# Commit: cf4d64328b6288ce2ec009f2c2fb15689de6a4cb (2026-08-15 21:12 -0700), the last
#   commit touching the .qmd before the Table S2 payload (built 2026-08-15 22:50:30 PDT,
#   Mac-Studio-3.local, forestsearch 0.2.0); render html committed in 43b051b6.
#   At HEAD the .qmd differs from cf4d6432 only at lines 144-145 (run_cv / run_loo
#   FALSE -> TRUE, 6ad5a724); no transplanted chunk below is affected.
# Chunks (name:fence lines, inclusive): bootstrap:931-963 fs-dina-bootstrap:1563-1597 fs-grf-bootstrap:1802-1832
# Chunk bodies are copied verbatim; only these header comments and the
# '## ---- chunk: <name> ----' separators were added. Transplanted 2026-10-28.
# ---------------------------------------------------------------------------

## ---- chunk: bootstrap (lines 932-962) ----
#| message: false
# message: false suppresses forestsearch_bootstrap_dofuture()'s per-iteration
# progress and its "Low events in NEW subgroup" notes (both message(), not
# cat()).  Roughly 200 lines of console noise in a NB = 1000 run; the
# low-event condition is summarised properly in the diagnostics chunk below.

t0 <- proc.time()

fs_bc <- forestsearch_bootstrap_dofuture(
  # Pinned to the package default.  The fit carries mr_inference = run_mr
  # (TRUE here), and each replicate/fold call copies the fit's argument
  # list, so without this the flag is inherited and MR runs inside every
  # replicate.  Nothing downstream consumes a per-replicate MR object.
  # No-op as written.
  mr_in_replicates = FALSE,
  fs.est = fs,
  nb_boots = NB,
  show_three = FALSE,
  details = FALSE,
  parallel_args = list(
    plan = "multisession",
    workers = n_cores,
    show_message = TRUE
  )
)

plan("sequential")
timings$bootstrap <- (proc.time() - t0)["elapsed"]

cat("\nBootstrap completed in",
    round(timings$bootstrap / 60, 1), "minutes\n")

## ---- chunk: fs-dina-bootstrap (lines 1564-1596) ----
#| message: false
# message: false suppresses forestsearch_bootstrap_dofuture()'s per-iteration
# progress and its "Low events in NEW subgroup" notes (both message(), not
# cat()).  Roughly 200 lines of console noise in a NB = 1000 run; the
# low-event condition is summarised properly in the diagnostics chunk below.
# Runs only when DINA identified a subgroup.  The bootstrap bias-corrects an
# identified subgroup's effect, so there is nothing to do otherwise.
t0 <- proc.time()

# details = TRUE keeps per-iteration output visible while validating DINA-mode
# bootstrap.  The plan is "multisession" (production speed); switch it to
# "sequential" if errors need to surface in loop order.

fs_dina_bc <- forestsearch_bootstrap_dofuture(
  # Pinned to the package default, as at the main bootstrap above -- see the
  # comment there for why this is a no-op.
  mr_in_replicates = FALSE,
  fs.est = fs_dina,
  nb_boots = NB,
  show_three = FALSE,
  details = TRUE,
  parallel_args = list(
    plan = "multisession",
    workers = n_cores,
    show_message = TRUE
  )
)

plan("sequential")
timings$fs_dina_bootstrap <- (proc.time() - t0)["elapsed"]

cat("\nDINA-mode bootstrap completed in",
    round(timings$fs_dina_bootstrap / 60, 1), "minutes\n")

## ---- chunk: fs-grf-bootstrap (lines 1803-1831) ----
#| message: false
# message: false suppresses forestsearch_bootstrap_dofuture()'s per-iteration
# progress and its "Low events in NEW subgroup" notes (both message(), not
# cat()).  Roughly 200 lines of console noise in a NB = 1000 run; the
# low-event condition is summarised properly in the diagnostics chunk below.
# Runs only when GRF identified a subgroup.  The bootstrap bias-corrects an
# identified subgroup's effect, so there is nothing to do otherwise.
t0 <- proc.time()

fs_grf_bc <- forestsearch_bootstrap_dofuture(
  # Pinned to the package default, as at the main bootstrap above -- see the
  # comment there for why this is a no-op.
  mr_in_replicates = FALSE,
  fs.est = fs_grf,
  nb_boots = NB,
  show_three = FALSE,
  details = TRUE,
  parallel_args = list(
    plan = "multisession",
    workers = n_cores,
    show_message = TRUE
  )
)

plan("sequential")
timings$fs_grf_bootstrap <- (proc.time() - t0)["elapsed"]

cat("\nGRF-mode bootstrap completed in",
    round(timings$fs_grf_bootstrap / 60, 1), "minutes\n")

