# =============================================================================
# run_fs_profile.R
#
# One thing to run to profile a single forestsearch() call from the
# actg175_binary_m1_harm_single-simtrial qmd.
#
# PREREQUISITE
#   In RStudio, open the qmd and run the chunks top-to-bottom THROUGH the
#   `single-trial-alt-sim` chunk (Run > Run All, or Ctrl+Alt+P up to it).
#   That attaches the package and defines df_sim, confounders, and fs_params
#   in your global session. Then run this file:
#
#       source("run_fs_profile.R")
#
#   Place fs_profile_harness.R and this file in the SAME folder as the qmd
#   (Quarto's working directory during/after a render is the qmd's folder).
#
# NOTE ON REGIME
#   The harness forces plan = "sequential" so Rprof can actually see the work
#   (multisession runs it in worker processes the profiler can't sample) and
#   so it matches the simulation/bootstrap inner-loop regime in fs3b_params.
#   This profiling run will therefore be SLOWER in wall-clock than the doc's
#   multisession single call -- that is expected. The phase PERCENTAGES are
#   the deliverable, not the absolute seconds of this run.
# =============================================================================

source("fs_profile_harness.R")

# -- Confirm the session has what we need (run the qmd prep first if not). ----
missing <- c("df_sim", "confounders", "fs_params")[
  !vapply(c("df_sim", "confounders", "fs_params"), exists, logical(1))
]
if (length(missing)) {
  stop("Missing object(s): ", paste(missing, collapse = ", "),
       ".\n  Run the qmd chunks through `single-trial-alt-sim` first, then ",
       "source() this file.")
}

# -- Match the cost-relevant setting of the real single-trial call. -----------
#    (The doc sets collapse_cuts = TRUE on fs_params_single; everything else
#    that differs -- multisession, show_candidate_summary, details -- the
#    harness deliberately overrides for a clean, sequential sample.)
fs_params_prof <- fs_params
fs_params_prof$collapse_cuts <- TRUE

# -- Profile. -----------------------------------------------------------------
#    reps = 2 averages out sampling noise. For a quick first look, set reps = 1
#    (a single multi-second call already yields thousands of samples, so the
#    phase shares are stable). nb_boots = 500 matches the doc's bootstrap and
#    is used only to extrapolate the per-trial budget -- it does NOT run it.
prof <- profile_forestsearch(
  df          = df_sim,
  confounders = confounders,
  fs_params   = fs_params_prof,
  reps        = 2L,
  interval    = 0.01,
  nb_boots    = 500L,
  top_n       = 20L
)

# `prof` holds the structured result (prof$phase_table, prof$call_seconds,
# prof$raw_total, prof$raw_self, prof$diagnostics) if you want to inspect or
# save it. The readable report has already printed above.
invisible(prof)
