# =============================================================================
# STAGE P -- PARALLEL BEHAVIOUR  (the measurements a small machine cannot make)
# =============================================================================
# Everything here uses the package's OWN public arguments. No code is patched.
#
# P2 is the headline. forestsearch() defaults stop_threshold = 0.95
# (forestsearch_main.R:887). With the default sg_focus = "hr",
# subgroup_consistency_main.R:837-840 sets batch_size_parallel <- 1L, so each
# future_lapply() call carries ONE candidate: workers idle, and the closure --
# which captures df, found.hrs, index.Z, confs_labels -- is serialized per
# candidate. That is a deliberate tradeoff for fine-grained early stopping.
# Whether it pays depends on where early stopping actually fires on YOUR data,
# which is exactly what P2 measures.
# =============================================================================

.p_row <- function(id, what, x, xlab, seconds, extra = "") {
  data.frame(stringsAsFactors = FALSE,
             id = id, what = what, x = x, xlab = xlab,
             seconds = seconds, extra = extra)
}

# Pull fields out of a forestsearch result without assuming the layout.
.fs_get <- function(obj, field) {
  if (is.null(obj)) return(NA)
  if (!is.null(obj[[field]])) return(obj[[field]])
  for (nm in names(obj)) {
    el <- obj[[nm]]
    if (is.list(el) && !is.null(el[[field]])) return(el[[field]])
  }
  NA
}

.fs_label <- function(obj) {
  sg <- .fs_get(obj, "sg.harm")
  if (is.null(sg) || (length(sg) == 1 && is.na(sg[1]))) return(NA_character_)
  paste(sort(as.character(sg)), collapse = " & ")
}


# -----------------------------------------------------------------------------
# P1 -- worker-pool spawn cost vs W
# -----------------------------------------------------------------------------
# subgroup_search.R:143/160 sets plan(multisession) then resets to
# plan(sequential), rather than the old_plan/on.exit pattern used in
# forestsearch_cross_validation.R:302 and subgroup_consistency_main.R:821.
# The reset tears the pool down, so the next call respawns it. Cost scales
# with W and is invisible on a small machine.
par_spawn_cost <- function(workers_sweep, reps = 3L) {
  out <- list()
  trivial <- function() future.apply::future_lapply(1:4, function(i) i)
  for (W in workers_sweep) {
    ts <- numeric(reps)
    for (r in seq_len(reps)) {
      future::plan(future::sequential)
      ts[r] <- system.time({
        future::plan(future::multisession, workers = W)
        invisible(trivial())
      })[["elapsed"]]
    }
    future::plan(future::sequential)
    out[[length(out) + 1L]] <- .p_row(
      "P1", "cold pool spawn + first dispatch", W, "workers",
      stats::median(ts), "cost paid per subgroup.search() call")
  }
  do.call(rbind, out)
}


# -----------------------------------------------------------------------------
# P2 -- batch_size sweep  (HEADLINE)
# -----------------------------------------------------------------------------
# Sweeps parallel_args$batch_size through the package's public interface with
# stop_threshold at its default. Records wall-clock, how many candidates were
# actually evaluated, whether early stopping fired, and -- critically -- the
# SELECTED SUBGROUP, so any speed difference can be checked against
# selection-identity rather than assumed harmless.
par_batch_sweep <- function(df, fs_args, W, batch_sizes, reps = 3L) {
  out <- list()
  for (bs in batch_sizes) {
    for (rep in seq_len(reps)) {
      args <- c(list(df.analysis = df), fs_args)
      args$parallel_args <- list(plan = "multisession", workers = W,
                                 batch_size = as.integer(bs))
      gc()
      tt <- NA_real_; res <- NULL
      tt <- tryCatch(system.time(
        res <- do.call(forestsearch::forestsearch, args))[["elapsed"]],
        error = function(e) { message("  batch_size=", bs, " rep=", rep,
                                      " failed: ", conditionMessage(e));
                              NA_real_ })
      n_eval <- .fs_get(res, "n_candidates_evaluated")
      es     <- .fs_get(res, "early_stop_triggered")
      out[[length(out) + 1L]] <- .p_row(
        "P2", "forestsearch() wall-clock by batch_size", bs, "batch_size", tt,
        sprintf("rep=%d evaluated=%s early_stop=%s selected={%s}",
                rep, format(n_eval), format(es), .fs_label(res)))
    }
  }
  future::plan(future::sequential)
  do.call(rbind, out)
}


# -----------------------------------------------------------------------------
# P3 -- consistency-stage scaling vs worker count
# -----------------------------------------------------------------------------
# Does the stage actually scale, or does serialization / idle-worker overhead
# flatten the curve? Run with a batch_size that permits parallelism (= W) so
# the curve reflects the stage rather than the batch_size = 1 degenerate case.
par_scaling <- function(df, fs_args, workers_sweep) {
  out <- list()
  for (W in workers_sweep) {
    args <- c(list(df.analysis = df), fs_args)
    args$parallel_args <- list(plan = "multisession", workers = W,
                               batch_size = as.integer(W))
    gc()
    tt <- tryCatch(system.time(
      res <- do.call(forestsearch::forestsearch, args))[["elapsed"]],
      error = function(e) NA_real_)
    out[[length(out) + 1L]] <- .p_row(
      "P3", "forestsearch() wall-clock by workers", W, "workers", tt,
      "batch_size = workers")
  }
  future::plan(future::sequential)
  do.call(rbind, out)
}


# -----------------------------------------------------------------------------
# P4 -- closure serialization cost
# -----------------------------------------------------------------------------
# The evaluation closure captures df and found.hrs. At batch_size = 1 it is
# shipped once per candidate. Measure the payload and the round-trip so the
# P2 result can be explained rather than merely observed.
par_serialization <- function(df, n_candidates = 200L) {
  payload <- list(df = df,
                  index.Z = matrix(0L, nrow = n_candidates, ncol = 8L),
                  names.Z = paste0("q", 1:8),
                  found.hrs = data.frame(HR = numeric(n_candidates),
                                         n = numeric(n_candidates),
                                         E = numeric(n_candidates),
                                         grp = numeric(n_candidates)))
  bytes <- length(serialize(payload, NULL))
  t_ser <- stats::median(replicate(20,
    system.time(serialize(payload, NULL))[["elapsed"]]))

  rbind(
    .p_row("P4a", "closure payload size", nrow(df), "rows in df",
           NA_real_, sprintf("%.2f MB per dispatch", bytes / 1024^2)),
    .p_row("P4b", "serialize() cost per dispatch", nrow(df), "rows in df",
           t_ser,
           sprintf("x %d candidates at batch_size=1 -> %.1f s of pure serialization",
                   n_candidates, t_ser * n_candidates))
  )
}


run_parallel <- function(cfg, df) {
  W_sweep <- cfg$workers_sweep
  if (is.null(W_sweep)) {
    mx <- tryCatch(unname(parallelly::availableCores()),
                   error = function(e) parallel::detectCores())
    W_sweep <- unique(pmax(1L, c(1L, mx %/% 16L, mx %/% 4L, mx %/% 2L, mx)))
  }
  W_main <- cfg$workers_main %||% max(W_sweep)

  bs <- cfg$batch_sweep
  if (is.null(bs)) {
    bs <- unique(pmax(1L, c(1L, W_main %/% 8L, W_main %/% 2L,
                            W_main, 2L * W_main)))
  }

  message("  P1: spawn cost over workers = ", paste(W_sweep, collapse = ", "))
  r1 <- par_spawn_cost(W_sweep)

  message("  P2: batch_size sweep = ", paste(bs, collapse = ", "),
          " at W = ", W_main)
  r2 <- par_batch_sweep(df, cfg$fs_args, W_main, bs)

  message("  P3: scaling over workers")
  r3 <- par_scaling(df, cfg$fs_args, W_sweep)

  message("  P4: serialization")
  r4 <- par_serialization(df, cfg$fs_args$max_subgroups_search %||% 200L)

  rbind(r1, r2, r3, r4)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
