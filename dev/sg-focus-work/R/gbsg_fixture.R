# =============================================================================
# Shared GBSG fixture for sg_focus acceptance checks (dev scratch)
# =============================================================================
# LASSO/GRF off and a sequential plan so runs are fast and deterministic;
# conf_force guarantees a non-empty candidate family. Mirrors the GBSG setup in
# vignettes/forestsearch.Rmd.

# NB: gbsg ships with `survival`, NOT with forestsearch -- despite the comment
# in vignettes/forestsearch.Rmd claiming otherwise.  In survival::gbsg `grade`
# is an integer, so compare numerically.
gbsg_df <- function() {
  d <- survival::gbsg
  d$id          <- seq_len(nrow(d))
  d$time_months <- d$rfstime / 30.4375
  d$grade3      <- as.integer(d$grade == 3)
  d
}

GBSG_CONFS <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")

#' One forestsearch() run on GBSG, with `...` overriding any default below.
#'
#' Overrides are applied with `args[nm] <- list(value)`, NOT
#' utils::modifyList(): modifyList treats a NULL value as "delete this element",
#' so `gbsg_fs(stop_threshold = NULL)` would silently be a no-op and the
#' package default (0.95) would still apply.  Passing NULL explicitly is
#' essential to the stop_threshold acceptance checks, so it must survive.
gbsg_fs <- function(...) {
  base <- list(
    df.analysis            = gbsg_df(),
    outcome.name           = "time_months",
    event.name             = "status",
    treat.name             = "hormon",
    id.name                = "id",
    confounders.name       = GBSG_CONFS,
    is.RCT                 = TRUE,
    seedit                 = 8316951,
    quiet                  = TRUE,
    use_lasso              = FALSE,
    use_grf                = FALSE,
    use_twostage           = FALSE,
    debias_gate            = FALSE,
    consistency_method     = "resample",
    hr.threshold           = 1.25,
    hr.consistency         = 1.00,
    pconsistency.threshold = 0.90,
    maxk                   = 2L,
    n.min                  = 60L,
    d0.min                 = 12L,
    d1.min                 = 12L,
    conf_force             = c("er <= 0", "pgr <= 0"),
    parallel_args          = list(plan = "sequential", workers = 1L)
  )
  over <- list(...)
  for (nm in names(over)) base[nm] <- list(over[[nm]])
  do.call(forestsearch::forestsearch, base)
}

#' The comparable part of a fit: the selected subgroup and its scored row.
#' Excludes timings and the echoed call, which differ by construction.
gbsg_selection <- function(fs) {
  gc_ <- fs$grp.consistency
  res <- gc_$out_sg$result
  list(
    sg_def   = if (is.null(fs$sg.harm)) NA_character_ else paste(fs$sg.harm, collapse = " & "),
    n_sel    = if (is.null(gc_$sg.harm.id)) NA_integer_ else sum(gc_$sg.harm.id == 1L),
    top_row  = if (is.null(res) || !nrow(res)) NULL else
                 as.list(as.data.frame(res)[1, c("Pcons", "hr", "N", "K", "m")]),
    n_eval   = gc_$n_candidates_evaluated,
    n_total  = gc_$n_candidates_total,
    n_passed = gc_$n_passed,
    early    = isTRUE(gc_$early_stop_triggered)
  )
}
