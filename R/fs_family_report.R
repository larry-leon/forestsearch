# =============================================================================
# fs_family_report(): what is, and is not, deterministic about the candidate
# family that forestsearch() will search, for a given set of arguments.
# -----------------------------------------------------------------------------
# A reporter, not a setter.  It computes nothing that enters any result and is
# called by nothing else in the package.  Every classification below quotes
# the governing line of the engine it mirrors; if the engine changes, the
# drift-guard test (tests/testthat/test-fs-family-report.R, test 2) fails.
# =============================================================================

#' Report which stages of the candidate family are data-dependent
#'
#' A user who sets \code{use_lasso = FALSE}, \code{use_grf = FALSE},
#' \code{use_dina = FALSE} and \code{vi.grf.min = NULL} may reasonably believe
#' the candidate family \code{\link{forestsearch}} searches is now fixed.  It
#' is not.  Continuous cuts are placed at \emph{sample} quantiles
#' (\code{get_FSdata()}); the prevalence, redundancy and size floors
#' (\code{minp}, \code{rmin}, \code{n.min}) are applied to \emph{sample}
#' counts (\code{subgroup.search()}); and the consistency stage removes
#' near-duplicate candidates keyed on their fitted \emph{statistics}, with no
#' argument that turns it off (\code{remove_near_duplicate_subgroups()}).
#' Those stages are the method, not screening.  This function says so, for a
#' given argument set, in a form that cannot overpromise.
#'
#' \strong{It reports and changes nothing.}  It never fits a model, never
#' runs the search, and is called by nothing else in the package.  \strong{No
#' combination of arguments makes the family deterministic while cuts are
#' placed at sample quantiles}; the best a caller can do is switch off the
#' disableable stages, and the printed footer lists what remains.
#'
#' The report mirrors \code{forestsearch()}'s own argument resolution rather
#' than reading arguments at face value: \code{sg_focus} aliases are
#' normalised; \code{sg_focus = "maxeff"} zeroes \code{minp}, \code{rmin} and
#' \code{pconsistency.threshold}, sets \code{stop_threshold} to \code{NULL},
#' \code{use_twostage} to \code{FALSE} and \code{max_subgroups_search} to
#' \code{Inf}, and disables the effect floor; \code{stop_threshold} is reset
#' to \code{NULL} for every focus other than \code{"maxeffCons"} (it is
#' meaningful only there); \code{max_n_confounders} is applied only inside the
#' GRF variable-importance block, so it is inert when \code{vi.grf.min} is
#' \code{NULL}; and the per-arm event floors \code{d0.min} / \code{d1.min} are
#' skipped for continuous and count outcomes.
#'
#' @param x Either a named list of \code{\link{forestsearch}} arguments (as
#'   one would pass to \code{do.call(forestsearch, x)}; unspecified arguments
#'   take \code{forestsearch()}'s formal defaults), or a fitted
#'   \code{forestsearch} object, in which case its \code{args_call_all} is
#'   read.
#' @param data Optional data frame.  When supplied, the report is grounded in
#'   counts: the number of cut columns \code{get_FSdata()} produces on this
#'   data under these arguments, and the number of up-to-\code{maxk}
#'   combinations that implies.  Nothing is fitted; the search is not run.
#' @param outcome_type One of \code{"survival"}, \code{"binary"},
#'   \code{"continuous"}, \code{"count"}.  Taken from \code{x} when present
#'   there; required otherwise, because the per-arm floor row depends on it.
#'
#' @return A data frame of class \code{c("fs_family_report", "data.frame")}
#'   with one row per stage and columns \code{stage}, \code{arguments} (the
#'   governing \code{forestsearch()} formals), \code{values} (their resolved
#'   values, formatted), \code{status} (one of \code{"deterministic"},
#'   \code{"disabled"}, \code{"inert"}, \code{"data-dependent"},
#'   \code{"data-dependent (not disableable)"}) and \code{note} (one line:
#'   why, and what would change it; \code{NA} where nothing would).
#'   Attributes: \code{verdict} (one sentence), \code{status_counts} (a
#'   named integer vector), \code{data_supplied} (logical), and -- when
#'   \code{data} is supplied -- \code{n_cut_columns} and
#'   \code{n_combinations}.
#'
#' @examples
#' rep <- fs_family_report(
#'   list(confounders.name = c("age", "preanti", "hemo"),
#'        conf.cont_jcuts = list(age = 10, preanti = 10),
#'        n.min = 60, maxk = 2, sg_focus = "maxeffCons",
#'        use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
#'        vi.grf.min = NULL),
#'   outcome_type = "continuous")
#' print(rep)
#' attr(rep, "verdict")
#'
#' # grounded in counts: cut columns and combinations on a data frame
#' set.seed(1)
#' df <- data.frame(age = round(rnorm(200, 40, 8)),
#'                  preanti = round(rexp(200, 1 / 500)),
#'                  hemo = rbinom(200, 1L, 0.1))
#' rep2 <- fs_family_report(
#'   list(confounders.name = c("age", "preanti", "hemo"),
#'        conf.cont_jcuts = list(age = 10, preanti = 10), maxk = 2,
#'        use_lasso = FALSE, use_grf = FALSE, vi.grf.min = NULL),
#'   data = df, outcome_type = "continuous")
#' attr(rep2, "n_cut_columns"); attr(rep2, "n_combinations")
#'
#' @seealso \code{\link{forestsearch}}, \code{\link{fs_oc_family_enumerate}}
#'   (the population-frame enumeration used by the operating-characteristic
#'   wrapper, which is deterministic in the DGM precisely because it does not
#'   cut at sample quantiles).
#' @export
fs_family_report <- function(x, data = NULL, outcome_type = NULL) {

  # ---- 1. what was passed: a fitted object or an argument list --------------
  if (inherits(x, "forestsearch") ||
      (is.list(x) && !is.null(x$args_call_all) && is.list(x$args_call_all))) {
    args <- x$args_call_all
    if (is.null(outcome_type)) outcome_type <- x$outcome_type %||% args$outcome_type
  } else if (is.list(x) && !is.data.frame(x) &&
             (length(x) == 0L || !is.null(names(x)))) {
    args <- x
    if (is.null(outcome_type)) outcome_type <- args$outcome_type
  } else {
    stop("`x` must be a named list of forestsearch() arguments or a fitted ",
         "forestsearch() object (one carrying `args_call_all`).", call. = FALSE)
  }
  unknown <- setdiff(names(args), names(formals(forestsearch)))
  unknown <- unknown[nzchar(unknown)]
  if (length(unknown)) {
    stop("`x` contains names that are not forestsearch() formals: ",
         paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  if (is.null(outcome_type)) {
    stop("`outcome_type` is required (it is not in `x`): the per-arm floor ",
         "row differs by outcome type.", call. = FALSE)
  }
  outcome_type <- match.arg(outcome_type,
                            c("survival", "binary", "continuous", "count"))

  # ---- 2. resolve every classified argument as forestsearch() does ---------
  fs_formals <- formals(forestsearch)
  .arg <- function(nm) {
    if (nm %in% names(args)) return(args[[nm]])
    d <- fs_formals[[nm]]
    if (is.symbol(d) && !nzchar(as.character(d))) return(NULL)   # no default
    tryCatch(eval(d, envir = list2env(args, parent = environment())),
             error = function(e) d)
  }
  # sg_focus aliases -> canonical (forestsearch_helpers.R: .normalize_sg_focus)
  sg_focus_raw <- .arg("sg_focus")
  sg_focus <- .normalize_sg_focus(sg_focus_raw)
  subgroup_method <- .arg("subgroup_method")
  if (is.null(subgroup_method)) subgroup_method <- "consistency"
  subgroup_method <- as.character(subgroup_method)[1L]
  is_maxeff <- identical(as.character(sg_focus), "maxeff") &&
               identical(subgroup_method, "consistency")

  use_lasso <- isTRUE(.arg("use_lasso")); use_grf <- isTRUE(.arg("use_grf"))
  use_dina  <- isTRUE(.arg("use_dina"))
  vi.grf.min <- .arg("vi.grf.min")
  max_n_confounders <- .arg("max_n_confounders")
  maxk <- .arg("maxk")
  minp <- .arg("minp")
  # rmin is not a forestsearch() formal: subgroup.search()'s default (5)
  # reaches the search unless sg_focus = "maxeff" (search_overrides$rmin <- 0)
  rmin <- eval(formals(subgroup.search)[["rmin"]])
  n.min <- .arg("n.min"); n.min.frac <- .arg("n.min.frac")
  d0.min <- .arg("d0.min"); d1.min <- .arg("d1.min")
  effect.threshold <- .arg("effect.threshold"); hr.threshold <- .arg("hr.threshold")
  consistency.threshold <- .arg("consistency.threshold"); hr.consistency <- .arg("hr.consistency")
  pconsistency.threshold <- .arg("pconsistency.threshold")
  stop_threshold <- if ("stop_threshold" %in% names(args)) args[["stop_threshold"]]
                    else pconsistency.threshold            # formal default
  use_twostage <- .arg("use_twostage"); twostage_args <- .arg("twostage_args")
  fs.splits <- .arg("fs.splits")
  consistency_method <- .arg("consistency_method")
  if (length(consistency_method) > 1L) consistency_method <- consistency_method[1L]
  max_subgroups_search <- .arg("max_subgroups_search")
  m1.threshold <- .arg("m1.threshold")
  selection_rule <- .arg("selection_rule"); effect_neighborhood <- .arg("effect_neighborhood")
  max.minutes <- .arg("max.minutes")
  cut_args <- c("conf.cont_jcuts", "cut_type", "cont.cutoff", "conf.cont_medians",
                "conf.cont_medians_force", "conf_force", "defaultcut_names",
                "exclude_cuts", "collapse_cuts", "collapse_cuts_args")
  cut_vals <- lapply(cut_args, .arg); names(cut_vals) <- cut_args

  # maxeff overrides (forestsearch_main.R L1476-1490; L2919, L2925)
  effect_floor_disabled <- FALSE
  if (is_maxeff) {
    pconsistency.threshold <- 0; stop_threshold <- NULL; use_twostage <- FALSE
    max_subgroups_search <- Inf; minp <- 0; rmin <- 0
    effect_floor_disabled <- TRUE
  }
  # stop_threshold is meaningful for maxeffCons only (forestsearch_main.R
  # L1554-1574; reset at L1575-1608 for the hr / hrMaxSG / hrMinSG / maxSG /
  # minSG family, and by the maxeff block above)
  if (!is.null(stop_threshold) &&
      as.character(sg_focus) %in% c("hrMaxSG", "hrMinSG", "hr", "maxSG", "minSG")) {
    stop_threshold <- NULL
  }
  # n.min: SECTION 1A2 -- supplied value, or max(60, ceiling(n.min.frac * N))
  n.min_note <- NULL
  if (is.null(n.min)) {
    if (!is.null(data)) {
      n.min <- max(60L, as.integer(ceiling(n.min.frac * nrow(data))))
      n.min_note <- sprintf("n.min = NULL resolved to max(60, ceiling(%s * %d)) = %d",
                            format(n.min.frac), nrow(data), n.min)
    } else {
      n.min_note <- sprintf("n.min = NULL: resolved at run time to max(60, ceiling(%s * N))",
                            format(n.min.frac))
    }
  }
  glm_like <- outcome_type %in% c("continuous", "count")

  # ---- 3. optional grounding: cut columns and combinations on `data` --------
  n_cut_columns <- NA_integer_; n_combinations <- NA_real_
  if (!is.null(data)) {
    if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
    conf <- .arg("confounders.name")
    if (is.null(conf)) stop("`confounders.name` is needed to count cut columns on `data`.", call. = FALSE)
    missing_conf <- setdiff(conf, names(data))
    if (length(missing_conf)) {
      stop("`data` lacks confounders: ", paste(missing_conf, collapse = ", "), ".", call. = FALSE)
    }
    # get_FSdata() validates an outcome and an event column; neither is used
    # with LASSO off, so two constant numeric columns satisfy the contract
    # (the same device fs_oc_family_enumerate() uses).
    df_cut <- data[, conf, drop = FALSE]
    df_cut[[".fs_fr_y"]] <- 0; df_cut[[".fs_fr_event"]] <- 1
    FSdata <- get_FSdata(
      df.analysis = df_cut, use_lasso = FALSE, use_grf = FALSE, grf_cuts = NULL,
      confounders.name = conf,
      cont.cutoff = cut_vals$cont.cutoff, conf_force = cut_vals$conf_force,
      conf.cont_medians = cut_vals$conf.cont_medians,
      conf.cont_medians_force = cut_vals$conf.cont_medians_force,
      conf.cont_jcuts = cut_vals$conf.cont_jcuts, dina_cuts = NULL,
      collapse_cuts = cut_vals$collapse_cuts, collapse_cuts_args = cut_vals$collapse_cuts_args,
      defaultcut_names = cut_vals$defaultcut_names, cut_type = cut_vals$cut_type,
      exclude_cuts = cut_vals$exclude_cuts,
      outcome.name = ".fs_fr_y", event.name = ".fs_fr_event",
      details = FALSE, outcome_type = outcome_type)
    Zdf <- dummy(FSdata$df[, FSdata$confs_names, drop = FALSE])
    n_cut_columns  <- ncol(Zdf)
    n_combinations <- calculate_max_combinations(n_cut_columns, as.integer(maxk))
  }

  # ---- 4. the stage table ----------------------------------------------------
  fmt <- function(v) {
    if (is.null(v)) return("NULL")
    if (is.list(v)) {
      if (!length(v)) return("list()")
      return(paste0("list(", paste(sprintf("%s = %s", names(v),
             vapply(v, function(e) paste(deparse(e), collapse = ""), "")), collapse = ", "), ")"))
    }
    paste(deparse(v, width.cutoff = 500L), collapse = "")
  }
  DD  <- "data-dependent"; DDN <- "data-dependent (not disableable)"
  DET <- "deterministic";  DIS <- "disabled"; INE <- "inert"
  rows <- list()
  add <- function(stage, arguments, values, status, note = NA_character_) {
    rows[[length(rows) + 1L]] <<- data.frame(
      stage = stage, arguments = paste(arguments, collapse = ", "),
      values = values, status = status, note = note, stringsAsFactors = FALSE)
  }

  add("cut construction", cut_args,
      paste(sprintf("%s = %s", cut_args, vapply(cut_vals, fmt, "")), collapse = "; "),
      DDN,
      paste0("continuous cuts are placed at SAMPLE quantiles (get_FSdata(): get_conf_force() / ",
             "cut_var_jq(), medians for cut_type = \"median\"); only a full set of cut VALUES ",
             "supplied through conf_force with the defaults excluded would fix them, and the ",
             "floors and near-duplicate removal below still vary"))
  add("LASSO screen", "use_lasso", sprintf("use_lasso = %s", fmt(use_lasso)),
      if (use_lasso) DD else DIS,
      if (use_lasso) "Cox/GLM-LASSO selection on the sample decides which covariates are cut (get_FSdata() L362-395); use_lasso = FALSE disables it"
      else "use_lasso = FALSE")
  add("DINA / GRF subgroup paths", c("use_dina", "use_grf", "subgroup_method"),
      sprintf("use_dina = %s; use_grf = %s; subgroup_method = %s", fmt(use_dina), fmt(use_grf), fmt(subgroup_method)),
      if (use_dina || use_grf || !identical(subgroup_method, "consistency")) DD else DIS,
      if (use_dina || use_grf || !identical(subgroup_method, "consistency"))
        "model-identified cuts (GRF section 3A) / DINA candidates are fitted on the sample and enter the pool; use_grf = use_dina = FALSE with subgroup_method = \"consistency\" disables them"
      else "use_grf = use_dina = FALSE, subgroup_method = \"consistency\"")
  vi_on <- !is.null(vi.grf.min)
  add("GRF variable-importance ordering", "vi.grf.min", sprintf("vi.grf.min = %s", fmt(vi.grf.min)),
      if (vi_on) DD else DIS,
      if (!vi_on) "vi.grf.min = NULL skips Section 5 entirely (forestsearch_main.R L2752)"
      else if (is.numeric(vi.grf.min) && vi.grf.min <= 0)
        "a causal forest is fitted per call and ORDERS the cut columns by importance; at values <= 0 nothing is filtered (variable_importance() is non-negative, block guarded by vi_max > 0), but the order reaches rmin (column-walk) and the clause order of the winner; vi.grf.min = NULL disables it"
      else "a causal forest is fitted per call; columns with vi_ratio <= vi.grf.min are dropped and the rest reordered; vi.grf.min = NULL disables it")
  cap_status <- if (!vi_on) INE
                else if (!is.na(n_cut_columns) && is.finite(max_n_confounders) && max_n_confounders >= n_cut_columns) INE
                else DD
  add("confounder cap", "max_n_confounders", sprintf("max_n_confounders = %s", fmt(max_n_confounders)),
      cap_status,
      if (!vi_on) "applied only inside the variable-importance block (forestsearch_main.R L2820, under vi_max > 0, under !is.null(vi.grf.min)); inert while vi.grf.min = NULL"
      else if (cap_status == INE) sprintf("cap %s >= %d cut columns on this data: never binds", fmt(max_n_confounders), n_cut_columns)
      else "truncates the importance-ORDERED column list, so which columns survive depends on the fitted forest")
  add("dummy expansion", character(0), "each 2-level factor -> two indicator columns (both directions)", DET,
      "deterministic given the cut columns")
  add("combination enumeration", "maxk", sprintf("maxk = %s%s", fmt(maxk),
      if (!is.na(n_cut_columns)) sprintf(" (%d cut columns -> %s combinations)", n_cut_columns, format(n_combinations, big.mark = ",")) else ""),
      DET, "all combinations of up to maxk indicator columns; deterministic given the columns")
  add("per-factor prevalence floor", "minp", sprintf("minp = %s%s", fmt(minp), if (is_maxeff) " (maxeff: -> 0)" else ""),
      if (isTRUE(minp > 0)) DD else DIS,
      if (isTRUE(minp > 0)) "all(colMeans(x) >= minp) on the SAMPLE (subgroup.search(): meets_prevalence_threshold); minp = 0 or sg_focus = \"maxeff\" disables it"
      else "minp = 0: never binds")
  add("redundancy (rmin)", "sg_focus", sprintf("rmin = %s (subgroup.search() default%s)", fmt(rmin), if (is_maxeff) "; maxeff -> 0" else ""),
      DDN,
      "each added factor must shrink the SAMPLE membership by more than rmin subjects, walking the columns in their current order (extract_idx_flagredundancy()); rmin = 0 (maxeff) still drops exact-membership prefixes, so it cannot be switched off, and rmin is not a forestsearch() formal")
  add("subgroup size floor", c("n.min", "n.min.frac"),
      sprintf("n.min = %s; n.min.frac = %s%s", fmt(n.min), fmt(n.min.frac), if (!is.null(n.min_note)) paste0(" (", n.min_note, ")") else ""),
      DDN,
      "sum(id.x) > n.min on the SAMPLE (subgroup.search() status 4); a floor on a sample count with no off switch")
  add("per-arm floors", c("d0.min", "d1.min"), sprintf("d0.min = %s; d1.min = %s", fmt(d0.min), fmt(d1.min)),
      if (glm_like) INE else DD,
      if (glm_like) sprintf("skipped entirely for outcome_type = \"%s\" (subgroup.search() L593-609: continuous and count rely on n.min)", outcome_type)
      else if (outcome_type == "binary") "minimum EVENTS (Y = 1) per arm within the subgroup, on the sample"
      else "minimum events per arm within the subgroup, on the sample")
  eff_val <- if (outcome_type == "survival") sprintf("hr.threshold = %s", fmt(hr.threshold))
             else sprintf("effect.threshold = %s", fmt(effect.threshold))
  add("effect screen", c("effect.threshold", "hr.threshold"), paste0(eff_val, if (effect_floor_disabled) " (maxeff: floor disabled)" else ""),
      if (effect_floor_disabled) DIS else DD,
      if (effect_floor_disabled) "sg_focus = \"maxeff\" disables the search-stage effect floor (disable_effect_floor)"
      else "per-candidate fit on the sample must clear the threshold to enter the consistency stage")
  add("near-duplicate removal", "sg_focus", sprintf("sg_focus = %s -> %s", fmt(sg_focus_raw),
      if (is_maxeff) "exact-membership dedup" else "statistics-keyed dedup"),
      DDN,
      if (is_maxeff) "identical SAMPLE memberships collapse to the fewest-cut representative (.maxeff_membership_dedup()); no argument disables it"
      else "candidates whose (K, n, E, d1, m1, m0, HR, L, U) agree to 0.001 collapse to the first (remove_near_duplicate_subgroups(), subgroup_consistency_main.R L579); statistics are sample-fitted; NO argument disables it")
  add("candidate cap", "max_subgroups_search", sprintf("max_subgroups_search = %s%s", fmt(max_subgroups_search), if (is_maxeff) " (maxeff: -> Inf)" else ""),
      if (is.finite(max_subgroups_search)) DD else INE,
      if (is.finite(max_subgroups_search)) "truncates the preview-ORDERED pool to the top max_subgroups_search before consistency (subgroup_consistency_main.R L605-625); the order is sample-fitted, so what is evaluated depends on the data"
      else "Inf: never binds")
  add("m1 filter", "m1.threshold", sprintf("m1.threshold = %s", fmt(m1.threshold)),
      if (is.finite(m1.threshold)) DD else INE,
      if (is.finite(m1.threshold)) "candidates with treatment-arm count m1 above the threshold are dropped before consistency (subgroup_consistency_main.R L516-531)"
      else "Inf: never binds")
  cons_args <- c("consistency_method", "pconsistency.threshold", "consistency.threshold", "hr.consistency", "fs.splits", "use_twostage", "twostage_args")
  add("consistency screen", cons_args,
      sprintf("consistency_method = %s; pconsistency.threshold = %s%s; %s; fs.splits = %s; use_twostage = %s; twostage_args = %s",
              fmt(consistency_method), fmt(pconsistency.threshold), if (is_maxeff) " (maxeff: -> 0)" else "",
              if (outcome_type == "survival") sprintf("hr.consistency = %s", fmt(hr.consistency)) else sprintf("consistency.threshold = %s", fmt(consistency.threshold)),
              fmt(fs.splits), fmt(use_twostage), fmt(twostage_args)),
      if (is_maxeff) DIS else DD,
      if (is_maxeff) "sg_focus = \"maxeff\" sets pconsistency.threshold to 0: every candidate passes"
      else "the gate: per-candidate consistency rate on resampled / split halves of the sample")
  add("early stopping", "stop_threshold", sprintf("stop_threshold = %s%s", fmt(stop_threshold),
      if (is.null(stop_threshold) && !is.null(.arg("stop_threshold")) && !identical(as.character(sg_focus), "maxeffCons")) sprintf(" (reset to NULL for sg_focus = %s)", fmt(sg_focus_raw)) else ""),
      if (is.null(stop_threshold)) DIS else DD,
      if (is.null(stop_threshold)) "NULL: every qualifying candidate is evaluated"
      else "meaningful for sg_focus = \"maxeffCons\" only: truncates what is evaluated once the running Pcons reaches stop_threshold, but the prefix winner is the global winner (forestsearch_main.R L1554-1574), so the answer is unchanged")
  add("winner selection", c("sg_focus", "selection_rule", "effect_neighborhood"),
      sprintf("sg_focus = %s (canonical %s); selection_rule = %s; effect_neighborhood = %s", fmt(sg_focus_raw), fmt(sg_focus), fmt(selection_rule), fmt(effect_neighborhood)),
      DET, "deterministic given the qualifying set")
  add("time cap", "max.minutes", sprintf("max.minutes = %s", fmt(max.minutes)), INE,
      "no path consults it (subgroup.search() forwards it and never compares it)")

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  counts <- table(factor(out$status, levels = c(DET, DIS, INE, DD, DDN)))
  counts <- stats::setNames(as.integer(counts), names(counts))
  n_dd <- counts[[DD]] + counts[[DDN]]
  verdict <- if (n_dd == 0L) {
    sprintf("Candidate family is deterministic given the data columns: 0 of %d stages vary with the sample.", nrow(out))
  } else {
    sprintf("Candidate family is data-dependent: %d of %d stages vary with the sample; %d %s disableable, %d %s intrinsic to the method.",
            n_dd, nrow(out), counts[[DD]], if (counts[[DD]] == 1L) "is" else "are",
            counts[[DDN]], if (counts[[DDN]] == 1L) "is" else "are")
  }
  attr(out, "verdict") <- verdict
  attr(out, "status_counts") <- counts
  attr(out, "data_supplied") <- !is.null(data)
  attr(out, "outcome_type") <- outcome_type
  if (!is.null(data)) {
    attr(out, "n_cut_columns") <- n_cut_columns
    attr(out, "n_combinations") <- n_combinations
  }
  class(out) <- c("fs_family_report", "data.frame")
  out
}

#' @rdname fs_family_report
#' @param ... Ignored.
#' @export
print.fs_family_report <- function(x, ...) {
  cat(attr(x, "verdict"), "\n", sep = "")
  if (isTRUE(attr(x, "data_supplied"))) {
    cat(sprintf("On the supplied data: %d cut columns, %s combinations of up to maxk.\n",
                attr(x, "n_cut_columns"), format(attr(x, "n_combinations"), big.mark = ",")))
  }
  cat("\n")
  tab <- data.frame(stage = x$stage, status = x$status,
                    values = ifelse(nchar(x$values) > 60L, paste0(substr(x$values, 1L, 57L), "..."), x$values),
                    stringsAsFactors = FALSE)
  print(tab, row.names = FALSE, right = FALSE)
  intrinsic <- x[x$status == "data-dependent (not disableable)", , drop = FALSE]
  if (nrow(intrinsic)) {
    cat("\nIntrinsic to the method -- no argument switches these off:\n")
    for (i in seq_len(nrow(intrinsic))) {
      cat(sprintf("  * %s: %s\n", intrinsic$stage[i], intrinsic$note[i]))
    }
  }
  invisible(x)
}
