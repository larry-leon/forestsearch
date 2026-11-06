# =============================================================================
# MRCT trial analysis: the shared core and the print / plot workhorses
# =============================================================================
#
# One analysis, two public entries: mrct_case_study() redraws a simulated trial,
# mrct_trial_analysis() takes an observed one, and both call .mrct_analyze_trial().
# mrct_four_sets() and the print / plot methods of all three classes share the
# workhorses below. The code was moved from R/mrct_case_study.R
# (TASK_mrct_trial_analysis_2026-11-06).
#
# =============================================================================


#' The analysis of one trial: the identifier, H-hat, the sets, the trigger
#'
#' Runs forestsearch() on the non-region patients with the region as the test
#' set, applies the declared rule to the whole trial, and fits the sets.
#'
#' @param dfs The trial (data.table).
#' @param region_var Region indicator column.
#' @param truth_var The column of H* (1 in H*), or `NULL` for none.
#' @param fs_spec Named list: the identifier's arguments.
#' @param mr_inference,mr_args Multiplier resampling, and its arguments.
#' @param nm Named list of column names: `outcome`, `treat`, `event`, `id`, `po`
#'   (the last two `NULL` when the trial has none).
#' @param caller The public entry's name, for its messages.
#' @return List: `trial`, `fs`, `sets`, `in_Hhat`, `confusion`, `Hhat`, `mr`,
#'   `trigger`, `warnings`, `fs_call_args`.
#' @noRd
.mrct_analyze_trial <- function(dfs, region_var, truth_var, fs_spec, mr_inference, mr_args,
                                nm, caller) {
  df_non <- dfs[dfs[[region_var]] == 0]
  df_reg <- dfs[dfs[[region_var]] == 1]

  fs_base <- list(
    df.analysis = as.data.frame(df_non), df.test = as.data.frame(df_reg),
    outcome.name = nm$outcome, treat.name = nm$treat, event.name = nm$event,
    id.name = nm$id,
    potentialOutcome.name = nm$po,
    show_candidate_summary = FALSE, details = FALSE, plot.sg = FALSE,
    parallel_args = list(plan = "sequential"), quiet = TRUE,
    mr_inference = isTRUE(mr_inference),
    mr_inference_args = if (isTRUE(mr_inference)) mr_args else list())
  fs_call_args <- utils::modifyList(fs_spec, fs_base)
  fs_warn <- character(0)
  fs <- withCallingHandlers(
    do.call(forestsearch, fs_call_args),
    warning = function(w) {
      if (!grepl("parallel_args missing required elements", conditionMessage(w)))
        fs_warn <<- c(fs_warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (length(fs_warn))
    message(caller, "(): forestsearch() raised ", length(fs_warn),
            " warning(s), kept in $warnings.")

  trial <- as.data.frame(dfs)
  declared <- !is.null(fs$sg.harm)
  in_Hhat <- if (declared) {
    get_dfpred(trial, sg.harm = fs$sg.harm, version = 2)$treat.recommend == 0
  } else {
    rep(FALSE, nrow(trial))
  }
  built <- .mrct_four_set_table(trial, region_var, truth_var, in_Hhat, nm)

  list(trial = trial, fs = fs, sets = built$sets, in_Hhat = in_Hhat,
       confusion = built$confusion,
       Hhat = if (declared) paste(fs$sg.harm, collapse = " & ") else "none",
       mr = fs$mr_inference, trigger = .mrct_trigger(built$sets),
       warnings = fs_warn, fs_call_args = fs_call_args)
}


#' The set table and the 2 x 2 tables of H-hat against H*
#'
#' @param trial The trial (data.frame).
#' @param region_var Region indicator column.
#' @param truth_var The column of H* (1 in H*), or `NULL`: no region(G*) row and
#'   no 2 x 2 tables.
#' @param in_Hhat Logical per row, or `NULL` before the identifier has run: no
#'   region(G-hat) row and no 2 x 2 tables.
#' @param nm Named list of column names: `outcome`, `treat`, `event`.
#' @return List: `sets`, `confusion`.
#' @noRd
.mrct_four_set_table <- function(trial, region_var, truth_var, in_Hhat, nm) {
  reg <- trial[[region_var]] == 1
  has_band <- !is.null(truth_var)
  band <- if (has_band) trial[[truth_var]] == 1 else NULL
  set_idx <- list(itt_all = rep(TRUE, nrow(trial)), itt_region = reg)
  if (has_band) set_idx$region_Gstar <- reg & !band
  if (!is.null(in_Hhat)) set_idx$region_Ghat <- reg & !in_Hhat
  labels <- c(itt_all = "All patients, ITT", itt_region = "Region, ITT",
              region_Gstar = "Region minus H* (G*)", region_Ghat = "Region minus H-hat (G-hat)")
  # .mrct_record_fit() reads these three column names
  dt <- data.table::data.table(y_sim = trial[[nm$outcome]], event_sim = trial[[nm$event]],
                               treat_sim = trial[[nm$treat]])
  sets <- do.call(rbind, lapply(names(set_idx), function(s) {
    f <- as.data.frame(.mrct_record_fit(dt, set_idx[[s]]))
    data.frame(set = s, label = unname(labels[s]), n = f$n, n_exp = f$n_exp,
               n_ctl = f$n - f$n_exp, events_exp = f$events_exp, events_ctl = f$events_ctl,
               hr = f$hr, lo = f$lo, hi = f$hi, se_loghr = f$se_loghr,
               stringsAsFactors = FALSE)
  }))
  rownames(sets) <- NULL

  confusion <- NULL
  if (has_band && !is.null(in_Hhat)) {
    cc <- function(w) c(TP = sum(in_Hhat[w] & band[w]), FP = sum(in_Hhat[w] & !band[w]),
                        FN = sum(!in_Hhat[w] & band[w]), TN = sum(!in_Hhat[w] & !band[w]))
    confusion <- list(nonregion = cc(!reg), region = cc(reg))
  }
  list(sets = sets, confusion = confusion)
}


#' The trigger: all-patient upper limit < 1.0 and region hazard ratio > 1.0
#' @noRd
.mrct_trigger <- function(sets) {
  itt <- sets[sets$set == "itt_all", ]; rgn <- sets[sets$set == "itt_region", ]
  !is.na(itt$hi) && itt$hi < 1 && !is.na(rgn$hr) && rgn$hr > 1
}


#' Rows of the trial in one of the sets
#' @noRd
.mrct_set_rows <- function(trial, set, region_var, truth_var, in_Hhat) {
  reg <- trial[[region_var]] == 1
  switch(set,
         itt_all = rep(TRUE, nrow(trial)),
         itt_region = reg,
         region_Gstar = reg & trial[[truth_var]] != 1,
         region_Ghat = reg & !in_Hhat)
}


#' The number formatter of the print methods
#' @noRd
.mrct_fm <- function(digits) {
  function(v) ifelse(is.na(v), "NA", formatC(v, format = "f", digits = digits))
}


#' Print the set table
#' @noRd
.mrct_print_sets <- function(s, fm) {
  tab <- data.frame(Set = s$label, n = s$n,
                    `exp n/events` = sprintf("%d/%d", s$n_exp, s$events_exp),
                    `ctl n/events` = sprintf("%d/%d", s$n_ctl, s$events_ctl),
                    HR = fm(s$hr), `95% CI` = sprintf("%s, %s", fm(s$lo), fm(s$hi)),
                    check.names = FALSE)
  print(tab, row.names = FALSE, right = FALSE)
}


#' Print the 2 x 2 tables of H-hat against H*, when there are any
#' @noRd
.mrct_print_confusion <- function(confusion) {
  if (!is.null(confusion)) {
    cat("\nH-hat against H* (TP, FP, FN, TN):\n")
    for (p in names(confusion))
      cat(sprintf("  %-10s %s\n", p, paste(confusion[[p]], collapse = ", ")))
  }
}


#' Print the multiplier-resampling block
#'
#' `evaluated = FALSE` is mrct_case_study()'s block, unchanged. `evaluated = TRUE`
#' reads mr_values(): the naive estimate, the two-term IJ interval, the field and
#' field-s one-sided bounds and the Bonferroni pairs, each bound shown with the side
#' it limits.
#' @noRd
.mrct_print_mr <- function(mr, fm, evaluated = FALSE) {
  if (is.null(mr) || is.null(mr$debiased)) {
    cat("\nMultiplier resampling: not run, or nothing declared.\n")
    return(invisible(NULL))
  }
  cat(sprintf("\nMultiplier resampling (%s, %s draws):\n", mr$ci_method,
              if (!is.null(mr$settings$draws)) mr$settings$draws else "?"))
  if (evaluated) {
    v <- mr_values(mr)
    is_H <- v$row == "H"
    side <- function(op, x) ifelse(is.na(x), "-", paste(op, fm(x)))
    rows <- data.frame(
      set = ifelse(is_H, "H-hat", "H-hat complement (non-region)"),
      naive = fm(v$naive_est),
      debiased = fm(v$debiased_est),
      `IJ 95% CI` = ifelse(is.na(v$ij_lower), "-", sprintf("%s, %s", fm(v$ij_lower), fm(v$ij_upper))),
      `field 1-sided 95%` = ifelse(is_H, side(">=", v$field_lower_1s), side("<=", v$field_upper_1s)),
      `field-s 1-sided 95%` = side("<=", v$field_upper_1s_s),
      Bonferroni = ifelse(is_H, side(">=", v$bonf_lower_H), side("<=", v$bonf_upper_Hc)),
      `Bonferroni, field-s` = ifelse(is_H, side(">=", v$bonf_lower_H_s), side("<=", v$bonf_upper_Hc_s)),
      check.names = FALSE)
    print(rows, row.names = FALSE, right = FALSE)
    return(invisible(NULL))
  }
  g <- function(...) { v <- mr; for (k in c(...)) { if (!is.list(v) || is.null(v[[k]])) return(NA_real_); v <- v[[k]] }; v }
  rows <- data.frame(
    set = c("H-hat", "H-hat complement (non-region)"),
    naive = fm(c(g("naive", "est"), g("complement", "naive", "est"))),
    debiased = fm(c(g("debiased", "est"), g("complement", "debiased", "est"))),
    `IJ 95% CI` = sprintf("%s, %s", fm(c(g("debiased", "lower"), g("complement", "debiased", "lower"))),
                          fm(c(g("debiased", "upper"), g("complement", "debiased", "upper")))),
    `field 1-sided 95%` = c(paste(">=", fm(g("field", "lower_1s"))),
                            paste("<=", fm(g("field", "complement", "upper_1s")))),
    Bonferroni = c(paste(">=", fm(g("field", "joint", "bonf_lower_H"))),
                   paste("<=", fm(g("field", "joint", "bonf_upper_Hc")))),
    check.names = FALSE)
  print(rows, row.names = FALSE, right = FALSE)
  invisible(NULL)
}


#' The forest plot and the Kaplan-Meier panels of the sets
#'
#' @param s The set table.
#' @param trial The trial (data.frame).
#' @param region_var,truth_var,in_Hhat As in `.mrct_set_rows()`.
#' @param which Any of `"forest"`, `"km"`, resolved by the method.
#' @param xlim,time_unit As in the plot methods.
#' @param analysis_time The Kaplan-Meier axis's upper limit.
#' @param nm Named list of column names: `outcome`, `treat`, `event`.
#' @param caller The method's name, for its message.
#' @noRd
.mrct_plot_sets <- function(s, trial, region_var, truth_var, in_Hhat, which, xlim, time_unit,
                            analysis_time, nm, caller) {
  if (is.null(time_unit)) time_unit <- "time"
  if ("forest" %in% which) {
    ord <- intersect(c("itt_all", "itt_region", "region_Ghat", "region_Gstar"), s$set)
    lab <- c(itt_all = "All patients, ITT", itt_region = "Region, ITT",
             region_Ghat = "Region minus the identified subgroup",
             region_Gstar = "Region minus the true harmed subgroup")
    cols <- c(itt_all = "black", itt_region = "firebrick", region_Ghat = "steelblue",
              region_Gstar = "grey40")[ord]
    f <- s[match(ord, s$set), ]; y <- rev(seq_along(ord))
    op <- graphics::par(mar = c(4, 17, 1, 9)); on.exit(graphics::par(op), add = TRUE)
    graphics::plot(pmin(pmax(f$hr, xlim[1]), xlim[2]), y, log = "x", xlim = xlim,
                   ylim = c(0.5, length(ord) + 0.5), pch = 15, cex = 1.3, yaxt = "n", xaxt = "n",
                   xlab = "hazard ratio, experimental vs control (log scale, 95% CI)", ylab = "",
                   col = cols)
    graphics::segments(pmax(f$lo, xlim[1]), y, pmin(f$hi, xlim[2]), y, lwd = 2, col = cols)
    graphics::abline(v = 1, lty = 2)
    at <- c(0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5)
    at <- at[at >= xlim[1] & at <= xlim[2]]
    graphics::axis(1, at = at, labels = format(at))
    graphics::axis(2, at = y, labels = lab[ord], las = 1, tick = FALSE)
    graphics::mtext(sprintf("%.3f (%.3f, %.3f)", f$hr, f$lo, f$hi), side = 4, at = y,
                    las = 1, line = 0.5)
    graphics::par(op)
  }
  if ("km" %in% which) {
    if (!requireNamespace("weightedsurv", quietly = TRUE)) {
      message(caller, "(): weightedsurv is not installed; the Kaplan-Meier ",
              "figure is skipped.")
      return(invisible(NULL))
    }
    ta <- analysis_time
    br <- ta / 4
    op <- graphics::par(mfrow = c(2, 2), mar = c(4.5, 4.2, 3, 1))
    on.exit(graphics::par(op), add = TRUE)
    for (i in seq_len(nrow(s))) {
      d <- trial[.mrct_set_rows(trial, s$set[i], region_var, truth_var, in_Hhat), ]
      dfc <- suppressWarnings(weightedsurv::df_counting(
        df = d, by.risk = br, tte.name = nm$outcome, event.name = nm$event,
        treat.name = nm$treat, arms = c("experimental", "control")))
      weightedsurv::plot_weighted_km(
        dfc, conf.int = FALSE, show.cox = FALSE, show.logrank = FALSE, show.med = FALSE,
        show_arm_legend = TRUE, arms = c("experimental", "control"), put.legend.arms = "left",
        by.risk = br, xmin = 0, xmax = ta, ymin = 0, ymax = 1.12, Xlab = time_unit,
        Ylab = "proportion event-free", col.0 = "black", col.1 = "blue", lwds = c(1.5, 1.5),
        censor.mark.all = FALSE, risk.cex = 0.8, arm.cex = 0.85)
      graphics::axis(1, at = seq(0, ta, br))
      graphics::text(ta, 1.12, adj = c(1, 1), cex = 0.85, labels = paste(
        sprintf("HR %.3f (95%% CI %.3f, %.3f)", s$hr[i], s$lo[i], s$hi[i]),
        sprintf("experimental n = %d, events %d", s$n_exp[i], s$events_exp[i]),
        sprintf("control n = %d, events %d", s$n_ctl[i], s$events_ctl[i]), sep = "\n"))
      graphics::title(sprintf("%d. %s", i, s$label[i]), cex.main = 1.05)
    }
  }
  invisible(NULL)
}
