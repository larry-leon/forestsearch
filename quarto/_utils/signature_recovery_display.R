# =============================================================================
# signature_recovery_display.R --- reshape signature-recovery agg -> OC display
# =============================================================================
# Location: quarto/_utils/signature_recovery_display.R
# Dependencies: data.table (load in the consuming .qmd's setup chunk).
#
# Maps the per-method aggregate table produced by the subgroup-method
# signature-recovery vignettes (subgroup_method_signature_recovery.qmd and its
# _2factor companion) onto the operating-characteristic display convention used
# by plot_oc_metrics() / compute_oc_pareto_flag() in this _utils directory, so
# the same bubble-chart and Pareto machinery applies without modification.
#
# The signature-recovery agg has columns:
#   label, n, tau_in, method,
#   detection_rate, covariate_set_hit, signature_correct,
#   sensitivity, specificity, ppv, npv, mean_sg_frac, mean_n_factors
#
# plot_oc_metrics() expects numeric metric columns (defaults Spec / NPV /
# Detection) plus categorical color (group) and shape (Analysis) columns.
# This helper renames to that convention and adds:
#   group     - the method label (color aesthetic)
#   Analysis  - method family: "FS-consistency", "FS-grf", or "DINA" (shape)
#   Detection - detection_rate;  Sen - sensitivity;  Spec - specificity;
#   PPV       - ppv;  NPV - npv
# and carries the recovery-specific metrics under tidy names:
#   CovHit (covariate_set_hit), SigCorrect (signature_correct),
#   Shat_frac (mean_sg_frac), NFactors (mean_n_factors).
#
# This file is NOT part of the forestsearch package; it is an exploratory
# utility shared across signature-recovery evaluation documents.
# =============================================================================


# .sigrec_method_family() ---------------------------------------------------
# Map a method label to its family, used as the plot shape aesthetic.
.sigrec_method_family <- function(method) {
  m <- as.character(method)
  out <- rep("other", length(m))
  out[grepl("^FS-grf", m)]              <- "FS-grf"
  out[grepl("^FS-(eff|effMaxSG)$", m)]  <- "FS-consistency"
  out[grepl("^DINA", m)]                <- "DINA"
  out
}


# sigrec_to_display() -------------------------------------------------------
#
# Arguments
#   agg   data.table/data.frame of per-method aggregate rows from the
#         signature-recovery vignette (one row per label x n x tau_in x method).
#
# Returns: a data.table in the plot_oc_metrics() display convention, one row
#   per original agg row, with both the OC column names (Detection/Sen/Spec/PPV/NPV)
#   and the recovery-specific metrics (CovHit/SigCorrect/Shat_frac/
#   NFactors).  The scenario keys (label, n, tau_in) and the original `method`
#   are retained for faceting and labelling.
sigrec_to_display <- function(agg) {
  dt <- data.table::as.data.table(data.table::copy(agg))

  req <- c("method", "n", "tau_in", "detection_rate", "sensitivity",
           "specificity", "ppv", "npv", "covariate_set_hit",
           "signature_correct", "mean_sg_frac")
  miss <- setdiff(req, names(dt))
  if (length(miss)) {
    stop("sigrec_to_display(): agg is missing column(s): ",
         paste(miss, collapse = ", "), call. = FALSE)
  }

  out <- data.table::data.table(
    group      = factor(dt$method),                    # color aesthetic
    Analysis   = factor(.sigrec_method_family(dt$method),
                        levels = c("FS-consistency", "FS-grf", "DINA", "other")),
    method     = dt$method,
    n          = dt$n,
    tau_in     = dt$tau_in,
    label      = if ("label" %in% names(dt)) dt$label else NA_character_,
    # OC-convention numeric columns (higher-is-better, on [0, 1]):
    Detection  = dt$detection_rate,
    Sen        = dt$sensitivity,
    Spec       = dt$specificity,
    PPV        = dt$ppv,
    NPV        = dt$npv,
    # recovery-specific metrics carried alongside:
    CovHit     = dt$covariate_set_hit,
    SigCorrect = dt$signature_correct,
    Shat_frac  = dt$mean_sg_frac,
    NFactors   = if ("mean_n_factors" %in% names(dt)) dt$mean_n_factors else NA_real_
  )
  out[]
}


# sigrec_fdr_to_display() ---------------------------------------------------
# Reshape the null-arm false-discovery-rate table for plotting.  The FDR
# story is categorical (a couple of methods near 1, the rest near 0), so it
# is shown as a bar/lollipop panel rather than a bubble; this helper just
# adds the group/Analysis aesthetics and an `FDR` alias.
#
# Input `fdr`: columns n, method, false_discovery_rate (the vignette's
#   null-arm aggregate).
# Returns a data.table with the originals plus group, Analysis, FDR.
sigrec_fdr_to_display <- function(fdr) {
  dt <- data.table::as.data.table(fdr)
  req <- c("method", "n", "false_discovery_rate")
  miss <- setdiff(req, names(dt))
  if (length(miss))
    stop("sigrec_fdr_to_display(): `fdr` missing columns: ",
         paste(miss, collapse = ", "), call. = FALSE)
  out <- data.table::data.table(
    method   = dt$method,
    n        = dt$n,
    group    = factor(dt$method),
    Analysis = factor(.sigrec_method_family(dt$method),
                      levels = c("FS-consistency", "FS-grf", "DINA", "other")),
    FDR      = dt$false_discovery_rate
  )
  out[]
}
