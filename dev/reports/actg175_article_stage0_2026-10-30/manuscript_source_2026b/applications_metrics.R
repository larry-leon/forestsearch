#!/usr/bin/env Rscript
# applications_metrics.R -- produces the applications extract under dev/payloads/applications/
# from the six committed GBSG and ACTG 175 analysis payloads, and COLUMNS.md beside it.
#
# It mirrors dev/payloads/R/identifiers_metrics.R: a transformation of committed payloads and
# nothing else. No trial data is read, no analysis is re-run, no model is fitted. Every CSV
# carries, on every row, the file the row was read from (`source_payload`), that file's
# md5 (`source_md5`) and the time this extract was built (`extract_built_at`).
#
# Output, in dev/payloads/applications/:
#   applications_metrics.csv               ten rows, one per analysis x region (H, Hc), wide
#   applications_funnels.csv               long, one row per analysis x stage
#   gbsg_agreement_pairwise.csv            the trio summary's three pairwise rows
#   gbsg_agreement_partition.csv           the trio summary's eight partition rows
#   gbsg_memberships.csv                   the trio summary's three membership-count rows
#   gbsg_rules_overlap.csv                 max N and max beta against max N(eps), from the FS payload
#   actg175_binary_published_overlap.csv   the published GRF and DINA regions, from the binary payload
#   COLUMNS.md                             every column of every file, with its source and derivation
#
# Analysis keys are fixed by file name: gbsg_fs, gbsg_dina, gbsg_grf, actg175_binary,
# actg175_continuous.
#
# Gates, checked here and reported by the task:
#   G2  region_def on H equals the trio summary's memberships$definitions entry (the three GBSG
#       analyses); region_label appears verbatim in that analysis's REPORT; n on H equals the
#       report's N for the selection.
#   G3  n(H) + n(Hc) == n_total on every analysis, and events(H) + events(Hc) == events_total
#       wherever events exist.
#
# Usage:  Rscript dev/payloads/R/applications_metrics.R [out_dir]
#         (default dev/payloads/applications)
# Base R only. Run from the repository root.

args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1L]] else "dev/payloads/applications"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

BUILT_AT <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
md5 <- function(f) unname(tools::md5sum(f))
stop_if <- function(cond, ...) if (!isTRUE(cond)) stop(..., call. = FALSE)

# ---- the five analyses, fixed by file name ---------------------------------------------------
SPEC <- data.frame(
  analysis   = c("gbsg_fs", "gbsg_dina", "gbsg_grf", "actg175_binary", "actg175_continuous"),
  trial      = c("GBSG", "GBSG", "GBSG", "ACTG175", "ACTG175"),
  identifier = c("FS", "DINA", "GRF", "FS", "FS"),
  file = c("quarto/gbsg/analysis_gbsg_mr_payload.rds",
           "quarto/gbsg/analysis_gbsg_dina_mr_payload.rds",
           "quarto/gbsg/analysis_gbsg_grf_mr_payload.rds",
           "quarto/actg175/analysis_actg175_binary_mr_payload.rds",
           "quarto/actg175/analysis_actg175_continuous_mr_payload.rds"),
  report = c("quarto/gbsg/REPORT_gbsg_mr_2026-09-17.md",
             "quarto/gbsg/REPORT_gbsg_dina_grf_mr_2026-09-17.md",
             "quarto/gbsg/REPORT_gbsg_dina_grf_mr_2026-09-17.md",
             "quarto/actg175/REPORT_actg175_binary_mr_2026-09-17.md",
             "quarto/actg175/REPORT_actg175_continuous_mr_2026-09-17.md"),
  # the plain-language description each analysis's report prints for its selection; asserted
  # verbatim against that report by Gate G2, never invented here
  region_label = c("estrogen receptor 0, progesterone receptor ≤ 26 fmol/l",
                   "grade 3 and progesterone receptor ≤ 20",
                   "estrogen receptor 0",
                   "weight > 86 kg and baseline CD4 count > 380 per mm³",
                   "baseline CD4 count > 507 per mm³ and male"),
  # the label element carrying the selected region's coded definition
  label_field = c("fs", "dina", "grf", "fs", "fs"),
  # forest search applies the splitting-consistency screen; DINA and GRF do not
  screen_applies = c(TRUE, FALSE, FALSE, TRUE, TRUE),
  stringsAsFactors = FALSE)
TRIO <- "quarto/gbsg/analysis_gbsg_trio_summary_payload.rds"

P <- lapply(SPEC$file, readRDS); names(P) <- SPEC$analysis
TR <- readRDS(TRIO)

# ---- helpers ---------------------------------------------------------------------------------
gv <- function(l, nm) { v <- l[[nm]]; if (is.null(v) || length(v) != 1L) NA else v }
# a setting, taken from meta$settings and, failing that, from meta$thresholds_received
setg <- function(p, nm) {
  v <- gv(p$meta$settings, nm)
  if (is.na(v) && is.list(p$meta$thresholds_received)) v <- gv(p$meta$thresholds_received, nm)
  if (is.na(v) && is.numeric(p$meta$thresholds_received) && nm %in% names(p$meta$thresholds_received))
    v <- unname(p$meta$thresholds_received[[nm]])
  v
}
# the whole-trial reference row: reference_hr / reference_or / reference_md, whichever the payload carries
refdf <- function(p) {
  for (nm in c("reference_hr", "reference_or", "reference_md"))
    if (!is.null(p$extras[[nm]])) return(p$extras[[nm]])
  NULL
}
ref_whole <- function(p) { d <- refdf(p); d[d$set == "whole trial", , drop = FALSE] }
est_col <- function(d) intersect(c("hr", "or", "md"), names(d))[1L]
# provenance, added to every file
prov <- function(d, file) {
  d$source_payload   <- file
  d$source_md5       <- md5(file)
  d$extract_built_at <- BUILT_AT
  d
}
wr <- function(d, name) {
  f <- file.path(out_dir, name)
  write.csv(d, f, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  cat(sprintf("wrote %-38s %3d rows x %2d cols\n", name, nrow(d), ncol(d)))
  invisible(f)
}

# ---- 1. applications_metrics.csv --------------------------------------------------------------
rows <- list()
for (i in seq_len(nrow(SPEC))) {
  a <- SPEC$analysis[i]; p <- P[[a]]; tb <- p$table
  stop_if(identical(sort(tb$region), c("H", "Hc")), a, ": table$region is not H and Hc")
  H <- tb[tb$region == "H", ]; C <- tb[tb$region == "Hc", ]
  rw <- ref_whole(p); ec <- est_col(rw)
  scale <- p$est_scale
  is_ratio <- scale %in% c("hr", "or")

  # region definition, in the search's own coding
  def <- p$labels[[SPEC$label_field[i]]]
  stop_if(is.character(def) && !is.na(def) && nzchar(def), a, ": no coded region definition")

  # arm counts and events on H; the complement's are the whole trial's less H's
  hr_ <- p$extras$harm_rates; am <- p$extras$arm_means
  trt_n_H  <- if (!is.null(hr_)) hr_$treated_n else if (!is.null(am)) am$treated_n else NA
  ctrl_n_H <- if (!is.null(hr_)) hr_$control_n else if (!is.null(am)) am$control_n else NA
  trt_e_H  <- if (!is.null(hr_)) hr_$treated_events else NA
  ctrl_e_H <- if (!is.null(hr_)) hr_$control_events else NA
  n_total      <- p$meta$n_total
  trt_n_tot    <- rw$n_trt; ctrl_n_tot <- rw$n_ctrl
  has_events   <- all(c("ev_trt", "ev_ctrl") %in% names(rw))
  events_total <- if (has_events) rw$ev_trt + rw$ev_ctrl else NA
  events_H     <- if (!is.na(trt_e_H)) trt_e_H + ctrl_e_H else NA
  events_C     <- if (!is.na(events_H) && !is.na(events_total)) events_total - events_H else NA
  trt_n_C  <- trt_n_tot - trt_n_H; ctrl_n_C <- ctrl_n_tot - ctrl_n_H
  trt_e_C  <- if (has_events && !is.na(trt_e_H)) rw$ev_trt - trt_e_H else NA
  ctrl_e_C <- if (has_events && !is.na(ctrl_e_H)) rw$ev_ctrl - ctrl_e_H else NA
  pc <- function(e, n) if (is.na(e) || is.na(n) || n == 0) NA else 100 * e / n
  # arm means: H from arm_means, the complement derived from the whole trial's arm means
  mean_trt_H <- if (!is.null(am)) am$treated_mean_decline else NA
  mean_ctl_H <- if (!is.null(am)) am$control_mean_decline else NA
  dmean <- function(n_tot, m_tot, n_h, m_h) if (any(is.na(c(n_tot, m_tot, n_h, m_h)))) NA else
    (n_tot * m_tot - n_h * m_h) / (n_tot - n_h)
  mean_trt_C <- if (!is.null(am)) dmean(trt_n_tot, rw$mean_trt, trt_n_H, mean_trt_H) else NA
  mean_ctl_C <- if (!is.null(am)) dmean(ctrl_n_tot, rw$mean_ctrl, ctrl_n_H, mean_ctl_H) else NA

  b <- p$bounds
  fc <- function(nv, mv) if (is_ratio) log(nv) - log(mv) else nv - mv
  ord_H <- b$field_lower_H <= b$field_est2_H && b$field_est2_H <= H$mr_est
  ord_C <- C$mr_est <= b$field_s_est2_Hc && b$field_s_est2_Hc <= b$field_s_upper_Hc

  common <- list(
    analysis = a, trial = SPEC$trial[i], identifier = SPEC$identifier[i],
    outcome_type = setg(p, "outcome_type"), est_scale = scale,
    adverse_outcome = setg(p, "adverse_outcome"),
    endpoint = p$meta$endpoint, outcome = p$meta$outcome,
    n_total = n_total, events_total = events_total,
    whole_trial_est = rw[[ec]], whole_trial_lo = rw$lo, whole_trial_hi = rw$hi,
    bonf_joint_prob = b$bonf_joint_prob, calib_gamma = b$calib_gamma,
    calib_joint_prob = b$calib_joint_prob, calib_lower_H = b$calib_lower_H,
    calib_upper_Hc = b$calib_upper_Hc, joint_corr = b$joint_corr,
    rho_c = if (is.null(b$rho_c)) NA else b$rho_c, alpha = b$alpha,
    selection_bias = p$meta$timings$fs$selection_bias, n_family = p$meta$timings$fs$n_family,
    sg_focus = setg(p, "sg_focus"), eps = setg(p, "effect_neighborhood"),
    selection_rule = setg(p, "selection_rule"),
    c1 = { v <- setg(p, "effect.threshold"); if (is.na(v)) setg(p, "hr.threshold") else v },
    c2 = { v <- setg(p, "consistency.threshold"); if (is.na(v)) setg(p, "hr.consistency") else v },
    pstar = setg(p, "pconsistency.threshold"), screen_applies = SPEC$screen_applies[i],
    consistency_method = setg(p, "consistency_method"), n_min = setg(p, "n.min"),
    d_min = setg(p, "d0.min"), maxk = setg(p, "maxk"), mr_draws = p$meta$mr_draws,
    field_R_out = setg(p, "field_R_out"), field_R_in = setg(p, "field_R_in"),
    seed = setg(p, "seedit"), n_cores = p$meta$cores,
    source_payload = SPEC$file[i], source_md5 = md5(SPEC$file[i]),
    payload_built_at = format(p$built_at), forestsearch_version = p$forestsearch_version,
    extract_built_at = BUILT_AT)

  side <- function(reg) {
    isH <- reg == "H"; r <- if (isH) H else C
    c(common[c("analysis", "trial", "identifier", "outcome_type", "est_scale", "adverse_outcome",
               "endpoint", "outcome")],
      list(region = reg,
           region_def   = if (isH) def else paste("complement of", def),
           region_label = if (isH) SPEC$region_label[i] else paste("complement of", SPEC$region_label[i]),
           n = r$n, pct = r$pct,
           events   = if (isH) events_H else events_C,
           trt_n    = if (isH) trt_n_H else trt_n_C,
           ctrl_n   = if (isH) ctrl_n_H else ctrl_n_C,
           trt_events  = if (isH) trt_e_H else trt_e_C,
           ctrl_events = if (isH) ctrl_e_H else ctrl_e_C,
           trt_pct  = if (isH) (if (is.null(hr_)) NA else hr_$treated_pct) else pc(trt_e_C, trt_n_C),
           ctrl_pct = if (isH) (if (is.null(hr_)) NA else hr_$control_pct) else pc(ctrl_e_C, ctrl_n_C),
           trt_mean  = if (isH) mean_trt_H else mean_trt_C,
           ctrl_mean = if (isH) mean_ctl_H else mean_ctl_C),
      common[c("n_total", "events_total", "whole_trial_est", "whole_trial_lo", "whole_trial_hi")],
      list(naive_est = r$naive_est, naive_lo = r$naive_lo, naive_hi = r$naive_hi,
           mr_est = r$mr_est, mr_lo = r$mr_lo, mr_hi = r$mr_hi,
           est2        = if (isH) b$field_est2_H else b$field_s_est2_Hc,
           lambda_mean = if (isH) b$field_lambda_mean_H else b$field_s_lambda_mean_Hc,
           est2_unstudentized        = if (isH) NA else b$field_est2_Hc,
           lambda_mean_unstudentized = if (isH) NA else b$field_lambda_mean_Hc,
           first_correction = fc(r$naive_est, r$mr_est),
           bound_side = if (isH) "lower" else "upper",
           field_bound = if (isH) b$field_lower_H else b$field_s_upper_Hc,
           bonf_bound  = if (isH) b$bonf_lower_H else b$bonf_upper_Hc),
      common[c("bonf_joint_prob", "calib_gamma", "calib_joint_prob", "calib_lower_H",
               "calib_upper_Hc", "joint_corr", "rho_c", "alpha")],
      list(n_out = if (isH) b$n_out_H else b$n_out_Hc,
           n_in_mean = if (isH) b$n_in_H else NA,
           ordering_holds = if (isH) ord_H else ord_C),
      common[c("selection_bias", "n_family", "sg_focus", "eps", "selection_rule", "c1", "c2",
               "pstar", "screen_applies", "consistency_method", "n_min", "d_min", "maxk",
               "mr_draws", "field_R_out", "field_R_in", "seed", "n_cores",
               "source_payload", "source_md5", "payload_built_at", "forestsearch_version",
               "extract_built_at")])
  }
  rows[[length(rows) + 1L]] <- as.data.frame(side("H"), stringsAsFactors = FALSE)
  rows[[length(rows) + 1L]] <- as.data.frame(side("Hc"), stringsAsFactors = FALSE)
}
MET <- do.call(rbind, rows)

# ---- Gate G2: definitions, labels and N ------------------------------------------------------
cat("\nGate G2 -- region definitions, labels and N\n")
defs <- TR$memberships$definitions
for (i in seq_len(nrow(SPEC))) {
  a <- SPEC$analysis[i]; h <- MET[MET$analysis == a & MET$region == "H", ]
  if (SPEC$trial[i] == "GBSG") {
    key <- SPEC$identifier[i]
    stop_if(identical(unname(defs[[key]]), h$region_def),
            "G2 FAIL: ", a, ": region_def '", h$region_def, "' is not the trio summary's '", defs[[key]], "'")
    cat(sprintf("  %-19s def == trio summary's %-4s : %s\n", a, key, h$region_def))
  }
  rpt <- readLines(SPEC$report[i], warn = FALSE)
  k <- grep(h$region_label, rpt, fixed = TRUE)
  stop_if(length(k) >= 1L, "G2 FAIL: ", a, ": region_label '", h$region_label,
          "' does not appear verbatim in ", SPEC$report[i])
  # N for the selection, as the report prints it
  npat <- sprintf("(N|n) %s\\b", format(h$n, big.mark = ""))
  kn <- grep(npat, rpt)
  stop_if(length(kn) >= 1L, "G2 FAIL: ", a, ": the report does not print N ", h$n, " for the selection")
  cat(sprintf("  %-19s label verbatim in %s line(s) %s; N %d printed on line(s) %s\n",
              a, basename(SPEC$report[i]), paste(head(k, 3), collapse = ","), h$n,
              paste(head(kn, 3), collapse = ",")))
}

# ---- Gate G3: identities ---------------------------------------------------------------------
cat("\nGate G3 -- n and events add up\n")
for (a in SPEC$analysis) {
  d <- MET[MET$analysis == a, ]
  nh <- d$n[d$region == "H"]; nc <- d$n[d$region == "Hc"]; nt <- d$n_total[1]
  stop_if(nh + nc == nt, "G3 FAIL: ", a, ": n(H) ", nh, " + n(Hc) ", nc, " != n_total ", nt)
  eh <- d$events[d$region == "H"]; ec <- d$events[d$region == "Hc"]; et <- d$events_total[1]
  if (!is.na(eh) && !is.na(ec) && !is.na(et))
    stop_if(eh + ec == et, "G3 FAIL: ", a, ": events(H) ", eh, " + events(Hc) ", ec, " != ", et)
  cat(sprintf("  %-19s %d + %d = %d%s\n", a, nh, nc, nt,
              if (!is.na(eh)) sprintf("; events %d + %d = %d", eh, ec, et) else "; no events (continuous)"))
}
wr(MET, "applications_metrics.csv")

# ---- 2. applications_funnels.csv --------------------------------------------------------------
fun <- list()
addf <- function(analysis, stage, count, source, note = NA, file) {
  fun[[length(fun) + 1L]] <<- data.frame(analysis = analysis, stage = stage, count = count,
                                         source = source, note = note,
                                         source_payload = file, source_md5 = md5(file),
                                         extract_built_at = BUILT_AT, stringsAsFactors = FALSE)
}
# the forest-search GBSG funnel is in neither its payload nor the criteria log
clog  <- "dev/payloads/criteria/gbsg/criteria_log.txt"
ltxt  <- readLines(clog, warn = FALSE)
# it records the funnel only if it carries a count against each of the four stage words
log_has <- all(vapply(c("enumerat", "fittable", "clear|screening threshold", "passing"),
                      function(w) any(grepl(w, ltxt, ignore.case = TRUE) & grepl("[0-9]", ltxt)), NA))
if (log_has) {
  stop("criteria_log.txt now records the funnel: take the counts from it (source = criteria_log)", call. = FALSE)
} else {
  rpt <- readLines(SPEC$report[SPEC$analysis == "gbsg_fs"], warn = FALSE)
  kf <- grep("The funnel\\*\\* at this fit", rpt)
  stop_if(length(kf) == 1L, "the FS GBSG funnel line is not in the report exactly once")
  m <- regmatches(rpt[kf], regexec(
    "([\\d,]+) enumerated, ([\\d,]+) fittable, ([\\d,]+) at (HR .[^,]*), \\*\\*([\\d,]+) passing", rpt[kf], perl = TRUE))[[1]]
  stop_if(length(m) == 6L, "the FS GBSG funnel line did not parse")
  n4 <- as.integer(gsub(",", "", m[c(2, 3, 4, 6)]))
  fsrep <- SPEC$report[SPEC$analysis == "gbsg_fs"]
  nt <- sprintf("%s line %d", fsrep, kf)
  for (j in seq_along(n4))
    addf("gbsg_fs", c("enumerated", "fittable", paste0("at ", m[5]), "passing")[j], n4[j], "report", nt, fsrep)
}
for (i in seq_len(nrow(SPEC))) {
  a <- SPEC$analysis[i]; p <- P[[a]]; f <- SPEC$file[i]
  if (!is.null(p$meta$filter_counts)) {
    d <- p$meta$filter_counts
    for (j in seq_len(nrow(d))) addf(a, d$stage[j], d$n[j], "meta$filter_counts", NA, f)
  }
  if (!is.null(p$meta$funnel)) {
    d <- p$meta$funnel
    for (j in seq_len(nrow(d))) addf(a, d$stage[j], d$n[j], "meta$funnel", NA, f)
  }
  for (nm in c("dina_surface", "grf_detail")) if (!is.null(p$extras[[nm]])) {
    s <- p$extras[[nm]]
    for (k in names(s)) addf(a, k, s[[k]], paste0("extras$", nm), NA, f)
  }
  addf(a, "mr_family", p$meta$timings$fs$n_family, "meta$timings$fs$n_family", NA, f)
}
FUN <- do.call(rbind, fun)
FUN <- FUN[order(match(FUN$analysis, SPEC$analysis)), ]
wr(FUN, "applications_funnels.csv")

# ---- 3. the trio summary's agreement, partition and memberships -------------------------------
wr(prov(TR$agreement$pairwise,  TRIO), "gbsg_agreement_pairwise.csv")
wr(prov(TR$agreement$partition, TRIO), "gbsg_agreement_partition.csv")
wr(prov(TR$memberships$counts,  TRIO), "gbsg_memberships.csv")

# ---- 4/5. the overlap files, from extras$agreement --------------------------------------------
overlap <- function(p, file, comps, defs_want, extra = list()) {
  o <- list()
  for (cmp in comps) {
    s <- p$extras$agreement[[cmp]]; d <- s$table
    d$comparison <- cmp
    for (nm in c("n_total", "n_A", "n_B", "proportion_agreeing", "kappa", "jaccard")) d[[nm]] <- s[[nm]]
    dd <- p$extras$agreement$definitions
    for (nm in defs_want) d[[nm]] <- unname(dd[[nm]])
    for (nm in names(extra)) d[[nm]] <- extra[[nm]]
    o[[length(o) + 1L]] <- d
  }
  prov(do.call(rbind, o), file)
}
fs_file <- SPEC$file[SPEC$analysis == "gbsg_fs"]
wr(overlap(P$gbsg_fs, fs_file, c("maxNeps_vs_maxN", "maxNeps_vs_maxbeta"),
           c("max_N_eps", "max_N", "max_beta")), "gbsg_rules_overlap.csv")
bn_file <- SPEC$file[SPEC$analysis == "actg175_binary"]
ab <- P$actg175_binary$extras$agreement
wr(overlap(P$actg175_binary, bn_file, c("maxNeps_vs_grf", "maxNeps_vs_dina"),
           c("max_N_eps", "published_grf", "published_dina"),
           extra = list(published_n_fs = unname(ab$published_n[["fs"]]),
                        published_n_grf = unname(ab$published_n[["grf"]]),
                        published_n_dina = unname(ab$published_n[["dina"]]),
                        notes = paste(ab$notes, collapse = "; "))),
   "actg175_binary_published_overlap.csv")

# ---- 6. COLUMNS.md ----------------------------------------------------------------------------
cm <- c(
"# COLUMNS -- the applications extract, and what it makes of the six committed payloads",
"",
sprintf("*Generated by `dev/payloads/R/applications_metrics.R` together with the CSV files beside it, %s. Regenerate with the script; do not edit by hand.*", BUILT_AT),
"",
"**Five analyses** -- `gbsg_fs`, `gbsg_dina`, `gbsg_grf`, `actg175_binary`, `actg175_continuous` -- read from the six committed payloads (the five analysis payloads and the GBSG trio summary). The analysis key is fixed by file name. No trial data is read and no analysis is re-run: every value here is a field of a committed payload, a field of a committed report where the payload has none (the forest-search GBSG funnel), or a derivation stated below.",
"",
"## Provenance, on every row of every file",
"",
"| column | what it holds |",
"|---|---|",
"| `source_payload` | the repository path of the file the row was read from. Normally a payload; for the `gbsg_fs` funnel rows it is the report, because the count is in no payload (`source` says which) |",
"| `source_md5` | `tools::md5sum()` of that file. A change to any source payload shows here |",
"| `extract_built_at` | when this extract was built |",
"",
"## 1. `applications_metrics.csv` -- ten rows, one per analysis x region",
"",
"Region rows are the `table$region` values `H` and `Hc`, as the trio summary reads them. Wide: every quantity is a column. A cell is empty where the payload has no such field; the last column below says which payloads carry it.",
"",
"| column | what it holds | payload field, and any derivation |",
"|---|---|---|",
"| `analysis` | the analysis key | fixed by file name |",
"| `trial` | `GBSG` or `ACTG175` | fixed by file name |",
"| `identifier` | `FS`, `DINA` or `GRF` | fixed by file name |",
"| `outcome_type` | `binary` / `continuous` | `meta$settings$outcome_type`. **The three GBSG payloads do not record it**, so their cells are empty |",
"| `est_scale` | the effect scale the payload records | `est_scale`, **verbatim**. The payloads record `hr`, `hr`, `hr`, `or` and `md_decline_oriented` -- lower case, and the continuous one carries its orientation in the name |",
"| `adverse_outcome` | the orientation flag | `meta$settings$adverse_outcome`; **the three GBSG payloads do not record it** (empty) |",
"| `endpoint`, `outcome` | the endpoint and the outcome column | `meta$endpoint`, `meta$outcome` |",
"| `region` | `H` or `Hc` | `table$region` |",
"| `region_def` | the region in the search's own coding | `labels$fs` (the three forest-search analyses), `labels$dina`, `labels$grf`. **Derived on `Hc`:** `complement of <def>` |",
"| `region_label` | the plain-language description the analysis's report prints for the selection | held by the builder and **asserted verbatim against that report** (Gate G2). **Derived on `Hc`:** `complement of <label>` |",
"| `n`, `pct` | size of the region and its percentage of the trial | `table$n`, `table$pct` (the payload's own rounding) |",
"| `events` | events in the region | **`H`:** `extras$harm_rates$treated_events + control_events`. **`Hc`:** `events_total - events(H)`. Empty on `actg175_continuous`, which has no events |",
"| `trt_n`, `ctrl_n` | arm sizes | **`H`:** `extras$harm_rates$treated_n` / `control_n`, and on `actg175_continuous`, which has no `harm_rates`, `extras$arm_means$treated_n` / `control_n`. **`Hc`:** the whole trial's `n_trt` / `n_ctrl` less `H`'s |",
"| `trt_events`, `ctrl_events` | events by arm | **`H`:** `extras$harm_rates`. **`Hc`:** the whole trial's `ev_trt` / `ev_ctrl` less `H`'s. Empty on the continuous analysis |",
"| `trt_pct`, `ctrl_pct` | the event percentage in each arm | **`H`:** `extras$harm_rates$treated_pct` / `control_pct` as recorded. **`Hc`: recomputed from the derived counts**, `100 * events / n` |",
"| `trt_mean`, `ctrl_mean` | mean CD4 decline by arm, the continuous analysis only | **`H`:** `extras$arm_means$treated_mean_decline` / `control_mean_decline`. **`Hc` derived** from the whole trial's arm means and counts: `(n_arm * mean_arm - n_arm(H) * mean_arm(H)) / (n_arm - n_arm(H))`, the whole-trial arm means being `extras$reference_md$mean_trt` / `mean_ctrl` on the same oriented scale |",
"| `n_total`, `events_total` | the whole trial, repeated on both rows | `meta$n_total`; `ev_trt + ev_ctrl` of the `whole trial` row of `extras$reference_hr` / `reference_or`. Empty for the continuous analysis |",
"| `whole_trial_est`, `_lo`, `_hi` | the whole-trial effect, repeated on both rows | the `whole trial` row of `extras$reference_hr` / `reference_or` / `reference_md` |",
"| `naive_est`, `naive_lo`, `naive_hi` | the unadjusted analysis, two-sided Wald | `table` |",
"| `mr_est`, `mr_lo`, `mr_hi` | the two-term corrected estimate with its IJ interval | `table` |",
"| `est2`, `lambda_mean` | the field's second-order estimate and the mean of the simulated root | **`H`:** `bounds$field_est2_H`, `field_lambda_mean_H`. **`Hc`:** `field_s_est2_Hc`, `field_s_lambda_mean_Hc` -- the **field-s** pair, which is the reported complement bound's |",
"| `est2_unstudentized`, `lambda_mean_unstudentized` | the unstudentized complement pair, `Hc` only | `bounds$field_est2_Hc`, `field_lambda_mean_Hc` |",
"| `first_correction` | what the two-term correction removed | **derived:** `log(naive_est) - log(mr_est)` on the two ratio scales, `naive_est - mr_est` on the CD4 scale |",
"| `bound_side` | `lower` on `H`, `upper` on `Hc` | fixed |",
"| `field_bound` | the reported one-sided bound | `bounds$field_lower_H` / `field_s_upper_Hc` |",
"| `bonf_bound` | the Bonferroni pair's bound on this side | `bounds$bonf_lower_H` / `bonf_upper_Hc` |",
"| `bonf_joint_prob`, `calib_gamma`, `calib_joint_prob`, `calib_lower_H`, `calib_upper_Hc`, `joint_corr`, `rho_c`, `alpha` | the joint quantities, repeated on both rows | `bounds`. `rho_c` is **NA on `actg175_continuous`**, whose source specification does not pass `field_decompose` |",
"| `n_out`, `n_in_mean` | the field's Monte Carlo accounting | `bounds$n_out_H` / `n_out_Hc`; `n_in_H` on the `H` row only |",
"| `ordering_holds` | whether the second-order estimate sits between the bound and the two-term estimate | **derived:** on `H`, `field_bound <= est2 & est2 <= mr_est`; on `Hc`, `mr_est <= est2 & est2 <= field_bound` |",
"| `selection_bias`, `n_family` | the correction's size and the family it re-selects over, repeated on both rows | `meta$timings$fs$selection_bias`, `$n_family` |",
"| `sg_focus`, `eps`, `selection_rule` | the selection rule | `meta$settings$sg_focus`, `effect_neighborhood`, `selection_rule`. `gbsg_fs` does not record `selection_rule` (empty) |",
"| `c1`, `c2` | the screening and per-split thresholds | `effect.threshold` where present, else `hr.threshold`; `consistency.threshold` where present, else `hr.consistency`. **`c2` is `NA` for `gbsg_dina` and `gbsg_grf`**: consistency is forest-search-only in the package, and neither document passes it |",
"| `pstar` | the consistency proportion | `pconsistency.threshold`; `NA` for DINA and GRF, for the same reason |",
"| `screen_applies` | whether the splitting-consistency screen applies | fixed: `TRUE` for the three forest-search analyses, `FALSE` for DINA and GRF |",
"| `consistency_method`, `n_min`, `d_min`, `maxk` | the eligibility settings | `meta$settings`, and failing that `meta$thresholds_received`. `d_min` is `d0.min`, the per-arm minimum. `gbsg_fs` records none of these in either place (empty); `consistency_method` is recorded by the two ACTG 175 payloads only |",
"| `mr_draws`, `field_R_out`, `field_R_in`, `seed`, `n_cores` | the multiplier and field settings | `meta$mr_draws`, `meta$settings$field_R_out` / `field_R_in` / `seedit`, `meta$cores` |",
"| `payload_built_at`, `forestsearch_version` | when the payload was built, and under which package version | `built_at`, `forestsearch_version` |",
"",
"## 2. `applications_funnels.csv` -- long, one row per analysis x stage",
"",
"| column | what it holds |",
"|---|---|",
"| `analysis` | the analysis key |",
"| `stage` | the stage name **as the source records it**, never renamed |",
"| `count` | the count at that stage |",
"| `source` | `meta$filter_counts` (the two ACTG 175 analyses), `meta$funnel` (DINA and GRF), `extras$dina_surface` / `extras$grf_detail` (their own rows, `stage` being the field name), `meta$timings$fs$n_family` (the `mr_family` row every analysis gets), or `report` |",
"| `note` | for a `report` row, the file and line the count was read from |",
"",
sprintf("**The forest-search GBSG funnel is in no payload.** `%s` records the six fits' timings and their candidate-family comparison, and **not** the enumerated, fittable, screen-clearing and consistency-passing counts, so the four rows are read from the report instead, with `source = report` and the line in `note`.", clog),
"",
"## 3. `gbsg_agreement_pairwise.csv`, `gbsg_agreement_partition.csv`, `gbsg_memberships.csv`",
"",
"The trio summary payload's `agreement$pairwise` (three rows: `pair`, `both_H`, `first_only`, `second_only`, `both_complement`, `agree`, `kappa`, `jaccard`), `agreement$partition` (eight rows: `region`, `N`, `events`, `description`) and `memberships$counts` (three rows: `identifier`, `N`, `events`, `hormonal`, `none`), each written as it stands with the provenance columns added. Nothing is derived.",
"",
"## 4. `gbsg_rules_overlap.csv`",
"",
"The forest-search GBSG payload's `extras$agreement`, for `maxNeps_vs_maxN` and `maxNeps_vs_maxbeta`: each comparison's four-cell table (`cell`, `n`, `n_ctrl`, `n_trt`, `ev_ctrl`, `ev_trt`, `hr`, `lo`, `hi`) as rows, with `comparison`, `n_total`, `n_A`, `n_B`, `proportion_agreeing`, `kappa` and `jaccard` repeated on each row and the three `definitions` (`max_N_eps`, `max_N`, `max_beta`) as columns. A cell too small to fit carries `NA` for `hr`, `lo` and `hi`, as the payload records it.",
"",
"## 5. `actg175_binary_published_overlap.csv`",
"",
"The ACTG 175 binary payload's `extras$agreement`, for `maxNeps_vs_grf` and `maxNeps_vs_dina`, in the same shape, with `definitions` (`max_N_eps`, `published_grf`, `published_dina`), `published_n_fs` / `_grf` / `_dina` and `notes` carried. **This file serves the supplement's S5; the long form does not quote it.**",
"",
"## 6. What is not here",
"",
"- **Per-subject membership.** Only the counts and the partition; no payload records which subject is in which region, except through the trio summary's own cross-tabulation, which is what `gbsg_agreement_*.csv` carries.",
"- **The panitumumab analyses.** Not among the five: their data are licensed and their placement is undecided.",
"- **Any operating characteristic.** This extract is five single fits; coverage, bias and declaration rates live in `dev/payloads/identifiers/`.",
"- **The full bootstrap.** `table$fb_est` / `fb_lo` / `fb_hi` are `NA` in every payload -- it was never run -- so no column reads them.",
"")
writeLines(cm, file.path(out_dir, "COLUMNS.md"), useBytes = FALSE)
cat(sprintf("wrote %-38s %3d lines\n", "COLUMNS.md", length(cm)))
cat("\ndone.\n")
