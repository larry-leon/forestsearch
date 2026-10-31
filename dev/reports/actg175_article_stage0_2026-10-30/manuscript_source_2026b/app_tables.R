# app_tables.R -- the applications tables, main paper (Section 6) and Supplementary Section S5.
#
# Every table is built from the applications extract, dev/payloads/applications/*.csv (built by
# dev/payloads/R/applications_metrics.R), and returned as a character data frame; the section files
# print them with knitr::kable. No number is typed here. The two counts no payload carries -- the
# frontier sizes of the breast cancer forest search and of the ACTG 175 binary search -- are read at
# render from the committed criteria payload and from the binary analysis's report, the same sources
# dev/writeup/applications/R/verify_applications.R addresses.
#
#   A <- app_load(extract_dir)
#   app_tbl_gbsg(A); app_tbl_five(A); app_tbl_funnels(A, criteria_csv, binary_report)
#   app_tbl_pairwise(A); app_tbl_partition(A); app_tbl_diag(A); app_tbl_location(A)
#
# Base R only. dev/writeup/applications/R/verify_applications.R sources this file as well, to read the
# cells the documents draw.

app_load <- function(dir) {
  if (!dir.exists(dir)) stop("applications extract not found: ", dir, call. = FALSE)
  rd <- function(f) utils::read.csv(file.path(dir, f), stringsAsFactors = FALSE, check.names = FALSE,
                                    encoding = "UTF-8")
  list(M = rd("applications_metrics.csv"), FUN = rd("applications_funnels.csv"),
       PW = rd("gbsg_agreement_pairwise.csv"), PT = rd("gbsg_agreement_partition.csv"),
       order = c("gbsg_fs", "gbsg_dina", "gbsg_grf", "actg175_binary", "actg175_continuous"),
       gbsg = c("gbsg_fs", "gbsg_dina", "gbsg_grf"))
}

.app_MINUS <- "−"
.app_EM    <- "—"
.app_g     <- function(A, a, r, col) A$M[[col]][A$M$analysis == a & A$M$region == r]
.app_ratio <- function(A, a) .app_g(A, a, "H", "est_scale") %in% c("hr", "or")
.app_signed <- function(x, d) paste0(if (isTRUE(x < 0)) .app_MINUS else "",
                                     formatC(abs(x), format = "f", digits = d, big.mark = ","))
# three decimals on the ratio scales, one on the CD4 scale
.app_num <- function(A, x, a) if (is.na(x)) .app_EM else .app_signed(x, if (.app_ratio(A, a)) 3L else 1L)
.app_cnt <- function(x) if (is.na(x)) .app_EM else formatC(as.integer(round(x)), format = "d", big.mark = ",")
.app_n_ev <- function(A, a, r, n = "n", e = "events") {
  ev <- .app_g(A, a, r, e); nn <- .app_g(A, a, r, n)
  if (is.na(ev)) .app_cnt(nn) else sprintf("%s (%s)", .app_cnt(nn), .app_cnt(ev))
}
.app_est_ci <- function(A, a, r, p) sprintf("%s (%s, %s)", .app_num(A, .app_g(A, a, r, paste0(p, "_est")), a),
                                            .app_num(A, .app_g(A, a, r, paste0(p, "_lo")), a),
                                            .app_num(A, .app_g(A, a, r, paste0(p, "_hi")), a))
# the short region label: the extract's plain-language label with the receptors abbreviated
.app_short <- function(A, a) {
  s <- .app_g(A, a, "H", "region_label")
  s <- gsub("estrogen receptor", "ER", s, fixed = TRUE)
  s <- gsub("progesterone receptor", "PR", s, fixed = TRUE)
  s <- gsub("baseline CD4 count", "CD4", s, fixed = TRUE)
  s <- gsub(" fmol/l", "", s, fixed = TRUE); s <- gsub(" per mm³", "", s, fixed = TRUE)
  gsub(", ", " and ", s, fixed = TRUE)
}
.app_head <- function(A, a, scale = TRUE) {
  gb <- .app_g(A, a, "H", "trial") == "GBSG"
  wh <- if (gb) .app_g(A, a, "H", "identifier") else .app_g(A, a, "H", "outcome_type")
  sc <- switch(.app_g(A, a, "H", "est_scale"), hr = "HR", or = "OR", "MD, CD4 cells")
  paste0(if (gb) "GBSG, " else "ACTG 175, ", wh, if (scale) sprintf(" (%s)", sc) else "")
}
.app_frame <- function(rows, heads) {
  d <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(d) <- c(" ", heads); rownames(d) <- NULL
  d
}
.app_body <- function(A, an, whole = FALSE, ij = FALSE, detail = TRUE) {
  ac <- function(f) vapply(an, f, "")
  est <- function(r) if (ij) ac(function(a) .app_est_ci(A, a, r, "mr"))
                     else ac(function(a) .app_num(A, .app_g(A, a, r, "mr_est"), a))
  lab <- if (ij) "Corrected (two-term), IJ interval" else "Corrected"
  c(if (whole) list(
      c("Whole trial N (events)", ac(function(a) .app_n_ev(A, a, "H", "n_total", "events_total"))),
      c("Whole-trial estimate, 95% CI", ac(function(a) .app_est_ci(A, a, "H", "whole_trial")))),
    list(
      c("Region", ac(function(a) .app_short(A, a))),
      c("N (events)", ac(function(a) .app_n_ev(A, a, "H"))),
      c("Unadjusted, 95% CI", ac(function(a) .app_est_ci(A, a, "H", "naive"))),
      c(lab, est("H")),
      c("Field, lower bound", ac(function(a) .app_num(A, .app_g(A, a, "H", "field_bound"), a))),
      c("Complement N (events)", ac(function(a) .app_n_ev(A, a, "Hc"))),
      c("Unadjusted, 95% CI", ac(function(a) .app_est_ci(A, a, "Hc", "naive"))),
      c(lab, est("Hc")),
      c("Field-s, upper bound", ac(function(a) .app_num(A, .app_g(A, a, "Hc", "field_bound"), a))),
      c("Bonferroni pair", ac(function(a)
        sprintf("%s / %s", .app_num(A, .app_g(A, a, "H", "bonf_bound"), a),
                .app_num(A, .app_g(A, a, "Hc", "bonf_bound"), a))))),
    if (detail) list(
      c("Selection bias", ac(function(a) .app_num(A, .app_g(A, a, "H", "selection_bias"), a))),
      c("Re-selection family", ac(function(a) .app_cnt(.app_g(A, a, "H", "n_family"))))))
}

# Table 1 of the main paper: the three breast cancer analyses, ten rows; the whole trial is in the
# caption (app_gbsg_whole). The selection bias and the re-selection family are not shown here
# (Larry, 2026-09-25); the supplement's five-analysis table keeps both.
app_tbl_gbsg <- function(A)
  .app_frame(.app_body(A, A$gbsg, detail = FALSE),
             vapply(A$gbsg, function(a) .app_g(A, a, "H", "identifier"), ""))

# the whole-trial line of Table 2's caption
app_gbsg_whole <- function(A) sprintf("%s patients, %s events, hazard ratio %s",
  .app_cnt(.app_g(A, "gbsg_fs", "H", "n_total")), .app_cnt(.app_g(A, "gbsg_fs", "H", "events_total")),
  .app_est_ci(A, "gbsg_fs", "H", "whole_trial"))

# S5: the five analyses side by side, with the whole-trial rows and the IJ intervals
app_tbl_five <- function(A)
  .app_frame(.app_body(A, A$order, whole = TRUE, ij = TRUE), vapply(A$order, function(a) .app_head(A, a), ""))

# S5: the candidate funnels on common stage names. Each cell names the extract stage it reads.
app_tbl_funnels <- function(A, criteria_csv, binary_report) {
  fn <- function(a, stage) { d <- A$FUN[A$FUN$analysis == a & A$FUN$stage == stage, ]
    if (nrow(d) != 1L) stop("funnel stage not found exactly once: ", a, " / ", stage, call. = FALSE)
    .app_cnt(d$count) }
  if (!file.exists(criteria_csv)) stop("criteria payload not found: ", criteria_csv, call. = FALSE)
  cr <- utils::read.csv(criteria_csv, stringsAsFactors = FALSE)
  fr_gbsg <- .app_cnt(sum(as.logical(cr$on_frontier)))
  if (!file.exists(binary_report)) stop("binary analysis report not found: ", binary_report, call. = FALSE)
  rp <- paste(readLines(binary_report, warn = FALSE), collapse = " ")
  m <- regmatches(rp, regexec("frontier holds (\\d+) of the (\\d+) candidates", rp))[[1]]
  if (length(m) != 3L) stop("the binary frontier count was not found in ", binary_report, call. = FALSE)
  e <- .app_EM
  rows <- list(
    c("Enumerated / proposed", fn("gbsg_fs", "enumerated"), fn("gbsg_grf", "n_pool"),
      fn("gbsg_dina", "n_candidates_searched"), fn("actg175_binary", "candidates enumerated"),
      fn("actg175_continuous", "candidates enumerated")),
    c("Fittable", fn("gbsg_fs", "fittable"), e, e, fn("actg175_binary", "fittable"),
      fn("actg175_continuous", "passing the per-arm minimum (10 / 10)")),
    c("Qualifying on the surface", e, e, fn("gbsg_dina", "n_candidates_qualifying"), e, e),
    c("Clearing the screening threshold", fn("gbsg_fs", "at HR ≥ 1.00"), fn("gbsg_grf", "admitted_n"),
      fn("gbsg_dina", "admitted_n"), fn("actg175_binary", "clearing the OR screen (≥ 1.00)"),
      fn("actg175_continuous", "clearing the screen (MD ≥ 10, oriented)")),
    c("Evaluated for consistency", e, e, e, fn("actg175_binary", "evaluated for consistency"),
      fn("actg175_continuous", "evaluated for consistency")),
    c("Passing the consistency screen", fn("gbsg_fs", "passing"), e, e,
      fn("actg175_binary", "passing consistency (p* = 0.90)"),
      fn("actg175_continuous", "passing consistency (c2 = 5, p* = 0.90)")),
    c("On the frontier", fr_gbsg, e, e, .app_cnt(as.numeric(m[2])),
      fn("actg175_continuous", "on the Pareto frontier")),
    c("Re-selection family", vapply(c("gbsg_fs", "gbsg_grf", "gbsg_dina", "actg175_binary", "actg175_continuous"),
                                    function(a) .app_cnt(.app_g(A, a, "H", "n_family")), "")))
  d <- .app_frame(rows, c("GBSG, FS", "GBSG, GRF", "GBSG, DINA", "ACTG 175, binary", "ACTG 175, continuous"))
  names(d)[1] <- "Stage"
  d
}

app_tbl_pairwise <- function(A) {
  P <- A$PW
  data.frame(Pair = P$pair, `Both selected` = vapply(P$both_H, .app_cnt, ""),
             `First only` = vapply(P$first_only, .app_cnt, ""), `Second only` = vapply(P$second_only, .app_cnt, ""),
             `Both complement` = vapply(P$both_complement, .app_cnt, ""),
             Agreement = sprintf("%.1f%%", 100 * P$agree), Kappa = formatC(P$kappa, format = "f", digits = 3),
             Jaccard = formatC(P$jaccard, format = "f", digits = 3), check.names = FALSE, stringsAsFactors = FALSE)
}

app_tbl_partition <- function(A) {
  P <- A$PT
  data.frame(Region = P$region, N = vapply(P$N, .app_cnt, ""), Events = vapply(P$events, .app_cnt, ""),
             `Who these patients are` = P$description, check.names = FALSE, stringsAsFactors = FALSE)
}

# S5: the supporting quantities. Row labels carry their symbols as math, so no glyph depends on the font.
app_tbl_diag <- function(A) {
  an <- A$order; ac <- function(f) vapply(an, f, "")
  pm <- function(x, d) if (is.na(x)) .app_EM else paste0(if (x >= 0) "+" else .app_MINUS, formatC(abs(x), format = "f", digits = d))
  .app_frame(list(
    c("Mean of the simulated root, region", ac(function(a) pm(.app_g(A, a, "H", "lambda_mean"), if (.app_ratio(A, a)) 3L else 2L))),
    c("Second-order estimate, region", ac(function(a) .app_num(A, .app_g(A, a, "H", "est2"), a))),
    c("Ordering bound $\\le$ second-order $\\le$ two-term", ac(function(a)
      if (isTRUE(as.logical(.app_g(A, a, "H", "ordering_holds")))) "holds" else "does not hold")),
    c("Second-order estimate, complement (field-s)", ac(function(a) .app_num(A, .app_g(A, a, "Hc", "est2"), a))),
    c("Scale ratio $\\rho^c$", ac(function(a) { x <- .app_g(A, a, "H", "rho_c")
      if (is.na(x)) "not computed" else formatC(x, format = "f", digits = 3) })),
    c("Bonferroni pair, simulated joint coverage", ac(function(a) formatC(.app_g(A, a, "H", "bonf_joint_prob"), format = "f", digits = 4))),
    c("Calibrated split $\\gamma$", ac(function(a) formatC(.app_g(A, a, "H", "calib_gamma"), format = "f", digits = 3))),
    c("Correlation of the aligned roots", ac(function(a) pm(.app_g(A, a, "H", "joint_corr"), 3L)))),
    vapply(an, function(a) .app_head(A, a, scale = FALSE), ""))
}

app_tbl_location <- function(A) {
  an <- A$order; ac <- function(f) unname(vapply(an, f, ""))
  data.frame(` ` = ac(function(a) .app_head(A, a)),
             `Field lower, region` = ac(function(a) .app_num(A, .app_g(A, a, "H", "field_bound"), a)),
             `Bonferroni lower, region` = ac(function(a) .app_num(A, .app_g(A, a, "H", "bonf_bound"), a)),
             `Field-s upper, complement` = ac(function(a) .app_num(A, .app_g(A, a, "Hc", "field_bound"), a)),
             `Bonferroni upper, complement` = ac(function(a) .app_num(A, .app_g(A, a, "Hc", "bonf_bound"), a)),
             `Whole trial` = ac(function(a) .app_num(A, .app_g(A, a, "H", "whole_trial_est"), a)),
             check.names = FALSE, stringsAsFactors = FALSE)
}
