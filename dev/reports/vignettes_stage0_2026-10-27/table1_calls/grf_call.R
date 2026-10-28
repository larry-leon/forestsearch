# ---------------------------------------------------------------------------
# Table 1 of León & Anderson (2026b), GRF column -- verbatim transplant.
# Source: fs-glms-interpretable/quarto/gbsg/analysis_gbsg_grf_mr.qmd
# Sibling repo HEAD: 19c6ac41c1198dbc513f5c6978d0d6c5df1fc157
# Chunks (name:start-end lines of the fence, inclusive): label-map:98-166 setup:174-195 fit:241-277
# Data: survival::gbsg, prepared in the 'setup' chunk (no file read).
# Chunk bodies below are copied verbatim; only these header comments and the
# '## ---- chunk: <name> ----' separators were added. Transplanted 2026-10-27.
# ---------------------------------------------------------------------------

## ---- chunk: label-map (lines 99-165) ----
covariate_labels <- data.frame(
  covariate  = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  label      = c("age", "postmenopausal", "tumor size", "grade 3", "positive lymph nodes",
                 "progesterone receptor", "estrogen receptor"),
  unit       = c("years", "", "mm", "", "", "fmol/l", "fmol/l"),
  complement = c(NA, "premenopausal", NA, "grade 1–2", NA, NA, NA),
  binary     = c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE),   # 0/1 indicators
  stringsAsFactors = FALSE)
knitr::kable(transform(covariate_labels, unit = ifelse(nzchar(unit), unit, "—"),
                       complement = ifelse(is.na(complement), "—", complement),
                       binary = ifelse(binary, "yes", "—")),
             col.names = c("Covariate", "Label", "Unit", "Complement label", "0/1 indicator"), row.names = FALSE,
             caption = "Covariate label map (the translation table itself). A 0/1 indicator compared at 1 reads as its label, and at 0 as its complement label.")

.parse_term <- function(term) {
  term <- trimws(term)
  neg  <- startsWith(term, "!")
  body <- sub("^!?\\{(.*)\\}$", "\\1", term)
  m <- regmatches(body, regexec("^([A-Za-z0-9_.]+)\\s*(<=|>=|<|>|==)\\s*(\\S+)$", body))[[1]]
  if (length(m) == 4L) list(var = m[2], op = m[3], value = as.numeric(m[4]), neg = neg)
  else list(var = body, op = NA_character_, value = NA_real_, neg = neg)
}
.flip <- c("<=" = ">", ">=" = "<", "<" = ">=", ">" = "<=")
.sym  <- c("<=" = "≤", ">=" = "≥", "<" = "<", ">" = ">", "==" = "=")
describe_subgroup <- function(def) {
  terms <- strsplit(def, " & ", fixed = TRUE)[[1]]
  out <- vapply(terms, function(tm) {
    p <- .parse_term(tm); i <- match(p$var, covariate_labels$covariate)
    if (is.na(i)) stop("no label for covariate ", p$var)
    if (is.na(p$op)) {
      if (!p$neg) return(covariate_labels$label[i])
      if (is.na(covariate_labels$complement[i])) stop("no complement label for ", p$var)
      return(covariate_labels$complement[i])
    }
    op <- if (p$neg) .flip[[p$op]] else p$op
    # a 0/1 indicator compared at its own levels is a flag, not a magnitude: {grade3 >= 1}
    # reads "grade 3", and {meno <= 0} reads "premenopausal"
    if (isTRUE(covariate_labels$binary[i])) {
      if ((op %in% c(">=", "==") && p$value == 1) || (op == ">" && p$value == 0))
        return(covariate_labels$label[i])
      if ((op %in% c("<=", "==") && p$value == 0) || (op == "<" && p$value == 1)) {
        if (is.na(covariate_labels$complement[i])) stop("no complement label for ", p$var)
        return(covariate_labels$complement[i])
      }
    }
    trimws(sprintf("%s %s %s %s", covariate_labels$label[i], .sym[[op]], format(p$value), covariate_labels$unit[i]))
  }, character(1))
  paste(out, collapse = " and ")
}
subgroup_member <- function(def, d) {
  terms <- strsplit(def, " & ", fixed = TRUE)[[1]]
  Reduce(`&`, lapply(terms, function(tm) {
    p <- .parse_term(tm); x <- d[[p$var]]; stopifnot(!is.null(x))
    v <- if (is.na(p$op)) x == 1 else switch(p$op, "<=" = x <= p$value, ">=" = x >= p$value,
                                               "<" = x < p$value, ">" = x > p$value, "==" = x == p$value)
    if (p$neg) !v else v
  }))
}
# Gate 4.1: a description carries no coding characters and no column name, and denotes the payload's set (by N)
raw_names <- c("meno", "grade3", "pgr", "er", "hormon", "rfstime", "status", "time_months")
check_description <- function(def, N, d) {
  s <- describe_subgroup(def)
  bad <- grepl("[{}!]", s) || any(vapply(raw_names, function(k) grepl(paste0("\\b", k, "\\b"), s), TRUE))
  if (bad) stop("description still carries coding: ", s)
  if (sum(subgroup_member(def, d)) != N) stop("description of ", def, " does not give N = ", N)
  s
}

## ---- chunk: setup (lines 175-194) ----
suppressPackageStartupMessages({
  library(forestsearch)
  library(survival)
})
t_doc <- proc.time()[["elapsed"]]

df <- survival::gbsg
df <- within(df, {
  id <- seq_len(nrow(df))                  # row number, the fit's subject identifier
  time_months <- rfstime / 30.4375         # recurrence-free survival time, days to months
  grade3 <- ifelse(grade == "3", 1, 0)     # grade 3 against grade 1 or 2
  treat <- hormon })                       # 1 = hormonal therapy, 0 = none

confounders <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")
stopifnot(nrow(df) == 686, !anyNA(df), identical(df$treat, df$hormon),
          sum(df$treat == 1) == 246, sum(df$treat == 0) == 440,
          sum(df$status) == 299, sum(df$status == 0) == 387,
          all(df$er >= 0), all(df$pgr >= 0))
c(subjects = nrow(df), hormonal_therapy = sum(df$treat == 1), no_hormonal_therapy = sum(df$treat == 0),
  events = sum(df$status), censored = sum(df$status == 0))

## ---- chunk: fit (lines 242-276) ----
mr_draws <- 5000L
t0 <- proc.time()[["elapsed"]]
fit <- forestsearch(
  df.analysis = df,
  outcome.name = "time_months", event.name = "status", treat.name = "hormon", id.name = "id",
  confounders.name = confounders,
  is.RCT = TRUE, seedit = 8316951, est.scale = "hr",
  use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
  max_n_confounders = 1000, conf_force = c("er <= 0", "pgr <= 0"),
  cut_type = "default", cont.cutoff = 4, conf.cont_jcuts = list(er = 10, pgr = 10),
  collapse_cuts = TRUE,
  subgroup_method = "grf",                                      # GRF's own
  grf_selection = "frontier", grf_select_statistic = "effect",  # GRF's own
  grf_depth = 2L, dmin.grf = 0.0,                               # GRF's own
  sg_focus = "effMaxSG",
  selection_rule = "neighborhood", effect_neighborhood = 0.20,
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
  hr.threshold = 1.00,
  m1.threshold = Inf, minp = 0.025, fs.splits = 1000,
  stop_threshold = NULL,
  mr_inference = TRUE,
  mr_inference_args = list(ci_method = "field", draws = mr_draws, include_complement = TRUE,
                           confirm_rule = "point", field_uniform = FALSE,
                           field_complement = TRUE, ij_residual = "two_term",
                           field_decompose = TRUE, field_scale_complement = "selected",
                           field_recovery = FALSE, return_reselection = TRUE))
fit_sec <- proc.time()[["elapsed"]] - t0
if (fit_sec > 900) stop(sprintf("the fit took %.0f s, over the 15-minute cap", fit_sec))
g <- fit$mr_inference
o <- fit$grp.consistency$out_sg
label_fs <- paste(fit$sg.harm, collapse = " & ")
selected <- length(fit$sg.harm) > 0L && !is.null(g)
cat(sprintf("forestsearch %s; fit including MR: %.1f s. Selection: %s\n",
            packageVersion("forestsearch"), fit_sec,
            if (selected) describe_subgroup(label_fs) else "NONE - see the empty-selection section"))
