# ---- transplant (verbatim chunk bodies; fences removed) ----
# source: fs-glms-interpretable/quarto/actg175/analysis_actg175_continuous_mr.qmd
# fs-glms-interpretable HEAD 19c6ac41c1198dbc513f5c6978d0d6c5df1fc157 (file last changed d97ca1f / 935ae87, 2026-09-23)
#   chunk 'label-map': fence lines 124-191
#   chunk 'setup': fence lines 197-233
#   chunk 'reference-effects': fence lines 240-275
# ----------------------------------------------------------

## ==== chunk 'label-map' ====
covariate_labels <- data.frame(
  covariate  = c("age", "preanti", "wtkg", "karnof", "cd40", "cd80",
                 "hemo", "homo", "drugs", "race", "gender", "symptom"),
  label      = c("age", "prior antiretroviral therapy", "weight", "Karnofsky score",
                 "baseline CD4 count", "baseline CD8 count",
                 "hemophilia", "homosexual activity", "history of intravenous drug use",
                 "non-white", "male", "symptomatic"),
  unit       = c("years", "days", "kg", "", "per mm³", "per mm³", "", "", "", "", "", ""),
  complement = c(NA, NA, NA, NA, NA, NA,
                 "no hemophilia", "no homosexual activity", "no history of intravenous drug use",
                 "white", "female", "asymptomatic"),
  stringsAsFactors = FALSE)
knitr::kable(transform(covariate_labels, unit = ifelse(nzchar(unit), unit, "—"),
                       complement = ifelse(is.na(complement), "—", complement)),
             col.names = c("Covariate", "Label", "Unit", "Complement label"), row.names = FALSE,
             caption = "Covariate label map (the translation table itself), from the speff2trial::ACTG175 documentation.")

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
    trimws(sprintf("%s %s %s %s", covariate_labels$label[i], .sym[[op]], format(p$value), covariate_labels$unit[i]))
  }, character(1))
  paste(out, collapse = " and ")
}
# the same rule, evaluated. The binary covariates are factors, as the applied document builds
# them (line 67), so a flag term compares on the level.
subgroup_member <- function(def, d) {
  terms <- strsplit(def, " & ", fixed = TRUE)[[1]]
  Reduce(`&`, lapply(terms, function(tm) {
    p <- .parse_term(tm); x <- d[[p$var]]; stopifnot(!is.null(x))
    v <- if (is.na(p$op)) as.character(x) == "1" else {
      x <- as.numeric(as.character(x))
      switch(p$op, "<=" = x <= p$value, ">=" = x >= p$value,
             "<" = x < p$value, ">" = x > p$value, "==" = x == p$value)
    }
    if (p$neg) !v else v
  }))
}
# Gate 7.1: a description carries no coding characters and no column name, and denotes the payload's set (by N)
raw_names <- c("preanti", "wtkg", "karnof", "cd40", "cd80", "cd420", "cd4_change", "y_decline",
               "hemo", "homo", "drugs", "race", "gender", "symptom", "str2", "arms")
check_description <- function(def, N, d) {
  s <- describe_subgroup(def)
  bad <- grepl("[{}!]", s) || any(vapply(raw_names, function(k) grepl(paste0("\\b", k, "\\b"), s), TRUE))
  if (bad) stop("description still carries coding: ", s)
  if (sum(subgroup_member(def, d)) != N) stop("description of ", def, " does not give N = ", N)
  s
}

## ==== chunk 'setup' ====
suppressPackageStartupMessages({
  library(forestsearch)
  library(speff2trial)
})
t_doc <- proc.time()[["elapsed"]]

cont_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")     # applied document line 58
bin_vars  <- c("hemo", "homo", "drugs", "race", "gender", "symptom")   # line 59
confounders <- c(cont_vars, bin_vars)                                  # line 60

df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))                # line 62
n_arms13  <- nrow(df)
df$id     <- seq_len(nrow(df))                                         # line 63
df$treat  <- ifelse(df$arms == 1L, 1L, 0L)                             # line 64: ZDV + ddI = 1
n_dropped <- sum(is.na(df$cd420))
df <- df[!is.na(df$cd420), ]                                           # line 65
df$y_decline <- df$cd40 - df$cd420                                     # line 66: a decline, positive = harm
for (v in bin_vars) df[[v]] <- as.factor(df[[v]])                      # line 67

stopifnot(n_arms13 == 1083, nrow(df) == 1083, n_dropped == 0,
          !anyNA(df[c(confounders, "y_decline", "treat")]),
          sum(df$treat == 1) == sum(df$arms == 1L))                     # treat = 1 IS the combination arm

# The harm direction, asserted rather than commented: the oriented mean difference equals
# mean(y_decline | ZDV+ddI) - mean(y_decline | ddI), so a POSITIVE value means a larger CD4
# decline -- a worse outcome -- under ZDV + ddI. Exactly one flip is in effect: y_decline is
# already adverse and adverse_outcome = TRUE applies no further negation.
md_whole <- coef(lm(y_decline ~ treat, data = df))[[2]]
stopifnot(abs(md_whole - (mean(df$y_decline[df$treat == 1]) - mean(df$y_decline[df$treat == 0]))) < 1e-9,
          abs(md_whole + coef(lm(I(cd420 - cd40) ~ treat, data = df))[[2]]) < 1e-9)   # the raw scale is its sign flip
c(subjects = nrow(df), ZDV_ddI = sum(df$treat == 1), ddI = sum(df$treat == 0),
  dropped_missing_week20_CD4 = n_dropped,
  whole_trial_oriented_MD = round(md_whole, 3),
  mean_decline_ZDV_ddI = round(mean(df$y_decline[df$treat == 1]), 1),
  mean_decline_ddI = round(mean(df$y_decline[df$treat == 0]), 1))

## ==== chunk 'reference-effects' ====
# Oriented mean difference on y_decline: positive = worse CD4 outcome under ZDV + ddI.
md_est <- function(d, min_per_arm = 0) {
  n <- tapply(d$y_decline, factor(d$treat, levels = 0:1), length)
  n[is.na(n)] <- 0L
  out <- data.frame(n = nrow(d), n_ctrl = n[["0"]], n_trt = n[["1"]],
                    mean_ctrl = NA_real_, mean_trt = NA_real_,
                    md = NA_real_, lo = NA_real_, hi = NA_real_)
  if (min(n) < max(2, min_per_arm)) return(out)
  out$mean_ctrl <- mean(d$y_decline[d$treat == 0])
  out$mean_trt  <- mean(d$y_decline[d$treat == 1])
  f  <- lm(y_decline ~ treat, data = d)
  b  <- coef(f)[[2]]; se <- sqrt(diag(vcov(f)))[[2]]
  out$md <- b; out$lo <- b - qnorm(0.975) * se; out$hi <- b + qnorm(0.975) * se
  out
}
f1  <- function(x) ifelse(is.na(x), "—", sprintf("%.1f", x))
fci <- function(e, l, u) ifelse(is.na(e), "—", sprintf("%.1f (%.1f, %.1f)", e, l, u))

ref_sets <- list(
  "whole trial"                    = rep(TRUE, nrow(df)),
  "baseline CD4 > 507 per mm³"     = df$cd40 > 507,
  "baseline CD4 ≤ 507 per mm³"     = df$cd40 <= 507,
  "male"                           = as.character(df$gender) == "1",
  "female"                         = as.character(df$gender) == "0")
reference_md <- do.call(rbind, lapply(names(ref_sets), function(k)
  cbind(set = k, md_est(df[ref_sets[[k]], ]))))
rownames(reference_md) <- NULL
knitr::kable(data.frame(set = reference_md$set, N = reference_md$n,
                        `ZDV+ddI / ddI` = sprintf("%d / %d", reference_md$n_trt, reference_md$n_ctrl),
                        `mean CD4 decline (ZDV+ddI / ddI)` = sprintf("%s / %s", f1(reference_md$mean_trt), f1(reference_md$mean_ctrl)),
                        `oriented MD (95% CI)` = fci(reference_md$md, reference_md$lo, reference_md$hi),
                        check.names = FALSE),
             row.names = FALSE,
             caption = "Mean difference in CD4 decline to week 20, CD4 cells per mm³: positive means a worse CD4 outcome under ZDV + ddI. A negative value is benefit from the combination.")

