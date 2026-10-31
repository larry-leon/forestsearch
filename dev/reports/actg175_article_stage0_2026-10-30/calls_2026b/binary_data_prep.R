# ---- transplant (verbatim chunk bodies; fences removed) ----
# source: fs-glms-interpretable/quarto/actg175/analysis_actg175_binary_mr.qmd
# fs-glms-interpretable HEAD 19c6ac41c1198dbc513f5c6978d0d6c5df1fc157 (file last changed d97ca1f / 935ae87, 2026-09-23)
#   chunk 'label-map': fence lines 106-169
#   chunk 'setup': fence lines 177-212
#   chunk 'reference-effects': fence lines 219-255
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
# the same rule, evaluated: membership of each subject in the described subgroup
subgroup_member <- function(def, d) {
  terms <- strsplit(def, " & ", fixed = TRUE)[[1]]
  Reduce(`&`, lapply(terms, function(tm) {
    p <- .parse_term(tm); x <- d[[p$var]]; stopifnot(!is.null(x))
    v <- if (is.na(p$op)) x == 1 else switch(p$op, "<=" = x <= p$value, ">=" = x >= p$value,
                                               "<" = x < p$value, ">" = x > p$value, "==" = x == p$value)
    if (p$neg) !v else v
  }))
}
# Gate 4: a description carries no coding characters and no column name, and denotes the payload's set (by N)
raw_names <- c("preanti", "wtkg", "karnof", "cd40", "cd80", "cd420", "hemo", "homo", "drugs",
               "race", "gender", "symptom", "arms", "y_binary_pos", "y_binary_neg")
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

df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
n_arms13 <- nrow(df)
df$id    <- seq_len(nrow(df))
df$treat <- ifelse(df$arms == 1L, 1L, 0L)         # 1 = ZDV + ddI, 0 = ddI

# Drop subjects with no week-20 CD4 measurement
n_dropped <- sum(is.na(df$cd420))
df <- df[!is.na(df$cd420), ]

# The beneficial scale, then the adverse scale the search runs on
df$y_binary_pos <- as.integer(df$cd420 > df$cd40)   # 1 = CD4 improved
df$y_binary_neg <- 1L - df$y_binary_pos             # 1 = no CD4 improvement (the analysed outcome)

# every candidate covariate numeric: forestsearch builds "var <= value" cuts, and cont.cutoff = 4
# treats the 0/1 indicators as categorical automatically
bin_vars  <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
cont_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
for (v in bin_vars) df[[v]] <- as.numeric(df[[v]])
confounders <- c(cont_vars, bin_vars)

stopifnot(n_arms13 == 1083, n_dropped == 0, nrow(df) == 1083,
          !anyNA(df[c(confounders, "y_binary_neg", "treat")]),
          identical(df$y_binary_neg, 1L - df$y_binary_pos),
          all(df$y_binary_neg %in% 0:1), setequal(unique(df$arms), c(1L, 3L)),
          all(vapply(bin_vars, function(v) setequal(unique(df[[v]]), 0:1), TRUE)))
c(subjects = nrow(df), ZDV_ddI = sum(df$treat == 1), ddI = sum(df$treat == 0),
  no_improvement = sum(df$y_binary_neg), improved = sum(df$y_binary_pos),
  dropped_missing_week20_CD4 = n_dropped)

## ==== chunk 'reference-effects' ====
# Unadjusted odds ratio of the adverse outcome, treatment against control, with a Wald interval.
# A cell is fitted only where each arm has at least min_per_arm adverse outcomes and at least
# min_per_arm non-adverse ones; otherwise the logistic fit is separated or unstable.
or_est <- function(d, min_per_arm = 0, y = "y_binary_neg") {
  ev <- tapply(d[[y]], factor(d$treat, levels = 0:1), sum)
  n  <- tapply(d[[y]], factor(d$treat, levels = 0:1), length)
  ev[is.na(ev)] <- 0L; n[is.na(n)] <- 0L
  non <- n - ev
  out <- data.frame(n = nrow(d), n_ctrl = n[["0"]], n_trt = n[["1"]],
                    ev_ctrl = ev[["0"]], ev_trt = ev[["1"]],
                    or = NA_real_, lo = NA_real_, hi = NA_real_)
  if (min(ev, non) < max(1, min_per_arm)) return(out)
  f  <- glm(d[[y]] ~ d$treat, family = binomial)
  b  <- coef(f)[[2]]; se <- sqrt(diag(vcov(f)))[[2]]
  out$or <- exp(b); out$lo <- exp(b - qnorm(0.975) * se); out$hi <- exp(b + qnorm(0.975) * se)
  out
}
fci <- function(e, l, u, d = 3) ifelse(is.na(e), "—", sprintf("%.*f (%.*f, %.*f)", d, e, d, l, d, u))

ref_sets <- list(
  "whole trial"               = rep(TRUE, nrow(df)),
  "weight > 86 kg"            = df$wtkg > 86,
  "weight ≤ 86 kg"            = df$wtkg <= 86,
  "baseline CD4 > 380 per mm³" = df$cd40 > 380,
  "baseline CD4 ≤ 380 per mm³" = df$cd40 <= 380)
reference_or <- do.call(rbind, lapply(names(ref_sets), function(k)
  cbind(set = k, or_est(df[ref_sets[[k]], ]))))
rownames(reference_or) <- NULL
knitr::kable(data.frame(set = reference_or$set, N = reference_or$n,
                        `ZDV+ddI / ddI` = sprintf("%d / %d", reference_or$n_trt, reference_or$n_ctrl),
                        `no improvement (ZDV+ddI / ddI)` = sprintf("%d / %d", reference_or$ev_trt, reference_or$ev_ctrl),
                        `OR (95% CI)` = fci(reference_or$or, reference_or$lo, reference_or$hi),
                        check.names = FALSE),
             row.names = FALSE,
             caption = "Odds ratio of NO CD4 improvement, ZDV + ddI against ddI (logistic, Wald 95%). OR above 1 is harm.")

