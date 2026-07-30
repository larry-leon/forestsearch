# maxeff_tests.R  --  unified acceptance tests T1-T7 for sg_focus = "maxeff".
#
# Self-contained.  Uses the installed forestsearch package (library(), NOT
# load_all): the fits run consistency splits on multisession workers, which
# require the installed namespace.  Prerequisite: the six maxeff edits are
# applied, document()+install() have run, R restarted.
#
# Run:  Rscript quarto/GuoHe/maxeff_tests.R
# Saves quarto/GuoHe/maxeff_tests.rds and prints one PASS/FAIL line per test.

suppressMessages({
  library(forestsearch)
  library(survival)
  library(future)
  library(data.table)
})

.this_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) return(dirname(normalizePath(f[1])))
  normalizePath(getwd())
})
cat("script directory:", .this_dir, "\n")

# ---- data, identical to fit_sg_focus_baseline.R ---------------------------
df.analysis <- gbsg
df.analysis <- within(df.analysis, {
  id          <- seq_len(nrow(df.analysis))
  time_months <- rfstime / 30.4375
  grade3      <- ifelse(grade == "3", 1, 0)
  treat       <- hormon
})
confounders.name <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")
outcome.name <- "time_months"; event.name <- "status"
id.name      <- "id";          treat.name <- "hormon"

n_cores <- max(1L, floor(0.80 * parallel::detectCores(logical = FALSE)))
cat("cores:", n_cores, "\n\n")

# ---- one fit, config identical to the baseline except the named overrides --
fit_fs <- function(sg_focus, overrides = list()) {
  base <- list(
    df.analysis  = df.analysis,
    outcome.name = outcome.name, event.name = event.name,
    treat.name = treat.name, id.name = id.name,
    confounders.name = confounders.name,
    is.RCT = TRUE, seedit = 8316951, est.scale = "hr",
    use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
    conf_force = c("er <= 0", "pgr <= 0"),
    conf.cont_jcuts = list(er = 10, pgr = 10),
    collapse_cuts = TRUE, cont.cutoff = 4,
    subgroup_method = "consistency",
    sg_focus = sg_focus,
    selection_rule = "neighborhood", effect_neighborhood = 0.10,
    mr_inference = TRUE, mr_inference_args = list(draws = 5000L),
    n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
    max_subgroups_search = 2000,
    hr.threshold = 1.0, hr.consistency = 1.0,
    pconsistency.threshold = 0.90,
    consistency_method = "resample",
    fs.splits = 1000, stop_threshold = 0.95,
    use_twostage = TRUE, max.minutes = 3,
    details = FALSE, quiet = TRUE, plot.sg = FALSE,
    parallel_args = list(plan = "multisession", workers = n_cores,
                         show_message = FALSE)
  )
  args <- utils::modifyList(base, overrides)
  future::plan(future::multisession, workers = n_cores)
  warns <- character(0)
  t0 <- proc.time()
  fs <- withCallingHandlers(
    do.call(forestsearch, args),
    warning = function(w) {
      warns[[length(warns) + 1L]] <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )
  future::plan(future::sequential)
  list(fs = fs, warnings = warns, elapsed = unname((proc.time() - t0)["elapsed"]))
}

# ---- helpers --------------------------------------------------------------
sg_label_of <- function(fs) {
  sg <- fs$sg.harm
  if (is.null(sg) || !length(sg)) NA_character_
  else paste(sort(as.character(sg)), collapse = " & ")
}
sg_n_of <- function(fs) {
  if (!is.null(fs$grp.consistency$sg.harm.id))
    sum(fs$grp.consistency$sg.harm.id == 1L) else NA_integer_
}
# naive within-subgroup Cox HR of fs$sg.harm -- identical logic to the driver.
naive_hr_of <- function(fs) {
  sg <- fs$sg.harm
  if (is.null(sg) || !length(sg) ||
      is.null(fs$grp.consistency$sg.harm.id) ||
      is.null(fs$df.est) ||
      length(fs$grp.consistency$sg.harm.id) != nrow(fs$df.est))
    return(NA_real_)
  dat <- fs$df.est
  pick <- function(c) { h <- c[c %in% names(dat)]; if (length(h)) h[1L] else NA_character_ }
  tcol <- pick(c("Treat", treat.name, "treat", "hormon"))
  ycol <- pick(c("Y", outcome.name, "time_months"))
  ecol <- pick(c("Event", event.name, "status"))
  if (anyNA(c(tcol, ycol, ecol))) return(NA_real_)
  memb <- fs$grp.consistency$sg.harm.id == 1L
  d <- dat[memb, , drop = FALSE]
  tryCatch(exp(unname(coxph(Surv(d[[ycol]], d[[ecol]]) ~ d[[tcol]])$coefficients[[1L]])),
           error = function(e) NA_real_)
}
# does the override warning name all listed overrides?  Patterns are regex so
# the max_subgroups_search value (2000 / 5 / ...) is matched value-agnostically.
warning_names <- function(warns, needles) {
  hit <- warns[grepl("effect maximiser", warns)]
  txt <- paste(hit, collapse = " || ")
  vapply(needles, function(n) grepl(n, txt), logical(1))
}
# normalise a subgroup label for comparison: confs_labels store cuts WITHOUT
# the display braces that fs$sg.harm adds (e.g. "er <= 0" vs "{er <= 0}").
norm_label <- function(s) {
  if (is.na(s)) return(NA_character_)
  s <- gsub("[{}]", "", s)
  s <- gsub("\\s+", " ", trimws(s))
  s
}

results <- list()
PF <- function(ok) if (isTRUE(ok)) "PASS" else "FAIL"

# ===========================================================================
# T1 -- accepted, announces its regime.  Capture the candidate cut matrix Z
#       (and outcome vectors + labels) from forestsearch's own frame so T2 can
#       reconstruct the family independently of the selection machinery.
# ===========================================================================
cat("=== T1: fitting sg_focus = 'maxeff' (with capture) ===\n")
tracer_exit <- quote({
  fr <- NULL
  for (k in sys.nframe():0) {
    e <- tryCatch(sys.frame(k), error = function(err) NULL)
    if (!is.null(e) && exists("conf.screen", envir = e, inherits = FALSE) &&
        exists("df", envir = e, inherits = FALSE)) { fr <- e; break }
  }
  if (!is.null(fr)) {
    nms <- c("df", "conf.screen", "confs_labels", "FSconfounders.name",
             "Y", "Event", "Treat", "id")
    have <- nms[vapply(nms, function(n) exists(n, envir = fr, inherits = FALSE), logical(1))]
    assign(".maxeff_cap", mget(have, envir = fr), envir = globalenv())
  }
})
suppressMessages(trace(forestsearch, exit = tracer_exit, print = FALSE))
r1 <- fit_fs("maxeff")
suppressMessages(untrace(forestsearch))
fs1 <- r1$fs

t1_needles <- c("pconsistency\\.threshold -> 0", "stop_threshold -> NULL",
                "use_twostage -> FALSE", "max_subgroups_search [0-9]+ -> Inf")
t1_named <- warning_names(r1$warnings, t1_needles)
t1_ok <- !is.null(fs1$sg.harm) && length(fs1$sg.harm) > 0 && all(t1_named)
t1_warn_text <- paste(r1$warnings[grepl("effect maximiser", r1$warnings)], collapse = "\n")
cat(sprintf("T1 %s -- sg.harm = %s ; overrides named: %s\n",
            PF(t1_ok), sg_label_of(fs1),
            paste(names(t1_named)[t1_named], collapse = " | ")))
results$T1 <- list(ok = t1_ok, label = sg_label_of(fs1), n = sg_n_of(fs1),
                   named = t1_named, warning = t1_warn_text, warnings = r1$warnings)

# ===========================================================================
# T2 -- maxeff's selection is the global argmax of the enumerated family.
#       Reconstruct the family independently: dummy cut matrix, all combos of
#       size <= maxk (=2), keep only candidates meeting n.min AND per-arm event
#       minima, fit each within-subgroup Cox HR, take the argmax (tie-break K).
# ===========================================================================
cat("\n=== T2: independent argmax reconstruction ===\n")
t2_ok <- NA; t2_report <- list()
if (!exists(".maxeff_cap") || is.null(.maxeff_cap$df) || is.null(.maxeff_cap$conf.screen)) {
  cat("T2 FAIL -- capture of Z/df failed; cannot build independent reference.\n")
  t2_ok <- FALSE
} else {
  cap <- .maxeff_cap
  Zdf <- forestsearch:::dummy(cap$df[, cap$conf.screen, drop = FALSE])
  Z <- as.matrix(Zdf); colnames(Z) <- names(Zdf)
  yy <- as.numeric(cap$Y); dd <- as.numeric(cap$Event); tt <- as.numeric(cap$Treat)
  L <- ncol(Z)
  n.min <- 60L; d0.min <- 10L; d1.min <- 10L

  eval_cand <- function(cols) {
    idx <- if (length(cols) == 1L) Z[, cols] else apply(Z[, cols, drop = FALSE], 1L, prod)
    inx <- idx == 1
    nx <- sum(inx)
    if (nx <= n.min) return(NULL)                          # Status 4: nx > n.min
    d0 <- sum(dd[inx & tt == 0]); d1 <- sum(dd[inx & tt == 1])
    if (d0 < d0.min || d1 < d1.min) return(NULL)           # Status 3: per-arm events
    fit <- try(coxph(Surv(yy[inx], dd[inx]) ~ tt[inx]), silent = TRUE)  # same fitter
    if (inherits(fit, "try-error")) return(NULL)
    cf <- fit$coefficients
    if (length(cf) < 1L || is.na(cf[1L])) return(NULL)
    data.frame(cols = paste(colnames(Z)[cols], collapse = " & "),
               K = length(cols), n = nx, d0 = d0, d1 = d1,
               hr = exp(unname(cf[1L])), stringsAsFactors = FALSE)
  }

  rows <- list()
  for (j in seq_len(L)) { rr <- eval_cand(j); if (!is.null(rr)) rows[[length(rows) + 1L]] <- rr }
  if (L >= 2L) {
    pr <- utils::combn(L, 2L)
    for (p in seq_len(ncol(pr))) { rr <- eval_cand(pr[, p]); if (!is.null(rr)) rows[[length(rows) + 1L]] <- rr }
  }
  fam <- do.call(rbind, rows)
  # maxeff ordering: (-hr, K).
  fam <- fam[order(-fam$hr, fam$K), , drop = FALSE]

  # map internal cut columns (e.g. "q1.1") to user labels via confs_labels.
  lab_of_cols <- function(colstr) {
    cols <- strsplit(colstr, " & ", fixed = TRUE)[[1]]
    labs <- vapply(cols, function(cc) {
      base <- sub("\\.(0|1)$", "", cc)
      dir  <- sub(".*\\.", "", cc)
      m <- match(base, cap$FSconfounders.name)
      lab <- if (!is.na(m)) as.character(cap$confs_labels[m]) else base
      if (identical(dir, "1")) lab else paste0("NOT(", lab, ")")
    }, character(1))
    paste(sort(labs), collapse = " & ")
  }

  top <- fam[1L, ]
  top_label <- lab_of_cols(top$cols)
  # runner-up = first row whose membership label differs from the winner's
  ru_hr <- NA_real_
  for (i in 2:nrow(fam)) {
    if (!identical(lab_of_cols(fam$cols[i]), top_label)) { ru_hr <- fam$hr[i]; break }
  }
  margin <- top$hr - ru_hr

  fs_label <- sg_label_of(fs1)
  t2_ok <- identical(norm_label(top_label), norm_label(fs_label))
  cat(sprintf("T2 %s -- argmax label: %s | maxeff label: %s\n",
              PF(t2_ok), top_label, fs_label))
  cat(sprintf("   argmax HR = %.6f ; runner-up HR = %.6f ; margin = %.6g\n",
              top$hr, ru_hr, margin))
  t2_report <- list(argmax_label = top_label, argmax_cols = top$cols,
                    maxeff_label = fs_label, argmax_hr = top$hr,
                    runnerup_hr = ru_hr, margin = margin,
                    n_family = nrow(fam), L = L,
                    top10 = utils::head(fam, 10))
}
results$T2 <- list(ok = t2_ok, report = t2_report)

# ===========================================================================
# T3 -- hostile gate settings do not move the selection.
# ===========================================================================
cat("\n=== T3: hostile settings ===\n")
r3 <- fit_fs("maxeff", list(
  pconsistency.threshold = 0.99,
  hr.threshold           = 5,
  stop_threshold         = 0.5,
  use_twostage           = TRUE,
  max_subgroups_search   = 5
))
fs3 <- r3$fs
t3_named <- warning_names(r3$warnings, t1_needles)
t3_same  <- identical(sg_label_of(fs3), sg_label_of(fs1))
t3_ok <- t3_same && all(t3_named)
cat(sprintf("T3 %s -- sg.harm = %s (T1 = %s) ; overrides named: %s\n",
            PF(t3_ok), sg_label_of(fs3), sg_label_of(fs1),
            paste(names(t3_named)[t3_named], collapse = " | ")))
if (!t3_same) cat("   *** SELECTION MOVED -- possible seventh gate. ***\n")
results$T3 <- list(ok = t3_ok, label = sg_label_of(fs3), same = t3_same,
                   named = t3_named, warnings = r3$warnings)

# ===========================================================================
# T4 -- regression: hr / effMaxSG / maxSG unchanged vs recorded baseline.
# ===========================================================================
cat("\n=== T4: regression on existing foci ===\n")
base_rds <- file.path(.this_dir, "sg_focus_baseline.rds")
baseline <- if (file.exists(base_rds)) readRDS(base_rds) else NULL
t4_rows <- list(); t4_ok <- TRUE
for (foc in c("hr", "effMaxSG", "maxSG")) {
  rr <- fit_fs(foc)
  lab <- sg_label_of(rr$fs); nn <- sg_n_of(rr$fs); hr <- naive_hr_of(rr$fs)
  exp_lab <- exp_n <- exp_hr <- NA
  if (!is.null(baseline)) {
    b <- baseline$summary[baseline$summary$sg_focus == foc, , drop = FALSE]
    if (nrow(b)) { exp_lab <- b$selected[1]; exp_n <- b$subgroup_n[1]; exp_hr <- b$naive_HR[1] }
  }
  ok <- identical(lab, exp_lab) && isTRUE(nn == exp_n) &&
        (is.na(exp_hr) || isTRUE(abs(hr - exp_hr) < 1e-3))
  t4_ok <- t4_ok && ok
  cat(sprintf("   %-9s %s -- %s (n=%s, HR=%.4f) vs baseline %s (n=%s, HR=%.4f)\n",
              foc, PF(ok), lab, nn, hr, exp_lab, exp_n, exp_hr))
  t4_rows[[foc]] <- list(ok = ok, label = lab, n = nn, hr = hr,
                         exp_label = exp_lab, exp_n = exp_n, exp_hr = exp_hr)
}
cat(sprintf("T4 %s\n", PF(t4_ok)))
results$T4 <- list(ok = t4_ok, rows = t4_rows)

# ===========================================================================
# T5 -- Pcons computed & reported but not filtering: the evaluated-candidate
#       table contains candidates with Pcons below pconsistency.threshold (0.90).
# ===========================================================================
cat("\n=== T5: Pcons retained but not filtering ===\n")
find_pcons_tbl <- function(x, depth = 0) {
  if (depth > 6) return(NULL)
  if (is.data.frame(x) && "Pcons" %in% names(x)) return(x)
  if (is.list(x)) for (el in x) { r <- find_pcons_tbl(el, depth + 1); if (!is.null(r)) return(r) }
  NULL
}
pt <- find_pcons_tbl(fs1$grp.consistency)
t5_ok <- FALSE; t5_min <- NA_real_; t5_nbelow <- NA_integer_; t5_ncand <- NA_integer_
if (!is.null(pt)) {
  pc <- suppressWarnings(as.numeric(pt$Pcons))
  pc <- pc[!is.na(pc)]
  t5_min <- if (length(pc)) min(pc) else NA_real_
  t5_nbelow <- sum(pc < 0.90)
  t5_ncand <- length(pc)
  t5_ok <- isTRUE(t5_nbelow > 0)
}
gc_eval   <- fs1$grp.consistency$n_candidates_evaluated
gc_total  <- fs1$grp.consistency$n_candidates_total
gc_passed <- fs1$grp.consistency$n_passed
cat(sprintf("T5 %s -- table rows = %s ; with Pcons < 0.90 = %s ; min Pcons = %s\n",
            PF(t5_ok), t5_ncand, t5_nbelow,
            ifelse(is.na(t5_min), "NA", sprintf("%.4f", t5_min))))
cat(sprintf("   grp.consistency: evaluated = %s, total = %s, passed = %s\n",
            gc_eval, gc_total, gc_passed))
if (!isTRUE(t5_ok))
  cat("   *** No sub-0.90 rows are STORED: the per-candidate evaluator dropped\n",
      "       them because pconsistency.threshold -> 0 did NOT propagate to\n",
      "       subgroup.consistency (read from args_call_all, never re-synced).\n", sep = "")
results$T5 <- list(ok = t5_ok, n_rows = t5_ncand, n_below = t5_nbelow,
                   min_pcons = t5_min, n_evaluated = gc_eval,
                   n_total = gc_total, n_passed = gc_passed)

# ===========================================================================
# T6 -- gate mapping: the Tier-2 gate re-selects under rule "maxeff".
# ===========================================================================
cat("\n=== T6: gate mapping ===\n")
mr_rule <- tryCatch(fs1$mr_inference$settings$reselection, error = function(e) NULL)
t6_ok <- identical(as.character(mr_rule), "maxeff")
cat(sprintf("T6 %s -- Re-selection rule = %s\n", PF(t6_ok),
            ifelse(is.null(mr_rule), "NULL", mr_rule)))
results$T6 <- list(ok = t6_ok, reselection = mr_rule)

# ===========================================================================
# T7 -- truncation is inert: max_subgroups_search = 5 does not move selection.
# ===========================================================================
cat("\n=== T7: truncation inert ===\n")
r7 <- fit_fs("maxeff", list(max_subgroups_search = 5))
fs7 <- r7$fs
t7_named <- warning_names(r7$warnings, "max_subgroups_search [0-9]+ -> Inf")
t7_same  <- identical(sg_label_of(fs7), sg_label_of(fs1))
t7_ok <- t7_same && all(t7_named)
cat(sprintf("T7 %s -- sg.harm = %s (T1 = %s) ; truncation override named: %s\n",
            PF(t7_ok), sg_label_of(fs7), sg_label_of(fs1), all(t7_named)))
if (!t7_same) cat("   *** SELECTION MOVED under truncation. ***\n")
results$T7 <- list(ok = t7_ok, label = sg_label_of(fs7), same = t7_same, named = t7_named)

# ===========================================================================
# Reporting
# ===========================================================================
saveRDS(results, file.path(.this_dir, "maxeff_tests.rds"))

cat("\n================= SUMMARY =================\n")
for (tt_ in paste0("T", 1:7)) cat(sprintf("%s : %s\n", tt_, PF(results[[tt_]]$ok)))

cat("\n--- T2 detail ---\n")
if (length(results$T2$report)) {
  rep <- results$T2$report
  cat(sprintf("argmax label : %s\n", rep$argmax_label))
  cat(sprintf("maxeff label : %s\n", rep$maxeff_label))
  cat(sprintf("runner-up margin (HR): %.6g  (argmax %.6f vs runner-up %.6f)\n",
              rep$margin, rep$argmax_hr, rep$runnerup_hr))
  cat(sprintf("family size evaluated: %d over L=%d cut columns\n", rep$n_family, rep$L))
}

cat("\n--- T1 override warning text ---\n")
cat(results$T1$warning, "\n")

cat("\nsaved:", file.path(.this_dir, "maxeff_tests.rds"), "\n")
