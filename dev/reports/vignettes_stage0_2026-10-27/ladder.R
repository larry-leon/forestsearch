# Timing ladder, vignettes stage 0 (2026-10-27) -- one rung per invocation.
#   Rscript ladder.R <rung>      (run from this directory; see run_ladder.sh)
# Rungs: S0 S1 S1r S2 S3 (FS), D1 D3 (DINA), G1 G3 (GRF).
# The transplanted call (table1_calls/<id>_call.R) is parsed, every expression
# before `fit <- forestsearch(...)` is evaluated verbatim (label-map + setup
# chunks, and the fit chunk's `mr_draws` / `t0` lines), and the forestsearch()
# call is evaluated with exactly the task's overrides applied to the call
# object. Nothing after the call in the fit chunk is evaluated.

rung <- commandArgs(trailingOnly = TRUE)[1]
spec <- list(
  S0  = list(id = "fs",   mr = FALSE, draws = NULL),
  S1  = list(id = "fs",   mr = TRUE,  draws = 500L),
  S1r = list(id = "fs",   mr = TRUE,  draws = 500L),
  S2  = list(id = "fs",   mr = TRUE,  draws = 2000L),
  S3  = list(id = "fs",   mr = TRUE,  draws = NULL),   # transplanted call's own
  D1  = list(id = "dina", mr = TRUE,  draws = 500L),
  G1  = list(id = "grf",  mr = TRUE,  draws = 500L),
  D3  = list(id = "dina", mr = TRUE,  draws = NULL),
  G3  = list(id = "grf",  mr = TRUE,  draws = NULL))[[rung]]
if (is.null(spec)) stop("unknown rung: ", rung)
dir.create("rungs", showWarnings = FALSE)
dir.create("fits", showWarnings = FALSE)

suppressPackageStartupMessages(library(forestsearch))
cat(sprintf("rung %s | forestsearch %s from %s | built %s\n", rung,
            packageVersion("forestsearch"), find.package("forestsearch"),
            packageDescription("forestsearch")$Built))
stopifnot(packageVersion("forestsearch") == "0.4.0.9000")

# ---- the transplanted call --------------------------------------------------
src <- file.path("table1_calls", paste0(spec$id, "_call.R"))
exprs <- parse(src, keep.source = FALSE)
is_fit <- vapply(exprs, function(e)
  is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("fit")) &&
    is.call(e[[3]]) && identical(e[[3]][[1]], as.name("forestsearch")), logical(1))
stopifnot(sum(is_fit) == 1L)
idx <- which(is_fit)
for (i in seq_len(idx - 1L)) eval(exprs[[i]], envir = globalenv())

cl <- exprs[[idx]][[3]]
call_verbatim <- cl
cl$parallel_args <- quote(list(plan = "sequential"))
cl$quiet <- TRUE
cl$details <- FALSE
cl$plot.sg <- FALSE
if ("plot.grf" %in% names(formals(forestsearch::forestsearch))) cl$plot.grf <- FALSE
if (!spec$mr) cl$mr_inference <- FALSE
if (!is.null(spec$draws)) cl$mr_inference_args$draws <- spec$draws
cat("---- call as evaluated ----\n"); cat(deparse(cl, width.cutoff = 100L), sep = "\n")

arg_vals <- lapply(as.list(cl)[-1], eval, envir = globalenv())
arg_vals$df.analysis <- NULL
saveRDS(list(rung = rung, source = src, call_verbatim = call_verbatim, call_evaluated = cl,
             args_resolved = arg_vals), file.path("rungs", paste0("args_", rung, ".rds")))

# ---- per-candidate p_cons capture (consistency stage only) ------------------
# subgroup.consistency() sets a flag; consistency_resample() and
# .consistency_via_splits() record their result while the flag is up. A
# .consistency_via_splits() call inside the stage is a fallback from the
# closed form (or the "split" method). Observation only: tracing does not
# change any return value.
.ladder_pc <- new.env()
.ladder_pc$in_cons <- FALSE; .ladder_pc$closed <- list(); .ladder_pc$split <- list()
.lbl <- function(pf) {
  l <- tryCatch(get("this.m_label", envir = pf), error = function(e) NA_character_)
  m <- tryCatch(get("m", envir = pf), error = function(e) NA)
  list(label = paste(l, collapse = " & "), m = suppressWarnings(as.integer(m[1])))
}
ns <- asNamespace("forestsearch")
suppressMessages({
  trace("subgroup.consistency", where = ns, print = FALSE,
        tracer = quote(assign("in_cons", TRUE, envir = .ladder_pc)),
        exit   = quote(assign("in_cons", FALSE, envir = .ladder_pc)))
  trace("consistency_resample", where = ns, print = FALSE,
        exit = quote(if (isTRUE(.ladder_pc$in_cons)) {
          .rv <- returnValue(); .id <- .lbl(parent.frame())
          .ladder_pc$closed[[length(.ladder_pc$closed) + 1L]] <- data.frame(
            m = .id$m, label = .id$label, n = .rv$n, d = .rv$d, beta_hat = .rv$beta_hat,
            sigma_D = .rv$sigma_D, rate_closed = .rv$rate_closed)
        }))
  trace(".consistency_via_splits", where = ns, print = FALSE,
        exit = quote(if (isTRUE(.ladder_pc$in_cons)) {
          .id <- .lbl(parent.frame())
          .ladder_pc$split[[length(.ladder_pc$split) + 1L]] <- data.frame(
            m = .id$m, label = .id$label, p_split = as.numeric(returnValue()))
        }))
})

# ---- the fit ----------------------------------------------------------------
warns <- character(0)
gc(reset = TRUE)
t0 <- proc.time()
fit <- withCallingHandlers(eval(cl, envir = globalenv()),
  warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
tt <- proc.time() - t0
mem <- gc()
fit_sec <- unname(tt[["elapsed"]])
cat(sprintf("---- fit elapsed %.2f s (user %.2f, sys %.2f); R max used %.1f Mb\n",
            fit_sec, tt[["user.self"]], tt[["sys.self"]], sum(mem[, ncol(mem)])))
cat(sprintf("---- warnings (%d):\n", length(warns))); if (length(warns)) cat(paste0("  W: ", warns), sep = "\n")
suppressMessages(for (f in c("subgroup.consistency", "consistency_resample", ".consistency_via_splits"))
  untrace(f, where = ns))

# ---- extraction --------------------------------------------------------------
pc_closed <- do.call(rbind, .ladder_pc$closed)
pc_split  <- do.call(rbind, .ladder_pc$split)
label <- paste(fit$sg.harm, collapse = " & ")
est <- fit$df.est
in_H <- if (!is.null(est) && length(fit$sg.harm)) est$treat.recommend[match(df$id, est$id)] == 0 else rep(FALSE, nrow(df))
nH <- sum(in_H); eH <- sum(df$status[in_H]); nHc <- sum(!in_H); eHc <- sum(df$status[!in_H])
gcs <- fit$grp.consistency
g <- fit$mr_inference
q <- list()
if (!is.null(g)) {
  cc <- g$complement; fb <- g$field; fc <- fb$complement; js <- fb$joint_s
  q <- list(naive_H = g$naive$est, naive_lo_H = g$naive$lower, naive_hi_H = g$naive$upper,
            corrected_H = g$debiased$est, field_lower_H = fb$lower_1s,
            naive_Hc = cc$naive$est, naive_lo_Hc = cc$naive$lower, naive_hi_Hc = cc$naive$upper,
            corrected_Hc = cc$debiased$est, field_s_upper_Hc = fc$upper_1s_s,
            bonf_lower_H = js$bonf_lower_H, bonf_upper_Hc = js$bonf_upper_Hc,
            n_family = g$n_family, mr_draws = g$settings$draws, mr_seed = g$settings$seed,
            field_R_out = fb$R_out, field_R_in = fb$R_in,
            field_n_out_used_H = fb$n_out_used, field_n_in_used_mean_H = fb$n_in_used_mean,
            field_n_out_used_Hc = fc$n_out_used, mr_sec = g$timing_seconds,
            field_sec_H = fb$timing_seconds, field_sec_Hc = fc$timing_seconds,
            harm_flag = g$harm_flag)
}
rec <- list(rung = rung, identifier = spec$id, fit_sec = fit_sec, warnings = warns,
            sg.harm = fit$sg.harm, label = label, N_H = nH, E_H = eH, N_Hc = nHc, E_Hc = eHc,
            in_H = in_H,
            filter_counts = fit$find.grps$filter_counts,
            n_candidates = if (!is.null(gcs)) gcs[grep("^n_candidates_|^n_passed$|^algorithm$", names(gcs))] else NULL,
            n_hr_subgroups = if (!is.null(fit$find.grps$out.found$hr.subgroups)) nrow(fit$find.grps$out.found$hr.subgroups) else NA_integer_,
            out_sg_result = gcs$out_sg$result,
            pcons_closed = pc_closed, pcons_split = pc_split,
            n_closed_calls = if (is.null(pc_closed)) 0L else nrow(pc_closed),
            n_fallback = if (is.null(pc_split)) 0L else nrow(pc_split),
            n_closed_na = if (is.null(pc_closed)) 0L else sum(is.na(pc_closed$rate_closed)),
            mr = q, consistency_algorithm = fit$consistency_algorithm,
            family_status = fit$family_status)
saveRDS(rec, file.path("rungs", paste0("rung_", rung, ".rds")))
saveRDS(fit, file.path("fits", paste0("fit_", rung, ".rds")))

f3 <- function(x) if (is.null(x) || !length(x)) NA_character_ else sprintf("%.6f", x)
row <- data.frame(rung = rung, identifier = spec$id, fit_sec = round(fit_sec, 2), region = label,
                  N_H = nH, E_H = eH, N_Hc = nHc, E_Hc = eHc,
                  n_screened_closed = rec$n_closed_calls, n_fallback = rec$n_fallback,
                  n_family = if (length(q)) q$n_family else NA, mr_draws = if (length(q)) q$mr_draws else NA,
                  field_R_out = if (length(q)) q$field_R_out else NA, field_R_in = if (length(q)) q$field_R_in else NA,
                  naive_H = f3(q$naive_H), naive_lo_H = f3(q$naive_lo_H), naive_hi_H = f3(q$naive_hi_H),
                  corrected_H = f3(q$corrected_H), field_lower_H = f3(q$field_lower_H),
                  naive_Hc = f3(q$naive_Hc), naive_lo_Hc = f3(q$naive_lo_Hc), naive_hi_Hc = f3(q$naive_hi_Hc),
                  corrected_Hc = f3(q$corrected_Hc), field_s_upper_Hc = f3(q$field_s_upper_Hc),
                  bonf_lower_H = f3(q$bonf_lower_H), bonf_upper_Hc = f3(q$bonf_upper_Hc),
                  n_warnings = length(warns), stringsAsFactors = FALSE)
qf <- "ladder_quantities.csv"
write.table(row, qf, sep = ",", row.names = FALSE, col.names = !file.exists(qf), append = file.exists(qf))
cat("---- record ----\n"); print(t(row))
cat("filter_counts:\n"); str(rec$filter_counts)
cat("consistency counts:\n"); str(rec$n_candidates)
if (!is.null(pc_closed)) { cat("p_cons (closed form, consistency stage):\n"); print(pc_closed[order(-pc_closed$rate_closed), ], row.names = FALSE) }
if (!is.null(pc_split)) { cat("fallback splits:\n"); print(pc_split, row.names = FALSE) }

# ---- GATE A (S0 only) --------------------------------------------------------
if (rung == "S0") {
  ok <- identical(label, "{er <= 0} & {pgr <= 26}") && nH == 75L && eH == 41L
  cat(sprintf("GATE A: %s (region %s, N %d, events %d)\n", if (ok) "PASS" else "FAIL", label, nH, eH))
  if (!ok) quit(status = 3L)
}
