# ---------------------------------------------------------------------------
# Table S2 of León & Anderson (2026a), table assembly (naive / full bootstrap / multiplier) -- verbatim transplant.
# Source: forestsearch quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd
# Commit: cf4d64328b6288ce2ec009f2c2fb15689de6a4cb (2026-08-15 21:12 -0700), the last
#   commit touching the .qmd before the Table S2 payload (built 2026-08-15 22:50:30 PDT,
#   Mac-Studio-3.local, forestsearch 0.2.0); render html committed in 43b051b6.
#   At HEAD the .qmd differs from cf4d6432 only at lines 144-145 (run_cv / run_loo
#   FALSE -> TRUE, 6ad5a724); no transplanted chunk below is affected.
# Chunks (name:fence lines, inclusive): bootstrap-summary:967-978 debias-MR-table:1089-1091 fs-dina-bootstrap-summary:1610-1619 fs-dina-debias-MR-table:1706-1708 fs-grf-bootstrap-summary:1843-1852 fs-grf-debias-MR-table:1942-1944 export-section6-2-payload:2434-2641 (the payload's $table is what Table S2 renders)
# Chunk bodies are copied verbatim; only these header comments and the
# '## ---- chunk: <name> ----' separators were added. Transplanted 2026-10-28.
# ---------------------------------------------------------------------------

## ---- chunk: bootstrap-summary (lines 968-977) ----
# Comprehensive summary with diagnostics
summaries <- summarize_bootstrap_results(
  sgharm = fs$sg.harm,
  boot_results = fs_bc,
  create_plots = TRUE,
  est.scale = "hr"
)

# Display bias-corrected estimates table
summaries$table

## ---- chunk: debias-MR-table (lines 1090-1090) ----
mr_estimates_table(fs, boot_results = fs_bc, est.scale = "hr", digits = 2)

## ---- chunk: fs-dina-bootstrap-summary (lines 1611-1618) ----
dina_summaries <- summarize_bootstrap_results(
  sgharm = fs_dina$sg.harm,
  boot_results = fs_dina_bc,
  create_plots = FALSE,
  est.scale = "hr"
)

dina_summaries$table

## ---- chunk: fs-dina-debias-MR-table (lines 1707-1707) ----
mr_estimates_table(fs_dina, boot_results = fs_dina_bc, est.scale = "hr", digits = 2)

## ---- chunk: fs-grf-bootstrap-summary (lines 1844-1851) ----
grf_summaries <- summarize_bootstrap_results(
  sgharm = fs_grf$sg.harm,
  boot_results = fs_grf_bc,
  create_plots = FALSE,
  est.scale = "hr"
)

grf_summaries$table

## ---- chunk: fs-grf-debias-MR-table (lines 1943-1943) ----
mr_estimates_table(fs_grf, boot_results = fs_grf_bc, est.scale = "hr", digits = 2)

## ---- chunk: export-section6-2-payload (lines 2435-2640) ----
#| echo: false
#| message: false

.qmd_dir  <- tryCatch(dirname(knitr::current_input(dir = TRUE)),
                      error = function(e) getwd())
.qmd_stem <- tryCatch(tools::file_path_sans_ext(basename(knitr::current_input())),
                      error = function(e) "analysis")
.dirout   <- if (is.null(dirout)) .qmd_stem else dirout
.base     <- if (is.null(results_dir)) file.path(.qmd_dir, "_payloads") else results_dir
.out_dir  <- file.path(.base, .dirout)
dir.create(.out_dir, recursive = TRUE, showWarnings = FALSE)

.get <- function(nm) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else NULL

.triple <- function(x) {
  if (is.null(x) || is.null(x$est)) c(NA_real_, NA_real_, NA_real_) else c(x$est, x$lower, x$upper)
}
.fb_triple <- function(est_list) {
  if (is.null(est_list)) c(NA_real_, NA_real_, NA_real_) else c(est_list$H2, est_list$H2_lower, est_list$H2_upper)
}

.engine_rows <- function(method, fit, boot) {
  if (is.null(fit) || is.null(fit$mr_inference)) return(NULL)
  g  <- fit$mr_inference
  cc <- g$complement
  tr <- fit$df.est$treat.recommend
  n_tot <- length(tr); n_H <- sum(tr == 0, na.rm = TRUE)
  nv_H  <- .triple(g$naive)
  nv_Hc <- .triple(if (!is.null(cc)) cc$naive)
  mr_H  <- .triple(g$debiased)
  mr_Hc <- .triple(if (!is.null(cc)) cc$debiased)
  fb_H  <- .fb_triple(if (!is.null(boot)) boot$H_estimates)
  fb_Hc <- .fb_triple(if (!is.null(boot)) boot$Hc_estimates)
  data.frame(
    method = method, region = c("H", "Hc"),
    n = c(n_H, n_tot - n_H),
    pct = round(100 * c(n_H, n_tot - n_H) / n_tot),
    naive_est = c(nv_H[1], nv_Hc[1]), naive_lo = c(nv_H[2], nv_Hc[2]), naive_hi = c(nv_H[3], nv_Hc[3]),
    fb_est = c(fb_H[1], fb_Hc[1]), fb_lo = c(fb_H[2], fb_Hc[2]), fb_hi = c(fb_H[3], fb_Hc[3]),
    mr_est = c(mr_H[1], mr_Hc[1]), mr_lo = c(mr_H[2], mr_Hc[2]), mr_hi = c(mr_H[3], mr_Hc[3]),
    stringsAsFactors = FALSE)
}

.engine_meta <- function(fit, boot_seconds) {
  if (is.null(fit) || is.null(fit$mr_inference))
    return(list(label = NA_character_, n_family = NA_integer_, selection_bias = NA_real_,
                boot_sec = NA_real_, mr_sec = NA_real_, speedup = NA_real_))
  g <- fit$mr_inference
  mr_sec <- if (!is.null(g$timing_seconds)) as.numeric(g$timing_seconds) else NA_real_
  list(label = paste(fit$sg.harm, collapse = " & "),
       n_family = g$n_family, selection_bias = g$selection_bias,
       boot_sec = boot_seconds, mr_sec = mr_sec,
       speedup = if (is.finite(boot_seconds) && is.finite(mr_sec) && mr_sec > 0) boot_seconds / mr_sec else NA_real_)
}

.tm <- .get("timings")
.boot_sec <- function(key) {
  v <- if (!is.null(.tm) && !is.null(.tm[[key]])) suppressWarnings(as.numeric(.tm[[key]])) else NA_real_
  if (length(v) == 1L && is.finite(v)) v else NA_real_
}

# FS vs GRF agreement on the harm flag (treat.recommend == 0), aligned by id, with Cohen's kappa.
.concordance <- function(a, b, idn) {
  if (is.null(a) || is.null(b)) return(NULL)
  da <- a$df.est; db <- b$df.est
  if (is.null(da) || is.null(db) || !(idn %in% names(da)) || !(idn %in% names(db))) return(NULL)
  m <- merge(data.frame(.id = da[[idn]], ha = as.integer(da$treat.recommend == 0)),
             data.frame(.id = db[[idn]], hb = as.integer(db$treat.recommend == 0)), by = ".id")
  if (nrow(m) == 0L) return(NULL)
  po <- mean(m$ha == m$hb); pa <- mean(m$ha); pb <- mean(m$hb)
  pe <- pa * pb + (1 - pa) * (1 - pb)
  list(n_common = nrow(m), n_agree = sum(m$ha == m$hb), pct_agree = 100 * po,
       kappa = if (pe < 1) (po - pe) / (1 - pe) else NA_real_)
}

# Observed event rate by arm inside an engine's harm region (treat.recommend == 0).
# Cox analog of the binary version: the "rate" is the proportion experiencing the
# event (status == 1), reported alongside n and event counts so a reader can
# recompute incidence later if person-time is added.
.harm_rates <- function(fit, ename, tname) {
  if (is.null(fit) || is.na(ename) || is.na(tname)) return(NULL)
  d <- fit$df.est
  if (is.null(d) || !all(c(ename, tname, "treat.recommend") %in% names(d))) return(NULL)
  h <- d[d$treat.recommend == 0, , drop = FALSE]
  if (nrow(h) == 0L) return(NULL)
  et <- h[[ename]][h[[tname]] == 1]; ec <- h[[ename]][h[[tname]] == 0]
  list(treated_n      = length(et),
       control_n      = length(ec),
       treated_events = if (length(et)) sum(et == 1, na.rm = TRUE) else NA_integer_,
       control_events = if (length(ec)) sum(ec == 1, na.rm = TRUE) else NA_integer_,
       treated_pct    = if (length(et)) 100 * mean(et == 1, na.rm = TRUE) else NA_real_,
       control_pct    = if (length(ec)) 100 * mean(ec == 1, na.rm = TRUE) else NA_real_)
}

# Machine / hardware fingerprint of whatever computer ran this render, so that
# timings pooled across different machines stay attributable. Pure base R plus
# guarded platform calls (Linux / macOS / Windows); any field that cannot be
# determined returns NA rather than erroring.
.machine_info <- function(n_workers = NULL) {
  si   <- Sys.info()
  sn   <- unname(si[["sysname"]])               # "Linux" | "Darwin" | "Windows"
  sess <- tryCatch(sessionInfo(), error = function(e) NULL)

  .sys <- function(cmd) tryCatch(
    trimws(paste(system(cmd, intern = TRUE, ignore.stderr = TRUE), collapse = " ")),
    error = function(e) NA_character_)

  cpu_model <- tryCatch({
    if (sn == "Linux") {
      ln <- grep("^model name", readLines("/proc/cpuinfo", warn = FALSE), value = TRUE)[1]
      if (length(ln) && !is.na(ln)) trimws(sub("^model name\\s*:\\s*", "", ln)) else NA_character_
    } else if (sn == "Darwin") {
      .sys("sysctl -n machdep.cpu.brand_string")
    } else if (sn == "Windows") {
      out <- strsplit(.sys("wmic cpu get name"), "  +")[[1]]
      out <- out[nzchar(out) & out != "Name"]
      if (length(out)) out[length(out)] else NA_character_
    } else NA_character_
  }, error = function(e) NA_character_)

  ram_gb <- tryCatch({
    if (sn == "Linux") {
      kb <- as.numeric(sub("\\D*(\\d+).*", "\\1",
              grep("^MemTotal", readLines("/proc/meminfo", warn = FALSE), value = TRUE)[1]))
      round(kb / 1024^2, 1)
    } else if (sn == "Darwin") {
      round(as.numeric(.sys("sysctl -n hw.memsize")) / 1024^3, 1)
    } else if (sn == "Windows") {
      round(as.numeric(gsub("\\D", "", .sys("wmic ComputerSystem get TotalPhysicalMemory"))) / 1024^3, 1)
    } else NA_real_
  }, error = function(e) NA_real_)

  list(
    hostname           = unname(si[["nodename"]]),
    user               = unname(si[["user"]]),
    os                 = if (!is.null(sess) && !is.null(sess$running)) sess$running else NA_character_,
    sysname            = sn,
    os_release         = unname(si[["release"]]),
    arch               = unname(si[["machine"]]),
    cpu_model          = if (length(cpu_model) == 1L) cpu_model else NA_character_,
    cores_physical     = tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_),
    cores_logical      = tryCatch(parallel::detectCores(logical = TRUE),  error = function(e) NA_integer_),
    cores_used         = if (!is.null(n_workers)) as.integer(n_workers) else NA_integer_,
    workers_registered = tryCatch(future::nbrOfWorkers(), error = function(e) NA_integer_),
    future_plan        = tryCatch(class(future::plan())[1], error = function(e) NA_character_),
    ram_gb             = if (length(ram_gb) == 1L) ram_gb else NA_real_,
    r_version          = R.version.string,
    r_platform         = R.version$platform,
    blas               = if (!is.null(sess) && !is.null(sess$BLAS)) sess$BLAS else NA_character_,
    lapack             = if (!is.null(sess) && !is.null(sess$LAPACK)) sess$LAPACK else NA_character_,
    timestamp          = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
}

.fs      <- .get("fs");      .fs_bc      <- .get("fs_bc")
.fs_dina <- .get("fs_dina"); .fs_dina_bc <- .get("fs_dina_bc")
.fs_grf  <- .get("fs_grf");  .fs_grf_bc  <- .get("fs_grf_bc")

.meta_fs   <- .engine_meta(.fs,      .boot_sec("bootstrap"))
.meta_dina <- .engine_meta(.fs_dina, .boot_sec("fs_dina_bootstrap"))
.meta_grf  <- .engine_meta(.fs_grf,  .boot_sec("fs_grf_bootstrap"))

table1 <- do.call(rbind, list(
  .engine_rows("FS",   .fs,      .fs_bc),
  .engine_rows("DINA", .fs_dina, .fs_dina_bc),
  .engine_rows("GRF",  .fs_grf,  .fs_grf_bc)))

.dsel <- .get("dina_sel")
dina_standalone <- if (!is.null(.dsel) && isTRUE(.dsel$found)) {
  op <- if (identical(.dsel$direction, "left")) "<=" else ">="
  list(found = TRUE, label = sprintf("%s %s %s", .dsel$covariate, op, format(.dsel$threshold)), n = .dsel$n_subgroup)
} else {
  list(found = FALSE, label = NA_character_, n = NA_integer_)
}

.idname <- if (exists("id.name")) id.name else "id"
.oname  <- if (exists("outcome.name")) outcome.name else NA_character_
.ename  <- if (exists("event.name")) event.name else NA_character_
.tname  <- if (exists("treat.name")) treat.name else NA_character_

payload <- list(
  table  = table1,
  labels = list(fs = .meta_fs$label, dina = .meta_dina$label, grf = .meta_grf$label,
                dina_standalone = dina_standalone),
  meta   = list(
    B          = if (exists("NB")) NB else NA_integer_,
    mr_draws = if (exists("mr_draws")) mr_draws else NA_integer_,
    n_total    = if (!is.null(.fs)) nrow(.fs$df.est) else NA_integer_,
    cores      = tryCatch(unname(future::availableCores()), error = function(e) NA_integer_),
    machine    = .machine_info(n_workers = if (exists("n_cores")) n_cores else NULL),
    outcome    = .oname,
    event      = .ename,
    covariates = if (exists("confounders.name")) confounders.name else NA_character_,
    timings    = list(fs = .meta_fs, dina = .meta_dina, grf = .meta_grf)),
  extras = list(
    concordance = .concordance(.fs, .fs_grf, .idname),
    harm_rates  = .harm_rates(.fs, .ename, .tname)),
  est_scale = "hr",
  built_at  = Sys.time(),
  forestsearch_version = tryCatch(as.character(utils::packageVersion("forestsearch")), error = function(e) NA_character_))

# Payload filename -- analog of actg175's "actg175_table_payload.rds".
# (If your gbsg reader/manuscript doc expects "gbsg_table1_payload.rds", change it here.)
.payload_file <- file.path(.out_dir, paste0(.dirout, "_payload.rds"))
saveRDS(payload, .payload_file)
cat("Wrote Section 6.2 payload to:", .payload_file, "\n")

