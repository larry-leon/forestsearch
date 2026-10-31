# ---- transplant (verbatim chunk bodies; fences removed) ----
# source: forestsearch/quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd
# forestsearch commit 43b051b6 (2026-08-16 12:07; the render that wrote the payload at 2026-08-16 16:50:40 UTC ran this text — its HTML is committed in the same commit; identical to HEAD; differs from fd1ce74f, the last commit before the build, in fs_selection_rule both->neighborhood)
#   chunk 'export-section6-2-payload': fence lines 1929-2124
# ----------------------------------------------------------

## ==== chunk 'export-section6-2-payload' ====
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

# Adverse-outcome rate by arm inside an engine's harm region (treat.recommend == 0).
.harm_rates <- function(fit, oname, tname) {
  if (is.null(fit) || is.na(oname) || is.na(tname)) return(NULL)
  d <- fit$df.est
  if (is.null(d) || !all(c(oname, tname, "treat.recommend") %in% names(d))) return(NULL)
  h <- d[d$treat.recommend == 0, , drop = FALSE]
  if (nrow(h) == 0L) return(NULL)
  yt <- h[[oname]][h[[tname]] == 1]; yc <- h[[oname]][h[[tname]] == 0]
  list(treated_pct = if (length(yt)) 100 * mean(yt) else NA_real_,
       control_pct = if (length(yc)) 100 * mean(yc) else NA_real_)
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
    covariates = if (exists("confounders.name")) confounders.name else NA_character_,
    timings    = list(fs = .meta_fs, dina = .meta_dina, grf = .meta_grf)),
  extras = list(
    concordance = .concordance(.fs, .fs_grf, .idname),
    harm_rates  = .harm_rates(.fs, .oname, .tname)),
  est_scale = "or",
  built_at  = Sys.time(),
  forestsearch_version = tryCatch(as.character(utils::packageVersion("forestsearch")), error = function(e) NA_character_))

.payload_file <- file.path(.out_dir, paste0(.dirout, "_payload.rds"))
saveRDS(payload, .payload_file)
cat("Wrote Section 6.2 payload to:", .payload_file, "\n")

