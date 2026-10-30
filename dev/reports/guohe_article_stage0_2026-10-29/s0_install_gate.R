# Stage 0 v2 sections 1.6-1.7: installed-build gate (fresh Rscript, repo root)
kv_file <- "dev/reports/guohe_article_stage0_2026-10-29/results/keyvalues.txt"
kv <- readLines(kv_file)
start <- trimws(sub("^task_start_utc = ", "", grep("^task_start_utc", kv, value = TRUE)))
pd <- utils::packageDescription("forestsearch")
built_dt <- trimws(strsplit(pd$Built, ";")[[1]][3])
lib <- dirname(dirname(dirname(attr(pd, "file"))))
cores <- parallelly::availableCores()
si <- Sys.info()[c("sysname", "machine")]
t_built <- as.POSIXct(sub(" UTC$", "", built_dt), tz = "UTC")
t_start <- as.POSIXct(start, tz = "UTC")
g_version <- identical(pd$Version, "0.4.1")
g_built <- isTRUE(t_built > t_start)
g_exp <- all(c("fs_mr_inference", "guohe_algorithm3") %in% getNamespaceExports("forestsearch"))
out <- c(
  sprintf("pd_Version = %s", pd$Version),
  sprintf("pd_Built = %s", pd$Built),
  sprintf("pd_Packaged = %s", if (is.null(pd$Packaged)) "<absent>" else pd$Packaged),
  sprintf("library_path = %s", lib),
  sprintf("survival_version = %s", format(packageVersion("survival"))),
  sprintf("R_version = %s", R.version.string),
  sprintf("sysname = %s", si[["sysname"]]),
  sprintf("machine = %s", si[["machine"]]),
  sprintf("availableCores = %d", as.integer(cores)),
  sprintf("gate_1_6_version_0.4.1 = %s", g_version),
  sprintf("gate_1_6_built_after_start = %s (Built %s vs start %s)", g_built, built_dt, start),
  sprintf("gate_1_7_exports = %s", g_exp)
)
write(out, kv_file, append = TRUE)
cat(out, sep = "\n")
if (!(g_version && g_built && g_exp)) quit(status = 2)
