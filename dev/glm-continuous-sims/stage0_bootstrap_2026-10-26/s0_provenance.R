# s0_provenance.R -- provenance record and the two gates that must hold before
# anything is measured.  Exits non-zero on a failed gate (the driver stops).
#
# GATE 1  installed forestsearch is version 0.4.0 (the CRAN release this audit
#         is based on).
# GATE 2  multicore is usable in this session (candidate C cannot be measured
#         otherwise; relaunch from a plain terminal).
# RECORD  repo R/ files under study vs their CRAN 0.4.0 MD5 sums (from the
#         CRAN tarball's MD5 file).  A mismatch is NOT a stop: the measurement
#         runs against the installed package, and GATE 1 covers that.  Every
#         mismatch is listed with its git log since the release so the report
#         can say whether the survey's line references still hold.

source("stage0_common.R")
unlink(out_path("keyvalues.txt"))

# --- installed package -------------------------------------------------------
pd <- utils::packageDescription("forestsearch")
kv(installed_version   = as.character(pd$Version),
   installed_packaged  = if (is.null(pd$Packaged)) "<absent>" else pd$Packaged,
   installed_built     = if (is.null(pd$Built)) "<absent>" else pd$Built,
   installed_repository = if (is.null(pd$Repository)) "<absent>" else pd$Repository,
   installed_remotetype = if (is.null(pd$RemoteType)) "<absent>" else pd$RemoteType,
   installed_remotesha  = if (is.null(pd$RemoteSha)) "<absent>" else pd$RemoteSha,
   installed_libpath   = dirname(dirname(dirname(attr(pd, "file")))),
   survival_version    = as.character(utils::packageVersion("survival")),
   future_version      = as.character(utils::packageVersion("future")),
   doFuture_version    = as.character(utils::packageVersion("doFuture")),
   r_version           = R.version.string,
   sysname             = unname(Sys.info()[["sysname"]]),
   machine             = unname(Sys.info()[["machine"]]),
   available_cores     = parallelly::availableCores(),
   supports_multicore  = parallelly::supportsMulticore(),
   in_rstudio          = nzchar(Sys.getenv("RSTUDIO")))

if (as.character(pd$Version) != "0.4.0") {
  cat("GATE 1 FAILED: installed forestsearch is", as.character(pd$Version),
      "- this audit is based on CRAN 0.4.0.\n")
  quit(status = 2)
}
if (!isTRUE(parallelly::supportsMulticore())) {
  cat("GATE 2 FAILED: multicore futures are not supported in this session",
      "(RStudio or Windows).  Relaunch run_stage0.sh from a plain terminal.\n")
  quit(status = 3)
}

# --- repo files vs CRAN 0.4.0 ------------------------------------------------
cran_md5 <- c(
  "DESCRIPTION"                           = "19b9459a03ed5c4a0c12943eb57d82b6",
  "NAMESPACE"                             = "00257d44d6f21d10f9af68b9f0dc227c",
  "R/Cox_estimation_helpers.R"            = "ffa9744bf4b3264041aa4c787cd5958f",
  "R/bootstrap_analysis_dofuture.R"       = "a891d66e7337b9a1ebdf524463102978",
  "R/bootstrap_calculations_helpers.R"    = "2c25b9b63263653c599f7143df6b854c",
  "R/bootstrap_dofuture_main.R"           = "65cb60f5b23303578362172bc3cb9291",
  "R/bootstrap_dofuture_setup_helpers.R"  = "97cac83ad17bff5bf16efe7f437e164a",
  "R/consistency_resample.R"              = "861c6b4e7d00c8176bbab6b7f863f724",
  "R/forestsearch_main.R"                 = "11d80e74c4fcdba43fe8bbf9a1287d3f",
  "R/get_fsdata.R"                        = "a7f71d659558706ae121d72e82099181",
  "R/guohe_algorithm3.R"                  = "5a6e842f871c61ddef41dbfcc93bf6c5",
  "R/subgroup_consistency_helpers.R"      = "f32b55b62d683337f6e1b36b3aa3bc17",
  "R/subgroup_consistency_main.R"         = "1285150ed4891883468d01d5ee981b78",
  "R/subgroup_search.R"                   = "f06d6bba8135f41a57962a659ebbe601"
)
root <- tryCatch(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE, stderr = FALSE),
                 error = function(e) character(0), warning = function(w) character(0))
if (length(root) != 1L || !nzchar(root)) {
  cat("not inside a git repository; repo-vs-CRAN comparison skipped\n")
  kv(repo_root = "<not a git repo>")
} else {
  head_sha <- system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE)
  branch   <- system2("git", c("-C", shQuote(root), "rev-parse", "--abbrev-ref", "HEAD"), stdout = TRUE)
  kv(repo_root = root, repo_head = head_sha, repo_branch = branch)
  files <- file.path(root, names(cran_md5))
  have  <- file.exists(files)
  got   <- ifelse(have, unname(tools::md5sum(files)), NA_character_)
  tab <- data.frame(file = names(cran_md5), cran_0_4_0 = unname(cran_md5),
                    repo_head = got,
                    match = ifelse(is.na(got), "MISSING", ifelse(got == unname(cran_md5), "match", "DIFFERS")),
                    stringsAsFactors = FALSE)
  utils::write.csv(tab, out_path("provenance_md5.csv"), row.names = FALSE)
  print(tab[, c("file", "match")], row.names = FALSE)
  kv(repo_files_matching_cran = sum(tab$match == "match"),
     repo_files_differing     = sum(tab$match != "match"))
  diff_files <- tab$file[tab$match != "match"]
  if (length(diff_files)) {
    cat("\nFiles differing from CRAN 0.4.0 -- git log since the release (2026-10-09):\n")
    for (f in diff_files) {
      cat("--", f, "\n")
      lg <- system2("git", c("-C", shQuote(root), "log", "--oneline", "--since=2026-10-09", "--", shQuote(f)),
                    stdout = TRUE)
      cat(if (length(lg)) paste0("   ", lg) else "   (no commits since 2026-10-09; file predates or is untracked)", sep = "\n")
    }
    cat("NOTE: measurements run against the INSTALLED package (GATE 1 = 0.4.0).",
        "If the installed build is not the CRAN build (installed_repository != CRAN),",
        "say so in the report's provenance section.\n")
  }
}
cat("\ns0 complete\n")
