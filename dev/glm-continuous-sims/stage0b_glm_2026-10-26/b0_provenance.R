# b0_provenance.R -- provenance record and gates for Stage 0b.
#   GATE 1  installed forestsearch is 0.4.0 (exit 2)
#   GATE 2  speff2trial (ACTG175) is installed (exit 3)
#   RECORD  repo files under study vs CRAN 0.4.0 MD5; fixture sources

source("stage0b_common.R")
unlink(out_path("keyvalues.txt"))

pd <- utils::packageDescription("forestsearch")
kv(installed_version    = as.character(pd$Version),
   installed_packaged   = if (is.null(pd$Packaged)) "<absent>" else pd$Packaged,
   installed_repository = if (is.null(pd$Repository)) "<absent>" else pd$Repository,
   installed_libpath    = dirname(dirname(dirname(attr(pd, "file")))),
   survival_version     = as.character(utils::packageVersion("survival")),
   future_version       = as.character(utils::packageVersion("future")),
   data_table_version   = as.character(utils::packageVersion("data.table")),
   r_version            = R.version.string,
   sysname              = unname(Sys.info()[["sysname"]]),
   available_cores      = parallelly::availableCores(),
   speff2trial_installed = requireNamespace("speff2trial", quietly = TRUE))
if (as.character(pd$Version) != "0.4.0") {
  cat("GATE 1 FAILED: installed forestsearch is", as.character(pd$Version), "- this audit is based on CRAN 0.4.0.\n")
  quit(status = 2)
}
if (!requireNamespace("speff2trial", quietly = TRUE)) {
  cat("GATE 2 FAILED: speff2trial (ACTG175) is not installed.\n")
  quit(status = 3)
}

cran_md5 <- c(
  "DESCRIPTION"                          = "19b9459a03ed5c4a0c12943eb57d82b6",
  "NAMESPACE"                            = "00257d44d6f21d10f9af68b9f0dc227c",
  "R/Cox_estimation_helpers.R"           = "ffa9744bf4b3264041aa4c787cd5958f",
  "R/consistency_resample.R"             = "861c6b4e7d00c8176bbab6b7f863f724",
  "R/forestsearch_main.R"                = "11d80e74c4fcdba43fe8bbf9a1287d3f",
  "R/glm_effect_estimators.R"            = "8e8239f94e5f832af33e4d9cc0ff64d3",
  "R/subgroup_consistency_helpers.R"     = "f32b55b62d683337f6e1b36b3aa3bc17",
  "R/subgroup_consistency_main.R"        = "1285150ed4891883468d01d5ee981b78",
  "R/subgroup_search.R"                  = "f06d6bba8135f41a57962a659ebbe601",
  "R/bootstrap_analysis_dofuture.R"      = "a891d66e7337b9a1ebdf524463102978",
  "R/forestsearch_cross_validation.R"    = "359955f9ba9bd850c978230aee3b4429"
)
root <- tryCatch(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE, stderr = FALSE),
                 error = function(e) character(0), warning = function(w) character(0))
if (length(root) != 1L || !nzchar(root)) {
  kv(repo_root = "<not a git repo>")
} else {
  kv(repo_root = root,
     repo_head = system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE),
     repo_branch = system2("git", c("-C", shQuote(root), "rev-parse", "--abbrev-ref", "HEAD"), stdout = TRUE))
  files <- file.path(root, names(cran_md5)); have <- file.exists(files)
  got <- ifelse(have, unname(tools::md5sum(files)), NA_character_)
  tab <- data.frame(file = names(cran_md5), cran_0_4_0 = unname(cran_md5), repo_head = got,
                    match = ifelse(is.na(got), "MISSING", ifelse(got == unname(cran_md5), "match", "DIFFERS")),
                    stringsAsFactors = FALSE)
  utils::write.csv(tab, out_path("provenance_md5.csv"), row.names = FALSE)
  print(tab[, c("file", "match")], row.names = FALSE)
  kv(repo_files_matching_cran = sum(tab$match == "match"), repo_files_differing = sum(tab$match != "match"))
  for (f in tab$file[tab$match != "match"]) {
    cat("--", f, "\n")
    lg <- system2("git", c("-C", shQuote(root), "log", "--oneline", "--since=2026-10-09", "--", shQuote(f)), stdout = TRUE)
    cat(if (length(lg)) paste0("   ", lg) else "   (no commits since 2026-10-09)", sep = "\n")
  }
}

# fixture sources (recorded; the task document says what they must be)
for (nm in names(FIXTURES)) {
  fx <- FIXTURES[[nm]]()
  do.call(kv, setNames(list(fx$label, fx$source, nrow(fx$df), paste(names(fx$args), collapse = ",")),
                       paste0("fixture_", nm, c("_label", "_source", "_n", "_argnames"))))
}
cat("\nb0 complete\n")
