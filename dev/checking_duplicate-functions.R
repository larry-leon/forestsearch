setwd("dev/R_revised")

# Use double backslashes for regex in R strings
extract_function_names <- function(file) {
  lines <- readLines(file, warn = FALSE)
  # Find lines that look like function definitions
  fun_lines <- grep("^[a-zA-Z0-9_.]+\\s*<-\\s*function", lines, value = TRUE)
  # Extract function names
  fun_names <- sub("\\s*<-.*", "", fun_lines)
  fun_names
}

# List of uploaded files
files <- c(
  "forestsearch_cross-validation.R",
  "crossvalidation_helpers.R",
  "forestsearch_bootstrap_dofuture.R",
  "grf_functions.R",
  "subgroup_consistency.R",
  "subgroup_search.R",
  "forest_search.R",
  "get_FSdata.R",
  "forestsearch_bootstrap_helpers.R",
  "get_FSdata_helpers.R",
  "summary_utility_functions.R"
)

# Extract function names from each file
fun_list <- lapply(files, extract_function_names)
names(fun_list) <- files

# Combine all function names with file info
df_fun <- data.frame(
  file = rep(names(fun_list), lengths(fun_list)),
  fun = unlist(fun_list, use.names = FALSE),
  stringsAsFactors = FALSE
)

# Find duplicates
fun_dups <- df_fun$fun[duplicated(df_fun$fun) | duplicated(df_fun$fun, fromLast = TRUE)]
df_dups <- df_fun[df_fun$fun %in% fun_dups, ]

nrow(df_dups)



