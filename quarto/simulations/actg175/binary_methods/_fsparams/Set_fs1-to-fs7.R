# First, remove clutter from prior use
# list.files("quarto/simulations/actg175/binary",
#            pattern = "\\.qmd(\\.bak)?$", full.names = TRUE)
# # de-clutter
# bak_files <- Sys.glob("quarto/simulations/actg175/binary/*.qmd.bak")
# bak_files                       # inspect first
# file.remove(bak_files)          # then delete

source("~/Documents/GitHub/forestsearch/quarto/_utils/gsub_in_files.R", echo = TRUE)

# Note: do not use rename = TRUE as that will over-write "A = fs2"

files <- Sys.glob("quarto/simulations/actg175/binary/*_fs1.qmd")

gsub_in_files(files, "fs1", "fs7",
                 copy = TRUE)        # PREVIEW


gsub_in_files(files, "fs1", "fs7",
             copy = TRUE,
             dry_run = FALSE)  # COMMIT
