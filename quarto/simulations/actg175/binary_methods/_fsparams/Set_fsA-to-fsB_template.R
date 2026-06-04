# First, remove clutter from prior use
# list.files("quarto/simulations/actg175/binary",
#            pattern = "\\.qmd(\\.bak)?$", full.names = TRUE)
# # de-clutter
# bak_files <- Sys.glob("quarto/simulations/actg175/binary/*.qmd.bak")
# bak_files                       # inspect first
# file.remove(bak_files)          # then delete


# Note: do not use rename = TRUE as that will over-write "A = fs2"

files <- Sys.glob("quarto/simulations/actg175/binary/*_fs2.qmd")

gsub_in_files(files, "fs2", "fs3",
                 copy = TRUE)        # PREVIEW


gsub_in_files(files, "fs2", "fs3",
             copy = TRUE,
             dry_run = FALSE)  # COMMIT
