
files <- Sys.glob("quarto/simulations/actg175/binary/*_fs1.qmd")

gsub_in_files(files, "fs1", "fs2",
                 copy = TRUE)        # PREVIEW


gsub_in_files(files, "fs1", "fs2",
             dry_run = FALSE, copy = TRUE)  # COMMIT
