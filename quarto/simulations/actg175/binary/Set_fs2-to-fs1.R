
files <- Sys.glob("quarto/simulations/actg175/binary/*_fs2.qmd")

gsub_in_files(files, "fs2", "fs1",
                 rename = TRUE)        # PREVIEW


gsub_in_files(files, "fs2", "fs1",
             rename = TRUE,
             dry_run = FALSE, backup = TRUE)  # COMMIT
