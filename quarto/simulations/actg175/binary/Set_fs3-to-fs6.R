
files <- Sys.glob("quarto/simulations/actg175/binary/*_fs3.qmd")

gsub_in_files(files, "fs3", "fs6",
                 rename = TRUE)        # PREVIEW


gsub_in_files(files, "fs3", "fs6",
             rename = TRUE,
             dry_run = FALSE, backup = TRUE)  # COMMIT
