files <- Sys.glob("quarto/simulations/actg175/binary/*_fs2.qmd")

# full to demo
# First check what is done
gsub_in_files(files,
              find    = 'sim_mode\\s*<-\\s*"full"',
              replace = 'sim_mode <- "demo"',
              fixed   = FALSE,
              dry_run = TRUE)

# Commit

gsub_in_files(files,
              find    = 'sim_mode\\s*<-\\s*"full"',
              replace = 'sim_mode <- "demo"',
              fixed   = FALSE,
              dry_run = FALSE, backup = TRUE)


# backup = TRUE allows for quick recovery
#file.copy("foo.qmd.bak", "foo.qmd", overwrite = TRUE)
#baks <- Sys.glob("quarto/simulations/*_fs2.qmd.bak")
#for (b in baks) file.copy(b, sub("\\.bak$", "", b), overwrite = TRUE)


# demo to full
# First check what is done
gsub_in_files(files,
              find    = 'sim_mode\\s*<-\\s*"demo"',
              replace = 'sim_mode <- "full"',
              fixed   = FALSE,
              dry_run = TRUE)

# Commit

gsub_in_files(files,
              find    = 'sim_mode\\s*<-\\s*"demo"',
              replace = 'sim_mode <- "full"',
              fixed   = FALSE,
              dry_run = FALSE, backup = TRUE)
