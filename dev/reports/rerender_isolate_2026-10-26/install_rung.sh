#!/usr/bin/env bash
# usage: install_rung.sh <rung> <source>
#   <source> = a directory (worktree or the main tree) -> devtools::install(dependencies = FALSE, upgrade = FALSE, quick = TRUE)
#            = "CRAN"                                   -> install.packages("forestsearch", repos = cloud.r-project.org)
# Appends to install_log.txt; exit status is the install's.
REC=/home/larryleon/Documents/GitHub/forestsearch/dev/reports/rerender_isolate_2026-10-26
rung=$1; src=$2
{
  echo "=== install rung=$rung source=$src  $(date -Is)"
  if [ "$src" = "CRAN" ]; then
    Rscript -e 'install.packages("forestsearch", repos = "https://cloud.r-project.org", dependencies = FALSE)' 2>&1
  else
    echo "source SHA: $(git -C "$src" rev-parse HEAD)"
    Rscript -e 'devtools::install(commandArgs(TRUE)[1], dependencies = FALSE, upgrade = FALSE, quick = TRUE)' "$src" 2>&1
  fi
  rc=$?
  echo "install exit=$rc"
  Rscript -e 'print(packageDescription("forestsearch")[c("Version","Packaged","Built","Repository")]); cat("lean fn exists:", exists(".fs_logit_two_group", asNamespace("forestsearch")), "\n")' 2>&1
  exit $rc
} >> "$REC/install_log.txt" 2>&1
