#!/usr/bin/env bash
# strip.sh -- Part R step 5 of TASK_classification_outside_package_2026-11-11_v2.
# The precedent's method (dev/reports/articles_finish_2026-10-30/strip_and_gate.sh), for one article:
# results/paper_simulations.txt is pandoc plain text of the whole pkgdown page; every <pre ...>...</pre>
# element is deleted (R gsub, perl, dot-all "(?s)<pre[ >].*?</pre>") and the rest converted with
# pandoc --wrap=none -f html -t plain -> results/paper_simulations_output.txt. The precedent's
# paper_simulations gates (5, 6, 12a, 12b) are repeated in results/keyvalues.txt. Run from the repo root.
set -u
SRC=dev/reports/paper_sims_rerun_2026-11-11
R=$SRC/results
KV=$R/keyvalues.txt
: > "$KV"
n=paper_simulations
pandoc --wrap=none -t plain "$SRC/$n.html" > "$R/$n.txt"
Rscript -e 'a <- commandArgs(TRUE); h <- paste(readLines(a[1], warn = FALSE), collapse = "\n");
  n <- lengths(regmatches(h, gregexpr("<pre[ >]", h)));
  h2 <- gsub("(?s)<pre[ >].*?</pre>", "", h, perl = TRUE);
  stopifnot(!grepl("<pre[ >]", h2)); writeLines(h2, a[2]); cat(a[3], "pre_elements_removed =", n, "\n")' \
  "$SRC/$n.html" "$R/${n}_nopre.html" "$n" | tee -a "$KV"
pandoc --wrap=none -f html -t plain "$R/${n}_nopre.html" > "$R/${n}_output.txt"
rm "$R/${n}_nopre.html"
g() { local c; c=$(grep -oF -- "$2" "$R/${n}_output.txt" | wc -l); echo "$1 [$n] '$2' = $c ; expect $3 ; $(awk -v n="$c" "BEGIN{print ($3) ? \"PASS\" : \"FAIL\"}")" | tee -a "$KV"; }
g 5   "44 simulated trials under the alternative and 44 under the null" "n==1"
g 6   "Identifying Patient Subgroups" "n==0"
g 12a "NaN" "n==0"
g 12b "?@" "n==0"
