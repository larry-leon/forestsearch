#!/usr/bin/env bash
# strip_and_gate.sh -- task sec. 1.7 text gates, output only. Run from the repo root.
# Adapted from dev/reports/articles_alignment_2026-10-30/strip_and_gate.sh: same strip method.
# Every <pre ...>...</pre> element is deleted from the pkgdown build docs/articles/<name>.html
# (R gsub, perl, dot-all "(?s)<pre[ >].*?</pre>"), then pandoc --wrap=none -f html -t plain
# -> results/<name>_output.txt. results/<name>.txt is the whole page, for the record.
# Prints one line per gate (also to results/keyvalues.txt); exit 1 on any FAIL.
set -u
SRC=docs/articles
R=dev/reports/articles_finish_2026-10-30/results
KV=$R/keyvalues.txt
: > "$KV"
fail=0
strip() {
  local n=$1
  pandoc --wrap=none -t plain "$SRC/$n.html" > "$R/$n.txt"
  Rscript -e 'a <- commandArgs(TRUE); h <- paste(readLines(a[1], warn = FALSE), collapse = "\n");
    n <- lengths(regmatches(h, gregexpr("<pre[ >]", h)));
    h2 <- gsub("(?s)<pre[ >].*?</pre>", "", h, perl = TRUE);
    stopifnot(!grepl("<pre[ >]", h2)); writeLines(h2, a[2]); cat(a[3], "pre_elements_removed =", n, "\n")' \
    "$SRC/$n.html" "$R/${n}_nopre.html" "$n" | tee -a "$KV" || { echo "STRIP FAILED $n"; exit 2; }
  pandoc --wrap=none -f html -t plain "$R/${n}_nopre.html" > "$R/${n}_output.txt"
  rm "$R/${n}_nopre.html"
}
# g id article pattern condition   (fixed strings)
g() {
  local c ok
  c=$(grep -oF -- "$3" "$R/$2_output.txt" | wc -l)
  ok=$(awk -v n="$c" "BEGIN{print ($4) ? \"PASS\" : \"FAIL\"}")
  [ "$ok" = FAIL ] && fail=1
  echo "$1 [$2] '$3' = $c ; expect $4 ; $ok" | tee -a "$KV"
}
ALL="methodology paper_simulations treatment_effect_definitions biomarker_effects causal_effects_brief_review survival-post-selection"
for n in $ALL; do strip $n; done
g 1   methodology "2026 note" "n==13"
g 2a  methodology "(eq-power?)" "n==0"
g 2b  methodology "?@eq-power" "n==0"
g 3   methodology "Equation" "n>=1"
g 4   methodology "These are the 2024 guidelines" "n==1"
g 5   paper_simulations "44 simulated trials under the alternative and 44 under the null" "n==1"
g 6   paper_simulations "Identifying Patient Subgroups" "n==0"
g 7   treatment_effect_definitions "@ref(" "n==0"
g 8   treatment_effect_definitions "The Post-Selection Target: the Standard Analysis of the Selected Region" "n>=1"
g 9   biomarker_effects "generate_aft_dgm_flex" "n>=2"
g 10a causal_effects_brief_review "What forestsearch Reports, Read Against This Review" "n>=1"
g 10b causal_effects_brief_review "not a claim of causal optimality" "n==2"
g 11  survival-post-selection "which this call reproduces" "n==1"
for n in $ALL; do
  g 12a $n "NaN" "n==0"
  g 12b $n "?@" "n==0"
done
[ $fail = 0 ] && echo "TEXT GATES: ALL PASS" | tee -a "$KV" || echo "TEXT GATES: FAIL" | tee -a "$KV"
exit $fail
