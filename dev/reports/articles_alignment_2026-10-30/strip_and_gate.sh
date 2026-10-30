#!/usr/bin/env bash
# strip_and_gate.sh -- task sec. 4 text gates, output only. Run from the repo root.
# Method (as dev/reports/guohe_article_2026-10-30/s2_text_gates.sh): every <pre ...>...</pre>
# element is deleted from the rendered HTML (R gsub, perl, dot-all "(?s)<pre[ >].*?</pre>"),
# then pandoc --wrap=none -f html -t plain -> results/<name>_output.txt.
# <name>.txt is the whole page, for the record. Prints one line per gate; exit 1 on any FAIL.
set -u
OUT=/tmp/fs_align_render_2026-10-30
R=dev/reports/articles_alignment_2026-10-30/results
fail=0
strip() {
  local n=$1
  pandoc --wrap=none -t plain "$OUT/$n.html" > "$R/$n.txt"
  Rscript -e 'a <- commandArgs(TRUE); h <- paste(readLines(a[1], warn = FALSE), collapse = "\n");
    n <- lengths(regmatches(h, gregexpr("<pre[ >]", h)));
    h2 <- gsub("(?s)<pre[ >].*?</pre>", "", h, perl = TRUE);
    stopifnot(!grepl("<pre[ >]", h2)); writeLines(h2, a[2]); cat(a[3], "pre_elements_removed =", n, "\n")' \
    "$OUT/$n.html" "$R/${n}_nopre.html" "$n" || { echo "STRIP FAILED $n"; exit 2; }
  pandoc --wrap=none -f html -t plain "$R/${n}_nopre.html" > "$R/${n}_output.txt"
  rm "$R/${n}_nopre.html"
}
# g id article mode(F fixed|E regex) pattern condition
g() {
  local c
  if [ "$3" = F ]; then c=$(grep -oF -- "$4" "$R/$2_output.txt" | wc -l); else c=$(grep -oE -- "$4" "$R/$2_output.txt" | wc -l); fi
  local ok; ok=$(awk -v n="$c" "BEGIN{print ($5) ? \"PASS\" : \"FAIL\"}")
  [ "$ok" = FAIL ] && fail=1
  echo "$1 [$2] '$4' = $c ; expect $5 ; $ok"
}
for n in methodology paper_simulations treatment_effect_definitions biomarker_effects causal_effects_brief_review; do strip $n; done
g g1  methodology F "2026 note" "n==13"
g g2a methodology F "consistency_method" "n>=1"
g g2b methodology F "fs_declaration_calibration" "n>=1"
g g2c methodology F "forestsearch_bootstrap_dofuture" "n>=1"
g g2d methodology F "hrMaxSG" "n>=1"
g g3a methodology F "1.276" "n>=1"
g g3b methodology F "0.634" "n>=1"
g g3c methodology F "0.645" "n>=1"
g g3d methodology F "0.804" "n>=1"
g g4  methodology F "These are the 2024 guidelines" "n==1"
g g5a paper_simulations F "Identifying Patient Subgroups" "n==0"
g g5b paper_simulations F "Exploratory subgroup identification in the heterogeneous Cox model" "n==1"
g g6a paper_simulations F "arXiv:2609.38361" "n>=2"
g g6b paper_simulations F "arXiv:2610.28379" "n>=2"
g g7a paper_simulations F "44 simulated trials under the alternative and 44 under the null" "n==1"
g g7b paper_simulations F "n_sims = 1000" "n==0"
g g8a paper_simulations F "random-split consistency" "n==0"
g g8b paper_simulations F "consistency_method" "n>=1"
g g9a treatment_effect_definitions F "The Post-Selection Target: the Standard Analysis of the Selected Region" "n==2"
g g9b treatment_effect_definitions F "identification factor" "n==1"
g g10a treatment_effect_definitions F "default estimand for ForestSearch simulation" "n==0"
g g10b treatment_effect_definitions F "arXiv:2609.38361" "n>=1"
g g10c treatment_effect_definitions F "arXiv:2610.28379" "n>=1"
g g11a biomarker_effects F "create_gbsg_dgm" "n>=1"
g g11b biomarker_effects F "generate_aft_dgm_flex" "n>=2"
g g11c biomarker_effects F "arXiv:2609.38361" "n>=1"
g g11d biomarker_effects F "arXiv:2610.28379" "n>=1"
g g12a causal_effects_brief_review F "What forestsearch Reports, Read Against This Review" "n==2"
g g12b causal_effects_brief_review F "not a claim of causal optimality" "n==2"
g g13a causal_effects_brief_review F "GRF/LASSO" "n==0"
g g13b causal_effects_brief_review F "CHR alongside AHR" "n==0"
g g13c causal_effects_brief_review F "biased comparator" "n==0"
g g13d causal_effects_brief_review F "GRF baseline covariate conditioning" "n==0"
g g14a causal_effects_brief_review F "arXiv:2609.38361" "n>=1"
g g14b causal_effects_brief_review F "arXiv:2610.28379" "n>=1"
for n in methodology paper_simulations treatment_effect_definitions biomarker_effects causal_effects_brief_review; do
  g g15 $n F "NaN" "n==0"
done
# g16: numbered-list hazard on the source vs master c2967362
for n in methodology paper_simulations treatment_effect_definitions biomarker_effects causal_effects_brief_review; do
  git show c2967362:vignettes/articles/$n.Rmd > "$R/.master_$n.Rmd"
  new=0
  while IFS= read -r line; do
    if ! grep -qxF -- "$line" "$R/.master_$n.Rmd"; then
      echo "  g16 [$n] not on master: $line"; new=$((new + 1))
      if [ "$n" = paper_simulations ] && echo "$line" | grep -qE '^(1|4|5)\. '; then new=$((new - 1)); echo "    (allowed: B6/B7 reference item)"; fi
    fi
  done < <(grep -E '^\s*[0-9]+[.)]( |$)' "vignettes/articles/$n.Rmd")
  rm "$R/.master_$n.Rmd"
  ok=PASS; [ $new -ne 0 ] && { ok=FAIL; fail=1; }
  echo "g16 [$n] numbered lines not on master (excluding allowed) = $new ; expect 0 ; $ok"
done
exit $fail
