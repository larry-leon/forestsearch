#!/usr/bin/env bash
# render_all.sh -- task sec. 3: render the five edited articles out of tree, in the stated order.
# Run from the repo root. Stops at the first non-zero exit (a timeout exits 124).
D=dev/reports/articles_alignment_2026-10-30
OUT=/tmp/fs_align_render_2026-10-30
for spec in methodology:10m treatment_effect_definitions:10m biomarker_effects:10m \
            causal_effects_brief_review:10m paper_simulations:45m; do
  name=${spec%%:*}; cap=${spec##*:}
  s=$(date +%s); echo "start $(date -Is)" > $D/logs/time_$name.txt
  timeout $cap Rscript -e "rmarkdown::render('vignettes/articles/$name.Rmd', output_dir = '$OUT', intermediates_dir = '$OUT/intermediates_$name', quiet = FALSE)" > $D/logs/render_$name.log 2>&1
  rc=$?; e=$(date +%s)
  echo "end $(date -Is)" >> $D/logs/time_$name.txt
  echo "$name exit=$rc cap=$cap wall_s=$((e - s))" | tee -a $D/logs/render_summary.txt
  [ $rc -ne 0 ] && { echo "STOP at $name"; exit $rc; }
done
echo "ALL RENDERS OK"
