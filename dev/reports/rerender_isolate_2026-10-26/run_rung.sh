#!/usr/bin/env bash
# usage: run_rung.sh <D: C|B|G> <rung> <cap, e.g. 20m> [worktree SHA or "-"]
# Renders a copy of document D in its own run directory (§4), writes build.txt, removes the copied .qmd,
# then runs the legacy-tree MD5 gate. Exit: 0 ok; 3 render failed/capped (ladder continues); 9 GATE FAILED (stop).
REPO=/home/larryleon/Documents/GitHub/forestsearch
REC=$REPO/dev/reports/rerender_isolate_2026-10-26
D=$1; rung=$2; cap=$3; sha=${4:--}
case $D in
  C) app=actg175; stem=analysis_actg175_continuous_intervals ;;
  B) app=actg175; stem=analysis_actg175_binary_multimethod_fixed_family ;;
  G) app=gbsg;    stem=analysis_gbsg_survival_effMaxSG ;;
esac
run=$REPO/quarto/applications/$app/rerun_2026-10-26/${D}_${rung}
mkdir -p "$run"
cp "$REPO/quarto/applications/$app/$stem.qmd" "$run/"
if [ "$D" = G ]; then mkdir -p "$run/_payloads"; export GH_DIR="$run/_payloads"; fi
export RSTUDIO_PANDOC=/usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64
{
  echo "document: $D ($stem.qmd at main-tree HEAD $(git -C $REPO rev-parse HEAD))"
  echo "rung: $rung"
  echo "worktree SHA: $sha"
  echo "date: $(date -Is)"
  [ "$D" = G ] && echo "env: GH_DIR=$GH_DIR (LOO cache redirected into the run directory; see inventory.txt)"
  Rscript -e 'pd <- packageDescription("forestsearch"); for (k in c("Version","Packaged","Built","Repository")) cat(k, ": ", if (is.null(pd[[k]])) NA else pd[[k]], "\n", sep = ""); cat("exists(.fs_logit_two_group):", exists(".fs_logit_two_group", asNamespace("forestsearch")), "\n")'
} > "$run/build.txt" 2>&1
cp "$run/build.txt" "$run/render.log"
echo "----- quarto render -----" >> "$run/render.log"
cd "$run"
t0=$(date +%s)
timeout "$cap" quarto render "$stem.qmd" >> render.log 2>&1
rc=$?
t1=$(date +%s)
echo "render: exit=$rc wall_s=$((t1-t0)) cap=$cap" >> build.txt
echo "$D $rung exit=$rc wall_s=$((t1-t0)) cap=$cap sha=$sha $(date -Is)" >> "$REC/timing.txt"
rm -f "$run/$stem.qmd"
rm -rf "$run/${stem}_files"
cd "$REPO"
# GATE: legacy tree unchanged
find quarto/applications/gbsg quarto/applications/actg175 -type f -not -path '*/rerun_2026-10-26/*' -exec md5sum {} + | sort -k2 \
  | diff - "$REC/legacy_tree_md5_before.txt" > "$REC/gate_last.diff"
g=$?
gh=$(md5sum quarto/GuoHe/cv_out_analysis_gbsg_survival_effMaxSG_hrMaxSG_neighborhood.rds | cut -d' ' -f1)
echo "$D $rung legacy-tree-diff-rc=$g guohe_cache_md5=$gh" >> "$REC/gate_log.txt"
if [ $g -ne 0 ] || [ "$gh" != 4844ae280990185655bed5a29b2d99bb ]; then echo "GATE FAILED after $D $rung"; cat "$REC/gate_last.diff"; exit 9; fi
[ $rc -ne 0 ] && { echo "render $D $rung exit=$rc"; tail -30 "$run/render.log"; exit 3; }
echo "OK $D $rung wall_s=$((t1-t0))"
