#!/bin/bash
cd ~/Documents/GitHub/forestsearch/quarto/simulations/actg175/continuous || exit 2
SP=/private/tmp/claude-501/-Users-larryleon-Downloads/476ca7fd-7bca-4851-bd29-048be32a2d48/scratchpad
run() { # md n workers nsims tag
  t0=$(date +%s)
  FS_MD_MD=$1 FS_MD_N=$2 FS_MD_WORKERS=$3 FS_MD_NSIMS=$4 FS_MD_CI=field FS_MD_CAMPAIGN=$5 \
    quarto render sim_fs_maxeffCons_mr_field_md_template.qmd --output cal_$5.html > $SP/cal_$5.log 2>&1
  echo "cal md=$1 n=$2 workers=$3 nsims=$4 tag=$5 -> exit $? in $(( $(date +%s) - t0 )) s total render"
  grep -E '^Batch sim_id|^This run' $SP/cal_$5.log | head -2
}
run 40 500 13 26 cal13n500
run 40 500 10 20 cal10n500
run 40 700 13 26 cal13n700
