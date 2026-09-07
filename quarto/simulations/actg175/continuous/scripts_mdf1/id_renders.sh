#!/bin/bash
cd "$(dirname "$0")" >/dev/null; cd ~/Documents/GitHub/forestsearch/quarto/simulations/actg175/continuous || exit 2
SP=/private/tmp/claude-501/-Users-larryleon-Downloads/476ca7fd-7bca-4851-bd29-048be32a2d48/scratchpad
for cell in "40 500" "40 700" "120 500" "null 500"; do
  set -- $cell
  for ci in ij field; do
    t0=$(date +%s)
    FS_MD_MD=$1 FS_MD_N=$2 FS_MD_CI=$ci FS_MD_NSIMS=5 FS_MD_CAMPAIGN=id$ci \
      quarto render sim_fs_maxeffCons_mr_field_md_template.qmd --output id_md$1_n$2_$ci.html > $SP/id_md$1_n$2_$ci.log 2>&1
    echo "md=$1 n=$2 ci=$ci -> exit $? in $(( $(date +%s) - t0 )) s"
  done
done
ls mr_md_harm | grep -E '_id(ij|field)_'
