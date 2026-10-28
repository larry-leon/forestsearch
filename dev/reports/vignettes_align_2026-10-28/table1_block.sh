#!/bin/bash
# §1.3 gate: Table 1 block (from the "Table 1" section heading to the next section) of two rendered HTMLs
# usage: table1_block.sh <pandoc> <old.html> <new.html> <outdir>
set -e
for t in old new; do
  f=$2; [ $t = new ] && f=$3
  "$1" -f html -t plain --wrap=none "$f" -o "$4/plain_$t.txt"
  awk '/^Table 1$/{on=1} /^Declaration at a protected level$/{on=0} on' "$4/plain_$t.txt" > "$4/table1_$t.txt"
  echo "$t: $(wc -l < "$4/table1_$t.txt") lines, md5 $(md5sum < "$4/table1_$t.txt" | cut -d' ' -f1)"
done
diff "$4/table1_old.txt" "$4/table1_new.txt" && echo "Table 1 block identical: TRUE"
