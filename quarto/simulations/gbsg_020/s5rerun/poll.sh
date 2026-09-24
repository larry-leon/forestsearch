#!/bin/bash
# quiet poll: wait up to ~9.5 min or until the run ends; print one status line
cd /Users/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020
L=s5rerun/logs/driver.log; N0=$(grep -c "CELL DONE\|NOT COMPLETE" $L)
for i in $(seq 1 57); do grep -q "RUN END" $L && break; [ $(grep -c "CELL DONE\|NOT COMPLETE" $L) -gt $N0 ] && break; sleep 10; done
echo "$(date +%T) | done=$(grep -c 'CELL DONE' $L)/18 | last: $(grep -E 'CELL|RUN|MEMORY|OVERRUN|FAIL' $L | tail -1 | cut -c1-110) | mem: $(tail -1 s5rerun/logs/s5rerun_mem.log | cut -d' ' -f2-)"
