#!/bin/bash
cd "$(dirname "$0")"
for s in 1 2 3 5 4; do for N in 5000 20000 100000; do
  echo "$(date +%T) start s$s N$N"; Rscript sim_study.R $s $N 100 4
done; done
echo ALLDONE
