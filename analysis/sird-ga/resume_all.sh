#!/bin/bash
# Runs the profile of alpha and the SIRD-GA simulation grid, skipping finished combinations.
cd "$(dirname "$0")"
echo "$(date +%T) start profile_alpha_sim"
Rscript profile_alpha_sim.R > results/profile_alpha_sim.log 2>&1
echo "$(date +%T) done profile_alpha_sim (exit $?)"
for s in 1 2 3 5 4; do for N in 5000 20000 100000; do
  if [ -f results/sim_s${s}_N${N}.rds ]; then echo "skip s$s N$N"; continue; fi
  echo "$(date +%T) start s$s N$N"; Rscript sim_study.R $s $N 100 4
done; done
echo ALLDONE
