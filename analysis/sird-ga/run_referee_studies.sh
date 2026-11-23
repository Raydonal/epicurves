#!/bin/bash
# Referee studies for the SIRD-GA article, run sequentially (each stage resumable: finished
# chunk/phase files are skipped). Progress and DONE markers in results/referee_run.log; each stage
# also writes its own log. Uses at most 3 cores (machine shared).
#   stage 2  paraiba_blocks.R            -> results/paraiba_blocks.log            BLOCKSDONE
#   stage 2b paraiba_blocks_warm.R       -> results/paraiba_blocks_warm.log       BLOCKSWARMDONE
#   stage 3  paraiba_floor.R             -> results/paraiba_floor.log             FLOORDONE
#   stage P  paraiba_profile_p_refined.R -> results/paraiba_profile_p_refined.log PROFILEDONE
#   stage 1  sim_paraiba_regime.R A,B,C  -> results/sim_paraiba_regime_<sc>.log   REGIMEDONE <sc>
#            summarise_paraiba_regime.R  -> results/paraiba_regime_summary.log
#   stage 4  sim_generic_starts.R s1,s4 x N 20000,100000 -> results/sim_generic_starts.log  GENERICDONE
#            summarise_generic_starts.R  -> results/sim_generic_summary.log
# Launch: setsid nohup ./run_referee_studies.sh >> results/referee_run.log 2>&1 < /dev/null &
cd "$(dirname "$0")"
CORES=3
echo "$(date '+%F %T') start referee studies"
Rscript paraiba_blocks.R $CORES >> results/paraiba_blocks.log 2>&1
grep -q BLOCKSDONE results/paraiba_blocks.log && echo "$(date '+%F %T') BLOCKSDONE"
Rscript paraiba_floor.R $CORES >> results/paraiba_floor.log 2>&1
grep -q FLOORDONE results/paraiba_floor.log && echo "$(date '+%F %T') FLOORDONE"
Rscript paraiba_blocks_warm.R $CORES >> results/paraiba_blocks_warm.log 2>&1
grep -q BLOCKSWARMDONE results/paraiba_blocks_warm.log && echo "$(date '+%F %T') BLOCKSWARMDONE"
Rscript paraiba_profile_p_refined.R $CORES >> results/paraiba_profile_p_refined.log 2>&1
grep -q PROFILEDONE results/paraiba_profile_p_refined.log && echo "$(date '+%F %T') PROFILEDONE"
for sc in A B C; do
  Rscript sim_paraiba_regime.R $sc "" $CORES >> results/sim_paraiba_regime_$sc.log 2>&1
  grep -q "REGIMEDONE $sc" results/sim_paraiba_regime_$sc.log && echo "$(date '+%F %T') REGIMEDONE $sc"
done
Rscript summarise_paraiba_regime.R > results/paraiba_regime_summary.log 2>&1
echo "$(date '+%F %T') REGIMESUMMARYDONE"
for s in 1 4; do for N in 20000 100000; do
  Rscript sim_generic_starts.R $s $N 100 $CORES >> results/sim_generic_starts.log 2>&1
  grep -q "GENERICDONE $s $N" results/sim_generic_starts.log && echo "$(date '+%F %T') GENERICDONE $s $N"
done; done
Rscript summarise_generic_starts.R > results/sim_generic_summary.log 2>&1
echo "$(date '+%F %T') ALLREFEREEDONE"
