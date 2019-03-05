#! /bin/bash

# Wrapper script to loop our comparison bash script

source $1

REP_TEMPFILE=/tmp/rep_temp.tmp
echo 0 > $REP_TEMPFILE

for ((n=0;n<$2;n++)); do

  REP_COUNTER=$[$(cat $REP_TEMPFILE) + 1]

  ./comparison_commands.sh $1 > comparison_run_log.txt 2>&1

  wait

  # move all the valuable files into a folder that shows which run its a part of

  mv ./combined_outputs ./combined_output-$REP_COUNTER
  mv $outdir ./$outdir-$REP_COUNTER
  mv ./phycorder_results ./phycorder_results-$REP_COUNTER
  mv ./gon_phy_runs_dir ./gon_phy_runs_dir-$REP_COUNTER
  mv ./gon_phy_results ./gon_phy_results-$REP_COUNTER

  echo $REP_COUNTER > $REP_TEMPFILE

done

unlink $REP_TEMPFILE
