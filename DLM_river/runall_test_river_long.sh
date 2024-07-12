export OMP_NUM_THREADS=1

dte=240708

num_of_trials=$(Rscript scripts/num_of_trials.R river_long test)
echo $num_of_trials

for data in river_long
do
  for round in test
  do
    for ii in $( seq 1 $num_of_trials )
    do
      echo "Rscript --no-restore 2_run_dhr_server.R $data $round $ii $dte &> logs/$data\_$round\_job$ii\_dte$dte.txt"
    done
  done
done | parallel --jobs 10
wait

echo "finished"
