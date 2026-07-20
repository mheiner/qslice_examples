# Normalized Weights

## Matthew Heiner

Benchmarking and performance comparisons among samplers in the ``normalized weights" or reweighted Multinomial model example.

The R scripts in this folder are numbered by their order in the workflow sequence. 
The general workflow is to schedule runs for the desired targets and samplers.

Begin with script `1_schedule_all.R` to produce an 
initial schedule for selected samplers. 
The schedule is saved in the `input` folder. 

The script `2_run_trials_all.R` executes the sampler for each scheduled job separately. 
It can also be run in a single instance to explore sampler output. 
The shell script `runallTrials.sh` can be used to run all jobs in the schedule with GNU Parallel. 
Summaries for each run are saved in the `output` folder. 

Once all jobs in a round are completed, scripts `3_combine_output.R` and 
`4_summarize_combined.R` summarize sampler performance.
