# DLM-River

## Joshua Christensen and Matthew Heiner

Benchmarking and performance comparisons among samplers on the non-Gaussian DLM 
applied to time series of nitrate concentrations in a river in France.

French water chemistry time series were provided by Naïades, 
extraction date 24 Nov. 2018, url: http://www.naiades.eaufrance.fr/france-entiere#/
The time series used in this example are from station id 3052134. 

The R scripts in this folder are numbered by their order in the workflow sequence. 
The general workflow is to schedule runs for the desired targets and samplers by round 
(round "tune" is for initial testing of settings within samplers; round "test" is for full
randomized comparison at selected settings).

Begin with script `1_schedule_river.R` to produce an 
initial schedule for selected samplers and round. 
To create a schedule for the simulated time series, change the variable `datatype` to "sim".
The resulting simulation schedule is saved in the `data` folder. 

The script `2_run_dhr_server.R` executes the sampler for each scheduled job. 
It can also be run in a single interactive instance and explored with `2_validate_results_sim.R` 
(uncomment and run the two sections labeled "for manual interactive run" in `2_run_dhr_server.R`). 
The shell script `runallTrials.sh` can be used to run all jobs in the schedule with GNU Parallel. 
Summaries for each run are saved in the `output` folder. 

Once all jobs in a round are completed, scripts `3_combine_output.R` and 
`4_summary.R` summarize sampler performance.
