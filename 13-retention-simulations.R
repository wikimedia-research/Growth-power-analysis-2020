## Run all simulations to calculate power for newcomer retention

## Assumes that 06-stat-load-libraries.R as well as the functions
## from 07-activation-simulations.R are all loaded into memory.

## 1: Load in the parameters
load('simulation_parameters/first4_retention_params.txt', verbose = TRUE)
## 2: Run simulations
first4_retention_results = run_power_simulations(act_counts_first4, ret_coeff_first4)

## 3: save output as a TSV
write.table(first4_retention_results, file = 'datasets/first4-retention-results.tsv',
            sep = '\t', quote = FALSE, row.names = FALSE)

