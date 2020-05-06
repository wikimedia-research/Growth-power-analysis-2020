## Run all simulations to calculate power for newcomer retention

## Assumes that 06-stat-load-libraries.R as well as the functions
## from 07-activation-simulations.R are all loaded into memory.

run_retention_simulations = function(act_counts, coeffs) {
  ## Run simulations of user retention with the given activation counts
  ## and coefficients. Returns a flattened data.table with results for
  ## 10%, 5%, and 2% effect size.

  ## We calculate retention power on a monthly basis from 1 to 6 months,
  ## in half-month increments
  multipliers <- seq(1, 6, by = 0.5)
  # So the output is also named, with the estimated power for each multiplier:
  multipliers <- setNames(multipliers, multipliers)
  # Then use the ~ syntax to define an anonymous function applied to each multiplier:
  retention_powers = list()
  retention_powers[['10%']] <- purrr::map_dbl(
    multipliers,
    ~ estimate_power(floor(act_counts * .x), coeffs)
  )

  ## What if the effect is 5%?
  coeffs[['beta']] = 0.2
  coeffs[['sigma_beta']] = 0.1

  retention_powers[['5%']] <- purrr::map_dbl(
    multipliers,
    ~ estimate_power(floor(act_counts * .x), coeffs)
  )

  ## What if the effect is 2%?
  coeffs[['beta']] = 0.08
  coeffs[['sigma_beta']] = 0.04

  retention_powers[['2%']] <- purrr::map_dbl(
    multipliers,
    ~ estimate_power(floor(act_counts * .x), coeffs)
  )
  
  ## Flatten the list into a data.frame
  x = data.table::rbindlist(lapply(retention_powers, as.data.frame.list))
  x[, `Effect size` := ordered(c("10%", "5%", "2%"), c("10%", "5%", "2%"))]
  return(x)
}


## 1: Load in the parameters
load('simulation_parameters/first4_retention_params.txt', verbose = TRUE)
## 2: Run simulations
## first4_retention_results = run_retention_simulations(act_counts_first4, ret_coeff_first4)

## 3: save output as a TSV
##write.table(first4_retention_results, file = 'datasets/first4-retention-results.tsv',
##            sep = '\t', quote = FALSE, row.names = FALSE)

## Current set of wikis
load('simulation_parameters/current_retention_params.txt', verbose = TRUE)

## current_retention_results = run_retention_simulations(act_counts_current, ret_coeff_current)

write.table(current_retention_results, file = 'datasets/current-retention-results.tsv',
            sep = '\t', quote = FALSE, row.names = FALSE)

## Add French
load('simulation_parameters/addfr_retention_params.txt', verbose = TRUE)

## addfr_retention_results = run_retention_simulations(act_counts_addfr, ret_coeff_addfr)

write.table(addfr_retention_results, file = 'datasets/addfr-retention-results.tsv',
            sep = '\t', quote = FALSE, row.names = FALSE)

## Add all
load('simulation_parameters/addall_retention_params.txt', verbose = TRUE)

## addall_retention_results = run_retention_simulations(act_counts_addall, ret_coeff_addall)

write.table(addall_retention_results, file = 'datasets/addall-retention-results.tsv',
            sep = '\t', quote = FALSE, row.names = FALSE)
