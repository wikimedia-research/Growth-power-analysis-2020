## Loading all necessary libraries to run this on a stat-machine

# for %>%
library(dplyr)

## Set BLAS threads to 1 because we parallelize using furrr
library(RhpcBLASctl)
blas_set_num_threads(1)

## parallelization
library(furrr) # future + purrr
options(mc.cores = 10)
plan(multiprocess)
