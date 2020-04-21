## Loading all necessary libraries to run this on a stat-machine

library(arm)
library(data.table)
library(magrittr)
library(purrr)
library(tidyr)
library(ggplot2)
library(pals)

## parallelization
library(furrr) # future + purrr
options(mc.cores = 8)
plan(multiprocess)
