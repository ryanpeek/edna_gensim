# generate_ms_call
library(tidyverse)
library(purrr)
# options(scipen = 999)
library(glue)


# SIM 005  -----------------------------------------------------------------

# this simulation uses just gamma and varies coverage, ll, and loci

resultsDir <- "sim005"

fs::dir_create(glue::glue("results/{resultsDir}"))

# SIMULATION 005 (10% error, haplo=< 2 threshold, distrib)
sim <- list(
  nInd = c(seq(5,40,5)),
  nLoci = c(100, 250, 500, 1000),
  ll = c(1000, 2500, 5000, 10000),
  coverage=c(10, 25, 50, 100),
  distrib=c("gamma")
  #reps=c(100) # number of simulations/repetitions
)

# make all combos
args <- sim %>% purrr::cross_df() %>% 
  mutate(call_list=paste(nInd, distrib, sep = " "))

# set up shell call
args <- args %>% 
  mutate(farm_runs = glue('sh ms_simulate_haplotypes_v3.sh {nInd} {nLoci} {ll} {coverage} {distrib} {resultsDir}'),
         local_runs = glue('sh code/ms_simulate_haplotypes_v3.sh {nInd} {nLoci} {ll} {coverage} {distrib} {resultsDir}'))

# preview
head(args$local_runs)
#head(args$farm_runs)

# WRITE FOR LOCAL RUNS
# make bash header for shell file
cat('#!/bin/bash -l\n\n', file = glue::glue("code/{resultsDir}_sbatch_local.sh"))
readr::write_lines(args$local_runs, path = glue::glue("code/{resultsDir}_sbatch_local.sh"), append = T)


# WRITE FOR FARM RUNS
# make bash header for shell file
cat('#!/bin/bash -l\n\n', file = glue::glue("code/{resultsDir}_sbatch_farm.sh"))
readr::write_lines(args$farm_runs, path = glue::glue("code/{resultsDir}_sbatch_farm.sh"), append = T)


# RUN SIMULATIONS ---------------------------------------------------------

## Put scripts on server:
# sftp put ms_sim002_sbatch.sh
## SBATCH RUN OF ALL FILES: 
# sbatch -t 2800 -p high ms_sim002_sbatch.sh

## RUN LOCAL/SINGLE RUN (nInd, distrib, outdir)
# sh code/ms_simulate_haplotypes_v3.sh 5 norm sim003


