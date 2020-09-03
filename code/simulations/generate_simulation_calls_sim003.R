# generate_ms_call
library(tidyverse)
library(purrr)
# options(scipen = 999)
library(glue)


# SIM 003 & 004  -----------------------------------------------------------------

# this simulation adds distributions for each individual, to simulate different
# eDNA shedding rates

resultsDir <- "sim003"

# SIMULATION 003 (10% error, haplo=< 2 threshold, distrib)
sim003 <- list(
  nInd = c(seq(5,35,5)),
  #nLoci = c(100),
  #ll = c(10000),
  #coverage=c(100),
  distrib=c("norm", "unif", "gamma", "beta", "betaU", "gammaX", "equal")
  #distrib=c("betaU", "gammaX", "equal")
  #reps=c(100) # number of simulations/repetitions
)

# make all combos
args <- sim003 %>% purrr::cross_df() %>% 
  mutate(call_list=paste(nInd, distrib, sep = " "))

# make the sh call with glue: "sh code/ms_simulate_haplotypes_local.sh 9 norm sim004"


args <- args %>% 
  mutate(farm_runs = glue('sh ms_simulate_haplotypes_v3.sh {nInd} {distrib} {resultsDir}'),
         local_runs = glue('sh code/ms_simulate_haplotypes_v3.sh {nInd} {distrib} {resultsDir}'))

# preview
head(args$local_runs)
#head(args$farm_runs)

# WRITE FOR LOCAL RUNS
# make bash header for shell file
cat('#!/bin/bash -l\n\n', file = "code/sim003_sbatch_local.sh")
readr::write_lines(args$local_runs, path = "code/sim003_sbatch_local.sh", append = T)


# WRITE FOR FARM RUNS
# make bash header for shell file
cat('#!/bin/bash -l\n\n', file = "code/sim003_sbatch_farm.sh")
readr::write_lines(args$farm_runs, path = "code/sim003_sbatch_farm.sh", append = T)



# A SINGLE LIST OF JUST PARAMS --------------------------------------------

# write out as single list
#readr::write_lines(args$call_list, path = "data/ms_params_50inds_2020_01_22.txt")

# append to existing file?
# readr::write_lines(call, path = "data/ms_params_multi", append = T)

# RUN SIMULATIONS ---------------------------------------------------------

## Put scripts on server:
# sftp put ms_sim002_sbatch.sh
## SBATCH RUN OF ALL FILES: 
# sbatch -t 2800 -p high ms_sim002_sbatch.sh

## RUN LOCAL/SINGLE RUN (nInd, distrib, outdir)
# sh code/ms_simulate_haplotypes_v3.sh 5 norm sim003


