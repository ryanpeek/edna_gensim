# generate_ms_call
library(tidyverse)
library(purrr)
# options(scipen = 999)
library(glue)

# SIM 001 -----------------------------------------------------------------

# this simulation has no error or haplotype requirements in the perl script.
# requires perl script ms_sample.pl

# now generate the actual call using purrr
# need the following: nInd=$1, nLoci=$2, ll=$3, cov=$4

# SIMULATION 001 (no error rate)
# sim001 <- list(
#   nInd = c(1:50),
#   nLoci = c(100, 1000),
#   ll = c(100, 1000, 10000),
#   coverage=c(10, 100)
# )
# 
# # make all combos
# args <- sim001 %>% purrr::cross_df() %>% 
#   mutate(call_list=paste(nInd, nLoci, ll, coverage, sep = " "))
# 
# # make the sh call with glue: "sh simulate_haplotypes.sh 1 100 100 10"
# args <- args %>% 
#   mutate(local_runs = glue('sh ms_simulate_haplotypes.sh {nInd} {nLoci} {ll} {coverage}'))
# 
# # make bash header for shell file
# cat('#!/bin/bash -l\n\n', file = "code/ms_sim001_sbatch.sh")
# readr::write_lines(args$local_runs, path = "code/ms_sim001_sbatch.sh", append = T)



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


