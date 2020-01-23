# generate_ms_call
library(tidyverse)
library(purrr)
# options(scipen = 999)
library(glue)

# Simple ms call ----------------------------------------------------------

nInd <- seq(1:100)
nLoci <- 100 # mean number of loci analyzed (100 or 1000)
ploidy <- 2 # diploid individuals
gc<- nInd * ploidy # gene copies (for ms)

mu <- 1e-8 # mutation rate
Ne <- 25000 # coalescent effect pop size
ll <- 100 # loci length (100, 1000, 10000)
theta <- 4*Ne*ll*mu # scaled mutation rate

coverage <- 10 # mean coverage per locus analyzed (10, 100)

# generate single combo list
#call <- paste(nInd, nLoci, ll, coverage, sep = ",")
#call # double check

# Generate function call --------------------------------------------------

# now generate the actual call using purrr
# need the following:
# nInd=$1, nLoci=$2, ll=$3, cov=$4

# create a list of combos
data <- list(
  nInd = c(1:50),
  nLoci = c(100, 1000),
  ll = c(100, 1000, 10000),
  coverage=c(10, 100)
)

args <- data %>% purrr::cross_df() %>% 
  mutate(call_list=paste(nInd, nLoci, ll, coverage, sep = " "))


# A SINGLE LIST OF JUST PARAMS --------------------------------------------

# write out as single list
readr::write_lines(args$call_list, path = "data/ms_params_multi")

# append to existing file?
# readr::write_lines(call, path = "data/ms_params_multi.txt", append = T)



# Make List Call for local runs -------------------------------------------

# want to generate sh calls
# "sh simulate_haplotypes.sh 1 100 100 10"

args <- args %>% 
  mutate(local_runs = glue('sh ms_simulate_haplotypes.sh {nInd} {nLoci} {ll} {coverage}'))

readr::write_lines(args$local_runs, path = "data/sh_runs_multi")



# NOTES -------------------------------------------------------------------

### maybe just need this: (can run locally with ~/Downloads/msdir/ms 4 2 -t 1)
### srun -t 20 simulate_haplotypes.sh 10 100 100 10

