# make haplotype plots


# Libraries ---------------------------------------------------------------

library(tidyverse)
#library(tidylog)


# Get Data ----------------------------------------------------------------

dat <- read_tsv("data/run.out", col_names = FALSE)


# Make Plots --------------------------------------------------------------

# plot # of sims by mean # of haplotypes per locus
dat_means <- dat %>% 
  summarize_all(mean, na.rm=TRUE) %>% 
  pivot_longer(cols=starts_with("X"), names_to = "sims")

hist(dat_means$value)


# from one simulation
# number of loci vs. # of haplotypes
hist(dat$X77)
