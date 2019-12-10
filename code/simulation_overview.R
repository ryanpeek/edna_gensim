# packages for sim:
# https://github.com/statgenlmu/coala
# https://academic.oup.com/bioinformatics/article/32/12/1903/1744397/
# whittaker plot Dunn et al 2007 Myrmar

library(purrr)
library(tidyverse)
library(coala)
library(phyclust) # to use ms sim
activate_ms(priority = 500)

# list simulators
list_simulators()

# check if model can be run
model <- coal_model(sample_size = 10, loci_number = 50, 
                    loci_length = 300, ploidy = 2) +
  feat_mutation(rate = 5, model = "IFS") + # rate= 4*N*mu (mutation rate per locus)
  feat_mutation(par_prior("theta", sample.int(100, 1))) +
  sumstat_nucleotide_div()
check_model(model)
model # see model
stats <- simulate(model, seed = 15, nsim=1000)

# plot
mean_pi <- map_dbl(stats, ~mean(.x$pi))
theta <- map_dbl(stats, "pars", .default="theta")
plot(theta, mean_pi, pch = 19, col = "orange")
abline(lm(mean_pi ~ theta), col = "red2", lty = 3)


# 20 samples from one pop, 
model <- coal_model(sample_size = 10, loci_number = 100, 
                    loci_length = 300, ploidy = 2) +
  feat_mutation(.1, model="IFS") +
  sumstat_tajimas_d()
check_model(model)
model
stats <- simulate(model, seed = 15, nsim=1000)


mean_tajD <- map_dbl(stats, ~mean(.x$tajimas_d, na.rm=T))
plot(density(mean_tajD, na.rm = TRUE), 
     main = "Neutral Distribution of Tajima's D")

# SFS
model <- coal_model(sample_size = 10, loci_number = 100, 
                    loci_length = 300, ploidy = 2) +
  feat_mutation(1, model="IFS") +
  sumstat_sfs()
check_model(model)
model

# option 1: 1 run
stats1 <- simulate(model, seed = 15)
barplot(stats1$sfs / sum(stats1$sfs), 
        names.arg = seq_along(stats1$sfs), 
        col = 3)
# ggplot option
sfs_out <- map(stats1[1], ~.x / sum(.x)) %>% purrr::as_vector() %>% 
  tibble(sfs=., numb=seq_along(stats1[[1]]))
ggplot() + geom_col(data=sfs_out, aes(x=numb, y=sfs), fill="cadetblue4")


# option 2: 1000 runs:
stats2 <- simulate(model, seed = 15, nsim=1000)

# extract col wise sfs from each sim
mean_sfs_col <- map_dfc(stats2, ~data_frame("sfs"=.x$sfs)) %>% 
  mutate(sfs_mean = rowMeans(.)) %>% 
  select(sfs_mean)
  
# ggplot option
(sfs_out <- mean_sfs_col %>% 
  mutate(sfs_no = as_factor(rownames(mean_sfs_col)),
         sfs_calc = (sfs_mean / sum(sfs_mean))) %>% 
  ggplot() + ggdark::dark_theme_classic(base_family = "Roboto Condensed") +
  geom_col(aes(x=sfs_no, y=sfs_calc), fill="cadetblue4") + 
    labs(title = "SFS", y="Site Freq.", x="Chromosomes (derived allele freq)"))


# Another Pop Gen in R Package --------------------------------------------

#https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.5412