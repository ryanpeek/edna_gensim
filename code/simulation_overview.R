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
model <- coal_model(10, 50, loci_length = 300, ploidy = 2) +
  feat_mutation(5, model = "IFS") +
  feat_mutation(par_prior("theta", sample.int(100, 1))) +
  sumstat_nucleotide_div()
check_model(model)
model # see model
stats <- simulate(model, seed = 15, nsim=50)

# plot
mean_pi <- map_dbl(stats, ~mean(.x$pi))
theta <- map_dbl(stats, "pars", .default="theta")
plot(theta, mean_pi, pch = 19, col = "orange")
abline(lm(mean_pi ~ theta), col = "red2", lty = 3)


# 20 samples from one pop, 
model <- coal_model(sample_size = 10, loci_number = 100, loci_length = 100, ploidy = 2) +
  feat_mutation(.1, model="IFS") +
  sumstat_tajimas_d()
check_model(model)
model
stats <- simulate(model, seed = 15)

plot(density(stats$tajimas_d, na.rm = TRUE), 
     main = "Neutral Distribution of Tajiam's D")

# SFS
model <- coal_model(sample_size = 10, loci_number = 100, loci_length = 300, ploidy = 2) +
  feat_mutation(1, model="IFS") +
  sumstat_sfs()
check_model(model)
model
stats <- simulate(model, seed = 15)

# base option
barplot(stats$sfs / sum(stats$sfs), 
        names.arg = seq_along(stats$sfs), 
        col = 3)

# ggplot option
sfs_out <- map(stats[1], ~.x / sum(.x)) %>% purrr::as_vector() %>% 
  tibble(sfs=., numb=seq_along(stats[[1]]))
ggplot() + geom_col(data=sfs_out, aes(x=numb, y=sfs), fill="cadetblue4")




# Another Pop Gen in R Package --------------------------------------------

#https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.5412