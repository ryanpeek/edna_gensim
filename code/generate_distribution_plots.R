# Generate distributions for simulations
# 2020, R. Peek

## distributions need to contain no negative numbers
## re-scaled in perl script to be between 0-1

# Generate distributions --------------------------------------------------

# number of individuals to simulate distrib
N <- 10

# number of samples to draw
S <- 10000

# Generate Normal Gaussian Distribution -----------------------------------

x_norm <- rnorm(S, mean = 10, sd = 2) # gen distrib
x_norm_s <- sample(x_norm, N) # take a sample
#cat(paste(x_norm_s, collapse=","), file = glue::glue("data/dist_norm_ind{N}.txt"))

# visualize
summary(x_norm_s)
plot(density(x_norm), lwd=2) # look at shape
hist(x_norm, breaks = c(seq(0,20,.5))) # histogram
hist(rep(x_norm_s,250), col = "orange", add=T)
abline(v=10, lwd=5)

# Generate Uniform Distribution -------------------------------------------

x_unif <- runif(S)
x_unif_s <- sample(x_unif, N) # take a sample
#cat(paste(x_unif, collapse=","), file = glue::glue("data/dist_unif_ind{N}.txt"))

# visualize
summary(x_unif)
plot(density(x_unif), lwd=2) # look at shape
hist(x_unif)
hist(rep(x_unif_s, 250), col="orange", add=T)

# Generate Gamm Distrib (Left Skew) ---------------------------------------

# Most commonly used as the distribution of the amount of time until the rth occurrence of an event in a Poisson process. If so, the following apply:
# 1. The number of possible occurrences in any unit of measurement is not limited to a fixed number.
# 2. The occurrences are independent. The number of occurrences in one unit of measurement does not affect the number of occurrences in other units.
# 3. The average number of occurrences must remain the same from unit to unit.

x_gamma <- rgamma(S, shape = 2, rate = 0.5) # left skewed
x_gamma_s <- sample(x_gamma, N) # take a sample
cat(paste(x_gamma_s, collapse=","), file = glue::glue("data/dist_gamma_ind{N}.txt"))

#visualize
summary(x_gamma)
hist(x_gamma, breaks = 100)
hist(rep(x_gamma_s, 250), col="orange", add=T)

# Generate Beta Distrib (Left Skew) ---------------------------------------

# beta distribution is the conjugate prior probability distribution for the Bernoulli, binomial, negative binomial and geometric distributions. The beta distribution is a suitable model for the random behavior of percentages and proportions.
# beta distribution is very flexible and is commonly used to represent variability over a fixed range. beta distribution is used to represent the uncertainty in the probability of occurrence of an event. It is also used to describe empirical data and predict the random behavior of percentages and fractions, as the range of outcomes is typically between 0 and 1.
# alpha (shape1) and beta (shape2): If the parameters are equal, the distribution is symmetrical. If either parameter is 1 and the other parameter is greater than 1, the distribution is J-shaped. If alpha is less than beta, the distribution is said to be positively skewed (most of the values are near the minimum value). If alpha is greater than beta, the distribution is negatively skewed (most of the values are near the maximum value).

x_beta <- rbeta(S, shape1 = 2, shape2 = 10) # left skewed
x_beta_s <- sample(x_beta, N) # take a sample
# x_beta <- rbeta(S, shape1 = 2, shape2 = 2) # normal without tails
# x_beta <- rbeta(S, shape1 = .5, shape2 = .5) # U shaped

# write as single line with no quotes:
#cat(paste(x_beta_s, collapse=","), file = glue::glue("data/dist_beta_ind{N}.txt"))

#visualize
plot(density(x_beta), lwd=2) # look at shape
hist(x_beta, breaks = 100, col="orange")

# Alternate Plotting Methods ----------------------------------------------

summary(x_norm)
hist(x_norm, breaks=100, xlim = c(0,20))
summary(x_gamma)
hist(x_gamma, col="orange", breaks = 100, xlim = c(0,20), add=T)
summary(x_beta)
hist(x_beta, col="purple", breaks = 100, xlim=c(0,20), add=T)

df_norm <- data.frame(x=density(x_norm)$x, y=density(x_norm)$y)
df_unif <- data.frame(x=density(x_unif)$x, y=density(x_unif)$y)
df_gamm <- data.frame(x=density(x_gamma)$x, y=density(x_gamma)$y)
df_beta <- data.frame(x=density(x_beta)$x, y=density(x_beta)$y)


# Make Distribution Quad Plot -----------------------------------------------------------

# make dataframe and combine
library(dplyr)
dfh_norm <- data.frame(x=x_norm, dist="norm")
dfh_unif <- data.frame(x=x_unif, dist="unif")
dfh_gamm <- data.frame(x=x_gamma, dist="gamma")
dfh_beta <- data.frame(x=x_beta, dist="beta")
df_hist <- bind_rows(dfh_norm, dfh_gamm, dfh_unif, dfh_beta)
dim(df_hist)
 
# PLOTS
library(ggplot2)
library(viridis)
library(cowplot)
 
# faceted for gamma/norm
dist1 <- ggplot() +
  geom_histogram(data = df_hist %>% filter(dist %in% c("gamma", "norm")),
                 aes(x = x, fill=dist), binwidth = .5, show.legend = FALSE) +
  scale_fill_manual("Distrib", values = c(viridis(4)[1], viridis(4)[2])) +
  facet_wrap(vars(dist), scales = "free_x") + theme_bw()
 
# faceted for unif/beta
dist2 <- ggplot() +
  geom_histogram(data = df_hist %>% filter(dist %in% c("unif", "beta")),
                 aes(x = x, fill=dist),
                 binwidth = .01, show.legend = FALSE) +
  scale_fill_manual("Distrib.", values = c(viridis(4)[3], viridis(4)[4])) +
  facet_wrap(vars(dist), scales = "free_x") + theme_bw()
 
# merge
pAll <- plot_grid(dist1, dist2, nrow=2)
 
ggsave(pAll, filename = "figs/distributions_faceted.png", width = 8, height = 7, units = "in", dpi = 300)

# NOTES -------------------------------------------------------------------

# add attributes for haplotypes beyond just mean
 # plot other attributes like variance, max, min, CV
 # now we can pick parameters that we know are reasonable:
 # USE:
 # no. loci: 100
 # theta: 10
 # error: 10%

# plot all 3 distributions on one plot for mean
# plot all 3 for max, etc