#!/usr/bin/Rscript

# Generate distributions for simulations
# 2020, R. Peek

## distributions need to contain no negative numbers
## re-scaled in perl script to be between 0-1

# Generate distributions --------------------------------------------------

# get commandArgs from bash script/command line
args <- commandArgs(trailingOnly = TRUE)

# number of individuals to simulate distrib * 2 =gene copies
N <- args[1]
distrib <- args[2]
resultsDir <- args[3]
repX <- args[4]
#print(args)

Ngc <- as.integer(N)
distrib <- as.character(distrib)

#print(args)
#N <- 10

# number of samples to draw
# S <- 10000

# Generate Gaussian Distribution -----------------------------------

if(distrib=="norm"){
        x_norm <- rnorm(Ngc, mean = 10, sd = 2) # gen distrib
        cat(paste(x_norm, collapse=","), file = glue::glue("results/{resultsDir}/dist_norm_ind{N}_{repX}.txt"))
}

# Generate Uniform Distribution -------------------------------------------

if(distrib=="unif"){
        x_unif <- runif(Ngc)
        cat(paste(x_unif, collapse=","), file = glue::glue("results/{resultsDir}/dist_unif_ind{N}_{repX}.txt"))
}

# Generate Equal Distribution -------------------------------------------

if(distrib=="equal"){
        x_equal <- rep(1, Ngc)
        cat(paste(x_equal, collapse=","), file = glue::glue("results/{resultsDir}/dist_equal_ind{N}_{repX}.txt"))
}

# Generate Gamma Distrib (Left Skew) ---------------------------------------

# Most commonly used as the distribution of the amount of time until the rth occurrence of an event in a Poisson process. If so, the following apply:
# 1. The number of possible occurrences in any unit of measurement is not limited to a fixed number.
# 2. The occurrences are independent. The number of occurrences in one unit of measurement does not affect the number of occurrences in other units.
# 3. The average number of occurrences must remain the same from unit to unit.

if(distrib=="gamma"){
        x_gamma <- rgamma(Ngc, shape = 2, rate = 0.5) # left skewed
        cat(paste(x_gamma, collapse=","), file = glue::glue("results/{resultsDir}/dist_gamma_ind{N}_{repX}.txt"))
}

# Generate Extreme Gamma Distrib (Left Skew) ------------------------------

# Extreme left skew largely peaked around 0-1

if(distrib=="gammaX"){
        x_gammaX <- rgamma(Ngc, shape = 0.1, rate = 0.1) # left skewed
        cat(paste(x_gammaX, collapse=","), file = glue::glue("results/{resultsDir}/dist_gammaX_ind{N}_{repX}.txt"))
}

# Generate Beta Distrib (Left Skew) ---------------------------------------

# beta distribution is the conjugate prior probability distribution for the Bernoulli, binomial, negative binomial and geometric distributions. The beta distribution is a suitable model for the random behavior of percentages and proportions.
# beta distribution is very flexible and is commonly used to represent variability over a fixed range. beta distribution is used to represent the uncertainty in the probability of occurrence of an event. It is also used to describe empirical data and predict the random behavior of percentages and fractions, as the range of outcomes is typically between 0 and 1.
# alpha (shape1) and beta (shape2): If the parameters are equal, the distribution is symmetrical. If either parameter is 1 and the other parameter is greater than 1, the distribution is J-shaped. If alpha is less than beta, the distribution is said to be positively skewed (most of the values are near the minimum value). If alpha is greater than beta, the distribution is negatively skewed (most of the values are near the maximum value).

if(distrib=="beta"){
        x_beta <- rbeta(Ngc, shape1 = 2, shape2 = 10) # left skewed
        cat(paste(x_beta, collapse=","), file = glue::glue("results/{resultsDir}/dist_beta_ind{N}_{repX}.txt"))
}


# Generate U-shaped Beta Distrib ---------------------------------------

if(distrib=="betaU"){
        x_betaU <- rbeta(Ngc, shape1 = 0.1, shape2 = 0.1) # left skewed
        cat(paste(x_betaU, collapse=","), file = glue::glue("results/{resultsDir}/dist_betaU_ind{N}_{repX}.txt"))
}

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