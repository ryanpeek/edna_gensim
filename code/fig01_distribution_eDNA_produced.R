# make figures of distributions

library(tidyverse)
library(ggfortify)
library(ggthemes)
library(cowplot)
library(scales)
# get colorblind palette
# scales::show_col(ggthemes::colorblind_pal()(8))
colblind_pal <- ggthemes::colorblind_pal()(4)


# Make Distributions ------------------------------------------------------

nnorm <- rnorm(1e6, mean = 5, sd=1.5)
nunif <- runif(1e6, min=-5, max=11)
ngam <- rgamma(1e6, shape = 1.7)

# make dataframe
df <- tibble(norm = nnorm, unif = nunif, gamma=ngam)


# Histogram Plots ---------------------------------------------------------

# geom histogram 
ggplot() + 
  geom_histogram(data=df, aes(gamma), bins = 100) +
  geom_histogram(data=df, aes(unif), fill=colblind_pal[2], bins=100, alpha=0.8) +
  geom_histogram(data=df, aes(norm), fill=colblind_pal[3], bins=100, alpha=0.8) +
  xlim(c(0,10)) +
  scale_y_continuous(labels = percent_format(scale =1e-3)) +
  theme_cowplot() +
  labs(
       y="Proportion of Individuals in Population",
       x="Amount of eDNA in sample over time")

ggsave(filename = "figs/edna_distributions_in_sample_histogram.png", 
       width = 10, height = 8, dpi=300)

# Density Plots -----------------------------------------------------------

ggplot() + 
  geom_density(data=df, aes(gamma), fill="gray50") +
  geom_density(data=df, aes(unif), color=colblind_pal[2],
               fill=colblind_pal[2], alpha=0.8) +
  geom_density(data=df, aes(norm), color=colblind_pal[3],
               fill=colblind_pal[3], alpha=0.8) +
  xlim(c(0,10)) +
  scale_y_continuous(labels = percent_format(scale =100)) +
  theme_cowplot() +
  labs( 
       y="Proportion of Individuals in Population",
       x="Amount of eDNA in sample over time")
ggsave(filename = "figs/edna_distributions_in_sample_density.png", 
       width = 10, height = 8, dpi=300)

# ggdistribution Plots ----------------------------------------------------


# using ggdistribution

# gamma left shifted
pgam <- ggdistribution(dgamma, seq(0, 10, 0.01), shape=1.5, 
                    colour = colblind_pal[2], 
                    fill = colblind_pal[2], alpha=0.4)
# normal dist
pnorm <- ggdistribution(dnorm, seq(0, 10, 0.01), mean = 5, sd=2,
                    colour = colblind_pal[3], linetype = 4,
                    fill = colblind_pal[3], alpha=0.4, p = pgam)
# log normal so left skewed like gamma
# p <- ggdistribution(dlnorm, seq(0, 10, 0.01), 
#                     meanlog = 0.5, colour = 'cyan4', p = p)
# uniform
punif <- ggdistribution(dunif, seq(0, 10, 0.01), min=0, max=10,
                    colour = colblind_pal[4],linetype = 2,
                    fill = colblind_pal[4], alpha=0.2, p=pnorm)

# plot together
(tst <- punif + theme_cowplot() +
    ylab("Proportion of Individuals in Population") +
    xlab("Amount of eDNA in sample over time")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())
ggsave(filename = "figs/edna_distributions_in_sample_ggdistrib_blank.png", 
       width = 10, height = 8, dpi=300)

