# PLOT

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(ggtext)
#library(here)

# get data dir
simdir <- "sim005"
reps <- 1000


# Read in Data ------------------------------------------------------------

# read in the by replicate
load(glue::glue("results/{simdir}/{simdir}_rep_df_{reps}reps.rda"))

# read in the by locus
load(glue::glue("results/{simdir}/{simdir}_locus_df_{reps}reps.rda"))

# Tidy Data For Plotting --------------------------------------------------

# sim_loc_ci <- sim_loc_df %>%
#   group_by(sampleID, stat) %>%
#   summarize("h_sd"= sd(haplos),
#             "CI"= list(mean_cl_normal(haplos) %>% rename(h_mean=y, ci_lo=ymin, ci_up=ymax))) %>%
#   unnest(CI) %>%
#   separate(col = sampleID, into = c("nInd","nLoci", "theta","coverage", "distrib"), sep = "_", remove = F) %>%
#   ungroup() %>%
#   mutate_at(c("nInd","nLoci","theta"), as.numeric) %>%
#   mutate_at(c("distrib", "coverage"), as.factor)
# 
# sim_rep_ci <- sim_rep_df %>%
#   group_by(sampleID, stat) %>%
#   summarize("h_sd"= sd(haplos),
#             "CI"= list(mean_cl_normal(haplos) %>% rename(h_mean=y, ci_lo=ymin, ci_up=ymax))) %>%
#   unnest(CI) %>%
#   separate(col = sampleID, into = c("nInd","nLoci", "theta","coverage", "distrib"), sep = "_", remove = F) %>%
#   ungroup() %>%
#   mutate_at(c("nInd","nLoci","theta"), as.numeric) %>%
#   mutate_at(c("distrib", "coverage"), as.factor)


# NOTCHED BOXPLOTS ------------------------------------------------------------

#notched box plot, the notches extend 1.58 * IQR / sqrt(n). This gives a roughly 95% confidence interval for comparing medians. See McGill et al. (1978) 

# make nice label names:
lociNames <- c(`100`="No. Loci = 100",
               `250`="No. Loci = 250",
               `500`="No. Loci = 500",
               `1000`="No. Loci = 1,000")
thetaNames <- c(`1`="theta=1 (1000)",`2.5`="theta=2.5",
                `5`="theta=5", `10`="theta=10  (10000)")



ggplot() + 
  #geom_jitter(data=sim_loc_df %>% filter(stat=="mean"), 
  #            aes(x=nInd, y=haplos, group=nInd, fill=stat), size=0.2, alpha = 0.2) +
  geom_boxplot(data=sim_loc_df,
               aes(x=nInd, y=haplos, group=nInd, color=as.factor(coverage)), 
               notch=TRUE, outlier.size = 0.5, outlier.alpha = 0.2) +
  
  facet_grid(theta~nLoci+coverage) +
  
              labeller= labeller(nLoci = as_labeller(lociNames),
                                 theta = as_labeller(thetaNames, label_parsed)) +
  theme_bw(base_family = "Roboto") +
  ggthemes::scale_fill_colorblind("Coverage") +
  theme(
    plot.caption = ggtext::element_markdown()) +
  #scale_x_continuous(minor_breaks = seq(0,50,2)) +
  #scale_y_continuous(breaks=seq(0,18,2)) +
  guides(color=FALSE)+
  labs(y="Number of Haplotypes", x="Number of Individuals",
       title="By Loci: Simulations of Haplotypes",
       caption="Based on 1000 replicate simulations for each parameter combination,<br>performed in the program <i>ms</i>")

# save
ggsave(glue::glue("figs/{simdir}_haplos_{coverage}x_100loci_boxplot_faceted.png"), width = 11, height = 8.5, units = "in", dpi=300)
ggsave(glue::glue("figs/{simdir}_haplos_{coverage}x_100loci_boxplot_faceted.pdf"), device = cairo_pdf, width = 11, height = 8.5, units = "in", dpi=300) 


# Notched Box of Variance -------------------------------------------------


ggplot() + 
  # geom_jitter(data=sim_loc_df %>% filter(stat=="var"), 
  #             aes(x=nInd, y=haplos, group=nInd, fill=stat), size=0.2, alpha = 0.4) +
  geom_boxplot(data=sim_loc_df %>% filter(stat=="var"), width=1.8,
               aes(x=nInd, y=haplos, group=nInd, fill=stat), notch = TRUE,
               outlier.size = 0.5, outlier.alpha = 0.2) +
  facet_grid(~distrib) +
  theme_bw(base_family = "Roboto") +
  scale_color_manual("Stat", values = "seagreen") +
  scale_fill_manual("Stat", values = "seagreen") +
  theme(plot.caption = ggtext::element_markdown()) +
  #scale_y_continuous(breaks=seq(0,30,3)) +
  guides(color=FALSE) +
  labs(y="Haplotype Variance", x="Number of Individuals",
       title="By Loci: Simulations of Haplotypes at 100x Coverage",
       caption="Based on 1000 replicate simulations for each parameter combination,<br>performed in the program *ms*, using base parameters: **\u0398=10**, **loci-length=10kb**, **100 loci**")

# save
ggsave(glue::glue("figs/{simdir}_haplos_{coverage}x_100loci_boxplot_variance.png"), width = 11, height = 8.5, units = "in", dpi=300)
ggsave(glue::glue("figs/{simdir}_haplos_{coverage}x_100loci_boxplot_variance.pdf"), device = cairo_pdf, width = 11, height = 8.5, units = "in", dpi=300) 
