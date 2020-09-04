# PLOT

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(ggtext)
library(glue)

# get data dir
simdir <- "sim005"
reps <- 1000


# Read in Data ------------------------------------------------------------

# read in the by locus
load(glue::glue("results/{simdir}/{simdir}_locus_df_{reps}reps.rda"))

# filter to a stat
sim_loc_df <- filter(sim_loc_df, stat=="mean") %>% 
  # make cover an ordered factor:
  mutate(coverage=forcats::fct_relevel(coverage, "10","25","50","100"))

# Tidy Data For Plotting --------------------------------------------------

sim_loc_ci <- sim_loc_df %>%  
  group_by(sampleID, stat) %>%
  summarize("h_sd"= sd(haplos),
            "h_max"= max(haplos),
            "h_min"= min(haplos),
            # calc 95% CI (based on Hmisc::smean.cl.normal())
            "CI"= list(mean_cl_normal(haplos) %>% rename(ci_mean=y, ci_lo=ymin, ci_up=ymax)),
            
            # calc 95% CI (based on Hmisc::smean.cl.normal())
            "CI_med"= list(median_hilow(haplos) %>% rename(ci_med=y, ci_mlo=ymin, ci_mup=ymax))) %>%
  unnest(c(CI, CI_med)) %>%
  separate(col = sampleID, into = c("nInd","nLoci", "theta","coverage", "distrib"), sep = "_", remove = F) %>%
  ungroup() %>%
  mutate_at(c("nInd","nLoci","theta"), as.numeric) %>%
  mutate_at(c("distrib", "coverage"), as.factor) %>% 
  mutate(coverage=forcats::fct_relevel(coverage, "10","25","50","100"))

# NOTCHED BOXPLOTS ------------------------------------------------------------

#notched box plot, the notches extend 1.58 * IQR / sqrt(n). This gives a roughly 95% confidence interval for comparing medians. See McGill et al. (1978) 

# make nice label names:
lociNames <- c(`100`="No. Loci = 100",
               `250`="No. Loci = 250",
               `500`="No. Loci = 500",
               `1000`="No. Loci = 1,000")

thetaNames <- c(`1`="theta==1 (1000)",`2.5`="theta==2.5",
                `5`="theta==5", `10`="theta==10  (10000)")


ggplot() + 
  # geom_boxplot(data=sim_loc_df %>% filter(coverage==10),
  #              aes(x=nInd, y=haplos, group=nInd, fill=coverage),
  #              notch=TRUE, outlier.size = 0.5, outlier.alpha = 0.2) +
  geom_boxplot(data=sim_loc_df %>% filter(coverage!=10),
               aes(x=factor(nInd), y=haplos, fill=coverage),
               notch=TRUE, outlier.size = 0.5, outlier.alpha = 0.2) +
  # geom_boxplot(data=sim_loc_df %>% filter(coverage==50),
  #              aes(x=nInd, y=haplos, group=nInd, fill=coverage), 
  #              notch=TRUE, outlier.size = 0.5, outlier.alpha = 0.2) +
  # geom_boxplot(data=sim_loc_df %>% filter(coverage==100),
  #              aes(x=nInd, y=haplos, group=nInd, fill=coverage),
  #              notch=TRUE, outlier.size = 0.5, outlier.alpha = 0.2) +
  facet_grid(theta~nLoci, scales = "free_y",
             labeller= labeller(nLoci = as_labeller(lociNames),
                                theta = as_labeller(thetaNames, label_parsed))) +
  theme_bw(base_family = "Roboto") +
  ggthemes::scale_fill_colorblind("Coverage") +
  theme(
    plot.caption = ggtext::element_markdown()) +
  #guides(color=FALSE)+
  labs(y="Mean No. Haplotypes", x="Number of Individuals",
       title="eDNA Simulations: Haplotype detection across sequencing parameters",
       subtitle="Using a gamma distribution for individual shedding probabilities",
       caption="Based on 1000 replicate simulations for each parameter combination,<br>performed in the program <i>ms</i> with an error rate of 0.1 and minimum haplotype detected=2")



# save
ggsave(glue::glue("figs/{simdir}_haplos_boxplot_faceted.png"), width = 11, height = 8.5, units = "in", dpi=300)
ggsave(glue::glue("figs/{simdir}_haplos_boxplot_faceted.pdf"), device = cairo_pdf, width = 11, height = 8.5, units = "in", dpi=300) 


# mean w 95 CI ------------------------------------------------------------

# make nice label names:
lociNames <- c(`100`="No. Loci = 100",
               `250`="No. Loci = 250",
               `500`="No. Loci = 500",
               `1000`="No. Loci = 1,000")
thetaNames <- c(`1`="theta==1 (1000)",`2.5`="theta==2.5",
                `5`="theta==5", `10`="theta==10  (10000)")


ggplot() + 

  # MEAN
  geom_ribbon(data=sim_loc_ci, aes(x=nInd, ymin=ci_mlo, ymax=ci_mup, fill=coverage), alpha = 0.5) +
  geom_line(data=sim_loc_ci, aes(x=nInd, y=ci_med, color=coverage), alpha = 0.9, lwd=0.5) +
  geom_point(data=sim_loc_ci, aes(x=nInd, y=ci_med, group=nInd, fill=coverage),
             pch=21, alpha = 0.7) +

  # formatting
  facet_grid(theta~nLoci, scales = "free_y",
             labeller= labeller(nLoci = as_labeller(lociNames),
                                theta = as_labeller(thetaNames, label_parsed))) +
  theme_bw(base_family = "Roboto") +
  ggthemes::scale_fill_colorblind("Coverage") +
  ggthemes::scale_color_colorblind("Coverage") +
  theme(
    plot.caption = ggtext::element_markdown()) +
  #guides(color=FALSE)+
  labs(y="mean no. of haplotypes", x="Number of Individuals",
       title="eDNA Simulations: Haplotype detection across sequencing parameters",
       subtitle="Using a gamma distribution for individual shedding probabilities, with 95% CI",
       caption="Based on 1000 replicate simulations for each parameter combination,<br>performed in the program <i>ms</i> with an error rate of 0.1 and minimum haplotype detected=2")

# save
ggsave(glue::glue("figs/{simdir}_haplos_ci_lines_faceted.png"), width = 11, height = 8.5, units = "in", dpi=300)
ggsave(glue::glue("figs/{simdir}_haplos_ci_lines_faceted.pdf"), device = cairo_pdf, width = 11, height = 8.5, units = "in", dpi=300) 
