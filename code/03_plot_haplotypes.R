# PLOT

# Libraries ---------------------------------------------------------------

library(tidyverse)
#library(tidylog)
library(here)

# get data dir
simdir <- "sim003"
#distribution <- "norm"
coverage <- 100


# Read in Data ------------------------------------------------------------

load(glue::glue("results/{simdir}/{simdir}_df_{coverage}x.rda"))

# Tidy Data For Plotting --------------------------------------------------

sim_all_df2 <- sim_all_df %>% 
  group_by(ID, stat) %>% 
  summarize("h_sd"= sd(haplos),
            "CI"= list(mean_cl_normal(haplos) %>% rename(h_mean=y, ci_lo=ymin, ci_up=ymax))) %>% 
  unnest(CI) %>% 
  separate(col = ID, into = c("nInd","nLoci", "theta","coverage", "distrib"), sep = "_", remove = F) %>%
  ungroup() %>% 
  mutate_at(c("nInd","nLoci","theta"), as.numeric) %>% 
  mutate_at(c("distrib", "coverage"), as.factor)

# GGPLOT DOTS -------------------------------------------------------------

# make nice label names:
# lociNames <- c(`100`="No. Loci = 100",
#                `1000`="No. Loci = 1,000")
# thetaNames <- c(`0.1`="theta==0.1 (100)",
#                 `1`="theta==1  (1000)",
#                 `10`="theta==10  (10000)")


(gg1 <- ggplot() + 
    geom_pointrange(data=sim_all_df2 %>% filter(stat!="var"), aes(x=nInd, y=h_mean, ymax=ci_up, ymin=ci_lo, color=stat, shape=distrib), size=1) +
    #geom_line(data=sim_all_df2 %>% filter(stat!="var"), 
    #          aes(x=nInd, y=h_mean, color=stat, group=stat), size=0.5, alpha=0.4) +
    ggthemes::scale_color_colorblind("Metric") +
    facet_grid(distrib ~ .) +
    theme_bw(base_family = "Roboto Condensed") +
    scale_x_continuous(minor_breaks = seq(0,50,2))+
    scale_y_continuous(breaks=seq(0,18,2))+
    labs(y="Number of Haplotypes", x="Number of Individuals",
         title="Simulations of Haplotypes for 100x Coverage",
         caption="based on 1000 replicate simulations for each \n parameter combination, performed in the program *ms*"))

plotly::ggplotly(gg1)

# SAVE OUT
ggsave(glue::glue("figs/{simdir}_haplos_{coverage}x_100loci_faceted.png"), width = 11, height = 8.5, units = "in", dpi=300)
ggsave(glue::glue("figs/{simdir}_haplos_{coverage}x_100loci_faceted.pdf"), 
       device = cairo_pdf, width = 11, height = 8.5, units = "in", dpi=300)  


# BOXPLOTS ----------------------------------------------------------------

ggplot() + 
  geom_boxplot(data=sim_all_df, #%>% filter(stat!="var"), 
               aes(x=nInd, y=haplos, group=nInd, color=stat), outlier.size = 0.5, outlier.alpha = 0.2) +
  facet_grid(distrib~stat) +
  #             labeller= labeller(nLoci = as_labeller(lociNames), 
  #                                theta = as_labeller(thetaNames, label_parsed))) +
  theme_bw(base_family = "Roboto Condensed") +
  ggthemes::scale_color_colorblind("Metric") +
  scale_x_continuous(minor_breaks = seq(0,50,2)) +
  scale_y_continuous(breaks=seq(0,18,2)) +
  guides(fill=FALSE)+
  labs(y="Number of Haplotypes", x="Number of Individuals",
       title="Simulations of Haplotypes at 100x Coverage",
       caption="based on 100 replicate simulations for each \n parameter combination, performed in the program *ms*")

ggsave(glue::glue("figs/{simdir}_haplos_{coverage}x_100loci_faceted.png"), width = 11, height = 8.5, units = "in", dpi=300)
ggsave(glue::glue("figs/{simdir}_haplos_{coverage}x_100loci_faceted.pdf"), device = cairo_pdf, width = 11, height = 8.5, units = "in", dpi=300) 


# Save Data ---------------------------------------------------------------

#dataAll2 <- dataAll[1,] %>% tidyr::pivot_wider(dat)

save(dataAll, file = "data/sim002_err0.1_dataAll_for_10x_100x.rda")

save(mean_haps_df, file= "data/sim002_err0.1_mean_haps_df_10x_100x.rda")



