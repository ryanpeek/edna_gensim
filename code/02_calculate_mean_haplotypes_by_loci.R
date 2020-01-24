# Make haplotype plots

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(tidylog)
library(here)

# Get Data ----------------------------------------------------------------

# need to read in all the data and collapse:

# get filenames
(filenames <- list.files(path = "data/simulations", pattern = ".out",full.names = TRUE))
# using fs fs::dir_ls(path = "data/simulations/", regexp = "\\.out")


# do some tidying
files_df <- filenames %>% tibble::enframe(name = "fileno", value = "filename") %>%
  separate(col = filename, 
           into = c("run","nInd","nLoci", "theta","coverage", "reps"), 
           sep = "_", remove = F) %>%
  select(-fileno, -run, -reps) %>% 
  mutate(theta = as.numeric(gsub(x=theta, pattern = "^t", replacement = "")),
         nInd = as.numeric(nInd),
         nLoci = as.numeric(nLoci),
         coverage = as.factor(coverage),
         sampleID = glue::glue("{nInd}_{nLoci}_{theta}_{coverage}")) %>% ## need to add coverage in now
  # arrange by 
  arrange(nLoci, theta, nInd)

# can check against the original run list, should be 600 runs
# ms_params <- read_delim(file = "data/ms_params_50inds_2020_01_22.txt", delim = " ", 
#                         col_names = c("nInd", "nLoci", "ll", "coverage")) 

# Read a Single File ------------------------------------------------------


dat <- read_tsv(file = files_df$filename[[8]], col_names = FALSE)

# plot # of sims by mean # of haplotypes per locus
dat_means <- dat %>% 
  summarize_all(mean, na.rm=TRUE) %>% 
  pivot_longer(cols=starts_with("X"), names_to = "sims")

hist(dat_means$value)

# from one simulation
# number of loci vs. # of haplotypes
hist(dat$X77)


# Read Multiple Files with Purrr ------------------------------------------

# reads all in at once, keeps data as list col
dataAll <- files_df %>% 
  mutate(dat = map(files_df$filename, ~read_tsv(.x, col_names = FALSE)))


# Calculate Mean Number of Haplotypes per Locus ---------------------------

# make custom function to calculate and pivot
meanHaplos_by_locus <- function(x) {
  summarize_all(x, mean) %>% 
    pivot_longer(cols=starts_with("X"), names_to="reps")
}

# Calculate Mean number of Haplotypes per individual/params
mean_haps <- map(dataAll$dat, ~meanHaplos_by_locus(.x)) %>% 
  map2_df(., files_df$sampleID, ~mutate(.x, ID=.y))

# clean a bit
mean_haps_df <- mean_haps %>% 
  separate(col = ID, into = c("nInd","nLoci", "theta","coverage"), sep = "_", remove = F) %>%
  rename(mean_haplos=value) %>% 
  mutate_at(c("nInd","nLoci","theta"), as.numeric) %>%
  mutate(coverage=as.factor(coverage)) %>% 
  mutate(reps = as.integer(gsub(reps, pattern = "^X",replacement = "")))


# PLOT --------------------------------------------------------------------

# number of individuals on X

# make nice label names:
lociNames <- c(`100`="No. Loci = 100",
               `1000`="No. Loci = 1,000")
thetaNames <- c(`0.1`="theta==0.1 (100)",
                `1`="theta==1  (1000)",
                `10`="theta==10  (10000)")

ggplot() + geom_boxplot(data=mean_haps_df, aes(x=nInd, y=mean_haplos, group=nInd), 
                        #outlier.shape = NA, 
                        outlier.size = 0.5, outlier.alpha = 0.2) +
  facet_grid(theta ~ nLoci, labeller= labeller(nLoci = as_labeller(lociNames),
                                               theta = as_labeller(thetaNames, label_parsed))) +
  theme_bw(base_family = "Roboto Condensed") +
  scale_x_continuous(minor_breaks = seq(0,50,2))+
  scale_y_continuous(breaks=seq(0,18,2))+
  labs(y="Mean Number of Haplotypes", x="Number of Individuals",
       title="Simulations of Mean Haplotypes at 100x Coverage",
       caption="based on 1000 replicate simulations for each \n parameter combination, performed in the program *ms*")

ggsave("figs/ms_mean_haplotypes_100x_50nInd.png", width = 11, height = 8.5, units = "in", dpi=300)
ggsave("figs/ms_mean_haplotypes_100x_50nInd.pdf", device = cairo_pdf,
       width = 11, height = 8.5, units = "in", dpi=300)  

