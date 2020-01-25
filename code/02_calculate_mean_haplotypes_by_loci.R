# Make haplotype plots

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(tidylog)
library(here)

# LOAD EXISTING DATA IF AVAIL. --------------------------------------------------------------

load("data/sim001_mean_haps_df_10x_100x.rda") # mean_haps_df

#load("data/sim002_err0.1_mean_haps_df_10x_100x.rda")
#sim002_haps_df <- mean_haps_df
#rm(mean_haps_df)

# Get Data ----------------------------------------------------------------

# get data dir
simdir <- "sim002_data"

# get filenames
(filenames <- list.files(path = paste0("data/",simdir), pattern = ".out",full.names = F))

# using fs fs::dir_ls(path = paste0("data/",simdir, "/"), regexp = "\\.out")

# do some tidying
files_df <- filenames %>% tibble::enframe(name = "fileno", value = "filename") %>%
  separate(col = filename, 
           into = c("run","nInd","nLoci", "theta",
                    "coverage", "error", "hapThresh", "reps"), 
           sep = "_", remove = F) %>%
  select(-fileno, -run, -reps, -hapThresh, -error) %>% 
  mutate(theta = as.numeric(gsub(x=theta, pattern = "^t", replacement = "")),
         nInd = as.numeric(nInd),
         nLoci = as.numeric(nLoci),
         coverage = as.factor(coverage),
         sampleID = glue::glue("{nInd}_{nLoci}_{theta}_{coverage}")) %>% 
  arrange(nLoci, theta, nInd)

# Read Multiple Files with Purrr ------------------------------------------

# reads all in at once, keeps data as list col

library(tictoc) # to time stuff
library(furrr) # use furrr and work in parallel
#plan(multiprocess) this does nothing

# now run in parallel
tic()
dataAll <- files_df %>% 
  #slice(1:100) %>% 
  mutate(dat = furrr::future_imap(paste0("data/", simdir,"/",filename),
                                  ~read_tsv(.x, col_names = FALSE),
                                  .progress = TRUE))
toc()


### try without parallel
# this is normal method (>2x slower than furrr)
# tic()
# dataAll <- files_df %>% 
#   slice(1:100) %>% 
#   mutate(dat = map(filename, ~read_tsv(.x, col_names = FALSE)))
# toc()


# Calculate Mean Number of Haplotypes per Locus ---------------------------

# make custom function to calculate and pivot
meanHaplos_by_locus <- function(x) {
  summarize_all(x, mean) %>% 
    pivot_longer(cols=starts_with("X"), names_to="reps")
}

# Calculate Mean number of Haplotypes per individual/params
mean_haps <- map(dataAll$dat, ~meanHaplos_by_locus(.x)) %>% 
  map2_df(., files_df$sampleID, ~mutate(.x, ID=.y))

# clean a bit and refactor
mean_haps_df <- mean_haps %>% 
  separate(col = ID, into = c("nInd","nLoci", "theta","coverage"), sep = "_", remove = F) %>%
  rename(mean_haplos=value) %>% 
  mutate_at(c("nInd","nLoci","theta"), as.numeric) %>%
  mutate(coverage=as.factor(coverage)) %>% 
  mutate(reps = as.integer(gsub(reps, pattern = "^X",replacement = "")))



# MEAN NUMBER HAPLOTYPES: ALL ----------------------------------------------------------------

# make nice label names:
lociNames <- c(`100`="No. Loci = 100",
               `1000`="No. Loci = 1,000")
thetaNames <- c(`0.1`="theta==0.1 (100)",
                `1`="theta==1  (1000)",
                `10`="theta==10  (10000)")

ggplot() + 
  geom_boxplot(data=mean_haps_df %>% 
                 filter(coverage=="10"), 
               aes(x=nInd, y=mean_haplos, group=nInd, color=coverage),
               #color="mediumpurple4",
               outlier.shape = NA, # make points null
               outlier.size = 0.5, outlier.alpha = 0.2) +
  geom_boxplot(data=mean_haps_df %>% 
                 filter(coverage=="100"), 
               aes(x=nInd, y=mean_haplos, group=nInd, color=coverage),
               #color="seagreen3",
               outlier.shape = NA, # make points null
               outlier.size = 0.5, outlier.alpha = 0.2) +
  scale_color_viridis_d("Coverage")+
  facet_grid(theta ~ nLoci, 
             labeller= labeller(nLoci = as_labeller(lociNames), theta = as_labeller(thetaNames, label_parsed))) +
  theme_bw(base_family = "Roboto Condensed") +
  scale_x_continuous(minor_breaks = seq(0,50,2))+
  scale_y_continuous(breaks=seq(0,18,2))+
  labs(y="Mean Number of Haplotypes", x="Number of Individuals",
       title="Simulations of Mean Haplotypes for 10x & 100x Coverage",
       caption="based on 1000 replicate simulations for each \n parameter combination, performed in the program *ms*")

ggsave("figs/ms_mean_haplotypes_10x100x_50nInd.png", width = 11, height = 8.5, units = "in", dpi=300)
ggsave("figs/ms_mean_haplotypes_10x100x_50nInd.pdf", device = cairo_pdf,
       width = 11, height = 8.5, units = "in", dpi=300)  


# MEAN NUMBER HAPLOTYPES: 10X ----------------------------------------------------------------

# filter to 10 X
#mean_haps_df <- mean_haps_df %>% filter(coverage=="10")

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
       title="Simulations of Mean Haplotypes at 10x Coverage",
       caption="based on 1000 replicate simulations for each \n parameter combination, performed in the program *ms*")

ggsave("figs/ms_mean_haplotypes_10x_50nInd.png", width = 11, height = 8.5, units = "in", dpi=300)
ggsave("figs/ms_mean_haplotypes_10x_50nInd.pdf", device = cairo_pdf,
       width = 11, height = 8.5, units = "in", dpi=300)  


# MEAN NUMBER HAPLOTYPES: 100X ----------------------------------------------------------------

# filter to 100 X
#mean_haps_df <- mean_haps_df %>% filter(coverage=="100")

# make nice label names:
lociNames <- c(`100`="No. Loci = 100",
               `1000`="No. Loci = 1,000")
thetaNames <- c(`0.1`="theta==0.1 (100)",
                `1`="theta==1  (1000)",
                `10`="theta==10  (10000)")

# make the plots
ggplot() + geom_boxplot(data=mean_haps_df, aes(x=nInd, y=mean_haplos, group=nInd), 
                        #outlier.shape = NA, 
                        outlier.size = 0.5, outlier.alpha = 0.2) +
  facet_grid(theta ~ nLoci,
             labeller= labeller(nLoci = as_labeller(lociNames), 
                                theta = as_labeller(thetaNames, label_parsed))) +
  theme_bw(base_family = "Roboto Condensed") +
  scale_x_continuous(minor_breaks = seq(0,50,2))+
  scale_y_continuous(breaks=seq(0,18,2))+
  labs(y="Mean Number of Haplotypes", x="Number of Individuals",
       title="Simulations of Mean Haplotypes at 100x Coverage",
       caption="based on 1000 replicate simulations for each \n parameter combination, performed in the program *ms*")

ggsave("figs/ms_mean_haplotypes_100x_50nInd_w_0.1_error.png", width = 11, height = 8.5, units = "in", dpi=300)
ggsave("figs/ms_mean_haplotypes_100x_50nInd_w_0.1_error.pdf", device = cairo_pdf,
       width = 11, height = 8.5, units = "in", dpi=300)  


# Save Data ---------------------------------------------------------------

#dataAll2 <- dataAll[1,] %>% tidyr::pivot_wider(dat)

save(dataAll, file = "data/sim002_err0.1_dataAll_for_10x_100x.rda")

save(mean_haps_df, file= "data/sim002_err0.1_mean_haps_df_10x_100x.rda")


# Read a Single File ------------------------------------------------------


# dat <- read_tsv(file = files_df$filename[[8]], col_names = FALSE)
# 
# # plot # of sims by mean # of haplotypes per locus
# dat_means <- dat %>% 
#   summarize_all(mean, na.rm=TRUE) %>% 
#   pivot_longer(cols=starts_with("X"), names_to = "sims")
# 
# hist(dat_means$value)
# 
# # from one simulation
# # number of loci vs. # of haplotypes
# hist(dat$X77)

