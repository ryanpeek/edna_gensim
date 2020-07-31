# Make haplotype plots

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(tidylog)
library(here)

# 00: LOAD EXISTING DATA IF AVAIL. --------------------------------------------------------------

# mean haplo list
load("data/sim001_mean_haps_df_10x_100x.rda") # mean_haps_df
sim001_haps_df <- mean_haps_df
load("data/sim002_err0.1_mean_haps_df_10x_100x.rda")
sim002_haps_df <- mean_haps_df; rm(mean_haps_df)

# raw data list (large)
# load("data/sim001_dataAll_for_10x_100x.rda")
# sim001_dataAll <- dataAll; rm(dataAll)
# load("data/sim002_err0.1_dataAll_for_10x_100x.rda")
# sim002_dataAll <- dataAll; rm(dataAll)

# raw file list
load("data/sim001_files_list.rda")
sim001_files_df <- files_df; rm(files_df)
load("data/sim002_files_list.rda")
sim002_files_df <- files_df; rm(files_df)

# tidy files
load("data/sim_all_mean.rda")
load("data/sim_all_median.rda")

# 01: GET RAW FILES LIST ----------------------------------------------------------


# get data dir
simdir <- "sim001_data"

# get filenames
filenames <- list.files(path = paste0("data/",simdir), pattern = ".out",full.names = F)
head(filenames)

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
  arrange(nLoci, theta, coverage, nInd)

# save this out
save(files_df, file = "data/sim001_files_list.rda")


# get data dir
simdir <- "sim002_data"

# get filenames
filenames <- list.files(path = paste0("data/",simdir), pattern = ".out",full.names = F)
head(filenames)

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
  arrange(nLoci, theta, coverage, nInd)

# save this out
save(files_df, file = "data/sim002_files_list.rda")


# 02: IMPORT FILES  ------------------------------------------

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

# 97.729 sec elapsed

# 03: CALCULATE MEAN HAPLOTYPES per LOCUS -----------------

# make custom function to calculate and pivot
meanHaplos_by_locus <- function(x) {
  dplyr::summarize_all(x, mean) %>% 
    pivot_longer(cols=starts_with("X"), names_to="reps")
}

mediHaplos_by_locus <- function(x) {
  dplyr::summarize_all(x, median) %>% 
    pivot_longer(cols=starts_with("X"), names_to="reps")
}


# ** Calc SIM001 -----------------------------------------------

# Calculate Mean number of Haplotypes per individual/params
sim001_mean_haps <- map(sim001_dataAll$dat, ~meanHaplos_by_locus(.x)) %>% 
  map2_df(., sim001_files_df$sampleID, ~mutate(.x, ID=.y))

# clean a bit and refactor
sim001_mean_haps_df <- sim001_mean_haps %>% 
  separate(col = ID, into = c("nInd","nLoci", "theta","coverage"), sep = "_", remove = F) %>%
  rename(mean_haplos=value) %>% 
  mutate_at(c("nInd","nLoci","theta"), as.numeric) %>%
  mutate(coverage=as.factor(coverage)) %>% 
  mutate(reps = as.integer(gsub(reps, pattern = "^X",replacement = "")))

# Calculate MeDian number of Haplotypes per individual/params
sim001_med_haps <- map(sim001_dataAll$dat, ~mediHaplos_by_locus(.x)) %>% 
  map2_df(., sim001_files_df$sampleID, ~mutate(.x, ID=.y))

# clean a bit and refactor
sim001_med_haps_df <- sim001_med_haps %>% 
  separate(col = ID, into = c("nInd","nLoci", "theta","coverage"), sep = "_", remove = F) %>%
  rename(med_haplos=value) %>% 
  mutate_at(c("nInd","nLoci","theta"), as.numeric) %>%
  mutate(coverage=as.factor(coverage)) %>% 
  mutate(reps = as.integer(gsub(reps, pattern = "^X",replacement = "")))

# ** Calc SIM002 -----------------------------------------------

# Calculate Mean number of Haplotypes per individual/params
sim002_mean_haps <- map(sim002_dataAll$dat, ~meanHaplos_by_locus(.x)) %>% 
  map2_df(., sim002_files_df$sampleID, ~mutate(.x, ID=.y))

# clean a bit and refactor
sim002_mean_haps_df <- sim002_mean_haps %>% 
  separate(col = ID, into = c("nInd","nLoci", "theta","coverage"), sep = "_", remove = F) %>%
  rename(mean_haplos=value) %>% 
  mutate_at(c("nInd","nLoci","theta"), as.numeric) %>%
  mutate(coverage=as.factor(coverage)) %>% 
  mutate(reps = as.integer(gsub(reps, pattern = "^X",replacement = "")))

# Calculate MeDian number of Haplotypes per individual/params
sim002_med_haps <- map(sim002_dataAll$dat, ~mediHaplos_by_locus(.x)) %>% 
  map2_df(., sim002_files_df$sampleID, ~mutate(.x, ID=.y))

# clean a bit and refactor
sim002_med_haps_df <- sim002_med_haps %>% 
  separate(col = ID, into = c("nInd","nLoci", "theta","coverage"), sep = "_", remove = F) %>%
  rename(med_haplos=value) %>% 
  mutate_at(c("nInd","nLoci","theta"), as.numeric) %>%
  mutate(coverage=as.factor(coverage)) %>% 
  mutate(reps = as.integer(gsub(reps, pattern = "^X",replacement = "")))

# 04: COMBINE INTO SINGLE DATASET ------------------------------

## MEAN
sim001_mean_haps_df <- sim001_mean_haps_df %>% 
  mutate(modrun = "error=0")

sim002_mean_haps_df <- sim002_mean_haps_df %>% 
  mutate(modrun = "error=0.10")

## MEDIAN
sim001_med_haps_df <- sim001_med_haps_df %>% 
  mutate(modrun = "error=0")

sim002_med_haps_df <- sim002_med_haps_df %>% 
  mutate(modrun = "error=0.10")

# COMBINE ALL: MEAN
sim_all_mean <- bind_rows(sim001_mean_haps_df, sim002_mean_haps_df) %>% 
  mutate("Model"=as.factor(glue::glue("{coverage}x {modrun}")))
save(sim_all_mean, file = "data/sim_all_mean.rda")

# COMBINE ALL: MEDIAN
sim_all_median <- bind_rows(sim001_med_haps_df, sim002_med_haps_df) %>% 
  mutate("Model"=as.factor(glue::glue("{coverage}x {modrun}")))
save(sim_all_median, file = "data/sim_all_median.rda")

# 05: PLOTS ---------

#  ** Mean Haplos All ------------------------------

sim_all_mean2 <- sim_all_mean %>% 
  group_by(ID, modrun) %>% 
  summarize_at(.vars = vars(mean_haplos), 
               .funs = c("hMax"=max, "hMin"=min, "hMean"=mean)) %>% 
  separate(col = ID, into = c("nInd","nLoci", "theta","coverage"), sep = "_", remove = F) %>%
  ungroup() %>% 
  mutate_at(c("nInd","nLoci","theta"), as.numeric)

# make nice label names:
lociNames <- c(`100`="No. Loci = 100",
               `1000`="No. Loci = 1,000")
thetaNames <- c(`0.1`="theta==0.1 (100)",
                `1`="theta==1  (1000)",
                `10`="theta==10  (10000)")

# now plot
plotly::ggplotly(ggplot() + 
  geom_pointrange(data=sim_all_mean2, aes(x=nInd, y=hMean, ymax=hMax, ymin=hMin, color=modrun, shape=coverage), size=0.5) +
  ggthemes::scale_color_colorblind("Model")+
  scale_shape_discrete("Coverage")+
  facet_grid(theta ~ nLoci, 
             labeller= labeller(nLoci = as_labeller(lociNames), theta = as_labeller(thetaNames, label_parsed))) +
  theme_bw(base_family = "Roboto Condensed") +
  scale_x_continuous(minor_breaks = seq(0,50,2))+
  scale_y_continuous(breaks=seq(0,18,2))+
  labs(y="Mean Number of Haplotypes", x="Number of Individuals",
       title="Simulations of Mean Haplotypes for 10x & 100x Coverage",
       caption="based on 1000 replicate simulations for each \n parameter combination, performed in the program *ms*"))

# SAVE OUT
ggsave("figs/ms_mean_haplotypes_10x100x_50nInd_w_error.png", width = 11, height = 8.5, units = "in", dpi=300)
ggsave("figs/ms_mean_haplotypes_10x100x_50nInd_w_error.pdf", device = cairo_pdf,
       width = 11, height = 8.5, units = "in", dpi=300)  


# ** Mean Haplos 10x ----------------------------------------------------------------

# filter to 10 X
#mean_haps_df <- mean_haps_df %>% filter(coverage=="10")

# make nice label names:
lociNames <- c(`100`="No. Loci = 100",
               `1000`="No. Loci = 1,000")
thetaNames <- c(`0.1`="theta==0.1 (100)",
                `1`="theta==1  (1000)",
                `10`="theta==10  (10000)")

ggplot() + geom_boxplot(data=sim001_haps_df, aes(x=nInd, y=mean_haplos, group=nInd, color=coverage), 
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


# ** Mean Haplos 100x ----------------------------------------------------------------

# filter to 100 X
#mean_haps_df <- mean_haps_df %>% filter(coverage=="100")

# make nice label names:
lociNames <- c(`100`="No. Loci = 100",
               `1000`="No. Loci = 1,000")
thetaNames <- c(`0.1`="theta==0.1 (100)",
                `1`="theta==1  (1000)",
                `10`="theta==10  (10000)")

# make the plots
ggplot() + geom_boxplot(data=sim002_haps_df, aes(x=nInd, y=mean_haplos, group=nInd, color=coverage), 
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



