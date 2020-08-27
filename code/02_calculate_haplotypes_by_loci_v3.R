# Make haplotype plots

# Libraries ---------------------------------------------------------------

library(tidyverse)
#library(tidylog)
library(here)

# 01: GET RAW FILES LIST ----------------------------------------------------------

# get data dir
simdir <- "sim003"
#distribution <- "norm"
coverage <- 100

# get filenames
filenames <- list.files(path = paste0("results/",simdir), pattern = ".out",full.names = F)
#head(filenames)

# do some tidying
files_df <- filenames %>% tibble::enframe(name = "fileno", value = "filename") %>%
  separate(col = filename, 
           into = c("run","nInd","nLoci", "theta",
                    "coverage", "error", "hapThresh",
                    "distrib", "reps"), 
           sep = "_", remove = F) %>%
  select(-fileno, -run, -reps, -hapThresh, -error) %>% 
  mutate(theta = as.numeric(gsub(x=theta, pattern = "^t", replacement = "")),
         nInd = as.numeric(nInd),
         nLoci = as.numeric(nLoci),
         coverage = as.factor(coverage),
         distrib = as.factor(distrib),
         sampleID = glue::glue("{nInd}_{nLoci}_{theta}_{coverage}_{distrib}")) %>% 
  arrange(nLoci, theta, coverage, nInd, distrib)

#View(files_df)

# save this out
save(files_df, file = glue::glue("results/{simdir}/{simdir}_files_list_{coverage}x.rda"))


# 02: IMPORT FILES  ------------------------------------------

# reads all in at once, keeps data as list col
library(tictoc) # to time stuff
library(furrr) # use furrr and work in parallel

# now run in parallel
#tic()
dataAll <- files_df %>% 
  #slice(1:100) %>% 
  mutate(dat = furrr::future_imap(glue::glue("results/{simdir}/{filename}"),
                                  ~read_tsv(.x, col_names = FALSE),
                                  .progress = TRUE))
#toc()

# 03: CALCULATE MEAN HAPLOTYPES per LOCUS -----------------


# make custom function to aggregate and pivot
f_agg_haps_by_locus <- function(x, stat) {
  dplyr::summarize_all(x, stat) %>% 
    pivot_longer(cols=starts_with("X"), names_to="reps")
}

# Function to make a dataframe and refactor
f_haps_df <- function(statname, df){
  df %>% 
    separate(col = ID, into = c("nInd","nLoci", "theta","coverage", "distrib"), sep = "_", remove = F) %>%
    rename(haplos=value) %>%
    mutate(stat=statname) %>% 
    mutate_at(c("nInd","nLoci","theta"), as.numeric) %>%
    mutate_at(c("distrib", "coverage"), as.factor) %>%  
    mutate(reps = as.integer(gsub(reps, pattern = "^X",replacement = "")))
}


# ** Calculate Mean -------------------------------------------

# Calculate Mean number of Haplotypes per individual/params
sim003_mean_df <- map(dataAll$dat, ~f_agg_haps_by_locus(.x, "mean")) %>% 
  map2_df(., files_df$sampleID, ~mutate(.x, ID=.y))

# apply function: MEAN
sim003_mean_df <- f_haps_df("mean", sim003_mean_df)

# ** Calculate Var -------------------------------------------

# Calculate Variance of number of Haplotypes per individual/params
sim003_var_df <- map(dataAll$dat, ~f_agg_haps_by_locus(.x, "var")) %>% 
  map2_df(., files_df$sampleID, ~mutate(.x, ID=.y))
sim003_var_df <- f_haps_df("var", sim003_var_df)

# ** Calculate Max -------------------------------------------

# Calculate Variance of number of Haplotypes per individual/params
sim003_max_df <- map(dataAll$dat, ~f_agg_haps_by_locus(.x, "max")) %>% 
  map2_df(., files_df$sampleID, ~mutate(.x, ID=.y))
sim003_max_df <- f_haps_df("max", sim003_max_df)

# ** Calculate Min -------------------------------------------

# Calculate Variance of number of Haplotypes per individual/params
sim003_min_df <- map(dataAll$dat, ~f_agg_haps_by_locus(.x, "min")) %>% 
  map2_df(., files_df$sampleID, ~mutate(.x, ID=.y))
sim003_min_df <- f_haps_df("min", sim003_min_df)

# 04: COMBINE INTO SINGLE DATASET ------------------------------

# COMBINE ALL:
sim_all_df <- bind_rows(sim003_mean_df, sim003_max_df, sim003_min_df, sim003_var_df) %>% 
  mutate("Model"=as.factor(glue::glue("{coverage}x {distrib}")))

# save out
save(sim_all_df, file = glue::glue("results/{simdir}/{simdir}_df_{coverage}x.rda"))


