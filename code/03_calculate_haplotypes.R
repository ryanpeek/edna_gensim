# Make haplotype plots

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(glue)
#library(tidylog)
#library(here)


# Get Data ----------------------------------------------------------------

# pull all recent runs from the server (do once)
system("sh code/02_rclone_get_sims.sh")

# 01: GET RAW FILES LIST ----------------------------------------------------------

# get data dir
simdir <- "sim005"
distribution <- "gamma"
#coverage <- 100

# get filenames
filepaths <- fs::dir_ls(path=glue::glue("~/Box Sync/PROJECTS/SERDP/{simdir}"), glob = "*run*out")

# filenames <- list.files(path = glue::glue("~/Box Sync/PROJECTS/SERDP/{simdir}"), pattern = ".out",full.names = F)

# do some tidying
files_df <- filepaths %>% tibble::enframe(name = "fileno", value = "filename") %>%
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
  arrange(nLoci, theta, coverage, nInd, distrib) %>% 
  mutate(files = basename(filename))

View(files_df)

# save this out
save(files_df, file = glue::glue("results/{simdir}/{simdir}_files_list_{distribution}.rda"))


# 02: IMPORT FILES  ------------------------------------------

# reads all in at once, keeps data as list col
library(tictoc) # to time stuff
library(furrr) # use furrr and work in parallel

# now run in parallel
tic()
dataAll <- files_df %>% 
  #slice(1:100) %>% 
  mutate(dat = furrr::future_imap(glue::glue("{filename}"),
                                  ~read_tsv(.x, col_names = FALSE),
                                  .progress = TRUE))
toc()

#save(dataAll, file = glue("results/{simdir}/{simdir}_dataAll.rda"))
# this is big = 114 MB
library(pryr)
object_size(dataAll)
mem_used()

# 03: CALCULATE HAPLOTYPES  -------------------------------------

# calculate by LOCI (columnwise, all replicates for one loci)
f_agg_haps_by_locus <- function(x, stat) {
  dplyr::summarize_all(x, stat) %>% 
    pivot_longer(cols=starts_with("X"), names_to="reps")
}

# calculate by replicate (so rowwise, all loci for one replicate)
f_agg_haps_by_rep <- function(x, stat) {
  x %>% unnest(dat) %>% rowwise %>% 
    group_by(sampleID) %>% summarize_at(vars(starts_with("X")), stat) %>% 
    pivot_longer(cols=starts_with("X"), names_to="reps")
}

# Function to make a dataframe and refactor
f_haps_df <- function(statname, df){
  df %>% 
    separate(col = sampleID, into = c("nInd","nLoci", "theta","coverage", "distrib"), sep = "_", remove = F) %>%
    rename(haplos=value) %>%
    mutate(stat=statname) %>% 
    mutate_at(c("nInd","nLoci","theta"), as.numeric) %>%
    mutate_at(c("distrib", "coverage"), as.factor) %>%  
    mutate(reps = as.integer(gsub(reps, pattern = "^X",replacement = "")))
}


# ** Calculate Mean -------------------------------------------

# BY LOCI: calculate haplos per loci (over all replicates)
tic()
loc_mean_df <- map(dataAll$dat, ~f_agg_haps_by_locus(.x, "mean")) %>% 
  map2_df(., files_df$sampleID, ~mutate(.x, sampleID=.y))
toc()

loc_mean_df <- f_haps_df("mean", loc_mean_df) # split out wide

# BY REPLICATE: calculate mean haplos per replicate (over all X loci)
#rep_mean_df <- f_agg_haps_by_rep(dataAll, "mean")
#rep_mean_df <- f_haps_df("mean", rep_mean_df)

# ** Calculate Median -------------------------------------------

# BY LOCI: calculate haplos per loci (over all replicates)
# tic()
# loc_median_df <- map(dataAll$dat, ~f_agg_haps_by_locus(.x, "median")) %>% 
#   map2_df(., files_df$sampleID, ~mutate(.x, sampleID=.y))
# toc()
# 
# loc_median_df <- f_haps_df("median", loc_median_df) # split out wide

# BY REPLICATE: calculate mean haplos per replicate (over all X loci)
#rep_median_df <- f_agg_haps_by_rep(dataAll, "median")
#rep_median_df <- f_haps_df("median", rep_median_df)

# # ** Calculate Var -------------------------------------------

# # BY LOCI: calculate haplos per loci (over all replicates)
# loc_var_df <- map(dataAll$dat, ~f_agg_haps_by_locus(.x, "var")) %>% 
#   map2_df(., files_df$sampleID, ~mutate(.x, sampleID=.y))
# loc_var_df <- f_haps_df("var", loc_var_df)
# 
# # BY REPLICATE: calculate haplos per replicate (over all X loci)
# rep_var_df <- f_agg_haps_by_rep(dataAll, "var")
# rep_var_df <- f_haps_df("var", rep_var_df)
# 
# 
# # ** Calculate Max -------------------------------------------
# 
# # BY LOCI: calculate haplos per loci (over all replicates)
# loc_max_df <- map(dataAll$dat, ~f_agg_haps_by_locus(.x, "max")) %>% 
#   map2_df(., files_df$sampleID, ~mutate(.x, sampleID=.y))
# loc_max_df <- f_haps_df("max", loc_max_df)
# 
# # BY REPLICATE: calculate haplos per replicate (over all X loci)
# rep_max_df <- f_agg_haps_by_rep(dataAll, "max")
# rep_max_df <- f_haps_df("max", rep_max_df)
# 
# # ** Calculate Min -------------------------------------------
# 
# # BY LOCI: calculate haplos per loci (over all replicates)
# loc_min_df <- map(dataAll$dat, ~f_agg_haps_by_locus(.x, "min")) %>% 
#   map2_df(., files_df$sampleID, ~mutate(.x, sampleID=.y))
# loc_min_df <- f_haps_df("min", loc_min_df)
# 
# # BY REPLICATE: calculate haplos per replicate (over all X loci)
# rep_min_df <- f_agg_haps_by_rep(dataAll, "min")
# rep_min_df <- f_haps_df("min", rep_min_df)



# 04: COMBINE INTO SINGLE DATASET ------------------------------

# COMBINE ALL BY LOCUS:
sim_loc_df <- bind_rows(loc_mean_df) %>% 
                        #loc_max_df, loc_min_df, loc_var_df) %>% 
  mutate("Model"=as.factor(glue::glue("t{theta}_{coverage}x {distrib} by Locus")))

object_size(sim_loc_df)

# save out
save(sim_loc_df, file = glue::glue("results/{simdir}/{simdir}_locus_df_1000reps.rda"))



# COMBINE ALL BY REPS:
# sim_rep_df <- bind_rows(rep_mean_df) %>% 
#                         #rep_max_df, rep_min_df, rep_var_df) %>% 
#   mutate("Model"=as.factor(glue::glue("t{theta}_{coverage}x {distrib} by Rep")))
# 
# # save out
# save(sim_rep_df, file = glue::glue("results/{simdir}/{simdir}_rep_df_1000reps.rda"))


