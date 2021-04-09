# sessionInfo()

# Load libraries
require(tidyverse)
require(here)

# Declare path 
here::i_am("invasive_spread.Rproj")
# here::here()

# Helper functions
'%ni%' <- function(x,y)!('%in%'(x,y))

numextract <- function(string){ 
  temp <- gsub("[^[:digit:].]", "", string)
  as.numeric(str_extract(temp, "\\-*\\d+\\.*\\d*"))
} 

webs <- (here("data/webs/bodymass_*.csv") %>% Sys.glob) %>% basename %>% 
  str_split(.,"_",simplify=TRUE) %>% as_tibble %>% pull(V2) %>% unique

landscapes <- (here("data/landscapes/landscape_*.csv") %>% Sys.glob) %>% basename %>%
  str_split(.,"_",simplify=TRUE) %>% as_tibble %>% pull(V2) %>% unique %>% numextract %>% sort

nutsupplys <- c(0.1, 1, 1000) 

# CONSTATNS
JOINT_BASH <- 0 # if 1: joint scenario (NOT RELEVANT HERE!)
ISOLATED_BASH <- 0 # if 1: isolated scenario (NOT RELEVANT HERE!)
TIMESERIES_BASH <- 0 # if 1: time series is written to file (instead of output file)

TEND_BASH <- 5000 # number of time steps for simulation run
TEVAL_BASH <- 1000 # time steps for evaluation 

DIR <- here("data/")
OUTPUT_DIR <- here("data/output/")
INVASION_DIR <- here("data/invasion/")
TS_DIR <- here("data/timeseries/")
TS_FILE <- here("data/timeseries/timeseries_")

tab <- expand_grid(WEB=webs, LANDSCAPE=landscapes, NUTSUPPLY=nutsupplys) %>%
    as_tibble %>% # make into tibble
    arrange(desc(WEB)) %>%  # sort by web 
    expand_grid(REMSPP = 1:21) %>% 
    add_column(SEED = .$WEB %>% as.integer()) %>%
    select(SEED, WEB, LANDSCAPE, REMSPP, NUTSUPPLY) %>% # change column order
    add_column(DIR, OUTPUT_DIR, INVASION_DIR, TS_DIR, TS_FILE,  
               JOINT_BASH, ISOLATED_BASH, TIMESERIES_BASH, TEND_BASH, TEVAL_BASH)

for (r in 1:nrow(tab)) {
    jobname <- sprintf(paste0(tab$WEB[r],"_%03d.sh"),r)
    f <- file(jobname)
    writeLines(c(
      paste("#!/bin/bash"),
      paste("#"),
      paste("# Input variables"),
      paste("# $SEED $WEB $LANDSCAPE $REMSPP $NUTSUPPLY $INVASION_BASH $DIR $OUTPUT_DIR $INVASION_DIR $TS_DIR $TS_FILE $JOINT_BASH $ISOLATED_BASH $TIMESERIES_BASH $TEND_BASH $TEVAL_BASH"), 
      paste("# Run simiulation"),
      paste("bash invasion_process.sh", paste(tab[r,], collapse=" ")), 
      paste("#"),
      paste("# Script ends here")
    ), f)
    close(f)
}
