
# Load libraries
require(tidyverse)
require(here)

# Declare path 
here::i_am("invasive_spread.Rproj")
# here::here()

infiles <- here("data/webs/bodymass_*_99.csv") %>% Sys.glob

for(infile in infiles){
  
  bodymass <- infile %>% read_csv 
  S <- bodymass %>% nrow
  
  for(i in 1:S)
    bodymass %>% slice(-i) %>% 
      write_csv(infile %>% str_replace("_99.csv", paste0("_", (i-1), ".csv")))
}



