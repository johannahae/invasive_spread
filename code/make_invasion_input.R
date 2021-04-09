
rm(list=ls())

# Load libraries
library(tidyverse)
library(here)

# Declare path 
here::i_am("invasive_spread.Rproj")
# here::here()

globals <- list.files(
  path = here("data/output"), 
  pattern = "^[global].*._0.csv$", 
  recursive = TRUE,
  full.names = TRUE
  ) 

for(q in globals){
  global <- q %>% read_csv(col_types = cols())
  name <- basename(q) %>% str_split(.,"global_", simplify=TRUE) %>% max ## extract file name
  mass <- paste0(dirname(q), "/mass_", name) %>% read_csv(col_types = cols()) ## read in results/mass.csv from simulation with S species 

  bodymass <- here("data/webs", paste0("bodymass_", global$web, "_99.csv")) %>% # invaded web
    read_csv(col_types = cols()) %>% #read in webs/bodymass.csv for S+1 species 
    mutate(species = 0:(nrow(.)-1))
  
  Sci <- bodymass %>% filter(if.basal == 0) %>% nrow # number of consumer species 
  Sbi <- bodymass %>% filter(if.basal == 1) %>% nrow # number of basal species 
  Si <- Sbi + Sci

  # testing...
  if(Si - global$number.of.spp != 1) "There is an issue with the number of species..." 
  
  remspp.is.basal <- ifelse((Sbi - global$number.of.plants)==1, 1, 0) # invader basal or consumer (1/0)?
  remspp.bodymass <- bodymass %>% filter(species == global$rem.spp) %>% pull(bodymass) ## invader bodymass
  
  if.invspp <- rep(0, nrow(bodymass)) # empty vector 
  for(s in 1:nrow(bodymass))
    if.invspp[s] <- ifelse(bodymass$species[s]==global$rem.spp,1,0) ## invasive species 
  
  #select 1st patch in the upper left corner in which the invasive species will be initialized
  landscape <- here("data/landscapes", paste0("landscape_", global$landscape,".csv")) %>% 
    read_csv(col_types = cols()) %>% 
    add_column(patch.number = 0:(nrow(.)-1))

  Z <- nrow(landscape) # number of patches
  if(Z != global$number.of.patch) "There is an issue with the number of patchess..."
  distdiff <- rep(0, Z)
  
  for(i in 1:Z){
    minx <- landscape %>% slice(i) %>% pull(x) %>% min
    maxy <- landscape %>% slice(i) %>% pull(y) %>% max
    distdiff[i] <- (maxy-minx)
  }
  
  pi <- landscape$patch.number[which(distdiff == max(distdiff))] # invasion start patch 
 
  # generate a tibble with the bodymasses etc. of S+1 species and save it in invasion/inputBodymss_*csv
  inputBodymass <- bodymass %>% add_column(if.invspp = if.invspp) %>% # add invasive species identifier
    select(bodymass, if.basal, if.invspp, species, dispersal.dist) %>% # change order of columns
    write_csv(paste0((dirname(q) %>% str_replace("output", "invasion")), "/inputBodymass_", name)) # save to invasion/
  
  mveci <- rep(0, length(1:Z)) # initizialize vector of length Z
  m <- 5.0 # mean initialized biomass density over all species 

  mveci[pi+1] <- m # initialize invader biomass on invasion start patch 
  # must be pi+1 because in C++ we start with 0 (and not 1)

  # generate a tibble with biomassess to initialize S+1 species and N nutrients on Z patches 
  # order in vector has to bee the same as in C++ B[] => D+DN => (S+1)*Z + N*Z !!! 
  inputB <- mass %>% select(Biomassess_tend, bodymass, patch, species) %>%
    add_row(patch=0:(Z-1), bodymass = remspp.bodymass, Biomassess_tend = mveci, species = global$rem.spp,.before=1) %>% ## add invasive species
    filter(bodymass>0) %>% arrange(desc(-species)) %>% arrange(desc(-patch)) %>% #sort vector 
    rbind(., mass %>% filter(bodymass==0) %>% select(Biomassess_tend, bodymass, patch, species) %>% 
            arrange(desc(-patch))) %>% # move nutrients to the tail 
    write_csv(paste0((dirname(q) %>% str_replace("output", "invasion")), "/inputBiomass_", name)) # save to invasion/
}

