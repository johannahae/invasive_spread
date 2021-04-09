
rm(list=ls())

# Load libraries
require(NetIndices)
require(graph)
require(tidyverse)
require(here)

# Declare path 
here::i_am("invasive_spread.Rproj")
# here::here()

input_dir <- here("data/output")

# Helper functions
'%ni%' <- function(x,y)!('%in%'(x,y))

isEmpty <- function(x)
  return(length(x)==0)

# Functcion to add the invasion start patch to each global file
add_invstartpatch <- function(inputfile){
  outdir <- dirname(inputfile)
  invasionstartpatch <- inputfile %>% str_replace("global", "mass") %>%
    read_csv(col_types = cols()) %>% 
    filter(if.inv.spp == 1, init.biomass > 0) %>% 
    pull(patch)
  inputfile %>% read_csv(col_types = cols()) %>%
    mutate(inv.start.patch = invasionstartpatch) %>%
    write_csv(inputfile)
}

# Summarize output files for analysis
global_files <- list.files(path=input_dir, pattern = "^global(.*)_1.csv$", all.files = TRUE, 
                           full.names = TRUE, no.. = TRUE, include.dirs = TRUE, recursive=TRUE) 

for(global_file in global_files) add_invstartpatch(global_files)

mass_files <- list.files(path=input_dir, pattern = "^mass(.*)_1.csv$", all.files = TRUE, 
                         full.names = TRUE, no.. = TRUE, include.dirs = TRUE, recursive=TRUE) 

EXTINCT <- 10^-20 # extinction threshold   

r <- 1
for(mass_file in mass_files){
     
      mass <- mass_file %>% read_csv(col_types = cols()) # read in MASS file 
    
      mass0 <- mass_file %>% str_replace("_1.csv","_0.csv") %>% 
        read_csv(col_types = cols()) # read in corresponding MASS file PRE-INVASION
      
      global <- mass_file %>% str_replace("mass", "global") %>% 
        read_csv(col_types = cols()) # read in corresponding GLOBAL file 
  
      params <- mass_file %>% str_replace("mass", "params") %>% 
        read_csv(col_types = cols()) # read in corresponding PARAMS file 
      
      Amat <- here("data/webs", paste0("web_", global$web, "_99.csv")) %>%
        read.table %>% # read in full adjacency matrix 
        as.matrix # transform into regular matrix 
      
      species_traits <- here("data/webs", paste0("bodymass_", global$web,"_99.csv")) %>%
        read_csv(col_types = cols()) # read in species traits (bodymass, if.basal.spp, D_Max)
      
      A <- mass_file %>% str_replace("mass", "web") %>%
        read.table %>% as.matrix # read in Amat
      A01 <- ifelse(A[,] > 0, 1, 0) # make binary
      
      A0 <- mass_file %>% str_replace("mass", "web") %>% str_replace("_1.csv", "_0.csv") %>%
        read.table %>% as.matrix # read in Amat PRE-INVASION
      A001 <- ifelse(A0[,] > 0, 1, 0) # make binary
    
      gamma.spp <- mass0 %>% filter(init.biomass > EXTINCT, species<1000) %>%
        pull(species) %>% unique %>% sort #regional species pool
      spp <- 0:(global$number.of.spp-1)
      
      gamma <- gamma.spp %>% length #regional species richness 
      mean.alpha <- mass0 %>% filter(species<1000, Biomassess_tend>EXTINCT) %>% group_by(patch) %>%
        summarise(alpha = species %>% length) %>% pull(alpha) %>% mean(., na.rm=TRUE) #local species richness (mean over patches)
      beta <- gamma/mean.alpha #differences in community composition)
  
      gamma.init.spp <- mass %>% filter(init.biomass > EXTINCT, species<1000) %>%
        pull(species) %>% unique %>% sort #regional species pool 
      gamma.init <- gamma.init.spp %>% length #regional species richness
      mean.alpha.init <- mass %>% filter(species<1000, init.biomass>EXTINCT) %>% group_by(patch) %>%
        summarise(alpha = species %>% length) %>% pull(alpha) %>% mean(., na.rm=TRUE) #local species richness (mean over patches)
      beta.init <- gamma.init/mean.alpha.init #beta-diversity (community composition)
      
      gamma.post.spp <- mass %>% filter(Biomassess_tend > EXTINCT, species<1000) %>%
        pull(species) %>% unique %>% sort #regional species pool 
      gamma.post <- gamma.post.spp %>% length #regional speceis richness
      mean.alpha.post <- mass %>% filter(species<1000, Biomassess_tend>EXTINCT) %>% group_by(patch) %>%
        summarise(alpha = species %>% length) %>% pull(alpha) %>% mean(., na.rm=TRUE) #local species richness (mean over patches)
      beta.post <- gamma.post/mean.alpha.post #differences in community composition)
      
      sumb.native.t0 <- mass0 %>% filter(init.biomass > EXTINCT, species<1000, if.inv.spp==0) %>%
        summarise(sumb = sum(init.biomass, na.rm=TRUE)) %>% ungroup %>% pull(sumb)
      sumb.invspp.t0 <- mass0 %>% filter(init.biomass > EXTINCT, species<1000, if.inv.spp==1) %>%
        summarise(sumb = sum(init.biomass, na.rm=TRUE)) %>% ungroup %>% pull(sumb)
      if(sumb.invspp.t0>0) print("Something went wrong...")
      
      cor.coeff <- cor(mass0 %>% pull(init.biomass), mass0 %>% pull(Biomassess_tend), method="pearson")
      cor.coeff.spp <- cor(mass0 %>% filter(species<100) %>% pull(init.biomass),
                           mass0 %>% filter(species<100) %>% pull(Biomassess_tend), method="pearson")
      cor.coeff.nut <- cor(mass0 %>% filter(species>100) %>% pull(init.biomass),
                           mass0 %>% filter(species>100) %>% pull(Biomassess_tend), method="pearson")
      
      meanb.native.t0 <- mass0 %>% filter(init.biomass > EXTINCT, species<1000, if.inv.spp==0) %>%
        summarise(sumb = mean(init.biomass, na.rm=TRUE)) %>% ungroup %>% pull(sumb)
      
      sumb.native.init <- mass %>% filter(init.biomass > EXTINCT, species<1000, if.inv.spp==0) %>%
        summarise(sumb = sum(init.biomass, na.rm=TRUE)) %>% ungroup %>% pull(sumb)
      sumb.invspp.init <- mass %>% filter(init.biomass > EXTINCT, species<1000, if.inv.spp==1) %>%
        summarise(sumb = sum(init.biomass, na.rm=TRUE)) %>% ungroup %>% pull(sumb)
      
      meanb.native.init <- mass %>% filter(init.biomass > EXTINCT, species<1000, if.inv.spp==0) %>%
        summarise(sumb = mean(init.biomass, na.rm=TRUE)) %>% ungroup %>% pull(sumb)
      meanb.invspp.init <- mass %>% filter(init.biomass > EXTINCT, species<1000, if.inv.spp==1) %>%
        summarise(sumb = mean(init.biomass, na.rm=TRUE)) %>% ungroup %>% pull(sumb)
      
      sumb.native.post <- mass %>% filter(Biomassess_tend > EXTINCT, species<1000, if.inv.spp==0) %>%
        summarise(sumb = sum(Biomassess_tend, na.rm=TRUE)) %>% ungroup %>% pull(sumb)
      sumb.invspp.post <- mass %>% filter(Biomassess_tend > EXTINCT, species<1000, if.inv.spp==1) %>%
        summarise(sumb = sum(Biomassess_tend, na.rm=TRUE)) %>% ungroup %>% pull(sumb)
  
      meanb.native.post <- mass %>% filter(Biomassess_tend > EXTINCT, species<1000, if.inv.spp==0) %>%
        summarise(sumb = mean(Biomassess_tend, na.rm=TRUE)) %>% ungroup %>% pull(sumb)
      meanb.invspp.post <- mass %>% filter(Biomassess_tend > EXTINCT, species<1000, if.inv.spp==1) %>%
        summarise(sumb = mean(Biomassess_tend, na.rm=TRUE)) %>% ungroup %>% pull(sumb)
      
      invaded.patches <- mass %>% filter(if.inv.spp == 1, Biomassess_tend >= EXTINCT, init.biomass==0) %>%
        pull(patch) %>% unique #"invaded" patches (without the invasion start patch)
    
      invspp.established.patch <- mass %>% 
        filter(if.inv.spp == 1, Biomassess_tend >= EXTINCT, init.biomass>0) %>% pull(patch) 
      invspp.established <- ifelse(length(invaded.patches)>0 &&
                                     !is_empty(invspp.established.patch), 1, 0) #invasive species established in invasion start patch 

      coords <-  params %>% filter(variable == "patch.location.XY.alternating") %>% #extract coords
        pull(value) #vector with XY alternating

      x <- coords[seq(1, length(coords), by = 2)] #extract X coordinates
      y <- coords[seq(2, length(coords), by = 2)] #extract Y coordinates
      landscape <- cbind(x,y) #generate landscape matrix
      distmat <- dist(landscape) %>% as.matrix #distance matrix containing the euclidean distances between all patches
      distvec <- distmat[global$inv.start.patch, ] #extract distances from invasion start patch to all other patches
      maxdist <- ifelse(invaded.patches %>% length>0, #invaded patch with maximum distance 
                        distvec[invaded.patches+1] %>% max, 0) #to the invasion start patch
    
      no.consumers.surv.init <- mass0 %>% filter(species<100, if.basal.spp==0) %>% group_by(species) %>% #globally extinct species
        summarize(surv=sum(ifelse(Biomassess_tend >= EXTINCT, 1, 0))) %>% 
        filter(surv>0) %>% pull(species) %>% length  
      no.basals.surv.init <- mass0 %>% filter(species<100, if.basal.spp==1) %>% group_by(species) %>% #globally extinct species
        summarize(surv=sum(ifelse(Biomassess_tend >= EXTINCT, 1, 0))) %>% 
        filter(surv>0) %>% pull(species) %>% length 
      
      no.consumers.surv.post <- mass %>% filter(species<100, if.basal.spp==0) %>% group_by(species) %>% #globally extinct species
        summarize(surv=sum(ifelse(Biomassess_tend >= EXTINCT, 1, 0))) %>% 
        filter(surv>0) %>% pull(species) %>% length  
      no.basals.surv.post <- mass %>% filter(species<100, if.basal.spp==1) %>% group_by(species) %>% #globally extinct species
        summarize(surv=sum(ifelse(Biomassess_tend >= EXTINCT, 1, 0))) %>% 
        filter(surv>0) %>% pull(species) %>% length 
    
      occ.patches <- mass %>% filter(species<1000) %>% group_by(patch) %>% #occupied patches
        summarize(surv=sum(ifelse(Biomassess_tend >= EXTINCT, 1, 0))) %>% 
        filter(surv>0) %>% pull(patch)
      
      spp <- 0:(nrow(Amat)-1)
      colnames(A01) <- rownames(A01) <- spp
      
      extinct.spp.init <- spp[-which(spp %in% gamma.init.spp)]
      if(is_empty(extinct.spp.init)) extinct.spp.init <- NA
      
      extinct.spp.post <- spp[-which(spp %in% gamma.post.spp)]
      if(is_empty(extinct.spp.post)) extinct.spp.post <- NA
      
      spp0 <- spp[-which(spp==global$rem.spp)]
      colnames(A001) <- rownames(A001) <- spp0
      
      #Native web properties initialized
      TL.native <- NetIndices::TrophInd(t(A001))$TL #vector with trophic levels (length S)
      OI.native <- NetIndices::TrophInd(t(A001))$OI #vector with omnivory indices (length S)
      G.native <- A001 %>% rowSums %>% as.integer #vector with generalities (length S)
      V.native <- A001 %>% colSums %>% as.integer #vector with vulnerabilities (length S)
      C.native <- sum(A001)/(nrow(A001)^2) #connectance C = L/S²
      L.native <- sum(A001)/nrow(A001) #link density L/S
      no.of.links.native <- G.native+V.native  #total number of links per speices (length S)

      #Native web properties post-simulation
      extinct.spp.native <- spp[-which(spp0 %in% gamma.init.spp)]
      if(is_empty(extinct.spp.native)) extinct.spp.native <- NA
      missing.spp.native <- which((colnames(A001) %in% extinct.spp.native)==TRUE)
      if(is_empty(missing.spp.native)) A.native <- A001
      if(!is_empty(missing.spp.native)) A.native <- A001[-missing.spp.native, -missing.spp.native]
      
      TL.native.post <- OI.native.post <- G.native.post <- V.native.post <- C.native.post <- L.native.post <- no.of.links.native.post <- NA
      if(is.matrix(A.native) && nrow(A.native>0)){
        TL.native.post <- NetIndices::TrophInd(t(A.native))$TL #vector with trophic levels (length S)
        OI.native.post <- NetIndices::TrophInd(t(A.native))$OI #vector with omnivory indices (length S)
        G.native.post <- A.native %>% rowSums %>% as.integer #vector with generalities (length S)
        V.native.post <- A.native %>% colSums %>% as.integer #vector with vulnerabilities (length S)
        C.native.post <- sum(A.native)/(nrow(A.native)^2) #connectance C = L/S²
        L.native.post <- sum(A.native)/nrow(A.native) #link density L/S
        no.of.links.native.post <- G.native.post+V.native.post  #total number of links per speices (length S)
      }
      
      #Potential web properties 
      TL <- NetIndices::TrophInd(t(A01))$TL #vector with trophic levels (length S)
      OI <- NetIndices::TrophInd(t(A01))$OI #vector with omnivory indices (length S)
      G <- A01 %>% rowSums %>% as.integer #vector with generalities (length S)
      V <- A01 %>% colSums %>% as.integer #vector with vulnerabilities (length S)
      C <- sum(A01)/(nrow(A01)^2) #connectance C = L/S²
      L <- sum(A01)/nrow(A01) #link density L/S
      no.of.links <- G+V  #total number of links per speices (length S)
      remspp.id <- global$rem.spp + 1
      
      #Web properties at time of introduction 
      missing.spp.init <- which((colnames(A01) %in% extinct.spp.init)==TRUE)
      if(is_empty(missing.spp.init)) A.init <- A01
      if(!is_empty(missing.spp.init)) A.init <- A01[-missing.spp.init, -missing.spp.init]
      
      TL.init <- OI.init <- G.init <- V.init <- C.init <- L.init <- no.of.links.init <- NA
      if(is.matrix(A.init) && nrow(A.init>0)){
        TL.init <- TrophInd(t(A.init))$TL #vector with trophic levels (length S)
        OI.init <- TrophInd(t(A.init))$OI #vector with omnivory indices (length S)
        G.init <- A.init %>% rowSums %>% as.integer #vector with generalities (length S)
        V.init <- A.init %>% colSums %>% as.integer #vector with vulnerabilities (length S)
        C.init <- sum(A.init)/(nrow(A.init)^2) #connectance C = L/S²
        L.init <- sum(A.init)/nrow(A.init) #link density L/S
        no.of.links.init <- G.init + V.init
      } 

      remspp.id.init <- which(colnames(A.init)==global$rem.spp)
      if(global$rem.spp %ni% colnames(A.init)) remspp.id.init <- 1

      #Web properties post invasion
      missing.spp.post <- which((colnames(A01) %in% extinct.spp.post)==TRUE)
      if(is_empty(missing.spp.post)) A.post <- A01
      if(!is_empty(missing.spp.post)) A.post <- A01[-missing.spp.post, -missing.spp.post]
      
      TL.post <- OI.post <- G.post <- V.post <- C.post <- L.post <- no.of.links.post <- NA
      is.matrix(A.post)
      if(is.matrix(A.post) && nrow(A.post>0)){
        TL.post <- TrophInd(t(A.post))$TL #vector with trophic levels (length S)
        OI.post <- TrophInd(t(A.post))$OI #vector with omnivory indices (length S)
        G.post <-  A.post %>% rowSums %>% as.integer #vector with generalities (length S)
        V.post <-  A.post %>% colSums %>% as.integer #vector with vulnerabilities (length S)
        C.post <- sum(A.post)/(nrow(A.post)^2) #connectance C = L/S²
        L.post <- sum(A.post)/nrow(A.post) #link density L/S
        no.of.links.post <- G.post + V.post
      } 
      
      remspp.id.post <- which(colnames(A.post)==global$rem.spp)
      if(global$rem.spp %ni% colnames(A.post)) remspp.id.post <- 1
      
      tmp <- global$rem.spp*(global$number.of.patch*global$number.of.patch) 
      tmp2 <- (global$rem.spp+1)*(global$number.of.patch*global$number.of.patch)
      tmp3 <- ifelse(params %>% filter(variable=="SW.dispersal.matrix") %>% 
                       slice((tmp+1):tmp2) %>% pull(value) > 0, 1, 0) %>% sum 
      rgg.con.invspp <- tmp3/(global$number.of.patch^2-global$number.of.patch) #no.realized.links/no.possible.links (Z^2-Z)

      bodymass.invspp <- params %>% filter(variable=="bodymasses") %>% slice(global$rem.spp+1) %>% pull(value)
      D_Max <- params %>% filter(variable=="D_Max") %>% pull(value)
      D_Max.invspp <- params %>% filter(variable=="D_Max") %>% slice(global$rem.spp+1) %>% pull(value)
      
      tmp4 <- global %>% mutate(landscape.type = ifelse(.$landscape > 10, "clustered", "random"),
                               invspp.type = ifelse(.$rem.spp < .$number.of.plants, "basal species", "consumer species"),
                               tot.biomass.invspp.post = sumb.invspp.post, 
                               tot.biomass.native.post = sumb.native.post,
                               tot.biomass.invspp.int = sumb.invspp.init, 
                               tot.biomass.native.int = sumb.native.init,
                               tot.biomass.native.t0 = sumb.native.t0,
                               mean.biomass.invspp.post = meanb.invspp.post, 
                               mean.biomass.native.post = meanb.native.post,
                               mean.biomass.invspp.int = meanb.invspp.init, 
                               mean.biomass.native.int = meanb.native.init,
                               mean.biomass.native.t0 = meanb.native.t0,
                               gamma = gamma, 
                               mean.alpha = mean.alpha, 
                               beta = beta, 
                               gamma.int = gamma.init, 
                               mean.alpha.int = mean.alpha.init, 
                               beta.int = beta.init, 
                               gamma.post = gamma.post, 
                               mean.alpha.post = mean.alpha.post, 
                               beta.post = beta.post, 
                               number.of.invaded.patches = invaded.patches %>% length, 
                               fraction.of.invaded.patches = (invaded.patches %>% length) / (number.of.patch-1),
                               invasion.distance = maxdist,
                               invspp.established = invspp.established, 
                               bodymass.invspp = bodymass.invspp,
                               dispersal.dist.invspp = D_Max.invspp, 
                               rgg.con.invspp = rgg.con.invspp, 
                               trophic.level.invspp = TL[global$rem.spp+1],
                               trophic.level.invspp.int = TL.init[remspp.id.init],
                               trophic.level.invspp.post = TL.post[remspp.id.post],
                               omnivory.index.invspp = OI[global$rem.spp+1],
                               omnivory.index.invspp.int = OI.init[remspp.id.init],
                               omnivory.index.invspp.post = OI.post[remspp.id.post],
                               generality.invspp = G[global$rem.spp+1],
                               generality.invspp.int = G.init[remspp.id.init],
                               generality.invspp.post = G.post[remspp.id.post],
                               vulnerability.invspp = V[global$rem.spp+1], 
                               vulnerability.invspp.int = V.init[remspp.id.init], 
                               vulnerability.invspp.post = V.post[remspp.id.post], 
                               no.of.links.invspp = no.of.links[global$rem.spp+1], 
                               no.of.links.invspp.int = no.of.links.init[remspp.id.init], 
                               no.of.links.invspp.post = no.of.links.post[remspp.id.post], 
                               T_mean = mean(TL), #mean trophic level
                               T_max = max(TL),  #maximum trophic level
                               G_mean = mean(G),  #mean generality
                               G_sd = sd(G), #sd generality
                               V_mean = mean(V),  #mean vulnerability
                               V_sd = sd(V), #sd vulnerability
                               C = C, #connectance 
                               L = L, #link density
                               T_mean.native = mean(TL.native), #mean trophic level native web
                               T_max.native = max(TL.native),  #maximum trophic level native web
                               G_mean.native = mean(G.native),  #mean generality native web
                               G_sd.native = sd(G.native), #sd generality native web
                               V_mean.native = mean(V.native),  #mean vulnerability native web
                               V_sd.native = sd(V.native),  #sd vulnerability native web
                               C.native = C.native, #connectance native web
                               L.native = L.native, #link density native web
                               T_mean.native.post = mean(TL.native.post), #mean trophic level native web
                               T_max.native.post = max(TL.native.post),  #maximum trophic level native web
                               G_mean.native.post = mean(G.native.post),  #mean generality native web
                               G_sd.native.post = sd(G.native.post), #sd generality native web
                               V_mean.native.post = mean(V.native.post),  #mean vulnerability native web
                               V_sd.native.post = sd(V.native.post),  #sd vulnerability native web
                               C.native.post = C.native.post, #connectance native web
                               L.native.post = L.native.post, #link density native web
                               T_mean.int = mean(TL.init), #mean trophic level at time of introduction
                               T_max.int = max(TL.init),  #maximum trophic level at time of introduction
                               G_mean.int = mean(G.init),  #mean generality at time of introduction
                               G_sd.int = sd(G.init), #sd generality at time of introduction
                               V_mean.int = mean(V.init),  #mean vulnerability at time of introduction
                               V_sd.int = sd(V.init),  #sd vulnerability at time of introduction
                               C.int = C.init, #connectance at time of introduction
                               L.int = L.init, #link density at time of introduction
                               T_mean.post = mean(TL.post), #mean trophic level at time of introduction
                               T_max.post = max(TL.post),  #maximum trophic level at time of introduction
                               G_mean.post = mean(G.post),  #mean generality at time of introduction
                               G_sd.post = sd(G.post), #sd generality at time of introduction
                               V_mean.post = mean(V.post),  #mean vulnerability at time of introduction
                               V_sd.post = sd(V.post),  #sd vulnerability at time of introduction
                               C.post = C.post, #connectance at time of introduction
                               L.post = L.post, #link density at time of introduction
                               no.occupied.patches = occ.patches %>% length,
                               extant.consumers.int = no.consumers.surv.init, 
                               extant.consumers.post = no.consumers.surv.post, 
                               extant.plants.int = no.basals.surv.init, 
                               extant.plants.post = no.basals.surv.post, 
                               cor.coeff = cor.coeff, 
                               cor.coeff.spp.biomass = cor.coeff.spp, 
                               cor.coeff.nut.concentration = cor.coeff.nut)
                                 
      if(r==1) dat <- tmp4 #intialize tibble once
      if(r>1)  dat <- rbind(dat,tmp4) #append table
      r <- r+1
      
      dat <- dat %>% mutate(trophic.level.invspp = ifelse(
            (invspp.type == "consumer species" & trophic.level.invspp ==  1), NA, 
            trophic.level.invspp),
        trophic.level.invspp.int = ifelse(
          (invspp.type == "consumer species" & trophic.level.invspp.int ==  1), NA, 
          trophic.level.invspp.int), 
      trophic.level.invspp.post = ifelse(
        (invspp.type == "consumer species" & trophic.level.invspp.post ==  1), NA, 
        trophic.level.invspp.post))
}

dat %>% write_csv(here("data/results/summary.csv"))

