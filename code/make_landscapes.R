
# Load libraries
require(igraph)
require(tidyverse)
require(here)

# Declare path 
here::i_am("invasive_spread.Rproj")
# here::here()

# Helper functions
numextract <- function(string){ 
  temp <- gsub("[^[:digit:].]", "", string)
  as.numeric(str_extract(temp, "\\-*\\d+\\.*\\d*"))
} 

# Truncated normal sampler
# Input: n = number of observations; mean = vector of means; sd = vector of standard deviations; 
# a = lower boundary; b: upper boundary; 
# Output: - vector of n random values within the set constraints.
rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){  
  qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}

# Randomly draw patch coordinates. 
# Input: N: number of patches; Output: matrix with N rows and 2 columns. 
generate_coordinates <- function(N) {
    return(matrix(runif(N * 2, 0, 1), nrow=N, ncol=2))
}


# Make landscapes 

N <- 40 # number of patches
seed <- 704 

## 1: Generate random landscapes
for(q in 1:10){ # landscapes 1-10
  set.seed(seed + num)
  landscape <- generate_coordinates(N) %>% as_tibble %>% rename(x=V1, y=V2) %>% 
    write_csv(sprintf(here("data/landscapes/landscape_%d.csv"), q))
  q <- q+1
}

## 2: Generate clustered landscapes => small world scenario with spatial autocorrelation 
no.of.clusters <- 8 # number of clusters

for(q in 11:20){   # landscapes 11-20
  
  set.seed(seed + q)
  p <- 1
  x <- rep(0, N)
  y <- rep(0, N)
  
  posx <- runif(no.of.clusters)
  posy <- runif(no.of.clusters)
  mat <- cbind(posx, posy)
  
  while(max(dist(mat))> 0.294 & min(dist(mat))<0.2){ # set max and min distances for clusters
    posx <- runif(no.of.clusters)
    posy <- runif(no.of.clusters)
    mat <- cbind(posx,posy)
  }
  
  p <- 1
  for(n in 1:no.of.clusters){
    # truncated normale sampler; the smaller sd, the higher the clustering
    x[p:(n*(N/no.of.clusters))] <- rtnorm(N/no.of.clusters, mean=posx[n], sd=0.03, a=0, b=1) 
    y[p:(n*(N/no.of.clusters))] <- rtnorm(N/no.of.clusters, mean=posy[n], sd=0.03, a=0, b=1)  
    p <- n*(N/no.of.clusters)+1 
  }
  
  landscape <- tibble(x=x, y=y) %>%  
    write_csv(sprintf(here("data/landscapes/landscape_%d.csv"), q))
  
  q <- q + 1
}



