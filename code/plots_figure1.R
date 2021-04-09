

# sessionInfo() # R version 3.6.3

#### Load libraries ####
require(tidyverse)
require(ggplot2)
require(igraph)
require(cowplot)
require(here)

#### Declar path ####
here::i_am("invasive_spread.Rproj")
# here::here()

#### Helper functions ####
'%ni%' <- function(x, y)
  ! ('%in%'(x, y))

quantiles_95 <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) #95%CI
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

numextract <- function(string) {
  temp <- gsub("[^[:digit:].]", "", string)
  as.numeric(str_extract(temp, "\\-*\\d+\\.*\\d*"))
}

#### Figure 1A: Expecatations landscape properties ####
# pdf(here("figures/figure_1A.pdf"), width=18,height=10)
par(pty = "s")

plot(
  1:10,
  type = "n",
  col = "blue",
  lwd = 2,
  xlim = c(0, 10),
  ylim = c(0, 10),
  xlab = " ",
  ylab = "Invasion success\n",
  xaxt = "n",
  yaxt = "n",
  cex.lab = 4,
  font.lab = 1,
  axes = F
)

points(
  x = c(0.5, 7.5),
  y = c(2, 5),
  col = c("blue", "darkgreen"),
  lwd = 2,
  pch = 15,
  cex = 10
)
points(
  x = c(0.5, 7.5),
  y = c(7.5, 8.5),
  col = c("blue", "darkgreen"),
  lwd = 2,
  pch = 16,
  cex = 10
)

arrows(
  -1,
  0,-1,
  10,
  col = "black",
  lwd = 3,
  pch = 15,
  cex = 10,
  xpd = TRUE
)
arrows(
  1.4,
  7.6,
  6.6,
  8.6,
  col = "black",
  lwd = 2,
  pch = 15,
  cex = 10
)
arrows(
  1.4,
  2.1,
  6.6,
  5,
  col = "black",
  lwd = 2,
  pch = 15,
  cex = 10
)
arrows(
  0.5,
  2.9,
  0.5,
  6.67,
  col = "black",
  lwd = 2,
  pch = 15,
  cex = 10
)
arrows(
  7.5,
  5.9,
  7.5,
  7.58,
  col = "black",
  lwd = 2,
  pch = 15,
  cex = 10
)

mtext(
  side = 1,
  text = "random \n",
  at = 7.5,
  col = "darkgreen",
  line = 1,
  cex = 3
)
mtext(
  side = 1,
  text = "clustered \n",
  at = 0.9,
  col = "blue",
  line = 1,
  cex = 3
)

mtext(
  side = 4,
  text = "oligotrophic",
  at = 3.5,
  col = "black",
  line = 1,
  cex = 3,
  las = 1
)
mtext(
  side = 4,
  text = "eutrophic",
  at = 8.5,
  col = "black",
  line = 1,
  cex = 3,
  las = 1
)

points(
  x = c(9.85, 9.85),
  y = c(3.5, 8.5),
  lwd = 2,
  pch = c(22, 21),
  cex = 10
)
# dev.off()

#### Figure 1B: Expecatations invasive spread ####
# pdf(here("figures/figure_1B.pdf"), width=18,height=10)

# Plot settings
par(pty = "s")
par(mar = c(0, 0, 0, 0) + 0.1) # c(bottom, left, top, right). The default is c(5, 4, 4, 2) + 0.1.
par(oma = c(0, 0, 0, 0)) # all sides have 3 lines of space

EXTINCT <- 10 ^ -20 # extinction threshold

# Read in mass files
massfiles <- c(
  "../data/output/mass_8516_20_12_1_1.csv",
  "../data/output/mass_8516_2_12_1_1.csv",
  "../data/output/mass_8516_20_8_1_1.csv",
  "../data/output/mass_8516_2_8_1_1.csv"
)

fig <- 1
for (massfile in massfiles) {
  png(here("figures", paste0("figure1B_", fig, ".png")), bg = "transparent")
  fig <- fig + 1
  
  mass <- massfile %>% read_csv()
  global <-
    massfile %>% str_replace("mass", "global") %>% read_csv
  
  dispmat <- massfile %>% str_replace("mass", "params") %>%
    read_csv %>%
    filter(variable == "SW.dispersal.matrix") %>%
    pull(value)
  
  Z <- global$number.of.patch
  x <- here("data/landscapes/",
            paste0("landscape_", global$landscape, ".csv")) %>%
    read_csv() %>%  pull(x)
  y <- here("data/landscapes/",
            paste0("landscape_", global$landscape, ".csv")) %>%
    read_csv() %>%  pull(y)
  
  patchcol <-
    ifelse(str_detect(massfile, "_20_"), "blue", "darkgreen")
  
  colvec <- rep("white", Z) # initialize color vector
  invaded <-
    mass %>% filter(if.inv.spp == 1, Biomassess_tend > EXTINCT) %>% pull(patch)
  invaded <- invaded + 1
  notinvaded <- (1:Z)[-invaded]
  colvec[invaded] <- patchcol # change color of invaded patches
  colvec[global$inv.start.patch + 1] <-
    "black" # change color of invasion start patch
  
  invspp <- global$rem.spp
  
  plot(
    1,
    type = "n",
    ylim = c(0, 1),
    xlim = c(0, 1),
    xlab = "",
    ylab = "",
    xaxt = "n",
    yaxt = "n"
  )
  dat <-
    matrix(dispmat[(1 + invspp * Z * Z):(Z * Z + invspp * Z * Z)], ncol =
             Z, nrow = Z)
  for (i in 1:Z) {
    for (j in 1:Z) {
      if (dat[i, j] != 0) {
        x1 = x[i]
        y1 = y[i]
        x2 = x[j]
        y2 = y[j]
        segments(x1,
                 y1,
                 x2,
                 y2,
                 col = "lightgrey",
                 lty = 3,
                 lwd = 1.5)
      }
    }
  }
  
  points(x[invaded],
         y[invaded],
         pch = 21,
         cex = 5,
         bg = colvec[invaded])
  
  points(x[notinvaded],
         y[notinvaded],
         pch = 21,
         cex = 5,
         bg = colvec[notinvaded])
  dev.off()
}

plots <- here("figures/figure1B_*.png") %>% Sys.glob()

df <- data.frame(
  group = c('clustered', 'random'),
  name = c(
    'good dispersers \n large animals \n generalists',
    'poor dispersers \n small animals \n specialists'
  ),
  x = c(1, 2),
  y = c(2, 3)
) %>%
  mutate(name = factor(
    name,
    levels = c(
      'good dispersers \n large animals \n generalists',
      'poor dispersers \n small animals \n specialists'
    )
  ))

theme_set(theme_half_open())

emptyplot <- ggplot(df, aes(x, y)) + facet_grid(name ~ group) +
  ggpubr::theme_pubr() +
  theme(
    strip.text.x  = element_text(size = 40),
    strip.text.y = element_text(size = 32),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.grid  = element_blank()
  )

figure1B <- ggdraw(emptyplot) +
  draw_image(
    plots[1],
    x = -0.085,
    y = 0.97,
    hjust = 0,
    vjust = 1,
    width = 0.6,
    height = 0.6
  ) +
  draw_image(
    plots[2],
    x = 0.38,
    y = 0.97,
    hjust = 0,
    vjust = 1,
    width = 0.6,
    height = 0.6
  ) +
  draw_image(
    plots[3],
    x = -0.085,
    y = 0.52,
    hjust = 0,
    vjust = 1,
    width = 0.6,
    height = 0.6
  ) +
  draw_image(
    plots[4],
    x = 0.38,
    y = 0.52,
    hjust = 0,
    vjust = 1,
    width = 0.6,
    height = 0.6
  )

# ggsave(filename= here("figures/figure1B.pdf"), width=15, height=8, plot = figure1B)

#### Figure 1C: Clustering coefficient / Transitivity ####
# Measures the probability that the adjacent vertices of a vertex are connected.
clustercoeff <-
  function(infile = "data/landscapes/landscape_1.csv", D = 0.3) {
    distmat <-
      infile %>% read_csv %>% as.matrix() %>% # read landscape data
      dist(.,
           method = "euclidean",
           diag = TRUE,
           upper = TRUE) %>% # make distance matrix
      as.matrix() # transform into a regular matrix
    cc <- ifelse(distmat < D, 1, 0) %>% # make it binary
      graph_from_adjacency_matrix(., mode = "undirected") %>% # transfrom distmat into a graph object
      transitivity(type = "global") # calcuate clustering coefficient from graph opject
    return(cc)
  }

#### Prepare landscape data ####
landscapes <- here("data/landscapes/*.csv") %>% Sys.glob() %>%
  as_tibble() %>% # transform vector into tibble
  mutate(landscape = basename(value) %>% numextract()) %>%
  arrange(landscape) %>%
  landscapes$type <- c(rep("random", 10), rep("clustered", 10))

landscapes <- landscapes %>% rowwise() %>%
  mutate(cc = clustercoeff(value, D = 0.3))

## Plot boxplots of clustering coeffiecients by landscape types
landscapes %>%
  ggplot(aes(x = as.factor(type), y = cc, alpha = 0.4)) +
  stat_summary(
    aes(colour = as.factor(type), fill = as.factor(type)),
    alpha = 0.4,
    show.legend = FALSE,
    position = position_dodge(1),
    fun.data = quantiles_95,
    geom = "boxplot"
  ) +
  scale_y_continuous(
    name = "Cluster coefficient\n",
    breaks = seq(0.4, 1, 0.2),
    limits = c(0.4, 1)
  ) +
  scale_x_discrete(name = "", breaks = "") +
  scale_colour_manual(values = c("blue", "darkgreen")) +
  scale_fill_manual(values = c("blue", "darkgreen")) +
  facet_grid(. ~ type, labeller = label_value) +
  ggpubr::theme_pubr() +
  theme(
    text = element_text(size = 20),
    axis.text.y.left = element_text(size = 28),
    strip.text.x = element_text(size = 48),
    axis.title.y.left = element_text(size = 48),
    axis.ticks.y.left = element_line(size = 1),
    axis.ticks.x = element_line(size = 1)
  )

# ggsave(here("figures/figure_1C_clustercoefficient.pdf"), width=12, height=8)

#### Figure 1C: Nearest neighbor distance ####
nndist <- function(infile = "data/landscapes/landscape_1.csv") {
  distmat <-
    infile %>% read_csv %>% as.matrix() %>% # read landscape data
    dist(method = "euclidean",
         diag = TRUE,
         upper = TRUE) %>% # make dstance matrix
    as.matrix() # transform into regular matrix
  is.na(distmat) <- !distmat # replace 0s with NAs
  
  tib <-
    as.numeric(apply(distmat, 1, min, na.rm = TRUE)) %>% # extract nearest neighbor distance for each patch
    as.tibble() %>% # transform into tibble
    rename(nndist = value) %>% # rename variable
    nest() # nest vector of nndists
  return(tib)
}

#### Prepare landscape data ####
landscapes_long <- landscapes %>%
  rowwise() %>% # for each landscape
  mutate(nndist(value)) %>% # get nearest neighbor dsitances
  unnest_legacy(data) # unnest nearest neighbor distances

## Plot histograms of nearest neighbor distances
landscapes_long %>% group_by(type) %>%
  summarise(grpmean = mean(nndist)) %>% # group mean for each landscape type
  inner_join(landscapes_long, by = "type") %>%
  ggplot(aes(x = nndist)) +
  geom_histogram(
    binwidth = 0.001,
    show.legend = FALSE,
    aes(colour = as.factor(type), fill = as.factor(type)),
    alpha = 0.4
  ) +
  geom_vline(
    aes(xintercept = grpmean, color = as.factor(type)),
    linetype = "dashed",
    show.legend = FALSE
  ) +
  geom_density(alpha = 0.1,
               show.legend = FALSE,
               aes(colour = as.factor(type), fill = as.factor(type))) +
  scale_colour_manual(values = c("blue", "darkgreen")) +
  scale_fill_manual(values = c("blue", "darkgreen")) +
  facet_grid(. ~ type, labeller = label_value) +
  labs(y = "Count", x = "Nearest neighbour distance") +
  ggpubr::theme_pubr() +
  theme(
    text = element_text(size = 20),
    axis.text.x = element_text(size = 28),
    axis.title.x = element_text(size = 48),
    axis.title.y.left = element_text(size = 48),
    axis.text.y.left = element_text(size = 28),
    axis.ticks.y.left = element_line(size = 1),
    axis.ticks.x = element_line(size = 1),
    strip.text.x =  element_text(size = 48),
    legend.title = element_blank(),
    legend.key.height = unit(2, "cm"),
    legend.key.width  = unit(1, "cm")
  )

# ggsave(here("figures/figure_1C_nndist.pdf"), width=16, height=8)
