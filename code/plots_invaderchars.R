

#### Preface ####

# Declare libraries
library(here)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(viridis)
library(RColorBrewer)
library(lme4)
require(ggpmisc)
require(broom)
require(mblm)
library(broom.mixed)

# Declare path
here::i_am("invasive_spread.Rproj")
# here::here()

# Define and set color palettes
cbPalette <-
  c("#999999",
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#0072B2",
    "#D55E00",
    "#CC79A7") # palette with grey
gg_color_hue <- function(n)
  cbPalette[1:n]

pal <- gg_color_hue(7)
ramp <- colorRampPalette(pal, bias = 0.5)

# Define and set theme
theme_half_open(
  font_size = 14,
  font_family = "",
  line_size = 0.5,
  rel_small = 11 / 14,
  rel_tiny = 12 / 14,
  rel_large = 16 / 14
)

theme_set(theme_half_open() + background_grid())

# Import data and functions
dat_plot <- here("data/results/dat_plot.rds") %>% read_rds()
source(here("code/boxplots.R"))

## Figure 3: Script to plot effects of invader characteristics.
# Boxplots with rows indicating the landscape type (clustered/random);
# columns nutrient availability.

# Fig. 3A: Body mass ####
datin <- dat_plot %>% drop_na(bodymass.invspp)
datin$x <-
  datin %>% pull(bodymass.invspp) %>% log10 %>% plyr::round_any(., 2, ceiling)
xu <- datin$x %>% unique
pal <- gg_color_hue(xu %>% length)

w <-
  datin %>% # width of boxplots depending on number of observations in each group
  group_by(landscape.type, nutrient.scenario, x) %>%
  tally() %>% pull(n) %>% sqrt / 10

bodymass <-
  biotic.boxplot.left(
    dat = datin,
    x = x,
    xu = xu,
    pal = pal,
    barwidth = w,
    namex = "Inv - Body mass [log10]",
    xmin = 1.5,
    xmax = 13,
    mylabels = scales::math_format(10 ^ .x)
  )

a <- annotation_logticks(base = 10,
                         sides = "b",
                         scaled = TRUE)
a$data <- data.frame(x = NA, landscape.type = "random")
bodymass.log <- bodymass + a

# bodymass.log %>% ggsave(filename = here("figures/figure_3A.pdf"), width=16, height=12)

# Figure 3B: Dispersal distance ####
datin <- dat_plot %>% drop_na(dispersal.dist.invspp)
datin$x <- datin %>% pull(dispersal.dist.invspp) %>% round(1) * 10
datin <- datin %>% filter(x < 7)
xu <- datin$x %>% unique %>% sort
pal <- gg_color_hue(xu %>% length)

w <-
  datin %>% #width of boxplots depending on number of observations in each group
  group_by(landscape.type, nutrient.scenario, x) %>%
  tally() %>% pull(n) %>% sqrt / 20

dispersal <-
  biotic.boxplot.right(
    dat = datin,
    x = x,
    xu = xu,
    pal = pal,
    barwidth = w,
    namex = "Inv - Dispersal distance",
    xmin = -0.5,
    xmax = 5,
    mylabels = xu / 10
  )

# dispersal %>% ggsave(filename = here("figures/figurue_3B.pdf"), width=16, height=12)

# Figure 3C: Trophic level ####
datin <- dat_plot %>% drop_na(trophic.level.invspp.int.rounded)
datin$x <- datin %>% pull(trophic.level.invspp.int.rounded)
xu <- unique(datin$x)

w <-
  datin %>% # width of boxplots depending on number of observations in each group
  group_by(landscape.type, nutrient.scenario, x) %>%
  tally() %>% pull(n) %>% sqrt / 20

pal <- ramp(xu %>% length)
trophiclevel <-
  biotic.boxplot.left(
    dat = datin,
    x,
    xu = xu,
    pal = pal,
    barwidth = w,
    namex = "Inv - Trophic level",
    xmin = 0.5,
    xmax = 7,
    mylabels = xu
  )

# trophiclevel %>% ggsave(filename = here("figures/figure_3C.pdf"), width=16, height=12)

# Figure 3D: Omnivory ####
datin <- dat_plot %>% drop_na(omnivory.index.invspp.int)
datin$x <-  datin %>% pull(omnivory.index.invspp.int) %>%
  plyr::round_any(., 0.3, ceiling) * 10
xu <- datin$x %>% unique
xu %>% range
pal <- gg_color_hue(xu %>% length)

w <-
  datin %>% #width of boxplots depending on number of observations in each group
  group_by(landscape.type, nutrient.scenario, x) %>%
  tally() %>% pull(n) %>% sqrt / 10

omnivory <-
  biotic.boxplot.right(
    dat = datin,
    x = x,
    xu = xu,
    pal = pal,
    barwidth = w,
    namex = "Inv - Omnivory index",
    xmin = -1.5,
    xmax = 13,
    mylabels = xu / 10
  )

# omnivory %>% ggsave(filename = here("figures/figure_3D.pdf"), width=16, height=12)

# Figure 3E: Generality ####
datin <- dat_plot %>% drop_na(generality.invspp.int)
datin$x <- datin %>% pull(generality.invspp.int) %>%
  divide_by(datin %>% pull(gamma.int)) %>%
  round(1) * 10
datin <- datin %>% filter(x < 7)
xu <- datin$x %>% unique %>% sort
pal <- gg_color_hue(xu %>% length)

w <-
  datin %>% #width of boxplots depending on number of observations in each group
  group_by(landscape.type, nutrient.scenario, x) %>%
  tally() %>% pull(n) %>% sqrt / 20

generality <-
  biotic.boxplot.left(
    dat = datin,
    x = x,
    xu = xu,
    pal = pal,
    barwidth = w,
    namex = "Inv - Generality (prey counts)",
    xmin = -0.75,
    xmax = 7,
    mylabels = xu / 10
  )

# generality %>% ggsave(filename = here("figures/figure_3E.pdf"), width=16, height=12)

# Figure 3F: Vulnerability ####
datin <- dat_plot %>% drop_na(vulnerability.invspp.int)
datin$x <- datin %>% pull(vulnerability.invspp.int) %>%
  divide_by(datin %>% pull(gamma.int)) %>%
  round(1) * 10
datin <- datin %>% filter(x < 7)
xu <- datin$x %>% unique %>% sort
pal <- gg_color_hue(xu %>% length)

w <-
  datin %>% #width of boxplots depending on number of observations in each group
  group_by(landscape.type, nutrient.scenario, x) %>%
  tally() %>% pull(n) %>% sqrt / 20

vulnerability <-
  biotic.boxplot.right(
    dat = datin,
    x = x,
    xu = xu,
    pal = pal,
    barwidth = w,
    namex = "Inv - Vulnerability (predator counts)",
    xmin = -0.75,
    xmax = 6,
    mylabels = xu / 10
  )

# vulnerability %>% ggsave(filename = here("figures/figure_3F.pdf"), width=16, height=12)
