

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

bg <- c("blue", "darkgreen")

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

## Script to plot effects of abiotic factors for manuscript and SI.
# Boxplots with landscape types (clustered/random) on the x-axis;
# columns indicate nutrient level

# Figure 2A: Fraction of successful invasions #####

# 1st, calculate fraction of successful invasions per group
n <- dat_plot %>%
  group_by(landscape.type, nutrient.scenario, web, landscape) %>%
  tally() %>% pull(n) # number of simulations in each group

sum_tib <-
  dat_plot %>% group_by(landscape.type, nutrient.scenario, web, landscape) %>%
  summarize(sumsuccess = sum(success)) %>%
  add_column(frac = .$sumsuccess / n * 100) %>% # fraction of successful invasions in each group
  ungroup

x <- sum_tib %>% pull(landscape.type) # x-var.
xu <- levels(x) # levels x-var.
y <- sum_tib %>% pull(frac)

success <-
  abiotic.boxplot(
    dat = sum_tib,
    x = x,
    xu = xu,
    y = y,
    pal = bg,
    # make boxplot
    namey = "Fraction successful invasions"
  )
# success %>% ggsave(filename = here("figures/figure_2A.pdf"), width=16, height=10)

# make color legend for landscape types:
# landscape_legend <- get_legend(success) %>% as_ggplot() # extract legend
# to do this, we must set the show.legend=TRUE in the function "abiotic.boxplot"
# ggsave(filename = here("figures/landscape_legend.pdf"), width=6, height=2, legend)

## Figure S7A: Fraction of successful invasions over all nutrient levels ###
totmean <- dat_plot %>%
  mutate(nrows = nrow(.)) %>%
  summarize(sumsuccess = sum(success), nrow = unique(nrows)) %>%
  mutate(frac = divide_by(sumsuccess, nrow)) # total fraction of successful invasions

n <- dat_plot %>%
  group_by(landscape.type) %>% # groups
  tally() %>% pull(n) # number of simulations per landscape type

lansdcmeans <- dat_plot %>% group_by(landscape.type) %>%
  summarize(sumsuccess = sum(success)) %>%
  mutate(frac = divide_by(sumsuccess, n)) # fraction of successful invasions per landscape type

success.allnuts <-
  abiotic.boxplot.allnuts(
    dat = sum_tib,
    x = x,
    xu = xu,
    y = y,
    pal = bg,
    namey = "Fraction successful invasions"
  ) +
  stat_compare_means(method = "anova",
                     label.x = 0.5,
                     label.y = 105) # add p-value anova
# success.allnuts %>% ggsave(filename = here("figures/figure_S7A.pdf"), width=16, height=10)

## Figure 2B: Fraction of invaded patches #####
x <- dat_plot %>% pull(landscape.type)
xu <- x %>% levels()
y <- dat_plot %>% pull(fraction.of.invaded.patches)

spread <-
  abiotic.boxplot(
    dat = dat_plot,
    x = x,
    xu = xu,
    y = y,
    pal = bg,
    namey = "Fraction invaded patches"
  )
# spread %>% ggsave(filename = here("figures/figure_2B.pdf"), width=16, height=10)

## Figure S7B: Fraction of invaded patches over all nutrient levels ####
spread.allnuts <-
  abiotic.boxplot.allnuts(
    dat = dat_plot,
    x = x,
    xu = xu,
    y = y,
    pal = bg,
    namey = "Fraction invaded patches"
  ) +
  stat_compare_means(method = "anova",
                     label.x = 0.5,
                     label.y = 105) # add p-value anova

# spread.allnuts %>% ggsave(filename = here("figures/figure_S7B.pdf"), width=16, height=10)

# Figure 2C: Total biomass density ####
x <- dat_plot %>% pull(landscape.type)
xu <- levels(x)

ytot <-
  dat_plot %>% pull(tot.biomass.invspp.int) +  dat_plot %>% pull(tot.biomass.native.int)
y <- dat_plot %>% pull(tot.biomass.invspp.int) / ytot
y <- y + 1 %>% log10
ytot <- ytot + 1 %>% log10

totbiomass <-
  abiotic.boxplot.log(
    dat = dat_plot,
    x = x,
    xu = xu,
    y = ytot,
    pal = bg,
    namey = "Total biomass density"
  )
# totbiomass %>% ggsave(filename = here("figures/figure_2C.pdf"), width=16, height=10)

# Figure 2D: Relative invader biomass density ####
relbiomass <-
  abiotic.boxplot.log(
    dat = dat_plot,
    x = x,
    xu = xu,
    y = y,
    pal = bg,
    namey = "Relative Inv-biomass density"
  )
# relbiomass %>% ggsave(filename = here("figures/figure_2Dd.pdf"), width=16, height=10)
