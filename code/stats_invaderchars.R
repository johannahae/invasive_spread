

## Table S1: Statistitcs for figure 3 ####

# Load libraries
library(here)
library(lme4)
require(ggpmisc)
require(broom)
require(mblm)
library(broom.mixed)
library(tidyverse)

# Declare path
here::i_am("invasive_spread.Rproj")
# here::here()

# Import data
dat_plot <- here("data/results/dat_plot.rds") %>% read_rds()
dat_plot <- dat_plot %>% arrange(landscape.type)


# Vulnerability
glmer_per_frac <- function(df)
  glmer(
    data = df,
    formula = fraction.of.invaded.patches / 100 ~ scale(vulnerability.invspp.int) + (1 |
                                                                                       web),
    family = binomial(link = "logit")
  )

model <- glmer_per_frac(dat_plot)

tab_stats <- dat_plot %>%
  group_by(landscape.type, nutrient.scenario) %>%
  nest() %>%
  ungroup() %>%
  mutate(stats = map(data, glmer_per_frac)) %>%
  mutate(params = map(stats, broom.mixed::tidy)) %>%
  unnest(params) %>%
  mutate(spchar = "Inv - Vulnerability")

# Generality
glmer_per_frac <- function(df)
  glmer(
    data = df,
    formula = fraction.of.invaded.patches /
      100 ~ scale(generality.invspp.int) + (1 | web),
    family = binomial(link = "logit")
  )

tab_stats <- dat_plot %>%
  group_by(landscape.type, nutrient.scenario) %>%
  nest() %>%
  ungroup() %>%
  mutate(stats = map(data, glmer_per_frac)) %>%
  mutate(params = map(stats, broom.mixed::tidy)) %>%
  unnest(params) %>%
  mutate(spchar = "Inv - Generality") %>%
  bind_rows(tab_stats)


# Omnivory index
glmer_per_frac <- function(df)
  glmer(
    data = df,
    formula = fraction.of.invaded.patches /
      100 ~ scale(omnivory.index.invspp.int) + (1 | web),
    family = binomial(link = "logit")
  )

tab_stats <- dat_plot %>%
  group_by(landscape.type, nutrient.scenario) %>%
  nest() %>%
  ungroup() %>%
  mutate(stats = map(data, glmer_per_frac)) %>%
  mutate(params = map(stats, broom.mixed::tidy)) %>%
  unnest(params) %>%
  mutate(spchar = "Inv - Omnivory index") %>%
  bind_rows(tab_stats)

# Trophic level
glmer_per_frac <- function(df)
  glmer(
    data = df,
    formula = fraction.of.invaded.patches /
      100 ~ scale(trophic.level.invspp.int.rounded) + (1 | web),
    family = binomial(link = "logit")
  )

tab_stats <- dat_plot %>%
  group_by(landscape.type, nutrient.scenario) %>%
  nest() %>%
  ungroup() %>%
  mutate(stats = map(data, glmer_per_frac)) %>%
  mutate(params = map(stats, broom.mixed::tidy)) %>%
  unnest(params) %>%
  mutate(spchar = "Inv - Trophic level") %>%
  bind_rows(tab_stats)

# Dispersal distance
glmer_per_frac <- function(df)
  glmer(
    data = df,
    formula = fraction.of.invaded.patches /
      100 ~ scale(dispersal.dist.invspp) + (1 | web),
    family = binomial(link = "logit")
  )

tab_stats <- dat_plot %>%
  group_by(landscape.type, nutrient.scenario) %>%
  nest() %>%
  ungroup() %>%
  mutate(stats = map(data, glmer_per_frac)) %>%
  mutate(params = map(stats, broom.mixed::tidy)) %>%
  unnest(params) %>%
  mutate(spchar = "Inv - Dispersal distance") %>%
  bind_rows(tab_stats)

# Log10 Bodymass
glmer_per_frac <-
  function(df)
    glmer(
      data = df,
      formula = fraction.of.invaded.patches / 100 ~ scale(log10(bodymass.invspp)) + (1 |
                                                                                       web),
      family = binomial(link = "logit")
    )

df = dat_plot %>% filter(landscape.type == "random", nutrient.scenario == "oligotrophic")
glmer_per_frac(df) %>% summary

tab_stats <- dat_plot %>%
  group_by(landscape.type, nutrient.scenario) %>%
  nest() %>%
  ungroup() %>%
  mutate(stats = map(data, glmer_per_frac)) %>%
  mutate(params = map(stats, broom.mixed::tidy)) %>%
  unnest(params) %>%
  mutate(spchar = "Inv - Body mass [log10]") %>%
  bind_rows(tab_stats)

tab_stats_sub <- tab_stats %>%
  filter(effect == "fixed", term != "(Intercept)") %>%
  dplyr::select(
    c(
      "spchar",
      "landscape.type",
      "nutrient.scenario",
      "estimate",
      "statistic",
      "p.value"
    )
  )
# select("statistic")

# print(xtable::xtable(tab_stats_sub, digits = 5), include.rownames=FALSE) # latex table with stats
