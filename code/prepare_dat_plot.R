
# Declare libraries
require(tidyverse)
require(here)

# Import data
dat <-
  here("data/results/summary.csv") %>% read_csv(col_types = cols())

dat %>% nrow # total number of 6300 simulation runs (invasion processes)

# Prepare and modify data for statistical analysis and to generate boxplots
dat_plot <- dat %>%
  filter(INVASION == 1) %>%
  mutate(
    fraction.of.invaded.patches = fraction.of.invaded.patches * 100,
    rem.spp = rem.spp + 1,
    landscape.type = # transform factor
      factor(landscape.type, levels = c('clustered', 'random')),
    feeding.type.invspp = # transform to factor
      factor(
        feeding.type.invspp,
        levels = c('basal', 'herbivore', 'omnivore', 'carnivore')
      ),
    nutrient.scenario = # define nutrient scenarios
      case_when(
        .$nutrient.supply == 0.1 ~ "oligotrophic",
        .$nutrient.supply == 1 ~ "mesotrophic",
        .$nutrient.supply == 1000 ~ "eutrophic"
      ),
    nutrient.scenario = # transform to factor
      factor(
        nutrient.scenario,
        levels = c("oligotrophic", "mesotrophic", "eutrophic")
      ),
    trophic.level.invspp.int.rounded = trophic.level.invspp.int %>% round,
    success = case_when(.$number.of.invaded.patches > 0 ~ 1, TRUE ~ 0),
    # create binary success variable
    log10mass = bodymass.invspp %>% log10(),
    log10tott0 =  tot.biomass.native.t0 %>% log10(),
    log10totpost = tot.biomass.native.post %>% log10(),
    log10invsppint = tot.biomass.invspp.int %>% log10(),
    log10invspppost = tot.biomass.invspp.post %>% log10(),
    normG = generality.invspp.int %% gamma.int,
    normV = vulnerability.invspp.int %% gamma.int,
    webno = case_when(
      .$web == 7516 ~ "food web 3",
      .$web == 8516 ~ "food web 4",
      .$web == 41516 ~ "food web 1",
      .$web == 42516 ~ "food web 5",
      .$web == 88516 ~ "food web 2"
    )
  ) %>%
  rename(`Nutrient supply rate` = nutrient.supply)

dat_plot <- dat_plot %>%
  filter(invspp.type == "consumer species" &
           trophic.level.invspp.int >= 2) %>% # remove sims in which all basal sp. went extinct
  bind_rows(dat_plot %>% filter(invspp.type == "basal species")) %>%
  drop_na(trophic.level.invspp.int) %>% # remova NAs
  filter(gamma.int > 0) # remove sims in which all species went extinct

dat_plot %>% nrow # 5185 simulations remain after removing non-meaningful output

# dat_plot %>% saveRDS(here("data/results/dat_plot.rds")) # save as rds-file to preserve column types
