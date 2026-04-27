# Epidemiological model of interactions between two seasonal respiratory viruses
# Copyright (C) 2025 George Shirreff
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

library(tidyverse)
library(stringi)
library(ggh4x)
library(odin)

suma = function(x){
  if (all(is.na(x))){
    x[NA_integer_]
  } else {
    sum(x, na.rm = TRUE) 
  }
} 
season_week <- function(season, week){
  season_week0 <- paste0("20", stri_sub(season, from = -2), "-01-2") %>% as.Date("%Y-%W-%w") - 4
  return(season_week0 + week*7)
}

source("odin_m5.R")

param_description <- read_csv2("input/param_description.csv")

# read data
allcatch_inc_age <- read_tsv("data/allcatch_inc_age.tsv", show_col_types = F)
allcatch_coinc_age <- read_tsv("data/allcatch_coinc_age.tsv", show_col_types = F)
# spenc <- read_csv("data/Spencer_params.csv")
# age_contact_standardised <- read_tsv("data/age_contact_standardised.tsv")
ili <- read_tsv("data_public/allcatch_ili_age2.tsv", show_col_types = F) %>% 
  filter(!(season %in% c("20/21"))) %>% 
  rename(age2_group = age_group, t = iso_standard) %>% 
  select(-virus, -CatchmentPop)

census2022 <- read_csv("data_public/census2022.csv", show_col_types = F)



dat <- allcatch_inc_age %>% 
  filter(#virus == this_virus, 
    # age_group == this_age,
    # age_group == "6-12m", 
    !(season %in% c("20/21")), 
    # season == this_season
  )

codat <- allcatch_coinc_age %>% 
  filter(#virus == this_virus, 
    # age_group == this_age,
    # age_group == "6-12m", 
    !(season %in% c("20/21")), 
    # season == this_season
  )




age_vec = dat$age_group %>% unique
age2_vec = ili$age2_group %>% unique
season_vec = dat$season %>% unique
n_seasons = season_vec %>% length
n_ages = age_vec %>% length


virus1 = "Influenza"
virus2 = "RSV"

dat_2virus <- rbind(dat %>% filter(virus == virus1) %>% mutate(virus = 1), 
                    dat %>% filter(virus == virus2) %>% mutate(virus = 2), 
                    codat %>% filter(CoVirus == paste0(virus1, "+", virus2)) %>% mutate(CoVirus = 12) %>% rename(virus = CoVirus)) %>% 
  transmute(season, t = iso_standard, age_group, virus, Inc, CatchmentPop) %>% 
  arrange(season, t, match(age_group, age_vec), virus)


data_catchments_weights = dat_2virus %>% 
  select(age_group, season, CatchmentPop) %>% unique %>% 
  left_join(census2022, by = "age_group") %>% 
  mutate(age2_group = factor(age2_age(age_group), levels = age2_vec)) %>% 
  group_by(season, age2_group) %>% 
  mutate(age2_weight = CensusPopulation/sum(CensusPopulation)) %>% 
  ungroup


# allcatch_inc_age <- read_tsv("data/allcatch_inc_age.tsv")
# dat <- allcatch_inc_age %>% 
#   filter(#virus == this_virus, 
#     # age_group == this_age,
#     # age_group == "6-12m", 
#     !(season %in% c("20/21")), 
#     # season == this_season
#   )
# 
# ili <- read_tsv("data_public/allcatch_ili_age2.tsv") %>% 
#   filter(!(season %in% c("20/21"))) %>% 
#   rename(age2_group = age_group, t = iso_standard) %>% 
#   select(-virus, -CatchmentPop)
# 
# age_vec = dat$age_group %>% unique
# age2_vec = ili$age2_group %>% unique
# season_vec = dat$season %>% unique
# n_seasons = season_vec %>% length
# n_ages = age_vec %>% length
# 
# 
# dat_2virus <- rbind(dat %>% filter(virus == virus1) %>% mutate(virus = 1), 
#                     dat %>% filter(virus == virus2) %>% mutate(virus = 2), 
#                     codat %>% filter(CoVirus == paste0(virus1, "+", virus2)) %>% mutate(CoVirus = 12) %>% rename(virus = CoVirus)) %>% 
#   transmute(season, t = iso_standard, age_group, virus, Inc, CatchmentPop) %>% 
#   arrange(season, t, match(age_group, age_vec), virus)

#### FIGURE 1 ####

include_ci = T

dat %>% 
  filter(virus %in% c("Influenza", "RSV")) %>% 
  mutate(age2_group = age2_age(age = age_group)) %>% 
  group_by(virus, season, iso_standard, age2_group) %>% 
  summarise(Inc = suma(Inc), 
            CatchmentPop = sum(CatchmentPop)) %>% 
  mutate(IncRate1e5 = Inc/CatchmentPop*100000) %>%
  # ungroup %>% tidyr::complete(virus, season, iso_standard, age2_group, fill = list(Inc = NA, IncRate1e5 = NA, CatchmentPop = NA)) %>% 
  select(-Inc) %>% 
  rbind(ili %>% transmute(virus = "ILI", age2_group, season, iso_standard = t, IncRate1e5, CatchmentPop = 1e5)) %>% 
  mutate(age2_group = factor(age2_group, levels = age2_vec)) %>% 
  mutate(midweek = season_week(season, iso_standard)) %>% 
  # mutate(rowlabel = ifelse(virus == "ILI", "ILI", "Hospitalisation")) %>% 
  mutate(virus = factor(virus, levels = c("Influenza", "RSV", "ILI"))) %>% 
  ggplot(aes(x = midweek, y = IncRate1e5, colour = age2_group)) + 
  geom_line(linewidth = 0.8) + 
  facet_grid(virus ~ ., scales = "free_y") + 
  theme_bw() + 
  scale_colour_brewer(palette = "Spectral") + 
  scale_x_date(date_breaks = "year", date_labels = "%Y", 
               limits = as.Date(c("2010-07-01", "2020-07-01")), expand = c(0,0)) + 
  labs(y = "Incidence per 100,000 per week", x = "", colour = "Age group (years)")


# ggsave("figtab/Fig3.png", width = 25, height = 18, units = "cm")
ggsave("figtab/Fig1_incidences.png", width = 25, height = 18, units = "cm")

dat %>% 
  filter(virus %in% c("Influenza", "RSV")) %>% 
  mutate(IncRate1e5 = Inc/CatchmentPop*100000) %>%
  mutate(age_group = factor(age_group, levels = age_vec)) %>% 
  mutate(midweek = season_week(season, iso_standard)) %>% 
  ggplot(aes(x = midweek, y = IncRate1e5, colour = virus)) + 
  geom_line(linewidth = 0.8) + 
  facet_nested(age_group ~ ., scales = "free_y") + 
  theme_bw() + 
  scale_colour_manual(values = c(RSV="blue", Influenza="orange3")) + 
  # scale_colour_manual(values = RColorBrewer::brewer.pal(3, "Set1")[1:2] %>% setNames(c("RSV", "Influenza"))) + 
  
  scale_x_date(date_breaks = "year", date_labels = "%Y", 
               limits = as.Date(c("2010-07-01", "2020-07-01")), expand = c(0,0)) + 
  # scale_y_continuous(breaks = seq(10, 600, by = 10), minor_breaks = seq(5, 605, by = 10)) + 
  labs(y = "Incidence per 100,000 per week", x = "", colour = "Virus")
  
# ggsave("figtab/FigS3.png", width = 25, height = 18, units = "cm")
ggsave("figtab/FigS3_incidencesAge.png", width = 25, height = 18, units = "cm")




