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

#### FIGURE 2 ####

allcatch_inc_age <- read_tsv("data/allcatch_inc_age.tsv")

include_ci = T

dat <- allcatch_inc_age %>% 
  filter(#virus == this_virus, 
    # age_group == this_age,
    # age_group == "6-12m", 
    !(season %in% c("20/21")), 
    # season == this_season
  )

ili <- read_tsv("data/allcatch_ili_age2.tsv") %>% 
  filter(!(season %in% c("20/21"))) %>% 
  rename(age2_group = age_group, t = iso_standard) %>% 
  select(-virus, -CatchmentPop)

age_vec = dat$age_group %>% unique
age2_vec = ili$age2_group %>% unique
season_vec = dat$season %>% unique
n_seasons = season_vec %>% length
n_ages = age_vec %>% length


dat_2virus <- rbind(dat %>% filter(virus == virus1) %>% mutate(virus = 1), 
                    dat %>% filter(virus == virus2) %>% mutate(virus = 2), 
                    codat %>% filter(CoVirus == paste0(virus1, "+", virus2)) %>% mutate(CoVirus = 12) %>% rename(virus = CoVirus)) %>% 
  transmute(season, t = iso_standard, age_group, virus, Inc, CatchmentPop) %>% 
  arrange(season, t, match(age_group, age_vec), virus)



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
ggsave("figtab/Fig2_incidences.png", width = 25, height = 18, units = "cm")

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

#### Figure 4 ####

best_params <- read_csv(paste0("output/", id = "validateSimpleInter_", virus1, "+", virus2, "_jump", 3, "_rep", 1, "_bestParams.csv"))

model_out <- generate_ModelHosp_season(pars = best_params %>% 
                                         select(-Iteration, -Loglik) %>% 
                                         unlist, generator = m5_generator, dat_2virus, data_catchments_weights = data_catchments_weights)

tib_ModelHosp <- model_out$hosp

plot_data <- dat_2virus %>% 
  transmute(season, age_group = factor(age_group, levels = age_vec), t, 
            virus = c(virus1, virus2, "Codetection")[match(virus, c(1, 2, 12))], 
            DataHosp = Inc, 
            ModelHosp = tib_ModelHosp$Inc)
# %>% 
#   mutate(DataHosp = DataHosp/CatchmentPop*1e5, 
#          ModelHosp = ModelHosp/CatchmentPop*1e5)

plot_data_date <- plot_data %>% 
  mutate(W1 = ISOweek::ISOweek2date(paste0("20", gsub(".*/", "", season), "-W01-1"))) %>% 
  mutate(week_begins = W1 - 7 + t*7) %>% 
  group_by(age_group, virus, season, week_begins) %>% 
  summarise(across(ends_with("Hosp"), ~sum(.x, na.rm = T))) %>% 
  ungroup %>% 
  complete(age_group, virus, season, week_begins, fill = list(DataHosp = NA, ModelHosp = NA))


catchment_season <- dat_2virus %>% 
  filter(virus == 1, t == 0) %>% 
  group_by(season) %>% 
  summarise(CatchmentPop = sum(CatchmentPop))

# plot_obj <- plot_data %>% 

# group_by(season, t, virus) %>% 
plot_obj <- plot_data_date %>% 
  mutate(virus = factor(virus, levels = c("Influenza", "RSV", "Codetection"))) %>% 
  mutate(virus = fct_recode(virus, `Influenza only`="Influenza", `RSV only`="RSV")) %>% 
  group_by(week_begins, virus) %>% 
  summarise(across(ends_with("Hosp"), ~sum(.x, na.rm = T)), season = season[1]) %>% 
  mutate(model_loCI = qpois(0.025, ModelHosp, lower.tail = T), 
         model_hiCI = qpois(0.025, ModelHosp, lower.tail = F)) %>%  
  left_join(catchment_season) %>%
  mutate(across(c("DataHosp", "ModelHosp", "model_loCI", "model_hiCI"), ~.x/CatchmentPop*1e5)) %>%
  ggplot(aes(x = week_begins)) + 
  {if(include_ci) geom_ribbon(aes(ymin = model_loCI, ymax = model_hiCI, fill = "CI"), show.legend = FALSE)} + 
  
  geom_line(aes(y = ModelHosp,  linetype = "Best model"), colour = "blue") + 
  geom_point(aes(y = DataHosp, shape = "Data"), colour = "red", size = 0.8) +
  facet_grid(virus~., scales = "free_y") + 
  labs(x = "Week", y = "Cases per 100 000", colour = "", linetype = "", shape = "") + 
  theme_bw()+ 
  scale_fill_manual(values = c(CI = "lightblue")) + 
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0))

# ggsave(plot_obj, filename = "figtab/Fig4.png", width = 25, height = 15, units = "cm", dpi = 600)
ggsave(plot_obj, filename = "figtab/Fig4_dataBestmodel.png", width = 25, height = 15, units = "cm", dpi = 600)
  

#### FIGURE S4

plot_data_date %>% 
  filter(virus == "Codetection") %>% 
  group_by(week_begins, age_group) %>% 
  summarise(across(ends_with("Hosp"), ~sum(.x, na.rm = T)), season = season[1]) %>% 
  mutate(model_loCI = qpois(0.025, ModelHosp, lower.tail = T), 
         model_hiCI = qpois(0.025, ModelHosp, lower.tail = F)) %>%  
  left_join(catchment_season) %>%
  mutate(across(c("DataHosp", "ModelHosp", "model_loCI", "model_hiCI"), ~.x/CatchmentPop*1e5)) %>%
  ggplot(aes(x = week_begins)) + 
  {if(include_ci) geom_ribbon(aes(ymin = model_loCI, ymax = model_hiCI, fill = "CI"), show.legend = FALSE)} + 
  
  geom_line(aes(y = ModelHosp,  linetype = "Best model"), colour = "blue") + 
  geom_point(aes(y = DataHosp, shape = "Data"), colour = "red", size = 0.8) +
  facet_grid(age_group~.
             # , scales = "free_y"
             ) + 
  labs(x = "Week", y = "Cases per 100 000", colour = "", linetype = "", shape = "") + 
  theme_bw()+ 
  coord_cartesian(ylim = c(0, 0.22)) + 
  scale_fill_manual(values = c(CI = "lightblue")) + 
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0))

ggsave("figtab/FigS4_dataBestmodelAge.png", width = 25, height = 18, units = "cm")

#### FIGURE S5 ILI

model_out <- generate_ModelHosp_season(pars = best_params %>% 
                                         select(-Iteration, -Loglik) %>% 
                                         unlist, generator = m5_generator, dat_2virus, data_catchments_weights = data_catchments_weights)

tib_ModelILI <- model_out$ili


plot_ili <- ili %>% 
  transmute(season, age2_group = factor(age2_group, levels = age2_vec), t, 
            DataILI = IncRate1e5, 
            ModelILI = tib_ModelILI$IncRate1e5)

plot_ili_date <- plot_ili %>% 
  mutate(W1 = ISOweek::ISOweek2date(paste0("20", gsub(".*/", "", season), "-W01-1"))) %>% 
  mutate(week_begins = W1 - 7 + t*7) %>% 
  ungroup %>% 
  complete(age2_group, week_begins, fill = list(DataILI = NA, ModelILI = NA))


plot_ili_obj <- plot_ili_date %>% 
  mutate(model_loCI = qlnorm(0.025, meanlog = log(ModelILI + 1), sdlog = sd(log(ili$IncRate1e5 + 1), na.rm = T), lower.tail = T), 
         model_hiCI = qlnorm(0.025, meanlog = log(ModelILI + 1), sdlog = sd(log(ili$IncRate1e5 + 1), na.rm = T), lower.tail = F)) %>%  
  ggplot(aes(x = week_begins)) + 
  {if(include_ci) geom_ribbon(aes(ymin = model_loCI, ymax = model_hiCI, fill = "CI"), show.legend = FALSE)} +
  
  geom_line(aes(y = ModelILI,  linetype = "Best model"), colour = "blue") + 
  geom_point(aes(y = DataILI, shape = "Data"), colour = "red", size = 0.8) +
  facet_grid(age2_group~., scales = "free_y") + 
  labs(x = "Week", y = "ILI per 100 000", colour = "", linetype = "", shape = "") + 
  theme_bw()+ 
  scale_fill_manual(values = c(CI = "lightblue")) + 
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0)) + 
  coord_cartesian(ylim = c(0, 800))

ggsave(plot_ili_obj, filename = "figtab/FigS5.png", width = 25, height = 20, units = "cm", dpi = 600)


#### FIGURE S5 do it but maybe cut it???

best_params <- read_csv(paste0("output/", id = "validateSimpleInter_", virus1, "+", virus2, "_jump", 3, "_rep", 1, "_bestParams.csv"))
baseline_params <- read_csv(paste0("output/", id = "simpleInter_", virus1, "+", virus2, "_jump", 0, "_rep", 0, "_bestParams.csv"))
nextbest_params <- read_csv(paste0("output/", id = "simpleInter_", virus1, "+", virus2, "_jump",5, "_rep", 2, "_bestParams.csv"))

tib_ModelHosp_best <- generate_ModelHosp_season(pars = best_params %>% 
                                         select(-Iteration, -Loglik) %>% 
                                         unlist, generator = m5_generator, dat_2virus = dat_2virus, data_catchments_weights = data_catchments_weights)$hosp

tib_ModelHosp_baseline <- generate_ModelHosp_season(baseline_params %>% 
                                                  select(-Iteration, -Loglik) %>% 
                                                  unlist, m5_generator, dat_2virus, data_catchments_weights = data_catchments_weights)$hosp

tib_ModelHosp_nextbest <- generate_ModelHosp_season(nextbest_params %>% 
                                                      select(-Iteration, -Loglik) %>% 
                                                      unlist, m5_generator, dat_2virus, data_catchments_weights = data_catchments_weights)$hosp




plot_comp_data <- dat_2virus %>% 
  transmute(season, age_group = factor(age_group, levels = age_vec), t, 
            virus = c(virus1, virus2, "Codetection")[match(virus, c(1, 2, 12))], 
            DataHosp = Inc, 
            best_ModelHosp = tib_ModelHosp_best$Inc, 
            baseline_ModelHosp = tib_ModelHosp_baseline$Inc, 
            nextbest_ModelHosp = tib_ModelHosp_nextbest$Inc)
# %>% 
#   mutate(DataHosp = DataHosp/CatchmentPop*1e5, 
#          ModelHosp = ModelHosp/CatchmentPop*1e5)

plot_comp_data_date <- plot_comp_data %>% 
  mutate(W1 = ISOweek::ISOweek2date(paste0("20", gsub(".*/", "", season), "-W01-1"))) %>% 
  mutate(week_begins = W1 - 7 + t*7) %>% 
  group_by(age_group, virus, season, week_begins) %>% 
  summarise(across(ends_with("Hosp"), ~sum(.x, na.rm = T))) %>% 
  ungroup %>% 
  complete(age_group, virus, season, week_begins, fill = list(DataHosp = NA, best_ModelHosp = NA, baseline_ModelHosp = NA, nextbest_ModelHosp = NA)) %>% 
  pivot_longer(cols = ends_with("ModelHosp"), names_sep = "_", names_to = c("Model", ".value"))


catchment_season <- dat_2virus %>% 
  filter(virus == 1, t == 0) %>% 
  group_by(season) %>% 
  summarise(CatchmentPop = sum(CatchmentPop))

# plot_obj <- plot_data %>% 

# group_by(season, t, virus) %>% 
plot_comp_obj <- plot_comp_data_date %>% 
  mutate(virus = factor(virus, levels = c("Influenza", "RSV", "Codetection"))) %>% 
  mutate(virus = fct_recode(virus, `Influenza only`="Influenza", `RSV only`="RSV")) %>% 
  group_by(week_begins, Model, virus) %>% 
  summarise(across(ends_with("Hosp"), ~sum(.x, na.rm = T)), season = season[1]) %>% 
  mutate(model_loCI = qpois(0.025, ModelHosp, lower.tail = T), 
         model_hiCI = qpois(0.025, ModelHosp, lower.tail = F)) %>%  
  left_join(catchment_season) %>%
  mutate(across(c("DataHosp", "ModelHosp", "model_loCI", "model_hiCI"), ~.x/CatchmentPop*1e5)) %>%
  filter(virus == "Codetection") %>% 
  mutate(Model = c("Best model", "No interaction", "Next best model")[match(Model, c("best", "baseline", "nextbest"))]) %>% 
  mutate(Model = factor(Model, levels = c("Best model", "No interaction", "Next best model"))) %>% 
  ggplot(aes(x = week_begins)) + 
  # {if(include_ci) geom_ribbon(aes(ymin = model_loCI, ymax = model_hiCI, fill = "CI"), show.legend = FALSE)} + 
  
  geom_line(aes(y = ModelHosp,  colour = Model), linewidth = 1) + 
  geom_point(aes(y = DataHosp, shape = "Data"), colour = "red", size = 0.8) +
  scale_colour_manual(values = c(`Best model` = "blue", `No interaction` = "orange3", `Next best model` = "lightblue4")) + 
  # facet_grid(Model~., scales = "free_y") + 
  labs(x = "Week", y = "Co-detection cases per 100 000", colour = "", linetype = "", shape = "") + 
  theme_bw()+ 
  scale_fill_manual(values = c(CI = "lightblue")) + 
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0))
# plot_comp_obj
ggsave(plot_comp_obj, filename = "figtab/FigS6_data2goodmodels.png", width = 25, height = 12, units = "cm", dpi = 600)

