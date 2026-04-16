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

# plot size of difference

#' @UPDATE-FIGURE5 
id = "validateSimpleInter_"
best_params <- read_csv(paste0("output/", id, virus1, "+", virus2, "_jump", 3, "_rep", 1, "_bestParams.csv"))
post <- read_csv(paste0("output/", id, virus1, "+", virus2, "_jump", 3, "_rep", 1, "_posterior.csv"))


plot_model_data(pars = best_params %>%
                  select(-Loglik) %>%
                  unlist,
                generator = m5_generator, virus1 = virus1, virus2 = virus2, dat_2virus = dat_2virus, ili = ili, data_catchments_weights = data_catchments_weights, 
                # type = "ageseason"
                type = "age"
                # type = "season"
                # type = "ili_ageseason"
                # type = "ili_age"
                # type = "ili_season"
)

best_ModelHosp <- generate_ModelHosp_season(pars = best_params %>%
                                              select(-Loglik) %>%
                                              unlist, generator = m5_generator, dat_2virus = dat_2virus,
                                            data_catchments_weights = data_catchments_weights
)$hosp


best_ModelHosp <- generate_ModelHosp_season(pars = best_params %>%
                                              select(-Loglik) %>%
                                              unlist, generator = m5_generator, dat_2virus = dat_2virus,
                                            data_catchments_weights = data_catchments_weights
)$hosp %>% 
  group_by(season, t, virus) %>% 
  summarise(Inc = sum(Inc))

# make a subsample
N = 1000
set.seed(1)
sample_params = post[sample(x = nrow(post), size = N), ] %>% 
  select(-Season)

sample_ModelHosp <- map(split(sample_params, 1:nrow(sample_params)), 
  # split(sample_params[1:3, ], 1:3), 
    .f = function(p) generate_ModelHosp_season(pars = unlist(p), generator = m5_generator, dat_2virus = dat_2virus,
                                               data_catchments_weights = data_catchments_weights)$hosp %>% 
      group_by(season, t, virus) %>% 
      summarise(Inc = sum(Inc)), 
    .progress = list(
        type = "iterator", 
        format = "Calculating {cli::pb_bar} {cli::pb_percent}",
        clear = TRUE)) %>% 
  data.table::rbindlist(idcol = "rep") %>% 
  as_tibble


sample_ModelHosp_LoHi <- sample_ModelHosp %>% 
  group_by(season, t, virus) %>% 
  summarise(LoInc = quantile(Inc, 0.025), HiInc = quantile(Inc, 0.975))



nointer_ModelHosp <- generate_ModelHosp_season(pars = best_params %>%
                                                 mutate(across(starts_with(c("lromega", "lrpsi")), ~-3), 
                                                        across(starts_with(c("delta", "sigma", "theta")), ~1)) %>% 
                                                 select(-Loglik) %>%
                                                 unlist, generator = m5_generator, dat_2virus = dat_2virus,
                                               data_catchments_weights = data_catchments_weights
)$hosp %>% 
  group_by(season, t, virus) %>% 
  summarise(Inc = sum(Inc)) %>% 
  rename(BaselineInc = Inc)

datainc = dat_2virus %>% 
  group_by(season, t, virus) %>% 
  summarise(DataInc = sum(Inc, na.rm = T))

incdiff <- best_ModelHosp %>% 
  left_join(nointer_ModelHosp) %>% 
  left_join(datainc) %>% 
  left_join(catchment_season) %>%
  left_join(sample_ModelHosp_LoHi) %>% 
  # left_join(dat_2virus %>% 
  #             filter(t == 0, virus == 1) %>% 
  #             group_by(season) %>% 
  #             summarise(CatchmentPop = sum(CatchmentPop))) %>% 
  mutate(across(ends_with("Inc"), ~.x/CatchmentPop*1e5)) %>%
  mutate(IncDiff = (Inc - BaselineInc)) %>% 
  mutate(virus =fct_recode(factor(virus, levels = c(1, 2, 12)), `Influenza only`="1", `RSV only`="2", `Codetection`="12"))

# calculate the proportion of the data points fall inside the 95% credibility interval
incdiff %>% 
  ungroup %>% 
  mutate(inblue = ((DataInc < LoInc) | (DataInc > HiInc))) %>% 
  group_by(virus) %>% 
  summarise(codetection_inblue = mean(inblue))

best_ModelHosp %>% 
  filter(virus == 1, season == "11/12") %>% 
  pull(Inc) %>% sum
nointer_ModelHosp %>% 
  filter(virus == 1, season == "11/12") %>% 
  pull(BaselineInc) %>% sum

datainc %>% 
  filter(virus == 1, season == "11/12") %>% 
  pull(DataInc) %>% sum

incdiff %>% 
  group_by(season, virus) %>% 
  summarise(IncDiff = sum(IncDiff), Inc = sum(Inc), BaselineInc = sum(BaselineInc)) %>% 
  group_by(virus) %>% 
  summarise(IncDiff = mean(IncDiff), Inc = mean(Inc), BaselineInc = mean(BaselineInc)) %>% 
  mutate(relIncDiff = IncDiff/BaselineInc)

incdiff %>% 
  # left_join(catchment_season) %>% 
  group_by(virus, season) %>% 
  mutate(CumInc = cumsum(Inc), 
         CumBaselineInc = cumsum(BaselineInc), 
         CumDataInc = cumsum(DataInc), 
         CumLoInc = cumsum(LoInc), 
         CumHiInc = cumsum(HiInc)) %>% 
  # mutate(CumInc = cumsum(Inc/CatchmentPop*1e5), CumBaselineInc = cumsum(BaselineInc/CatchmentPop*1e5)) %>% 
  ggplot(aes(x = t)) + 
  geom_ribbon(aes(ymin = CumLoInc, ymax = CumHiInc, fill = "Best model")) + 
  geom_line(aes(y = CumInc, colour = "Best model"), linewidth = 1) + 
  # geom_line(aes(y = CumLoInc, colour = "Best model"), alpha = 0.5) +
  # geom_line(aes(y = CumHiInc, colour = "Best model"), alpha = 0.5) +
  geom_line(aes(y = CumBaselineInc, colour = "No interaction"), linewidth = 1) + 
  geom_line(aes(y = CumDataInc, colour = "Data"), linewidth = 1) +
  facet_grid(virus~season, scales = "free_y") + 
  labs(x = "Week relative to start of year", y = "Cumulative cases per 100 000", colour = "", fill = "Credibility interval") + 
  scale_colour_manual(values = c(`Best model` = "blue", `No interaction`="orange3", Data = "red")) + 
  scale_fill_manual(values = c(`Best model` = "lightblue")) + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 6))

ggsave(filename = "figtab/Fig5_cumulativeCases_CI.png", width = 25, height = 15, units = "cm", dpi = 600)

incdiff %>% 
  # left_join(catchment_season) %>% 
  group_by(virus, season) %>% 
  mutate(CumInc = cumsum(Inc), 
         CumBaselineInc = cumsum(BaselineInc), 
         CumDataInc = cumsum(DataInc), 
         CumLoInc = cumsum(LoInc), 
         CumHiInc = cumsum(HiInc)) %>% 
  # mutate(CumInc = cumsum(Inc/CatchmentPop*1e5), CumBaselineInc = cumsum(BaselineInc/CatchmentPop*1e5)) %>% 
  ggplot(aes(x = t)) + 
  geom_line(aes(y = CumInc, colour = "Best model")) + 
  geom_line(aes(y = CumBaselineInc, colour = "No interaction")) + 
  geom_line(aes(y = CumDataInc, colour = "Data")) +
  facet_grid(virus~season, scales = "free_y") + 
  labs(x = "Week relative to start of year", y = "Cumulative cases per 100 000", colour = "", fill = "Credibility interval") + 
  scale_colour_manual(values = c(`Best model` = "blue", `No interaction`="orange3", Data = "red")) + 
  # scale_fill_manual(values = c(`Best model` = "blue")) + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 6))

ggsave(filename = "figtab/Fig5_cumulativeCases.png", width = 25, height = 15, units = "cm", dpi = 600)


incdiff %>% 
  ggplot(aes(x = t, y = IncDiff)) + geom_line() + 
  facet_grid(virus~season, scales = "free_y") + 
  labs(x = "Week", y = "Increase in cases due to interaction") + 
  theme_bw()

ggsave(filename = "~/SanOdin/output/averted_plot.png", width = 30, height = 16, units = "cm")
