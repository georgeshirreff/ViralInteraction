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

baseline_params <- read_csv(paste0("output/", id = "simpleInter_", virus1, "+", virus2, "_jump", 0, "_rep", 0, "_bestParams.csv"))

best_post <- read_csv2(paste0("output/", id = "interOnly_", virus1, "+", virus2, "_jump", 5, "_rep", 0, "_posterior.csv")) %>% 
  select(-Accepted, -Season)
best_params = best_post %>% 
  arrange(-Loglik) %>% 
  slice(1)

nextbest_post <- read_csv2(paste0("output/", id = "interOnly_", virus1, "+", virus2, "_jump", 2, "_rep", 0, "_posterior.csv")) %>% 
  select(-Accepted, -Season)
nextbest_params = nextbest_post %>% 
  arrange(-Loglik) %>% 
  slice(1)


season_vec

nointer <- generate_ModelHosp_season(pars = baseline_params %>% 
                                       select(-Loglik) %>%
                                       unlist, generator = m5_generator, dat_2virus = dat_2virus,
                                     data_catchments_weights = data_catchments_weights)$hosp %>% 
  mutate(interaction = F, intLevel = NA, in_phase = F)

bestinter <- generate_ModelHosp_season(pars = best_params %>%
                            select(-Loglik) %>%
                            unlist, generator = m5_generator, dat_2virus = dat_2virus,
                          data_catchments_weights = data_catchments_weights)$hosp %>% 
  mutate(interaction = T, intLevel = "theta1", in_phase = F)


nextbestinter <- generate_ModelHosp_season(pars = nextbest_params %>% 
                                       select(-Loglik) %>%
                                       replace("theta1", 0) %>% 
                                       unlist, generator = m5_generator, dat_2virus = dat_2virus,
                                     data_catchments_weights = data_catchments_weights)$hosp %>% 
  mutate(interaction = T, intLevel = "psi1", in_phase = F)


# best_params["eta2_1213"]
# best_params["eta1_1213"]
phase_noninter <- generate_ModelHosp_season(pars = baseline_params %>%
                                     replace("eta1_1213", -17.0) %>% 
                                     
                                        select(-Loglik) %>%
                                        unlist, generator = m5_generator, dat_2virus = dat_2virus,
                                      data_catchments_weights = data_catchments_weights)$hosp %>% 
  mutate(interaction = T, intLevel = "best", in_phase = T)

phase %>% 
  filter(season == "12/13") %>% 
  mutate(virus =fct_recode(factor(virus, levels = c(1, 2, 12)), `Influenza only`="1", `RSV only`="2", `Codetection`="12")) %>% 
  group_by(t, virus) %>% 
  summarise(Inc = sum(Inc)) %>% 
  ggplot(aes(x = t, y = Inc, colour = virus)) + geom_line()

phase_bestinter <- generate_ModelHosp_season(pars = best_params %>%
                                               replace("eta1_1213", -17.0) %>% 
                                               select(-Loglik) %>%
                                         unlist, generator = m5_generator, dat_2virus = dat_2virus,
                                       data_catchments_weights = data_catchments_weights)$hosp %>% 
  mutate(interaction = T, intLevel = "theta1", in_phase = T)


phase_nextbestinter <- generate_ModelHosp_season(pars = nextbest_params %>% 
                                             replace("eta1_1213", -17.0) %>% 
                                             select(-Loglik) %>%
                                             unlist, generator = m5_generator, dat_2virus = dat_2virus,
                                           data_catchments_weights = data_catchments_weights)$hosp %>% 
  mutate(interaction = T, intLevel = "psi1", in_phase = T)



# comp_table <- rbind(original, nointer, phase, phase_nointer) %>% 
#   filter(season == "12/13") %>% 
#   group_by(virus, in_phase, interaction, t) %>% 
#   summarise(Inc = sum(Inc)) %>% 
#   group_by(virus, in_phase, interaction) %>% 
#   summarise(peak_time = t[Inc == max(Inc)], Inc = sum(Inc)) %>% 
#   arrange(virus, in_phase) %>% 
#   mutate(virus =fct_recode(factor(virus, levels = c(1, 2, 12)), `Influenza only`="1", `RSV only`="2", `Codetection`="12")) %>% 
#   mutate(interaction = ifelse(interaction, "Interaction", "No interaction")) %>% 
#   mutate(in_phase = ifelse(in_phase, "Simultaneous peaks", "Original"))

comp_table <- rbind(nointer, bestinter, nextbestinter, phase_nointer, phase_bestinter, phase_nextbestinter) %>% 
  filter(season == "12/13") %>% 
  group_by(virus, in_phase, interaction, intLevel, t) %>% 
  summarise(Inc = sum(Inc)) %>% 
  group_by(virus, in_phase, interaction, intLevel) %>% 
  summarise(peak_time = t[Inc == max(Inc)], Inc = sum(Inc)) %>% 
  arrange(virus, in_phase) %>% 
  mutate(virus =fct_recode(factor(virus, levels = c(1, 2, 12)), `Influenza only`="1", `RSV only`="2", `Codetection`="12")) %>% 
  mutate(interaction = ifelse(interaction, "Interaction", "No interaction")) %>% 
  mutate(in_phase = ifelse(in_phase, "Simultaneous peaks", "Original"))


comp_table_wide <- comp_table %>% 
  pivot_wider(id_cols = c("virus", "in_phase"), names_from = c("interaction", "intLevel"), values_from = c("Inc", "peak_time")) %>% 
  mutate(averted_theta1 = `Inc_No interaction_NA` - `Inc_Interaction_theta1`, 
         averted_psi1 = `Inc_No interaction_NA` - `Inc_Interaction_psi1`) %>% 
  mutate(pAverted_theta1 = averted_theta1/`Inc_No interaction_NA`, 
         pAverted_psi1 = averted_psi1/`Inc_No interaction_NA`) %>% 
  relocate(starts_with("peak_time"), .after = everything()) %>% 
  mutate(shift_theta1 = `peak_time_Interaction_theta1` - `peak_time_No interaction_NA`, 
         shift_psi1 = `peak_time_Interaction_psi1` - `peak_time_No interaction_NA`) %>% 
  mutate(averted_theta1 = paste0(round(averted_theta1), " (", round(pAverted_theta1*100), "%)"), 
         averted_psi1 = paste0(round(averted_psi1), " (", round(pAverted_psi1*100), "%)")) %>% 
  select(virus, in_phase, starts_with("averted"), starts_with("shift")) %>% 
  ungroup

comp_table_wide_long <-
  comp_table_wide %>% pivot_longer(-c("virus", "in_phase"), names_sep = "_", names_to = c(".value", "intLevel")) %>% 
  mutate(estimate = case_match(intLevel, "theta1" ~ "Best (theta1)", 
                               "psi1" ~ "2nd best (psi1)")) %>% 
  mutate(estimate = factor(estimate, levels = c("Best (theta1)", "2nd best (psi1)"))) %>% 
  arrange(virus, in_phase, estimate) %>% 
  group_by(virus, in_phase) %>% 
  summarise(across(c("estimate", "averted", "shift"), ~paste0(.x, collapse = "\n")))
  


tbs <- lapply(split(comp_table_wide_long, comp_table_wide_long$in_phase), "[", -2)
# tbs <- lapply(split(comp_table_wide, ~in_phase+virus), "[", -2)
df <- tibble(x   = rep(Inf, length(tbs)), 
             y   = rep(Inf, length(tbs)), 
             in_phase = levels(as.factor(comp_table_wide_long$in_phase)), 
             tbl = tbs)



write_csv2(comp_table_wide_long, file = "figtab/Tab3_simultaneouspeaks_lobesthi.csv")

vir_cols = c("Influenza only"="red", "RSV only" = "blue", Codetection ="violet")
simult_plot <- rbind(nointer, bestinter, nextbestinter, phase_nointer, phase_bestinter, phase_nextbestinter) %>% 
  filter(season == "12/13") %>% 
  mutate(intLevel = case_match(intLevel, "theta1" ~ "Best (theta1)", 
                               "psi1" ~ "2nd best (psi1)")) %>% 
  mutate(intLevel = factor(intLevel, levels = c("None", "Best (theta1)", "2nd best (psi1)"))) %>% 
  mutate(intLevel = replace_na(intLevel, "None")) %>% 
  # filter(!intLevel == "None") %>% 
  group_by(virus, in_phase, intLevel, t) %>% 
  summarise(Inc = sum(Inc)) %>% 
  arrange(virus, in_phase, intLevel) %>% 
  mutate(virus =fct_recode(factor(virus, levels = c(1, 2, 12)), `Influenza only`="1", `RSV only`="2", `Codetection`="12")) %>% 
  # mutate(interaction = ifelse(interaction, "Interaction", "No interaction")) %>% 
  mutate(in_phase = ifelse(in_phase, "Simultaneous peaks", "Original")) %>% 
  ggplot(aes(x = t, y = Inc, colour = virus, linetype = intLevel)) + 
  geom_line() + 
  # scale_colour_manual(values = vir_cols)+ 
  scale_linetype_manual(values = c(None = "solid", `Best (theta1)` = "dashed", `2nd best (psi1)` = "dotted", `Low estimate` = "dotdash")) + 
  facet_grid(. ~ in_phase) + 
  theme_bw() + 
  labs(x = "Weeks relative to start of year", 
       y = "Cases per 100 000", 
       colour = "Virus", linetype = "Interaction") + 
  coord_cartesian(ylim = c(0, 75))

simult_plot + 
  ggpp::geom_table(data = df, aes(x = x, y = y, label = tbs),
                   inherit.aes = FALSE,
                   table.colnames = T, table.hjust = 1,
                   table.theme = gridExtra::ttheme_default(base_size = 5, 
                     core = list(
                       padding=unit(c(3, 2), "mm"), 
                       core.just = "right", 
                       # fg_params = list(cex = 0.5, lineheight = 0.6, col = "black"),  # shrink body rows
                       # fg_params = list(col = "darkgrey"),
                       bg_params = list(cex = 0.5, fill = "white", col = "darkgrey")
                     )
                   ))


ggsave(filename = "figtab/Fig5_simultaneous_peaks.png", width = 20, height = 12, units = "cm", dpi = 600)




