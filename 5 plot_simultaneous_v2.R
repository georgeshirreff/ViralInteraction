season_vec
original <- generate_ModelHosp_season(pars = best_params %>%
                            select(-Loglik) %>%
                            unlist, generator = m5_generator, dat_2virus = dat_2virus,
                          data_catchments_weights = data_catchments_weights)$hosp %>% 
  mutate(interaction = T, intLevel = "best", in_phase = F)

nointer <- generate_ModelHosp_season(pars = best_params %>% 
                                        replace("delta1", 1) %>% 
                                        select(-Loglik) %>%
                                        unlist, generator = m5_generator, dat_2virus = dat_2virus,
                                      data_catchments_weights = data_catchments_weights)$hosp %>% 
  mutate(interaction = F, intLevel = NA, in_phase = F)

interHi <- generate_ModelHosp_season(pars = best_params %>% 
                                       replace("delta1", 0) %>% 
                                       select(-Loglik) %>%
                                       unlist, generator = m5_generator, dat_2virus = dat_2virus,
                                     data_catchments_weights = data_catchments_weights)$hosp %>% 
  mutate(interaction = T, intLevel = "hi", in_phase = F)

interLo <- generate_ModelHosp_season(pars = best_params %>% 
                                        replace("delta1", 1.8) %>% 
                                        select(-Loglik) %>%
                                        unlist, generator = m5_generator, dat_2virus = dat_2virus,
                                      data_catchments_weights = data_catchments_weights)$hosp %>% 
  mutate(interaction = T, intLevel = "lo", in_phase = F)




# best_params["eta2_1213"]
# best_params["eta1_1213"]
phase <- generate_ModelHosp_season(pars = best_params %>%
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

phase_nointer <- generate_ModelHosp_season(pars = best_params %>% 
                                             replace("eta1_1213", -17.0) %>% 
                                             replace("delta1", 1) %>% 
                                             select(-Loglik) %>%
                                             unlist, generator = m5_generator, dat_2virus = dat_2virus,
                                           data_catchments_weights = data_catchments_weights)$hosp %>% 
  mutate(interaction = F, intLevel = NA, in_phase = T)


phase_interHi <- generate_ModelHosp_season(pars = best_params %>% 
                                             replace("eta1_1213", -17.0) %>% 
                                             replace("delta1", 0) %>% 
                                             select(-Loglik) %>%
                                             unlist, generator = m5_generator, dat_2virus = dat_2virus,
                                           data_catchments_weights = data_catchments_weights)$hosp %>% 
  mutate(interaction = T, intLevel = "hi", in_phase = T)

phase_interLo <- generate_ModelHosp_season(pars = best_params %>% 
                                             replace("eta1_1213", -17.0) %>% 
                                             replace("delta1", 1.8) %>% 
                                             select(-Loglik) %>%
                                             unlist, generator = m5_generator, dat_2virus = dat_2virus,
                                           data_catchments_weights = data_catchments_weights)$hosp %>% 
  mutate(interaction = T, intLevel = "lo", in_phase = T)

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

comp_table <- rbind(original, nointer, interLo, interHi, phase, phase_interLo, phase_interHi, phase_nointer) %>% 
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
  mutate(averted_best = `Inc_No interaction_NA` - `Inc_Interaction_best`, 
         averted_hi = `Inc_No interaction_NA` - `Inc_Interaction_hi`, 
         averted_lo = `Inc_No interaction_NA` - `Inc_Interaction_lo`) %>% 
  mutate(pAverted_best = averted_best/`Inc_No interaction_NA`, 
         pAverted_hi = averted_hi/`Inc_No interaction_NA`,
         pAverted_lo = averted_lo/`Inc_No interaction_NA`) %>% 
  relocate(starts_with("peak_time"), .after = everything()) %>% 
  mutate(shift_best = `peak_time_Interaction_best` - `peak_time_No interaction_NA`, 
         shift_hi = `peak_time_Interaction_hi` - `peak_time_No interaction_NA`, 
         shift_lo = `peak_time_Interaction_lo` - `peak_time_No interaction_NA`) %>% 
  mutate(averted_best = paste0(round(averted_best), " (", round(pAverted_best*100), "%)"), 
         averted_hi = paste0(round(averted_hi), " (", round(pAverted_hi*100), "%)"), 
         averted_lo = paste0(round(averted_lo), " (", round(pAverted_lo*100), "%)")) %>% 
  select(virus, in_phase, starts_with("averted"), starts_with("shift")) %>% 
  ungroup

comp_table_wide_long <-
  comp_table_wide %>% pivot_longer(-c("virus", "in_phase"), names_sep = "_", names_to = c(".value", "intLevel")) %>% 
  mutate(estimate = case_match(intLevel, "best" ~ "Best", 
                               "hi" ~ "High", 
                               "lo" ~ "Low")) %>% 
  mutate(estimate = factor(estimate, levels = c("Low", "Best", "High"))) %>% 
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
simult_plot <- rbind(original, nointer, interLo, interHi, phase, phase_interLo, phase_interHi, phase_nointer) %>% 
  filter(season == "12/13") %>% 
  mutate(intLevel = case_match(intLevel, "best" ~ "Best estimate", 
                               "lo" ~ "Low estimate", 
                               "hi" ~ "High estimate")) %>% 
  mutate(intLevel = factor(intLevel, levels = c("None", "Low estimate", "Best estimate", "High estimate"))) %>% 
  mutate(intLevel = replace_na(intLevel, "None")) %>% 
  group_by(virus, in_phase, intLevel, t) %>% 
  summarise(Inc = sum(Inc)) %>% 
  arrange(virus, in_phase, intLevel) %>% 
  mutate(virus =fct_recode(factor(virus, levels = c(1, 2, 12)), `Influenza only`="1", `RSV only`="2", `Codetection`="12")) %>% 
  # mutate(interaction = ifelse(interaction, "Interaction", "No interaction")) %>% 
  mutate(in_phase = ifelse(in_phase, "Simultaneous peaks", "Original")) %>% 
  ggplot(aes(x = t, y = Inc, colour = virus, linetype = intLevel)) + 
  geom_line() + 
  # scale_colour_manual(values = vir_cols)+ 
  scale_linetype_manual(values = c(None = "solid", `Best estimate` = "dashed", `High estimate` = "dotted", `Low estimate` = "dotdash")) + 
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


ggsave(filename = "figtab/Fig6_simultaneous_peaks.png", width = 20, height = 12, units = "cm", dpi = 600)




