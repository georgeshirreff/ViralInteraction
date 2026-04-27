#### FIGURE S5 do it but maybe cut it???

# # paper version v16
# best_params <- read_csv(paste0("output/", id = "validateSimpleInter_", virus1, "+", virus2, "_jump", 3, "_rep", 1, "_bestParams.csv"))
# baseline_params <- read_csv(paste0("output/", id = "simpleInter_", virus1, "+", virus2, "_jump", 0, "_rep", 0, "_bestParams.csv"))
# nextbest_params <- read_csv(paste0("output/", id = "simpleInter_", virus1, "+", virus2, "_jump",5, "_rep", 2, "_bestParams.csv"))

# paper version v17
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


tib_ModelHosp_baseline <- generate_ModelHosp_season(baseline_params %>% 
                                                      select(-Iteration, -Loglik) %>% 
                                                      unlist, m5_generator, dat_2virus, data_catchments_weights = data_catchments_weights)$hosp

tib_ModelHosp_best <- generate_ModelHosp_season(pars = best_params %>% 
                                                  select(-Iteration, -Loglik) %>% 
                                                  unlist, generator = m5_generator, dat_2virus = dat_2virus, data_catchments_weights = data_catchments_weights)$hosp

tib_ModelHosp_nextbest <- generate_ModelHosp_season(nextbest_params %>% 
                                                      select(-Iteration, -Loglik) %>% 
                                                      unlist, m5_generator, dat_2virus, data_catchments_weights = data_catchments_weights)$hosp

phase_tib_ModelHosp_baseline <- generate_ModelHosp_season(baseline_params %>% 
                                                            replace("eta1_1213", -17.0) %>% 
                                                            select(-Iteration, -Loglik) %>% 
                                                            unlist, m5_generator, dat_2virus, data_catchments_weights = data_catchments_weights)$hosp

# phase_tib_ModelHosp_baseline %>% 
#   filter(season == "12/13") %>% 
#   mutate(virus =fct_recode(factor(virus, levels = c(1, 2, 12)), `Influenza only`="1", `RSV only`="2", `Codetection`="12")) %>% 
#   group_by(t, virus) %>% 
#   summarise(Inc = sum(Inc)) %>% 
#   ggplot(aes(x = t, y = Inc, colour = virus)) + geom_line()

phase_tib_ModelHosp_best <- generate_ModelHosp_season(pars = best_params %>% 
                                                        replace("eta1_1213", -17.0) %>% 
                                                        select(-Iteration, -Loglik) %>% 
                                                        unlist, generator = m5_generator, dat_2virus = dat_2virus, data_catchments_weights = data_catchments_weights)$hosp

phase_tib_ModelHosp_nextbest <- generate_ModelHosp_season(nextbest_params %>% 
                                                            replace("eta1_1213", -17.0) %>% 
                                                            select(-Iteration, -Loglik) %>% 
                                                            unlist, m5_generator, dat_2virus, data_catchments_weights = data_catchments_weights)$hosp


nointer <- tib_ModelHosp_baseline %>% mutate(interaction = F, intLevel = NA, in_phase = F)
bestinter <- tib_ModelHosp_best %>% mutate(interaction = T, intLevel = "best", in_phase = F)
nextbestinter <- tib_ModelHosp_nextbest %>% mutate(interaction = T, intLevel = "nextbest", in_phase = F)

phase_nointer <- phase_tib_ModelHosp_baseline %>% mutate(interaction = F, intLevel = NA, in_phase = T)
phase_bestinter <- phase_tib_ModelHosp_best %>% mutate(interaction = T, intLevel = "best", in_phase = T)
phase_nextbestinter <- phase_tib_ModelHosp_nextbest %>% mutate(interaction = T, intLevel = "nextbest", in_phase = T)




# make Figure 6 (basic comparison of two bets models vs baseline)

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
  {if(F) geom_ribbon(aes(ymin = model_loCI, ymax = model_hiCI, fill = "CI"), show.legend = FALSE)} +
  
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





#####

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
  mutate(averted_best = `Inc_No interaction_NA` - `Inc_Interaction_best`, 
         averted_nextbest = `Inc_No interaction_NA` - `Inc_Interaction_nextbest`) %>% 
  mutate(pAverted_best = averted_best/`Inc_No interaction_NA`, 
         pAverted_nextbest = averted_nextbest/`Inc_No interaction_NA`) %>% 
  relocate(starts_with("peak_time"), .after = everything()) %>% 
  mutate(shift_best = `peak_time_Interaction_best` - `peak_time_No interaction_NA`, 
         shift_nextbest = `peak_time_Interaction_nextbest` - `peak_time_No interaction_NA`) %>% 
  mutate(averted_best = paste0(round(averted_best), " (", round(pAverted_best*100), "%)"), 
         averted_nextbest = paste0(round(averted_nextbest), " (", round(pAverted_nextbest*100), "%)")) %>% 
  select(virus, in_phase, starts_with("averted"), starts_with("shift")) %>% 
  ungroup

comp_table_wide_long <-
  comp_table_wide %>% pivot_longer(-c("virus", "in_phase"), names_sep = "_", names_to = c(".value", "intLevel")) %>% 
  mutate(estimate = case_match(intLevel, "best" ~ "Best (theta1)", 
                               "nextbest" ~ "2nd best (psi1)")) %>% 
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
  mutate(intLevel = case_match(intLevel, "best" ~ "Best (theta1)", 
                               "nextbest" ~ "2nd best (psi1)")) %>% 
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


