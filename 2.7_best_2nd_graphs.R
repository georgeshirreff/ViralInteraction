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

