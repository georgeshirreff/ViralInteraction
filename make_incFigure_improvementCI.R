best_post <- read_csv2(paste0("output/", id = "interOnly_", virus1, "+", virus2, "_jump", 5, "_rep", 0, "_posterior.csv")) %>% 
  select(-Accepted, -Season)
baseline_loglik <- read_csv(paste0("output/", id = "simpleInter_", virus1, "+", virus2, "_jump", 0, "_rep", 0, "_posterior.csv")) %>% 
  pull(Loglik) %>% 
  max(na.rm = T)
baseline_aic = 2*K - 2*baseline_loglik

best_params = best_post %>% 
  arrange(-Loglik) %>% 
  slice(1)


improvement_post = best_post %>% 
  mutate(AIC = 2*67 - 2*Loglik) %>% 
  filter(AIC < baseline_aic)

lo_params = improvement_post %>% 
  arrange(theta1) %>% 
  slice(1) 

hi_params = improvement_post %>% 
  arrange(-theta1) %>% 
  slice(1) 






model_best <- generate_ModelHosp_season(pars = best_params %>% 
                                         select(-Iteration, -Loglik) %>% 
                                         unlist, generator = m5_generator, dat_2virus = dat_2virus, data_catchments_weights = data_catchments_weights)

model_lo <- generate_ModelHosp_season(pars = lo_params %>% 
                                          select(-Iteration, -Loglik) %>% 
                                          unlist, generator = m5_generator, dat_2virus = dat_2virus, data_catchments_weights = data_catchments_weights)

model_hi <- generate_ModelHosp_season(pars = hi_params %>% 
                                        select(-Iteration, -Loglik) %>% 
                                        unlist, generator = m5_generator, dat_2virus = dat_2virus, data_catchments_weights = data_catchments_weights)

remove(tib_ModelHosp)
# tib_ModelHosp <- model_best$hosp

plot_data <- dat_2virus %>% 
  transmute(season, age_group = factor(age_group, levels = age_vec), t, 
            virus = c(virus1, virus2, "Codetection")[match(virus, c(1, 2, 12))], 
            DataHosp = Inc, 
            ModelHosp = model_best$hosp$Inc, 
            loModelHosp = model_lo$hosp$Inc, 
            hiModelHosp = model_hi$hosp$Inc
  )
# %>% 
#   mutate(DataHosp = DataHosp/CatchmentPop*1e5, 
#          ModelHosp = ModelHosp/CatchmentPop*1e5)

plot_data_date <- plot_data %>% 
  mutate(W1 = ISOweek::ISOweek2date(paste0("20", gsub(".*/", "", season), "-W01-1"))) %>% 
  mutate(week_begins = W1 - 7 + t*7) %>% 
  group_by(age_group, virus, season, week_begins) %>% 
  summarise(across(ends_with("Hosp"), ~sum(.x, na.rm = T))) %>% 
  ungroup %>% 
  complete(age_group, virus, season, week_begins, fill = list(DataHosp = NA, ModelHosp = NA, loModelHosp = NA, hiModelHosp = NA))


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
  # mutate(model_loCI = qpois(0.025, ModelHosp, lower.tail = T), 
  #        model_hiCI = qpois(0.025, ModelHosp, lower.tail = F)) %>%  
  left_join(catchment_season) %>%
  mutate(across(c("DataHosp", "ModelHosp", "loModelHosp", "hiModelHosp"), ~.x/CatchmentPop*1e5)) %>%
  ggplot(aes(x = week_begins)) + 
  {if(T) geom_ribbon(aes(ymin = loModelHosp, ymax = hiModelHosp, fill = "CI"), show.legend = FALSE)} + 
  
  geom_line(aes(y = ModelHosp,  linetype = "Best model"), colour = "blue") + 
  geom_point(aes(y = DataHosp, shape = "Data"), colour = "red", size = 0.8) +
  facet_grid(virus~., scales = "free_y") + 
  labs(x = "Week", y = "Cases per 100 000", colour = "", linetype = "", shape = "") + 
  theme_bw()+ 
  scale_fill_manual(values = c(CI = "lightblue")) + 
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0))


