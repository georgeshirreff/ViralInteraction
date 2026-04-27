

#### Figure 3 ####

# from paper v16
# best_params <- read_csv(paste0("output/", id = "validateSimpleInter_", virus1, "+", virus2, "_jump", 3, "_rep", 1, "_bestParams.csv"))

# for paper v17
# best_params <- read_csv(paste0("output/", id = "interOnly_", virus1, "+", virus2, "_jump", 5, "_rep", 0, "_bestParams.csv"))
best_post <- read_csv2(paste0("output/", id = "interOnly_", virus1, "+", virus2, "_jump", 5, "_rep", 0, "_posterior.csv")) %>% 
  select(-Accepted, -Season)
best_params = best_post %>% 
  arrange(-Loglik) %>% 
  slice(1)

model_out <- generate_ModelHosp_season(pars = best_params %>% 
                                         select(-Iteration, -Loglik) %>% 
                                         unlist, generator = m5_generator, dat_2virus = dat_2virus, data_catchments_weights = data_catchments_weights)

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
  {if(T) geom_ribbon(aes(ymin = model_loCI, ymax = model_hiCI, fill = "CI"), show.legend = FALSE)} + 
  
  geom_line(aes(y = ModelHosp,  linetype = "Best model"), colour = "blue") + 
  geom_point(aes(y = DataHosp, shape = "Data"), colour = "red", size = 0.8) +
  facet_grid(virus~., scales = "free_y") + 
  labs(x = "Week", y = "Cases per 100 000", colour = "", linetype = "", shape = "") + 
  theme_bw()+ 
  scale_fill_manual(values = c(CI = "lightblue")) + 
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0))

# ggsave(plot_obj, filename = "figtab/Fig4.png", width = 25, height = 15, units = "cm", dpi = 600)
# ggsave(plot_obj, filename = "figtab/Fig4_dataBestmodel.png", width = 25, height = 15, units = "cm", dpi = 600)
ggsave(plot_obj, filename = "figtab/Fig3_dataBestmodel.png", width = 25, height = 15, units = "cm", dpi = 600)

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


