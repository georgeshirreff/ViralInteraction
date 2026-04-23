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

# compare_reps <- function(stem = "output/", id = id, virus1 = virus1, virus2 = virus2, 
#                          jump = 0, reps = 0:5, burnin = 0, thin = 1, 
#                          exclude_params = c(paste0(rep(c("alpha", "gamma", "omega", "psi", "delta", "sigma", "theta"), each = 2), 1:2), paste0("pi_", age_vec))
#                          # exclude_params = c(paste0(rep(c("alpha", "gamma", "psi", "delta", "sigma", "theta"), each = 2), 1:2), paste0("pi_", age_vec))
#                          # exclude_params = c("alpha1", "gamma1", "alpha2", "gamma2", paste0("pi_", age_vec))
# ){

#' @UPDATE-FIGURE6 

  library(ggstance)
  library(ggpubr) 
  
  stem = "output/"
  
  # these values from v16 of the paper
  # id = "validateSimpleInter_"
  # post = read_csv(paste0(stem, id, virus1, "+", virus2, "_jump", 3, "_rep", 1, 
  #                            "_posterior.csv"))
  
  # these values for v17
  id = "simpleInter_"
  post = read_csv(paste0(stem, id, virus1, "+", virus2, "_jump", 0, "_rep", 0, 
                         "_posterior.csv"))
  
  burnin = 0  
  
  if(burnin > max(post$Iteration)){stop("burnin too long")}
  
  post_cropped <- post %>% 
    filter(Iteration >= burnin) %>% 
    filter((Iteration %% thin) == 0) %>% 
    select(-ends_with("ccepted"), -ends_with("eason")) %>% 
    select_if(~n_distinct(.x) > 1)
  
  post_long <- post_cropped %>% 
    pivot_longer(-c("Iteration", "Loglik"), 
                 names_pattern = "([a-z]+)([12])([_]?.*)", 
                 names_to = c("Parameter", "virus", "category")) %>% 
    mutate(category = gsub("_", "", category)) %>% 
    mutate(season = ifelse(grepl("^[0-9]{4}$", category), category, "common"), 
           age_group = ifelse((season == "common") & (!(category == "")), category, "all")) %>% 
    select(-category)
  
  
  K <- post_long %>% 
    select(Parameter, virus, age_group, season) %>% 
    unique %>% nrow
  
  post_ci <- post_long %>% 
    filter(!is.na(Loglik)) %>% 
    group_by(Parameter, virus, age_group, season) %>% 
    mutate(maxLoglik = max(Loglik, na.rm = T), 
           deltaLoglik = Loglik - maxLoglik, 
           LRTstat = -2*deltaLoglik, 
           LRTp = pchisq(q = LRTstat, df = K, lower.tail = F), 
           inCI = (0.05 <= LRTp)) %>% 
    summarise(best = median(value[Loglik == maxLoglik]), 
              loCI = min(value[inCI]), 
              hiCI = max(value[inCI]))
  
  custom_labels <- c(
    kappa = "kappa ~ '(Probability of community ILI)'",
    beta  = "beta ~ '(Transmission rate)'", 
    h  = "h ~ '(Probability of hospitalisation)'",
    eta  = "eta ~ '(Introduction date)'"
  )
  
  plot_kappa <- post_ci %>% 
    filter(Parameter == "kappa") %>% 
    # mutate(Parameter = "ILI rate") %>% 
    # mutate(Parameter = expression(kappa)) %>%
    mutate(virus = c("Influenza", "RSV")[as.numeric(virus)]) %>% 
    mutate(age_group = factor(age_group, levels = age2_vec)) %>% 
    ggplot(aes(y = age_group, colour = virus)) + 
    geom_point(aes(x = best), position = position_dodgev(height = 0.5)) + 
    geom_errorbarh(aes(xmin = loCI, xmax = hiCI), position = position_dodgev(height = 0.5)) + 
    facet_grid(~Parameter, scales = "free_x", labeller = labeller(Parameter = as_labeller(custom_labels, label_parsed))) +
    # facet_grid(~Parameter, scales = "free_x", labeller = label_parsed) +
    theme_bw() + 
    scale_x_continuous(labels = scales::percent) + 
    labs(x = "Proportion of incident infections", y = "Age group (b)")
  
  plot_lh <- post_ci %>% 
    filter(Parameter == "lh") %>% 
    mutate(Parameter = "h") %>% 
    mutate(virus = c("Influenza", "RSV")[as.numeric(virus)]) %>% 
    mutate(age_group = factor(age_group, levels = age_vec)) %>% 
    ggplot(aes(y = age_group, colour = virus)) + 
    geom_point(aes(x = 10^(best)), position = position_dodgev(height = 0.5)) + 
    geom_errorbarh(aes(xmin = 10^(loCI), xmax = 10^(hiCI)), position = position_dodgev(height = 0.5)) + 
    facet_grid(~Parameter, scales = "free_x", labeller = labeller(Parameter = as_labeller(custom_labels, label_parsed))) +
    # facet_grid(.~Parameter, scales = "free_x") +
    theme_bw() + 
    scale_x_continuous(labels = function(x) scales::percent(x, drop0trailing=TRUE), trans = "log10") + 
    labs(x = "Proportion of incident infections", y = "Age group (a)")
  
  
  
  plot_beta <- post_ci %>% 
    filter(Parameter %in% c("beta")) %>% 
    filter(age_group == "all") %>% 
    mutate(season = gsub("^([0-9][0-9])", "20\\1/", season)) %>%
    mutate(virus = c("Influenza", "RSV")[as.numeric(virus)]) %>% 
    ggplot(aes(y = season, colour = virus)) + 
    geom_point(aes(x = best), position = position_dodgev(height = 0.5)) + 
    geom_errorbarh(aes(xmin = loCI, xmax = hiCI), position = position_dodgev(height = 0.5)) + 
    facet_grid(~Parameter, scales = "free_x", labeller = labeller(Parameter = as_labeller(custom_labels, label_parsed))) +
    # facet_grid(.~Parameter, scales = "free_x", labeller = label_parsed)  +
    theme_bw() + 
    labs(x = "Transmissions per potentially infectious contact per week", y = "Season")

  Sys.setlocale("LC_ALL", "EN")
  plot_eta <- post_ci %>% 
    filter(Parameter %in% c("eta")) %>% 
    filter(age_group == "all") %>% 
    mutate(season = gsub("^([0-9][0-9])", "20\\1/", season)) %>%
    mutate(virus = c("Influenza", "RSV")[as.numeric(virus)]) %>% 
    mutate(across(c("best", "loCI", "hiCI"), ~as.Date("1969-12-31") + .x*7)) %>%
    ggplot(aes(y = season, colour = virus)) + 
    geom_point(aes(x = best), position = position_dodgev(height = 0.5)) + 
    geom_errorbarh(aes(xmin = loCI, xmax = hiCI), position = position_dodgev(height = 0.5)) + 
    facet_grid(~Parameter, scales = "free_x", labeller = labeller(Parameter = as_labeller(custom_labels, label_parsed))) +
    # facet_grid(.~Parameter, scales = "free_x", labeller = label_parsed)  +
    theme_bw() + 
    labs(x = "Week of first transmission", y = "Season")
  # 
  # plot_interaction <- post_ci %>% 
  #   filter(age_group == "all" & season == "common") %>% 
  #   ggplot(aes(y = Parameter, group = rep)) + 
  #   geom_point(aes(x = best), position = position_dodgev(height = 0.5)) + 
  #   geom_errorbarh(aes(xmin = loCI, xmax = hiCI), position = position_dodgev(height = 0.5)) + 
  #   facet_grid(virus~., scales = "free_x")  +
  #   theme_bw()

  plot_reps <- ggarrange(ggarrange(plot_beta, plot_eta, nrow = 1, common.legend = T, labels = c("a", "b"), align = "h"), 
                         ggarrange(plot_lh, plot_kappa, nrow = 1, legend = F, labels = c("c", "d")), 
                         nrow = 2, common.legend = T)
  
  # plot_reps <- ggarrange(plot_lh, plot_kappa, plot_betaeta,  nrow = 3, common.legend = T)
  
  plot_reps %>% ggsave(filename = "figtab/Fig2_noninteractionParams.png", width = 25, height = 25, units = "cm", dpi = 600)
  # plot_reps %>% ggsave(filename = "~/SanOdin/output/nonint_params.png", width = 15, height = 15, units = "cm")
                         
  