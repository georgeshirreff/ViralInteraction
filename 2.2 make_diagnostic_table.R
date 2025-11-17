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

library(coda)
library(data.table)
library(tidyverse)

param_jumps = read_csv(
  # file = paste0("input/odin2v_jumps.csv")
  # file = paste0("input/odin2v_jumps_etaonly.csv")
  # file = paste0("input/odin2v_jumps_etabeta.csv")
  # file = paste0("input/odin2v_jumps_etabetakappa.csv")
  # file = paste0("input/odin2v_jumps_etabetakappalh.csv")
  # file = paste0("input/odin2v_jumps_", "all", ".csv")
  file = paste0("input/lrodin2v_jumps_", "all", ".csv")
)
# this_jump = 3
# this_rep = 0:5
param_table <- function(stem = "output/", id, virus1, virus2, 
             this_jump = 0, this_rep = 0, burnin = 0, thin = 1, param_jumps){

  if(length(this_rep) == 1){
    post = read_csv(paste0(stem, id, virus1, "+", virus2, 
                           ifelse(is.null(this_jump), "", paste0("_jump", this_jump)), 
                           ifelse(is.null(this_rep), "", paste0("_rep", this_rep)), 
                           "_posterior.csv"), show_col_types = FALSE)
  } else {
    post = NULL
    for(r in this_rep){
      cat(r)
      post_piece = read_csv(paste0(stem, id, virus1, "+", virus2, 
                                   ifelse(is.null(this_jump), "", paste0("_jump", this_jump)), 
                                   ifelse(is.null(this_rep), "", paste0("_rep", r)), 
                                   "_posterior.csv"), show_col_types = FALSE)
      
      if(is.null(post)){
        post = post_piece
      } else {
        post = rbind(post, post_piece)
      }
    }
  }

  if(burnin > max(post$Iteration)){stop("burnin too long")}
  
  which_params = param_jumps %>% 
    filter(jump == this_jump) %>% 
    select(-jump) %>% unlist %>% 
    {.[. > 0]} %>% 
    names
  
  which_interaction = param_jumps %>% 
    filter(jump == this_jump) %>% 
    select(-starts_with("beta"), -starts_with("eta"), -starts_with("alpha"), -starts_with("gamma"), -starts_with("pi"), -starts_with("lh"), -starts_with("kappa")) %>% 
    select(-jump) %>% unlist %>% 
    {.[. > 0]} %>% 
    names
  

  post_cropped <- post %>% 
    filter(Iteration >= burnin) %>% 
    filter((Iteration %% thin) == 0) %>% 
    select(-ends_with("ccepted"), -ends_with("eason"))
  
  # print("reciprocal of process parameters")
  # post_cropped <- post_cropped %>% 
  #   mutate(across(starts_with(c("omega", "psi")), ~1/.x))
  
  post_long <- post_cropped %>% 
    select(Iteration, which_params, Loglik) %>% 
    pivot_longer(-c("Iteration", "Loglik"), 
                 names_pattern = "([a-z]+)([12])([_]?.*)", 
                 names_to = c("Parameter", "virus", "category")) %>% 
    mutate(category = gsub("_", "", category)) %>% 
    mutate(season = ifelse(grepl("^[0-9]{4}$", category), category, "common"), 
           age_group = ifelse((season == "common") & (!(category == "")), category, "all")) %>% 
    select(-category)
  
  if(F){
    post_cropped %>% 
      select(which_interaction, Loglik) %>% 
      ggplot(aes(x = get(which_interaction), y = Loglik)) + geom_point()
    
    post_cropped %>% 
      select(Iteration, which_interaction) %>% 
      ggplot(aes(x = Iteration, y = get(which_interaction))) + geom_point()
    
    
  }
  
  
  K = length(which_params)
  
  maxLoglik = post_cropped$Loglik %>% max(na.rm = T)
  
  post_return = post_long %>% 
    filter(!is.na(Loglik)) %>% 
    filter(Parameter != "alpha", Parameter != "gamma", Parameter != "pi") %>% 
    group_by(Parameter, virus, season, age_group) %>% 
    mutate(              deltaLoglik = Loglik - maxLoglik, 
                         LRTstat = -2*deltaLoglik, 
                         LRTp = pchisq(q = LRTstat, df = K, lower.tail = F), 
                         inCI = (0.05 <= LRTp)) %>% 
    summarise(best = last(value[Loglik == max(Loglik, na.rm = T)]),
              loCI = min(na.rm = T, value[inCI]), 
              hiCI = max(na.rm = T, value[inCI]), 
              
              # best = median(value),
              # loCI = quantile(value, 0.025), 
              # hiCI = quantile(value, 0.975),
              
              effSize = effectiveSize(value), 
              K = K, 
              maxLoglik = maxLoglik, 
              AIC = 2*K - 2*maxLoglik
    ) %>% 
    mutate(jump = this_jump, rep = -1) %>% 
    select(jump, rep, everything())
  
  post_return
    
  
}


# View(classif)
aic_table = NULL
aic_comp = NULL
for(j in 0:10){
  print(j)
  r_vec = 0:5
  
  if(j == 0){
    this_id = "simpleInter_"
  } else {
    this_id = "validateSimpleInter_"
  }
  this_burnin = 0

  aic_piece = param_table(stem = "output/", id = this_id, virus1 = "Influenza", virus2 = "RSV", 
              this_jump = j, this_rep = r_vec, burnin = this_burnin, thin = 1, param_jumps)
  
  if(j == 0){
    aic_comp_piece <- aic_piece %>% 
      transmute(jump, rep = -1, Parameter = "BASELINE", virus = NA, season = NA, age_group = NA, best = NA, loCI = NA, hiCI = NA, effSize = NA, 
                K, maxLoglik, AIC) %>% 
      unique
  } else {
    aic_comp_piece <- aic_piece %>% 
      filter(Parameter %in% c("lromega", "lrpsi", "delta", "sigma", "theta"))
    
  }
      
  if(is.null(aic_table)){
    aic_table <- aic_piece 
    aic_comp <- aic_comp_piece
  } else {
    aic_table <- rbind(aic_table, aic_piece)
    aic_comp <- rbind(aic_comp, aic_comp_piece)
  }
}

# whole table
aic_comp %>% 
  filter(is.na(effSize) | effSize > 50) %>%
  arrange(AIC) %>% 
  View

#' @UPDATE-TABLE2 
param_description <- read_csv2("input/param_description.csv", col_types = c(virus = "c"))

table2 <- aic_comp %>%
  ungroup %>% 
  mutate(across(c("best", "loCI", "hiCI"), ~ifelse(grepl("^lr", Parameter), 10^.x, .x))) %>% 
  mutate(across("Parameter", ~ifelse(grepl("^lr", .x), gsub("lr", "1/", .x), .x))) %>% 
  # filter(is.na(effSize) | effSize > 50) %>% 
  left_join(param_description %>% transmute(jump, Parameter, virus, description)) %>% 
  transmute(jump, rep, 
            Parameter = paste0(Parameter, replace_na(virus, "")), 
            description, 
            # Description = NA, 
            K, effSize, maxLoglik, 
            deltaAIC = (aic_comp %>% filter(Parameter == "BASELINE") %>% pull(AIC) %>% min) - AIC, 
            Baseline = case_when(grepl("omega", Parameter) ~ 0.001, 
                                 grepl("psi", Parameter) ~ 0.001, 
                                 T ~ 1), 
            estimate = paste0(round(best, 1), " (", round(loCI, 1), "-", round(hiCI, 1), ")"))# %>% 
  # arrange(-deltaAIC)
 
table2 %>% 
  arrange(-deltaAIC) %>% 
  write_csv2("figtab/Tab2_interactionFits.csv")




aic_table %>% 
  filter(jump %in% c(0, 3)) %>% 
  filter(Parameter == "beta") %>%
  # filter(Parameter == "eta") %>% 
  ggplot(aes(y = as.factor(jump), colour = virus)) + geom_point(aes(x = best)) + 
  geom_errorbarh(aes(xmin = loCI, xmax = hiCI)) +
  facet_grid(.~season, scales = "free_x") + 
  theme_bw()

##### TABLE S4

tableS4_comprehensive <- aic_table %>% 
  filter(jump == 3) %>% 
  mutate(season = ifelse(season == "common", season, gsub("^([0-9][0-9])", "20\\1/", season))) %>%
  mutate(W1 = ifelse(season == "common", NA, as.character(ISOweek::ISOweek2date(paste0("20", gsub(".*/", "", season), "-W01-1"))))) %>% 
  # mutate(across(c("best", "loCI", "hiCI"), ~ifelse(Parameter == "lh", 10^.x, .x))) %>% 
  # mutate(across(c("best", "loCI", "hiCI"), ~ifelse(Parameter == "eta", 
  #                                                  as.character(format.Date(as.Date(W1) - 7 + .x*7, "%b %d")), 
  #                                                  as.character(format(signif(.x, 3), scientific = F, drop0trailing = T))))) %>% 
  mutate(across(c("best", "loCI", "hiCI"), ~ifelse(Parameter == "lh", 10^.x, .x))) %>%
  mutate(Parameter = ifelse(Parameter == "lh", "h", Parameter)) %>% 
  mutate(across(c("best", "loCI", "hiCI"), ~case_when(Parameter == "eta" ~ as.character(format.Date(as.Date(W1) - 7 + .x*7, "%b %d")),
                                                      Parameter == "beta" ~ as.character(format(signif(.x, 3), scientific = F, drop0trailing = T)),
                                                      Parameter %in% c("h", "kappa") ~ paste0(format(signif(.x*100, 3), scientific = F, drop0trailing = T), "%"), 
                                                      T ~ NA)
                                                   )) %>% 
  mutate(Parameter = factor(Parameter, levels = c("beta", "eta", "h", "kappa"))) %>% 
  mutate(age_group = factor(age_group, levels = c(age_vec, age2_vec, "all"))) %>% 
  arrange(Parameter, virus, season, age_group) %>% 
  filter(!is.na(Parameter)) %>% 
  mutate(virus = c("Influenza", "RSV")[as.numeric(virus)]) %>% 
  relocate(age_group, .before = "season")


tableS4_comprehensive %>% 
  transmute(Parameter, virus, age_group, season, estimate = paste0(best, " (", loCI, "-", hiCI, ")")) %>% 
  write_csv2("figtab/TabS4_allparamFits.csv")  

