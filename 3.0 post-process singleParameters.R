library(tidyverse)

virus1 = "Influenza"
virus2 = "RSV"
id = "interOnly_" #
read_id = "simpleInter_"
jumpfile = "interactionOnly"

# param_jumpsWRONG = read_csv(
#   file = paste0("input/lrodin2v_jumps_", "interactiononlyDoubleWRONG", ".csv"), 
#   show_col_types = F
# ) 

param_jumps = read_csv(
  file = paste0("input/lrodin2v_jumps_", "interactiononlyDouble", ".csv"), 
  show_col_types = F
) 

baseline_post %>% 
  select(-Iteration, -Loglik, -Accepted, -Season) %>% 
  select_if(~length(unique(.x)) > 1) %>% 
  ncol

param_jumps %>% 
  select(jump, interaction_params) %>% 
  View

K = 66

this_rep = 3
this_jump = 3
for(this_rep in 0:5){
  baseline_stem = paste0("output/", read_id, virus1, "+", virus2, "_jump", 0, "_rep", this_rep)
  baseline_post <- read_csv(file = paste0(baseline_stem, "_posterior.csv"), show_col_types = F) %>% 
    mutate(across(starts_with(c("delta", "sigma", "theta")), ~1)) %>% 
    mutate(across(starts_with(c("omega", "psi")), ~-3)) %>% 
    rename_all(~gsub("^omega", "lromega", .x)) %>% 
    rename_all(~gsub("^psi", "lrpsi", .x))
  
  
  baseline_loglik <- baseline_post %>% 
    arrange(-Loglik) %>% 
    slice(1) %>% 
    pull(Loglik)
  
  baseline_aic = 2*K - 2*baseline_loglik
    
  for(this_jump in 55:1){
    print(c(this_rep, this_jump))
    
    full_stem = paste0("output/", id, virus1, "+", virus2, "_jump", this_jump, "_rep", this_rep)
    
    mcmc_post <- read_csv2(file = paste0(full_stem, "_posterior.csv"), show_col_types = F)

  
    # mcmc_post %>% 
    #   select(-Iteration, -Loglik, -Accepted, -Season) %>% 
    #   select_if(~length(unique(.x)) > 1) %>% 
    #   names
    
    fit_params <- param_jumps %>% 
      filter(jump == this_jump) %>% 
      select(-jump) %>% 
      select(-starts_with("beta"), -starts_with("eta"), -starts_with("pi"), -starts_with("lh"), -starts_with("kappa"), -starts_with("alpha"), -starts_with("gamma")) %>% 
      select_if(~.x > 0) %>% names
    
    baseline_paramDefault <- baseline_post[fit_params] %>% unique %>% unlist
    
    if(this_jump <= 10){
      AICimprovement_table_piece <- mcmc_post %>% 
        select(Iteration, any_of(fit_params), Loglik) %>% 
        rename(fit1 = !!(fit_params[1])) %>% 
        mutate(improvement = Loglik >= baseline_loglik) %>% 
        mutate(AIC = 2*(K+length(fit_params)) - 2*Loglik) %>% 
        mutate(AICimprovement = AIC <= baseline_aic) %>% 
        # arrange(-Loglik) %>% 
        filter(!is.na(Loglik)) %>% 
        summarise(rep = this_rep, 
                  jump = this_jump, 
                  param1 = fit_params[1], 
                  paramDefault1 = baseline_paramDefault[1], 
                  baseline_aic = baseline_aic, 
                  best_AIC = min(AIC), 
                  best1 = first(fit1[order(Loglik, decreasing = T)]), 
                  # lo_MCMC = min(fit1), 
                  # hi_MCMC = max(fit1), 
                  # lo_improvement = min(fit1[improvement]), 
                  # hi_improvement = max(fit1[improvement]),                 
                  lo1_AICimprovement = min(fit1[AICimprovement]) %>% ifelse(is.finite(.), ., NA), 
                  hi1_AICimprovement = max(fit1[AICimprovement]) %>% ifelse(is.finite(.), ., NA), 
                  ESS1 = coda::effectiveSize(fit1)
        )
      
      if(F){
        mcmc_post %>% 
          select(Iteration, any_of(fit_params), Loglik) %>% 
          # pull(sigma2) %>% coda::effectiveSize()
          # ggplot(aes(x = Iteration, y = sigma1)) + geom_line()
        ggplot(aes(x = sigma1, y = Loglik)) + geom_point()
        
      }
    } else {
      # if(length(fit_params) == 1){
      #   AICimprovement_table_piece = NULL
      # } else {
      #   this_jump = param_jumps %>% 
      #     rename(fit1 = !!(fit_params[1]), fit2 = !!(fit_params[2])) %>% 
      #     filter(fit1 > 0 & fit2 > 0) %>% 
      #     pull(jump)
        
      AICimprovement_table_piece <- mcmc_post %>% 
        select(Iteration, any_of(fit_params), Loglik) %>% 
        rename(fit1 = !!(fit_params[1]), fit2 = !!(fit_params[2])) %>% 
        mutate(improvement = Loglik >= baseline_loglik) %>% 
        mutate(AIC = 2*(K+length(fit_params)) - 2*Loglik) %>% 
        mutate(AICimprovement = AIC <= baseline_aic) %>% 
        # arrange(-Loglik) %>% 
        filter(!is.na(Loglik)) %>% 
        summarise(rep = this_rep, 
                  jump = this_jump, 
                  param1 = fit_params[1], 
                  paramDefault1 = baseline_paramDefault[1], 
                  param2 = fit_params[2], 
                  paramDefault2 = baseline_paramDefault[2], 
                  baseline_aic = baseline_aic, 
                  best_AIC = min(AIC), 
                  K = (K+length(fit_params)), 
                  best1 = first(fit1[order(Loglik, decreasing = T)]), 
                  best2 = first(fit2[order(Loglik, decreasing = T)]), 
                  # best1 = first(fit1), 
                  # best2 = first(fit2), 
                  # lo_MCMC = min(fit1), 
                  # hi_MCMC = max(fit1), 
                  # lo_improvement = min(fit1[improvement]), 
                  # hi_improvement = max(fit1[improvement]),                 
                  lo1_AICimprovement = min(fit1[AICimprovement]) %>% ifelse(is.finite(.), ., NA), 
                  hi1_AICimprovement = max(fit1[AICimprovement]) %>% ifelse(is.finite(.), ., NA), 
                  lo2_AICimprovement = min(fit2[AICimprovement]) %>% ifelse(is.finite(.), ., NA), 
                  hi2_AICimprovement = max(fit2[AICimprovement]) %>% ifelse(is.finite(.), ., NA), 
                  ESS1 = coda::effectiveSize(fit1), 
                  ESS2 = coda::effectiveSize(fit2))
      # }
      
    }
    
    
    if(this_rep == 0 & this_jump == 55){
      AICimprovement_table = AICimprovement_table_piece
    } else {
      AICimprovement_table = bind_rows(AICimprovement_table, AICimprovement_table_piece)
    }
      
    # mcmc_post %>% 
    #   select(Iteration, any_of(fit_params), Loglik) %>% 
    #   mutate(improvement = Loglik >= baseline_loglik) %>% 
    #   ggplot(aes(x = get(fit_params), y = Loglik, colour = improvement)) + geom_point()
    # 
    # mcmc_post %>% 
    #   select(Iteration, any_of(fit_params), Loglik) %>% 
    #   ggplot(aes(x = Iteration, y = get(fit_params))) + geom_line()
    # 
    # mcmc_post %>% 
    #   filter()
    #   
    
  }  
}


AICimprovement_table %>% 
  filter(rep == 0) %>% 
  filter(jump < 11) %>%
  # filter(rep == 0) %>% 
  arrange(best_AIC)

AICimprovement_table %>% 
  filter(rep == 0) %>% 
  filter(jump < 11) %>%
  # filter(rep == 0) %>% 
  arrange(best_AIC) %>% 
  mutate(across(c("paramDefault1", "best1", "lo1_AICimprovement", "hi1_AICimprovement"), ~ifelse(grepl("^lr", param1), 10^.x, .x))) %>%
  mutate(param1 = gsub("^lr", "1/", param1)) %>% 
  transmute(jump, param1, K = 67, deltaAIC = baseline_aic - best_AIC, Baseline = paramDefault1, Estimate = paste0(round(best1, 1), " (", round(lo1_AICimprovement, 1), "-", round(hi1_AICimprovement, 1), ")")) %>% 
  {rbind(., tibble(jump = 0, param1 = "BASELINE", K = 66, deltaAIC = 0, Baseline = NA, Estimate = NA))} %>% 
  left_join(param_description %>% transmute(param1 = paste0(Parameter, ifelse(is.na(virus), "", virus)), description)) %>% 
  arrange(-deltaAIC) %>% 
  relocate(description, .after = "param1") %>% 
  mutate(deltaAIC = case_when(deltaAIC > 0 ~ format(round(deltaAIC, 1), digits = 2), 
                              deltaAIC == 0 ~ "0", 
                              T ~ "<0")) %>% 
  mutate(Estimate = ifelse(deltaAIC == "<0", NA, Estimate)) %>% 
  write_csv2("figtab/Tab2_interactionFits_twoStage.csv")





AICimprovement_table %>% 
  filter(jump < 11) %>%
  mutate(across(c("paramDefault1", "best1", "lo1_AICimprovement", "hi1_AICimprovement"), ~ifelse(grepl("^lr", param1), 10^.x, .x))) %>%
  mutate(param1 = gsub("^lr", "1/", param1)) %>% 
  separate(col = "param1", into = c("param", "virus"), sep = -1, remove = F) %>% 
  mutate(virus = c("Influenza", "RSV")[as.numeric(virus)]) %>% 
  mutate(valid = ifelse(best_AIC < baseline_aic, "Improvement", "No improvement")) %>% 
  ggplot(aes(y = rep, colour = valid)) + 
  geom_point(aes(x = best1)) + 
  geom_errorbarh(aes(xmin = lo1_AICimprovement, xmax = hi1_AICimprovement)) + 
  facet_grid(virus~param, scales = "free_x") + 
  theme_bw() + 
  labs(x = "Parameter value", y = "Repeat", colour = "")

ggsave(filename = "figtab/compare_repetitions_improvement.png", width = 28, height = 12, units = "cm", dpi = 300)
  

AICimprovement_table %>% 
  arrange(best_AIC) %>% 
  # filter(jump < 11) %>% 
  filter(best_AIC < baseline_aic) %>% 
  ggplot(aes(colour = rep)) + 
  geom_point(aes(x = best1, y = best2)) + 
  geom_errorbarh(aes(xmin = lo1_AICimprovement, xmax = hi1_AICimprovement)) + 
  facet_wrap(.~param1, scales = "free_x")



# FOR THE TABLE, IDENTIFY PARAMETER RANGES WHICH ARE AN IMPROVEMENT ON BASELINE (FOR EACH REP?)
# USE THIS FOR THE CIs?
# OR MAKE THE CIs consistently poisson





