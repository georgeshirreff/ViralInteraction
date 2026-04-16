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

#odin
# library(odin2)
library(tidyverse) 
library(odin)
library(FME)


source("odin_m5.R")

# cmds = commandArgs(trailingOnly = T)

if(grepl("C:", getwd())){
  virus1 = "Influenza"
  virus2 = "RSV"
  include_ili = (virus1 == "Influenza") & (virus2 == "RSV") # should the ILI data be included?
  id = "estimateSingleInteractions_" #
  append = F #should the data be appended to an existing analysis?
  method = "cyc" #"basic" #"mod" #
  this_jump = 9 #which line in the proposal file should be used (baseline is 0)
  jumpfile = "interactiononly" #"all" #"all_BIGeta_BIGompsi" #"all" #"interactiononly" #"etabeta" #"kappalh" #
  this_rep = 0 #repeats, with the same parameters and proposals
  new_iter_total = 1000 #the number of new iterations to add
  thin = 1 #how many iterations should be saved (to save file space, increase this number to 10, 100, 1000 etc)
  save_every = 1000 #how often should the results be saved?
  
} else {
  cmds = commandArgs(trailingOnly = T)
  
  virus1 = cmds[1]
  virus2 = cmds[2]
  include_ili = (virus1 == "Influenza") & (virus2 == "RSV")
  id = cmds[3]
  append = cmds[4]=="append"
  method = cmds[5]
  this_jump = cmds[6]
  jumpfile = cmds[7]
  this_rep = cmds[8]
  new_iter_total = as.numeric(cmds[9])
  thin = as.numeric(cmds[10])
  save_every = 1000
  
}
eta_per_beta = 5

full_stem = paste0("output/", id, virus1, "+", virus2, "_jump", this_jump, "_rep", this_rep)
print(full_stem)
print(commandArgs(trailingOnly = T))
if(save_every < thin) stop(paste("save frequency", save_every, "can't be less than thinning factor", thin))

# read data
allcatch_inc_age <- read_tsv("data/allcatch_inc_age.tsv", show_col_types = F)
allcatch_coinc_age <- read_tsv("data/allcatch_coinc_age.tsv", show_col_types = F)
# spenc <- read_csv("data/Spencer_params.csv")
# age_contact_standardised <- read_tsv("data/age_contact_standardised.tsv")
ili <- read_tsv("data_public/allcatch_ili_age2.tsv", show_col_types = F) %>% 
  filter(!(season %in% c("20/21"))) %>% 
  rename(age2_group = age_group, t = iso_standard) %>% 
  select(-virus, -CatchmentPop)

census2022 <- read_csv("data_public/census2022.csv", show_col_types = F)

# this_virus = "Influenza"
param_lower = read_csv(file = paste0("input/lrodin2v_lower.csv"), show_col_types = F) %>% unlist
param_upper = read_csv(file = paste0("input/lrodin2v_upper.csv"), show_col_types = F) %>% unlist
param_jump = read_csv(
  # file = paste0("input/odin2v_jumps.csv")
  # file = paste0("input/odin2v_jumps_etaonly.csv")
  # file = paste0("input/odin2v_jumps_etabeta.csv")
  # file = paste0("input/odin2v_jumps_etabetakappa.csv")
  # file = paste0("input/odin2v_jumps_etabetakappalh.csv")
  file = paste0("input/lrodin2v_jumps_", jumpfile, ".csv"), 
  show_col_types = F
  
                      ) %>% 
  filter(jump == this_jump) %>% 
  select(-jump) %>% 
  unlist 

fit_param = read_csv(
  file = paste0("input/lrodin2v_jumps_", jumpfile, ".csv"), 
  show_col_types = F
) %>% 
  filter(jump == this_jump) %>% 
  select(-jump) %>% 
  select(-starts_with("beta"), -starts_with("eta"), -starts_with("pi"), -starts_with("lh"), -starts_with("kappa"), -starts_with("alpha"), -starts_with("gamma")) %>% 
  select_if(~.x > 0) %>% names


dat <- allcatch_inc_age %>% 
  filter(#virus == this_virus, 
         # age_group == this_age,
         # age_group == "6-12m", 
         !(season %in% c("20/21")), 
         # season == this_season
  )

codat <- allcatch_coinc_age %>% 
  filter(#virus == this_virus, 
    # age_group == this_age,
    # age_group == "6-12m", 
    !(season %in% c("20/21")), 
    # season == this_season
  )




age_vec = dat$age_group %>% unique
age2_vec = ili$age2_group %>% unique
season_vec = dat$season %>% unique
n_seasons = season_vec %>% length
n_ages = age_vec %>% length


dat_2virus <- rbind(dat %>% filter(virus == virus1) %>% mutate(virus = 1), 
                    dat %>% filter(virus == virus2) %>% mutate(virus = 2), 
                    codat %>% filter(CoVirus == paste0(virus1, "+", virus2)) %>% mutate(CoVirus = 12) %>% rename(virus = CoVirus)) %>% 
  transmute(season, t = iso_standard, age_group, virus, Inc, CatchmentPop) %>% 
  arrange(season, t, match(age_group, age_vec), virus)


data_catchments_weights = dat_2virus %>% 
  select(age_group, season, CatchmentPop) %>% unique %>% 
  left_join(census2022, by = "age_group") %>% 
  mutate(age2_group = factor(age2_age(age_group), levels = age2_vec)) %>% 
  group_by(season, age2_group) %>% 
  mutate(age2_weight = CensusPopulation/sum(CensusPopulation)) %>% 
  ungroup


if(!append){
  read_params <- read_csv(file =
                            paste0("output/", id, virus1, "+", virus2, "_jump", this_jump, "_rep", this_rep, "_bestParams.csv")
                            , show_col_types = F)
  
  mcmc_post = read_params %>% 
    mutate(Iteration = 0, Loglik = NA, Accepted = NA_real_, Season = "All") %>% 
    select(Iteration, everything())

  last_iter = 0
  
  param_start = read_params %>% 
    select(-contains("oglik"), -contains("teration")) %>% unlist

} else {
  mcmc_post = read_csv(file = paste0(full_stem, "_posterior.csv"), show_col_types = F)
  last_iter = mcmc_post$Iteration %>% last()

  mcmc_post %>%
    slice(n()) %>%
    select(Loglik)
  
  param_start <- mcmc_post %>% 
    slice(n()) %>% 
    select(-Iteration, -Loglik, -Accepted, -Season) %>%
    unlist
} 


if(method == "mod"){
  #### METHOD modMCMC ####
  # runs adaptive MCMC
  
  
  # iter_chunk = 200

  target_iter = last_iter + new_iter_total
  it = last_iter + save_every
  while(it <= target_iter){
    print(paste0(virus1, "+", virus2, ": ", last_iter + 1, " to ", it))
    
    #fit model
    mcmc_out <- modMCMC(niter = save_every, verbose = 50, updatecov = 50, ntrydr = 1, 
                        p = param_start, lower = param_lower, upper = param_upper, jump = param_jump, 
                        f = compute_ss, generator = m5_generator, dat_2virus = dat_2virus, ili, include_ili = include_ili, data_catchments_weights = data_catchments_weights)
  
    
    mcmc_lines = mcmc_out$pars %>% 
      as_tibble %>% 
      mutate(Iteration = (last_iter + 1):it, Loglik = mcmc_out$SS*-0.5, 
             Accepted = NA_real_, 
             Season = "All")
    
    mcmc_post <- rbind(mcmc_post, mcmc_lines) %>% 
      filter((Iteration %% thin) == 0)
  
  
    write_csv(mcmc_post, file = paste0(full_stem, "_posterior.csv"))
  
    best_params = mcmc_post[rev(which.max(mcmc_post$Loglik))[1], ] %>% select(-Iteration, -Accepted, -Season)
    best_params %>% write_csv(paste0(full_stem, "_bestParams.csv"))
    
    print(best_params["Loglik"])
  
    # prepare for next chunk
    last_iter = it
    it = last_iter + save_every
    param_start = mcmc_post[nrow(mcmc_post), ] %>% select(-Iteration, -Loglik, -Accepted, -Season) %>% unlist
  
  }
}

# set.seed(1)
if(method %in% c("cyc", "basic")){
  #### METHOD cyclic ####
  # "cyc" runs through a cycle of proposals by season
  # "basic" proposes all parameters simultaneously
  
  current_params = param_start
  
  if(length(current_params) != length(param_jump)) stop("starting parameters and jump params are not the same length")
  
  current_state <- generate_ModelHosp_season(pars = current_params, generator = m5_generator, dat_2virus, data_catchments_weights)
  current_ModelHosp <- current_state$hosp
  
  
  if(include_ili){
    current_ModelILI <- current_state$ili
    current_loglik = get_loglik_hosp(dat_2virus, current_ModelHosp) + 
      get_loglik_ili(ili, current_ModelILI)
  } else {
    current_loglik = get_loglik_hosp(dat_2virus, current_ModelHosp)
  }
  
  # current_loglik = sum(dpois(dat_2virus$Inc, current_ModelHosp$Inc + 1e-25, log = T), na.rm = T) + 
  #   # sum(dpois(round(ili$IncRate1e5*1e7), lambda = (current_ModelILI$IncRate1e5 + 1e-25)*1e7, log = T), na.rm = T)
  #   sum(dexp(ili$IncRate1e5*1e6, lambda = current_ModelILI$IncRate1e5*1e7, log = T), na.rm = T)
  

  
  # mcmc_post <- read_csv(file = paste0("output/", id, virus1, "+", virus2, "_posterior.csv"))
  # mcmc_post <- mcmc_post %>% mutate(Iteration = 1:n() - 1)
  # write_csv(mcmc_post, file = paste0("output/", id, virus1, "+", virus2, "_posterior.csv"))
  
  it = last_iter
  # season_index = 0
  # accepted = NA_real_
  target_iter = last_iter + new_iter_total
  
  if(method == "cyc"){
    marker_vec = 0:(n_seasons*eta_per_beta)
  } else if (method == "basic"){
    marker_vec = 0
  }
  while(it < target_iter){
    for(season_index_marker in marker_vec){
      
      season_index = ceiling(season_index_marker/eta_per_beta)
      last_iter = it
      it = it + 1
      
      if(it > target_iter) break
      
      if(season_index == 0){
        
        proposal_params = Vectorize(FUN = rnorm)(n = 1, mean = current_params, sd = param_jump) %>% 
          pmax(param_lower) %>% 
          pmin(param_upper) %>% 
          setNames(names(current_params))
        
        proposal_state <- generate_ModelHosp_season(pars = proposal_params, generator = m5_generator, dat_2virus = dat_2virus, data_catchments_weights)

        # proposal_loglik = sum(dpois(dat_2virus$Inc, proposal_ModelHosp$Inc + 1e-25, log = T), na.rm = T) + 
          # sum(dpois(ili$IncRate1e5*1e6, lambda = (proposal_ModelILI$IncRate1e5 + 1e-25)*1e6, log = T), na.rm = T)
        
      } else {
        param_jump_season = param_jump
        param_jump_season[case_when(grepl(gsub("/", "", season_vec[season_index]), names(param_jump)) ~ F, 
                                    grepl("^beta", names(param_jump)) ~ T, 
                                    grepl("^eta", names(param_jump)) ~ T, 
                                    T ~ F) | 
                            ((!grepl("eta", names(param_jump))) & (season_index_marker %% eta_per_beta != 0))] = 0
        
        proposal_params = Vectorize(FUN = rnorm)(n = 1, mean = current_params, sd = param_jump_season) %>% 
          pmax(param_lower) %>% 
          pmin(param_upper) %>% 
          setNames(names(current_params))
        
        proposal_state <- generate_ModelHosp_season(pars = proposal_params, generator = m5_generator, dat_2virus = dat_2virus, data_catchments_weights, input_season_index = season_index, current_state$carry)
      }
      
      proposal_ModelHosp <- proposal_state$hosp
      # proposal_ModelHosp <- generate_ModelHosp_season(pars = proposal_params, generator = m5_generator, dat_2virus = dat_2virus)
      # proposal_loglik = sum(dpois(dat_2virus$Inc, proposal_ModelHosp$Inc + 1e-25, log = T), na.rm = T)
      if(include_ili){
        proposal_ModelILI <- proposal_state$ili
        proposal_loglik = get_loglik_hosp(dat_2virus, proposal_ModelHosp) + 
          get_loglik_ili(ili, proposal_ModelILI)
      } else {
        proposal_loglik = get_loglik_hosp(dat_2virus, proposal_ModelHosp)
      }

            
      lalpha = proposal_loglik - current_loglik # log of acceptance ratio
      alpha = exp(lalpha)
      
      u = runif(1) # draw a uniform variable which will be less than alpha with probability min(1, alpha)
      if (u < alpha) { # then accept the candidate
        current_params <- proposal_params
        current_loglik <- proposal_loglik
        current_state <- proposal_state
        # current_ModelHosp <- proposal_ModelHosp
        
        accepted = 1
        
      } else {
        accepted = 0
      }
      cat(c(it, fit_param, current_params[fit_param], current_loglik))
      cat('\n')
      
      mcmc_lines = c(Iteration = it, current_params, Loglik = current_loglik, Accepted = accepted, Season = ifelse(season_index == 0, "All", season_vec[season_index])) %>% t %>% as_tibble %>% 
        mutate(across(-Season, as.numeric))
      
      mcmc_post <- rbind(mcmc_post, mcmc_lines) %>% 
        filter((Iteration %% thin) == 0)
      
      if((it %% save_every) == 0){
        write_csv(mcmc_post, file = paste0(full_stem, "_posterior.csv"))
        
        best_params = mcmc_post[rev(which.max(mcmc_post$Loglik))[1], ] %>% select(-Accepted, -Season)
        best_params %>% write_csv(paste0(full_stem, "_bestParams.csv"))
        
        # diagnostic step
        if(F){
          check_current_loglik = compute_ss(pars = current_params, generator = m5_generator, dat_2virus = dat_2virus, ili = ili, include_ili = T, data_catchments_weights = data_catchments_weights, make_loglik = T)
          check_current_state <- compute_ss(pars = current_params, generator = m5_generator, dat_2virus = dat_2virus, ili = ili, include_ili = T, data_catchments_weights = data_catchments_weights, make_loglik = T, return_solution = T)
          check_current_ModelHosp = check_current_state$hosp
          check_current_ModelILI = check_current_state$ili
          tibble(Loglik = check_current_loglik) %>% 
            write_csv(file = paste0(full_stem, "_currentLoglik.csv"))
          check_current_ModelHosp %>% 
            write_csv(file = paste0(full_stem, "_ModelHosp.csv"))
          check_current_ModelHosp %>% 
            write_csv(file = paste0(full_stem, "_ModelILI.csv"))
        }
        
        print(c(it, current_loglik, best_params["Loglik"]))
      }
      
      
      param_start = current_params #for transfer back to modMCMC method
    }
    
  }
}


#######################
## some basic processing ##

if(F){
  
  id = "estimateSingleInteractions_"
  plot_post_2v(stem = "output/", id = id, virus1 = "Influenza", virus2 = "RSV", 
               jump = 1, rep = 1,
               burnin = 0, thin = 1,
               exclude_params = c(paste0(rep(c("alpha", "gamma"
                                               #, "omega", "psi", "delta", "sigma", "theta"
                                               ), each = 2), 1:2), paste0("pi_", age_vec))
               # exclude_params = c(paste0(rep(c("alpha", "gamma", "omega", "psi", "delta", "sigma", "theta"), each = 2), 1:2), paste0("pi_", age_vec))
               # exclude_params = c(paste0(rep(c("alpha", "gamma", "psi", "delta", "sigma", "theta"), each = 2), 1:2), paste0("pi_", age_vec))
               # exclude_params = c("alpha1", "gamma1", "alpha2", "gamma2", paste0("pi_", age_vec))
  )
  
  
  check_jump = 1
  interaction_params = param_start %>% names %>% rev %>% {.[1:10]} %>% rev
  if(check_jump == 0){
    interaction_params_excl = NULL
  } else {
    interaction_params_excl = interaction_params[-check_jump]
  }
  compare_reps(stem = "output/", id = "estimateSingleInteractions_", virus1 = "Influenza", virus2 = "RSV", 
               jump = check_jump, reps = 0:5,
               thin = 10,
               # exclude_params = c(paste0(rep(c("alpha", "gamma", "omega", "psi", "delta", "sigma", "theta"), each = 2), 1:2), paste0("pi_", age_vec))
               exclude_params = c(interaction_params_excl, paste0(rep(c("alpha", "gamma"), each = 2), 1:2), paste0("pi_", age_vec))
               # exclude_params = c(paste0(rep(c("alpha", "gamma", "psi", "delta", "sigma", "theta"), each = 2), 1:2), paste0("pi_", age_vec))
               # exclude_params = c("alpha1", "gamma1", "alpha2", "gamma2", paste0("pi_", age_vec))
  )
  
  plot_post_2v(stem = "output/", id = "estimateSingleInteractions_", virus1 = "Influenza", virus2 = "RSV", 
               jump = check_jump, rep = 0,
               burnin = 0, thin = 1,
               exclude_params = c(interaction_params_excl, paste0(rep(c("alpha", "gamma"
                                               #, "omega", "psi", "delta", "sigma", "theta"
               ), each = 2), 1:2), paste0("pi_", age_vec))
               # exclude_params = c(paste0(rep(c("alpha", "gamma", "omega", "psi", "delta", "sigma", "theta"), each = 2), 1:2), paste0("pi_", age_vec))
               # exclude_params = c(paste0(rep(c("alpha", "gamma", "psi", "delta", "sigma", "theta"), each = 2), 1:2), paste0("pi_", age_vec))
               # exclude_params = c("alpha1", "gamma1", "alpha2", "gamma2", paste0("pi_", age_vec))
  )
  
  chosen_rep <- tibble(jump = 0:10, chosenRep = c(0, 
                                                  3, 1, 0, 1, 2, 
                                                  0, 1, 0, 1, 4), 
         burnin = c(2.5e5, 
                    1e5, 2e5, 1e5, 2e5, 1e5, 
                    1e5, 1e5, 1.5e5, 1e5, 1e5))
  
  compare_reps(stem = "output/", id = "estimateSingleInteractions_", virus1 = "Influenza", virus2 = "RSV", 
               jump = 6, reps = 0,
               thin = 10,
               # exclude_params = c(paste0(rep(c("alpha", "gamma", "omega", "psi", "delta", "sigma", "theta"), each = 2), 1:2), paste0("pi_", age_vec))
               exclude_params = c(paste0(rep(c("alpha", "gamma"), each = 2), 1:2), paste0("pi_", age_vec))
               # exclude_params = c(paste0(rep(c("alpha", "gamma", "psi", "delta", "sigma", "theta"), each = 2), 1:2), paste0("pi_", age_vec))
               # exclude_params = c("alpha1", "gamma1", "alpha2", "gamma2", paste0("pi_", age_vec))
  )
  
  plot_acceptance_2v(stem = "output/", id = id, virus1 = virus1, virus2 = virus2, 
                     jump = 0, rep = 0)

  post %>% 
    select(Iteration, omega1, Loglik) %>% 
    ggplot(aes(x = Iteration, y = 1/omega1))+ geom_line()
    ggplot(aes(x = 1/omega1, y = Loglik)) + geom_point()
  
  post %>% 
    select(Iteration, Loglik, starts_with("eta")) %>% 
    pivot_longer(-c("Iteration", "Loglik")) %>% 
    ggplot(aes(x = Iteration, y = value)) + geom_line() + 
    facet_wrap(.~name, scales = "free")
  
  post %>% 
    select(Iteration, Loglik, starts_with("eta"), , starts_with("beta")) %>% 
    pivot_longer(-c("Iteration", "Loglik"), names_pattern = "([a-z]+)([12])_([0-9]+)", names_to = c(".value", "virus", "season")) %>% 
    ggplot(aes(x = eta, y = beta, colour = Loglik)) + geom_point() + 
    facet_grid(season~virus, scales = "free")
    
  
  best_params <- read_csv(paste0("output/", id = "estimateSingleInteractions_", virus1, "+", virus2, "_jump", 3, "_rep", 1, "_bestParams.csv"))
  
  #' @UPDATE-FIGURE4 
  
  plot_obj <- plot_model_data(pars = best_params %>%
                    select(-Loglik) %>%
                    unlist,
                  generator = m5_generator, virus1 = virus1, virus2 = virus2, dat_2virus = dat_2virus, ili = ili, data_catchments_weights = data_catchments_weights, 
                  # type = "ageseason"
                  # type = "age"
                  type = "season"
                  # type = "ili_ageseason"
                  # type = "ili_age"
                  # type = "ili_season", 
                  , include_ci = T
  ) 
  
  ggsave(plot_obj, filename = "~/SanOdin/output/best_season_plot.png", width = 30, height = 16, units = "cm")
  # 
  
  # ggsave(paste0("output/", id, virus1, "+", virus2, "_solution.png"))
  
}



  
