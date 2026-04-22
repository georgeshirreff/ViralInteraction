library(tidyverse) 
library(odin)
library(FME)

source("odin_m5.R")

# cmds = commandArgs(trailingOnly = T)

if(grepl("C:", getwd())){
  virus1 = "Influenza"
  virus2 = "RSV"
  include_ili = (virus1 == "Influenza") & (virus2 == "RSV") # should the ILI data be included?
  id = "interOnly_" #
  append = F #should the data be appended to an existing analysis?
  method = "mod" #"cyc" #"basic" #
  this_jump = 65 #which line in the proposal file should be used (baseline is 0)
  jumpfile = "interactiononlyDouble" #"all" #"all_BIGeta_BIGompsi" #"all" #"interactiononly" #"etabeta" #"kappalh" #
  this_rep = 5 #repeats, with the same parameters and proposals
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
param_jumps = read_csv(
  file = paste0("input/lrodin2v_jumps_", jumpfile, ".csv"), 
  show_col_types = F
) 
param_jump = param_jumps %>% 
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
  read_id = "simpleInter_"
  read_params <- read_csv(file =
                            paste0("output/", read_id, virus1, "+", virus2, "_jump", 0, "_rep", this_rep, "_bestParams.csv")
                          , show_col_types = F) %>% 
    mutate(across(starts_with(c("delta", "sigma", "theta")), ~1)) %>% 
    mutate(across(starts_with(c("omega", "psi")), ~-3)) %>% 
    rename_all(~gsub("^omega", "lromega", .x)) %>% 
    rename_all(~gsub("^psi", "lrpsi", .x)) #%>% 
    # mutate(across(starts_with(c("lromega", "lrpsi")), ~log10(1/.x))) 
    
  
  mcmc_post = read_params %>% 
    mutate(Iteration = 0, Loglik = NA, Accepted = NA_real_, Season = "All") %>% 
    select(Iteration, everything())
  
  last_iter = 0
  
  param_start = read_params %>% 
    select(-contains("oglik"), -contains("teration")) %>% unlist
  
} else {
  mcmc_post = read_csv2(file = paste0(full_stem, "_posterior.csv"), show_col_types = F)
  last_iter = mcmc_post$Iteration %>% last()
  
  mcmc_post %>%
    slice(n()) %>%
    select(Loglik)
  
  param_start <- mcmc_post %>% 
    slice(n()) %>% 
    select(-Iteration, -Loglik, -Accepted, -Season) %>%
    unlist
} 


# for(rep in 0:5){
#   post = read_csv(paste0("output/", "simpleInter_",  virus1, "+", virus2, 
#                          "_jump", 0, 
#                          "_rep", rep, 
#                          "_posterior.csv")) %>% 
#     arrange(-Loglik) %>% 
#     slice(1)
#   
#   if(rep == 0){
#     baseline_rep = post %>% mutate(rep = rep)
#   } else {
#     baseline_rep = rbind(baseline_rep, post %>% mutate(rep = rep))
#   }
#   
# }


# baseline_rep <- baseline_rep %>% 
#   mutate(across(starts_with(c("delta", "sigma", "theta")), ~1), 
#          across(starts_with(c("omega", "psi")), list(lr=~log10(1/.x)), .names = "{.fn}{.col}"))


param_lower
param_upper = param_upper[names(param_lower)]
param_jump = param_jump[names(param_lower)]
param_start = param_start[names(param_lower)]

# # baseline_rep %>% t
# # baseline_rep %>% names
# # param_lower %>% names
# # param_upper %>% names
# # param_jumps %>% names
# # this_jump = 1
# # this_rep = 0
# # for(this_rep in 0:5){
# #   for(this_jump in 0:10){
#     
# 
#     # param_jump = param_jumps %>% 
#     #   filter(jump == this_jump) %>% 
#     #   select(-jump) %>% 
#     #   unlist 
#     
#     # param_start = baseline_rep %>% 
#     #   filter(rep == this_rep) %>% 
#     #   select(-rep) %>% 
#     #   unlist
#     
#     mcmc_out <- modMCMC(niter = 1000, verbose = 10, updatecov = 100, ntrydr = 1, 
#                         p = param_start[param_jump > 0], 
#                         
#                         lower = param_lower[param_jump > 0], 
#                         upper = param_upper[param_jump > 0], 
#                         jump = param_jump[param_jump > 0], 
#                         f = compute_ss, 
#                         fixed_pars = param_start[param_jump == 0], 
#                         generator = m5_generator, dat_2virus = dat_2virus, ili, include_ili = include_ili, data_catchments_weights = data_catchments_weights)
#     
#     plot(mcmc_out$pars[, "lromega1"], mcmc_out$SS*-0.5)
#     
#   }
# }

target_iter = last_iter + new_iter_total
it = last_iter + save_every
while(it <= target_iter){
  print(paste0(virus1, "+", virus2, ": ", last_iter + 1, " to ", it))
  
  #fit model
  mcmc_out <- modMCMC(niter = save_every, verbose = 50, updatecov = 50, ntrydr = 1, 
                      p = param_start[param_jump > 0], 
                      
                      lower = param_lower[param_jump > 0], 
                      upper = param_upper[param_jump > 0], 
                      jump = param_jump[param_jump > 0], 
                      fixed_pars = param_start[param_jump == 0], 
                      # p = param_start, lower = param_lower, upper = param_upper, jump = param_jump, 
                      f = compute_ss, generator = m5_generator, dat_2virus = dat_2virus, ili, include_ili = include_ili, data_catchments_weights = data_catchments_weights)
  
  
  mcmc_lines = mcmc_out$pars %>% 
    as_tibble %>% 
    cbind(t(param_start[param_jump == 0])) %>% 
    select(names(param_jump)) %>% 
    mutate(Iteration = (last_iter + 1):it, Loglik = mcmc_out$SS*-0.5, 
           Accepted = NA_real_, 
           Season = "All")
    
  mcmc_post <- rbind(mcmc_post, mcmc_lines) %>% 
    filter((Iteration %% thin) == 0)
  
  write_csv2(mcmc_post, file = paste0(full_stem, "_posterior.csv"))
  
  best_params = mcmc_post[rev(which.max(mcmc_post$Loglik))[1], ] %>% select(-Iteration, -Accepted, -Season)
  best_params %>% write_csv2(paste0(full_stem, "_bestParams.csv"))
  
  print(best_params["Loglik"])
  
  # prepare for next chunk
  last_iter = it
  it = last_iter + save_every
  param_start = mcmc_post[nrow(mcmc_post), ] %>% select(-Iteration, -Loglik, -Accepted, -Season) %>% unlist
  
# }
}


