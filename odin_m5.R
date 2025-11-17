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

#odin_m1
m5_generator <- odin(verbose = T, {
  

lambda1[] <- import_rate1[i] + beta1[i] * (IS[i] + IE[i] + II[i]*sigma1 + IRv[i] + IRf[i] + IR[i])/if(N[i]==0) 1 else N[i]
lambda2[] <- import_rate2[i] + beta2[i] * (SI[i] + EI[i] + II[i]*sigma2 + RvI[i] + RfI[i] + RI[i])/if(N[i]==0) 1 else N[i]
# deriv(S[]) <- - lambda[i] * S[i]
# deriv(E[]) <- lambda[i] * S[i] - alpha * E[i]
# deriv(I[]) <- alpha * E[i] - gamma * I[i]
# deriv(Rv[]) <- gamma * I[i] - omega * Rv[i]
# deriv(Rf[]) <- omega * Rv[i] - psi * Rf[i]
# deriv(R[]) <- psi * Rf[i]
deriv(SS[]) <- -SS[i]*(lambda1[i]+lambda2[i])
deriv(ES[]) <- SS[i]*lambda1[i]-ES[i]*(alpha1+lambda2[i])
deriv(IS[]) <- ES[i]*alpha1-IS[i]*(gamma1+lambda2[i]*delta2)
deriv(RvS[]) <- IS[i]*gamma1-RvS[i]*(omega1+lambda2[i])
deriv(RfS[]) <- RvS[i]*omega1-RfS[i]*(psi1+0)
deriv(RS[]) <- RfS[i]*psi1-RS[i]*(lambda2[i])

deriv(SE[]) <- SS[i]*lambda2[i]-SE[i]*(lambda1[i]+alpha2)
deriv(EE[]) <- SE[i]*lambda1[i]+ES[i]*lambda2[i]-EE[i]*(alpha1+alpha2)
deriv(IE[]) <- EE[i]*alpha1+IS[i]*lambda2[i]*delta2-IE[i]*(gamma1+alpha2)
deriv(RvE[]) <- IE[i]*gamma1+RvS[i]*lambda2[i]-RvE[i]*(omega1+alpha2)
deriv(RfE[]) <- RvE[i]*omega1+RfS[i]*0-RfE[i]*(psi1+alpha2)
deriv(RE[]) <- RfE[i]*psi1+RS[i]*lambda2[i]-RE[i]*(alpha2)

deriv(SI[]) <- SE[i]*alpha2-SI[i]*(lambda1[i]*delta1+gamma2)
deriv(EI[]) <- SI[i]*lambda1[i]*delta1+EE[i]*alpha2-EI[i]*(alpha1+gamma2)
deriv(II[]) <- EI[i]*alpha1+IE[i]*alpha2-II[i]*(gamma1+gamma2)
deriv(RvI[]) <- II[i]*gamma1+RvE[i]*alpha2-RvI[i]*(omega1+gamma2)
deriv(RfI[]) <- RvI[i]*omega1+RfE[i]*alpha2-RfI[i]*(psi1+gamma2)
deriv(RI[]) <- RfI[i]*psi1+RE[i]*alpha2-RI[i]*(gamma2)

deriv(SRv[]) <- SI[i]*gamma2-SRv[i]*(lambda1[i]+omega2)
deriv(ERv[]) <- SRv[i]*lambda1[i]+EI[i]*gamma2-ERv[i]*(alpha1+omega2)
deriv(IRv[]) <- ERv[i]*alpha1+II[i]*gamma2-IRv[i]*(gamma1+omega2)
deriv(RvRv[]) <- IRv[i]*gamma1+RvI[i]*gamma2-RvRv[i]*(omega1+omega2)
deriv(RfRv[]) <- RvRv[i]*omega1+RfI[i]*gamma2-RfRv[i]*(psi1+omega2)
deriv(RRv[]) <- RfRv[i]*psi1+RI[i]*gamma2-RRv[i]*(omega2)

deriv(SRf[]) <- SRv[i]*omega2-SRf[i]*(0+psi2)
deriv(ERf[]) <- SRf[i]*0+ERv[i]*omega2-ERf[i]*(alpha1+psi2)
deriv(IRf[]) <- ERf[i]*alpha1+IRv[i]*omega2-IRf[i]*(gamma1+psi2)
deriv(RvRf[]) <- IRf[i]*gamma1+RvRv[i]*omega2-RvRf[i]*(omega1+psi2)
deriv(RfRf[]) <- RvRf[i]*omega1+RfRv[i]*omega2-RfRf[i]*(psi1+psi2)
deriv(RRf[]) <- RfRf[i]*psi1+RRv[i]*omega2-RRf[i]*(psi2)

deriv(SR[]) <- SRf[i]*psi2-SR[i]*(lambda1[i])
deriv(ER[]) <- SR[i]*lambda1[i]+ERf[i]*psi2-ER[i]*(alpha1)
deriv(IR[]) <- ER[i]*alpha1+IRf[i]*psi2-IR[i]*(gamma1)
deriv(RvR[]) <- IR[i]*gamma1+RvRf[i]*psi2-RvR[i]*(omega1)
deriv(RfR[]) <- RvR[i]*omega1+RfRf[i]*psi2-RfR[i]*(psi1)
deriv(RR[]) <- RfR[i]*psi1+RRf[i]*psi2

deriv(CumInc1_simple[]) <- alpha1 * (ES[i] + EE[i] + ER[i]) #incidence unaffected by interaction
# deriv(CumInc1_refract[]) <- alpha1 * (ERf[i] * 0) #this is always zero so ignore
deriv(CumInc1_coinf[]) <- alpha1 * (EI[i] * delta1) #incidence of virus 1 with coinfection with virus 2
deriv(CumInc1_hitch[]) <- alpha1 * (ERv[i]) #incidence of virus 1 with codetection of virus 2

deriv(CumInc2_simple[]) <- alpha2 * (SE[i] + EE[i] + RE[i]) #incidence unaffected by interaction
# deriv(CumInc1_refract[]) <- alpha1 * (ERf[i] * 0) #this is always zero so ignore
deriv(CumInc2_coinf[]) <- alpha2 * (IE[i] * delta2) #incidence of virus 2 with coinfection with virus 1
deriv(CumInc2_hitch[]) <- alpha2 * (RvE[i]) #incidence of virus 2 with codetection of virus 1

import_rate1[] <- rho*if (t < eta1[i] - 0.75) 0 else if (t < eta1[i] - 0.25) (t - (eta1[i] - 0.75))/0.5 else if (t < eta1[i] + 0.25) 1 else if (t < eta1[i] + 0.75) (t - (eta1[i] + 0.75))/-0.5 else 0
import_rate2[] <- rho*if (t < eta2[i] - 0.75) 0 else if (t < eta2[i] - 0.25) (t - (eta2[i] - 0.75))/0.5 else if (t < eta2[i] + 0.25) 1 else if (t < eta2[i] + 0.75) (t - (eta2[i] + 0.75))/-0.5 else 0
output(import_rate1[]) = import_rate1[i]
output(import_rate2[]) = import_rate2[i]


# import_rate[] <- interpolate(import_t[i,], import_y[i,], "linear")
# import_rate[] <- interpolate(tt, import_y[i,], "linear")
# import_rate[] <- interpolate(tt, yy, "linear")
# import_rate[] = rho/50


beta1[] <- user()
beta2[] <- user()
dim(beta1) <- user()
dim(beta2) <- user()
eta1[] = user()
eta2[] = user()
dim(eta1) = user()
dim(eta2) = user()
N[] <- user()
dim(N) <- user()

alpha1 <- user()
alpha2 <- user()
gamma1 <- user()
gamma2 <- user()
omega1 <- user()
omega2 <- user()
psi1 <- user()
psi2 <- user()
delta1 <- user()
delta2 <- user()
sigma1 <- user()
sigma2 <- user()
rho <- 4e-4

n_seasons = length(beta1)

dim(import_rate1) = n_seasons
dim(import_rate2) = n_seasons

dim(lambda1) = n_seasons
dim(lambda2) = n_seasons
dim(SS) = n_seasons
dim(ES) = n_seasons
dim(IS) = n_seasons
dim(RvS) = n_seasons
dim(RfS) = n_seasons
dim(RS) = n_seasons

dim(SE) = n_seasons
dim(EE) = n_seasons
dim(IE) = n_seasons
dim(RvE) = n_seasons
dim(RfE) = n_seasons
dim(RE) = n_seasons

dim(SI) = n_seasons
dim(EI) = n_seasons
dim(II) = n_seasons
dim(RvI) = n_seasons
dim(RfI) = n_seasons
dim(RI) = n_seasons

dim(SRv) = n_seasons
dim(ERv) = n_seasons
dim(IRv) = n_seasons
dim(RvRv) = n_seasons
dim(RfRv) = n_seasons
dim(RRv) = n_seasons

dim(SRf) = n_seasons
dim(ERf) = n_seasons
dim(IRf) = n_seasons
dim(RvRf) = n_seasons
dim(RfRf) = n_seasons
dim(RRf) = n_seasons

dim(SR) = n_seasons
dim(ER) = n_seasons
dim(IR) = n_seasons
dim(RvR) = n_seasons
dim(RfR) = n_seasons
dim(RR) = n_seasons

dim(CumInc1_simple) = n_seasons
dim(CumInc1_coinf) = n_seasons
dim(CumInc1_hitch) = n_seasons

dim(CumInc2_simple) = n_seasons
dim(CumInc2_coinf) = n_seasons
dim(CumInc2_hitch) = n_seasons


initial(SS[]) = N[i] 
initial(ES[]) = 0
initial(IS[]) = 0
initial(RvS[]) = 0
initial(RfS[]) = 0
initial(RS[]) = 0

initial(SE[]) = 0
initial(EE[]) = 0
initial(IE[]) = 0
initial(RvE[]) = 0
initial(RfE[]) = 0
initial(RE[]) = 0

initial(SI[]) = 0
initial(EI[]) = 0
initial(II[]) = 0
initial(RvI[]) = 0
initial(RfI[]) = 0
initial(RI[]) = 0

initial(SRv[]) = 0
initial(ERv[]) = 0
initial(IRv[]) = 0
initial(RvRv[]) = 0
initial(RfRv[]) = 0
initial(RRv[]) = 0

initial(SRf[]) = 0
initial(ERf[]) = 0
initial(IRf[]) = 0
initial(RvRf[]) = 0
initial(RfRf[]) = 0
initial(RRf[]) = 0

initial(SR[]) = 0
initial(ER[]) = 0
initial(IR[]) = 0
initial(RvR[]) = 0
initial(RfR[]) = 0
initial(RR[]) = 0

initial(CumInc1_simple[]) = 0
initial(CumInc1_coinf[]) = 0
initial(CumInc1_hitch[]) = 0

initial(CumInc2_simple[]) = 0
initial(CumInc2_coinf[]) = 0
initial(CumInc2_hitch[]) = 0

})

options(dplyr.summarise.inform = F)

# pars = param_start
# generator = m5_generator
# dat_2virus

# generate_ModelHosp<- function(pars, generator, dat_2virus){
#   
#   # tic(msg = "initial variables")
#   model_CatchmentPop = 1000
#   age_vec = dat_2virus$age_group %>% unique
#   season_vec = dat_2virus$season %>% unique
#   data_catchments = dat_2virus %>% select(age_group, season, CatchmentPop) %>% unique
#   # toc(quiet = FALSE, func.toc = my.msg.toc)
#   
#   # tic("make age params")
#   age_params = tibble(incidentVirus = rep(1:2, each = length(age_vec)), 
#                       age_group = rep(age_vec, times = 2), pi = rep(pars[paste("pi", sep = "_", age_vec)], times = 2), 
#                       h = 10^pars[c(paste("lh1", sep = "_", age_vec), paste("lh2", sep = "_", age_vec))], 
#                       theta = rep(pars[paste0("theta", 1:2)], each = length(age_vec)))
#   # toc(quiet = FALSE, func.toc = my.msg.toc)
#   # age_params = tibble(age_group = age_vec, virus = 1, 
#   #                     h = 10^pars[paste("lh1", sep = "_", age_vec)],
#   #                     pi = pars[paste("pi", sep = "_", age_vec)])
#   # tic("make lpars")
#   process_params = c("alpha", "gamma", "omega", "psi", "delta", "sigma"
#                      # , "theta"
#                      )
#   lpars = c(list(N = rep(model_CatchmentPop, length(season_vec))), 
#             list(beta1 = pars[paste0("beta1_", gsub("/", "", season_vec))]), 
#             list(beta2 = pars[paste0("beta2_", gsub("/", "", season_vec))]), 
#             list(eta1 = pars[paste0("eta1_", gsub("/", "", season_vec))]), 
#             list(eta2 = pars[paste0("eta2_", gsub("/", "", season_vec))]), 
#             as.list(pars[c(paste0(process_params, 1), paste0(process_params, 2))])
#             # list(alpha1 = pars["alpha1"], gamma1 = pars["gamma1"], 
#             #      alpha2 = pars["alpha2"], gamma2 = pars["gamma2"])
#             )
#   # toc(quiet = FALSE, func.toc = my.msg.toc)
#   
#   # tic("generator")
#   mod = generator$new(user = lpars)
#   # toc(quiet = FALSE, func.toc = my.msg.toc)
#   # tic("execute")
#   out = mod$run(t = sort(unique(dat$iso_standard)), pretty = T)
#   # toc(quiet = FALSE, func.toc = my.msg.toc)
#   
#   # tic("tib_ModelInc_1000")
#   # EXTEND THIS TO Inc1, Inc2 and Inc12
#   # CORRESPOND WITH THE PAPER
#   tib_ModelInc_1000 <- out %>% as_tibble %>% mutate_all(as.double) %>% 
#     # select(fromHere) %>% 
#     select(t, starts_with("CumInc")) %>% 
#     mutate(across(starts_with("CumInc"), function(x) diff(c(0, x)))) %>% 
#     rename_all(~gsub("CumInc", "Inc", .x)) %>% 
#     # pivot_longer(-t, names_pattern = "([_A-Za-z]+)\\[([0-9]+)\\]", names_to = c(".value", "season_index")) %>% 
#     pivot_longer(-t, names_pattern = "([A-Za-z]+)([12])_([a-z]+)\\[([0-9]+)\\]", names_to = c(".value", "incidentVirus", "infection", "season_index")) %>%
#     mutate(season = season_vec[as.numeric(season_index)], 
#            incidentVirus = as.numeric(incidentVirus)) %>% 
#     select(-season_index)
#   # toc(quiet = FALSE, func.toc = my.msg.toc)
#   
#   # tic("tib_ModelInc")
#   tib_ModelInc <- tib_ModelInc_1000 %>% 
#     transmute(t, season, incidentVirus, infection, ModelInc_1000 = Inc, model_CatchmentPop) %>% 
#     left_join(age_params, by = "incidentVirus", relationship = "many-to-many") %>%  #expands this file according to age_params
#     left_join(data_catchments, by = c("season", "age_group")) %>% 
#     mutate(ModelInc = ModelInc_1000/1000*CatchmentPop) %>% 
#     select(-ModelInc_1000)
#   # toc(quiet = FALSE, func.toc = my.msg.toc)
#   
#   # tic("tib_ModelHosp")
#   tib_ModelHosp <- tib_ModelInc %>% 
#     mutate(ModelHosp1 = ModelInc*case_when(infection == "simple" & incidentVirus == 1 ~ pi*h,  
#                                            infection == "coinf" & incidentVirus == 1 ~ pi*h*(1-pi), 
#                                            infection == "hitch" & incidentVirus == 1 ~ pi*h*(1-pi),
#                                            T ~ 0), 
#            ModelHosp2 = ModelInc*case_when(infection == "simple" & incidentVirus == 2 ~ pi*h, 
#                                            infection == "coinf" & incidentVirus == 2 ~ pi*h*(1-pi), 
#                                            infection == "hitch" & incidentVirus == 2 ~ pi*h*(1-pi),
#                                            T ~ 0), 
#            ModelHosp12 = ModelInc*case_when(infection == "simple" & incidentVirus == 1 ~ 0,  
#                                             infection == "coinf" & incidentVirus == 1 ~ pi*h*pi*theta, 
#                                             infection == "hitch" & incidentVirus == 1 ~ pi*h*pi,
#                                             infection == "simple" & incidentVirus == 2 ~ 0, 
#                                             infection == "coinf" & incidentVirus == 2 ~ pi*h*pi*theta, 
#                                             infection == "hitch" & incidentVirus == 2 ~ pi*h*pi,
#                                             T ~ 0)) %>% 
#     group_by(season, t, age_group) %>% 
#     summarise(across(starts_with("ModelHosp"), sum), .group = "drop") %>% 
#     ungroup %>% 
#     pivot_longer(cols = starts_with("ModelHosp"), names_prefix = "ModelHosp", names_to = "virus", values_to = "Inc") %>% 
#     mutate(virus = as.numeric(virus)) %>% 
#     select(season, t, age_group, virus, Inc) %>% 
#     arrange(season, t, match(age_group, age_vec), virus)
#   # toc(quiet = FALSE, func.toc = my.msg.toc)
#   
#   tib_ModelHosp
# }
# 
# generate_ModelHosp_season<- function(pars, generator, dat_2virus, input_season_index = 0, tib_ModelHosp = NULL){
#   
#   if((input_season_index == 0) != is.null(tib_ModelHosp)) stop("either include both input_season_index and tib_ModelHosp, or neither")
#   # tic(msg = "initial variables")
#   model_CatchmentPop = 1000
#   age_vec = dat_2virus$age_group %>% unique
#   season_vec = dat_2virus$season %>% unique
#   data_catchments = dat_2virus %>% select(age_group, season, CatchmentPop) %>% unique
#   # toc(quiet = FALSE, func.toc = my.msg.toc)
#   
#   # tic("make age params")
#   age_params = tibble(incidentVirus = rep(1:2, each = length(age_vec)), 
#                       age_group = rep(age_vec, times = 2), pi = rep(pars[paste("pi", sep = "_", age_vec)], times = 2), 
#                       h = 10^pars[c(paste("lh1", sep = "_", age_vec), paste("lh2", sep = "_", age_vec))], 
#                       theta = rep(pars[paste0("theta", 1:2)], each = length(age_vec)))
#   # toc(quiet = FALSE, func.toc = my.msg.toc)
#   # age_params = tibble(age_group = age_vec, virus = 1, 
#   #                     h = 10^pars[paste("lh1", sep = "_", age_vec)],
#   #                     pi = pars[paste("pi", sep = "_", age_vec)])
#   # tic("make lpars")
#   process_params = c("alpha", "gamma", "omega", "psi", "delta", "sigma"
#                      # , "theta"
#   )
#   lpars = c(list(N = rep(model_CatchmentPop, length(season_vec))), 
#             list(beta1 = pars[paste0("beta1_", gsub("/", "", season_vec))]), 
#             list(beta2 = pars[paste0("beta2_", gsub("/", "", season_vec))]), 
#             list(eta1 = pars[paste0("eta1_", gsub("/", "", season_vec))]), 
#             list(eta2 = pars[paste0("eta2_", gsub("/", "", season_vec))]), 
#             as.list(pars[c(paste0(process_params, 1), paste0(process_params, 2))])
#             # list(alpha1 = pars["alpha1"], gamma1 = pars["gamma1"], 
#             #      alpha2 = pars["alpha2"], gamma2 = pars["gamma2"])
#   )
#   # toc(quiet = FALSE, func.toc = my.msg.toc)
#   
#   if(input_season_index != 0){
#     # take the parameter values only relevant to the season
#     lpars$N = lpars$N[input_season_index]
#     lpars$beta1 = lpars$beta1[input_season_index]
#     lpars$beta2 = lpars$beta2[input_season_index]
#     lpars$eta1 = lpars$eta1[input_season_index]
#     lpars$eta2 = lpars$eta2[input_season_index]
#   }
#   
#   # tic("generator")
#   mod = generator$new(user = lpars)
#   # toc(quiet = FALSE, func.toc = my.msg.toc)
#   # tic("execute")
#   out = mod$run(t = sort(unique(dat$iso_standard)), pretty = T)
#   # toc(quiet = FALSE, func.toc = my.msg.toc)
#   
#   # tic("tib_ModelInc_1000")
#   # EXTEND THIS TO Inc1, Inc2 and Inc12
#   # CORRESPOND WITH THE PAPER
#   stib_ModelInc_1000 <- out %>% as_tibble %>% mutate_all(as.double) %>% 
#     # select(fromHere) %>% 
#     select(t, starts_with("CumInc")) %>% 
#     mutate(across(starts_with("CumInc"), function(x) diff(c(0, x)))) %>% 
#     rename_all(~gsub("CumInc", "Inc", .x)) %>% 
#     # pivot_longer(-t, names_pattern = "([_A-Za-z]+)\\[([0-9]+)\\]", names_to = c(".value", "season_index")) %>% 
#     pivot_longer(-t, names_pattern = "([A-Za-z]+)([12])_([a-z]+)\\[([0-9]+)\\]", names_to = c(".value", "incidentVirus", "infection", "s_index")) %>% 
#     mutate(s_index = case_when(input_season_index == 0 ~ as.numeric(s_index), 
#                                T ~ input_season_index)) %>% 
#     mutate(season = season_vec[s_index]) %>% 
#     mutate(incidentVirus = as.numeric(incidentVirus)) %>% 
#     select(-s_index)
#   
#   
#   # toc(quiet = FALSE, func.toc = my.msg.toc)
#   
#   # tic("tib_ModelInc")
#   stib_ModelInc <- stib_ModelInc_1000 %>% 
#     transmute(t, season, incidentVirus, infection, ModelInc_1000 = Inc, model_CatchmentPop) %>% 
#     left_join(age_params, by = "incidentVirus", relationship = "many-to-many") %>%  #expands this file according to age_params
#     left_join(data_catchments, by = c("season", "age_group")) %>% 
#     mutate(ModelInc = ModelInc_1000/1000*CatchmentPop) %>% 
#     select(-ModelInc_1000)
#   # toc(quiet = FALSE, func.toc = my.msg.toc)
#   
#   # tic("tib_ModelHosp")
#   stib_ModelHosp <- stib_ModelInc %>% 
#     mutate(ModelHosp1 = ModelInc*case_when(infection == "simple" & incidentVirus == 1 ~ pi*h,  
#                                            infection == "coinf" & incidentVirus == 1 ~ pi*h*(1-pi), 
#                                            infection == "hitch" & incidentVirus == 1 ~ pi*h*(1-pi),
#                                            T ~ 0), 
#            ModelHosp2 = ModelInc*case_when(infection == "simple" & incidentVirus == 2 ~ pi*h, 
#                                            infection == "coinf" & incidentVirus == 2 ~ pi*h*(1-pi), 
#                                            infection == "hitch" & incidentVirus == 2 ~ pi*h*(1-pi),
#                                            T ~ 0), 
#            ModelHosp12 = ModelInc*case_when(infection == "simple" & incidentVirus == 1 ~ 0,  
#                                             infection == "coinf" & incidentVirus == 1 ~ pi*h*pi*theta, 
#                                             infection == "hitch" & incidentVirus == 1 ~ pi*h*pi,
#                                             infection == "simple" & incidentVirus == 2 ~ 0, 
#                                             infection == "coinf" & incidentVirus == 2 ~ pi*h*pi*theta, 
#                                             infection == "hitch" & incidentVirus == 2 ~ pi*h*pi,
#                                             T ~ 0)) %>% 
#     group_by(season, t, age_group) %>% 
#     summarise(across(starts_with("ModelHosp"), sum), .group = "drop") %>% 
#     ungroup %>% 
#     pivot_longer(cols = starts_with("ModelHosp"), names_prefix = "ModelHosp", names_to = "virus", values_to = "Inc") %>% 
#     mutate(virus = as.numeric(virus)) %>% 
#     select(season, t, age_group, virus, Inc) %>% 
#     arrange(season, t, match(age_group, age_vec), virus)
#   # toc(quiet = FALSE, func.toc = my.msg.toc)
#   
#   if(input_season_index == 0){
#     tib_ModelHosp = stib_ModelHosp
#   } else {
#     tib_ModelHosp <- tib_ModelHosp %>% rows_update(stib_ModelHosp, by = c("season", "t", "age_group", "virus"))
#   }
#   
#   tib_ModelHosp
# }
# 
# compute_ss <- function(pars, generator, dat_2virus, make_loglik = F){
#   tib_ModelHosp <- generate_ModelHosp(pars, generator, dat_2virus)
#   
#   
#   # # MAKE DAMN SURE WERE COMPARING THE RIGHT DATA
#   # all.equal(tib_ModelHosp %>% select(season, t, age_group, virus), dat_2virus %>% select(season, t, age_group, virus))
# 
#   loglik = sum(dpois(dat_2virus$Inc, tib_ModelHosp$Inc + 1e-25, log = T), na.rm = T)
#   
#   if(make_loglik) return(loglik)
#   -2*loglik
# }

generate_ModelHosp_season<- function(pars, generator, dat_2virus, data_catchments_weights, input_season_index = 0, tib_ModelInc_1000 = NULL){
  
  if((input_season_index == 0) != is.null(tib_ModelInc_1000)) stop("either include both input_season_index and tib_ModelHosp, or neither")
  # tic(msg = "initial variables")
  model_CatchmentPop = 1000
  age_vec = dat_2virus$age_group %>% unique
  season_vec = dat_2virus$season %>% unique
  
  # toc(quiet = FALSE, func.toc = my.msg.toc)
  
  # tic("make age params")
  age_params = tibble(incidentVirus = rep(1:2, each = length(age_vec)), 
                      age_group = rep(age_vec, times = 2), pi = rep(pars[paste("pi", sep = "_", age_vec)], times = 2), 
                      h = 10^pars[c(paste("lh1", sep = "_", age_vec), paste("lh2", sep = "_", age_vec))], 
                      theta = rep(pars[paste0("theta", 1:2)], each = length(age_vec)))
  
  age2_params = tibble(age2_group = rep(age2_vec, times = 2), 
                       incidentVirus = rep(1:2, each = length(age2_vec)), 
                       kappa = pars[c(paste("kappa1", sep = "_", age2_vec), paste("kappa2", sep = "_", age2_vec))])
  # toc(quiet = FALSE, func.toc = my.msg.toc)
  # age_params = tibble(age_group = age_vec, virus = 1, 
  #                     h = 10^pars[paste("lh1", sep = "_", age_vec)],
  #                     pi = pars[paste("pi", sep = "_", age_vec)])
  # tic("make lpars")
  
  recipR_flag = "lromega1" %in% names(pars)
  
  process_params = c("alpha", "gamma", "delta", "sigma"
                     # , "theta"
  )
  
  lpars = c(list(N = rep(model_CatchmentPop, length(season_vec))), 
            list(beta1 = pars[paste0("beta1_", gsub("/", "", season_vec))]), 
            list(beta2 = pars[paste0("beta2_", gsub("/", "", season_vec))]), 
            list(eta1 = pars[paste0("eta1_", gsub("/", "", season_vec))]), 
            list(eta2 = pars[paste0("eta2_", gsub("/", "", season_vec))]), 
            list(omega1 = ifelse(recipR_flag, 1/(10^pars["lromega1"]), pars["omega1"])),
            list(omega2 = ifelse(recipR_flag, 1/(10^pars["lromega2"]), pars["omega2"])),
            list(psi1 = ifelse(recipR_flag, 1/(10^pars["lrpsi1"]), pars["psi1"])),
            list(psi2 = ifelse(recipR_flag, 1/(10^pars["lrpsi2"]), pars["psi2"])),
            
            as.list(pars[c(paste0(process_params, 1), paste0(process_params, 2))])
            # list(alpha1 = pars["alpha1"], gamma1 = pars["gamma1"], 
            #      alpha2 = pars["alpha2"], gamma2 = pars["gamma2"])
  )
  # toc(quiet = FALSE, func.toc = my.msg.toc)
  
  if(input_season_index != 0){
    # take the parameter values only relevant to the season
    lpars$N = lpars$N[input_season_index]
    lpars$beta1 = lpars$beta1[input_season_index]
    lpars$beta2 = lpars$beta2[input_season_index]
    lpars$eta1 = lpars$eta1[input_season_index]
    lpars$eta2 = lpars$eta2[input_season_index]
  }
  
  # tic("generator")
  mod = generator$new(user = lpars)
  # toc(quiet = FALSE, func.toc = my.msg.toc)
  # tic("execute")
  # out = mod$run(t = sort(unique(dat_2virus$t)), pretty = T)
  out = mod$run(t = sort(unique(dat_2virus$t)))
  # toc(quiet = FALSE, func.toc = my.msg.toc)
  
  # tic("tib_ModelInc_1000")
  # EXTEND THIS TO Inc1, Inc2 and Inc12
  # CORRESPOND WITH THE PAPER
  stib_ModelInc_1000 <- out %>% as_tibble %>% mutate_all(as.double) %>% 
    # select(fromHere) %>% 
    select(t, starts_with("CumInc")) %>% 
    mutate(across(starts_with("CumInc"), function(x) diff(c(0, x)))) %>% 
    rename_all(~gsub("CumInc", "Inc", .x)) %>% 
    # pivot_longer(-t, names_pattern = "([_A-Za-z]+)\\[([0-9]+)\\]", names_to = c(".value", "season_index")) %>% 
    pivot_longer(-t, names_pattern = "([A-Za-z]+)([12])_([a-z]+)\\[([0-9]+)\\]", names_to = c(".value", "incidentVirus", "infection", "s_index")) %>% 
    mutate(s_index = case_when(input_season_index == 0 ~ as.numeric(s_index), 
                               T ~ input_season_index)) %>% 
    mutate(season = season_vec[s_index]) %>% 
    mutate(incidentVirus = as.numeric(incidentVirus)) %>% 
    select(-s_index)
  
  
  if(input_season_index == 0){
    tib_ModelInc_1000 = stib_ModelInc_1000
  } else {
    tib_ModelInc_1000 <- tib_ModelInc_1000 %>% rows_update(stib_ModelInc_1000, by = c("season", "t", "incidentVirus", "infection"))
  }
  
  # toc(quiet = FALSE, func.toc = my.msg.toc)
  
  # tic("tib_ModelInc")
  tib_ModelInc <- tib_ModelInc_1000 %>% 
    transmute(t, season, incidentVirus, infection, ModelInc_1000 = Inc, model_CatchmentPop) %>% 
    left_join(age_params, by = "incidentVirus", relationship = "many-to-many") %>%  #expands this file according to age_params
    left_join(data_catchments_weights %>% select(age_group, age2_group, season, CatchmentPop), by = c("season", "age_group")) %>% 
    mutate(ModelInc = ModelInc_1000/1000*CatchmentPop, 
           ModelInc_1e5 = ModelInc_1000*100) %>% 
    select(-ModelInc_1000)
  
  tib_ModelILI = tib_ModelInc %>% 
    left_join(data_catchments_weights %>% select(age_group, season, age2_group, age2_weight), by = c("age_group", "season", "age2_group")) %>% 
    left_join(age2_params, by = join_by(incidentVirus, age2_group)) %>% 
    group_by(age2_group, season, t) %>% 
    summarise(IncRate1e5 = sum(ModelInc_1e5*pi*age2_weight*kappa)) %>% 
    ungroup %>% 
    arrange(age2_group, season, t)
  # toc(quiet = FALSE, func.toc = my.msg.toc)
  
  # tic("tib_ModelHosp")
  tib_ModelHosp <- tib_ModelInc %>% 
    mutate(ModelHosp1 = ModelInc*case_when(infection == "simple" & incidentVirus == 1 ~ pi*h,  
                                           infection == "coinf" & incidentVirus == 1 ~ pi*h*(1-pi), 
                                           infection == "hitch" & incidentVirus == 1 ~ pi*h*(1-pi),
                                           T ~ 0), 
           ModelHosp2 = ModelInc*case_when(infection == "simple" & incidentVirus == 2 ~ pi*h, 
                                           infection == "coinf" & incidentVirus == 2 ~ pi*h*(1-pi), 
                                           infection == "hitch" & incidentVirus == 2 ~ pi*h*(1-pi),
                                           T ~ 0), 
           ModelHosp12 = ModelInc*case_when(infection == "simple" & incidentVirus == 1 ~ 0,  
                                            infection == "coinf" & incidentVirus == 1 ~ pi*h*pi*theta, 
                                            infection == "hitch" & incidentVirus == 1 ~ pi*h*pi,
                                            infection == "simple" & incidentVirus == 2 ~ 0, 
                                            infection == "coinf" & incidentVirus == 2 ~ pi*h*pi*theta, 
                                            infection == "hitch" & incidentVirus == 2 ~ pi*h*pi,
                                            T ~ 0)) %>% 
    group_by(season, t, age_group) %>% 
    summarise(across(starts_with("ModelHosp"), sum), .group = "drop") %>% 
    ungroup %>% 
    pivot_longer(cols = starts_with("ModelHosp"), names_prefix = "ModelHosp", names_to = "virus", values_to = "Inc") %>% 
    mutate(virus = as.numeric(virus)) %>% 
    select(season, t, age_group, virus, Inc) %>% 
    arrange(season, t, match(age_group, age_vec), virus)
  # toc(quiet = FALSE, func.toc = my.msg.toc)
  
  # the hosp part is used to calculate likelihood
  # the carry part can be carried over so that only a single season needs to be numerically integrated, if all other season-specific parameters remain fixed
  return_object = list(hosp = tib_ModelHosp, ili = tib_ModelILI, carry = tib_ModelInc_1000)
  return_object
}

get_loglik_hosp <- function(dat_2virus, tib_ModelHosp, return_vector = F){
  if(nrow(dat_2virus) != nrow(tib_ModelHosp)) stop("hospital data and model are different lengths")
  loglik_vector = dpois(dat_2virus$Inc, tib_ModelHosp$Inc + 1e-10, log = T)
  if(return_vector){
    return(loglik_vector)
  } else {
    return(sum(loglik_vector, na.rm = T))
  }
  # sum(dpois(dat_2virus$Inc, tib_ModelHosp$Inc + 1e-10, log = T), na.rm = T)
}
get_loglik_ili <- function(ili, tib_ModelILI, return_vector = F){
  if(nrow(ili) != nrow(tib_ModelILI)) stop("ILI data and model are different lengths")
  
  epsilon = 1
  loglik_vector = dlnorm(x = ili$IncRate1e5 + epsilon, meanlog = log(tib_ModelILI$IncRate1e5 + epsilon), sdlog = ((ili$IncRate1e5 + epsilon) %>% log %>% sd(na.rm = T)), log = T)
  if(return_vector){
    return(loglik_vector)
  } else {
    return(sum(loglik_vector, na.rm = T))
  }
  # sum(dlnorm(x = ili$IncRate1e5 + epsilon, meanlog = log(tib_ModelILI$IncRate1e5 + epsilon), sdlog = ((ili$IncRate1e5 + epsilon) %>% log %>% sd(na.rm = T)), log = T), na.rm = T)
}


# compute_ss <- function(pars, generator, dat_2virus, ili, include_ili = T, data_catchments_weights, make_loglik = F){
#   model_out <- generate_ModelHosp_season(pars, generator, dat_2virus, data_catchments_weights)
#   tib_ModelHosp <- model_out$hosp
#   tib_ModelILI <- model_out$ili
#   # tib_ModelHosp <- generate_ModelHosp_season(pars, generator, dat_2virus)$hosp
#   
#   # # MAKE DAMN SURE WERE COMPARING THE RIGHT DATA
#   # all.equal(tib_ModelHosp %>% select(season, t, age_group, virus), dat_2virus %>% select(season, t, age_group, virus))
#   
#   loglik_hosp = get_loglik_hosp(dat_2virus, tib_ModelHosp)
#   
#   # loglik_hosp = sum(dpois(dat_2virus$Inc, tib_ModelHosp$Inc + 1e-25, log = T), na.rm = T)
#   # loglik_ili = sum(dpois(ili$IncRate1e5*1e6, lambda = (tib_ModelILI$IncRate1e5 + 1e-25)*1e6, log = T), na.rm = T)
#   if(include_ili){
#     loglik_ili = get_loglik_ili(ili, tib_ModelILI)
#     loglik = loglik_hosp + loglik_ili
#   } else {
#     loglik = loglik_hosp
#   }
#   
#   if(make_loglik) return(loglik)
#   -2*loglik
# }

# pars = best_params %>% select(-Loglik) %>% unlist
# generator = m5_generator
# dat_2vi
compute_ss <- function(pars, generator, dat_2virus, ili, include_ili = T, data_catchments_weights, make_loglik = F, return_solution = F){
  model_out <- generate_ModelHosp_season(pars, generator, dat_2virus, data_catchments_weights)
  tib_ModelHosp <- model_out$hosp
  tib_ModelILI <- model_out$ili
  
  if(nrow(tib_ModelHosp) != nrow(dat_2virus)) stop("Hospital data and model are not the same length")
  if(nrow(tib_ModelILI) != nrow(ili)) stop("ILI data and model are not the same length")
  # tib_ModelHosp <- generate_ModelHosp_season(pars, generator, dat_2virus)$hosp
  
  # # MAKE DAMN SURE WERE COMPARING THE RIGHT DATA
  # all.equal(tib_ModelHosp %>% select(season, t, age_group, virus), dat_2virus %>% select(season, t, age_group, virus))
  if(return_solution){
    loglik_hosp_vec = get_loglik_hosp(dat_2virus, tib_ModelHosp, return_vector = T)
    tib_ModelHosp$Loglik = loglik_hosp_vec
    
    if(include_ili){
      loglik_ili_vec = get_loglik_ili(ili, tib_ModelILI, return_vector = T)
      tib_ModelILI$Loglik = loglik_ili_vec
      return(list(hosp = tib_ModelHosp, ili = tib_ModelILI))
    }
    
    return(tib_ModelHosp)
    
  } else {
    loglik_hosp = get_loglik_hosp(dat_2virus, tib_ModelHosp)
    if(include_ili){
      loglik_ili = get_loglik_ili(ili = ili, tib_ModelILI = tib_ModelILI)
      loglik = loglik_hosp + loglik_ili
    } else {
      loglik = loglik_hosp
    }
    
    if(make_loglik){
      return(loglik)
    } else {
      return(-2*loglik)
    }
  }
}


plot_model_data <- function(pars, generator, dat_2virus, ili = NULL, virus1, virus2, type = "ageseason", data_catchments_weights, include_ci = F){
  
  model_out <- generate_ModelHosp_season(pars, generator, dat_2virus, data_catchments_weights = data_catchments_weights)
  
  if(grepl("ili", type)){
    tib_ModelILI <- model_out$ili
    if(is.null(ili)) stop("no ILI data provided.")
    
    plot_data <- ili %>% 
      transmute(season, age2_group = factor(age2_group, levels = age2_vec), t, 
                DataILI = IncRate1e5, 
                ModelILI = tib_ModelILI$IncRate1e5)
    
  } else  {
    tib_ModelHosp <- model_out$hosp
    
    plot_data <- dat_2virus %>% 
      transmute(season, age_group = factor(age_group, levels = age_vec), t, 
                virus = c(virus1, virus2, "Codetection")[match(virus, c(1, 2, 12))], 
                DataHosp = Inc, 
                ModelHosp = tib_ModelHosp$Inc)
  }
  
  
  # # MAKE DAMN SURE WERE COMPARING THE RIGHT DATA
  # all.equal(tib_ModelHosp %>% select(season, t, age_group, virus), dat_2virus %>% select(season, t, age_group, virus))
  

  censusweights = data_catchments_weights %>% 
    group_by(season, age2_group) %>% 
    summarise(CensusPopulation = sum(CensusPopulation)) %>% 
    group_by(season) %>% 
    mutate(CensusWeight = CensusPopulation/sum(CensusPopulation)) %>% 
    ungroup
  
  
  if(type == "ili_ageseason"){
    plot_obj <- plot_data %>% 
      mutate(model_loCI = qpois(0.025, ModelILI, lower.tail = T), 
             model_hiCI = qpois(0.025, ModelILI, lower.tail = F)) %>%  
      ggplot(aes(x = t)) + 
      {if(include_ci) geom_ribbon(aes(ymin = model_loCI, ymax = model_hiCI, fill = "CI"))} + 
      geom_line(aes(y = ModelILI,  colour = "Model")) + 
      geom_point(aes(y = DataILI, colour = "Data")) + 
      facet_grid(age2_group~season, scales = "free_y") + 
      labs(x = "Week", y = "ILI incidence rate", colour = "", linetype = "", shape = "") + 
      theme_bw() + 
      scale_fill_manual(values = c(CI = "grey"))
  } else if (type == "ili_season"){
    plot_obj <- plot_data %>% 
      left_join(censusweights) %>% 
      group_by(season, t) %>% 
      summarise(across(ends_with("ILI"), ~sum(.x*CensusWeight, na.rm = T))) %>% 
      mutate(model_loCI = qpois(0.025, ModelILI, lower.tail = T), 
             model_hiCI = qpois(0.025, ModelILI, lower.tail = F)) %>%  
      ggplot(aes(x = t)) + 
      {if(include_ci) geom_ribbon(aes(ymin = model_loCI, ymax = model_hiCI, fill = "CI"))} + 
      geom_line(aes(y = ModelILI,  colour = "Model")) + 
      geom_point(aes(y = DataILI, colour = "Data")) + 
      facet_grid(.~season, scales = "free_y") + 
      labs(x = "Week", y = "ILI incidence rate", colour = "", linetype = "", shape = "") + 
      theme_bw() + 
      scale_fill_manual(values = c(CI = "grey"))
  } else if (type == "ili_age"){
    plot_obj <- plot_data %>% 
      group_by(age2_group, t) %>% 
      summarise(across(ends_with("ILI"), ~mean(.x, na.rm = T))) %>% 
      mutate(model_loCI = qpois(0.025, ModelILI, lower.tail = T), 
             model_hiCI = qpois(0.025, ModelILI, lower.tail = F)) %>%  
      ggplot(aes(x = t)) + 
      {if(include_ci) geom_ribbon(aes(ymin = model_loCI, ymax = model_hiCI, fill = "CI"))} + 
      geom_line(aes(y = ModelILI,  colour = "Model")) + 
      geom_point(aes(y = DataILI, colour = "Data")) + 
      facet_grid(.~age2_group, scales = "free_y") + 
      labs(x = "Week", y = "ILI incidence rate", colour = "", linetype = "", shape = "") + 
      theme_bw() + 
      scale_fill_manual(values = c(CI = "grey"))
  } else if(type == "ageseason"){
    plot_obj <- plot_data %>% 
      mutate(model_loCI = qpois(0.025, ModelHosp, lower.tail = T), 
             model_hiCI = qpois(0.025, ModelHosp, lower.tail = F)) %>%  
      ggplot(aes(x = t)) + 
      {if(include_ci) geom_ribbon(aes(ymin = model_loCI, ymax = model_hiCI, fill = "CI"))} + 
      
      geom_line(aes(y = ModelHosp,  linetype = "Model", colour = virus)) + 
      geom_point(aes(y = DataHosp, shape = "Data", colour = virus)) +
      facet_grid(age_group~season, scales = "free_y") + 
      labs(x = "Week", y = "Cases", colour = "", linetype = "", shape = "") + 
      theme_bw() + 
      scale_fill_manual(values = c(CI = "grey"))
    
  } else if (type == "season") {
    
    plot_obj <- plot_data %>% 
      group_by(season, t, virus) %>% 
      summarise(across(ends_with("Hosp"), ~sum(.x, na.rm = T))) %>% 
      mutate(model_loCI = qpois(0.025, ModelHosp, lower.tail = T), 
             model_hiCI = qpois(0.025, ModelHosp, lower.tail = F)) %>%  
      ggplot(aes(x = t)) + 
      {if(include_ci) geom_ribbon(aes(ymin = model_loCI, ymax = model_hiCI, fill = "CI"))} + 
      
      geom_line(aes(y = ModelHosp,  linetype = "Model", colour = virus)) + 
      geom_point(aes(y = DataHosp, shape = "Data", colour = virus)) +
      facet_grid(virus~season, scales = "free_y") + 
      labs(x = "Week", y = "Cases", colour = "", linetype = "", shape = "") + 
      theme_bw()+ 
      scale_fill_manual(values = c(CI = "grey"))
    
    
  } else if (type == "age"){
    plot_obj <- plot_data %>% 
      group_by(age_group, t, virus) %>% 
      summarise(across(ends_with("Hosp"), ~sum(.x, na.rm = T))) %>% 
      mutate(model_loCI = qpois(0.025, ModelHosp, lower.tail = T), 
             model_hiCI = qpois(0.025, ModelHosp, lower.tail = F)) %>%  
      ggplot(aes(x = t)) + 
      {if(include_ci) geom_ribbon(aes(ymin = model_loCI, ymax = model_hiCI, fill = "CI"))} + 
      
      geom_line(aes(y = ModelHosp,  linetype = "Model", colour = virus)) + 
      geom_point(aes(y = DataHosp, shape = "Data", colour = virus)) + 
      facet_grid(virus~age_group, scales = "free_y") + 
      labs(x = "Week", y = "Cases", colour = "", linetype = "", shape = "") + 
      theme_bw()+ 
      scale_fill_manual(values = c(CI = "grey"))
    
    
  } else {
    print("unknown type")
  }
  
  
  return(plot_obj)
}

# stem = "~/SanOdin/output/"
# id = "odin1v_"
# virus = "Influenza"
# exclude_params = c("alpha", "gamma", paste0("pi_", age_vec))
# burnin = 0
# thin = 1

# stem = "output/"
# id = "simpleInter_"
# virus1 = "Influenza"
# virus2 = "RSV"
# burnin = 0
# thin = 1
# exclude_params = c(paste0(rep(c("alpha", "gamma", "omega", "psi", "delta", "sigma", "theta"), each = 2), 1:2), paste0("pi_", age_vec))
plot_post_2v <- function(stem, id, virus1, virus2, jump = NULL, rep = NULL, exclude_params = NULL, thin = 1, burnin = 0, save = F, w = 20, h = 12){
  
  
  post = read_csv(paste0(stem, id, virus1, "+", virus2, 
                         ifelse(is.null(jump), "", paste0("_jump", jump)), 
                         ifelse(is.null(rep), "", paste0("_rep", rep)), 
                         "_posterior.csv"))
  
  if(burnin > max(post$Iteration)){stop("burnin too long")}
  
  post_cropped <- post %>% 
    filter(Iteration >= burnin) %>% 
    filter((Iteration %% thin) == 0) %>% 
    select(-exclude_params, -ends_with("ccepted"), -ends_with("eason"))
  
  post_cropped %>% 
    select(-Iteration, -ends_with("oglik")) %>% 
    summarise_all(effectiveSize) %>% 
    t %>% 
    print
  
  post_long <- post_cropped %>% 
    pivot_longer(-c("Iteration"), names_sep = "_", names_to = c("Parameter", "category")) %>% 
    mutate(category = replace_na(category, "common"), 
           age_group = ifelse(!grepl("^[0-9]{4}$", category) & (category != "common"), category, "all"), 
           season = ifelse(grepl("^[0-9]{4}$", category), category, "common")) %>% 
    select(-category)
  
  post_plotobj <- post_long %>% 
    ggplot(aes(x = Iteration, y = value, alpha = season, colour = age_group
               )) + 
    geom_line() + 
    facet_wrap(.~Parameter, scales = "free_y") + 
    theme_bw() + 
    labs(title = paste0(virus1, "+", virus2), x = "Iterations")
  
  print(post_plotobj)
  
  if(save){
    ggsave(post_plotobj, filename = paste0(stem, id, virus, "_trace.png"), width = w, height = h, units = "cm")
  }
}

plot_acceptance_2v <- function(stem, id, virus1, virus2, jump = NULL, rep = NULL){
  
  
  post = read_csv(paste0(stem, id, virus1, "+", virus2, 
                         ifelse(is.null(jump), "", paste0("_jump", jump)), 
                         ifelse(is.null(rep), "", paste0("_rep", rep)), 
                         "_posterior.csv"))
  acc_plot <- post %>% 
    filter(!is.na(Accepted)) %>% 
    mutate(CumAccepted = cumsum(Accepted)/(1:n())) %>% 
    ggplot(aes(x = Iteration, y = CumAccepted)) + geom_line()
  
  print(acc_plot)
}

compare_reps <- function(stem = "output/", id = id, virus1 = virus1, virus2 = virus2, 
                         jump = 0, reps = 0:5, burnin = 0, thin = 1, 
                         exclude_params = c(paste0(rep(c("alpha", "gamma", "omega", "psi", "delta", "sigma", "theta"), each = 2), 1:2), paste0("pi_", age_vec))
                         # exclude_params = c(paste0(rep(c("alpha", "gamma", "psi", "delta", "sigma", "theta"), each = 2), 1:2), paste0("pi_", age_vec))
                         # exclude_params = c("alpha1", "gamma1", "alpha2", "gamma2", paste0("pi_", age_vec))
){
  library(ggstance)
  library(ggpubr)
  
  post = NULL
  for(rep in reps){
    post_tmp = read_csv(paste0(stem, id, virus1, "+", virus2, 
                               ifelse(is.null(jump), "", paste0("_jump", jump)), 
                               ifelse(is.null(rep), "", paste0("_rep", rep)), 
                               "_posterior.csv")) %>% 
      mutate(rep = rep)
    
    if(is.null(post)){
      post <- post_tmp
    } else {
      post <- rbind(post, post_tmp)
    }
  }
  
  if(burnin > max(post$Iteration)){stop("burnin too long")}
  
  post_cropped <- post %>% 
    filter(Iteration >= burnin) %>% 
    filter((Iteration %% thin) == 0) %>% 
    select(-exclude_params, -ends_with("ccepted"), -ends_with("eason"))
  
  post_long <- post_cropped %>% 
    pivot_longer(-c("Iteration", "rep", "Loglik"), 
                 names_pattern = "([a-z]+)([12])([_]?.*)", 
                 names_to = c("Parameter", "virus", "category")) %>% 
    mutate(category = gsub("_", "", category)) %>% 
    mutate(season = ifelse(grepl("^[0-9]{4}$", category), category, "common"), 
           age_group = ifelse((season == "common") & (!(category == "")), category, "all")) %>% 
    select(-category)
  
  
  post_long %>% 
    group_by(rep) %>% 
    summarise(Loglik = max(Loglik, na.rm = T)) %>% 
    print
  
  K <- post_long %>% 
    select(Parameter, virus, age_group, season) %>% 
    unique %>% nrow
  
  post_ci <- post_long %>% 
    filter(!is.na(Loglik)) %>% 
    group_by(rep, Parameter, virus, age_group, season) %>% 
    mutate(maxLoglik = max(Loglik, na.rm = T), 
           deltaLoglik = Loglik - maxLoglik, 
           LRTstat = -2*deltaLoglik, 
           LRTp = pchisq(q = LRTstat, df = K, lower.tail = F), 
           inCI = (0.05 <= LRTp)) %>% 
    summarise(best = median(value[Loglik == maxLoglik]), 
              loCI = min(value[inCI]), 
              hiCI = max(value[inCI])) %>% 
    mutate(rep = as.factor(rep))
  
  plot_kappa <- post_ci %>% 
    filter(Parameter == "kappa") %>% 
    ggplot(aes(y = age_group, colour = rep, group = rep)) + 
    geom_point(aes(x = best), position = position_dodgev(height = 0.5)) + 
    geom_errorbarh(aes(xmin = loCI, xmax = hiCI), position = position_dodgev(height = 0.5)) + 
    facet_grid(virus~Parameter, scales = "free_x") +
    theme_bw()
  
  plot_lh <- post_ci %>% 
    filter(Parameter == "lh") %>% 
    mutate(Parameter == "h") %>% 
    ggplot(aes(y = age_group, colour = rep, group = rep)) + 
    geom_point(aes(x = 10^(best)), position = position_dodgev(height = 0.5)) + 
    geom_errorbarh(aes(xmin = 10^(loCI), xmax = 10^(hiCI)), position = position_dodgev(height = 0.5)) + 
    scale_x_log10() + 
    facet_grid(virus~Parameter, scales = "free_x") +
    theme_bw()
  
  plot_betaeta <- post_ci %>% 
    filter(age_group == "all") %>%
    filter(Parameter %in% c("beta", "eta")) %>% 
    ggplot(aes(y = season, colour = rep, group = rep)) + 
    geom_point(aes(x = best), position = position_dodgev(height = 0.5)) + 
    geom_errorbarh(aes(xmin = loCI, xmax = hiCI), position = position_dodgev(height = 0.5)) + 
    facet_grid(virus~Parameter, scales = "free_x")  +
    theme_bw()
  
  
  plot_interaction <- post_ci %>% 
    filter(age_group == "all" & season == "common") %>% 
    mutate(process = ifelse(Parameter %in% c("omega", "psi"), "Process", "Modifier")) %>% 
    ggplot(aes(y = Parameter, colour = rep, group = rep)) + 
    geom_point(aes(x = best), position = position_dodgev(height = 0.5)) + 
    geom_errorbarh(aes(xmin = loCI, xmax = hiCI), position = position_dodgev(height = 0.5)) + 
    facet_grid(virus~process, scales = "free_x")  +
    theme_bw()
  
  plot_reps <- ggarrange(ggarrange(plot_lh, plot_kappa, nrow = 1), 
                         plot_betaeta, 
                         plot_interaction, nrow = 3)
  
  print(plot_reps)
}



param_estimates <- function(stem, id, virus){
  post = read_csv(paste0(stem, id, virus, "_posterior.csv"))
  
  post_long <- post %>% 
    select(-alpha, -gamma) %>% 
    pivot_longer(-c("Iteration", "Loglik"), names_sep = "_", names_to = c("Parameter", "category")) %>% 
    mutate(category = replace_na(category, "common"), 
           age_group = ifelse(!grepl("^[0-9]{4}$", category) & (category != "common"), category, "all"), 
           season = ifelse(grepl("^[0-9]{4}$", category), category, "common")) %>% 
    select(-category)
  
  K = post %>% select(starts_with("beta"), starts_with("eta"), starts_with("lh")) %>% ncol
  
  post_ci <- post_long %>% 
    mutate(age_group = factor(age_group, levels = c("all", age_vec))) %>% 
    group_by(age_group, Parameter, season) %>% 
    mutate(maxLoglik = max(Loglik), 
           deltaLoglik = Loglik - max(Loglik), 
           LRTstat = -2*deltaLoglik, 
           LRTp = pchisq(q = LRTstat, df = K, lower.tail = F), 
           inCI = (0.05 <= LRTp))  %>% 
    # filter(inCI, Parameter == "omega1") %>% pull(value)
    # filter(!in_ci) %>% View
    summarise(across(value, list(#median = median, 
      best = ~median(.x[Loglik == maxLoglik]), 
      # lo = ~quantile(.x, 0.025), 
      # hi = ~quantile(.x, 0.975), 
      lo = ~min(.x[inCI]), 
      hi = ~max(.x[inCI]), 
      ESS = ~effectiveSize(.x)), .names = "{.fn}"), 
      bestLoglik = max(Loglik), 
      K = K[1])
  
  post_ci
}

week_to_date_num = function(week) {
  ymd = as.Date("0000-12-31") + week*7
  gsub("000[01]-", "", as.character(ymd))
} 

week_to_date_name = function(week) {
  Sys.setlocale(locale = "en")
  ymd = as.Date("0000-12-31") + week*7
  ymd = format(ymd, "%Y-%b-%d")
  gsub("-[0]?", " ", gsub("000[01]-", "", as.character(ymd)))
} 


week_to_date = function(week) {
  if(!is.numeric(week)) return(week)
  Sys.setlocale(locale = "en")
  ymd = as.Date("0000-12-31") + week*7
  ymd
} 

age2_age <- function(age){
  age2_vec[c(1, 1, 1, 1, 2, 3, 3, 4, 4)[match(age, age_vec)]]
}
