get_best_rep <- function(stem = "output/", id, virus1, virus2, 
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
                                   "_posterior.csv"), show_col_types = FALSE) %>% 
        mutate(rep = r)
      
      # print(post_piece)
      
      if(is.null(post)){
        post = post_piece
      } else {
        post = rbind(post, post_piece)
      }
    }
  }
  

  post %>% 
    group_by(rep) %>% 
    summarise(maxLoglik = max(Loglik, na.rm = T)) %>% 
    return
}




# jump_rep_performance = NULL
for(j in 0:10){
  print(j)
  r_vec = 0:5
  
  if(j == 0){
    this_id = "simpleInter_"
  } else {
    this_id = "validateSimpleInter_"
  }
  this_burnin = 0
  
  # aic_piece = param_table(stem = "output/", id = this_id, virus1 = "Influenza", virus2 = "RSV", 
  #                         this_jump = j, this_rep = r_vec, burnin = this_burnin, thin = 1, param_jumps)
  
  rep_performance = get_best_rep(stem = "output/", id = this_id, virus1 = "Influenza", virus2 = "RSV", 
                          this_jump = j, this_rep = r_vec, burnin = this_burnin, thin = 1, param_jumps)
  
  best_rep = rep_performance %>% arrange(-maxLoglik) %>% slice(1) %>% pull(rep)
  
  print(c(j, best_rep))
  
  
  best_params <- read_csv(paste0("output/", id, virus1, "+", virus2, "_jump", 3, "_rep", best_rep, "_bestParams.csv"), show_col_types = F)
  post <- read_csv(paste0("output/", id, virus1, "+", virus2, "_jump", 3, "_rep", best_rep, "_posterior.csv"), show_col_types = F)
  
  best_ModelHosp <- generate_ModelHosp_season(pars = best_params %>%
                                                select(-Loglik) %>%
                                                unlist, generator = m5_generator, dat_2virus = dat_2virus,
                                              data_catchments_weights = data_catchments_weights
  )$hosp %>% 
    group_by(season, t, virus) %>% 
    summarise(Inc = sum(Inc))
  
  # make a subsample
  N = 20
  set.seed(1)
  sample_params = post[sample(x = nrow(post), size = N), ] %>% 
    select(-Season)
  
  sample_ModelHosp <- map(split(sample_params, 1:nrow(sample_params)), 
                          #split(sample_params[1:3, ], 1:3), 
                          .f = function(p) generate_ModelHosp_season(pars = unlist(p), generator = m5_generator, dat_2virus = dat_2virus,
                                                                     data_catchments_weights = data_catchments_weights)$hosp %>% 
                            group_by(season, t, virus) %>% 
                            summarise(Inc = sum(Inc)), 
                          .progress = list(
                            type = "iterator", 
                            format = "Calculating {cli::pb_bar} {cli::pb_percent}",
                            clear = TRUE)) %>% 
    data.table::rbindlist(idcol = "rep") %>% 
    as_tibble
  
  
  sample_ModelHosp_LoHi <- sample_ModelHosp %>% 
    group_by(season, t, virus) %>% 
    summarise(LoInc = quantile(Inc, 0.025), HiInc = quantile(Inc, 0.975))
  
  
  incdiff <- best_ModelHosp %>% 
    # left_join(nointer_ModelHosp) %>% 
    left_join(datainc, by = c("season", "t", "virus")) %>% 
    left_join(catchment_season, by = c("season")) %>%
    left_join(sample_ModelHosp_LoHi, by = c("season", "t", "virus")) %>% 
    # left_join(dat_2virus %>% 
    #             filter(t == 0, virus == 1) %>% 
    #             group_by(season) %>% 
    #             summarise(CatchmentPop = sum(CatchmentPop))) %>% 
    mutate(across(ends_with("Inc"), ~.x/CatchmentPop*1e5)) %>%
    # mutate(IncDiff = (Inc - BaselineInc)) %>% 
    mutate(virus =fct_recode(factor(virus, levels = c(1, 2, 12)), `Influenza only`="1", `RSV only`="2", `Codetection`="12"))
  
  # calculate the proportion of the data points fall inside the 95% credibility interval
  inblue_piece = incdiff %>% 
    ungroup %>% 
    mutate(inblue = ((DataInc < LoInc) | (DataInc > HiInc))) %>% 
    group_by(virus) %>% 
    summarise(codetection_inblue = mean(inblue)) %>% 
    mutate(jump = j)
  
  if(j == 0){
    inblue = inblue_piece
  } else {
    inblue = rbind(inblue, inblue_piece)
  }
  
  print(inblue %>% 
          filter(virus == "Codetection"), n = 100)
  
  # if(j == 0){
  #   jump_rep_performance = rep_performance %>% mutate(jump = j)
  # } else {
  #   jump_rep_performance = rbind(jump_rep_performance, 
  #                                rep_performance %>% mutate(jump = j))
  # }
  
}
print(inblue, n = 100)
inblue %>% 
  filter(virus == "Codetection") %>% 
  left_join(param_description %>% select(-virus) %>% slice(1:11))
