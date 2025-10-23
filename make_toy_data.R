coinc <- read_tsv("~/ViralInteraction/data/allcatch_coinc_age.tsv")
coinc %>% 
  filter(CoVirus == "Influenza+RSV") %>% 
  group_by(age_group, season) %>% 
  filter(iso_standard %in% c(min(iso_standard), max(iso_standard), 0)) %>% 
  write_tsv("~/ViralInteraction/data_public/toydata_allcatch_coinc_age.tsv")


inc <- read_tsv("~/ViralInteraction/data/allcatch_inc_age.tsv")
inc %>% 
  filter(virus %in% c("Influenza", "RSV")) %>% 
  group_by(age_group, season) %>% 
  filter(iso_standard %in% c(min(iso_standard), max(iso_standard), 0)) %>% 
  write_tsv("~/ViralInteraction/data_public/toydata_allcatch_inc_age.tsv")
