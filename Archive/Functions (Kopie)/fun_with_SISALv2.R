####-------------------- FUN WITH SISAL v2 ---------------------------####

# load libraries
library(tidyverse)

### load and transform data --------------------------------------###
# load SISAL v2 csv files

load_sisal_data <- function(prefix = "", path = "/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/02_SISAL/SISAL_v2_CARLA/"){
  prefix = ""
  path = "/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/02_SISAL/SISAL_v2_CARLA/"
  composite_link_entity <- read.csv(paste(path, prefix,'composite_link_entity.csv',sep = ''), header = T,stringsAsFactors = F)
  d13C <- read.csv(paste(path, prefix,'d13C.csv',sep='') ,header = T, stringsAsFactors = F)
  d13C <- plyr::rename(d13C, c("iso_std" = "iso_std_d13C"))
  d18O <- read.csv(paste(path, prefix,'d18O.csv', sep =''),header = T, stringsAsFactors = F)
  d18O <- plyr::rename(d18O, c("iso_std" = "iso_std_d18O"))
  dating_lamina <- read.csv(paste(path, prefix,'dating_lamina.csv', sep = ''), header = T, stringsAsFactors = F)
  dating <- read.csv(paste(path, prefix,'dating.csv',sep = ''), header = T, stringsAsFactors = F)
  entity_link_reference <- read.csv(paste(path, prefix,'entity_link_reference.csv', sep = ''), header =T, stringsAsFactors = F)
  entity <- read.csv(paste(path, prefix,'entity.csv', sep = ''), header = T, stringsAsFactors = F)
  gap <- read.csv(paste(path, prefix,'gap.csv', sep = ''), header = T, stringsAsFactors = F)
  hiatus <- read.csv(paste(path, prefix,'hiatus.csv', sep =''), header = T, stringsAsFactors = F)
  notes <- read.csv(paste(path, prefix,'notes.csv', sep = ''), header = T, stringsAsFactors = F)
  original_chronology <- read.csv(paste(path, prefix,'original_chronology.csv', sep = ''), header = T, stringsAsFactors = F)
  reference <- read.csv(paste(path, prefix,'reference.csv', sep = ''), header = T, stringsAsFactors = F)
  sample <- read.csv(paste(path, prefix,'sample.csv', sep = ''), header = T, stringsAsFactors = F)
  sisal_chronology <- read.csv(paste(path, prefix,'sisal_chronology.csv', sep = ''), header = T, stringsAsFactors = F)
  site <- read.csv(paste(path, prefix,'site_countries.csv', sep = ''), header = T, stringsAsFactors = F)


  # build SISAL tables
  site_tb <- left_join(site, entity, by = 'site_id') %>% left_join(., entity_link_reference, by = 'entity_id') %>% 
    left_join(., reference, by = 'ref_id') %>% left_join(., notes, by = 'site_id') %>% mutate_at(vars(site_id, entity_id), as.numeric)
  dating_tb <- left_join(dating, entity) %>% group_by(entity_id) %>%mutate(laminar_dated = if_else((entity_id %in% dating_lamina$entity_id), 'yes', 'no')) %>% 
    mutate_at(vars(dating_id, depth_dating, dating_thickness, X14C_correction, corr_age, corr_age_uncert_pos, corr_age_uncert_neg), as.numeric) %>%ungroup()
  sample_tb <- plyr::join_all(list(sample,hiatus, gap, original_chronology, sisal_chronology, d13C, d18O), by = 'sample_id', type = 'full', match = 'all') %>% 
    mutate_at(vars(entity_id, sample_id, sample_thickness, depth_sample, interp_age, interp_age_uncert_pos, interp_age_uncert_neg, COPRA_age,
                   COPRA_age_uncert_pos, COPRA_age_uncert_neg, linear_age, linear_age_uncert_pos, linear_age_uncert_neg, d13C_measurement,
                   d13C_precision, d18O_measurement, d18O_precision), as.numeric)


  # filter for 'from base' dated entities
  entity_from_base <- site_tb %>% filter(depth_ref == 'from base') %>% distinct(entity_id)
  sample_from_base <- sample_tb %>% filter(entity_id %in% entity_from_base$entity_id) %>% 
    select(entity_id,depth_sample) %>% group_by(entity_id) %>% dplyr::summarise(max = max(depth_sample))

  # transform depths for 'from base' dated entities in dating file
  dating_tb_new <- full_join(dating_tb, sample_from_base, by = 'entity_id') %>% group_by(entity_id) %>% 
    mutate(depth_conv = if_else(entity_id %in% entity_from_base$entity_id, max-depth_dating, NA_real_)) %>% 
    mutate(depth_dating = if_else(!is.na(depth_conv), depth_conv, depth_dating)) %>%
    select(-depth_conv) %>% arrange(., depth_dating, .by_group = T)

  #transform depths for 'from base' dated entities in sample file
  sample_tb_new <- full_join(sample_tb, sample_from_base, by = 'entity_id') %>% group_by(entity_id) %>% 
    mutate(depth_conv = if_else(entity_id %in% entity_from_base$entity_id, max-depth_sample, NA_real_)) %>% 
    mutate(depth_sample = if_else(!is.na(depth_conv), depth_conv, depth_sample)) %>%
    select(-depth_conv) %>% arrange(., depth_sample, .by_group = T)

  ### entity objects ---------------------------------------------###
  # In the MCensemble file each object name is designed as follows: entity_id - entity_name
  # Each object is a named list containing information (site information, references, notes), the proxy values, the original chronology,
  # the available dating information, hiatus depths and the chronologies and ensembles for the 6 different AM.

  # load entity object and assign name for, e.g. entity_id = 1
  eID_1 <- get(load('/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/02_SISAL/SISAL_v2_CARLA/1-BT-1'))


  ### filter data ------------------------------------------------###

  ## filter for climate periods -----------------------------##
  # filter dating table
  dating_tb_filtered <- dating_tb_new %>% filter(entity_status == 'current') %>% mutate_at(vars(corr_age),as.numeric) %>% 
    filter(date_used == 'yes' & date_type != 'Event; hiatus')

  # define periods and determine age ranges for each period
  p1 <- dating_tb %>% group_by(entity_id) %>% filter(corr_age <= 10000) %>% 
    dplyr::summarise(diff1 = max(corr_age, na.rm = T) - min(corr_age, na.rm = T),n1 = n(), z1 = 1) %>% 
    filter(n1 > 2 & diff1 >= 4000)          # filter for at least 3 dates and a range of at least 4000 years
  p2 <- dating_tb %>% group_by(entity_id) %>% filter(corr_age <= 20000 & corr_age > 10000) %>% 
    dplyr::summarise(diff2 = max(corr_age, na.rm = T) - min(corr_age, na.rm = T), n2 = n(), z2 = 1) %>% 
    filter(n2 > 2 & diff2 >= 4000)
  p3 <- dating_tb %>% group_by(entity_id) %>% filter(corr_age <= 30000 & corr_age > 20000) %>% 
    dplyr::summarise(diff3 = max(corr_age, na.rm = T) - min(corr_age, na.rm = T),n3 = n(), z3 = 1) %>% 
    filter(n3 > 2 & diff3 >= 4000)

  # join climate period tables p1,p2 and p3 and filter for z > 1 to get those entity ids covering more than just one period with a high enough resolution
  p <- plyr::join_all(list(p1,p2,p3), by = 'entity_id', type = 'full') %>% group_by(entity_id) %>% 
    mutate(z = sum(c(z1,z2,z3), na.rm = T)) %>% 
    filter(z > 1) %>% 
    arrange(entity_id)

  # retrieve long and lat for the filtered entty ids
  p_sites <- site_tb %>% filter(entity_id %in% p$entity_id) %>% select(site_id, site_name, entity_id, longitude, latitude)



  ## filter for entities with certain dating properties -----------------------------##
  # entities with no published AM
  no_pub_AM <- left_join(sisal_chronology, sample) %>% distinct(entity_id)

  # entities containing uncertainties in original chronology
  ci_orig <- sample_tb %>% filter(interp_age_uncert_pos != 'NULL') %>% distinct(entity_id)

  # entities only containing hiatuses
  only_h <-  dating %>% group_by(entity_id) %>% filter(all(date_type == 'Event; hiatus')) %>% distinct(entity_id)

  # entities not included in sample table
  no_sample_info <- entity %>% filter(!(entity_id %in% sample_tb$entity_id))

  # entities containing hiatuses
  hiatus_tb <- sample_tb %>% filter(hiatus == 'H') %>% distinct(entity_id)

  # entities containing reversals
  reversals <- dating_tb %>% filter(date_used == 'yes') %>%
    dplyr::select(entity_id, depth_dating, corr_age, corr_age_uncert_pos, corr_age_uncert_neg) %>%
    group_by(entity_id) %>%
    arrange(depth_dating, .by_group = T) %>%
    mutate(diff = lead(corr_age)-corr_age) %>%
    mutate(reversal = if_else(!is.na(diff) & diff < 0, TRUE, FALSE)) %>%
    group_by(entity_id) %>%
    arrange(depth_dating, .by_group = TRUE) %>%
    mutate(tractable = if_else(reversal & (abs(corr_age-lead(corr_age)) < (corr_age_uncert_pos + lead(corr_age_uncert_neg))), TRUE, FALSE))

  # entities with tractable/nontractable reversals
  reversal_tractable <- reversals %>% dplyr::select(entity_id, tractable) %>% group_by(entity_id) %>% filter(tractable) %>% distinct(entity_id)
  reversal_nontractable <- reversals %>% filter(reversal) %>%filter(!(entity_id %in% reversal_tractable$entity_id)) %>% distinct(entity_id)

  # entities with more than 3 dates
  nr_dates <- dating_tb %>% filter(date_used == 'yes' & date_used != 'Event; hiatus') %>% dplyr::select(entity_id, corr_age) %>% group_by(entity_id) %>%count() %>% filter(n>=3)

  # entities with missing sample depths
  no_depth_sample <- sample_tb %>% group_by(entity_id) %>% dplyr::summarise(depths = if_else(all(is.na(depth_sample)), FALSE, TRUE)) %>% filter(!depths)

  # entities that are not only U/Th dated
  not_UTh_dates <- dating_tb %>% filter(date_used == 'yes') %>% filter(date_type == 'Event; start of laminations' | date_type == 'Event; end of laminations' | date_type == 'C14' | date_type =='Multiple methods' | date_type =='other') %>%
    distinct(entity_id) 

  # entities that are only U/TH dated
  UTh_dates <- dating_tb %>% filter(date_used == 'yes') %>% filter(date_type == "MC-ICP-MS U/Th" | date_type == "TIMS" | date_type ==  "ICP-MS U/Th Other" | date_type == "Alpha U/Th" | date_type ==  "U/Th unspecified")  %>% 
    distinct(entity_id) %>% filter(!(entity_id %in% not_UTh_dates$entity_id))

  # entities containing only U/Th dates, enough dates, sample depths; 523
  run <- dating_tb_filtered %>% distinct(entity_id) %>%
    filter(entity_id %in% nr_dates$entity_id) %>%
    filter(!(entity_id %in% no_depth_sample$entity_id)) %>%
    filter(!(entity_id %in% not_UTh_dates$entity_id)) %>%
    filter(!(entity_id %in% no_sample_info$entity_id)) %>%
    filter(!(entity_id %in% only_h$entity_id)) %>%
    arrange(., entity_id)

  # entities that are not only U/Th dated, might have hiatuses, have enough dates 
  #run <- not_UTh_dates %>% filter(!(entity_id %in% no_depth_sample$entity_id)) %>% filter(entity_id %in% nr_dates$entity_id) %>% arrange(., entity_id)


  # entities where no AM can be computed (not enough dates)
  #run <- dating_tb %>% distinct(entity_id) %>%
  #  filter(!(entity_id %in% nr_dates$entity_id) | entity_id %in% no_depth_sample$entity_id | entity_id %in% no_sample_info$entity_id | entity_id %in% only_h$entity_id)
  
  ###############################################
  
  dating_tb_filtered_2 <- dating_tb_filtered %>% filter(entity_id %in% run$entity_id)
  sample_tb_new_filtered <- sample_tb_new %>% filter(entity_id %in% run$entity_id)
  site_tb_filtered <- site_tb %>% filter(entity_id %in% run$entity_id)


#return(list(site_tb, dating_tb, dating_tb_2, sample_tb))

return(list(site_tb_filtered, dating_tb_filtered_2, sample_tb_new_filtered, run))  
  
}











