load_sisal_data <- function(prefix, wd) {
  composite_link_entity <- read.csv(paste(wd, prefix, 'composite_link_entity.csv',sep = ''), header = T,stringsAsFactors = F)
  d13C <- read.csv(paste(wd, prefix, 'd13C.csv',sep='') ,header = T, stringsAsFactors = F)
  d13C <- rename(d13C, iso_std_d13C = iso_std )
  d18O <- read.csv(paste(wd, prefix, 'd18O.csv', sep =''),header = T, stringsAsFactors = F)
  d18O <- rename(d18O, iso_std_d18O = iso_std)
  dating_lamina <- read.csv(paste(wd, prefix, 'dating_lamina.csv', sep = ''), header = T, stringsAsFactors = F)
  dating <- read.csv(paste(wd, prefix, 'dating.csv',sep = ''), header = T, stringsAsFactors = F)
  entity_link_reference <- read.csv(paste(wd, prefix, 'entity_link_reference.csv', sep = ''), header =T, stringsAsFactors = F)
  entity <- read.csv(paste(wd, prefix, 'entity.csv', sep = ''), header = T, stringsAsFactors = F)
  gap <- read.csv(paste(wd, prefix, 'gap.csv', sep = ''), header = T, stringsAsFactors = F)
  hiatus <- read.csv(paste(wd, prefix, 'hiatus.csv', sep =''), header = T, stringsAsFactors = F)
  notes <- read.csv(paste(wd, prefix, 'notes.csv', sep = ''), header = T, stringsAsFactors = F)
  original_chronology <- read.csv(paste(wd, prefix, 'original_chronology.csv', sep = ''), header = T, stringsAsFactors = F)
  reference <- read.csv(paste(wd, prefix, 'reference.csv', sep = ''), header = T, stringsAsFactors = F)
  sample <- read.csv(paste(wd, prefix, 'sample.csv', sep = ''), header = T, stringsAsFactors = F)
  sisal_chronology <- read.csv(paste(wd, prefix, 'sisal_chronology.csv', sep = ''), header = T, stringsAsFactors = F)
  site <- read.csv(paste(wd, prefix, 'site.csv', sep = ''), header = T, stringsAsFactors = F)
  
  site_tb <- left_join(site, entity, by = 'site_id') %>% left_join(., entity_link_reference, by = 'entity_id') %>%
    left_join(., reference, by = 'ref_id') %>% left_join(., notes, by = 'site_id') %>% mutate_at(vars(site_id, entity_id), as.numeric)
  dating_tb <- dating %>% group_by(entity_id) %>%mutate(laminar_dated = if_else((entity_id %in% dating_lamina$entity_id), 'yes', 'no')) %>%
    mutate_at(vars(dating_id, depth_dating, dating_thickness, X14C_correction, corr_age, corr_age_uncert_pos, corr_age_uncert_neg), as.numeric) %>%ungroup()
  dating_tb_2 <- dating %>% left_join(.,entity, by = "entity_id") %>% filter(entity_status == "current") %>%
    mutate_at(vars(dating_id, depth_dating, dating_thickness, X14C_correction, corr_age, corr_age_uncert_pos, corr_age_uncert_neg), as.numeric)
  sample_tb <- join_all(list(sample,hiatus, gap, original_chronology, sisal_chronology, d13C, d18O), by = 'sample_id', type = 'left', match = 'all') %>%
    mutate_at(vars(entity_id, sample_id, sample_thickness, depth_sample, interp_age, interp_age_uncert_pos, interp_age_uncert_neg, COPRA_age,
                   COPRA_age_uncert_pos, COPRA_age_uncert_neg, linear_age, linear_age_uncert_pos, linear_age_uncert_neg, d13C_measurement,
                   d13C_precision, d18O_measurement, d18O_precision), as.numeric)
  return(list(site_tb, dating_tb, dating_tb_2, sample_tb))
}