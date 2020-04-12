#################################################
## DATA EXPORT ##################################
#################################################

DATA_EXPORT_SITE <- list()
DATA_EXPORT_ENTITY <- list()

#################################################

# Tabelle mit Site - year - temp - prec - isot
# Tabelle mit Entity - year - temp - prec - isot

## SITE 

for(site in DATA_past1000$CAVES$site_info$site_id){
  data_new = array(dim = c(diff(DATA_past1000$time),6))
  colnames(data_new) = c("site_id", "year(BP)", "mean_temp(°C)", "mean_prec(kg m^-2 s^-1)", "mean_d18O (permil)", "d18Opw (permil)")
  data_new[,1] = site
  data_new[,2] = seq(from = 1950 - DATA_past1000$time[1],to =  1950-DATA_past1000$time[2]+1, by = -1)
  data_new[,3] = DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$TEMP_a
  data_new[,4] = DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$PREC_a
  data_new[,5] = DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ISOT_a
  data_new[,6] = DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ITPC_a
  
  if(site == 1){ data = data_new}
  else{
    data = rbind(data, data_new)
  }
}

DATA_EXPORT_SITE <- data

##ENTITY

for(entity in DATA_past1000$CAVES$entity_info$entity_id){
  data_new = array(dim = c(length(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$TEMP_a),6))
  colnames(data_new) = c("site_id", "year(BP)", "mean_temp(°C)", "mean_prec(kg m^-2 s^-1)", "mean_d18O (permil)", "d18Opw (permil)")
  data_new[,1] = entity
  data_new[,2] = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age
  data_new[,3] = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$TEMP_a
  data_new[,4] = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$PREC_a
  data_new[,5] = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ISOT_a
  data_new[,6] = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ITPC_a
  
  if(entity == 9){ data = data_new}
  else{
    data = rbind(data, data_new)
  }
}

DATA_EXPORT_ENTITY <- data

write.csv(DATA_EXPORT_SITE, file = "SISAL_HadCM3xnapa_PMIL_yearly.csv", row.names = F)
write.csv(DATA_EXPORT_ENTITY, file = "SISAL_HadCM3xnapa_PMIL_downsampled.csv", row.names = F)
