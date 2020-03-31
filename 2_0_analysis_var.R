#################################################
## ANALYSIS VAR #################################
#################################################


# neccessary analysis that is too big for the analysis in the plotting skripts

#ANALYSIS <- list()

#################################################

## MEAN #########################################

ANALYSIS$MEAN <- list()

mean_bias_full <- list()
mean_bias_ds <- list()
for(entity in DATA_past1000$CAVES$entity_info$entity_id[mask_mean]){
  site = DATA_past1000$CAVES$entity_info$site_id[match(entity, DATA_past1000$CAVES$entity_info$entity_id)]
  mean_bias_full <- c(mean_bias_full, mean(c(mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ITPC_a, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_dw_eq_a, na.rm = T),
                                             mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ITPC_b, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_dw_eq_b, na.rm = T),
                                             mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ITPC_c, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_dw_eq_c, na.rm = T)), na.rm = T))
  mean_bias_ds <- c(mean_bias_ds, mean(c(mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ITPC_a, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_dw_eq_a, na.rm = T),
                                         mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ITPC_b, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_dw_eq_b, na.rm = T),
                                         mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ITPC_c, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_dw_eq_c, na.rm = T)), na.rm = T))
}

ANALYSIS$MEAN$global_mean_ds <- mean(as.numeric(mean_bias_ds), na.rm = T)
ANALYSIS$MEAN$global_mean_full <- mean(as.numeric(mean_bias_full), na.rm = T)

# CLUSTER

cluster_mean <- list()

for(cluster in 1:9){
  print(cluster)
  entity_list <- ANALYSIS$NETWORK$entity_meta %>% filter(cluster_id == cluster)
  entity_list <- entity_list$entity_id
  
  mean_bias_full <- list()
  mean_bias_ds <- list()
  for(entity in entity_list){
    site = DATA_past1000$CAVES$entity_info$site_id[match(entity, DATA_past1000$CAVES$entity_info$entity_id)]
    mean_bias_full <- c(mean_bias_full, mean(c(mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ITPC_a, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_dw_eq_a, na.rm = T),
                                               mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ITPC_b, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_dw_eq_b, na.rm = T),
                                               mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ITPC_c, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_dw_eq_c, na.rm = T)), na.rm = T))
    mean_bias_ds <- c(mean_bias_ds, mean(c(mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ITPC_a, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_dw_eq_a, na.rm = T),
                                           mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ITPC_b, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_dw_eq_b, na.rm = T),
                                           mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ITPC_c, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_dw_eq_c, na.rm = T)), na.rm = T))
  }
  print(quantile(as.numeric(mean_bias_ds), prob = seq(0,1,0.05)))
  cluster_mean[paste0("CLUSTER", cluster, "_ds")] <- mean(as.numeric(mean_bias_ds), na.rm = T)
  cluster_mean[paste0("CLUSTER", cluster, "_full")] <- mean(as.numeric(mean_bias_full), na.rm = T)
}

ANALYSIS$MEAN$CLUSTER <- cluster_mean

#################################################

source("Functions/aw_mean.R")

## VARIANCE #####################################

ANALYSIS$VARIANCE <- list()

for(run in c("a", "b", "c")){
  ANALYSIS$VARIANCE$FIELDS[[run]]$ISOTlyr <- array(dim = c(96,73))
  ANALYSIS$VARIANCE$FIELDS[[run]]$ITPClyr <- array(dim = c(96,73))
  ANALYSIS$VARIANCE$FIELDS[[run]]$TEMPlyr <- array(dim = c(96,73))
  ANALYSIS$VARIANCE$FIELDS[[run]]$PREClyr <- array(dim = c(96,73))
}

for(run in c("a", "b", "c")){
  for(lon in 1:96){
    for(lat in 1:73){
      ANALYSIS$VARIANCE$FIELDS[[run]]$ISOTlyr[lon,lat] = var(DATA_past1000[[paste0("SIM_yearly_", run)]]$ISOT[lon,lat,])
      ANALYSIS$VARIANCE$FIELDS[[run]]$ITPClyr[lon,lat] = var(DATA_past1000[[paste0("SIM_yearly_", run)]]$ITPC[lon,lat,])
      ANALYSIS$VARIANCE$FIELDS[[run]]$TEMPlyr[lon,lat] = var(DATA_past1000[[paste0("SIM_yearly_", run)]]$TEMP[lon,lat,])
      ANALYSIS$VARIANCE$FIELDS[[run]]$PREClyr[lon,lat] = var(DATA_past1000[[paste0("SIM_yearly_", run)]]$PREC[lon,lat,])
    }
  }
  
}


remove(lat,lon, run)


for(run in c("a", "b", "c")){
  ANALYSIS$VARIANCE$FIELDS[[run]]$ISOTlyr <- rbind(ANALYSIS$VARIANCE$FIELDS[[run]]$ISOT[49:96,1:73], ANALYSIS$VARIANCE$FIELDS[[run]]$ISOT[1:48, 1:73])
  ANALYSIS$VARIANCE$FIELDS[[run]]$ITPClyr <- rbind(ANALYSIS$VARIANCE$FIELDS[[run]]$ITPC[49:96,1:73], ANALYSIS$VARIANCE$FIELDS[[run]]$ITPC[1:48, 1:73])
  ANALYSIS$VARIANCE$FIELDS[[run]]$TEMPlyr <- rbind(ANALYSIS$VARIANCE$FIELDS[[run]]$TEMP[49:96,1:73], ANALYSIS$VARIANCE$FIELDS[[run]]$TEMP[1:48, 1:73])
  ANALYSIS$VARIANCE$FIELDS[[run]]$PREClyr <- rbind(ANALYSIS$VARIANCE$FIELDS[[run]]$PREC[49:96,1:73], ANALYSIS$VARIANCE$FIELDS[[run]]$PREC[1:48, 1:73])
  
}

ANALYSIS$VARIANCE$POINTS$CAVElyr <- list(
  lon = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  lat = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  value_record = numeric(length(DATA_past1000$CAVES$entity_info$entity_id))
)

for(run in c("a", "b", "c")){
  ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_temp = numeric(length(DATA_past1000$CAVES$entity_info$entity_id))
  ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_temp_ds = numeric(length(DATA_past1000$CAVES$entity_info$entity_id))
  ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_prec = numeric(length(DATA_past1000$CAVES$entity_info$entity_id))
  ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_prec_ds = numeric(length(DATA_past1000$CAVES$entity_info$entity_id))
  ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_isot = numeric(length(DATA_past1000$CAVES$entity_info$entity_id))
  ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_isot_ds = numeric(length(DATA_past1000$CAVES$entity_info$entity_id))
  ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_itpc = numeric(length(DATA_past1000$CAVES$entity_info$entity_id))
  ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_itpc_ds = numeric(length(DATA_past1000$CAVES$entity_info$entity_id))
}

for(run in c("a", "b", "c")){
  Var_temp_norm <- simpleawmean(ANALYSIS$VARIANCE$FIELDS[[run]]$TEMPlyr, seq(from = -90, to = 90, length.out = 73))
  Var_prec_norm <- simpleawmean(ANALYSIS$VARIANCE$FIELDS[[run]]$PREClyr, seq(from = -90, to = 90, length.out = 73))
  Var_isot_norm <- simpleawmean(ANALYSIS$VARIANCE$FIELDS[[run]]$ISOTlyr, seq(from = -90, to = 90, length.out = 73))
  Var_itpc_norm <- simpleawmean(ANALYSIS$VARIANCE$FIELDS[[run]]$ITPClyr, seq(from = -90, to = 90, length.out = 73))
  
  for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
    site = DATA_past1000$CAVES$entity_info$site_id[ii]
    entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
    ANALYSIS$VARIANCE$POINTS$CAVElyr$lon[ii]   = DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site]
    ANALYSIS$VARIANCE$POINTS$CAVElyr$lat[ii]   = DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site]
    ANALYSIS$VARIANCE$POINTS$CAVElyr$value_record[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]][[paste0("d18O_dw_eq_",run)]], na.rm = T)
    #Temp
    # Take normalized Variances
    ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_temp[ii]    = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]][[paste0("d18O_dw_eq_",run)]], na.rm = T)/
      var(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]][[paste0("TEMP_", run)]], na.rm = T)*Var_temp_norm/Var_isot_norm
    ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_temp_ds[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]][[paste0("d18O_dw_eq_",run)]], na.rm = T)/
      var(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]][[paste0("TEMP_", run)]], na.rm = T) * Var_temp_norm/Var_isot_norm
    #Prec
    # Take normalized Variances
    ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_prec[ii]    = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]][[paste0("d18O_dw_eq_",run)]], na.rm = T)/
      var(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]][[paste0("PREC_", run)]], na.rm = T)*Var_prec_norm/Var_isot_norm
    ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_prec_ds[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]][[paste0("d18O_dw_eq_",run)]], na.rm = T)/
      var(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]][[paste0("PREC_", run)]], na.rm = T) * Var_prec_norm/Var_isot_norm
    # ISOT
    ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_isot[ii]    = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]][[paste0("d18O_dw_eq_",run)]], na.rm = T)/
      var(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]][[paste0("ISOT_", run)]], na.rm = T)
    ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_isot_ds[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]][[paste0("d18O_dw_eq_",run)]], na.rm = T)/
      var(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]][[paste0("ISOT_", run)]], na.rm = T)
    # ITPC
    ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_itpc[ii]    = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]][[paste0("d18O_dw_eq_",run)]], na.rm = T)/
      var(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]][[paste0("ITPC_", run)]], na.rm = T)
    ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_itpc_ds[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]][[paste0("d18O_dw_eq_",run)]], na.rm = T)/
      var(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]][[paste0("ITPC_", run)]], na.rm = T)
  }  
}

vr_tracker <- data.frame(
  entity_id = DATA_past1000$CAVES$entity_info$entity_id[mask_spec],
  var_itpc_ds = ANALYSIS$VARIANCE$POINTS$CAVElyr$a$value_VR_isot_ds[mask_spec],
  cluster_id = ANALYSIS$NETWORK$entity_meta$cluster_id
)

cluster_mean <- numeric(9)
for(cluster in 1:9){
  cluster_var = vr_tracker %>% filter(cluster_id == cluster)
  cluster_mean[cluster] = mean(cluster_var$var_itpc_ds, na.rm = T) 
}

ANALYSIS$VARIANCE$CLUSTER_MEAN <- cluster_mean

remove(entity, ii, site, Var_isot_norm, Var_itpc_norm, Var_prec_norm, Var_temp_norm, run)
