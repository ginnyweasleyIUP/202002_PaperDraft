#################################################
## ANALYSIS VAR #################################
#################################################


# neccessary analysis that is too big for the analysis in the plotting skripts

ANALYSIS <- list()

#################################################

source("Functions/aw_mean.R")

## VARIANCE #####################################

ANALYSIS$VARIANCE <- list()

ANALYSIS$VARIANCE$FIELDS$ISOTlyr <- array(dim = c(96,73))
ANALYSIS$VARIANCE$FIELDS$ITPClyr <- array(dim = c(96,73))
ANALYSIS$VARIANCE$FIELDS$TEMPlyr <- array(dim = c(96,73))
ANALYSIS$VARIANCE$FIELDS$PREClyr <- array(dim = c(96,73))

for(lon in 1:96){
  for(lat in 1:73){
    ANALYSIS$VARIANCE$FIELDS$ISOTlyr[lon,lat] = var(DATA_past1000$SIM_yearly$ISOT[lon,lat,])
    ANALYSIS$VARIANCE$FIELDS$ITPClyr[lon,lat] = var(DATA_past1000$SIM_yearly$ITPC[lon,lat,])
    ANALYSIS$VARIANCE$FIELDS$TEMPlyr[lon,lat] = var(DATA_past1000$SIM_yearly$TEMP[lon,lat,])
    ANALYSIS$VARIANCE$FIELDS$PREClyr[lon,lat] = var(DATA_past1000$SIM_yearly$PREC[lon,lat,])
  }
}

remove(lat,lon)



ANALYSIS$VARIANCE$FIELDS$ISOTlyr <- rbind(ANALYSIS$VARIANCE$FIELDS$ISOT[49:96,1:73],
                                          ANALYSIS$VARIANCE$FIELDS$ISOT[1:48, 1:73])

ANALYSIS$VARIANCE$FIELDS$ITPClyr <- rbind(ANALYSIS$VARIANCE$FIELDS$ITPC[49:96,1:73],
                                          ANALYSIS$VARIANCE$FIELDS$ITPC[1:48, 1:73])

ANALYSIS$VARIANCE$FIELDS$TEMPlyr <- rbind(ANALYSIS$VARIANCE$FIELDS$TEMP[49:96,1:73],
                                          ANALYSIS$VARIANCE$FIELDS$TEMP[1:48, 1:73])

ANALYSIS$VARIANCE$FIELDS$PREClyr <- rbind(ANALYSIS$VARIANCE$FIELDS$PREC[49:96,1:73],
                                          ANALYSIS$VARIANCE$FIELDS$PREC[1:48, 1:73])

ANALYSIS$VARIANCE$POINTS$CAVElyr <- data.frame(
  lon = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  lat = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)), 
  value_record = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  value_VR_temp = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  value_VR_temp_ds = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  value_VR_prec = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  value_VR_prec_ds = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  value_VR_isot = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  value_VR_isot_ds = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  value_VR_itpc = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  value_VR_itpc_ds = numeric(length(DATA_past1000$CAVES$entity_info$entity_id))
)

Var_temp_norm <- simpleawmean(ANALYSIS$VARIANCE$FIELDS$TEMPlyr, seq(from = -90, to = 90, length.out = 73))
Var_prec_norm <- simpleawmean(ANALYSIS$VARIANCE$FIELDS$PREClyr, seq(from = -90, to = 90, length.out = 73))
Var_isot_norm <- simpleawmean(ANALYSIS$VARIANCE$FIELDS$ISOTlyr, seq(from = -90, to = 90, length.out = 73))
Var_itpc_norm <- simpleawmean(ANALYSIS$VARIANCE$FIELDS$ITPClyr, seq(from = -90, to = 90, length.out = 73))

for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
  site = DATA_past1000$CAVES$entity_info$site_id[ii]
  entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
  ANALYSIS$VARIANCE$POINTS$CAVElyr$lon[ii]   = DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site]
  ANALYSIS$VARIANCE$POINTS$CAVElyr$lat[ii]   = DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site]
  ANALYSIS$VARIANCE$POINTS$CAVElyr$value[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement, na.rm = T)
  #Temp
  # Take normalized Variances
  ANALYSIS$VARIANCE$POINTS$CAVElyr$value_VR_temp[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement, na.rm = T)/
    var(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$TEMP, na.rm = T)*Var_temp_norm/Var_isot_norm
  ANALYSIS$VARIANCE$POINTS$CAVElyr$value_VR_temp_ds[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement, na.rm = T)/
    var(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$TEMP, na.rm = T) * Var_temp_norm/Var_isot_norm
  #Prec
  #Take normalized variances
  ANALYSIS$VARIANCE$POINTS$CAVElyr$value_VR_prec[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement, na.rm = T)/
    var(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$PREC, na.rm = T) * Var_prec_norm/Var_isot_norm
  ANALYSIS$VARIANCE$POINTS$CAVElyr$value_VR_prec_ds[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement, na.rm = T)/
    var(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$PREC, na.rm = T) * Var_prec_norm/Var_isot_norm
  #Isot "d18O"
  ANALYSIS$VARIANCE$POINTS$CAVElyr$value_VR_isot[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement, na.rm = T)/
    var(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ISOT, na.rm = T)
  ANALYSIS$VARIANCE$POINTS$CAVElyr$value_VR_isot_ds[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement, na.rm = T)/
    var(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ISOT, na.rm = T)
  #ITPC  "d18O_pw"
  ANALYSIS$VARIANCE$POINTS$CAVElyr$value_VR_itpc[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement, na.rm = T)/
    var(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ITPC, na.rm = T)
  ANALYSIS$VARIANCE$POINTS$CAVElyr$value_VR_itpc_ds[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement, na.rm = T)/
    var(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ITPC, na.rm = T)
}

remove(entity, ii, site, Var_isot_norm, Var_itpc_norm, Var_prec_norm, Var_temp_norm)
