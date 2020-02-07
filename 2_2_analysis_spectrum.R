#################################################
## Paper Figure SPECTRUM ########################
#################################################

library(plyr)
library(dplyr)
library(tidyverse)
library(zoo)
library(PaleoSpec)
library(nest)
library(latex2exp)

#################################################

## CALC

# We also use d18O measurement as it is not influenced by simulation

# We need 
# [ ] full yearly spectrum    TEMP, PREC, TEMP (prec weighted), ISOT, ITPC
# [ ] downsamples spectrum    TEMP, PREC, TEMP (prec weighted), ISOT, ITPC
# [ ] record spectrum
# [ ] gridbox weighing

ANALYSIS$SPECTRA <- list()
ANALYSIS$SPECTRA$RECORDS <- list()
ANALYSIS$SPECTRA$SIM_ds <- list("a" <- list(),"b" <- list(),"c" <- list())
ANALYSIS$SPECTRA$SIM_full <- list("a" <- list(),"b" <- list(),"c" <- list())
ANALYSIS$SPECTRA$MEAN_SPEC <- list()
ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH <- list()

## RECORDS:

length_cave <-length(DATA_past1000$CAVES$entity_info$entity_id)
entities_spec <- list()

# PaleoSpec::SpecMTM needs equally spaced data...

for(ii in 1:length_cave){
  print(ii)
  entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
  name = paste0("ENTITY", entity)
  if(length(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age) > 8 & ii != 95 & ii != 53 & ii != 109 & mask_spec[ii]){
    entities_spec = c(entities_spec, entity)
    start_ts = ceiling(head(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age, n = 1))
    length = length(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age)
    stop_ts = floor(tail(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age, n = 1))
    if(length > (stop_ts-start_ts)){length = (stop_ts-start_ts)}
    stop_ts = floor((stop_ts-start_ts)/length)*length+start_ts
    
    record <- PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age,DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement,
                                         time.target = seq(from = start_ts, to = stop_ts, by = floor((stop_ts-start_ts)/length)))
    record = na.omit(record)
    data <- ts(data = record, start = start_ts, end = stop_ts, deltat = floor((stop_ts-start_ts)/length))
    
    ANALYSIS$SPECTRA$RECORDS[[name]] = SpecMTM(data)
    
  }
  
}

ANALYSIS$SPECTRA$entities_spec_rec <- as.numeric(entities_spec)

ANALYSIS$SPECTRA$MEAN_SPEC[["Record"]] <- MeanSpectrum(ANALYSIS$SPECTRA$RECORDS)


## WEIGHING

entity_gridbox = array(dim = c(length(ANALYSIS$SPECTRA$entities_spec_rec), 4))
counter = 1

for(entity in ANALYSIS$SPECTRA$entities_spec_rec){
  entity_gridbox[counter, 1] = entity 
  entity_gridbox[counter, 2] = ANALYSIS$NETWORK$entity_meta$gridbox_id[ANALYSIS$NETWORK$entity_meta$entity_id == entity]
  entity_gridbox[counter, 3] = ANALYSIS$NETWORK$entity_meta$gridbox_lat[ANALYSIS$NETWORK$entity_meta$entity_id == entity]
  counter = counter+1
}

entity_gridbox <- as.data.frame(entity_gridbox)
colnames(entity_gridbox) <- c("entity_id", "gridbox_id", "gridbox_lat", "weighing")
total_gb <- entity_gridbox %>% select(gridbox_id, gridbox_lat) %>% group_by(gridbox_id) %>% count()
total_lat <- entity_gridbox %>% select(gridbox_lat) %>% group_by(gridbox_lat) %>% count()


for(entity in entity_gridbox$entity_id){
  gridbox = entity_gridbox$gridbox_id[entity_gridbox$entity_id == entity]
  lat = entity_gridbox$gridbox_lat[entity_gridbox$entity_id == entity]
  lat_weight = cos((90-lat*2.5)*pi/180)/sum(cos(seq(from = -90, to = 90, length.out = 73)*pi/180))/total_lat$n[total_lat$gridbox_lat == lat]
  entity_gridbox$weighing[entity_gridbox$entity_id == entity] = 1/dim(total_gb)[1]/total_gb$n[total_gb$gridbox_id == gridbox]*lat_weight
}

entity_gridbox$weighing <- entity_gridbox$weighing/sum(entity_gridbox$weighing)



ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[["Record"]] <- MeanSpectrum(ANALYSIS$SPECTRA$RECORDS, weights = entity_gridbox$weighing)


## SIMULATION DOWNSAMPLED
for(run in c("a","b","c")){
  ANALYSIS$SPECTRA$SIM_ds[[run]] <- list(
    TEMP = list(),
    PREC = list(),
    ISOT = list(), 
    ITPC = list()
  )
  # colnames(data) = c("interp_age", 
  #                    "TEMP_a", "PREC_a","ISOT_a", "ITPC_a", 
  #                    "TEMP_b", "PREC_b","ISOT_b", "ITPC_b", 
  #                    "TEMP_c", "PREC_c","ISOT_c", "ITPC_c")
  # 
  # DATA_past1000$CAVES$sim_data_downsampled[[nameE]] <- as.tibble(data)
  
  
  for(ii in 1:length(ANALYSIS$SPECTRA$entities_spec_rec)){
    print(ii)
    entity = ANALYSIS$SPECTRA$entities_spec_rec[ii]
    name = paste0("ENTITY", entity)
    
    start_ts = ceiling(head(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, n = 1))
    length = length(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age)
    stop_ts = floor(tail(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, n = 1))
    if(length > (stop_ts-start_ts)){length = (stop_ts-start_ts)}
    stop_ts = floor((stop_ts-start_ts)/length)*length+start_ts
    #ISOT
    record <- PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age,
                                         DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]][[paste0("ISOT", run)]],
                                         time.target = seq(from = start_ts, to = stop_ts, by = floor((stop_ts-start_ts)/length)))
    record = na.omit(record)
    data  = ts(data = record, start = start_ts, end   = stop_ts, deltat = floor((stop_ts-start_ts)/length))
    
    ANALYSIS$SPECTRA$SIM_ds[[run]]$ISOT[[name]] = SpecMTM(data)
    
    #ITPC
    record <- PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age,
                                         DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]][[paste0("ITPC", run)]],
                                         time.target = seq(from = start_ts, to = stop_ts, by = floor((stop_ts-start_ts)/length)))
    record = na.omit(record)
    data  = ts(data = record, start = start_ts, end   = stop_ts, deltat = floor((stop_ts-start_ts)/length))
    
    ANALYSIS$SPECTRA$SIM_ds[[run]]$ITPC[[name]] = SpecMTM(data)
    
    #TEMP
    record <- PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age,
                                         DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]][[paste0("TEMP", run)]],
                                         time.target = seq(from = start_ts, to = stop_ts, by = floor((stop_ts-start_ts)/length)))
    record = na.omit(record)
    data  = ts(data = record, start = start_ts, end   = stop_ts, deltat = floor((stop_ts-start_ts)/length))
    
    ANALYSIS$SPECTRA$SIM_ds[[run]]$TEMP[[name]] = SpecMTM(data)
    
    #PREC
    record <- PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age,
                                         DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]][[paste0("PREC", run)]],
                                         time.target = seq(from = start_ts, to = stop_ts, by = floor((stop_ts-start_ts)/length)))
    record = na.omit(record)
    data  = ts(data = record, start = start_ts, end   = stop_ts, deltat = floor((stop_ts-start_ts)/length))
    
    ANALYSIS$SPECTRA$SIM_ds[[run]]$PREC[[name]] = SpecMTM(data)
    
  }
  
  ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_ds_temp_", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_ds[[run]]$TEMP)
  ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_ds_prec_", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_ds[[run]]$PREC)
  ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_ds_isot_", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_ds[[run]]$ISOT)
  ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_ds_itpc_", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_ds[[run]]$ITPC)
  
  ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_ds_temp", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_ds[[run]]$TEMP, weights = entity_gridbox$weighing)
  ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_ds_prec", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_ds[[run]]$PREC, weights = entity_gridbox$weighing)
  ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_ds_isot", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_ds[[run]]$ISOT, weights = entity_gridbox$weighing)
  ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_ds_itpc", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_ds[[run]]$ITPC, weights = entity_gridbox$weighing)
  
  ## SIMULATION FULL 
  
  ANALYSIS$SPECTRA$SIM_full[[run]] <- list(
    TEMP = list(),
    PREC = list(),
    ISOT = list(), 
    ITPC = list()
  )
  
  for(ii in 1:length(ANALYSIS$SPECTRA$entities_spec_rec)){
    print(ii)
    entity = ANALYSIS$SPECTRA$entities_spec_rec[ii]
    site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
    name = paste0("ENTITY", entity)
    
    record = na.omit(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]][[paste0("ISOT", run)]])
    data  = ts(data = rev(record), 
               start = -49, 
               end   = 1100, 
               deltat = 1)
    
    ANALYSIS$SPECTRA$SIM_full[[run]]$ISOT[[name]] = SpecMTM(data)
    
    record = na.omit(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]][[paste0("ITPC", run)]])
    data  = ts(data = rev(record), 
               start = -49, 
               end   = 1100, 
               deltat = 1)
    
    ANALYSIS$SPECTRA$SIM_full[[run]]$ITPC[[name]] = SpecMTM(data)
    
    record = na.omit(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]][[paste0("TEMP", run)]])
    data  = ts(data = rev(record), 
               start = -49, 
               end   = 1100, 
               deltat = 1)
    
    ANALYSIS$SPECTRA$SIM_full[[run]]$TEMP[[name]] = SpecMTM(data)
    
    record = na.omit(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]][[paste0("PREC", run)]])
    data  = ts(data = rev(record), 
               start = -49, 
               end   = 1100, 
               deltat = 1)
    
    ANALYSIS$SPECTRA$SIM_full[[run]]$PREC[[name]] = SpecMTM(data)
    
  }
  
  ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_temp", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_full[[run]]$TEMP)
  ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_prec", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_full[[run]]$PREC)
  ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_isot", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_full[[run]]$ISOT)
  ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_itpc", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_full[[run]]$ITPC)
  
  ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_full_temp", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_full[[run]]$TEMP, weights = entity_gridbox$weighing)
  ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_full_prec", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_full[[run]]$PREC, weights = entity_gridbox$weighing)
  ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_full_isot", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_full[[run]]$ISOT, weights = entity_gridbox$weighing)
  ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_full_itpc", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_full[[run]]$ITPC, weights = entity_gridbox$weighing)
}


#################################################
## Sim full FILTER Spectra ######################
# filtere full Sim signal auf down sim signal
# filtere full sim signal auf rec signal
# filtere down sim signal auf rec signal

source('Functions/Filter/EASY_Sensor_WM4.R')
source('Functions/Filter/filter_function3.R')


ANALYSIS$SPECTRA$SIM_filter_full_down <- list(TEMP = list(), PREC = list(), ISOT = list(), ITPC = list())
ANALYSIS$SPECTRA$SIM_filter_full_rec <- list(TEMP = list(), PREC = list(), ISOT = list(), ITPC = list())
ANALYSIS$SPECTRA$SIM_filter_down_rec <- list(TEMP = list(), PREC = list(), ISOT = list(), ITPC = list())

filter1 = 2.0
filter2 = 20.0
filter3 = 10.0

for(ii in 1:length(ANALYSIS$SPECTRA$entities_spec_rec)){
  print(ii)
  entity = ANALYSIS$SPECTRA$entities_spec_rec[ii]
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  name = paste0("ENTITY", entity)
  
  #TEMP
  Results <- easy_sensor_wm4(1.0, na.omit(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$TEMP), filter1)
  data  = ts(data = rev(Results), start = (-49-10.0), end   = 1100, deltat = 1)
  
  ANALYSIS$SPECTRA$SIM_filter_full_down$TEMP[[name]] = SpecMTM(data)
  
  #PREC
  Results <- easy_sensor_wm4(1.0, na.omit(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$PREC), filter1)
  data  = ts(data = rev(Results), start = (-49-10.0), end   = 1100, deltat = 1)
  
  ANALYSIS$SPECTRA$SIM_filter_full_down$PREC[[name]] = SpecMTM(data)
  
  #ISOT
  Results <- easy_sensor_wm4(1.0, na.omit(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ISOT),filter1)
  data  = ts(data = rev(Results), start = (-49-10.0), end   = 1100, deltat = 1)
  
  ANALYSIS$SPECTRA$SIM_filter_full_down$ISOT[[name]] = SpecMTM(data)
  
  #ITPC
  Results <- easy_sensor_wm4(1.0, na.omit(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ITPC),filter1)
  data  = ts(data = rev(Results), start = (-49-10.0), end   = 1100, deltat = 1)
  
  ANALYSIS$SPECTRA$SIM_filter_full_down$ITPC[[name]] = SpecMTM(data)
    
  
  #full-< rec
  #TEMP
  Results <- easy_sensor_wm4(1.0, na.omit(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$TEMP), filter2)
  data  = ts(data = rev(Results), start = (-49-10.0), end   = 1100, deltat = 1)
  
  ANALYSIS$SPECTRA$SIM_filter_full_rec$TEMP[[name]] = SpecMTM(data)
  
  #PREC
  Results <- easy_sensor_wm4(1.0, na.omit(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$PREC), filter2)
  data  = ts(data = rev(Results), start = (-49-10.0), end   = 1100, deltat = 1)
  
  ANALYSIS$SPECTRA$SIM_filter_full_rec$PREC[[name]] = SpecMTM(data)
  
  #ISOT
  Results <- easy_sensor_wm4(1.0, na.omit(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ISOT), filter2)
  data  = ts(data = rev(Results), start = (-49-10.0), end   = 1100, deltat = 1)
  
  ANALYSIS$SPECTRA$SIM_filter_full_rec$ISOT[[name]] = SpecMTM(data)
  
  #ITPC
  Results <- easy_sensor_wm4(1.0, na.omit(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ITPC), filter2)
  data  = ts(data = rev(Results), start = (-49-10.0), end   = 1100, deltat = 1)
  
  ANALYSIS$SPECTRA$SIM_filter_full_rec$ITPC[[name]] = SpecMTM(data)
  
  #down-<rec
  start_ts = ceiling(head(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, n = 1))
  length = length(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age)
  stop_ts = floor(tail(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, n = 1))
  diff = floor((stop_ts-start_ts)/length)
  if(diff<1){diff = 1}
  
  #TEMP
  record <- PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age,DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$TEMP,
                                       time.target = seq(from = start_ts, to = stop_ts, by = diff))
  Results <- easy_sensor_wm4(diff, na.omit(rev(record)), filter3)
  data = ts(data = Results, start = start_ts-10, deltat = diff)
  ANALYSIS$SPECTRA$SIM_filter_down_rec$TEMP[[name]] = SpecMTM(data)
  
  #PREC
  record <- PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age,DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$PREC,
                                       time.target = seq(from = start_ts, to = stop_ts, by = diff))
  Results <- easy_sensor_wm4(diff, na.omit(rev(record)), filter3)
  data = ts(data = Results, start = start_ts-10, deltat = diff)
  ANALYSIS$SPECTRA$SIM_filter_down_rec$PREC[[name]] = SpecMTM(data)
  
  #ISOT
  record <- PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age,DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ISOT,
                                       time.target = seq(from = start_ts, to = stop_ts, by = diff))
  Results <- easy_sensor_wm4(diff, na.omit(rev(record)), filter3)
  data = ts(data = Results, start = start_ts-10, deltat = diff)
  ANALYSIS$SPECTRA$SIM_filter_down_rec$ISOT[[name]] = SpecMTM(data)
  
  #ITPC
  record <- PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age,DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ITPC,
                                       time.target = seq(from = start_ts, to = stop_ts, by = diff))
  Results <- easy_sensor_wm4(diff, na.omit(rev(record)), filter3)
  data = ts(data = Results, start = start_ts-10, deltat = diff)
  ANALYSIS$SPECTRA$SIM_filter_down_rec$ITPC[[name]] = SpecMTM(data)
  
}

ANALYSIS$SPECTRA$MEAN_SPEC[["SIM_full_down_temp"]]<- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_full_down$TEMP)
ANALYSIS$SPECTRA$MEAN_SPEC[["SIM_full_down_prec"]]<- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_full_down$PREC)
ANALYSIS$SPECTRA$MEAN_SPEC[["SIM_full_down_isot"]]<- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_full_down$ISOT)
ANALYSIS$SPECTRA$MEAN_SPEC[["SIM_full_down_itpc"]]<- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_full_down$ITPC)

ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[["SIM_full_down_temp"]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_full_down$TEMP, weights = entity_gridbox$weighing)
ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[["SIM_full_down_prec"]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_full_down$PREC, weights = entity_gridbox$weighing)
ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[["SIM_full_down_isot"]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_full_down$ISOT, weights = entity_gridbox$weighing)
ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[["SIM_full_down_itpc"]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_full_down$ITPC, weights = entity_gridbox$weighing)


ANALYSIS$SPECTRA$MEAN_SPEC[["SIM_full_rec_temp"]]<- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_full_rec$TEMP)
ANALYSIS$SPECTRA$MEAN_SPEC[["SIM_full_rec_prec"]]<- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_full_rec$PREC)
ANALYSIS$SPECTRA$MEAN_SPEC[["SIM_full_rec_isot"]]<- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_full_rec$ISOT)
ANALYSIS$SPECTRA$MEAN_SPEC[["SIM_full_rec_itpc"]]<- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_full_rec$ITPC)

ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[["SIM_full_rec_temp"]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_full_rec$TEMP, weights = entity_gridbox$weighing)
ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[["SIM_full_rec_prec"]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_full_rec$PREC, weights = entity_gridbox$weighing)
ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[["SIM_full_rec_isot"]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_full_rec$ISOT, weights = entity_gridbox$weighing)
ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[["SIM_full_rec_itpc"]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_full_rec$ITPC, weights = entity_gridbox$weighing)


ANALYSIS$SPECTRA$MEAN_SPEC[["SIM_down_rec_temp"]]<- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_down_rec$TEMP)
ANALYSIS$SPECTRA$MEAN_SPEC[["SIM_down_rec_prec"]]<- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_down_rec$PREC)
ANALYSIS$SPECTRA$MEAN_SPEC[["SIM_down_rec_isot"]]<- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_down_rec$ISOT)
ANALYSIS$SPECTRA$MEAN_SPEC[["SIM_down_rec_itpc"]]<- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_down_rec$ITPC)

ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[["SIM_down_rec_temp"]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_down_rec$TEMP, weights = entity_gridbox$weighing)
ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[["SIM_down_rec_prec"]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_down_rec$PREC, weights = entity_gridbox$weighing)
ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[["SIM_down_rec_isot"]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_down_rec$ISOT, weights = entity_gridbox$weighing)
ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[["SIM_down_rec_itpc"]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_down_rec$ITPC, weights = entity_gridbox$weighing)

