#################################################
## SPECTRAL ANALYSIS ############################
#################################################

## ToDo:
# [ ] Mean Spectrum of Sim at Cave Sites vs Rec at Cave Site
# [ ] Mean Spectrum of down sampled Sim at Cave Site vs. Rec at Cave Site
# [ ]
# [ ]
# [ ]
# [ ]
# [ ]

SPEC_ANALYSIS <- list()

#################################################
## Record Spectra ###############################

library(zoo)
library(PaleoSpec)
library(nest)
library(latex2exp)

SPEC_ANALYSIS$RECORDS <- list()

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
    
    SPEC_ANALYSIS$RECORDS[[name]] = SpecMTM(data)
    
  }
    
}

SPEC_ANALYSIS$GENERAL$entities_spec_rec <- entities_spec

SPEC_ANALYSIS$MEAN_SPEC[["Record"]] <- MeanSpectrum(SPEC_ANALYSIS$RECORDS)

#################################################
## Sim DS Spectra ###############################

SPEC_ANALYSIS$SIM_ds <- list()

length_cave <-length(DATA_past1000$CAVES$entity_info$entity_id)
entities_spec <- list()

# PaleoSpec::SpecMTM needs equally spaced data...

for(ii in 1:length_cave){
  entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
  name = paste0("ENTITY", entity)
  if(length(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age) > 8 & ii != 95& ii != 53 & ii != 109 & mask_spec[ii]){
    entities_spec = c(entities_spec, entity)
    start_ts = ceiling(head(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, n = 1))
    length = length(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age)
    stop_ts = floor(tail(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, n = 1))
    if(length > (stop_ts-start_ts)){length = (stop_ts-start_ts)}
    stop_ts = floor((stop_ts-start_ts)/length)*length+start_ts
    record <- PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age,DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$d18O,
                                         time.target = seq(from = start_ts, to = stop_ts, by = floor((stop_ts-start_ts)/length)))
    record = na.omit(record)
    data  = ts(data = record, start = start_ts, end   = stop_ts, deltat = floor((stop_ts-start_ts)/length))
    
    SPEC_ANALYSIS$SIM_ds[[name]] = SpecMTM(data)
    
  }
  
}

SPEC_ANALYSIS$GENERAL$entities_spec_sim_ds <- entities_spec

SPEC_ANALYSIS$MEAN_SPEC[["SIM_ds"]] <- MeanSpectrum(SPEC_ANALYSIS$SIM_ds)


#################################################
## Sim full Spectra #############################

SPEC_ANALYSIS$SIM_full <- list()

length_cave <-length(DATA_past1000$CAVES$entity_info$entity_id)
entities_spec <- list()

for(ii in 1:length_cave){
  entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
  site = DATA_past1000$CAVES$entity_info$site_id[ii]
  name = paste0("ENTITY", entity)
  if(length(DATA_past1000$CAVES$sim_data_raw[[paste0("CAVE", site)]]$ISOT) > 8 & mask_spec[ii]){
    entities_spec = c(entities_spec, entity)
    record = na.omit(DATA_past1000$CAVES$sim_data_raw[[paste0("CAVE", site)]]$ISOT)
    data  = ts(data = rev(record), 
               start = -49, 
               end   = 1100, 
               deltat = 1)
    
    SPEC_ANALYSIS$SIM_full[[name]] = SpecMTM(data)
    
  }
  
}

SPEC_ANALYSIS$GENERAL$entities_spec_sim_full <- entities_spec

SPEC_ANALYSIS$MEAN_SPEC[["SIM_full"]] <- MeanSpectrum(SPEC_ANALYSIS$SIM_full)

#################################################
## Sim full FILTER Spectra ######################
# filtere full Sim signal auf down sim signal
# filtere full sim signal auf rec signal
# filtere down sim signal auf rec signal

source('Functions/Filter/bwf_filter.R')
source('Functions/Filter/EASY_Sensor_WM4.R')
source('Functions/Filter/filter_function3.R')


SPEC_ANALYSIS$SIM_filter_full_down <- list()
SPEC_ANALYSIS$SIM_filter_full_rec <- list()
SPEC_ANALYSIS$SIM_filter_down_rec <- list()

length_cave <-length(DATA_past1000$CAVES$entity_info$entity_id)
entities_spec <- list()

for(ii in 1:length_cave){
  entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
  site = DATA_past1000$CAVES$entity_info$site_id[ii]
  name = paste0("ENTITY", entity)
  if(length(DATA_past1000$CAVES$sim_data_raw[[paste0("CAVE", site)]]$ISOT) > 8 & mask_spec[ii]){
    entities_spec = c(entities_spec, entity)
    Results <- easy_sensor_wm4(1.0, 
                               na.omit(DATA_past1000$CAVES$sim_data_raw[[paste0("CAVE", site)]]$ISOT),
                               5.0, na.omit(DATA_past1000$CAVES$sim_data_raw[[paste0("CAVE", site)]]$TEMP))
    data  = ts(data = rev(Results), 
               start = (-49-10.0), 
               end   = 1100, 
               deltat = 1)
    
    SPEC_ANALYSIS$SIM_filter_full_down[[name]] = SpecMTM(data)
    
    Results <- easy_sensor_wm4(1.0, 
                               na.omit(DATA_past1000$CAVES$sim_data_raw[[paste0("CAVE", site)]]$ISOT), 
                               50.0, 
                               na.omit(DATA_past1000$CAVES$sim_data_raw[[paste0("CAVE", site)]]$TEMP))
    data  = ts(data = rev(Results), 
               start = (-49-10.0), 
               end   = 1100, 
               deltat = 1)
    
    SPEC_ANALYSIS$SIM_filter_full_rec[[name]] = SpecMTM(data)
    
    
    start_ts = ceiling(head(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, n = 1))
    length = length(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age)
    stop_ts = floor(tail(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, n = 1))
    diff = floor((stop_ts-start_ts)/length)
    if(diff<1){diff = 1}
    record <- PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age,DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$d18O,
                                         time.target = seq(from = start_ts, to = stop_ts, by = diff))
    
    Results <- easy_sensor_wm4(diff, 
                               na.omit(rev(record)), 
                               10.0, 
                               na.omit(record))
    
    data = ts(data = Results, start = start_ts-10, deltat = diff)
    
    SPEC_ANALYSIS$SIM_filter_down_rec[[name]] = SpecMTM(data)
    
  }
  
}

SPEC_ANALYSIS$GENERAL$entities_spec_sim_full_down <- entities_spec

SPEC_ANALYSIS$MEAN_SPEC[["SIM_full_down"]] <- MeanSpectrum(SPEC_ANALYSIS$SIM_filter_full_down)
SPEC_ANALYSIS$MEAN_SPEC[["SIM_full_rec"]] <- MeanSpectrum(SPEC_ANALYSIS$SIM_filter_full_rec)
SPEC_ANALYSIS$MEAN_SPEC[["SIM_down_rec"]] <- MeanSpectrum(SPEC_ANALYSIS$SIM_filter_down_rec)


#################################################
## AREA WEIGHING ################################
#################################################

# Hier müssen alle Speleotheme innerhalb einer Grid box zu einer Box geaveraged werden, bevor es dann global geaveraged wird. 

# 1) Finde heraus, welche Speleotheme innerhalb einer Grid Box sind
# 2) Average dort

list_spectra_gridboxes <- list()

for(lon in 1:96){
  for(lat in 1:72){
    mask_entities <- logical(length = DATA_past1000$CAVES$entity_info$entity_id[mask_spec])
    long_start = 360/96*lon - 
    long_stop =
    lat_tart = 
    lat_stop =
    for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])){
      entity = DATA_past1000$CAVES$entity_info$entity_id[mask_spec][ii]
      
    }
  }
}


#################################################
## Plotting #####################################
#################################################

pdf(file = "Plots/Spectra/Global_Spectra_compare.pdf", width = 8, height = 5)
LPlot(SPEC_ANALYSIS$MEAN_SPEC$SIM_full$spec, col = "#074893", 
      ylim = c(0.00001,500), xlim = c(1/500, 0.5),
      xaxt = 'n',
      xlab = "",
      main = TeX("Mean Spectra from cave locations (res>8)"))
mtext("Periode (years)", side = 1, line= 2)
LLines(SPEC_ANALYSIS$MEAN_SPEC$SIM_full_down$spec, col = "#074893", lty = 3 )
LLines(SPEC_ANALYSIS$MEAN_SPEC$SIM_full_rec$spec, col = "#074893", lty = 3)
text(0.2, 2e2, "HadCM3 yearly res.", col = "#074893")
text(0.3, 4e0, "5y filter", col = "#074893")
text(0.3, 3e-2, "50y filter", col = "#074893")

axis(side = 1, at = c(0.002, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5), 
     labels = c(1/0.001, 1/0.005, 1/0.01, 1/0.02, 1/0.05, 1/0.2, 1/0.5))

LLines(SPEC_ANALYSIS$MEAN_SPEC$SIM_ds$spec, col = "#91002B")
LLines(SPEC_ANALYSIS$MEAN_SPEC$SIM_down_rec$spec, col = "#91002B", lty = 3)
text(0.02, 0.01, "HadCM3 down-sampled", col = "#91002B")
text(0.02, 0.001, "10y filter", col = "#91002B")

LLines(SPEC_ANALYSIS$MEAN_SPEC$Record$spec, col = "black", lw = 2)
text(0.02, 0.0005, "Records", col = "black")



legend("bottomleft", legend = c("HadCM3 yearly res.", "HadCM3 down-sampled to record res.", "Records"), 
       col = c("#074893","#91002B","black"), lwd = c(2,2,2))
dev.off()


