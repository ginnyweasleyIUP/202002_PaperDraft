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
  entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
  name = paste0("ENTITY", entity)
  if(length(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age) > 8){
    entities_spec = c(entities_spec, entity)
    diff_dt = mean(diff(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age), na.rm = T)
    record <- PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age,DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement,
                                         time.target = seq(from = head(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age, n = 1),
                                                           to = tail(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age, n = 1),
                                                           by = diff_dt))
    record = na.omit(record)
    data  = ts(data = record, 
               start = head(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age, n = 1), 
               end   = tail(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age, n = 1), 
               deltat = diff_dt)
    
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
  if(length(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age) > 8){
    entities_spec = c(entities_spec, entity)
    diff_dt = mean(diff(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age), na.rm = T)
    record <- PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age,DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$d18O,
                                         time.target = seq(from = head(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, n = 1),
                                                           to = tail(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, n = 1),
                                                           by = diff_dt))
    record = na.omit(record)
    data  = ts(data = record, 
               start = head(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, n = 1), 
               end   = tail(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, n = 1), 
               deltat = diff_dt)
    
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

# PaleoSpec::SpecMTM needs equally spaced data...

for(ii in 1:length_cave){
  entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
  site = DATA_past1000$CAVES$entity_info$site_id[ii]
  name = paste0("ENTITY", entity)
  if(length(DATA_past1000$CAVES$sim_data_raw[[paste0("CAVE", site)]]$ISOT) > 8){
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
## Plotting #####################################
#################################################

pdf(file = "Plots/Spectra/Global_Spectra_compare.pdf", width = 8, height = 5)
LPlot(SPEC_ANALYSIS$MEAN_SPEC$SIM_full$spec, col = "#074893", 
      ylim = c(0.000001,500), xlim = c(1/500, 5),
      xaxt = 'n',
      xlab = "",
      main = TeX("Mean Spectra from cave locations (res>8)"))
mtext("Periode (years)", side = 1, line= 2)
axis(side = 1, at = c(0.002, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5, 1.0, 2.0, 5.0), 
     labels = c(1/0.001, 1/0.005, 1/0.01, 1/0.02, 1/0.05, 1/0.2, 1/0.5, 1/1.0, 1/2.0, 1/5.0))
abline(v=0.5)
LLines(SPEC_ANALYSIS$MEAN_SPEC$SIM_ds$spec, col = "#91002B")
LLines(SPEC_ANALYSIS$MEAN_SPEC$Record$spec, col = "black")
legend("bottomleft", legend = c("HadCM3 yearly res.", "HadCM3 down-sampled to record res.", "Records"), 
       col = c("#074893","#91002B","black"), lwd = c(2,2,2))
dev.off()

#49, 110, 144, 148?, 203, 209, 277, 278



