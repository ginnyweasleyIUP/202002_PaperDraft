#################################################
## Paper Figure SPECTRUM ########################
#################################################

library(zoo)
library(PaleoSpec)
library(nest)
library(latex2exp)

#################################################

## CALC

# We need 
# [ ] full yearly spectrum    TEMP, PREC, TEMP (prec weighted), ISOT, ITPC
# [ ] downsamples spectrum    TEMP, PREC, TEMP (prec weighted), ISOT, ITPC
# [ ] record spectrum
# [ ] gridbox weighing

ANALYSIS$SPECTRA <- list()
ANALYSIS$SPECTRA$RECORDS <- list()
ANALYSIS$SPECTRA$SIM_ds <- list()
ANALYSIS$SPECTRA$SIM_full <- list()

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
    
    SPEC_ANALYSIS$RECORDS[[name]] = SpecMTM(data)
    
  }
  
}

SPEC_ANALYSIS$GENERAL$entities_spec_rec <- entities_spec

SPEC_ANALYSIS$MEAN_SPEC[["Record"]] <- MeanSpectrum(SPEC_ANALYSIS$RECORDS)


## Weighing sollte ja wohl ohne Internet gehen...





## PLOT

pdf(file = "Plots/Paper/Paper_Plot_7_Spectra.pdf", width = PLOTTING_VARIABLES$WIDTH/2, height = PLOTTING_VARIABLES$HEIGHT/2*1.2)
LPlot(SPEC_ANALYSIS$MEAN_SPEC$SIM_full$spec, col = "#074893", 
      ylim = c(0.00001,1000), xlim = c(1/500, 0.5),
      xaxt = 'n',
      xlab = "")#,
#main = TeX("Mean Spectra from cave locations (res>8)"))
mtext("Periode (years)", side = 1, line= 2)
LLines(SPEC_ANALYSIS$MEAN_SPEC$SIM_full_down$spec, col = "#074893", lty = 3 )
LLines(SPEC_ANALYSIS$MEAN_SPEC$SIM_full_rec$spec, col = "#074893", lty = 3)
#text(0.2, 8e2, "HadCM3 yearly res.", col = "#074893")
#text(0.3, 3e2, "5y filter", col = "#074893")
#text(0.3, 1e2, "50y filter", col = "#074893")

axis(side = 1, at = c(0.002, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5), 
     labels = c(1/0.001, 1/0.005, 1/0.01, 1/0.02, 1/0.05, 1/0.2, 1/0.5))

LLines(SPEC_ANALYSIS$MEAN_SPEC$SIM_ds$spec, col = "#91002B")
LLines(SPEC_ANALYSIS$MEAN_SPEC$SIM_down_rec$spec, col = "#91002B", lty = 3)
#text(0.02, 0.01, "HadCM3 down-sampled", col = "#91002B")
#text(0.02, 0.001, "10y filter", col = "#91002B")

LLines(SPEC_ANALYSIS$MEAN_SPEC$Record$spec, col = "black", lw = 2)
#text(0.02, 0.0005, "Records", col = "black")



legend("bottomleft", legend = c("HadCM3 yearly res.", " ... 5y filter", " ... 50y filter", "HadCM3 down-sampled to record res.", " ... 10y filter", "Records"), 
       col = c("#074893","#074893","#074893","#91002B","#91002B","black"), lwd = c(1,1,1,1,1,2), lty = c(1,3,3,1,3,1), bty = "n")
dev.off()