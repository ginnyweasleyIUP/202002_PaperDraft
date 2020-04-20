#################################################
## DISCUSSION Time Scale Variance ###############
#################################################

library(plyr)
library(dplyr)
library(tidyverse)
library(PaleoSpec)
library(nest)

DISCUSSION$VARIANCE <- list(Rec_short = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                            Rec_long = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                            Sim_full_short = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                            Sim_full_long = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                            Sim_down_short = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                            Sim_down_long = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])))

for(entity_number in 1:length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])){
  entity = DATA_past1000$CAVES$entity_info$entity_id[mask_var][entity_number]
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  TS_rec = zoo(x = DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_measurement, order.by = DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$interp_age)
  TS_sim_full = zoo(x = DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ITPC_a, order.by = seq(1950-DATA_past1000$time[1],1950-DATA_past1000$time[2]+1, by = -1))
  TS_sim_down = zoo(x = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$ITPC_a, order.by = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$interp_age)

  DISCUSSION$VARIANCE$Rec_short[entity_number] = tsc_dep_var(timser = TS_rec, tsc.in = c(10,50))$var.tsc
  DISCUSSION$VARIANCE$Rec_long[entity_number] = tsc_dep_var(timser = TS_rec, tsc.in = c(50,200))$var.tsc
  DISCUSSION$VARIANCE$Sim_full_short[entity_number] = tsc_dep_var(timser = TS_sim_full, tsc.in = c(10,50))$var.tsc
  DISCUSSION$VARIANCE$Sim_full_long[entity_number] = tsc_dep_var(timser = TS_sim_full, tsc.in = c(50,200))$var.tsc
  DISCUSSION$VARIANCE$Sim_down_short[entity_number] = tsc_dep_var(timser = TS_sim_down, tsc.in = c(10,50))$var.tsc
  DISCUSSION$VARIANCE$Sim_down_long[entity_number] = tsc_dep_var(timser = TS_sim_down, tsc.in = c(50,200))$var.tsc
  
}

pdf(file = paste0("Plots/Discussion/Variance_TimeScaleDep.pdf"), width = 1.3*6, height = 3*PLOTTING_VARIABLES$HEIGHT/1.5)
par(mfrow=c(3,1),oma = c(1,3,0,0) + 0.1,mar = c(3,1,0,1) + 0.1)
hist(log10(DISCUSSION$VARIANCE$Rec_short), 
     breaks = 9, border = "white", prob = TRUE, 
     ylim = c(0,2), xlim = c(-3,1), xlab = "",xaxt = 'n', ylab = "",
     main = "", cex.main = 1.5, cex.axis = 1.5)
axis(side = 1, at = c(log10(0.001), log10(0.01), log10(0.1), 0, log10(10), log10(100), log10(1000)), 
     labels = c(0.001,0.01, 0.1, 1, 10, 100, 1000), cex.axis = 1.5)
lines(density(log10(DISCUSSION$VARIANCE$Rec_short), na.rm = T),
      lwd = 2, col = "black")
lines(density(log10(DISCUSSION$VARIANCE$Rec_long), na.rm = T),
      lwd = 2, col = "#B2182B")
mtext(text = "density",side = 2,line = 2.5, cex = 1.5)
text(log10(1), 1.3, "long (50-200y)", col = "#B2182B", cex = 1.5)
text(log10(1), 1.0, "short (10-50y)", col = "black", cex = 1.5)
mtext(text = "Record", side = 3, line = -2, adj = 1,col = "black", cex = 1.5, at = log10(10))

hist(log10(DISCUSSION$VARIANCE$Sim_full_short), 
     breaks = 9, border = "white", prob = TRUE, 
     ylim = c(0,2), xlim = c(-3,1), xlab = "",xaxt = 'n', ylab = "",
     main = "", cex.main = 1.5, cex.axis = 1.5)
axis(side = 1, at = c(log10(0.001), log10(0.01), log10(0.1), 0, log10(10), log10(100), log10(1000)), 
     labels = c(0.001,0.01, 0.1, 1, 10, 100, 1000), cex.axis = 1.5)
lines(density(log10(DISCUSSION$VARIANCE$Sim_full_short), na.rm = T),
      lwd = 2, col = "black")
lines(density(log10(DISCUSSION$VARIANCE$Sim_full_long), na.rm = T),
      lwd = 2, col = "#B2182B")
mtext(text = "density",side = 2,line = 2.5, cex = 1.5)
mtext(text = "Sim full", side = 3, line = -2, adj = 1,col = "black", cex = 1.5, at = log10(10))

hist(log10(DISCUSSION$VARIANCE$Sim_down_short), 
     breaks = 9, border = "white", prob = TRUE, 
     ylim = c(0,2), xlim = c(-3,1), xlab = "",xaxt = 'n', ylab = "",
     main = "", cex.main = 1.5, cex.axis = 1.5)
axis(side = 1, at = c(log10(0.001), log10(0.01), log10(0.1), 0, log10(10), log10(100), log10(1000)), 
     labels = c(0.001,0.01, 0.1, 1, 10, 100, 1000), cex.axis = 1.5)
lines(density(log10(DISCUSSION$VARIANCE$Sim_down_short), na.rm = T),
      lwd = 2, col = "black")
lines(density(log10(DISCUSSION$VARIANCE$Sim_down_long), na.rm = T),
      lwd = 2, col = "#B2182B")
mtext(text = "density",side = 2,line = 2.5, cex = 1.5)
mtext(text = "Sim down", side = 3, line = -2, adj = 1,col = "black", cex = 1.5, at = log10(10))
mtext(text = "var.tsc(d18Opw)",side = 1,line = 2.7, cex = 1.5)


dev.off()
