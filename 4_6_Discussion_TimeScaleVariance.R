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
                            Sim_down_long = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                            Power_Law_1_short = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                            Power_Law_1_long = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                            Power_Law_2_short = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                            Power_Law_2_long = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                            WhiteNoise_short = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                            WhiteNoise_long = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])))

for(entity_number in 1:length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])){
  entity = DATA_past1000$CAVES$entity_info$entity_id[mask_var][entity_number]
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  TS_rec = zoo(x = DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_measurement, order.by = DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$interp_age)
  TS_sim_full = zoo(x = DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ITPC_a, order.by = seq(1950-DATA_past1000$time[1],1950-DATA_past1000$time[2]+1, by = -1))
  TS_sim_down = zoo(x = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$ITPC_a, order.by = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$interp_age)
  
  TS_PL1 = PaleoSpec::SimPowerlaw(1, 1110)
  TS_PL2 = PaleoSpec::SimPowerlaw(0.5, 1110)
  TS_WN = runif(1110)
  
  #Record
  if(is.na(tsc_dep_var(timser = TS_rec, tsc.in = c(10,50))$var.tsc)){
    DISCUSSION$VARIANCE$Rec_short[entity_number] = NA
  }else{
    DISCUSSION$VARIANCE$Rec_short[entity_number] = tsc_dep_var(timser = TS_rec, tsc.in = c(10,50))$var.tsc/tsc_dep_var(timser = TS_rec, tsc.in = c(10,50))$var.tot
  }
  if(is.na(tsc_dep_var(timser = TS_rec, tsc.in = c(50,200))$var.tsc)){
    DISCUSSION$VARIANCE$Rec_long[entity_number] = NA
  }else{
    DISCUSSION$VARIANCE$Rec_long[entity_number] = tsc_dep_var(timser = TS_rec, tsc.in = c(50,200))$var.tsc/tsc_dep_var(timser = TS_rec, tsc.in = c(50,200))$var.tot
  }
  
  
  #SIMfull
  DISCUSSION$VARIANCE$Sim_full_short[entity_number] = tsc_dep_var(timser = TS_sim_full, tsc.in = c(10,50))$var.tsc/tsc_dep_var(timser = TS_sim_full, tsc.in = c(10,50))$var.tot
  DISCUSSION$VARIANCE$Sim_full_long[entity_number] = tsc_dep_var(timser = TS_sim_full, tsc.in = c(50,200))$var.tsc/tsc_dep_var(timser = TS_sim_full, tsc.in = c(50,200))$var.tot
  
  if(is.na(tsc_dep_var(timser = TS_sim_down, tsc.in = c(10,50))$var.tsc)){
    DISCUSSION$VARIANCE$Sim_down_short[entity_number] = NA
  }else{
    DISCUSSION$VARIANCE$Sim_down_short[entity_number] =tsc_dep_var(timser = TS_sim_down, tsc.in = c(10,50))$var.tsc/tsc_dep_var(timser = TS_sim_down, tsc.in = c(10,50))$var.tot
  }
  if(is.na(tsc_dep_var(timser = TS_sim_down, tsc.in = c(50,200))$var.tsc)){
    DISCUSSION$VARIANCE$Sim_down_long[entity_number] = NA
  }else{
    DISCUSSION$VARIANCE$Sim_down_long[entity_number] = tsc_dep_var(timser = TS_sim_down, tsc.in = c(50,200))$var.tsc/tsc_dep_var(timser = TS_sim_down, tsc.in = c(50,200))$var.tot
  }
  
  
  DISCUSSION$VARIANCE$Power_Law_1_short[entity_number] = tsc_dep_var(timser = TS_PL1, tsc.in = c(10,50))$var.tsc/tsc_dep_var(timser = TS_PL1, tsc.in = c(10,50))$var.tot
  DISCUSSION$VARIANCE$Power_Law_1_long[entity_number] = tsc_dep_var(timser = TS_PL1, tsc.in = c(50,200))$var.tsc/tsc_dep_var(timser = TS_PL1, tsc.in = c(50,200))$var.tot
  DISCUSSION$VARIANCE$Power_Law_2_short[entity_number] = tsc_dep_var(timser = TS_PL2, tsc.in = c(10,50))$var.tsc/tsc_dep_var(timser = TS_PL2, tsc.in = c(10,50))$var.tot
  DISCUSSION$VARIANCE$Power_Law_2_long[entity_number] = tsc_dep_var(timser = TS_PL2, tsc.in = c(50,200))$var.tsc/tsc_dep_var(timser = TS_PL2, tsc.in = c(50,200))$var.tot
  DISCUSSION$VARIANCE$WhiteNoise_short[entity_number] = tsc_dep_var(timser = TS_WN, tsc.in = c(10,50))$var.tsc/tsc_dep_var(timser = TS_WN, tsc.in = c(10,50))$var.tot
  DISCUSSION$VARIANCE$WhiteNoise_long[entity_number] = tsc_dep_var(timser = TS_WN, tsc.in = c(50,200))$var.tsc/tsc_dep_var(timser = TS_WN, tsc.in = c(50,200))$var.tot
  
}

pdf(file = paste0("Plots/Discussion/Variance_TimeScaleDep.pdf"), width = 1.3*6, height = 3*PLOTTING_VARIABLES$HEIGHT/1.5)
par(mfrow=c(3,1),oma = c(1,3,0,0) + 0.1,mar = c(3,1,0,1) + 0.1)
hist(log10(DISCUSSION$VARIANCE$Rec_short), 
     breaks = 9, border = "white", prob = TRUE, 
     ylim = c(0,5), xlim = c(-2,0), xlab = "",xaxt = 'n', ylab = "",
     main = "", cex.main = 1.5, cex.axis = 1.5)
axis(side = 1, at = c(log10(0.001), log10(0.01), log10(0.1), 0), 
     labels = c(0.001,0.01, 0.1, 1), cex.axis = 1.5)
lines(density(log10(DISCUSSION$VARIANCE$Rec_short), na.rm = T),
      lwd = 2, col = "black")
lines(density(log10(DISCUSSION$VARIANCE$Rec_long), na.rm = T),
      lwd = 2, col = "#B2182B")
mtext(text = "density",side = 2,line = 2.5, cex = 1.5)
text(log10(0.1), 3, "long (50-200y)", col = "#B2182B", cex = 1.5)
text(log10(0.1), 2, "short (10-50y)", col = "black", cex = 1.5)
mtext(text = "Record", side = 3, line = -2, adj = 1,col = "black", cex = 1.5, at = log10(1))
mtext(text = "(a)", side = 3, line = -2, adj = 0,col = "black", cex = 1.5, at = log10(0.001))

hist(log10(DISCUSSION$VARIANCE$Sim_full_short), 
     breaks = 9, border = "white", prob = TRUE, 
     ylim = c(0,5), xlim = c(-2,0), xlab = "",xaxt = 'n', ylab = "",
     main = "", cex.main = 1.5, cex.axis = 1.5)
axis(side = 1, at = c(log10(0.001), log10(0.01), log10(0.1), 0), 
     labels = c(0.001,0.01, 0.1, 1), cex.axis = 1.5)
lines(density(log10(DISCUSSION$VARIANCE$Sim_full_short), na.rm = T),
      lwd = 2, col = "black")
lines(density(log10(DISCUSSION$VARIANCE$Sim_full_long), na.rm = T),
      lwd = 2, col = "#B2182B")
mtext(text = "density",side = 2,line = 2.5, cex = 1.5)
mtext(text = "Sim full", side = 3, line = -2, adj = 1,col = "black", cex = 1.5, at = log10(1))
mtext(text = "(b)", side = 3, line = -2, adj = 0,col = "black", cex = 1.5, at = log10(0.001))

hist(log10(DISCUSSION$VARIANCE$Sim_down_short), 
     breaks = 9, border = "white", prob = TRUE, 
     ylim = c(0,5), xlim = c(-2,0), xlab = "",xaxt = 'n', ylab = "",
     main = "", cex.main = 1.5, cex.axis = 1.5)
axis(side = 1, at = c(log10(0.001), log10(0.01), log10(0.1), 0), 
     labels = c(0.001,0.01, 0.1, 1), cex.axis = 1.5)
lines(density(log10(DISCUSSION$VARIANCE$Sim_down_short), na.rm = T),
      lwd = 2, col = "black")
lines(density(log10(DISCUSSION$VARIANCE$Sim_down_long), na.rm = T),
      lwd = 2, col = "#B2182B")
mtext(text = "density",side = 2,line = 2.5, cex = 1.5)
mtext(text = "Sim down", side = 3, line = -2, adj = 1,col = "black", cex = 1.5, at = log10(0.2))
mtext(text = "var.tsc(d18Opw)/var.tot(d18Opw)",side = 1,line = 2.7, cex = 1.5)
mtext(text = "(c)", side = 3, line = -2, adj = 0,col = "black", cex = 1.5, at = log10(0.001))


dev.off()

pdf(file = paste0("Plots/Discussion/Variance_TimeScaleDep_Example.pdf"), width = 1.3*6, height = 3*PLOTTING_VARIABLES$HEIGHT/1.5)
par(mfrow=c(3,1),oma = c(1,3,0,0) + 0.1,mar = c(3,1,0,1) + 0.1)
hist(log10(DISCUSSION$VARIANCE$Power_Law_1_short), 
     breaks = 9, border = "white", prob = TRUE, 
     ylim = c(0,5), xlim = c(-2,0), xlab = "",xaxt = 'n', ylab = "",
     main = "", cex.main = 1.5, cex.axis = 1.5)
axis(side = 1, at = c(log10(0.001), log10(0.01), log10(0.1), 0), 
     labels = c(0.001,0.01, 0.1, 1), cex.axis = 1.5)
lines(density(log10(DISCUSSION$VARIANCE$Power_Law_1_short), na.rm = T),
      lwd = 2, col = "black")
lines(density(log10(DISCUSSION$VARIANCE$Power_Law_1_long), na.rm = T),
      lwd = 2, col = "#B2182B")
mtext(text = "density",side = 2,line = 2.5, cex = 1.5)
text(log10(0.1), 3, "long (50-200y)", col = "#B2182B", cex = 1.5)
text(log10(0.1), 2, "short (10-50y)", col = "black", cex = 1.5)
mtext(text = "Power Law beta=1", side = 3, line = -2, adj = 1,col = "black", cex = 1.5, at = log10(1))
mtext(text = "(a)", side = 3, line = -2, adj = 0,col = "black", cex = 1.5, at = log10(0.001))

hist(log10(DISCUSSION$VARIANCE$Power_Law_2_long), 
     breaks = 9, border = "white", prob = TRUE, 
     ylim = c(0,5), xlim = c(-2,0), xlab = "",xaxt = 'n', ylab = "",
     main = "", cex.main = 1.5, cex.axis = 1.5)
axis(side = 1, at = c(log10(0.001), log10(0.01), log10(0.1), 0), 
     labels = c(0.001,0.01, 0.1, 1), cex.axis = 1.5)
lines(density(log10(DISCUSSION$VARIANCE$Power_Law_2_short), na.rm = T),
      lwd = 2, col = "black")
lines(density(log10(DISCUSSION$VARIANCE$Power_Law_2_long), na.rm = T),
      lwd = 2, col = "#B2182B")
mtext(text = "density",side = 2,line = 2.5, cex = 1.5)
mtext(text = "Power Law beta=0.5", side = 3, line = -2, adj = 1,col = "black", cex = 1.5, at = log10(1))
mtext(text = "(b)", side = 3, line = -2, adj = 0,col = "black", cex = 1.5, at = log10(0.001))

hist(log10(DISCUSSION$VARIANCE$WhiteNoise_short), 
     breaks = 9, border = "white", prob = TRUE, 
     ylim = c(0,5), xlim = c(-2,0), xlab = "",xaxt = 'n', ylab = "",
     main = "", cex.main = 1.5, cex.axis = 1.5)
axis(side = 1, at = c(log10(0.001), log10(0.01), log10(0.1), 0), 
     labels = c(0.001,0.01, 0.1, 1), cex.axis = 1.5)
lines(density(log10(DISCUSSION$VARIANCE$WhiteNoise_short), na.rm = T),
      lwd = 2, col = "black")
lines(density(log10(DISCUSSION$VARIANCE$WhiteNoise_long), na.rm = T),
      lwd = 2, col = "#B2182B")
mtext(text = "density",side = 2,line = 2.5, cex = 1.5)
mtext(text = "White Noise", side = 3, line = -2, adj = 1,col = "black", cex = 1.5, at = log10(1))
mtext(text = "var.tsc(d18Opw)/var.tot(d18Opw)",side = 1,line = 2.7, cex = 1.5)
mtext(text = "(c)", side = 3, line = -2, adj = 0,col = "black", cex = 1.5, at = log10(0.001))

dev.off()
