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

## PLOT

COLZ <- c("#1A2254", "#0A4296", "#278BCE","#91002B", "#BD6B73", "black")

for(var in c("ISOT", "ITPC")){
  cairo_pdf(file = paste0("Plots/Discussion/TS_1-48_",var,".pdf"), width = 21, height = 29.7)
  par(mfrow=c(8,6),oma = c(1,3,0,0) + 0.1,mar = c(3,0,1,3) + 0.1)
  for(entity in DATA_past1000$CAVES$entity_info$entity_id[mask_mean][1:48]){
    site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
    yrange = range(c(range(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a, na.rm = T), 
                     range(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]][[paste0(var, "_a")]], na.rm = T)))
    yrange[1] = yrange[1]-1
    yrange[2] = yrange[2]+1
    yrange = c(-15,0)
    plot(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$interp_age, DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a, 
         col = adjustcolor(COLZ[1], alpha = 0.5), ylim = yrange, xlim = c(0,1150), type = "l",
          ylab = "",
          xlab = "", lwd = 2)
    #gaussdetr(zoo(x = Timeseries$pages2k$value[850:2000],order.by = Timeseries$pages2k$time[850:2000]), tsc.in = 100)$Xsmooth
    #gaussdetr()
    lines(gaussdetr(zoo(x=DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a, 
                        order.by = DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$interp_age), tsc.in = 100)$Xsmooth, col = COLZ[1], lw = 3)#,
    #main = TeX("Mean Spectra from cave locations (res>8)"))
    mtext("years BP", side = 1, line= 2, cex = 0.8)
    mtext("d18O in [‰]", side = 2, line= 2, cex = 0.8)
    lines(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$interp_age, 
          DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]][[paste0(var, "_a")]], 
          col = adjustcolor(COLZ[4], alpha.f = 0.5), lw = 2)
    lines(gaussdetr(zoo(x = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]][[paste0(var, "_a")]], 
                        order.by = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$interp_age), tsc.in = 100)$Xsmooth, 
          col = COLZ[4], lw = 2)
    
    if(var == "ITPC"){
      legend("bottomleft", legend = c(paste0("Entity ", entity," dweq"), "HadCM3 down-sampled"), 
             col = c(COLZ[1],COLZ[4]), lwd = c(2,2), lty = c(1,1), bty = "n", cex = cex_text)
    }else{
      legend("bottomleft", legend = c(paste0("Entity ", entity," dweq"), "HadCM3 yearly"), 
             col = c(COLZ[1],COLZ[4]), lwd = c(2,2), lty = c(1,1), bty = "n", cex = cex_text)
    }
    
  }
  
  
  
  dev.off()
  
  cairo_pdf(file = paste0("Plots/Discussion/TS_49-76_",var,".pdf"), width = 21, height = 29.7)
  par(mfrow=c(8,6),oma = c(1,3,0,0) + 0.1,mar = c(3,0,1,3) + 0.1)
  for(entity in DATA_past1000$CAVES$entity_info$entity_id[mask_mean][49:96]){
    site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
    yrange = range(c(range(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a, na.rm = T), 
                     range(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]][[paste0(var, "_a")]], na.rm = T)))
    yrange[1] = yrange[1]-1
    yrange[2] = yrange[2]+1
    yrange = c(-15,0)
    plot(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$interp_age, DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a, 
         col = adjustcolor(COLZ[1], alpha = 0.5), ylim = yrange, xlim = c(0,1150), type = "l",
         ylab = "",
         xlab = "", lwd = 2)
    #gaussdetr(zoo(x = Timeseries$pages2k$value[850:2000],order.by = Timeseries$pages2k$time[850:2000]), tsc.in = 100)$Xsmooth
    #gaussdetr()
    lines(gaussdetr(zoo(x=DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a, 
                        order.by = DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$interp_age), tsc.in = 100)$Xsmooth, col = COLZ[1], lw = 3)#,
    #main = TeX("Mean Spectra from cave locations (res>8)"))
    mtext("years BP", side = 1, line= 2, cex = 0.8)
    mtext("d18O in [‰]", side = 2, line= 2, cex = 0.8)
    lines(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$interp_age, DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]][[paste0(var, "_a")]], 
          col = adjustcolor(COLZ[4], alpha.f = 0.5), lw = 2)
    lines(gaussdetr(zoo(x = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]][[paste0(var, "_a")]], 
                        order.by = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$interp_age), tsc.in = 100)$Xsmooth, 
          col = COLZ[4], lw = 2)
    
    
    if(var == "ITPC"){
      legend("bottomleft", legend = c(paste0("Entity ", entity," dweq"), "HadCM3 down-sampled"), 
             col = c(COLZ[1],COLZ[4]), lwd = c(2,2), lty = c(1,1), bty = "n", cex = cex_text)
    }else{
      legend("bottomleft", legend = c(paste0("Entity ", entity," dweq"), "HadCM3 yearly"), 
             col = c(COLZ[1],COLZ[4]), lwd = c(2,2), lty = c(1,1), bty = "n", cex = cex_text)
    }
  }
  
  
  
  dev.off()
  
  cairo_pdf(file = paste0("Plots/Discussion/TS_77-144_",var,".pdf"), width = 21, height = 29.7)
  par(mfrow=c(8,6),oma = c(1,3,0,0) + 0.1,mar = c(3,0,1,3) + 0.1)
  for(entity in DATA_past1000$CAVES$entity_info$entity_id[mask_mean][97:108]){
    site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
    yrange = range(c(range(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a, na.rm = T), 
                     range(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]][[paste0(var, "_a")]], na.rm = T)))
    yrange[1] = yrange[1]-1
    yrange[2] = yrange[2]+1
    yrange = c(-15,0)
    plot(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$interp_age, DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a, 
         col = adjustcolor(COLZ[1], alpha = 0.5), ylim = yrange, xlim = c(0,1150), type = "l",
         ylab = "",
         xlab = "", lwd = 2)
    #gaussdetr(zoo(x = Timeseries$pages2k$value[850:2000],order.by = Timeseries$pages2k$time[850:2000]), tsc.in = 100)$Xsmooth
    #gaussdetr()
    lines(gaussdetr(zoo(x=DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a, 
                        order.by = DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$interp_age), tsc.in = 100)$Xsmooth, col = COLZ[1], lw = 3)#,
    #main = TeX("Mean Spectra from cave locations (res>8)"))
    mtext("years BP", side = 1, line= 2, cex = 0.8)
    mtext("d18O in [‰]", side = 2, line= 2, cex = 0.8)
    lines(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$interp_age, DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]][[paste0(var, "_a")]], 
          col = adjustcolor(COLZ[4], alpha.f = 0.5), lw = 2)
    lines(gaussdetr(zoo(x = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]][[paste0(var, "_a")]], 
                        order.by = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$interp_age), tsc.in = 100)$Xsmooth, 
          col = COLZ[4], lw = 2)
    
    
    if(var == "ITPC"){
      legend("bottomleft", legend = c(paste0("Entity ", entity," dweq"), "HadCM3 down-sampled"), 
             col = c(COLZ[1],COLZ[4]), lwd = c(2,2), lty = c(1,1), bty = "n", cex = cex_text)
    }else{
      legend("bottomleft", legend = c(paste0("Entity ", entity," dweq"), "HadCM3 yearly"), 
             col = c(COLZ[1],COLZ[4]), lwd = c(2,2), lty = c(1,1), bty = "n", cex = cex_text)
    }
  }
  
  
  
  dev.off()
  
}


