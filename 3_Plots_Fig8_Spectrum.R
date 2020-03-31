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

filter <- list(
  ISOT = c(2.0, 16.0, 8.0),
  ITPC = c(3.0, 9.0, 3.0)
)

COLZ <- c("#1A2254", "#0A4296", "#278BCE","#91002B", "#BD6B73", "black")

for(var in c("ISOT", "ITPC")){
  pdf(file = paste0("Plots/Paper_Plot_7_Spectra_",var,".pdf"), width = PLOTTING_VARIABLES$WIDTH*2/3, height = PLOTTING_VARIABLES$HEIGHT*1.2*2/3)
  #png(file = paste0("Plots/Paper_Plot_7_Spectra_",var,".png"), width = 70*PLOTTING_VARIABLES$WIDTH*2/3, height = 70*PLOTTING_VARIABLES$HEIGHT*1.2*2/3)
  LPlot(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec, col = adjustcolor(COLZ[1], alpha = 0.5), 
        ylim = c(0.00005,1000), xlim = c(1/300, 0.5),
        ylab = "",
        xaxt = 'n',
        yaxt = "n",
        xlab = "", lwd = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec), col = COLZ[1], lw = 3)#,
  #main = TeX("Mean Spectra from cave locations (res>8)"))
  mtext("Period (y)", side = 1, line= 2)
  mtext("Power spectral sensity", side = 2, line= 2)
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_down_",var)]]$spec, col = COLZ[2], lty = 3, lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_down_",var)]]$spec), col = COLZ[2], lw = 2)
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_rec_",var)]]$spec, col = COLZ[3], lty = 3, lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_rec_",var)]]$spec), col = COLZ[3], lw = 2)
  #text(0.2, 8e2, "HadCM3 yearly res.", col = "#074893")
  #text(0.3, 3e2, "5y filter", col = "#074893")
  #text(0.3, 1e2, "50y filter", col = "#074893")
  
  axis(side = 1, at = c(0.002,0.0033333333, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5), 
       labels = c(1/0.002, 300, 1/0.005, 1/0.01, 1/0.02, 1/0.05, 1/0.2, 1/0.5))
  axis(side = 2, at = c(1e-3, 1e-1, 1e1, 1e3), 
       labels = c(1e-3, 0.1, 10, 1000))
  
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_ds_",var)]]$spec, col = adjustcolor(COLZ[4], alpha = 0.5), lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_ds_",var)]]$spec), col = COLZ[4], lw = 2)
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_down_rec_",var)]]$spec, col = COLZ[5], lty = 3, lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_down_rec_",var)]]$spec), col = COLZ[5], lw = 2)
  #text(0.02, 0.01, "HadCM3 down-sampled", col = COLZ[4])
  #text(0.02, 0.001, "10y filter", col = COLZ[4])
  
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC$Record$spec, col = adjustcolor(COLZ[6], alpha = 0.5), lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC$Record$spec), col = COLZ[6], lw = 3)
  #text(0.02, 0.0005, "Records", col = "black")
  legend("bottomleft", legend = c("HadCM3 (full)", paste0(" ...  ",filter[[var]][1],"y filter"), paste0(" ...  ",filter[[var]][2],"y filter"), "HadCM3 (down-sampled)", paste0(" ...  ",filter[[var]][3],"y filter"), "Records"), 
         col = c(COLZ[1],COLZ[2],COLZ[3],COLZ[4],COLZ[5],COLZ[6]), lwd = c(1,1,1,1,1,2), lty = c(1,1,1,1,1,1), bty = "n", cex = cex_text)
  
  if(var == "ISOT"){
    text(0.5, 500, TeX("HadCM3: $\\delta^{18}O$ in prec"), adj = 1)
  }
  if(var == "ITPC"){
    text(0.5, 500, TeX("HadCM3: prec-weighted $\\delta^{18}O$"), adj = 1)
  }
  
  
  dev.off()
  
}

for(var in c("ISOT", "ITPC")){
  cairo_pdf(file = paste0("Plots/Appendix/Paper_Plot_7_Spectra_Compare",var,".pdf"), width = PLOTTING_VARIABLES$WIDTH*2/3, height = PLOTTING_VARIABLES$HEIGHT*1.2*2/3)
  #par(mfrow=c(1,2))
  cex_text = 1
  cex_axis = 0.6
  cex_axis_text = 0.7
  par(mar=c(1,0,0,1),oma=c(5,5,2,3),xaxs="i",yaxs="i",cex=1,lwd=2)
  layout(matrix(c(1,2,3,4,4,4), 3, 2, byrow = F))
  
  plot(DATA_past1000$CAVES$record_data$ENTITY240$interp_age, DATA_past1000$CAVES$record_data$ENTITY240$d18O_dw_eq_a,
       col = "black", type = "l", ylim = c(-10,-5.5), xlim = c(0,1100), xaxt = "n", yaxt = "n",ylab = "", xlab = "")
  axis(2,at=seq(-10,-6,by=1),labels=FALSE,col=COLZ[6])
  mtext(side=2,at=seq(-10,-6,by=1),line = 1, seq(-10,-6,by=1),las=1,col=COLZ[6], cex = cex_axis)
  text(1000, -6, TeX("Example: eID 240 (Bunker Cave)"), col = COLZ[6], adj = 1, cex = cex_text)
  mtext("a)", side = 3, line = -1.3, adj = 0.025, cex = cex_axis_text)
  
  plot(seq(from = 1950 - DATA_past1000$time[[1]], to = 1950 - DATA_past1000$time[[2]]+1, by = -1), DATA_past1000$CAVES$sim_data_yearly$CAVE117[[paste0(var,"_a")]], 
       col = COLZ[1], type = "l", ylim = c(-10,-5.5), xlim = c(0,1100), xaxt = "n", yaxt = "n",ylab = "", xlab = "")
  axis(2,at=seq(-10,-6,by=1),labels=FALSE,col=COLZ[1])
  mtext(side=2,at=seq(-10,-6,by=1),line = 1, seq(-10,-6,by=1),col=COLZ[1], cex = cex_axis)
  text(1000, -6, TeX("HadCM3 at eID 240"), col = COLZ[1], adj = 1, cex = cex_text)
  mtext(TeX("$\\delta^{18}O$"),side=2, at = -8, line = 2.5, las = 1, col = COLZ[1], cex = cex_axis_text)
  mtext(text = "[‰]", side = 2, at = -9, line = 2.5, las = 1, col = COLZ[1], cex = cex_axis_text)
  mtext("b)", side = 3, line = -1.3, adj = 0.025, cex = cex_axis_text)
  
  plot(DATA_past1000$CAVES$sim_data_downsampled$ENTITY240$interp_age, DATA_past1000$CAVES$sim_data_downsampled$ENTITY240[[paste0(var,"_a")]],
        col = COLZ[4], type = "l", ylim = c(-10,-5.5), xlim = c(0,1100), xaxt = "n", yaxt = "n",ylab = "", xlab = "")
  axis(2,at=seq(-10,-6,by=1),labels=FALSE,col=COLZ[4])
  mtext(side=2,at=seq(-10,-6,by=1),line = 1, seq(-10,-6,by=1),las=1,col=COLZ[4], cex = cex_axis)
  text(1000, -6, TeX("HadCM3 at eID 240 down-sampled"), col = COLZ[4], adj = 1, cex = cex_text)
  mtext("c)", side = 3, line = -1.3, adj = 0.025, cex = cex_axis_text)
  
  
  
  axis(1,at=seq(0,1000,by=200),labels=FALSE,col="black")
  mtext(side=1,at=seq(0,1000,by=200),line = 1, seq(0,1000,by=200),las=1,col="black", cex = cex_axis)
  mtext("years BP", side = 1, line = 3, at = 500, cex = cex_axis_text)
  
  #png(file = paste0("Plots/Paper_Plot_7_Spectra_",var,".png"), width = 70*PLOTTING_VARIABLES$WIDTH*2/3, height = 70*PLOTTING_VARIABLES$HEIGHT*1.2*2/3)
  LPlot(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec, col = adjustcolor(COLZ[1], alpha = 0.5), 
        ylim = c(0.00005,1000), xlim = c(1/300, 0.5),
        ylab = "",
        xaxt = 'n',
        yaxt = "n",
        xlab = "", lwd = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec), col = COLZ[1], lw = 3)#,
  #main = TeX("Mean Spectra from cave locations (res>8)"))
  mtext("Period (y)", side = 1, line= 3, cex = cex_axis_text)
  mtext("PSD", side = 4, line= 1.5, at = 0.5, las = 1, cex = cex_axis_text)
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_down_",var)]]$spec, col = COLZ[2], lty = 3, lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_down_",var)]]$spec), col = COLZ[2], lw = 2)
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_rec_",var)]]$spec, col = COLZ[3], lty = 3, lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_rec_",var)]]$spec), col = COLZ[3], lw = 2)
  #text(0.2, 8e2, "HadCM3 yearly res.", col = "#074893")
  #text(0.3, 3e2, "5y filter", col = "#074893")
  #text(0.3, 1e2, "50y filter", col = "#074893")
  
  axis(side = 1, at = c(0.002,0.0033333333, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5), labels = FALSE)
  mtext(side=1,at=c(0.0033333333, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5),line = 1, 
        c(300, 1/0.005, 1/0.01, 1/0.02, 1/0.05, 1/0.2, 1/0.5),las=1,col="black", cex = cex_axis)
  axis(side = 4, at = c(1e-3, 1e-1, 1e1, 1e3), labels = FALSE)
  mtext(side=4,at=c(1e-3, 1e-1, 1e1, 1e3),line = 1, c(1e-3, 0.1, 10, 1000),las=1,col="black", cex = cex_axis)
  
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_ds_",var)]]$spec, col = adjustcolor(COLZ[4], alpha = 0.5), lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_ds_",var)]]$spec), col = COLZ[4], lw = 2)
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_down_rec_",var)]]$spec, col = COLZ[5], lty = 3, lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_down_rec_",var)]]$spec), col = COLZ[5], lw = 2)
  #text(0.02, 0.01, "HadCM3 down-sampled", col = COLZ[4])
  #text(0.02, 0.001, "10y filter", col = COLZ[4])
  
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC$Record$spec, col = adjustcolor(COLZ[6], alpha = 0.5), lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC$Record$spec), col = COLZ[6], lw = 3)
  #text(0.02, 0.0005, "Records", col = "black")
  
  legend("bottomleft", legend = c("HadCM3 (full)", paste0(" ...  ",filter[[var]][1],"y filter"), paste0(" ...  ",filter[[var]][2],"y filter"), "HadCM3 (down-sampled)", paste0(" ...  ",filter[[var]][3],"y filter"), "Records"), 
         col = c(COLZ[1],COLZ[2],COLZ[3],COLZ[4],COLZ[5],COLZ[6]), lwd = c(1,1,1,1,1,2), lty = c(1,1,1,1,1,1), bty = "n", cex = cex_text)
  
  if(var == "ISOT"){
    text(0.4, 500, TeX("Spectrum HadCM3: $\\delta^{18}O$"), adj = 1, cex = cex_text)
  }
  if(var == "ITPC"){
    text(0.4, 500, TeX("Spectrum HadCM3: $\\delta^{18}O_{pw}$"), adj = 1, cex = cex_text)
  }
  mtext("d)", side = 3, line = -1.3, adj = 0.025, cex = cex_axis_text)
  
  
  dev.off()
  
}

#################################################
## WITH FILTER ##################################
#################################################

source('Functions/Filter/EASY_Sensor_WM4.R')
source('Functions/Filter/filter_function3.R')

entity = 240
site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
name = paste0("ENTITY", entity)

for(var in c("ISOT", "ITPC")){
  #Filter 
  #full -> down
  Results <- easy_sensor_wm4(1.0, na.omit(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]][[paste0(var,"_a")]]), filter[[var]][1])
  TS_full_down  = ts(data = rev(Results), start = 1950-DATA_past1000$time[2]-29, end   = 1950-DATA_past1000$time[1], deltat = 1)
  #full-> rec
  Results <- easy_sensor_wm4(1.0, na.omit(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]][[paste0(var,"_a")]]), filter[[var]][2])
  TS_full_rec  = ts(data = rev(Results), start = 1950-DATA_past1000$time[2]-144, end   = 1950-DATA_past1000$time[1], deltat = 1)
  #down->rec
  start_ts = ceiling(head(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, n = 1))
  length = length(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age)
  stop_ts = floor(tail(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, n = 1))
  diff = floor((stop_ts-start_ts)/length)
  if(diff<1){diff = 1}
  record <- PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age,
                                       DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]][[paste0(var,"_a")]],
                                       time.target = seq(from = start_ts, to = stop_ts, by = diff))
  Results <- easy_sensor_wm4(diff, na.omit(rev(record)), filter[[var]][3])
  TS_down_rec = ts(data = Results, start = start_ts-length(Results)+length(record), deltat = diff)
  
  cairo_pdf(file = paste0("Plots/Appendix/Paper_Plot_7_Spectra_Compare_Filter_",var,".pdf"), width = PLOTTING_VARIABLES$WIDTH*2/3, height = PLOTTING_VARIABLES$HEIGHT*1.2*2/3)
  #par(mfrow=c(1,2))
  cex_text = 1
  cex_axis = 0.6
  cex_axis_text = 0.7
  par(mar=c(1,0,0,1),oma=c(3,5,2,3),xaxs="i",yaxs="i",cex=1,lwd=2)
  layout(matrix(c(1,2,3,4,4,4), 3, 2, byrow = F))
  
  plot(DATA_past1000$CAVES$record_data$ENTITY240$interp_age, DATA_past1000$CAVES$record_data$ENTITY240$d18O_dw_eq_a,
       col = "black", type = "l", ylim = c(-10,-5.5), xlim = c(0,1100), xaxt = "n", yaxt = "n",ylab = "", xlab = "")
  axis(2,at=seq(-10,-6,by=1),labels=FALSE,col=COLZ[6])
  mtext(side=2,at=seq(-10,-6,by=1),line = 1, seq(-10,-6,by=1),las=1,col=COLZ[6], cex = cex_axis)
  text(1050, -6, TeX("Example: eID 240 (Bunker Cave)"), col = COLZ[6], adj = 1, cex = cex_text)
  mtext("a)", side = 3, line = -1.3, adj = 0.025, cex = cex_axis_text)
  
  plot(seq(from = 1950 - DATA_past1000$time[[1]], to = 1950 - DATA_past1000$time[[2]]+1, by = -1), DATA_past1000$CAVES$sim_data_yearly$CAVE117[[paste0(var,"_a")]], 
       col = COLZ[1], type = "l", ylim = c(-10,-5.5), xlim = c(0,1100), xaxt = "n", yaxt = "n",ylab = "", xlab = "")
  lines(window(TS_full_down, start = 30, end = 1100), col = COLZ[2])
  lines(window(TS_full_rec, start = 30, end = 1030), col = COLZ[3])
  axis(2,at=seq(-10,-6,by=1),labels=FALSE,col=COLZ[1])
  mtext(side=2,at=seq(-10,-6,by=1),line = 1, seq(-10,-6,by=1),col=COLZ[1], cex = cex_axis)
  text(1050, -6, TeX("HadCM3 at eID 240"), col = COLZ[1], adj = 1, cex = cex_text)
  text(1050, -6.5, paste0(filter[[var]][1], "y filter"), col = COLZ[2], adj = 1, cex = cex_text)
  text(1050, -9.5, paste0(filter[[var]][2], "y filter"), col = COLZ[3], adj = 1, cex = cex_text)
  mtext(TeX("$\\delta^{18}O$"),side=2, at = -8, line = 2.5, las = 1, col = COLZ[1], cex = cex_axis_text)
  mtext(text = "[‰]", side = 2, at = -9, line = 2.5, las = 1, col = COLZ[1], cex = cex_axis_text)
  mtext("b)", side = 3, line = -1.3, adj = 0.025, cex = cex_axis_text)
  
  plot(DATA_past1000$CAVES$sim_data_downsampled$ENTITY240$interp_age, DATA_past1000$CAVES$sim_data_downsampled$ENTITY240[[paste0(var,"_a")]],
       col = COLZ[4], type = "l", ylim = c(-10,-5.5), xlim = c(0,1100), xaxt = "n", yaxt = "n",ylab = "", xlab = "")
  lines(window(TS_down_rec, start = 150, end = 1100), col = COLZ[5])
  axis(2,at=seq(-10,-6,by=1),labels=FALSE,col=COLZ[4])
  mtext(side=2,at=seq(-10,-6,by=1),line = 1, seq(-10,-6,by=1),las=1,col=COLZ[4], cex = cex_axis)
  text(1050, -6, TeX("HadCM3 at eID 240 down-sampled"), col = COLZ[4], adj = 1, cex = cex_text)
  text(1050, -6.5, paste0(filter[[var]][3], "y filter"), col = COLZ[5], adj = 1, cex = cex_text)
  mtext("c)", side = 3, line = -1.3, adj = 0.025, cex = cex_axis_text)
  
  
  
  axis(1,at=seq(0,1000,by=200),labels=FALSE,col="black")
  mtext(side=1,at=seq(0,1000,by=200),line = 1, seq(0,1000,by=200),las=1,col="black", cex = cex_axis)
  mtext("years BP", side = 1, line = 2, at = 500, cex = cex_axis_text)
  
  #png(file = paste0("Plots/Paper_Plot_7_Spectra_",var,".png"), width = 70*PLOTTING_VARIABLES$WIDTH*2/3, height = 70*PLOTTING_VARIABLES$HEIGHT*1.2*2/3)
  LPlot(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec, col = adjustcolor(COLZ[1], alpha = 0.5), 
        ylim = c(0.00005,1000), xlim = c(1/300, 0.5),
        ylab = "",
        xaxt = 'n',
        yaxt = "n",
        xlab = "", lwd = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec), col = COLZ[1], lw = 3)#,
  #main = TeX("Mean Spectra from cave locations (res>8)"))
  mtext("Period (y)", side = 1, line= 2, cex = cex_axis_text)
  mtext("PSD", side = 4, line= 1.5, at = 0.5, las = 1, cex = cex_axis_text)
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_down_",var)]]$spec, col = COLZ[2], lty = 3, lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_down_",var)]]$spec), col = COLZ[2], lw = 2)
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_rec_",var)]]$spec, col = COLZ[3], lty = 3, lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_rec_",var)]]$spec), col = COLZ[3], lw = 2)
  #text(0.2, 8e2, "HadCM3 yearly res.", col = "#074893")
  #text(0.3, 3e2, "5y filter", col = "#074893")
  #text(0.3, 1e2, "50y filter", col = "#074893")
  
  axis(side = 1, at = c(0.002,0.0033333333, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5), labels = FALSE)
  mtext(side=1,at=c(0.0033333333, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5),line = 1, 
        c(300, 1/0.005, 1/0.01, 1/0.02, 1/0.05, 1/0.2, 1/0.5),las=1,col="black", cex = cex_axis)
  axis(side = 4, at = c(1e-3, 1e-1, 1e1, 1e3), labels = FALSE)
  mtext(side=4,at=c(1e-3, 1e-1, 1e1, 1e3),line = 1, c(1e-3, 0.1, 10, 1000),las=1,col="black", cex = cex_axis)
  
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_ds_",var)]]$spec, col = adjustcolor(COLZ[4], alpha = 0.5), lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_ds_",var)]]$spec), col = COLZ[4], lw = 2)
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_down_rec_",var)]]$spec, col = COLZ[5], lty = 3, lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_down_rec_",var)]]$spec), col = COLZ[5], lw = 2)
  #text(0.02, 0.01, "HadCM3 down-sampled", col = COLZ[4])
  #text(0.02, 0.001, "10y filter", col = COLZ[4])
  
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC$Record$spec, col = adjustcolor(COLZ[6], alpha = 0.5), lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC$Record$spec), col = COLZ[6], lw = 3)
  #text(0.02, 0.0005, "Records", col = "black")
  
  legend("bottomleft", legend = c("HadCM3 (full)", paste0(" ...  ",filter[[var]][1],"y filter"), paste0(" ...  ",filter[[var]][2],"y filter"), "HadCM3 (down-sampled)", paste0(" ...  ",filter[[var]][3],"y filter"), "Records"), 
         col = c(COLZ[1],COLZ[2],COLZ[3],COLZ[4],COLZ[5],COLZ[6]), lwd = c(1,1,1,1,1,2), lty = c(1,1,1,1,1,1), bty = "n", cex = cex_text)
  
  if(var == "ISOT"){
    text(0.4, 500, TeX("Spectrum HadCM3: $\\delta^{18}O$"), adj = 1, cex = cex_text)
  }
  if(var == "ITPC"){
    text(0.4, 500, TeX("Spectrum HadCM3: $\\delta^{18}O_{pw}$"), adj = 1, cex = cex_text)
  }
  mtext("d)", side = 3, line = -1.3, adj = 0.025, cex = cex_axis_text)
  
  
  dev.off()
  
}

