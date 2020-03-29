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

for(var in c("ISOT", "ITPC")){
  pdf(file = paste0("Plots/Paper_Plot_7_Spectra_",var,".pdf"), width = PLOTTING_VARIABLES$WIDTH*2/3, height = PLOTTING_VARIABLES$HEIGHT*1.2*2/3)
  #png(file = paste0("Plots/Paper_Plot_7_Spectra_",var,".png"), width = 70*PLOTTING_VARIABLES$WIDTH*2/3, height = 70*PLOTTING_VARIABLES$HEIGHT*1.2*2/3)
  LPlot(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec, col = adjustcolor("#08306b", alpha = 0.5), 
        ylim = c(0.00005,1000), xlim = c(1/300, 0.5),
        ylab = "",
        xaxt = 'n',
        yaxt = "n",
        xlab = "", lwd = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec), col = "#08306b", lw = 3)#,
  #main = TeX("Mean Spectra from cave locations (res>8)"))
  mtext("Period (y)", side = 1, line= 2)
  mtext("Power spectral sensity", side = 2, line= 2)
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_down_",var)]]$spec, col = "#4292c6", lty = 3, lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_down_",var)]]$spec), col = "#4292c6", lw = 2)
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_rec_",var)]]$spec, col = "#c6dbef", lty = 3, lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_rec_",var)]]$spec), col = "#c6dbef", lw = 2)
  #text(0.2, 8e2, "HadCM3 yearly res.", col = "#074893")
  #text(0.3, 3e2, "5y filter", col = "#074893")
  #text(0.3, 1e2, "50y filter", col = "#074893")
  
  axis(side = 1, at = c(0.002,0.0033333333, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5), 
       labels = c(1/0.002, 300, 1/0.005, 1/0.01, 1/0.02, 1/0.05, 1/0.2, 1/0.5))
  axis(side = 2, at = c(1e-3, 1e-1, 1e1, 1e3), 
       labels = c(1e-3, 0.1, 10, 1000))
  
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_ds_",var)]]$spec, col = adjustcolor("#91002B", alpha = 0.5), lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_ds_",var)]]$spec), col = "#91002B", lw = 2)
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_down_rec_",var)]]$spec, col = "#91002B", lty = 3, lw = 2)
  #text(0.02, 0.01, "HadCM3 down-sampled", col = "#91002B")
  #text(0.02, 0.001, "10y filter", col = "#91002B")
  
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC$Record$spec, col = adjustcolor("black", alpha = 0.5), lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC$Record$spec), col = "black", lw = 3)
  #text(0.02, 0.0005, "Records", col = "black")
  
  if(var == "ISOT"){
    legend("bottomleft", legend = c("HadCM3 (full)", " ...  2y filter", " ... 20y filter", "HadCM3 (down-sampled)", " ... 10y filter", "Records"), 
           col = c("#08306b","#4292c6","#c6dbef","#91002B","#91002B","black"), lwd = c(1,1,1,1,1,2), lty = c(1,1,1,1,3,1), bty = "n") 
    text(0.5, 500, TeX("HadCM3: $\\delta^{18}O$ in prec"), adj = 1)
  }
  if(var == "ITPC"){
    legend("bottomleft", legend = c("HadCM3 (full)", " ...  3y filter", " ... 12y filter", "HadCM3 (down-sampled)", " ... 4y filter", "Records"), 
           col = c("#08306b","#4292c6","#c6dbef","#91002B","#91002B","black"), lwd = c(1,1,1,1,1,2), lty = c(1,1,1,1,3,1), bty = "n") 
    text(0.5, 500, TeX("HadCM3: prec-weighted $\\delta^{18}O$"), adj = 1)
  }
  
  
  dev.off()
  
}

for(var in c("ISOT", "ITPC")){
  pdf(file = paste0("Plots/Appendix/Paper_Plot_7_Spectra_Compare",var,".pdf"), width = PLOTTING_VARIABLES$WIDTH*2/3, height = PLOTTING_VARIABLES$HEIGHT*1.2*2/3)
  #par(mfrow=c(1,2))
  cex_text = 1
  cex_axis = 0.6
  cex_axis_text = 0.7
  par(mar=c(1,0,0,1),oma=c(5,5,2,3),xaxs="i",yaxs="i",cex=1,lwd=2)
  layout(matrix(c(1,2,3,4,4,4), 3, 2, byrow = F))
  plot(seq(from = 1950 - DATA_past1000$time[[1]], to = 1950 - DATA_past1000$time[[2]]+1, by = -1), DATA_past1000$CAVES$sim_data_yearly$CAVE117[[paste0(var,"_a")]], 
       col = "#08306b", type = "l", ylim = c(-10,-5.5), xlim = c(0,1100), xaxt = "n", yaxt = "n",ylab = "", xlab = "")
  axis(2,at=seq(-10,-6,by=1),labels=FALSE,col="#08306b")
  mtext(side=2,at=seq(-10,-6,by=1),line = 1, seq(-10,-6,by=1),col="#08306b", cex = cex_axis)
  text(1000, -6, TeX("Bunker eID 240 full"), col = "#08306b", adj = 1, cex = cex_text)
  
  plot(DATA_past1000$CAVES$sim_data_downsampled$ENTITY240$interp_age, DATA_past1000$CAVES$sim_data_downsampled$ENTITY240[[paste0(var,"_a")]],
        col = "#91002B", type = "l", ylim = c(-10,-5.5), xlim = c(0,1100), xaxt = "n", yaxt = "n",ylab = "", xlab = "")
  axis(2,at=seq(-10,-6,by=1),labels=FALSE,col="#91002B")
  mtext(side=2,at=seq(-10,-6,by=1),line = 1, seq(-10,-6,by=1),las=1,col="#91002B", cex = cex_axis)
  text(1000, -6, TeX("down-sampled"), col = "#91002B", adj = 1, cex = cex_text)
  mtext(TeX("$\\delta^{18}O$"),side=2, at = -8, line = 2.5, las = 1, col = "#91002B", cex = cex_axis_text)
  
  plot(DATA_past1000$CAVES$record_data$ENTITY240$interp_age, DATA_past1000$CAVES$record_data$ENTITY240$d18O_dw_eq_a,
       col = "black", type = "l", ylim = c(-10,-5.5), xlim = c(0,1100), xaxt = "n", yaxt = "n",ylab = "", xlab = "")
  axis(2,at=seq(-10,-6,by=1),labels=FALSE,col="black")
  mtext(side=2,at=seq(-10,-6,by=1),line = 1, seq(-10,-6,by=1),las=1,col="black", cex = cex_axis)
  text(1000, -6, TeX("record"), col = "black", adj = 1, cex = cex_text)
  
  axis(1,at=seq(0,1000,by=200),labels=FALSE,col="black")
  mtext(side=1,at=seq(0,1000,by=200),line = 1, seq(0,1000,by=200),las=1,col="black", cex = cex_axis)
  mtext("years BP", side = 1, line = 3, at = 500, cex = cex_axis_text)
  
  #png(file = paste0("Plots/Paper_Plot_7_Spectra_",var,".png"), width = 70*PLOTTING_VARIABLES$WIDTH*2/3, height = 70*PLOTTING_VARIABLES$HEIGHT*1.2*2/3)
  LPlot(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec, col = adjustcolor("#08306b", alpha = 0.5), 
        ylim = c(0.00005,1000), xlim = c(1/300, 0.5),
        ylab = "",
        xaxt = 'n',
        yaxt = "n",
        xlab = "", lwd = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec), col = "#08306b", lw = 3)#,
  #main = TeX("Mean Spectra from cave locations (res>8)"))
  mtext("Period (y)", side = 1, line= 3, cex = cex_axis_text)
  mtext("Power spectral sensity", side = 4, line= 2.5, cex = cex_axis_text)
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_down_",var)]]$spec, col = "#4292c6", lty = 3, lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_down_",var)]]$spec), col = "#4292c6", lw = 2)
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_rec_",var)]]$spec, col = "#c6dbef", lty = 3, lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_rec_",var)]]$spec), col = "#c6dbef", lw = 2)
  #text(0.2, 8e2, "HadCM3 yearly res.", col = "#074893")
  #text(0.3, 3e2, "5y filter", col = "#074893")
  #text(0.3, 1e2, "50y filter", col = "#074893")
  
  axis(side = 1, at = c(0.002,0.0033333333, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5), labels = FALSE)
  mtext(side=1,at=c(0.0033333333, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5),line = 1, 
        c(300, 1/0.005, 1/0.01, 1/0.02, 1/0.05, 1/0.2, 1/0.5),las=1,col="black", cex = cex_axis)
  axis(side = 4, at = c(1e-3, 1e-1, 1e1, 1e3), labels = FALSE)
  mtext(side=4,at=c(1e-3, 1e-1, 1e1, 1e3),line = 1, c(1e-3, 0.1, 10, 1000),las=1,col="black", cex = cex_axis)
  
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_ds_",var)]]$spec, col = adjustcolor("#91002B", alpha = 0.5), lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_ds_",var)]]$spec), col = "#91002B", lw = 2)
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_down_rec_",var)]]$spec, col = "#91002B", lty = 3, lw = 2)
  #text(0.02, 0.01, "HadCM3 down-sampled", col = "#91002B")
  #text(0.02, 0.001, "10y filter", col = "#91002B")
  
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC$Record$spec, col = adjustcolor("black", alpha = 0.5), lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC$Record$spec), col = "black", lw = 3)
  #text(0.02, 0.0005, "Records", col = "black")
  
  if(var == "ISOT"){
    legend("bottomleft", legend = c("HadCM3 (full)", " ...  2y filter", " ... 20y filter", "HadCM3 (down-sampled)", " ... 10y filter", "Records"), 
           col = c("#08306b","#4292c6","#c6dbef","#91002B","#91002B","black"), lwd = c(1,1,1,1,1,2), lty = c(1,1,1,1,3,1), bty = "n", cex = cex_text) 
    text(0.2, 500, TeX("HadCM3: $\\delta^{18}O$ in prec"), adj = 1, cex = cex_text)
  }
  if(var == "ITPC"){
    legend("bottomleft", legend = c("HadCM3 (full)", " ...  3y filter", " ... 12y filter", "HadCM3 (down-sampled)", " ... 4y filter", "Records"), 
           col = c("#08306b","#4292c6","#c6dbef","#91002B","#91002B","black"), lwd = c(1,1,1,1,1,2), lty = c(1,1,1,1,3,1), bty = "n", cex = cex_text) 
    text(0.5, 500, TeX("Spectrum HadCM3: prec-weighted $\\delta^{18}O$"), adj = 1, cex = cex_text)
  }
  
  
  dev.off()
  
}

