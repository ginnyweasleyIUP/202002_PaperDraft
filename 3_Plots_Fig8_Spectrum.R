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
  pdf(file = paste0("Plots/Paper_Plot_7_Spectra_",var,".pdf"), width = PLOTTING_VARIABLES$WIDTH, height = PLOTTING_VARIABLES$HEIGHT*1.2)
  #png(file = paste0("Plots/Paper_Plot_7_Spectra_",var,".png"), width = 70*PLOTTING_VARIABLES$WIDTH, height = 70*PLOTTING_VARIABLES$HEIGHT*1.2)
  LPlot(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec, col = adjustcolor("#074893", alpha = 0.5), 
        ylim = c(0.00001,1000), xlim = c(1/300, 0.5),
        xaxt = 'n',
        xlab = "")
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec), col = "#074893", lw = 2)#,
  #main = TeX("Mean Spectra from cave locations (res>8)"))
  mtext("Periode (years)", side = 1, line= 2)
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_down_",var)]]$spec, col = "#074893", lty = 3 )
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_rec_",var)]]$spec, col = "#074893", lty = 3)
  #text(0.2, 8e2, "HadCM3 yearly res.", col = "#074893")
  #text(0.3, 3e2, "5y filter", col = "#074893")
  #text(0.3, 1e2, "50y filter", col = "#074893")
  
  axis(side = 1, at = c(0.002,0.0033333333, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5), 
       labels = c(1/0.002, 300, 1/0.005, 1/0.01, 1/0.02, 1/0.05, 1/0.2, 1/0.5))
  
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_ds_",var)]]$spec, col = adjustcolor("#91002B", alpha = 0.5))
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_ds_",var)]]$spec), col = "#91002B")
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_down_rec_",var)]]$spec, col = "#91002B", lty = 3)
  #text(0.02, 0.01, "HadCM3 down-sampled", col = "#91002B")
  #text(0.02, 0.001, "10y filter", col = "#91002B")
  
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC$Record$spec, col = adjustcolor("black", alpha = 0.5), lw = 1)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC$Record$spec), col = "black", lw = 2)
  #text(0.02, 0.0005, "Records", col = "black")
  
  if(var == "ISOT"){
    legend("bottomleft", legend = c("HadCM3 yearly res.", " ...  2y filter", " ... 20y filter", "HadCM3 down-sampled to record res.", " ... 10y filter", "Records"), 
           col = c("#074893","#074893","#074893","#91002B","#91002B","black"), lwd = c(1,1,1,1,1,2), lty = c(1,3,3,1,3,1), bty = "n")  
  }
  if(var == "ITPC"){
    legend("bottomleft", legend = c("HadCM3 yearly res.", " ...  3y filter", " ... 12y filter", "HadCM3 down-sampled to record res.", " ... 4y filter", "Records"), 
           col = c("#074893","#074893","#074893","#91002B","#91002B","black"), lwd = c(1,1,1,1,1,2), lty = c(1,3,3,1,3,1), bty = "n")
  }
  
  
  dev.off()
  
}

