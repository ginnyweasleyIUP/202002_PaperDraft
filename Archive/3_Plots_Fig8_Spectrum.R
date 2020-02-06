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

pdf(file = "Plots/Paper_Plot_7_Spectra.pdf", width = PLOTTING_VARIABLES$WIDTH, height = PLOTTING_VARIABLES$HEIGHT*1.2)
LPlot(ANALYSIS$SPECTRA$MEAN_SPEC$SIM_full_isot$spec, col = "#074893", 
      ylim = c(0.00001,1000), xlim = c(1/500, 0.5),
      xaxt = 'n',
      xlab = "")#,
#main = TeX("Mean Spectra from cave locations (res>8)"))
mtext("Periode (years)", side = 1, line= 2)
LLines(ANALYSIS$SPECTRA$MEAN_SPEC$SIM_full_down_isot$spec, col = "#074893", lty = 3 )
LLines(ANALYSIS$SPECTRA$MEAN_SPEC$SIM_full_rec_isot$spec, col = "#074893", lty = 3)
#text(0.2, 8e2, "HadCM3 yearly res.", col = "#074893")
#text(0.3, 3e2, "5y filter", col = "#074893")
#text(0.3, 1e2, "50y filter", col = "#074893")

axis(side = 1, at = c(0.002, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5), 
     labels = c(1/0.001, 1/0.005, 1/0.01, 1/0.02, 1/0.05, 1/0.2, 1/0.5))

LLines(ANALYSIS$SPECTRA$MEAN_SPEC$SIM_ds_isot$spec, col = "#91002B")
LLines(ANALYSIS$SPECTRA$MEAN_SPEC$SIM_down_rec_isot$spec, col = "#91002B", lty = 3)
#text(0.02, 0.01, "HadCM3 down-sampled", col = "#91002B")
#text(0.02, 0.001, "10y filter", col = "#91002B")

LLines(ANALYSIS$SPECTRA$MEAN_SPEC$Record$spec, col = "black", lw = 2)
#text(0.02, 0.0005, "Records", col = "black")



legend("bottomleft", legend = c("HadCM3 yearly res.", " ... 2y filter", " ... 20y filter", "HadCM3 down-sampled to record res.", " ... 10y filter", "Records"), 
       col = c("#074893","#074893","#074893","#91002B","#91002B","black"), lwd = c(1,1,1,1,1,2), lty = c(1,3,3,1,3,1), bty = "n")
dev.off()


## APPENDIX #####################################

pdf(file = "Plots/Appendix/Paper_Plot_7_Spectra_isot_weighted.pdf", width = PLOTTING_VARIABLES$WIDTH, height = PLOTTING_VARIABLES$HEIGHT*1.2)
LPlot(ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH$SIM_full_isot$spec, col = "#074893", 
      ylim = c(0.00001,1000), xlim = c(1/500, 0.5),
      xaxt = 'n',
      xlab = "")#,
#main = TeX("Mean Spectra from cave locations (res>8)"))
mtext("Periode (years)", side = 1, line= 2)
LLines(ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH$SIM_full_down_isot$spec, col = "#074893", lty = 3 )
LLines(ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH$SIM_full_rec_isot$spec, col = "#074893", lty = 3)
#text(0.2, 8e2, "HadCM3 yearly res.", col = "#074893")
#text(0.3, 3e2, "5y filter", col = "#074893")
#text(0.3, 1e2, "50y filter", col = "#074893")

axis(side = 1, at = c(0.002, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5), 
     labels = c(1/0.001, 1/0.005, 1/0.01, 1/0.02, 1/0.05, 1/0.2, 1/0.5))

LLines(ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH$SIM_ds_isot$spec, col = "#91002B")
LLines(ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH$SIM_down_rec_isot$spec, col = "#91002B", lty = 3)
#text(0.02, 0.01, "HadCM3 down-sampled", col = "#91002B")
#text(0.02, 0.001, "10y filter", col = "#91002B")

LLines(ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH$Record$spec, col = "black", lw = 2)
#text(0.02, 0.0005, "Records", col = "black")



legend("bottomleft", legend = c("HadCM3 yearly res.", " ... 2y filter", " ... 20y filter", "HadCM3 down-sampled to record res.", " ... 10y filter", "Records"), 
       col = c("#074893","#074893","#074893","#91002B","#91002B","black"), lwd = c(1,1,1,1,1,2), lty = c(1,3,3,1,3,1), bty = "n")
dev.off()


pdf(file = "Plots/Appendix/Paper_Plot_7_Spectra_temp.pdf", width = PLOTTING_VARIABLES$WIDTH, height = PLOTTING_VARIABLES$HEIGHT*1.2)
LPlot(ANALYSIS$SPECTRA$MEAN_SPEC$SIM_full_temp$spec, col = "#074893", 
      ylim = c(0.00001,1000), 
      xlim = c(1/500, 0.5),
      xaxt = 'n',
      xlab = "")#,
#main = TeX("Mean Spectra from cave locations (res>8)"))
mtext("Periode (years)", side = 1, line= 2)
LLines(ANALYSIS$SPECTRA$MEAN_SPEC$SIM_full_down_temp$spec, col = "#074893", lty = 3 )
LLines(ANALYSIS$SPECTRA$MEAN_SPEC$SIM_full_rec_temp$spec, col = "#074893", lty = 3)
#text(0.2, 8e2, "HadCM3 yearly res.", col = "#074893")
#text(0.3, 3e2, "5y filter", col = "#074893")
#text(0.3, 1e2, "50y filter", col = "#074893")

axis(side = 1, at = c(0.002, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5), 
     labels = c(1/0.001, 1/0.005, 1/0.01, 1/0.02, 1/0.05, 1/0.2, 1/0.5))

LLines(ANALYSIS$SPECTRA$MEAN_SPEC$SIM_ds_temp$spec, col = "#91002B")
LLines(ANALYSIS$SPECTRA$MEAN_SPEC$SIM_down_rec_temp$spec, col = "#91002B", lty = 3)
#text(0.02, 0.01, "HadCM3 down-sampled", col = "#91002B")
#text(0.02, 0.001, "10y filter", col = "#91002B")

LLines(ANALYSIS$SPECTRA$MEAN_SPEC$Record$spec, col = "black", lw = 2)
#text(0.02, 0.0005, "Records", col = "black")



legend("bottomleft", legend = c("HadCM3 yearly res.", " ... 2y filter", " ... 20y filter", "HadCM3 down-sampled to record res.", " ... 10y filter", "Records"), 
       col = c("#074893","#074893","#074893","#91002B","#91002B","black"), lwd = c(1,1,1,1,1,2), lty = c(1,3,3,1,3,1), bty = "n")
dev.off()
