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
  pdf(file = paste0("Plots/Discussion/Spectra_1-48_",var,".pdf"), width = 21, height = 29.7)
  par(mfrow=c(8,6),oma = c(1,3,0,0) + 0.1,mar = c(3,0,1,3) + 0.1)
  for(entity in ANALYSIS$SPECTRA$entities_spec_rec[1:48]){
    LPlot(ANALYSIS$SPECTRA$SIM_full$a[[var]][[paste0("ENTITY",entity)]], col = adjustcolor(COLZ[1], alpha = 0.5), 
          ylim = c(0.00005,1000), xlim = c(1/300, 0.5),
          ylab = "",
          xaxt = 'n',
          yaxt = "n",
          xlab = "", lwd = 2)
    LLines(LogSmooth(ANALYSIS$SPECTRA$SIM_full$a[[var]][[paste0("ENTITY",entity)]]), col = COLZ[1], lw = 3)#,
    #main = TeX("Mean Spectra from cave locations (res>8)"))
    mtext("Period (y)", side = 1, line= 2, cex = 0.8)
    mtext("Power spectral sensity", side = 2, line= 2, cex = 0.8)
    LLines(ANALYSIS$SPECTRA$SIM_filter_full_down[[var]]$a[[paste0("ENTITY",entity)]], col = COLZ[2], lty = 3, lw = 2)
    LLines(LogSmooth(ANALYSIS$SPECTRA$SIM_filter_full_down[[var]]$a[[paste0("ENTITY",entity)]]), col = COLZ[2], lw = 2)
    LLines(ANALYSIS$SPECTRA$SIM_filter_full_rec[[var]]$a[[paste0("ENTITY",entity)]], col = COLZ[3], lty = 3, lw = 2)
    LLines(LogSmooth(ANALYSIS$SPECTRA$SIM_filter_full_rec[[var]]$a[[paste0("ENTITY",entity)]]), col = COLZ[3], lw = 2)
    #text(0.2, 8e2, "HadCM3 yearly res.", col = "#074893")
    #text(0.3, 3e2, "5y filter", col = "#074893")
    #text(0.3, 1e2, "50y filter", col = "#074893")
    
    axis(side = 1, at = c(0.002, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5), 
         labels = c(1/0.002, 1/0.005, 1/0.01, 1/0.02, 1/0.05, 1/0.2, 1/0.5))
    axis(side = 2, at = c(1e-3, 1e-1, 1e1, 1e3), 
         labels = c(1e-3, 0.1, 10, 1000))
    
    LLines(ANALYSIS$SPECTRA$SIM_ds$a[[var]][[paste0("ENTITY",entity)]], col = adjustcolor(COLZ[4], alpha = 0.5), lw = 2)
    LLines(LogSmooth(ANALYSIS$SPECTRA$SIM_ds$a[[var]][[paste0("ENTITY",entity)]]), col = COLZ[4], lw = 2)
    LLines(ANALYSIS$SPECTRA$SIM_filter_down_rec[[var]]$a[[paste0("ENTITY",entity)]], col = COLZ[5], lty = 3, lw = 2)
    LLines(LogSmooth(ANALYSIS$SPECTRA$SIM_filter_down_rec[[var]]$a[[paste0("ENTITY",entity)]]), col = COLZ[5], lw = 2)
    #text(0.02, 0.01, "HadCM3 down-sampled", col = COLZ[4])
    #text(0.02, 0.001, "10y filter", col = COLZ[4])
    
    LLines(ANALYSIS$SPECTRA$RECORDS[[paste0("ENTITY",entity)]], col = adjustcolor(COLZ[6], alpha = 0.5), lw = 2)
    LLines(LogSmooth(ANALYSIS$SPECTRA$RECORDS[[paste0("ENTITY",entity)]]), col = COLZ[6], lw = 3)
    #text(0.02, 0.0005, "Records", col = "black")
    legend("bottomleft", legend = c("HadCM3 (full)", paste0(" ...  ",filter[[var]][1],"y filter"), paste0(" ...  ",filter[[var]][2],"y filter"), "HadCM3 (down-sampled)", paste0(" ...  ",filter[[var]][3],"y filter"), "Record"), 
           col = c(COLZ[1],COLZ[2],COLZ[3],COLZ[4],COLZ[5],COLZ[6]), lwd = c(1,1,1,1,1,2), lty = c(1,1,1,1,1,1), bty = "n", cex = cex_text)
    
    text(0.5,200, paste0("ENTITY ", entity), adj = 1)
    if(var == "ISOT"){
      text(0.5, 500, TeX("HadCM3: $\\delta^{18}O$ in prec"), adj = 1)
    }
    if(var == "ITPC"){
      text(0.5, 500, TeX("HadCM3: prec-weighted $\\delta^{18}O$"), adj = 1)
    }
  }
  
  
  
  dev.off()
  
  pdf(file = paste0("Plots/Discussion/Spectra_49-87_",var,".pdf"), width = 21, height = 29.7)
  par(mfrow=c(8,6),oma = c(1,3,0,0) + 0.1,mar = c(3,0,1,3) + 0.1)
  for(entity in ANALYSIS$SPECTRA$entities_spec_rec[49:87]){
    if(entity ==351 | entity == 390 ){next}
    LPlot(ANALYSIS$SPECTRA$SIM_full$a[[var]][[paste0("ENTITY",entity)]], col = adjustcolor(COLZ[1], alpha = 0.5), 
          ylim = c(0.00005,1000), xlim = c(1/300, 0.5),
          ylab = "",
          xaxt = 'n',
          yaxt = "n",
          xlab = "", lwd = 2)
    LLines(LogSmooth(ANALYSIS$SPECTRA$SIM_full$a[[var]][[paste0("ENTITY",entity)]]), col = COLZ[1], lw = 3)#,
    #main = TeX("Mean Spectra from cave locations (res>8)"))
    mtext("Period (y)", side = 1, line= 2, cex = 0.8)
    mtext("Power spectral sensity", side = 2, line= 2, cex = 0.8)
    LLines(ANALYSIS$SPECTRA$SIM_filter_full_down[[var]]$a[[paste0("ENTITY",entity)]], col = COLZ[2], lty = 3, lw = 2)
    LLines(LogSmooth(ANALYSIS$SPECTRA$SIM_filter_full_down[[var]]$a[[paste0("ENTITY",entity)]]), col = COLZ[2], lw = 2)
    LLines(ANALYSIS$SPECTRA$SIM_filter_full_rec[[var]]$a[[paste0("ENTITY",entity)]], col = COLZ[3], lty = 3, lw = 2)
    LLines(LogSmooth(ANALYSIS$SPECTRA$SIM_filter_full_rec[[var]]$a[[paste0("ENTITY",entity)]]), col = COLZ[3], lw = 2)
    #text(0.2, 8e2, "HadCM3 yearly res.", col = "#074893")
    #text(0.3, 3e2, "5y filter", col = "#074893")
    #text(0.3, 1e2, "50y filter", col = "#074893")
    
    axis(side = 1, at = c(0.002, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5), 
         labels = c(1/0.002, 1/0.005, 1/0.01, 1/0.02, 1/0.05, 1/0.2, 1/0.5))
    axis(side = 2, at = c(1e-3, 1e-1, 1e1, 1e3), 
         labels = c(1e-3, 0.1, 10, 1000))
    
    LLines(ANALYSIS$SPECTRA$SIM_ds$a[[var]][[paste0("ENTITY",entity)]], col = adjustcolor(COLZ[4], alpha = 0.5), lw = 2)
    LLines(LogSmooth(ANALYSIS$SPECTRA$SIM_ds$a[[var]][[paste0("ENTITY",entity)]]), col = COLZ[4], lw = 2)
    LLines(ANALYSIS$SPECTRA$SIM_filter_down_rec[[var]]$a[[paste0("ENTITY",entity)]], col = COLZ[5], lty = 3, lw = 2)
    LLines(LogSmooth(ANALYSIS$SPECTRA$SIM_filter_down_rec[[var]]$a[[paste0("ENTITY",entity)]]), col = COLZ[5], lw = 2)
    #text(0.02, 0.01, "HadCM3 down-sampled", col = COLZ[4])
    #text(0.02, 0.001, "10y filter", col = COLZ[4])
    
    LLines(ANALYSIS$SPECTRA$RECORDS[[paste0("ENTITY",entity)]], col = adjustcolor(COLZ[6], alpha = 0.5), lw = 2)
    LLines(LogSmooth(ANALYSIS$SPECTRA$RECORDS[[paste0("ENTITY",entity)]]), col = COLZ[6], lw = 3)
    #text(0.02, 0.0005, "Records", col = "black")
    legend("bottomleft", legend = c("HadCM3 (full)", paste0(" ...  ",filter[[var]][1],"y filter"), paste0(" ...  ",filter[[var]][2],"y filter"), "HadCM3 (down-sampled)", paste0(" ...  ",filter[[var]][3],"y filter"), "Record"), 
           col = c(COLZ[1],COLZ[2],COLZ[3],COLZ[4],COLZ[5],COLZ[6]), lwd = c(1,1,1,1,1,2), lty = c(1,1,1,1,1,1), bty = "n", cex = cex_text)
    
    text(0.5,200, paste0("ENTITY ", entity), adj = 1)
    if(var == "ISOT"){
      text(0.5, 500, TeX("HadCM3: $\\delta^{18}O$ in prec"), adj = 1)
    }
    if(var == "ITPC"){
      text(0.5, 500, TeX("HadCM3: prec-weighted $\\delta^{18}O$"), adj = 1)
    }
  }
  
  
  
  dev.off()
  
}


