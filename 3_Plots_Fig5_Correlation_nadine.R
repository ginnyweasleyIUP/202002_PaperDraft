#################################################
## Paper Figure 1 ###############################
#################################################

## Here analysis and Plotting

## Einführungsplot (GENERAL --> von Kira kopieren)

library(dplyr)
library(latex2exp)
source("Functions/Plotting/STACYmap_5.R")
source("Functions/Plotting/STACYmap_5_1_NAgrid.R")

# this is correlation with downsampled temp and prec

#Ersetze einfach die Plot_lyr_temp mit der Martix inder für jede gridbox der Wert der Korrelation zw temp und d18O steht
# Plot_lyr_temp_p ist dann der p-value (damit nachher die unsignifikanten rausgefiltert werden können)
# Same for prec

# run a,b,c musst du wahrscheinlich raus machen (das ist nur, weil ich mit allen dreien arbeite)

#################################################
## PLOTS ########################################
#################################################
for(run in c("a","b","c")){
  print(run)
  Plot_lyr_temp <- ANALYSIS$CORR$FIELD[[run]]$CORR_TEMP_ISOT
  Plot_lyr_temp_p <- ANALYSIS$CORR$FIELD[[run]]$CORR_TEMP_ISOT_P
  Plot_lyr_temp[Plot_lyr_temp_p > 0.1] <- NA
  #Plot_lyr_temp[abs(Plot_lyr_temp) < 0.2] <- NA
  Plot_lyr_prec <- ANALYSIS$CORR$FIELD[[run]]$CORR_PREC_ISOT
  Plot_lyr_prec_p <- ANALYSIS$CORR$FIELD[[run]]$CORR_PREC_ISOT_P
  Plot_lyr_prec[Plot_lyr_prec_p > 0.1] <- NA
  #Plot_lyr_prec[abs(Plot_lyr_prec) < 0.2] <- NA
  
  
  #Das hier ist zum anpassen weil hadCM3 von 0 bis 360 und nicht von -180 bis 180 geht. 
  #Falls du das für dich schon woanders geändert hast, dann musst du diesen Schritt einfach auskommentieren 
  Plot_lyr_temp <- rbind(Plot_lyr_temp[49:96,1:73],
                         Plot_lyr_temp[1:48,1:73])
  Plot_lyr_prec <- rbind(Plot_lyr_prec[49:96,1:73],
                         Plot_lyr_prec[1:48,1:73])
  
  
  #NA Layer
  NA_plot_lyr <- Plot_lyr_temp
  NA_plot_lyr[!is.na(NA_plot_lyr)] <- 0
  NA_plot_lyr[is.na(NA_plot_lyr)] <- 1
  
  #Das ist meine extra NA plto Funktion die für alle NA leere Kästchen plottet
  plot_temp <- STACYmap_NA(gridlyr = Plot_lyr_temp, centercolor = 0, graticules = T,
                           NA_gridlyr = NA_plot_lyr, NA_color = "grey", legend_names = list(grid = TeX("$\\rho (T, \\delta^{18}O)$"))) +
    theme(panel.border = element_blank(),
          legend.background = element_blank(),
          axis.text = element_blank(),
          text = element_text(size = 12),
          legend.title = element_text(size = 12))
  
  plot_temp
  
  NA_plot_lyr <- Plot_lyr_prec
  NA_plot_lyr[!is.na(NA_plot_lyr)] <- 0
  NA_plot_lyr[is.na(NA_plot_lyr)] <- 1
  
  plot_prec <- STACYmap_NA(gridlyr = Plot_lyr_prec, centercolor = 0, graticules = T,
                           NA_gridlyr = NA_plot_lyr, NA_color = "grey",
                           legend_names = list(grid = TeX("$\\rho (P, \\delta^{18}O)$"))) + 
    theme(panel.border = element_blank(),
          legend.background = element_blank(),
          axis.text = element_blank(),
          text = element_text(size = 12),
          legend.title = element_text(size = 12))
  
  #plot_prec
  
  
  #zusammenbauen von zwei 
  library(ggpubr)
  plot <- ggarrange(plot_temp, plot_prec,
                    labels = c("(a)", "(b)"),
                    ncol = 2, nrow = 1)
  
  #Pfade ändern!!!
  plot  %>% ggsave(filename = paste0('Paper_Plot_5_Correlation_xnap',run, '.pdf'), plot = ., path = 'Plots', 
                   width = 2*12, height = 8, units = 'cm', dpi = 'print', device = "pdf")
  plot  %>% ggsave(filename = paste0('Paper_Plot_5_Correlation_xnap',run, '.png'), plot = ., path = 'Plots', 
                   width = 2*12, height = 8, units = 'cm', dpi = 'print', device = "png")
}

remove(COR, double_time, plot, Plot_lyr_prec, Plot_lyr_prec_p, Plot_lyr_temp, Plot_lyr_temp_p, plot_prec, plot_temp, Point_Lyr_prec, Point_Lyr_temp, s, test)
remove(allmax, allmax_real, diff_dt, entity, ii, length_cave, record, run, sim, site, var)
