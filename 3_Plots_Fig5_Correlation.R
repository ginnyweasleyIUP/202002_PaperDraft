#################################################
## Paper Figure 1 ###############################
#################################################

## Here analysis and Plotting

## Einführungsplot (GENERAL --> von Kira kopieren)

library(dplyr)
library(latex2exp)
source("Functions/STACYmap_5.R")


#################################################
## PLOTS ########################################
#################################################
for(run in c("a","b","c")){
  Plot_lyr_temp <- ANALYSIS$CORR$FIELD[[run]]$CORR_TEMP_ISOT
  Plot_lyr_temp_p <- ANALYSIS$CORR$FIELD[[run]]$CORR_TEMP_ISOT_P
  Plot_lyr_temp[Plot_lyr_temp_p > 0.1] <- NA
  Plot_lyr_temp[abs(Plot_lyr_temp) < 0.2] <- NA
  Plot_lyr_prec <- ANALYSIS$CORR$FIELD[[run]]$CORR_PREC_ISOT
  Plot_lyr_prec_p <- ANALYSIS$CORR$FIELD[[run]]$CORR_PREC_ISOT_P
  Plot_lyr_prec[Plot_lyr_prec_p > 0.1] <- NA
  Plot_lyr_prec[abs(Plot_lyr_prec) < 0.2] <- NA
  
  Plot_lyr_temp <- rbind(Plot_lyr_temp[49:96,1:73],
                         Plot_lyr_temp[1:48,1:73])
  Plot_lyr_prec <- rbind(Plot_lyr_prec[49:96,1:73],
                         Plot_lyr_prec[1:48,1:73])
  
  ##### Point Layer
  
  Point_Lyr_temp <- list(lon = list(), lat = list(), value = list())
  Point_Lyr_prec <- list(lon = list(), lat = list(), value = list())
  
  length_cave = length(DATA_past1000$CAVES$entity_info$site_id)
  
  for(ii in 1:length_cave){
    site <- DATA_past1000$CAVES$entity_info$site_id[ii]
    print(ii)
    if(is.na(ANALYSIS$CORR$POINTS[[run]]$p_TEMP[ii])){next}
    # 1) sortiert aus, was nicht signifikant ist
    if(ANALYSIS$CORR$POINTS[[run]]$p_TEMP[ii] > 0.1){
      Point_Lyr_temp$lon = c(Point_Lyr_temp$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
      Point_Lyr_temp$lat = c(Point_Lyr_temp$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
      Point_Lyr_temp$value = c(Point_Lyr_temp$value, ANALYSIS$CORR$POINTS[[run]]$CORR_TEMP[ii])
      # 2) betrachte signifikante Korrelationen:
    }
    if(is.na(ANALYSIS$CORR$POINTS[[run]]$p_PREC[ii])){next}
    # 1) sortiert aus, was nicht signifikant ist
    if(ANALYSIS$CORR$POINTS[[run]]$p_PREC[ii] > 0.1){
      Point_Lyr_prec$lon = c(Point_Lyr_prec$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
      Point_Lyr_prec$lat = c(Point_Lyr_prec$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
      Point_Lyr_prec$value = c(Point_Lyr_prec$value, ANALYSIS$CORR$POINTS[[run]]$CORR_PREC[ii])
      # 2) betrachte signifikante Korrelationen:
    }
  }
  
  
  
  Point_Lyr_temp$lon = as.numeric(Point_Lyr_temp$lon)
  Point_Lyr_temp$lat = as.numeric(Point_Lyr_temp$lat)
  Point_Lyr_temp$value = as.numeric(Point_Lyr_temp$value)
  
  Point_Lyr_prec$lon = as.numeric(Point_Lyr_prec$lon)
  Point_Lyr_prec$lat = as.numeric(Point_Lyr_prec$lat)
  Point_Lyr_prec$value = as.numeric(Point_Lyr_prec$value)
  
  
  GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE <- 3
  
  plot_temp <- STACYmap(gridlyr = Plot_lyr_temp, centercolor = 0, graticules = T,
                        ptlyr = as.data.frame(Point_Lyr_temp), legend_names = list(grid = 'Temp.-Correlation (p<0.1)')) + 
    theme(panel.border = element_blank(),
          legend.background = element_blank(),
          axis.text = element_blank(),
          text = element_text(size = 12),
          legend.title = element_text(size = 12))
  
  plot_prec <- STACYmap(gridlyr = Plot_lyr_prec, centercolor = 0, graticules = T,
                        ptlyr = as.data.frame(Point_Lyr_prec), legend_names = list(grid = 'Prec.-Correlation (p<0.1)')) + 
    theme(panel.border = element_blank(),
          legend.background = element_blank(),
          axis.text = element_blank(),
          text = element_text(size = 12),
          legend.title = element_text(size = 12))
  
  library(ggpubr)
  plot <- ggarrange(plot_temp, plot_prec,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)
  
  plot  %>% ggsave(filename = paste('Paper_Plot_5_Correlation', 'pdf', sep = '.'), plot = ., path = 'Plots', 
                   width = 2*12, height = 12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")
}

