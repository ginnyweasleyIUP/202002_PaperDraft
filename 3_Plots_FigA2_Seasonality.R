#################################################
## Paper Figure 10 ###############################
#################################################

## SEASONALITY

#################################################

source("Functions/Plotting/STACYmap_5.R")
source("Functions/projection_ptlyr.R")

for(run in c("a","b","c")){
  PLOT <- list()
  for(var in c("TEMP", "PREC","ISOT")){
    Plot_Lyr <- data.frame(lon = ANALYSIS$SEASONS[[paste0("Plot_Lyr_", var, "_", run)]]$lon, 
                           lat = ANALYSIS$SEASONS[[paste0("Plot_Lyr_", var, "_", run)]]$lat, 
                           layer = ANALYSIS$SEASONS[[paste0("Plot_Lyr_", var, "_", run)]]$season, 
                           value = ANALYSIS$SEASONS[[paste0("Plot_Lyr_", var, "_", run)]]$value)
    
    mask_china <- logical(length(Plot_Lyr$lon))
    
    for(ii in 1:length(Plot_Lyr$lon)){
      if(is.na(Plot_Lyr$lon[ii])){next}
      if(Plot_Lyr$lon[ii] > 100 & Plot_Lyr$lon[ii] < 120){
        if(Plot_Lyr$lat[ii] < 35 & Plot_Lyr$lat[ii] > 22){
          mask_china[ii] = T}
      }
    }
    
    ptlyr_china <- data.frame(
      lon = Plot_Lyr$lon[mask_china],
      lat = Plot_Lyr$lat[mask_china],
      layer = Plot_Lyr$layer[mask_china],
      value = Plot_Lyr$value[mask_china]
    )
    
    ptlyr_rest <- data.frame(
      lon = Plot_Lyr$lon[!mask_china],
      lat = Plot_Lyr$lat[!mask_china],
      layer = Plot_Lyr$layer[!mask_china],
      value = Plot_Lyr$value[!mask_china]
    )
    
    
    ptlyr_china_p <- projection_ptlyr(ptlyr_china, as.character('+proj=robin +datum=WGS84'))#
    ptlyr_rest_p <- projection_ptlyr(ptlyr_rest, as.character('+proj=robin +datum=WGS84'))
    ptlyr_china_p$value <- abs(ptlyr_china$value)
    ptlyr_rest_p$value <- abs(ptlyr_rest$value)
    
    ptlyr_china_p$layer <- factor(ptlyr_china$layer)
    ptlyr_rest_p$layer <- factor(ptlyr_rest$layer)
    
    if(var =="TEMP"){title = "Temperature"}
    if(var =="PREC"){title = "Precipitation"}
    if(var =="ISOT"){title = "Isotopic composition"}
    
    
    plot <- STACYmap(coastline = T) +
      geom_point(data = ptlyr_china_p, aes(x = long, y = lat, color = layer, size = value), shape = 19, 
                 show.legend = c(color =TRUE, size = TRUE), position = position_jitter(width = 1000000, height = 500000)) +
      geom_point(data = ptlyr_rest_p, aes(x = long, y = lat, color = layer, size = value), shape = 19, 
                 show.legend = c(color =TRUE, size = TRUE)) +
      scale_color_manual(name = "Seasons", labels = c("DJF", "MAM", "JJA", "SON", "year"), 
                         
                         values = c("#0072b2", "#009e73", "#d55e00", "#f0e442", "#483D3F")) +
      scale_size(name = "abs. corr.") +
      ggtitle(title) +
      theme(panel.border = element_blank(),
            legend.background = element_blank(),
            axis.text = element_blank(),
            legend.direction = "horizontal",
            text = element_text(size = 10),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            plot.title = element_text(hjust = 0.5))+
      guides(color=guide_legend(nrow=2,byrow=TRUE), size = guide_legend(nrow = 2, byrow = T))
    
    #plot
    # plot  %>% ggsave(filename = paste('Paper_Plot_A2_Seasons_prec', 'pdf', sep = '.'), plot = ., path = 'Plots', 
    #                  width = 2*8.3, height =2*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")
    PLOT[[var]] <- plot
  }
  
  library(ggpubr)
  plot <- ggarrange(PLOT$TEMP, PLOT$PREC, PLOT$ISOT,
                    labels = c("(a)", "(b)", "(c)"),
                    ncol = 1, nrow = 3)
  plot
  
  
  plot  %>% ggsave(filename = paste0('Paper_Plot_A2_Seasons_xnap',run, '.pdf'), plot = ., path = 'Plots', 
                   width = 2*PLOTTING_VARIABLES$WIDTH, height = 4*12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")
  plot  %>% ggsave(filename = paste0('Paper_Plot_A2_Seasons_xnap',run, '.png'), plot = ., path = 'Plots', 
                   width = 2*PLOTTING_VARIABLES$WIDTH, height = 4*12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "png")
}
