#################################################
## Paper Figure 1 ###############################
#################################################

## Here analysis and Plotting

## Einführungsplot (GENERAL --> von Kira kopieren)

# ToDo:
# [X] auf downsampled umsteigen
# [X] colorbar breiter machen
# [ ] Schriftgrößen anpassen
# [ ] ABCDE hinzufügen

#PLOTTING_DATA <- list()

#################################################

PLOTTING_DATA$FIG4 <- list()
PLOTTING_DATA$FIG4$POINTS$Cavelyr <- ANALYSIS$VARIANCE$POINTS$CAVElyr
PLOTTING_DATA$FIG4$VAR_RATIOS <- list()
for(run in c("a","b","c")){
  PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_temp <- ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_temp[mask_var]
  PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_temp_ds <- ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_temp_ds[mask_var]
  PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_prec <- ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_prec[mask_var]
  PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_prec_ds <- ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_prec_ds[mask_var]
  PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_isot <- ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_isot[mask_var]
  PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_isot_ds <- ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_isot_ds[mask_var]
  PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_itpc <- ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_itpc[mask_var]
  PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_itpc_ds <- ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_itpc_ds[mask_var]
  
}



## PLOTTING #####################################

library(latex2exp)

## Total of 5 Plots that then have to be arranged
## 1) Variance Map
## 2) ISOT density
## 3) ITPC density
## 4) TEMP density
## 5) PREC density

## Variance Map #################################

source("Functions/Plotting/var_map_plot.R")

# mask_var_2 <- mask_var
# mask_var_2[85] <- FALSE
# mask_var_2[63] <- FALSE
# mask_var_2[64] <- FALSE
# mask_var_2[104] <- FALSE
# mask_var_2[38] <- FALSE

for(run in c("a","b","c")){
  Point_Lyr <- data.frame(
    lon = ANALYSIS$VARIANCE$POINTS$CAVElyr$lon[mask_var],
    lat = ANALYSIS$VARIANCE$POINTS$CAVElyr$lat[mask_var],
    value = log10(ANALYSIS$VARIANCE$POINTS$CAVElyr[[run]]$value_VR_itpc_ds[mask_var])
    #value = log(VAR_ANALYSIS$POINTS$CAVElyr$value_VR_isot_ds[mask_var])
  )
  
  Point_Lyr$lon[50] <- NA
  Point_Lyr$lat[50] <- NA
  Point_Lyr$value[50] <- NA
  
  GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE = 2.5
  
  # Points in China werden gejittered und sind nicht an ihrer richtigen Position!
  plot_var <- var_map_plot(Point_Lyr =  Point_Lyr, pt_size =  3, txt_size =  10)
  plot_var
  
  plot_var %>% ggsave(filename = paste0('Paper_Plot_4_Variance_1_map_xnap',run, '.pdf'), plot = ., path = 'Plots', 
                      width = 2*6, height = 2*PLOTTING_VARIABLES$HEIGHT/1.5, units = 'cm', dpi = 'print', device = "pdf")
  plot_var %>% ggsave(filename = paste0('Paper_Plot_4_Variance_1_map_xnap',run, '.png'), plot = ., path = 'Plots', 
                      width = 2*6, height = 2*PLOTTING_VARIABLES$HEIGHT/1.5, units = 'cm', dpi = 'print', device = "png")
}



#lines(c(median(log10(PLOTTING_DATA$FIG4$VAR_RATIOS$VR_sim_temp)),median(log10(PLOTTING_DATA$FIG4$VAR_RATIOS$VR_sim_temp))), c(0, ))
for(run in c("a","b","c")){
  #pdf(file = paste0("Plots/Paper_Plot_4_Variance_2_histo_xnap",run,".pdf"), width = 2*6, height = 2*PLOTTING_VARIABLES$HEIGHT/1.5)
  png(file = paste0("Plots/Paper_Plot_4_Variance_2_histo_xnap",run,".png"), width = 50*2*6, height = 50*2*PLOTTING_VARIABLES$HEIGHT/1.5)
  par(mfrow=c(2,2),oma = c(1,3,0,0) + 0.1,mar = c(3,1,0,1) + 0.1, new = FALSE)
  
  hist(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_isot), 
       breaks = 9, border = "white", prob = TRUE, 
       ylim = c(0,1), xlim = c(-2.5, 2.5), xlab = "",xaxt = 'n',
       main = "", cex.main = 1.5, cex.axis = 1.5)
  axis(side = 1, at = c(log10(0.01), log10(0.1), 0, log10(10), log10(100)), 
       labels = c(0.01, 0.1, 1, 10, 100), cex.axis = 1.5)
  lines(c(median(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_isot)),median(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_isot))), 
        c(0, max(density(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_isot), na.rm = T)$y)-0.01),
        lwd = 2, col = "black", lty = 2)
  lines(density(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_isot), na.rm = T),
        lwd = 2, col = "black")
  lines(c(median(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_isot_ds)),median(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_isot_ds))), 
        c(0, max(density(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_isot_ds), na.rm = T)$y)-0.01),
        lwd = 2, col = "#B2182B", lty = 2)
  lines(density(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_isot_ds), na.rm = T),
        lwd = 2, col = "#B2182B")
  abline(v=0, col = "grey60", lty = 3)
  #mtext(text = TeX("$$Var_{Rec}/Var_{Sim}$$"),side = 1,line = 2.5, cex = 1.5)
  mtext(text = "density",side = 2,line = 2.5, cex = 1.5)
  text(1.0, 0.6, "down-sampled", col = "#B2182B", cex = 1.5)
  text(-1.8, 0.6, "full", col = "black", cex = 1.5)
  text(0, 0.95, "d18O in precipitation", col = "black", cex = 1.5)
  
  
  hist(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_itpc), 
       breaks = 9, border = "white", prob = TRUE, 
       ylim = c(0,1), xlim = c(-2.5, 2.5), xlab = "",xaxt = 'n',
       main = "", cex.main = 1.5, cex.axis = 1.5)
  axis(side = 1, at = c(log10(0.01), log10(0.1), 0, log10(10), log10(100)), 
       labels = c(0.01, 0.1, 1, 10, 100), cex.axis = 1.5)
  lines(c(median(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_itpc)),median(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_itpc))), 
        c(0, max(density(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_itpc), na.rm = T)$y)-0.01),
        lwd = 2, col = "black", lty = 2)
  lines(density(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_itpc), na.rm = T),
        lwd = 2, col = "black")
  lines(c(median(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_itpc_ds)),median(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_itpc_ds))), 
        c(0, max(density(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_itpc_ds), na.rm = T)$y)-0.05),
        lwd = 2, col = "#B2182B", lty = 2)
  lines(density(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_itpc_ds), na.rm = T),
        lwd = 2, col = "#B2182B")
  abline(v=0, col = "grey60", lty = 3)
  text(0, 0.95, "prec-weighted d18O", col = "black", cex = 1.5)
  
  hist(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_temp), 
       breaks = 9, border = "white", prob = TRUE, 
       ylim = c(0,1), xlim = c(-2.5, 2.5), xlab = "",xaxt = 'n',
       cex.axis = 1.5, main = NULL)
  axis(side = 1, at = c(log10(0.01), log10(0.1), 0, log10(10), log10(100)), 
       labels = c(0.01, 0.1, 1, 10, 100), cex.axis = 1.5)
  lines(c(median(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_temp)),median(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_temp))), 
        c(0, max(density(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_temp), na.rm = T)$y)-0.01),
        lwd = 2, col = "black", lty = 2)
  lines(density(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_temp), na.rm = T),
        lwd = 2, col = "black")
  lines(c(median(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_temp_ds)),median(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_temp_ds))), 
        c(0, max(density(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_temp_ds), na.rm = T)$y)-0.01),
        lwd = 2, col = "#B2182B", lty = 2)
  lines(density(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_temp_ds), na.rm = T),
        lwd = 2, col = "#B2182B")
  abline(v=0, col = "grey60", lty = 3)
  #mtext(text = TeX("$$Var_{Rec}/Var_{Sim}$$"),side = 1,line = 2.5, cex = 1.5)
  mtext(text = "density",side = 2,line = 2.5, cex = 1.5)
  text(0, 0.95, "Temperature", col = "black", cex = 1.5)
  mtext(text = TeX("$$Var_{Rec}/Var_{Sim}$$"),side = 1,line = 2.7, cex = 1.5)
  
  hist(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_prec), 
       breaks = 9, border = "white", prob = TRUE, 
       ylim = c(0,1), xlim = c(-2.5, 2.5), xlab = "",xaxt = 'n',
       main = "", cex.main = 1.5, cex.axis = 1.5)
  axis(side = 1, at = c(log10(0.01), log10(0.1), 0, log10(10), log10(100)), 
       labels = c(0.01, 0.1, 1, 10, 100), cex.axis = 1.5)
  lines(c(median(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_prec)),median(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_prec))), 
        c(0, max(density(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_prec), na.rm = T)$y)-0.01),
        lwd = 2, col = "black", lty = 2)
  lines(density(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_prec), na.rm = T),
        lwd = 2, col = "black")
  lines(c(median(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_prec_ds)),median(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_prec_ds))), 
        c(0, max(density(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_prec_ds), na.rm = T)$y)-0.01),
        lwd = 2, col = "#B2182B", lty = 2)
  lines(density(log10(PLOTTING_DATA$FIG4$VAR_RATIOS[[run]]$VR_sim_prec_ds), na.rm = T),
        lwd = 2, col = "#B2182B")
  abline(v=0, col = "grey60", lty = 3)
  text(0, 0.95, "Precipitation", col = "black", cex = 1.5)
  mtext(text = TeX("$$Var_{Rec}/Var_{Sim}$$"),side = 1,line = 2.7, cex = 1.5)
  #mtext(text = TeX("$$Var_{Rec}/Var_{Sim}$$"),side = 1,line = 2.5, cex = 1.5)
  
  dev.off()
}



remove(plot_var)


