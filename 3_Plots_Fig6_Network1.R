#################################################
## Paper Figure 6 ###############################
#################################################

## Here analysis and Plotting

## Network Plot 1

#################################################

library(plyr)
library(dplyr)
library(tidyverse)
library(zoo)

## USE nest::network_links

lats = c()
longs = c()

#counter = 1
for (entity in DATA_past1000$CAVES$entity_info$entity_id[mask_spec]){
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  #prep_corr_matrix[,counter] = DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ISOT
  #counter = counter + 1
  lats = c(lats, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
  longs = c(longs, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
}

load("ENSEMBLE.RData")
c_ensemble <- ENSEMBLE$NETWORK$C
rm(ENSEMBLE)
load("C_ensemble_max.RData")
c_ensemble <- C


###################################################################################################

source("Functions/Plotting/networkmap_simple3.R")

plot_dist <- ANALYSIS$NETWORK$DIST
plot_dist[lower.tri(ANALYSIS$NETWORK$DIST)] = NA
lowess_dist <- as.vector(ANALYSIS$NETWORK$DIST[upper.tri(ANALYSIS$NETWORK$DIST)])
o <- order(lowess_dist)
lowess_dist_sorted <- lowess_dist[o]

for(run in c("a","b","c")){
  link_density = 0.05

  C_SIM_p <- ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$C
  C_REC_p <- ANALYSIS$NETWORK$GLOBAL$C
  
  C_SIM_p[ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$P>0.1] = NA
  C_REC_p[ANALYSIS$NETWORK$GLOBAL$P>0.1] = NA
  
  o_sim = order(abs(C_SIM_p), na.last = F)
  o_rec = order(abs(C_REC_p), na.last = F)
  C_SIM_p[o_sim[1:floor((length(o_sim)-link_density*length(o_sim)))]]<- NA
  C_REC_p[o_rec[1:floor((length(o_rec)-link_density*length(o_rec)))]]<- NA
  
  plot_c_sim <- ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$C
  plot_c_sim[lower.tri(ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$C, diag = FALSE)] = NA
  lowess_c_sim <- as.vector(ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$C[upper.tri(ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$C)])
  lowess_c_sim_sorted <- lowess_c_sim[o]
  plot_c_rec <- ANALYSIS$NETWORK$GLOBAL$C
  plot_c_rec[lower.tri(ANALYSIS$NETWORK$GLOBAL$C)] = NA
  lowess_c_rec <- as.vector(ANALYSIS$NETWORK$GLOBAL$C[upper.tri(ANALYSIS$NETWORK$GLOBAL$C)])
  lowess_c_rec_sorted <- lowess_c_rec[o]
  plot_c_rec_max <- ANALYSIS$NETWORK$GLOBAL_CHRONO$C_max
  plot_c_rec_max[lower.tri(ANALYSIS$NETWORK$GLOBAL_CHRONO$C_max)] = NA
  lowess_c_rec_max <- as.vector(ANALYSIS$NETWORK$GLOBAL_CHRONO$C_max[upper.tri(ANALYSIS$NETWORK$GLOBAL_CHRONO$C_max)])
  lowess_c_rec_max_sorted <- lowess_c_rec_max[o]
  
  plot_c_rec_max_e <- c_ensemble
  plot_c_rec_max_e[lower.tri(c_ensemble)] = NA
  lowess_c_rec_max_e <- as.vector(c_ensemble[upper.tri(c_ensemble)])
  lowess_c_rec_max_e_sorted <- lowess_c_rec_max_e[o]
  
  lo <- loess(lowess_c_rec_sorted ~ lowess_dist_sorted)
  
  boxes_sim <- list()
  boxes_rec <- list()
  
  for(ii in 1:20){
    boxes_sim[[paste0(ii*1000)]] <- na.omit(as.numeric(plot_c_sim[plot_dist<ii*1000 & plot_dist>(ii-1)*1000]))
    boxes_rec[[paste0(ii*1000)]] <- na.omit(as.numeric(plot_c_rec[plot_dist<ii*1000 & plot_dist>(ii-1)*1000]))
  }
  
  scaling = 1.5
  spacing = 0.7
  namcex = 1
  colorbar_length = 0.7

  # vielleicht erst die stärksten Links aussuchen udn dann unsignifikante raus schmeißen! Damit die gleiche Anzahl an Links da ist!
  rbPal <- colorRampPalette(c("#2166ac", "grey", "#b2182b"))
  COLZ <- array(rbPal(9))
  
  
  pdf(file = paste0("Plots/Paper_Plot_6_Network_a_xnap",run,".pdf"), height= PLOTTING_VARIABLES$HEIGHT, width = PLOTTING_VARIABLES$WIDTH)
  #png(file = paste0("Plots/Paper_Plot_6_Network_a_xnap",run,".png"), height= 100*PLOTTING_VARIABLES$HEIGHT, width = 100*PLOTTING_VARIABLES$WIDTH)
  par(mfrow=c(2,2), mai = c(rep(spacing, 4)), mar = c(1.5,2.5,0,0.5), oma = c(2,0.5,0.5,0.5))
  #SIM MAP
  networkmap_simple3(CMAT = C_SIM_p, 
                     lat = lats, 
                     lon = longs, 
                     title = "",
                     thresh = 0.1)
  fields::colorbar.plot(x = 190,y = 10, col = rev(COLZ), 
                        strip = c(1,0.75, 0.5, 0.25,0,-0.25, -0.5, -0.75,-1), horizontal = F, strip.length = colorbar_length, strip.width = 0.04)
  axis(4,at=seq(80,-60,by=-2*140/8),labels=FALSE)
  mtext("HadCM3", side = 1, cex = namcex, line = -2, font = 2)
  mtext("(a)", side = 3, adj = 0, cex = namcex, line = -1.5, at = -185)
  #SIM Cor-Dist
  plot(plot_dist[seq(1,length(plot_dist), by = 2)], plot_c_sim[seq(1,length(plot_c_sim), by = 2)], 
       ylim = c(-1,1),
       xlim = c(0,20000),
       ylab = "",
       xlab = "", 
       cex = 1, 
       lwd = 0.5, 
       panel.first = grid(), col = adjustcolor("grey", alpha.f = 0.7), xaxt = "n")
  abline(h=0)
  for(ii in 1:20){
    boxplot(boxes_sim[[paste0(ii*1000)]], add = TRUE, at = c(ii*1000-500),boxwex = 1000, names = "n", axes = F, outline = F)  
  }
  
  lo <- loess(lowess_c_sim_sorted ~ lowess_dist_sorted, span = 0.2)
  
  lines(lo$x, lo$fitted, lwd = 4, col = "#B2182B")
  #lines(lowess(lowess_dist_sorted,lowess_c_sim_sorted, f=0.1), lwd = 4, col = "#B2182B")
  mtext("(b)", side = 3, adj = 0, cex = namcex, line = -1.5, at = 1000)
  
  #SISAL MAP
  networkmap_simple3(CMAT = C_REC_p, 
                     lat = lats, 
                     lon = longs,
                     title = "",
                     thresh = 0.1)
  fields::colorbar.plot(x = 190,y = 10, col = rev(COLZ), 
                        strip = c(1,0.75, 0.5, 0.25,0,-0.25, -0.5, -0.75,-1), horizontal = F, strip.length = colorbar_length, strip.width = 0.04)
  axis(4,at=seq(80,-60,by=-2*140/8),labels=FALSE)
  mtext("SISAL", side = 1, cex = namcex, line = -2, font = 2)
  mtext("(c)", side = 3, adj = 0, cex = namcex, line = -1.5, at = -185)
  
  #SISAL Cor-Dist
  plot(plot_dist[seq(1,length(plot_dist), by = 2)], plot_c_rec[seq(1,length(plot_c_rec), by = 2)], 
       ylim = c(-1,1),
       xlim = c(0,20000),
       ylab = "",
       xlab = "", 
       cex = 1, 
       lwd = 1, 
       panel.first = grid(), col = adjustcolor("grey", alpha.f = 0.7))
  abline(h=0)
  for(ii in 1:20){
    boxplot(boxes_rec[[paste0(ii*1000)]], add = TRUE, at = c(ii*1000-500),boxwex = 1000, names = "n", axes = F, outline = F)  
  }
  lo <- loess(lowess_c_rec_sorted ~ lowess_dist_sorted, span = 0.2)
  lines(lo$x, lo$fitted, lwd = 4, col = "#B2182B")
  lo <- loess(lowess_c_rec_max_sorted ~ lowess_dist_sorted, span = 0.2)
  lines(lo$x, lo$fitted, lwd = 4, col = "#1D3461")
  lo <- loess(lowess_c_rec_max_e_sorted ~ lowess_dist_sorted, span = 0.2)
  lines(lo$x, lo$fitted, lwd = 4, col = "#0068C4")
  mtext("Distance between pairs (km)", side= 1, line = 2)
  mtext("(d)", side = 3, adj = 0, cex = namcex, line = -1.5, at = 0)
  text(20000, 0.9, "original chron.", col = "#B2182B", adj = 1, font = 2)
  text(20000, 0.77, "sisal chron.", col = "#1D3461", adj = 1, font = 2)
  text(20000, 0.64, "sisal ensemble", col = "#0068C4", adj = 1, font = 2)
  
  dev.off()
}

rm(C, C_rec, C_REC_p, C_sim, C_SIM_p, boxes_rec, boxes_sim, COLZ)
rm(lo, mat_comp, network_corplot, network_lyr_rec, network_lyr_sim, network_p, network_sim, order, P, P_rec, P_sim, PLOT, plot_c_rec, plot_c_rec_max)
rm(plot_dist, Plot_Lyr, point_lyr, ptlyr_shina, ptlyr_china_p, ptlyr_rest, ptlyr_rest_p, s, sim_p)
rm(site_list, temp_sim, TEST, TS, TS_rec, TS_sim, chronology, cluster, corr, counter, entity_list, gridbox, i, ii, index, j, lats, line_x, line_y, link_density)
rm(longs, lowess_c_rec, lowess_c_rec_max, lowess_c_rec_max_sorted, lowess_c_rec_sorted, lowess_c_sim, lowess_c_sim_sorted, lowess_dist, lowess_dist_sorted)
rm(namcex, name, noPts, o, o_rec, o_sim, plot, scaling, seasons, season_num, site_corr, spacing, test, title, var, year_start, year_stop)
rm(box_corr, box_map, data, dating_all, double_time, list, map_window, names_move_plot_c_sim, plot_data, ptlyr_china)
rm(plot_c_sim, names_move)
