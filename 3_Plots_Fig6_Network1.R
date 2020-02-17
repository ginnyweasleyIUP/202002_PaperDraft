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


source("Functions/Plotting/networkmap_simple2.R")

networkmap_simple3(CMAT = ANALYSIS$NETWORK$GLOBAL$C, 
                   lat = lats, 
                   lon = longs, 
                   title = "Correlation HadCM3, sig level = 0.1", 
                   thresh = 0.6)



###################################################################################################

source("Functions/Plotting/networkmap_simple3.R")

plot_dist <- ANALYSIS$NETWORK$DIST
plot_dist[lower.tri(ANALYSIS$NETWORK$DIST)] = NA
lowess_dist <- as.vector(ANALYSIS$NETWORK$DIST[upper.tri(ANALYSIS$NETWORK$DIST)])
o <- order(lowess_dist)
lowess_dist_sorted <- lowess_dist[o]

for(run in c("a","b","c")){
  link_density = 0.15
  C_SIM_p <- ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$C
  C_REC_p <- ANALYSIS$NETWORK$GLOBAL$C
  
  o_sim = order(C_SIM_p, na.last = NA)
  C_SIM_p[o_sim[1:floor((length(o_sim)-link_density*length(C_SIM_p)))]]<- NA
  o_rec = order(C_REC_p, na.last = NA)
  C_REC_p[o_sim[1:floor((length(o_rec)-link_density*length(C_REC_p)))]]<- NA
  
  C_SIM_p[ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$P>0.1] = NA
  C_REC_p[ANALYSIS$NETWORK$GLOBAL$P>0.1] = NA
  
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
  
  lo <- loess(lowess_c_rec_sorted ~ lowess_dist_sorted)
  
  boxes_sim <- list()
  boxes_rec <- list()
  
  for(ii in 1:10){
    boxes_sim[[paste0(ii*2000)]] <- na.omit(as.numeric(plot_c_sim[plot_dist<ii*2000 & plot_dist>(ii-1)*2000]))
    boxes_rec[[paste0(ii*2000)]] <- na.omit(as.numeric(plot_c_rec[plot_dist<ii*2000 & plot_dist>(ii-1)*2000]))
  }
  
  scaling = 1.5
  spacing = 0.7
  namcex = 1

  # vielleicht erst die stärksten Links aussuchen udn dann unsignifikante raus schmeißen! Damit die gleiche Anzahl an Links da ist!
  
  
  
  #pdf(file = paste0("Plots/Paper_Plot_6_Network_a_xnap",run,".pdf"), height= PLOTTING_VARIABLES$HEIGHT, width = PLOTTING_VARIABLES$WIDTH)
  png(file = paste0("Plots/Paper_Plot_6_Network_a_xnap",run,".png"), height= 100*PLOTTING_VARIABLES$HEIGHT, width = 100*PLOTTING_VARIABLES$WIDTH)
  par(mfrow=c(2,2), mai = c(rep(spacing, 4)), mar = c(3,3,2,0.5))
  #SIM MAP
  networkmap_simple3(CMAT = C_SIM_p, 
                     lat = lats, 
                     lon = longs, 
                     title = "",
                     thresh = 0.1)
  mtext("Corr HadCM3, p<0.1, c>0.1", side = 3, cex = namcex)
  mtext("A", side = 3, adj = 0, cex = namcex)
  #SIM Cor-Dist
  plot(plot_dist, plot_c_sim, 
       ylim = c(-1,1),
       xlim = c(0,20000),
       ylab = "",
       xlab = "", 
       cex = 1, 
       lwd = 0.5, 
       panel.first = grid(), col = "grey", type = "n", xaxt = "n")
  for(ii in 1:10){
    boxplot(boxes_sim[[paste0(ii*2000)]], add = TRUE, at = c(ii*2000-1000),boxwex = 1000, names = "n")  
  }
  
  lo <- loess(lowess_c_sim_sorted ~ lowess_dist_sorted, span = 0.2)
  
  lines(lo$x, lo$fitted, lwd = 4, col = "#B2182B")
  #lines(lowess(lowess_dist_sorted,lowess_c_sim_sorted, f=0.1), lwd = 4, col = "#B2182B")
  mtext("Distance between pairs (km)", side= 1, line = 2)
  mtext("B", side = 3, adj = 0, cex = namcex)
  
  #SISAL MAP
  networkmap_simple3(CMAT = C_REC_p, 
                     lat = lats, 
                     lon = longs,
                     title = "",
                     thresh = 0.1)
  mtext("Corr SISAL , p<0.1, c>0.1", side = 3, cex = namcex)
  mtext("C", side = 3, adj = 0, cex = namcex)
  
  #SISAL Cor-Dist
  plot(plot_dist, plot_c_rec, 
       ylim = c(-1,1),
       xlim = c(0,20000),
       ylab = "",
       xlab = "", 
       cex = 1, 
       lwd = 1, 
       panel.first = grid(), col = "grey", type = "n")
  for(ii in 1:10){
    boxplot(boxes_rec[[paste0(ii*2000)]], add = TRUE, at = c(ii*2000-1000),boxwex = 1000, names = "n")  
  }
  lo <- loess(lowess_c_rec_sorted ~ lowess_dist_sorted, span = 0.2)
  lines(lo$x, lo$fitted, lwd = 4, col = "#B2182B")
  lo <- loess(lowess_c_rec_max_sorted ~ lowess_dist_sorted, span = 0.2)
  lines(lo$x, lo$fitted, lwd = 4, col = "darkblue")
  mtext("Distance between pairs (km)", side= 1, line = 2)
  mtext("D", side = 3, adj = 0, cex = namcex)
  legend("topright", legend = c("original chron", "sisal chron."), col = c("#B2182B", "darkblue"), lty = c(1,1), lwd = c(4,4))
  
  
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
