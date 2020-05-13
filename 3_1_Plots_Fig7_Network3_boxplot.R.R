#################################################
## DISCUSSION NETWORK TUNING ####################
#################################################

library(plyr)
library(dplyr)
library(tidyverse)
library(zoo)
library(PaleoSpec)
library(nest)
library(latex2exp)

DISCUSSION$NETWORK_TUNING <- list()

#################################################

# was brauch ich?
#   - masken für 50% geringster offset
#   - masken für cluster
#   - masken innerhalb cluster mit nur 50% geringsten Abstand
#       --> dafür pro cluster eine distance-maske

# SITE
mask = matrix(logical(length = 89*89), ncol = 89)
sites = ANALYSIS$NETWORK$entity_meta %>% select(entity_id, site_id) %>% group_by(site_id) %>% count() %>% filter(n>1)

for(site in sites$site_id){
  entity_total <- ANALYSIS$NETWORK$entity_meta$entity_id
  entity_list = ANALYSIS$NETWORK$entity_meta$entity_id[ANALYSIS$NETWORK$entity_meta$site_id == site]
  
  index.c <- length(entity_list)
  for(ii in 1:length(entity_list)){
    index.c[ii] = which(entity_list[ii] == entity_total)
  }
  
  for(ii in 1:(length(index.c)-1)){
    for(jj in ii:length(index.c)){
      mask[index.c[ii],index.c[jj]] = T
      mask[index.c[jj],index.c[ii]] = T
      mask[index.c[ii],index.c[ii]] = T
      mask[index.c[jj],index.c[jj]] = T
    }
  }
}

DISCUSSION$NETWORK_TUNING$SITES <- list()
DISCUSSION$NETWORK_TUNING$SITES$C_rec <- matrix(ANALYSIS$NETWORK$GLOBAL$C[mask])
DISCUSSION$NETWORK_TUNING$SITES$P_rec <- matrix(ANALYSIS$NETWORK$GLOBAL$P[mask])
DISCUSSION$NETWORK_TUNING$SITES$C_rec_gauss <- matrix(ANALYSIS$NETWORK$GLOBAL$C_gauss[mask])
DISCUSSION$NETWORK_TUNING$SITES$P_rec_gauss <- matrix(ANALYSIS$NETWORK$GLOBAL$P_gauss[mask])

DISCUSSION$NETWORK_TUNING$SITES$C_rec_ensemble <- matrix(c_ensemble[mask])
DISCUSSION$NETWORK_TUNING$SITES$C_rec_chrono <- matrix(ANALYSIS$NETWORK$GLOBAL_CHRONO$C_max[mask])



# GRIDBOX 

mask = matrix(logical(length = 89*89), ncol = 89)
gridboxes = ANALYSIS$NETWORK$entity_meta %>% select(entity_id, gridbox_id) %>% group_by(gridbox_id) %>% count() %>% filter(n>1)

for(gridbox in gridboxes$gridbox_id){
  entity_total <- ANALYSIS$NETWORK$entity_meta$entity_id
  entity_list = ANALYSIS$NETWORK$entity_meta$entity_id[ANALYSIS$NETWORK$entity_meta$gridbox_id == gridbox]
  
  index.c <- length(entity_list)
  for(ii in 1:length(entity_list)){
    index.c[ii] = which(entity_list[ii] == entity_total)
  }
  
  for(ii in 1:(length(index.c)-1)){
    for(jj in ii:length(index.c)){
      mask[index.c[ii],index.c[jj]] = T
      mask[index.c[jj],index.c[ii]] = T
      mask[index.c[ii],index.c[ii]] = T
      mask[index.c[jj],index.c[jj]] = T
    }
  }
}

DISCUSSION$NETWORK_TUNING$GRIDBOX <- list()
DISCUSSION$NETWORK_TUNING$GRIDBOX$C_rec <- matrix(ANALYSIS$NETWORK$GLOBAL$C[mask])
DISCUSSION$NETWORK_TUNING$GRIDBOX$P_rec <- matrix(ANALYSIS$NETWORK$GLOBAL$P[mask])
DISCUSSION$NETWORK_TUNING$GRIDBOX$C_rec_gauss <- matrix(ANALYSIS$NETWORK$GLOBAL$C_gauss[mask])
DISCUSSION$NETWORK_TUNING$GRIDBOX$P_rec_gauss <- matrix(ANALYSIS$NETWORK$GLOBAL$P_gauss[mask])

DISCUSSION$NETWORK_TUNING$GRIDBOX$C_rec_ensemble <- matrix(c_ensemble[mask])
DISCUSSION$NETWORK_TUNING$GRIDBOX$C_rec_chrono <- matrix(ANALYSIS$NETWORK$GLOBAL_CHRONO$C_max[mask])

# CLUSTER

DISCUSSION$NETWORK_TUNING$CLUSTER <- list()

for(cluster in 1:9){
  
  DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER", cluster)]] <- list()
    
  entity_total <- ANALYSIS$NETWORK$entity_meta$entity_id
  cluster_list <- ANALYSIS$NETWORK$entity_meta$entity_id[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster]
  index.c <- length(cluster_list)
  for(ii in 1:length(cluster_list)){
    index.c[ii] = which(cluster_list[ii] == entity_total)
  }
  mask = matrix(logical(length = 89*89), ncol = 89)
  
  for(ii in 1:(length(index.c)-1)){
    for(jj in ii:length(index.c)){
      mask[index.c[ii],index.c[jj]] = T
      mask[index.c[jj],index.c[ii]] = T
      mask[index.c[ii],index.c[ii]] = T
      mask[index.c[jj],index.c[jj]] = T
    }
  }
  
  DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$C_rec <- matrix(ANALYSIS$NETWORK$GLOBAL$C[mask], ncol = length(cluster_list))
  DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$P_rec <- matrix(ANALYSIS$NETWORK$GLOBAL$P[mask], ncol = length(cluster_list))
  DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$C_rec_gauss <- matrix(ANALYSIS$NETWORK$GLOBAL$C_gauss[mask], ncol = length(cluster_list))
  DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$P_rec_gauss <- matrix(ANALYSIS$NETWORK$GLOBAL$P_gauss[mask], ncol = length(cluster_list))
  
  DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$C_rec_ensemble <- matrix(c_ensemble[mask], ncol = length(cluster_list))
  DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$C_rec_chrono <- matrix(ANALYSIS$NETWORK$GLOBAL_CHRONO$C_max[mask], ncol = length(cluster_list))
  
  ## mask 50% closest within cluster
  mask_small = matrix(logical(length = length(cluster_list)^2), ncol = length(cluster_list))
  
  dist = matrix(ANALYSIS$NETWORK$DIST[mask], ncol = length(cluster_list))
  diag(dist) = NA
  mask_small[dist<=median(dist,na.rm = T)] = TRUE
  
  DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$closest <- list()
  DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$closest$C_rec <- ANALYSIS$NETWORK$GLOBAL$C[mask][mask_small]
  DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$closest$P_rec <- ANALYSIS$NETWORK$GLOBAL$P[mask][mask_small]
  DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$closest$C_rec_gauss <- ANALYSIS$NETWORK$GLOBAL$C_gauss[mask][mask_small]
  DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$closest$P_rec_gauss <- ANALYSIS$NETWORK$GLOBAL$P_gauss[mask][mask_small]
  
  DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$closest$C_rec_ensemble <- c_ensemble[mask][mask_small]
  DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$closest$C_rec_chrono <- ANALYSIS$NETWORK$GLOBAL_CHRONO$C_max[mask][mask_small]
  
  ## mask 50% smallest offset
  entity_total <- ANALYSIS$NETWORK$entity_meta$entity_id
  cluster_list <- ANALYSIS$NETWORK$entity_meta$entity_id[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster]
  index.c <- length(cluster_list)
  offset = scatter_data %>% filter(entity_id %in% cluster_list)
  offset = offset$diff_down
  cluster_list[abs(offset)>median(abs(offset), na.rm = T)] = NA
  cluster_list = na.omit(cluster_list)
  for(ii in 1:length(cluster_list)){
    index.c[ii] = which(cluster_list[ii] == entity_total)
  }
  mask = matrix(logical(length = 89*89), ncol = 89)
  
  for(ii in 1:(length(index.c)-1)){
    for(jj in ii:length(index.c)){
      mask[index.c[ii],index.c[jj]] = T
      mask[index.c[jj],index.c[ii]] = T
      mask[index.c[ii],index.c[ii]] = T
      mask[index.c[jj],index.c[jj]] = T
    }
  }

  
  DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$smallest_offset <- list()
  DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$smallest_offset$C_rec <- ANALYSIS$NETWORK$GLOBAL$C[mask]
  DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$smallest_offset$P_rec <- ANALYSIS$NETWORK$GLOBAL$P[mask]
  DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$smallest_offset$C_rec_gauss <- ANALYSIS$NETWORK$GLOBAL$C_gauss[mask]
  DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$smallest_offset$P_rec_gauss <- ANALYSIS$NETWORK$GLOBAL$P_gauss[mask]
  
  DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$smallest_offset$C_rec_ensemble <- c_ensemble[mask]
  DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$smallest_offset$C_rec_chrono <- ANALYSIS$NETWORK$GLOBAL_CHRONO$C_max[mask]
  
}

#################################################
## EXTENDED BOXPLOT #############################
#################################################

#sites = 4
#gridbox = 4
#cluster = 7*6 = 42
#global = 6
 
#total = 56+3 = 59

line_names = 1
col_chrono = "#C47900"
col_ensemble = "#004F00"

col_closest <- c("#6A51A3", "#9E9AC8", "#CBC9E2", "#F2E6F7")
col_offset <- c("#CB181D", "#FB6A4A", "#FCAE91", "#FEE5D9")
for(plot in 1:1){
  pdf(file = paste0("Plots/Paper_Plot_6_Network_c_extended.pdf"), height= PLOTTING_VARIABLES$HEIGHT, width = PLOTTING_VARIABLES$WIDTH)
  par(mar = c(2,11,4,2))
  plot(c(-1,1), c(1,28), type = "n", axes = FALSE, xlab = "", ylab = "" )
  abline(v=0)
  ## SITES
  corr <- list(full     = DISCUSSION$NETWORK_TUNING$SITES$C_rec[DISCUSSION$NETWORK_TUNING$SITES$P_rec<0.1], 
               gauss    = DISCUSSION$NETWORK_TUNING$SITES$C_rec_gauss[DISCUSSION$NETWORK_TUNING$SITES$P_rec_gauss<0.1],
               ensemble = DISCUSSION$NETWORK_TUNING$SITES$C_rec_ensemble)
  
  boxplot(as.numeric(corr$full), add = TRUE, at = 28, boxwex = 1, names = "n", horizontal = T, outline = F) 
  boxplot(as.numeric(corr$gauss), add = TRUE, at = 27.5, boxwex = 1, names = "n", horizontal = T, col = "grey", axes = F, outline = F)
  boxplot(as.numeric(corr$ensemble), add = TRUE, at = 27, boxwex = 1, names = "n", horizontal = T, col = "#004F00", axes = F, outline = F)
  
  
  
  #mtext(side=2,"GMST",                cex = unitscex,    line = unitslinno, las = 1, col = "black", at = 1)
  mtext(side = 2, "sites (27/12)", cex = 1, line = line_names, las = 1, col = "black", at = 28)
  abline(h=26.5, lty = 3, lwd = 2)
  ##GRIGBOX
  
  corr <- list(full     = DISCUSSION$NETWORK_TUNING$GRIDBOX$C_rec[DISCUSSION$NETWORK_TUNING$GRIDBOX$P_rec<0.1], 
               gauss    = DISCUSSION$NETWORK_TUNING$GRIDBOX$C_rec_gauss[DISCUSSION$NETWORK_TUNING$GRIDBOX$P_rec_gauss<0.1],
               ensemble = DISCUSSION$NETWORK_TUNING$GRIDBOX$C_rec_ensemble)
  
  boxplot(as.numeric(corr$full), add = TRUE, at = 26, boxwex = 1, names = "n", horizontal = T, outline = F) 
  boxplot(as.numeric(corr$gauss), add = TRUE, at = 25.5, boxwex = 1, names = "n", horizontal = T, col = "grey", axes = F, outline = F)
  boxplot(as.numeric(corr$ensemble), add = TRUE, at = 25, boxwex = 1, names = "n", horizontal = T, col = "#004F00", axes = F, outline = F)
  
  #mtext(side=2,"GMST",                cex = unitscex,    line = unitslinno, las = 1, col = "black", at = 1)
  mtext(side = 2, "gridbox (45/18)", cex = 1, line = line_names, las = 1, col = "black", at = 26)
  abline(h=24.5, lty = 3, lwd = 2)
  
  
  position <- list(cluster6 = c(24, 23.5, 23, 22.5,22),
                   cluster2 = c(21, 20.5, 20, 19.5,19),
                   cluster3 = c(18, 17.5, 17, 16.5,16),
                   cluster7 = c(15, 14.5, 14, 13.5,13),
                   cluster1 = c(12, 11.5, 11, 10.5,10),
                   cluster5 = c( 9,  8.5,  8,  7.5,7),
                   cluster9 = c( 6, 5.5, 5, 4.5, 4))
  text <- list(cluster1 = "India + Centr.Asia", cluster2 = "SouthAm", cluster3 = "Europe", cluster4 = "Africa", 
               cluster5 = "China + East.Asia", cluster6 = "NorthAm", cluster7 = "Arabia", cluster8 = "NZ", cluster9 = "SE Asia")
  cluster_list <- ANALYSIS$NETWORK$entity_meta %>% select(cluster_id, entity_id) %>% group_by(cluster_id) %>% count() %>% filter(n>1)
  cluster_number <- c(6,2,3,4,7,1,5,9,8)
  
  ## CLUSTER
  #for(cluster in c(1,2,3,5,6,7)){
  for(cluster in c(6,2,3,7,1,5,9)){
    corr <- list(full = as.numeric(DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$C_rec[DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$P_rec<0.1]), 
                 gauss = as.numeric(DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$C_rec_gauss[DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$P_rec_gauss<0.1]), 
                 ensemble = as.numeric(DISCUSSION$NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$C_rec_ensemble),
                 
                 full_sim = list(as.numeric(c(ANALYSIS$NETWORK$CLUSTER_SIM_a[[paste0("CLUSTER",cluster)]]$C[ANALYSIS$NETWORK$CLUSTER_SIM_a[[paste0("CLUSTER",cluster)]]$P<0.1],
                                              ANALYSIS$NETWORK$CLUSTER_SIM_b[[paste0("CLUSTER",cluster)]]$C[ANALYSIS$NETWORK$CLUSTER_SIM_b[[paste0("CLUSTER",cluster)]]$P<0.1],
                                              ANALYSIS$NETWORK$CLUSTER_SIM_c[[paste0("CLUSTER",cluster)]]$C[ANALYSIS$NETWORK$CLUSTER_SIM_c[[paste0("CLUSTER",cluster)]]$P<0.1]))), 
                 gauss_sim = list(as.numeric(c(ANALYSIS$NETWORK$CLUSTER_SIM_a[[paste0("CLUSTER",cluster)]]$C_gauss[ANALYSIS$NETWORK$CLUSTER_SIM_a[[paste0("CLUSTER",cluster)]]$P_gauss<0.1],
                                               ANALYSIS$NETWORK$CLUSTER_SIM_b[[paste0("CLUSTER",cluster)]]$C_gauss[ANALYSIS$NETWORK$CLUSTER_SIM_b[[paste0("CLUSTER",cluster)]]$P_gauss<0.1],
                                               ANALYSIS$NETWORK$CLUSTER_SIM_c[[paste0("CLUSTER",cluster)]]$C_gauss[ANALYSIS$NETWORK$CLUSTER_SIM_c[[paste0("CLUSTER",cluster)]]$P_gauss<0.1]))))
    
    boxplot(corr$full, add = TRUE, at = position[[paste0("cluster",cluster)]][1] ,boxwex = 1, names = "n", horizontal = T, axes = F, outline = F) 
    boxplot(corr$gauss, add = TRUE, at = position[[paste0("cluster",cluster)]][2] , boxwex = 1, names = "n", horizontal = T, col = "grey", axes = F, outline = F)
    boxplot(corr$ensemble, add = TRUE, at = position[[paste0("cluster",cluster)]][3] , boxwex = 1, names = "n", horizontal = T, col = col_ensemble, axes = F, outline = F)
    
    boxplot(corr$full_sim, add = TRUE, at = position[[paste0("cluster",cluster)]][4]  ,boxwex = 1, names = "n", horizontal = T, col = "dodgerblue3", axes = F, outline = F) 
    boxplot(corr$gauss_sim, add = TRUE, at = position[[paste0("cluster",cluster)]][5] , boxwex = 1, names = "n", horizontal = T, col = adjustcolor("dodgerblue3", alpha.f = 0.5), axes = F, outline = F)
    
    mtext(side = 2, paste0("c",cluster_number[cluster],"/",text[[cluster]]," [",cluster_list$n[cluster],"]"), cex = 1, line = line_names, las = 1, col = "black", at = position[[paste0("cluster",cluster)]][1])
    
  }
  
  abline(h=3.5, lty =3, lwd = 2)
  
  ##Global
  
  corr$full = as.numeric(as.numeric(ANALYSIS$NETWORK$GLOBAL$C[ANALYSIS$NETWORK$GLOBAL$P<0.1]))
  corr$gauss = as.numeric(as.numeric(ANALYSIS$NETWORK$GLOBAL$C_gauss[ANALYSIS$NETWORK$GLOBAL$P_gauss<0.1]))
  corr$ensembls = as.numeric(c_ensemble)
  corr$full_sim = as.numeric(ANALYSIS$NETWORK$GLOBAL_SIM_a$C[ANALYSIS$NETWORK$GLOBAL_SIM_a$P<0.1], 
                             ANALYSIS$NETWORK$GLOBAL_SIM_b$C[ANALYSIS$NETWORK$GLOBAL_SIM_b$P<0.1],
                             ANALYSIS$NETWORK$GLOBAL_SIM_c$C[ANALYSIS$NETWORK$GLOBAL_SIM_c$P<0.1])
  corr$gauss_sim = as.numeric(ANALYSIS$NETWORK$GLOBAL_SIM_a$C_gauss[ANALYSIS$NETWORK$GLOBAL_SIM_a$P_gauss<0.1], 
                              ANALYSIS$NETWORK$GLOBAL_SIM_b$C_gauss[ANALYSIS$NETWORK$GLOBAL_SIM_b$P_gauss<0.1],
                              ANALYSIS$NETWORK$GLOBAL_SIM_c$C_gauss[ANALYSIS$NETWORK$GLOBAL_SIM_c$P_gauss<0.1])
  
  boxplot(corr$full, add = TRUE, at = 3 ,boxwex = 1, names = "n", horizontal = T, axes = F) 
  boxplot(corr$gauss, add = TRUE, at = 2.5 , boxwex = 1, names = "n", horizontal = T, col = "grey", axes = F)
  boxplot(corr$ensemble, add = TRUE, at = 2 , boxwex = 1, names = "n", horizontal = T, col = col_ensemble, axes = F)
  boxplot(corr$full_sim, add = TRUE, at = 1.5,boxwex = 1, names = "n", horizontal = T, col = "dodgerblue3", axes = F, outline = F) 
  boxplot(corr$gauss_sim, add = TRUE, at = 1 , boxwex = 1, names = "n", horizontal = T, col = adjustcolor("dodgerblue3", alpha.f = 0.5), axes = F, outline = F)
  
  
  mtext(side = 2, paste0("global (",dim(ANALYSIS$NETWORK$GLOBAL$C)[[1]],")"), cex = 1, line = line_names, las = 1, col = "black", at = 3)
  
  
  legend(-1.075,34,xpd = T,inset=-0.2, bty='n', x.intersp=0.5,text.width=c(1,0.3,0.4),
         c("record", "record 100gauss", "record best selection"), fill = c("white", "grey", col_ensemble), horiz = TRUE, cex = 1)
  legend(-1.075, 32,xpd = T,inset=-0.2, bty='n', x.intersp=0.5,text.width=c(1,0.3),
         c("simulation", "simulation 100gauss"), fill = c("dodgerblue3", adjustcolor("dodgerblue3", alpha.f = 0.5)), horiz=TRUE, cex=1)
  
  
  dev.off()
  
}

rm(corr, dist, dist_matrix, cluster_list, gridboxes, mask, mask_small, position, sites)
rm(d18O, data, DATA_EXPORT_ENTITY, DATA_EXPORT_SITE, data_new, double_time, ens_1000, filter, mask_cluster, mask_cluster1, org_chronology, PLOT, plot_c_rec_max_e, Plot_Lyr,
   ptlyr_china, ptlyr_china_p, ptlyr_rest, ptlyr_rest_p, rm_list, s, sample, temp, TEST, text, TS)
rm(cluster, cluster_number, col_chrono, col_closest, col_ensemble, col_offset, entity, entity_list, entity_total, gridbox, ii, index, index.c, jj, line_names, offset, site)
rm(cex_axis, cex_axis_text, cex_text, chronology_2, colorbar_length, COLZ, counter, diff, e_count, e_count2, entity_2, gmd18O,lats, length, longs, lowess_c_rec_max_e, lowess_c_rec_max_e_sorted,
   name,noE, plot, record, Results, run, season, season_num,start_ts, stop_ts, title, TS_full_down, TS_down_rec, TS_full_rec, var)
