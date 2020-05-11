#################################################
## DISCUSSION ANALYSIS ##########################
#################################################

library(plyr)
library(dplyr)
library(tidyverse)

DISCUSSION <- list()

#################################################
## VARIANCE SCATTER #############################

# Gibt es wirklich keinen Zusammenhang zwischen Offset und Varianz?

DISCUSSION$scatter_data_variance =array(dim = c(length(DATA_past1000$CAVES$entity_info$entity_id[mask_var]),19))

colnames(DISCUSSION$scatter_data_variance) = c("entity_id", "site_id", "diff_full", "diff_down", "vr_full_isot", "vr_full_itpc", "vr_ds_isot", "vr_ds_itpc",
                                               "elevation", "latitude", "mean_temp", "mean_prec", "winter_mean_prec", "elevation_diff", "dist_entrance", 
                                               "geology", "cover_thickness", "d18Oc", "var_record")

DISCUSSION$scatter_data_variance[,1] = DATA_past1000$CAVES$entity_info$entity_id[mask_var]
for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])){
  entity = DISCUSSION$scatter_data_variance[[ii,1]]
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  index = match(site, DATA_past1000$CAVES$site_info$site_id)
  DISCUSSION$scatter_data_variance[ii,2] = as.numeric(site)
  DISCUSSION$scatter_data_variance[ii,3] = as.numeric(mean(c(mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ISOT_a, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a, na.rm = T),
                                                             mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ISOT_b, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_b, na.rm = T),
                                                             mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ISOT_c, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_c, na.rm = T))))
  DISCUSSION$scatter_data_variance[ii,4] = as.numeric(mean(c(mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$ISOT_a, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a, na.rm = T),
                                                             mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$ISOT_b, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_b, na.rm = T),
                                                             mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$ISOT_c, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_c, na.rm = T))))
  #Var Ratio
  DISCUSSION$scatter_data_variance[ii,5] = as.numeric(mean(c(ANALYSIS$VARIANCE$POINTS$CAVElyr$a$value_VR_isot[ii],
                                                             ANALYSIS$VARIANCE$POINTS$CAVElyr$b$value_VR_isot[ii],
                                                             ANALYSIS$VARIANCE$POINTS$CAVElyr$c$value_VR_isot[ii]), na.rm = T))
  DISCUSSION$scatter_data_variance[ii,6] = as.numeric(median(c(ANALYSIS$VARIANCE$POINTS$CAVElyr$a$value_VR_itpc[ii],
                                                             ANALYSIS$VARIANCE$POINTS$CAVElyr$b$value_VR_itpc[ii],
                                                             ANALYSIS$VARIANCE$POINTS$CAVElyr$c$value_VR_itpc[ii]), na.rm = T))
  DISCUSSION$scatter_data_variance[ii,7] = as.numeric(median(c(ANALYSIS$VARIANCE$POINTS$CAVElyr$a$value_VR_isot_ds[ii],
                                                             ANALYSIS$VARIANCE$POINTS$CAVElyr$b$value_VR_isot_ds[ii],
                                                             ANALYSIS$VARIANCE$POINTS$CAVElyr$c$value_VR_isot_ds[ii]), na.rm = T))
  DISCUSSION$scatter_data_variance[ii,8] = as.numeric(median(c(ANALYSIS$VARIANCE$POINTS$CAVElyr$a$value_VR_itpc_ds[ii],
                                                             ANALYSIS$VARIANCE$POINTS$CAVElyr$b$value_VR_itpc_ds[ii],
                                                             ANALYSIS$VARIANCE$POINTS$CAVElyr$c$value_VR_itpc_ds[ii]), na.rm = T))
  
  #parameters
  DISCUSSION$scatter_data_variance[ii,9] = as.numeric(DATA_past1000$CAVES$site_info$elevation[index])
  DISCUSSION$scatter_data_variance[ii,10] = as.numeric(DATA_past1000$CAVES$site_info$latitude[index])
  DISCUSSION$scatter_data_variance[ii,11] = as.numeric(mean(c(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$TEMP_a,
                                                              DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$TEMP_b,
                                                              DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$TEMP_c), na.rm = T))
  DISCUSSION$scatter_data_variance[ii,12] = as.numeric(mean(c(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$PREC_a,
                                                              DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$PREC_b,
                                                              DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$PREC_c), na.rm = T))
  DISCUSSION$scatter_data_variance[ii,13] = as.numeric(mean(c(DATA_past1000$CAVES$sim_data_seasonal$a[[paste0("CAVE", site)]]$WINTER$prec_mean,
                                                              DATA_past1000$CAVES$sim_data_seasonal$b[[paste0("CAVE", site)]]$WINTER$prec_mean,
                                                              DATA_past1000$CAVES$sim_data_seasonal$c[[paste0("CAVE", site)]]$WINTER$prec_mean), na.rm = T))
  DISCUSSION$scatter_data_variance[ii,14] = as.numeric(DATA_past1000$CAVES$elevation_cave_sim$`sim-cave`[DATA_past1000$CAVES$elevation_cave_sim$entity_id == entity])
  DISCUSSION$scatter_data_variance[ii,15] = as.numeric(DATA_past1000$CAVES$entity_info$distance_entrance[match(entity, DATA_past1000$CAVES$entity_info$entity_id)])
  DISCUSSION$scatter_data_variance[ii,16] = DATA_past1000$CAVES$site_to_entity$geology[match(entity, DATA_past1000$CAVES$site_to_entity$entity_id)]
  DISCUSSION$scatter_data_variance[ii,17] = as.numeric(DATA_past1000$CAVES$site_to_entity$cover_thickness[match(entity, DATA_past1000$CAVES$site_to_entity$entity_id)])
  DISCUSSION$scatter_data_variance[ii,18] = as.numeric(mean(c(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a,
                                                              DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a,
                                                              DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a), na.rm = T))
  DISCUSSION$scatter_data_variance[ii,19] = as.numeric(var(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_measurement, na.rm = T))
}

DISCUSSION$scatter_data_variance = as.tibble(DISCUSSION$scatter_data_variance)
DISCUSSION$scatter_data_variance$site_id = as.numeric(DISCUSSION$scatter_data_variance$site_id)
DISCUSSION$scatter_data_variance$diff_full = as.numeric(DISCUSSION$scatter_data_variance$diff_full)
DISCUSSION$scatter_data_variance$diff_down = as.numeric(DISCUSSION$scatter_data_variance$diff_down)
DISCUSSION$scatter_data_variance$vr_full_isot = as.numeric(DISCUSSION$scatter_data_variance$vr_full_isot)
DISCUSSION$scatter_data_variance$vr_full_itpc = as.numeric(DISCUSSION$scatter_data_variance$vr_full_itpc)
DISCUSSION$scatter_data_variance$vr_ds_isot = as.numeric(DISCUSSION$scatter_data_variance$vr_ds_isot)
DISCUSSION$scatter_data_variance$vr_ds_itpc = as.numeric(DISCUSSION$scatter_data_variance$vr_ds_itpc)
DISCUSSION$scatter_data_variance$elevation = as.numeric(DISCUSSION$scatter_data_variance$elevation)
DISCUSSION$scatter_data_variance$latitude = as.numeric(DISCUSSION$scatter_data_variance$latitude)
DISCUSSION$scatter_data_variance$mean_temp = as.numeric(DISCUSSION$scatter_data_variance$mean_temp)
DISCUSSION$scatter_data_variance$mean_prec = as.numeric(DISCUSSION$scatter_data_variance$mean_prec)
DISCUSSION$scatter_data_variance$winter_mean_prec = as.numeric(DISCUSSION$scatter_data_variance$winter_mean_prec)
DISCUSSION$scatter_data_variance$elevation_diff = as.numeric(DISCUSSION$scatter_data_variance$elevation_diff)
DISCUSSION$scatter_data_variance$dist_entrance = as.numeric(DISCUSSION$scatter_data_variance$dist_entrance)
DISCUSSION$scatter_data_variance$cover_thickness = as.numeric(DISCUSSION$scatter_data_variance$cover_thickness)
DISCUSSION$scatter_data_variance$d18Oc = as.numeric(DISCUSSION$scatter_data_variance$d18Oc)
DISCUSSION$scatter_data_variance$var_record = as.numeric(DISCUSSION$scatter_data_variance$var_record)
#DISCUSSION$scatter_data_variance$winter_mean_prec[57] = NA

mask_mean_calcite = logical(length = length(DISCUSSION$scatter_data_variance$entity_id))
mask_mean_aragonite  = logical(length = length(DISCUSSION$scatter_data_variance$entity_id))

for(ii in 1:length(DISCUSSION$scatter_data_variance$entity_id)){
  entity = DISCUSSION$scatter_data_variance$entity_id[ii]
  if(DATA_past1000$CAVES$entity_info$mineralogy[DATA_past1000$CAVES$entity_info$entity_id == entity] == "calcite") {mask_mean_calcite[ii] = T}
  else{mask_mean_aragonite[ii] = T}
}


##down-sampled
for(plot in 1:1){
  pdf(file = "Plots/Discussion/Variance_Scatterplot.pdf", width = 8, height = 10)
  #png(file = "Plots/Appendix/A1_ScatterMean_diff_full.png", width = 50*8, height = 50*10)
  par(mfrow=c(5,2),oma = c(1,3,0,0) + 0.1,mar = c(3,0,1,1) + 0.1)
  #latitude d18Oc
  plot(DISCUSSION$scatter_data_variance$latitude[mask_mean_calcite], DISCUSSION$scatter_data_variance$vr_ds_isot[mask_mean_calcite],
       xlab = "", ylab = "", ylim = c(0.01,100), panel.first = grid(), log = "y",
       pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
  points(DISCUSSION$scatter_data_variance$latitude[mask_mean_aragonite], DISCUSSION$scatter_data_variance$vr_ds_isot[mask_mean_aragonite], 
         pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
  abline(h=1)
  lines(lowess(DISCUSSION$scatter_data_variance$latitude, DISCUSSION$scatter_data_variance$vr_ds_isot, f = 2/3, delta = 0.01*diff(range(DISCUSSION$scatter_data_variance$latitude, na.rm = T))), lwd = 4, col = "#B2182B")
  mtext(text = "latitude",side = 1,line = 2)
  mtext(text = "VR d18O ds",side = 2,line = 2)
  text(10, 10, "calcite", cex = 1.5)
  text(10, 0.05, "aragonite", col = "blue", cex = 1.5)
  mtext(text = "(a)", side = 3, line = -1.5, adj = 0.02, cex = 1)
  
  #temperature d18Oc
  plot(DISCUSSION$scatter_data_variance$mean_temp[mask_mean_calcite], DISCUSSION$scatter_data_variance$vr_ds_isot[mask_mean_calcite],
       yaxt = 'n', xlab = "", ylab = "", ylim = c(0.01,100), panel.first = grid(), log = "y",
       pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
  points(DISCUSSION$scatter_data_variance$mean_temp[mask_mean_aragonite], DISCUSSION$scatter_data_variance$vr_ds_isot[mask_mean_aragonite], 
         pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
  abline(h=1)
  lines(lowess(DISCUSSION$scatter_data_variance$mean_temp, DISCUSSION$scatter_data_variance$vr_ds_isot, f = 2/3, delta = 0.01*diff(range(DISCUSSION$scatter_data_variance$mean_temp, na.rm = T))), lwd = 4, col = "#B2182B")
  mtext(text = "mean temp",side = 1,line = 2)
  #mtext(text = "d18Oc",side = 2,line = 2)
  #text(10, 3, "calcite", cex = 1.5)
  #text(10, -13., "aragonite", col = "blue", cex = 1.5)
  mtext(text = "(b)", side = 3, line = -1.5, adj = 0.02, cex = 1)
  
  #latitude
  plot(DISCUSSION$scatter_data_variance$diff_down[mask_mean_calcite], DISCUSSION$scatter_data_variance$vr_ds_isot[mask_mean_calcite], 
       xlab = "", ylab = "",  ylim = c(0.01,100), panel.first = grid(), log = "y",
       pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
  points(DISCUSSION$scatter_data_variance$diff_down[mask_mean_aragonite], DISCUSSION$scatter_data_variance$vr_ds_isot[mask_mean_aragonite], 
         pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
  abline(h=1)
  lines(lowess(DISCUSSION$scatter_data_variance$diff_down, DISCUSSION$scatter_data_variance$vr_ds_isot, f = 2/3, delta = 0.01*diff(range(DISCUSSION$scatter_data_variance$diff_down, na.rm = T))), lwd = 4, col = "#B2182B")
  mtext(text = "d18O-d18Oc",side = 1,line = 2)
  mtext(text = "VR d18O ds",side = 2,line = 2)
  #text(10, 6.7, "calcite", cex = 1.5)
  #text(-5, -7., "aragonite", col = "blue", cex = 1.5)
  mtext(text = "(c)", side = 3, line = -1.5, adj = 0.02, cex = 1)
  
  #temperature
  plot(DISCUSSION$scatter_data_variance$d18Oc[mask_mean_calcite], DISCUSSION$scatter_data_variance$vr_ds_isot[mask_mean_calcite], 
       xlab = "", ylab = "",  ylim = c(0.01,100), panel.first = grid(), log = "y",
       pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
  points(DISCUSSION$scatter_data_variance$d18Oc[mask_mean_aragonite], DISCUSSION$scatter_data_variance$vr_ds_isot[mask_mean_aragonite], 
         pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
  abline(h=1)
  lines(lowess(DISCUSSION$scatter_data_variance$d18Oc, DISCUSSION$scatter_data_variance$vr_ds_isot, f = 2/3, delta = 0.01*diff(range(DISCUSSION$scatter_data_variance$d18Oc, na.rm = T))), lwd = 4, col = "#B2182B")
  mtext(text = "d18Oc",side = 1,line = 2)
  #text(10, 6.7, "calcite", cex = 1.5)
  #text(-5, -7., "aragonite", col = "blue", cex = 1.5)
  mtext(text = "(d)", side = 3, line = -1.5, adj = 0.02, cex = 1)
  
  #prec
  plot(DISCUSSION$scatter_data_variance$mean_prec[mask_mean_calcite]*8.6148e4, DISCUSSION$scatter_data_variance$vr_ds_isot[mask_mean_calcite], 
       xlab = "", ylab = "", log = "xy", xlim = c(0.2,15), ylim = c(0.01,100), panel.first = grid(equilogs = FALSE),
       pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
  abline(h=1)
  points(DISCUSSION$scatter_data_variance$mean_prec[mask_mean_aragonite]*8.6148e4, DISCUSSION$scatter_data_variance$vr_ds_isot[mask_mean_aragonite], 
         pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
  lines(lowess(DISCUSSION$scatter_data_variance$mean_prec*8.6148e4, DISCUSSION$scatter_data_variance$vr_ds_isot, f = 2/3, delta = 0.01*diff(range(DISCUSSION$scatter_data_variance$winter_mean_prec, na.rm = T))), lwd = 4, col = "#B2182B")
  mtext(text = "mean prec",side = 1,line = 2)
  mtext(text = "VR d18O ds",side = 2,line = 2)
  mtext(text = "(e)", side = 3, line = -1.5, adj = 0.02, cex = 1)
  
  #DJF prec
  plot(DISCUSSION$scatter_data_variance$winter_mean_prec[mask_mean_calcite]*8.6148e4, DISCUSSION$scatter_data_variance$vr_ds_isot[mask_mean_calcite], 
       yaxt = 'n', xlab = "", ylab = "", xlim = c(0.4,15), log = "xy", ylim = c(0.01,100), yaxt = "n", panel.first = grid(equilogs = FALSE),
       pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
  abline(h=1)
  points(DISCUSSION$scatter_data_variance$winter_mean_prec[mask_mean_aragonite]*8.6148e4, DISCUSSION$scatter_data_variance$vr_ds_isot[mask_mean_aragonite],
         pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
  lines(lowess(DISCUSSION$scatter_data_variance$winter_mean_prec*8.6148e4, DISCUSSION$scatter_data_variance$vr_ds_isot, f = 2/3, delta = 0.01*diff(range(DISCUSSION$scatter_data_variance$winter_mean_prec, na.rm = T))), lwd = 4, col = "#B2182B")
  mtext(text = "DJF prec",side = 1,line = 2)
  mtext(text = "(f)", side = 3, line = -1.5, adj = 0.02, cex = 1)
  
  # elevation
  plot(DISCUSSION$scatter_data_variance$elevation[mask_mean_calcite], DISCUSSION$scatter_data_variance$vr_ds_isot[mask_mean_calcite], 
       xlab = "", ylab = "", ylim = c(0.01,100), xlim = c(0,4000), panel.first = grid(), log = "y",
       pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
  abline(h=1)
  points(DISCUSSION$scatter_data_variance$elevation[mask_mean_aragonite], DISCUSSION$scatter_data_variance$vr_ds_isot[mask_mean_aragonite], 
         pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
  lines(lowess(na.omit(as.numeric(DISCUSSION$scatter_data_variance$elevation)), DISCUSSION$scatter_data_variance$vr_ds_isot[!is.na(as.numeric(DISCUSSION$scatter_data_variance$elevation))], f = 2/3, delta = 0.01*diff(range(DISCUSSION$scatter_data_variance$elevation, na.rm = T))), lwd = 4, col = "#B2182B")
  mtext(text = "elevation",side = 1,line = 2)
  mtext(text = "VR d18O ds",side = 2,line = 2)
  mtext(text = "(g)", side = 3, line = -1.5, adj = 0.02, cex = 1)
  
  plot(DISCUSSION$scatter_data_variance$elevation_diff[mask_mean_calcite], DISCUSSION$scatter_data_variance$vr_ds_isot[mask_mean_calcite], 
       yaxt = 'n', xlab = "", ylab = "", ylim = c(0.01,100), xlim = c(-1500,1500), panel.first = grid(), log = "y",
       pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
  abline(h=1)
  points(DISCUSSION$scatter_data_variance$elevation_diff[mask_mean_aragonite], DISCUSSION$scatter_data_variance$vr_ds_isot[mask_mean_aragonite], 
         pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
  lines(lowess(na.omit(DISCUSSION$scatter_data_variance$elevation_diff), DISCUSSION$scatter_data_variance$vr_ds_isot[!is.na(DISCUSSION$scatter_data_variance$diff_full)], f = 2/5, delta = 0.01*diff(range(DISCUSSION$scatter_data_variance$elevation, na.rm = T))), lwd = 4, col = "#B2182B")
  mtext(text = "elevation difference (Sim-Rec)",side = 1,line = 2)
  mtext(text = "(h)", side = 3, line = -1.5, adj = 0.02, cex = 1)
  
  #cover thickness
  plot(DISCUSSION$scatter_data_variance$cover_thickness[mask_mean_calcite], DISCUSSION$scatter_data_variance$vr_ds_isot[mask_mean_calcite], 
       xlab = "", ylab = "", ylim = c(0.01,100), panel.first = grid(), log = "y",
       pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
  abline(h=1)
  points(DISCUSSION$scatter_data_variance$cover_thickness[mask_mean_aragonite], DISCUSSION$scatter_data_variance$vr_ds_isot[mask_mean_aragonite], 
         pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
  lines(lowess(na.omit(DISCUSSION$scatter_data_variance$cover_thickness), DISCUSSION$scatter_data_variance$vr_ds_isot[!is.na(DISCUSSION$scatter_data_variance$cover_thickness)], f = 2/3, delta = 0.01*diff(range(DISCUSSION$scatter_data_variance$mean_temp, na.rm = T))), lwd = 4, col = "#B2182B")
  mtext(text = "cover thickness",side = 1,line = 2)
  mtext(text = "VR d18O ds",side = 2,line = 2)
  mtext(text = "(i)", side = 3, line = -1.5, adj = 0.02, cex = 1)
  
  plot(c(0,4), c(0.01,100), type = "n", xlab = "", ylab = "", yaxt = 'n', xaxt = 'n', log = "y")
  abline(h=1)
  boxplot(as.numeric(DISCUSSION$scatter_data_variance$vr_ds_isot[DISCUSSION$scatter_data_variance$geology == "limestone"]), add = T, at = 0.5, yaxt ="n", axes = F)
  boxplot(as.numeric(DISCUSSION$scatter_data_variance$vr_ds_isot[DISCUSSION$scatter_data_variance$geology == "dolomite"]), add = T, at = 1.5, yaxt ="n", axes = F)
  boxplot(as.numeric(DISCUSSION$scatter_data_variance$vr_ds_isot[DISCUSSION$scatter_data_variance$geology == "marble"]), add = T, at = 2.5, yaxt ="n", axes = F)
  boxplot(as.numeric(DISCUSSION$scatter_data_variance$vr_ds_isot[DISCUSSION$scatter_data_variance$geology == "unknown"]), add = T, at = 3.5, yaxt ="n", axes = F)
  mtext(text = "geology", side = 1, line = 2)
  axis(1,at=c(0.5,1.5,2.5,3.5),labels=c("limestone", "dolomite", "marble", "unknown"))
  mtext(text = "(j)", side = 3, line = -1.5, adj = 0.02, cex = 1)
  
  
  dev.off()
}


##variance
for(plot in 1:1){
  pdf(file = "Plots/Discussion/Variance_Scatterplot_RecordVar.pdf", width = 8, height = 10)
  #png(file = "Plots/Appendix/A1_ScatterMean_diff_full.png", width = 50*8, height = 50*10)
  par(mfrow=c(5,2),oma = c(1,3,0,0) + 0.1,mar = c(3,0,1,1) + 0.1)
  #latitude d18Oc
  plot(DISCUSSION$scatter_data_variance$latitude[mask_mean_calcite], DISCUSSION$scatter_data_variance$var_record[mask_mean_calcite],
       xlab = "", ylab = "", ylim = c(0.005,2.5), panel.first = grid(), log = "y",
       pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
  points(DISCUSSION$scatter_data_variance$latitude[mask_mean_aragonite], DISCUSSION$scatter_data_variance$var_record[mask_mean_aragonite], 
         pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
  abline(h=1)
  lines(lowess(DISCUSSION$scatter_data_variance$latitude, DISCUSSION$scatter_data_variance$var_record, f = 2/3, delta = 0.01*diff(range(DISCUSSION$scatter_data_variance$latitude, na.rm = T))), lwd = 4, col = "#B2182B")
  mtext(text = "latitude",side = 1,line = 2)
  mtext(text = "var d18O",side = 2,line = 2)
  text(10, 10, "calcite", cex = 1.5)
  text(10, 0.05, "aragonite", col = "blue", cex = 1.5)
  mtext(text = "(a)", side = 3, line = -1.5, adj = 0.02, cex = 1)
  
  #temperature d18Oc
  plot(DISCUSSION$scatter_data_variance$mean_temp[mask_mean_calcite], DISCUSSION$scatter_data_variance$var_record[mask_mean_calcite],
       yaxt = 'n', xlab = "", ylab = "", ylim = c(0.005,2.5), panel.first = grid(), log = "y",
       pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
  points(DISCUSSION$scatter_data_variance$mean_temp[mask_mean_aragonite], DISCUSSION$scatter_data_variance$var_record[mask_mean_aragonite], 
         pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
  abline(h=1)
  lines(lowess(DISCUSSION$scatter_data_variance$mean_temp, DISCUSSION$scatter_data_variance$var_record, f = 2/3, delta = 0.01*diff(range(DISCUSSION$scatter_data_variance$mean_temp, na.rm = T))), lwd = 4, col = "#B2182B")
  mtext(text = "mean temp",side = 1,line = 2)
  #mtext(text = "d18Oc",side = 2,line = 2)
  #text(10, 3, "calcite", cex = 1.5)
  #text(10, -13., "aragonite", col = "blue", cex = 1.5)
  mtext(text = "(b)", side = 3, line = -1.5, adj = 0.02, cex = 1)
  
  #latitude
  plot(DISCUSSION$scatter_data_variance$diff_down[mask_mean_calcite], DISCUSSION$scatter_data_variance$var_record[mask_mean_calcite], 
       xlab = "", ylab = "",  ylim = c(0.005,2.5), panel.first = grid(), log = "y",
       pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
  points(DISCUSSION$scatter_data_variance$diff_down[mask_mean_aragonite], DISCUSSION$scatter_data_variance$var_record[mask_mean_aragonite], 
         pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
  abline(h=1)
  lines(lowess(DISCUSSION$scatter_data_variance$diff_down, DISCUSSION$scatter_data_variance$var_record, f = 2/3, delta = 0.01*diff(range(DISCUSSION$scatter_data_variance$diff_down, na.rm = T))), lwd = 4, col = "#B2182B")
  mtext(text = "d18O-d18Oc",side = 1,line = 2)
  mtext(text = "var d18O",side = 2,line = 2)
  #text(10, 6.7, "calcite", cex = 1.5)
  #text(-5, -7., "aragonite", col = "blue", cex = 1.5)
  mtext(text = "(c)", side = 3, line = -1.5, adj = 0.02, cex = 1)
  
  #temperature
  plot(DISCUSSION$scatter_data_variance$d18Oc[mask_mean_calcite], DISCUSSION$scatter_data_variance$var_record[mask_mean_calcite], 
       xlab = "", ylab = "",  ylim = c(0.005,2.5), panel.first = grid(), log = "y", yaxt = "n",
       pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
  points(DISCUSSION$scatter_data_variance$d18Oc[mask_mean_aragonite], DISCUSSION$scatter_data_variance$var_record[mask_mean_aragonite], 
         pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
  abline(h=1)
  lines(lowess(DISCUSSION$scatter_data_variance$d18Oc, DISCUSSION$scatter_data_variance$var_record, f = 2/3, delta = 0.01*diff(range(DISCUSSION$scatter_data_variance$d18Oc, na.rm = T))), lwd = 4, col = "#B2182B")
  mtext(text = "d18Oc",side = 1,line = 2)
  #text(10, 6.7, "calcite", cex = 1.5)
  #text(-5, -7., "aragonite", col = "blue", cex = 1.5)
  mtext(text = "(d)", side = 3, line = -1.5, adj = 0.02, cex = 1)
  
  #prec
  plot(DISCUSSION$scatter_data_variance$mean_prec[mask_mean_calcite]*8.6148e4, DISCUSSION$scatter_data_variance$var_record[mask_mean_calcite], 
       xlab = "", ylab = "", log = "xy", xlim = c(0.2,15), ylim = c(0.005,2.5), panel.first = grid(equilogs = FALSE),
       pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
  abline(h=1)
  points(DISCUSSION$scatter_data_variance$mean_prec[mask_mean_aragonite]*8.6148e4, DISCUSSION$scatter_data_variance$var_record[mask_mean_aragonite], 
         pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
  lines(lowess(DISCUSSION$scatter_data_variance$mean_prec*8.6148e4, DISCUSSION$scatter_data_variance$var_record, f = 2/3, delta = 0.01*diff(range(DISCUSSION$scatter_data_variance$winter_mean_prec, na.rm = T))), lwd = 4, col = "#B2182B")
  mtext(text = "mean prec",side = 1,line = 2)
  mtext(text = "var d18O",side = 2,line = 2)  
  mtext(text = "(e)", side = 3, line = -1.5, adj = 0.02, cex = 1)
  
  #DJF prec
  plot(DISCUSSION$scatter_data_variance$winter_mean_prec[mask_mean_calcite]*8.6148e4, DISCUSSION$scatter_data_variance$var_record[mask_mean_calcite], 
       yaxt = 'n', xlab = "", ylab = "", xlim = c(0.4,15), log = "xy", ylim = c(0.005,2.5), yaxt = "n", panel.first = grid(equilogs = FALSE),
       pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
  abline(h=1)
  points(DISCUSSION$scatter_data_variance$winter_mean_prec[mask_mean_aragonite]*8.6148e4, DISCUSSION$scatter_data_variance$var_record[mask_mean_aragonite],
         pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
  lines(lowess(DISCUSSION$scatter_data_variance$winter_mean_prec*8.6148e4, DISCUSSION$scatter_data_variance$var_record, f = 2/3, delta = 0.01*diff(range(DISCUSSION$scatter_data_variance$winter_mean_prec, na.rm = T))), lwd = 4, col = "#B2182B")
  mtext(text = "DJF prec",side = 1,line = 2)
  mtext(text = "(f)", side = 3, line = -1.5, adj = 0.02, cex = 1)
  
  # elevation
  plot(DISCUSSION$scatter_data_variance$elevation[mask_mean_calcite], DISCUSSION$scatter_data_variance$var_record[mask_mean_calcite], 
       xlab = "", ylab = "", ylim = c(0.005,2.5), xlim = c(0,4000), panel.first = grid(), log = "y",
       pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
  abline(h=1)
  points(DISCUSSION$scatter_data_variance$elevation[mask_mean_aragonite], DISCUSSION$scatter_data_variance$var_record[mask_mean_aragonite], 
         pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
  lines(lowess(na.omit(as.numeric(DISCUSSION$scatter_data_variance$elevation)), DISCUSSION$scatter_data_variance$var_record[!is.na(as.numeric(DISCUSSION$scatter_data_variance$elevation))], f = 2/3, delta = 0.01*diff(range(DISCUSSION$scatter_data_variance$elevation, na.rm = T))), lwd = 4, col = "#B2182B")
  mtext(text = "elevation",side = 1,line = 2)
  mtext(text = "var d18O",side = 2,line = 2)  
  mtext(text = "(g)", side = 3, line = -1.5, adj = 0.02, cex = 1)
  
  plot(DISCUSSION$scatter_data_variance$elevation_diff[mask_mean_calcite], DISCUSSION$scatter_data_variance$var_record[mask_mean_calcite], 
       yaxt = 'n', xlab = "", ylab = "", ylim = c(0.005,2.5), xlim = c(-1500,1500), panel.first = grid(), log = "y",
       pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
  abline(h=1)
  points(DISCUSSION$scatter_data_variance$elevation_diff[mask_mean_aragonite], DISCUSSION$scatter_data_variance$var_record[mask_mean_aragonite], 
         pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
  lines(lowess(na.omit(DISCUSSION$scatter_data_variance$elevation_diff), DISCUSSION$scatter_data_variance$var_record[!is.na(DISCUSSION$scatter_data_variance$diff_full)], f = 2/5, delta = 0.01*diff(range(DISCUSSION$scatter_data_variance$elevation, na.rm = T))), lwd = 4, col = "#B2182B")
  mtext(text = "elevation difference (Sim-Rec)",side = 1,line = 2)
  mtext(text = "(h)", side = 3, line = -1.5, adj = 0.02, cex = 1)
  
  #cover thickness
  plot(DISCUSSION$scatter_data_variance$cover_thickness[mask_mean_calcite], DISCUSSION$scatter_data_variance$var_record[mask_mean_calcite], 
       xlab = "", ylab = "", ylim = c(0.005,2.5), panel.first = grid(), log = "y",
       pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
  abline(h=1)
  points(DISCUSSION$scatter_data_variance$cover_thickness[mask_mean_aragonite], DISCUSSION$scatter_data_variance$var_record[mask_mean_aragonite], 
         pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
  lines(lowess(na.omit(DISCUSSION$scatter_data_variance$cover_thickness), DISCUSSION$scatter_data_variance$var_record[!is.na(DISCUSSION$scatter_data_variance$cover_thickness)], f = 2/3, delta = 0.01*diff(range(DISCUSSION$scatter_data_variance$mean_temp, na.rm = T))), lwd = 4, col = "#B2182B")
  mtext(text = "cover thickness",side = 1,line = 2)
  mtext(text = "var d18O",side = 2,line = 2)
  mtext(text = "(i)", side = 3, line = -1.5, adj = 0.02, cex = 1)
  
  plot(c(0,4), c(0.01,100), type = "n", xlab = "", ylab = "", yaxt = 'n', xaxt = 'n', log = "y")
  abline(h=1)
  boxplot(as.numeric(DISCUSSION$scatter_data_variance$var_record[DISCUSSION$scatter_data_variance$geology == "limestone"]), add = T, at = 0.5, yaxt ="n", axes = F)
  boxplot(as.numeric(DISCUSSION$scatter_data_variance$var_record[DISCUSSION$scatter_data_variance$geology == "dolomite"]), add = T, at = 1.5, yaxt ="n", axes = F)
  boxplot(as.numeric(DISCUSSION$scatter_data_variance$var_record[DISCUSSION$scatter_data_variance$geology == "marble"]), add = T, at = 2.5, yaxt ="n", axes = F)
  boxplot(as.numeric(DISCUSSION$scatter_data_variance$var_record[DISCUSSION$scatter_data_variance$geology == "unknown"]), add = T, at = 3.5, yaxt ="n", axes = F)
  mtext(text = "geology", side = 1, line = 2)
  axis(1,at=c(0.5,1.5,2.5,3.5),labels=c("limestone", "dolomite", "marble", "unknown"))
  mtext(text = "(j)", side = 3, line = -1.5, adj = 0.02, cex = 1)
  
  
  dev.off()
}
