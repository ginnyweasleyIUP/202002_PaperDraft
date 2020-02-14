#################################################
## MEAN SCATTER PLOTS ###########################
#################################################

# 2 FRAGEN:
#           1) downsampeln oder full res?
#           2) 

scatter_data = array(dim = c(length(DATA_past1000$CAVES$entity_info$entity_id[mask_mean]),17))
colnames(scatter_data) = c("entity_id", "site_id", "diff_full_a", "diff_down_a","diff_full_b", "diff_down_b","diff_full_c", "diff_down_c",
                           "elevation", "latitude", "mean_temp", "mean_prec", "winter_mean_prec", "elevation_diff", 
                           "dist_entrance", "geology", "cover_thickness")

scatter_data[,1] = DATA_past1000$CAVES$entity_info$entity_id[mask_mean]
for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id[mask_mean])){
  entity = scatter_data[[ii,1]]
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  index = match(site, DATA_past1000$CAVES$site_info$site_id)
  scatter_data[ii,2] = site
  scatter_data[ii,3] = mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ISOT_a, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a, na.rm = T)
  scatter_data[ii,4] = mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$ISOT_a, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a, na.rm = T)
  scatter_data[ii,5] = mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ISOT_b, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_b, na.rm = T)
  scatter_data[ii,6] = mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$ISOT_b, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_b, na.rm = T)
  scatter_data[ii,7] = mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ISOT_c, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_c, na.rm = T)
  scatter_data[ii,8] = mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$ISOT_c, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_c, na.rm = T)
  
 
  scatter_data[ii,9] = DATA_past1000$CAVES$site_info$elevation[index]
  scatter_data[ii,10] = DATA_past1000$CAVES$site_info$latitude[index]
  scatter_data[ii,11] = mean(c(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$TEMP_a,
                               DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$TEMP_b,
                               DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$TEMP_c), na.rm = T)
  scatter_data[ii,12] = mean(c(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$PREC_a,
                               DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$PREC_b,
                               DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$PREC_c), na.rm = T)
  scatter_data[ii,13] = mean(c(DATA_past1000$CAVES$sim_data_seasonal$a[[paste0("CAVE", site)]]$WINTER$prec_mean,
                               DATA_past1000$CAVES$sim_data_seasonal$b[[paste0("CAVE", site)]]$WINTER$prec_mean,
                               DATA_past1000$CAVES$sim_data_seasonal$c[[paste0("CAVE", site)]]$WINTER$prec_mean), na.rm = T)
  scatter_data[ii,14] = as.numeric(DATA_past1000$CAVES$elevation_cave_sim$`sim-cave`[DATA_past1000$CAVES$elevation_cave_sim$entity_id == entity])
  scatter_data[ii,15] = as.numeric(DATA_past1000$CAVES$entity_info$distance_entrance[match(entity, DATA_past1000$CAVES$entity_info$entity_id)])
  scatter_data[ii,16] = as.numeric(DATA_past1000$CAVES$entity_info$geology[match(entity, DATA_past1000$CAVES$entity_info$entity_id)])
  scatter_data[ii,17] = NA
}

scatter_data = as.tibble(scatter_data)

#View(scatter_data)
# mask for aragonite and calcite


mask_mean_calcite = logical(length = length(MEAN_ANALYSIS$MEAN_DIFF_ISOT$entity_id[mask_mean]))
mask_mean_aragonite  = logical(length = length(MEAN_ANALYSIS$MEAN_DIFF_ISOT$entity_id[mask_mean]))

for(ii in 1:length(MEAN_ANALYSIS$MEAN_DIFF_ISOT$entity_id[mask_mean])){
  entity = MEAN_ANALYSIS$MEAN_DIFF_ISOT$entity_id[mask_mean][ii]
  if(DATA_past1000$CAVES$entity_info$mineralogy[DATA_past1000$CAVES$entity_info$entity_id == entity] == "calcite") {mask_mean_calcite[ii] = T}
  else{mask_mean_aragonite[ii] = T}
}


##diff_full

pdf(file = "Plots/Paper/Appendix_01_ScatterMean_diff-full.pdf", width = 8, height = 8)
par(mfrow=c(3,2),oma = c(1,3,0,0) + 0.1,mar = c(3,0,1,1) + 0.1)
# elevation
plot(scatter_data$elevation[mask_mean_calcite], scatter_data$diff_full[mask_mean_calcite], xlab = "", ylab = "", ylim = c(-10,10), xlim = c(0,4000), panel.first = grid())
abline(h=0)
points(scatter_data$elevation[mask_mean_aragonite], scatter_data$diff_full[mask_mean_aragonite], pch = 8, col = "blue")
lines(lowess(scatter_data$elevation, scatter_data$diff_full, f = 2/3, delta = 0.01*diff(range(scatter_data$elevation, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "elevation",side = 1,line = 2)
mtext(text = "d18O-d18Oc",side = 2,line = 2)
text(350, 5.7, "calcite")
text(1150, -7., "aragonite", col = "blue")

plot(scatter_data$elevation_diff[mask_mean_calcite], scatter_data$diff_full[mask_mean_calcite], yaxt = 'n', xlab = "", ylab = "", ylim = c(-10,10), panel.first = grid())
abline(h=0)
points(scatter_data$elevation_diff[mask_mean_aragonite], scatter_data$diff_full[mask_mean_aragonite], pch = 8, col = "blue")
lines(lowess(scatter_data$elevation_diff, scatter_data$diff_full, f = 2/5, delta = 0.01*diff(range(scatter_data$elevation, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "elevation diff. (sim-rec)",side = 1,line = 2)

plot(scatter_data$latitude[mask_mean_calcite], scatter_data$diff_full[mask_mean_calcite], xlab = "", ylab = "", ylim = c(-10,10), panel.first = grid())
points(scatter_data$latitude[mask_mean_aragonite], scatter_data$diff_full[mask_mean_aragonite], pch = 8, col = "blue")
abline(h=0)
lines(lowess(scatter_data$latitude, scatter_data$diff_full, f = 2/3, delta = 0.01*diff(range(scatter_data$latitude, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "latitude",side = 1,line = 2)
mtext(text = "d18O-d18Oc",side = 2,line = 2)

plot(scatter_data$mean_temp[mask_mean_calcite], scatter_data$diff_full[mask_mean_calcite], yaxt = 'n', xlab = "", ylab = "", ylim = c(-10,10), panel.first = grid())
abline(h=0)
points(scatter_data$mean_temp[mask_mean_aragonite], scatter_data$diff_full[mask_mean_aragonite], pch = 8, col = "blue")
lines(lowess(scatter_data$mean_temp, scatter_data$diff_full, f = 2/3, delta = 0.01*diff(range(scatter_data$mean_temp, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "mean temp",side = 1,line = 2)

plot(scatter_data$mean_prec[mask_mean_calcite]*8.6148e4, scatter_data$diff_full[mask_mean_calcite], xlab = "", ylab = "", log = "x", ylim = c(-10,10), xlim = c(0.2,15), panel.first = grid(equilogs = FALSE))
abline(h=0)
points(scatter_data$mean_prec[mask_mean_aragonite]*8.6148e4, scatter_data$diff_full[mask_mean_aragonite], pch = 8, col = "blue")
lines(lowess(scatter_data$mean_prec*8.6148e4, scatter_data$diff_full, f = 2/3, delta = 0.01*diff(range(scatter_data$winter_mean_prec, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "mean prec",side = 1,line = 2)
mtext(text = "d18O-d18Oc",side = 2,line = 2)

plot(scatter_data$winter_mean_prec[mask_mean_calcite]*8.6148e4, scatter_data$diff_full[mask_mean_calcite], yaxt = 'n', xlab = "", ylab = "", xlim = c(0.2,15), log = "x", ylim = c(-10,10), yaxt = "n", panel.first = grid(equilogs = FALSE))
abline(h=0)
points(scatter_data$winter_mean_prec[mask_mean_aragonite]*8.6148e4, scatter_data$diff_full[mask_mean_aragonite], pch = 8, col = "blue")
lines(lowess(scatter_data$winter_mean_prec*8.6148e4, scatter_data$diff_full, f = 2/3, delta = 0.01*diff(range(scatter_data$winter_mean_prec, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "winter prec",side = 1,line = 2)



dev.off()

##down_sampled

pdf(file = "Plots/Paper/Appendix_01_ScatterMean_diff-down.pdf", width = 8, height = PLOTTING_VARIABLES$HEIGHT/2)
par(mfrow=c(2,2),oma = c(1,3,0,0) + 0.1,mar = c(3,0,1,1) + 0.1)
plot(scatter_data$elevation[mask_mean_calcite], scatter_data$diff_down[mask_mean_calcite], xlab = "", ylab = "", ylim = c(-10,10), xlim = c(0,4000), panel.first = grid())
points(scatter_data$elevation[mask_mean_aragonite], scatter_data$diff_down[mask_mean_aragonite], pch = 8, col = "blue")
lines(lowess(scatter_data$elevation, scatter_data$diff_down, f = 2/3, delta = 0.01*diff(range(scatter_data$elevation, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "elevation",side = 1,line = 2)
mtext(text = "diff-d18O",side = 2,line = 2)
text(350, 5.7, "calcite")
text(1150, -7., "aragonite", col = "blue")
plot(scatter_data$latitude[mask_mean_calcite], scatter_data$diff_down[mask_mean_calcite], xlab = "", ylab = "", ylim = c(-10,10), panel.first = grid(), yaxt ="n")
points(scatter_data$latitude[mask_mean_aragonite], scatter_data$diff_down[mask_mean_aragonite], pch = 8, col = "blue")
lines(lowess(scatter_data$latitude, scatter_data$diff_down, f = 2/3, delta = 0.01*diff(range(scatter_data$latitude, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "latitude",side = 1,line = 2)

plot(scatter_data$mean_temp[mask_mean_calcite], scatter_data$diff_down[mask_mean_calcite], xlab = "", ylab = "", ylim = c(-10,10), panel.first = grid())
points(scatter_data$mean_temp[mask_mean_aragonite], scatter_data$diff_down[mask_mean_aragonite], pch = 8, col = "blue")
lines(lowess(scatter_data$mean_temp, scatter_data$diff_down, f = 2/3, delta = 0.01*diff(range(scatter_data$mean_temp, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "mean temp",side = 1,line = 2)
mtext(text = "diff-d18O",side = 2,line = 2)
plot(scatter_data$winter_mean_prec[mask_mean_calcite], scatter_data$diff_down[mask_mean_calcite], xlab = "", ylab = "", log = "x", ylim = c(-10,10), yaxt = "n", panel.first = grid(equilogs = FALSE))
points(scatter_data$winter_mean_prec[mask_mean_aragonite], scatter_data$diff_down[mask_mean_aragonite], pch = 8, col = "blue")
lines(lowess(scatter_data$winter_mean_prec, scatter_data$diff_down, f = 2/3, delta = 0.01*diff(range(scatter_data$winter_mean_prec, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "winter prec",side = 1,line = 2)
dev.off()

