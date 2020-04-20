#################################################
## MEAN SCATTER PLOTS ###########################
#################################################

# 2 FRAGEN:
#           1) downsampeln oder full res?
#           2) 

scatter_data = array(dim = c(length(DATA_past1000$CAVES$entity_info$entity_id[mask_mean]),20))
colnames(scatter_data) = c("entity_id", "site_id", "diff_full", "diff_down",
                           "diff_full_a", "diff_down_a","diff_full_b", "diff_down_b","diff_full_c", "diff_down_c",
                           "elevation", "latitude", "mean_temp", "mean_prec", "winter_mean_prec", "elevation_diff", 
                           "dist_entrance", "geology", "cover_thickness", "d18Oc")

scatter_data[,1] = DATA_past1000$CAVES$entity_info$entity_id[mask_mean]
for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id[mask_mean])){
  entity = scatter_data[[ii,1]]
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  index = match(site, DATA_past1000$CAVES$site_info$site_id)
  scatter_data[ii,2] = as.numeric(site)
  scatter_data[ii,3] = as.numeric(mean(c(mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ISOT_a, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a, na.rm = T),
                              mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ISOT_b, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_b, na.rm = T),
                              mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ISOT_c, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_c, na.rm = T))))
  scatter_data[ii,4] = as.numeric(mean(c(mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$ISOT_a, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a, na.rm = T),
                              mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$ISOT_b, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_b, na.rm = T),
                              mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$ISOT_c, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_c, na.rm = T))))
  scatter_data[ii,5] = as.numeric(mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ISOT_a, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a, na.rm = T))
  scatter_data[ii,6] = as.numeric(mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$ISOT_a, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a, na.rm = T))
  scatter_data[ii,7] = as.numeric(mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ISOT_b, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_b, na.rm = T))
  scatter_data[ii,8] = as.numeric(mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$ISOT_b, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_b, na.rm = T))
  scatter_data[ii,9] = as.numeric(mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ISOT_c, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_c, na.rm = T))
  scatter_data[ii,10] = as.numeric(mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$ISOT_c, na.rm = T) - mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_c, na.rm = T))
  
 
  scatter_data[ii,11] = as.numeric(DATA_past1000$CAVES$site_info$elevation[index])
  scatter_data[ii,12] = as.numeric(DATA_past1000$CAVES$site_info$latitude[index])
  scatter_data[ii,13] = as.numeric(mean(c(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$TEMP_a,
                                          DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$TEMP_b,
                                          DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$TEMP_c), na.rm = T))
  scatter_data[ii,14] = as.numeric(mean(c(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$PREC_a,
                                          DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$PREC_b,
                                          DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$PREC_c), na.rm = T))
  scatter_data[ii,15] = as.numeric(mean(c(DATA_past1000$CAVES$sim_data_seasonal$a[[paste0("CAVE", site)]]$WINTER$prec_mean,
                                          DATA_past1000$CAVES$sim_data_seasonal$b[[paste0("CAVE", site)]]$WINTER$prec_mean,
                                          DATA_past1000$CAVES$sim_data_seasonal$c[[paste0("CAVE", site)]]$WINTER$prec_mean), na.rm = T))
  scatter_data[ii,16] = as.numeric(DATA_past1000$CAVES$elevation_cave_sim$`sim-cave`[DATA_past1000$CAVES$elevation_cave_sim$entity_id == entity])
  scatter_data[ii,17] = as.numeric(DATA_past1000$CAVES$entity_info$distance_entrance[match(entity, DATA_past1000$CAVES$entity_info$entity_id)])
  scatter_data[ii,18] = DATA_past1000$CAVES$site_to_entity$geology[match(entity, DATA_past1000$CAVES$site_to_entity$entity_id)]
  scatter_data[ii,19] = as.numeric(DATA_past1000$CAVES$site_to_entity$cover_thickness[match(entity, DATA_past1000$CAVES$site_to_entity$entity_id)])
  scatter_data[ii,20] = as.numeric(mean(c(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a,
                                          DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a,
                                          DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq_a), na.rm = T))
}

scatter_data = as.tibble(scatter_data)
scatter_data$entity_id = as.numeric(scatter_data$entity_id)
scatter_data$site_id = as.numeric(scatter_data$site_id)
scatter_data$diff_full = as.numeric(scatter_data$diff_full)
scatter_data$diff_down = as.numeric(scatter_data$diff_down)
scatter_data$elevation = as.numeric(scatter_data$elevation)
scatter_data$latitude = as.numeric(scatter_data$latitude)
scatter_data$mean_temp = as.numeric(scatter_data$mean_temp)
scatter_data$mean_prec = as.numeric(scatter_data$mean_prec)
scatter_data$winter_mean_prec = as.numeric(scatter_data$winter_mean_prec)
scatter_data$elevation_diff = as.numeric(scatter_data$elevation_diff)
scatter_data$dist_entrance = as.numeric(scatter_data$dist_entrance)
scatter_data$cover_thickness = as.numeric(scatter_data$cover_thickness)
scatter_data$d18Oc = as.numeric(scatter_data$d18Oc)

scatter_data$winter_mean_prec[57] = NA

#View(scatter_data)
# mask for aragonite and calcite


mask_mean_calcite = logical(length = length(scatter_data$entity_id))
mask_mean_aragonite  = logical(length = length(scatter_data$entity_id))

for(ii in 1:length(scatter_data$entity_id)){
  entity = scatter_data$entity_id[ii]
  if(DATA_past1000$CAVES$entity_info$mineralogy[DATA_past1000$CAVES$entity_info$entity_id == entity] == "calcite") {mask_mean_calcite[ii] = T}
  else{mask_mean_aragonite[ii] = T}
}


##diff_full

pdf(file = "Plots/Appendix/A1_ScatterMean_diff_full.pdf", width = 8, height = 10)
#png(file = "Plots/Appendix/A1_ScatterMean_diff_full.png", width = 50*8, height = 50*10)
par(mfrow=c(5,2),oma = c(1,3,0,0) + 0.1,mar = c(3,0,1,1) + 0.1)
#latitude d18Oc
plot(scatter_data$latitude[mask_mean_calcite], scatter_data$d18Oc[mask_mean_calcite],
     xlab = "", ylab = "", ylim = c(-20,5), panel.first = grid(),
     pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
points(scatter_data$latitude[mask_mean_aragonite], scatter_data$d18Oc[mask_mean_aragonite], 
       pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
abline(h=0)
lines(lowess(scatter_data$latitude, scatter_data$d18Oc, f = 2/3, delta = 0.01*diff(range(scatter_data$latitude, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "latitude",side = 1,line = 2)
mtext(text = "d18Oc",side = 2,line = 2)
text(10, 3, "calcite", cex = 1.5)
text(10, -13., "aragonite", col = "blue", cex = 1.5)
mtext(text = "(a)", side = 3, line = -1.5, adj = 0.02, cex = 1)

#temperature d18Oc
plot(scatter_data$mean_temp[mask_mean_calcite], scatter_data$d18Oc[mask_mean_calcite],
     yaxt = 'n', xlab = "", ylab = "", ylim = c(-20,5), panel.first = grid(),
     pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
points(scatter_data$mean_temp[mask_mean_aragonite], scatter_data$d18Oc[mask_mean_aragonite], 
       pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
abline(h=0)
lines(lowess(scatter_data$mean_temp, scatter_data$d18Oc, f = 2/3, delta = 0.01*diff(range(scatter_data$mean_temp, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "mean temp",side = 1,line = 2)
#mtext(text = "d18Oc",side = 2,line = 2)
#text(10, 3, "calcite", cex = 1.5)
#text(10, -13., "aragonite", col = "blue", cex = 1.5)
mtext(text = "(b)", side = 3, line = -1.5, adj = 0.02, cex = 1)

#latitude
plot(scatter_data$latitude[mask_mean_calcite], scatter_data$diff_full[mask_mean_calcite], 
     xlab = "", ylab = "", ylim = c(-10,10), panel.first = grid(),
     pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
points(scatter_data$latitude[mask_mean_aragonite], scatter_data$diff_full[mask_mean_aragonite], 
       pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
abline(h=0)
lines(lowess(scatter_data$latitude, scatter_data$diff_full, f = 2/3, delta = 0.01*diff(range(scatter_data$latitude, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "latitude",side = 1,line = 2)
mtext(text = "d18O-d18Oc",side = 2,line = 2)
#text(10, 6.7, "calcite", cex = 1.5)
#text(-5, -7., "aragonite", col = "blue", cex = 1.5)
mtext(text = "(c)", side = 3, line = -1.5, adj = 0.02, cex = 1)

#temperature
plot(scatter_data$mean_temp[mask_mean_calcite], scatter_data$diff_full[mask_mean_calcite], 
     yaxt = 'n', xlab = "", ylab = "", ylim = c(-10,10), panel.first = grid(),
     pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
abline(h=0)
points(scatter_data$mean_temp[mask_mean_aragonite], scatter_data$diff_full[mask_mean_aragonite], 
       pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
lines(lowess(scatter_data$mean_temp, scatter_data$diff_full, f = 2/3, delta = 0.01*diff(range(scatter_data$mean_temp, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "mean temp",side = 1,line = 2)
mtext(text = "(d)", side = 3, line = -1.5, adj = 0.02, cex = 1)

#prec
plot(scatter_data$mean_prec[mask_mean_calcite]*8.6148e4, scatter_data$diff_full[mask_mean_calcite], 
     xlab = "", ylab = "", log = "x", ylim = c(-10,10), xlim = c(0.2,15), panel.first = grid(equilogs = FALSE),
     pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
abline(h=0)
points(scatter_data$mean_prec[mask_mean_aragonite]*8.6148e4, scatter_data$diff_full[mask_mean_aragonite], 
       pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
lines(lowess(scatter_data$mean_prec*8.6148e4, scatter_data$diff_full, f = 2/3, delta = 0.01*diff(range(scatter_data$winter_mean_prec, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "mean prec",side = 1,line = 2)
mtext(text = "d18O-d18Oc",side = 2,line = 2)
mtext(text = "(e)", side = 3, line = -1.5, adj = 0.02, cex = 1)

#DJF prec
plot(scatter_data$winter_mean_prec[mask_mean_calcite]*8.6148e4, scatter_data$diff_full[mask_mean_calcite], 
     yaxt = 'n', xlab = "", ylab = "", xlim = c(0.4,15), log = "x", ylim = c(-10,10), yaxt = "n", panel.first = grid(equilogs = FALSE),
     pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
abline(h=0)
points(scatter_data$winter_mean_prec[mask_mean_aragonite]*8.6148e4, scatter_data$diff_full[mask_mean_aragonite],
       pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
lines(lowess(scatter_data$winter_mean_prec*8.6148e4, scatter_data$diff_full, f = 2/3, delta = 0.01*diff(range(scatter_data$winter_mean_prec, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "DJF prec",side = 1,line = 2)
mtext(text = "(f)", side = 3, line = -1.5, adj = 0.02, cex = 1)

# elevation
plot(scatter_data$elevation[mask_mean_calcite], scatter_data$diff_full[mask_mean_calcite], 
     xlab = "", ylab = "", ylim = c(-10,10), xlim = c(0,4000), panel.first = grid(),
     pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
abline(h=0)
points(scatter_data$elevation[mask_mean_aragonite], scatter_data$diff_full[mask_mean_aragonite], 
       pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
lines(lowess(na.omit(as.numeric(scatter_data$elevation)), scatter_data$diff_full[!is.na(as.numeric(scatter_data$elevation))], f = 2/3, delta = 0.01*diff(range(scatter_data$elevation, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "elevation",side = 1,line = 2)
mtext(text = "d18O-d18Oc",side = 2,line = 2)
mtext(text = "(g)", side = 3, line = -1.5, adj = 0.02, cex = 1)

plot(scatter_data$elevation_diff[mask_mean_calcite], scatter_data$diff_full[mask_mean_calcite], 
     yaxt = 'n', xlab = "", ylab = "", ylim = c(-10,10), panel.first = grid(),
     pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
abline(h=0)
points(scatter_data$elevation_diff[mask_mean_aragonite], scatter_data$diff_full[mask_mean_aragonite], 
       pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
lines(lowess(na.omit(scatter_data$elevation_diff), scatter_data$diff_full[!is.na(scatter_data$diff_full)], f = 2/5, delta = 0.01*diff(range(scatter_data$elevation, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "elevation diff. (sim-rec)",side = 1,line = 2)
mtext(text = "(h)", side = 3, line = -1.5, adj = 0.02, cex = 1)

#cover thickness
plot(scatter_data$cover_thickness[mask_mean_calcite], scatter_data$diff_full[mask_mean_calcite], 
     xlab = "", ylab = "", ylim = c(-10,10), panel.first = grid(),
     pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
abline(h=0)
points(scatter_data$cover_thickness[mask_mean_aragonite], scatter_data$diff_full[mask_mean_aragonite], 
       pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
lines(lowess(na.omit(scatter_data$cover_thickness), scatter_data$diff_full[!is.na(scatter_data$cover_thickness)], f = 2/3, delta = 0.01*diff(range(scatter_data$mean_temp, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "cover thickness",side = 1,line = 2)
mtext(text = "d18O-d18Oc",side = 2,line = 2)
mtext(text = "(i)", side = 3, line = -1.5, adj = 0.02, cex = 1)

plot(c(0,4), c(-10,10), type = "n", xlab = "", ylab = "", yaxt = 'n', xaxt = 'n')
abline(h=0)
boxplot(as.numeric(scatter_data$diff_full[scatter_data$geology == "limestone"]), add = T, at = 0.5, yaxt ="n", axes = F)
boxplot(as.numeric(scatter_data$diff_full[scatter_data$geology == "dolomite"]), add = T, at = 1.5, yaxt ="n", axes = F)
boxplot(as.numeric(scatter_data$diff_full[scatter_data$geology == "marble"]), add = T, at = 2.5, yaxt ="n", axes = F)
boxplot(as.numeric(scatter_data$diff_full[scatter_data$geology == "unknown"]), add = T, at = 3.5, yaxt ="n", axes = F)
mtext(text = "geology", side = 1, line = 2)
axis(1,at=c(0.5,1.5,2.5,3.5),labels=c("limestone", "dolomite", "marble", "unknown"))
mtext(text = "(j)", side = 3, line = -1.5, adj = 0.02, cex = 1)


dev.off()

# ##diff_down
# 
# pdf(file = "Plots/Appendix/A1_ScatterMean_diff_down.pdf", width = 8, height = 10)
# #png(file = "Plots/Appendix/A1_ScatterMean_diff_down.png", width = 50*8, height = 50*10)
# par(mfrow=c(5,2),oma = c(1,3,0,0) + 0.1,mar = c(3,0,1,1) + 0.1)
# 
# #latitude d18Oc
# plot(scatter_data$latitude[mask_mean_calcite], scatter_data$d18Oc[mask_mean_calcite],
#      xlab = "", ylab = "", ylim = c(-20,5), panel.first = grid(),
#      pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
# points(scatter_data$latitude[mask_mean_aragonite], scatter_data$d18Oc[mask_mean_aragonite], 
#        pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
# abline(h=0)
# lines(lowess(scatter_data$latitude, scatter_data$d18Oc, f = 2/3, delta = 0.01*diff(range(scatter_data$latitude, na.rm = T))), lwd = 4, col = "#B2182B")
# mtext(text = "latitude",side = 1,line = 2)
# mtext(text = "d18Oc",side = 2,line = 2)
# text(10, 3, "calcite", cex = 1.5)
# text(10, -13., "aragonite", col = "blue", cex = 1.5)
# mtext(text = "(a)", side = 3, line = -1.5, adj = 0.02, cex = 1)
# 
# #temperature d18Oc
# plot(scatter_data$mean_temp[mask_mean_calcite], scatter_data$d18Oc[mask_mean_calcite],
#      xlab = "", ylab = "", ylim = c(-20,5), panel.first = grid(),
#      pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
# points(scatter_data$mean_temp[mask_mean_aragonite], scatter_data$d18Oc[mask_mean_aragonite], 
#        pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
# abline(h=0)
# lines(lowess(scatter_data$mean_temp, scatter_data$d18Oc, f = 2/3, delta = 0.01*diff(range(scatter_data$mean_temp, na.rm = T))), lwd = 4, col = "#B2182B")
# mtext(text = "mean_temp",side = 1,line = 2)
# mtext(text = "d18Oc",side = 2,line = 2)
# #text(10, 3, "calcite", cex = 1.5)
# #text(10, -13., "aragonite", col = "blue", cex = 1.5)
# mtext(text = "(b)", side = 3, line = -1.5, adj = 0.02, cex = 1)
# 
# #latitude
# plot(scatter_data$latitude[mask_mean_calcite], scatter_data$diff_down[mask_mean_calcite], 
#      xlab = "", ylab = "", ylim = c(-10,10), panel.first = grid(),
#      pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
# points(scatter_data$latitude[mask_mean_aragonite], scatter_data$diff_down[mask_mean_aragonite], 
#        pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
# abline(h=0)
# lines(lowess(scatter_data$latitude, scatter_data$diff_down, f = 2/3, delta = 0.01*diff(range(scatter_data$latitude, na.rm = T))), lwd = 4, col = "#B2182B")
# mtext(text = "latitude",side = 1,line = 2)
# mtext(text = "d18O-d18Oc",side = 2,line = 2)
# text(10, 6.7, "calcite", cex = 1.5)
# text(-5, -7., "aragonite", col = "blue", cex = 1.5)
# 
# plot(scatter_data$mean_temp[mask_mean_calcite], scatter_data$diff_down[mask_mean_calcite], 
#      yaxt = 'n', xlab = "", ylab = "", ylim = c(-10,10), panel.first = grid(),
#      pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
# abline(h=0)
# points(scatter_data$mean_temp[mask_mean_aragonite], scatter_data$diff_down[mask_mean_aragonite], 
#        pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
# lines(lowess(scatter_data$mean_temp, scatter_data$diff_down, f = 2/3, delta = 0.01*diff(range(scatter_data$mean_temp, na.rm = T))), lwd = 4, col = "#B2182B")
# mtext(text = "mean temp",side = 1,line = 2)
# 
# plot(scatter_data$mean_prec[mask_mean_calcite]*8.6148e4, scatter_data$diff_down[mask_mean_calcite], 
#      xlab = "", ylab = "", log = "x", ylim = c(-10,10), xlim = c(0.2,15), panel.first = grid(equilogs = FALSE),
#      pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
# abline(h=0)
# points(scatter_data$mean_prec[mask_mean_aragonite]*8.6148e4, scatter_data$diff_down[mask_mean_aragonite], 
#        pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
# lines(lowess(scatter_data$mean_prec*8.6148e4, scatter_data$diff_down, f = 2/3, delta = 0.01*diff(range(scatter_data$winter_mean_prec, na.rm = T))), lwd = 4, col = "#B2182B")
# mtext(text = "mean prec",side = 1,line = 2)
# mtext(text = "d18O-d18Oc",side = 2,line = 2)
# 
# plot(scatter_data$winter_mean_prec[mask_mean_calcite]*8.6148e4, scatter_data$diff_down[mask_mean_calcite], 
#      yaxt = 'n', xlab = "", ylab = "", xlim = c(0.4,15), log = "x", ylim = c(-10,10), yaxt = "n", panel.first = grid(equilogs = FALSE),
#      pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
# abline(h=0)
# points(scatter_data$winter_mean_prec[mask_mean_aragonite]*8.6148e4, scatter_data$diff_down[mask_mean_aragonite],
#        pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
# lines(lowess(scatter_data$winter_mean_prec*8.6148e4, scatter_data$diff_down, f = 2/3, delta = 0.01*diff(range(scatter_data$winter_mean_prec, na.rm = T))), lwd = 4, col = "#B2182B")
# mtext(text = "DJF prec",side = 1,line = 2)
# # elevation
# plot(scatter_data$elevation[mask_mean_calcite], scatter_data$diff_down[mask_mean_calcite], 
#      xlab = "", ylab = "", ylim = c(-10,10), xlim = c(0,4000), panel.first = grid(),
#      pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
# abline(h=0)
# points(scatter_data$elevation[mask_mean_aragonite], scatter_data$diff_down[mask_mean_aragonite], 
#        pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
# lines(lowess(na.omit(as.numeric(scatter_data$elevation)), scatter_data$diff_down[!is.na(as.numeric(scatter_data$elevation))], f = 2/3, delta = 0.01*diff(range(scatter_data$elevation, na.rm = T))), lwd = 4, col = "#B2182B")
# mtext(text = "elevation",side = 1,line = 2)
# mtext(text = "d18O-d18Oc",side = 2,line = 2)
# 
# plot(scatter_data$elevation_diff[mask_mean_calcite], scatter_data$diff_down[mask_mean_calcite], 
#      yaxt = 'n', xlab = "", ylab = "", ylim = c(-10,10), panel.first = grid(),
#      pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
# abline(h=0)
# points(scatter_data$elevation_diff[mask_mean_aragonite], scatter_data$diff_down[mask_mean_aragonite], 
#        pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
# lines(lowess(na.omit(scatter_data$elevation_diff), scatter_data$diff_down[!is.na(scatter_data$diff_down)], f = 2/5, delta = 0.01*diff(range(scatter_data$elevation, na.rm = T))), lwd = 4, col = "#B2182B")
# mtext(text = "elevation diff. (sim-rec)",side = 1,line = 2)
# plot(scatter_data$cover_thickness[mask_mean_calcite], scatter_data$diff_down[mask_mean_calcite], 
#      xlab = "", ylab = "", ylim = c(-10,10), panel.first = grid(),
#      pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
# abline(h=0)
# points(scatter_data$cover_thickness[mask_mean_aragonite], scatter_data$diff_down[mask_mean_aragonite], 
#        pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2)
# lines(lowess(na.omit(scatter_data$cover_thickness), scatter_data$diff_down[!is.na(scatter_data$cover_thickness)], f = 2/3, delta = 0.01*diff(range(scatter_data$mean_temp, na.rm = T))), lwd = 4, col = "#B2182B")
# mtext(text = "cover thickness",side = 1,line = 2)
# mtext(text = "d18O-d18Oc",side = 2,line = 2)
# 
# plot(c(0,4), c(-10,10), type = "n", xlab = "", ylab = "", yaxt = 'n', xaxt = 'n')
# abline(h=0)
# boxplot(as.numeric(scatter_data$diff_down[scatter_data$geology == "limestone"]), add = T, at = 0.5, yaxt ="n", axes = F)
# boxplot(as.numeric(scatter_data$diff_down[scatter_data$geology == "dolomite"]), add = T, at = 1.5, yaxt ="n", axes = F)
# boxplot(as.numeric(scatter_data$diff_down[scatter_data$geology == "marble"]), add = T, at = 2.5, yaxt ="n", axes = F)
# boxplot(as.numeric(scatter_data$diff_down[scatter_data$geology == "unknown"]), add = T, at = 3.5, yaxt ="n", axes = F)
# mtext(text = "geology", side = 1, line = 2)
# axis(1,at=c(0.5,1.5,2.5,3.5),labels=c("limestone", "dolomite", "marble", "unknown"))
# 
# 
# dev.off()

