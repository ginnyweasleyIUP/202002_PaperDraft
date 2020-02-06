###################################################################################################
## ANALYSIS OF MEAN VALUES --> MEAN BIAS ##########################################################
###################################################################################################
library(plyr)
library(dplyr)
library(latex2exp)
source("Functions/STACYmap_5.R")

# ToDo
# [ ] Mean-Diff over Time:
#      [ ] Cluster Analysis 
#      [ ] aufteilen auf latitudes/temp/prec
#      [ ] prec umrechnen in mm/day
# [ ] Plotting Limit mit Kira festlegen --> was passiert mit NAs? Schwarz machen und per Inkscape auf Skala?
# [ ] Fehlerbalken für Barplot
# [ ]
# [ ] Umrechnung  1 kg/m²/s = 8.6148e4 mm/day

PLOTTING_VARIABLES <- list()

MEAN_ANALYSIS <- list()

ls_mask <- ncdf4::nc_open("landseamask_preindustrial_hadcm3.nc")
plotting_ls_mask <- as.array(ncdf4::ncvar_get(ls_mask))
ncdf4::nc_close(ls_mask)
remove(ls_mask)
plotting_ls_mask[plotting_ls_mask == 0] <- NA

PLOTTING_VARIABLES$ls_mask <- plotting_ls_mask

remove(plotting_ls_mask)

## 1) Create map without any Power Test --> just include all ##
##
## 1.1) Use full simulation over 1000 years



MEAN_ANALYSIS$FIELDS$ISOTlyr <- rbind(DATA_past1000$SIM_mean$ISOT[49:96,1:73]*PLOTTING_VARIABLES$ls_mask[49:96,1:73],
                 DATA_past1000$SIM_mean$ISOT[1:48,1:73]*PLOTTING_VARIABLES$ls_mask[1:48,1:73])

MEAN_ANALYSIS$FIELDS$ITPClyr <- rbind(DATA_past1000$SIM_mean$ITPC[49:96,1:73]*PLOTTING_VARIABLES$ls_mask[49:96,1:73],
                                      DATA_past1000$SIM_mean$ITPC[1:48,1:73]*PLOTTING_VARIABLES$ls_mask[1:48,1:73])
MEAN_ANALYSIS$FIELDS$TEMPlyr <- rbind(DATA_past1000$SIM_mean$TEMP[49:96,1:73]*PLOTTING_VARIABLES$ls_mask[49:96,1:73],
                                      DATA_past1000$SIM_mean$TEMP[1:48,1:73]*PLOTTING_VARIABLES$ls_mask[1:48,1:73])
MEAN_ANALYSIS$FIELDS$PREClyr <- rbind(DATA_past1000$SIM_mean$PREC[49:96,1:73]*PLOTTING_VARIABLES$ls_mask[49:96,1:73],
                                      DATA_past1000$SIM_mean$PREC[1:48,1:73]*PLOTTING_VARIABLES$ls_mask[1:48,1:73])
MEAN_ANALYSIS$FIELDS$SLPRlyr <- rbind(DATA_past1000$SIM_mean$SLPR[49:96,1:73]*PLOTTING_VARIABLES$ls_mask[49:96,1:73],
                                      DATA_past1000$SIM_mean$SLPR[1:48,1:73]*PLOTTING_VARIABLES$ls_mask[1:48,1:73])


MEAN_ANALYSIS$POINTS$CAVElyr_isot <- data.frame(
  lon = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  lat = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)), 
  value = numeric(length(DATA_past1000$CAVES$entity_info$entity_id))
)
MEAN_ANALYSIS$POINTS$CAVElyr_itpc <- data.frame(
  lon = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  lat = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)), 
  value = numeric(length(DATA_past1000$CAVES$entity_info$entity_id))
)

plotting_limit = -10


# Aufpassen. Es ist natürlich nur sinnvoll die Messwerte mit der Simulation zu vergleichen und nicht die Simulation mit der Simulation!!!
for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
  site = DATA_past1000$CAVES$entity_info$site_id[ii]
  entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
  MEAN_ANALYSIS$POINTS$CAVElyr_isot$lon[ii]   = DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site]
  MEAN_ANALYSIS$POINTS$CAVElyr_isot$lat[ii]   = DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site]
  #MEAN_ANALYSIS$POINTS$CAVElyr_isot$value[ii] = mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ISOT, na.rm = T)
  MEAN_ANALYSIS$POINTS$CAVElyr_isot$value[ii] = mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_dw_eq, na.rm = T)
  
  MEAN_ANALYSIS$POINTS$CAVElyr_itpc$lon[ii]   = DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site]
  MEAN_ANALYSIS$POINTS$CAVElyr_itpc$lat[ii]   = DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site]
  #MEAN_ANALYSIS$POINTS$CAVElyr_itpc$value[ii] = mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ITPC, na.rm = T)
  MEAN_ANALYSIS$POINTS$CAVElyr_itpc$value[ii] = mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_dw_eq, na.rm = T)
  
}



# Plot_lyr <- MEAN_ANALYSIS$FIELDS$ISOTlyr
# Plot_lyr[Plot_lyr<plotting_limit] = NA
# 
# GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE = 1
# GLOBAL_STACY_OPTIONS$GLOBAL_FONT_SIZE = 4.5
# GLOBAL_STACY_OPTIONS$GLOBAL_FONT_FACE_TITLE = "plain"
# 
# plot_isot <- STACYmap(gridlyr = Plot_lyr,
#                       ptlyr = MEAN_ANALYSIS$POINTS$CAVElyr_itpc,
#                       zoom = c(-180, -60, 180, 73),
#                       legend_names = list(grid = "Mean d18O in \U2030"), 
#                       graticules = TRUE,
#                       colorscheme = RColorBrewer::brewer.pal(9, 'Reds')) +
#   theme(panel.border = element_blank()) +
#   ggtitle("Mean d18O of HadCM3 prec vs SISALv2 Caves Mean (last 1000y)") +
#   theme(plot.title = element_text(h = 0.5))
# 
# plot_isot
# 
# Plot_lyr <- MEAN_ANALYSIS$FIELDS$ITPClyr
# Plot_lyr[Plot_lyr<plotting_limit] = NA
# 
# plot_itpc <- STACYmap(gridlyr = Plot_lyr,
#                       ptlyr = MEAN_ANALYSIS$POINTS$CAVElyr_itpc,
#                       zoom = c(-180, -60, 180, 73),
#                       legend_names = list(grid = "Mean d18O in \U2030"),
#                                           graticules = TRUE,
#                       colorscheme = RColorBrewer::brewer.pal(9, 'Reds')) +
#   theme(panel.border = element_blank()) +
#   ggtitle("Prec weighted d18O of HadCM3 vs SISALv2 Caves Mean (last 1000y)") +
#   theme(plot.title = element_text(h = 0.5))
# 
# plot_itpc
# 
# plot_isot %>% ggsave(filename = paste('Map_d18O-prec_unfilteredcaves', 'pdf', sep = '.'), plot = ., path = 'Plots/Mean', 
#                        width = 9, height = 5.5, units = 'cm', dpi = 'print', device = "pdf")
# 
# plot_itpc %>% ggsave(filename = paste('Map_d18O-prec-weighted_unfilteredcaves', 'pdf', sep = '.'), plot = ., path = 'Plots/Mean', 
#                      width = 9, height = 5.5, units = 'cm', dpi = 'print', device = "pdf")


## 1.2) Use downsampled as Mean


###################################################################################################

## Masken die verwendet werden zum weiterrechnen!
##
## Da dies noch nicht gemacht ist, verwende ich hier zwei dummy-Masken

mask_mean = logical(length = length(DATA_past1000$CAVES$entity_info$entity_id))
mask_var  = logical(length = length(DATA_past1000$CAVES$entity_info$entity_id))
mask_spec = logical(length = length(DATA_past1000$CAVES$entity_info$entity_id))

for(entity in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
  if(DATA_past1000$CAVES$entity_info$n[entity] > 10 & DATA_past1000$CAVES$entity_info$period[entity] > 600){mask_mean[entity] = T}
  if(DATA_past1000$CAVES$entity_info$n[entity] > 20 & DATA_past1000$CAVES$entity_info$period[entity] > 600){mask_var[entity] = T}
  if(DATA_past1000$CAVES$entity_info$n[entity] > 30 & DATA_past1000$CAVES$entity_info$period[entity] > 600){mask_spec[entity] = T}
}

###################################################################################################

plotting_limit = -10

MEAN_ANALYSIS$POINTS$CAVElyr_isot_used <- data.frame(
  lon = MEAN_ANALYSIS$POINTS$CAVElyr_isot$lon[mask_mean],
  lat = MEAN_ANALYSIS$POINTS$CAVElyr_isot$lat[mask_mean],
  value = MEAN_ANALYSIS$POINTS$CAVElyr_isot$value[mask_mean]
)

MEAN_ANALYSIS$POINTS$CAVElyr_itpc_used <- data.frame(
  lon = MEAN_ANALYSIS$POINTS$CAVElyr_itpc$lon[mask_mean],
  lat = MEAN_ANALYSIS$POINTS$CAVElyr_itpc$lat[mask_mean],
  value = MEAN_ANALYSIS$POINTS$CAVElyr_itpc$value[mask_mean]
)



GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE = 1
GLOBAL_STACY_OPTIONS$GLOBAL_FONT_SIZE = 4.5

Plot_lyr <- MEAN_ANALYSIS$FIELDS$ISOTlyr
Plot_lyr[Plot_lyr<plotting_limit] = NA

plot_isot <- STACYmap(gridlyr = Plot_lyr,
                      ptlyr = MEAN_ANALYSIS$POINTS$CAVElyr_isot_used,
                      zoom = c(-180, -60, 180, 73),
                      legend_names = list(grid = "Mean d18O in \U2030"), 
                      graticules = TRUE,
                      colorscheme = RColorBrewer::brewer.pal(9, 'Reds')) +
  theme(panel.border = element_blank()) +
  ggtitle("Mean d18O of HadCM3 prec vs SISALv2 Caves Mean (last 1000y)") +
  theme(plot.title = element_text(h = 0.5))

#plot_isot

Plot_lyr <- MEAN_ANALYSIS$FIELDS$ITPClyr
Plot_lyr[Plot_lyr<plotting_limit] = NA

plot_itpc <- STACYmap(gridlyr = Plot_lyr,
                      ptlyr = MEAN_ANALYSIS$POINTS$CAVElyr_isot_used,
                      zoom = c(-180, -60, 180, 73),
                      legend_names = list(grid = "Mean d18O in \U2030"), 
                      graticules = TRUE,
                      colorscheme = RColorBrewer::brewer.pal(9, 'Reds')) +
  theme(panel.border = element_blank()) +
  ggtitle("Prec weighted d18O of HadCM3 vs SISALv2 Caves Mean (last 1000y)") +
  theme(plot.title = element_text(h = 0.5))

#plot_itpc

plot_isot %>% ggsave(filename = paste('Map_d18O-prec_filteredcaves', 'pdf', sep = '.'), plot = ., path = 'Plots/Mean', 
                                         width = 9, height = 5.5, units = 'cm', dpi = 'print', device = "pdf")

plot_itpc %>% ggsave(filename = paste('Map_d18O-prec-weighted_filteredcaves', 'pdf', sep = '.'), plot = ., path = 'Plots/Mean', 
                                         width = 9, height = 5.5, units = 'cm', dpi = 'print', device = "pdf")


## Now Calculate Mean and put it in histogram and in latitudinal boxes

matrix <- as.tibble(array(dim = c(length(DATA_past1000$CAVES$entity_info$entity_id), 6)))
colnames(matrix) <- c("entity_id", "mean_sim", "mean_down_sim", "mean_sisal", "diff_full", "diff_down")
MEAN_ANALYSIS$MEAN_DIFF_ISOT <- matrix
MEAN_ANALYSIS$MEAN_DIFF_ITPC <- matrix

for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
  entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
  site = DATA_past1000$CAVES$entity_info$site_id[ii]
  MEAN_ANALYSIS$MEAN_DIFF_ISOT[ii,1] <- DATA_past1000$CAVES$entity_info$entity_id[ii]
  MEAN_ANALYSIS$MEAN_DIFF_ISOT[ii,2] <- mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ISOT, na.rm = T)
  MEAN_ANALYSIS$MEAN_DIFF_ISOT[ii,3] <- mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$d18O, na.rm = T)
  MEAN_ANALYSIS$MEAN_DIFF_ISOT[ii,4] <- mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq, na.rm = T)
  
  MEAN_ANALYSIS$MEAN_DIFF_ISOT[ii,5] <- MEAN_ANALYSIS$MEAN_DIFF_ISOT[ii,2]-MEAN_ANALYSIS$MEAN_DIFF_ISOT[ii,4]
  MEAN_ANALYSIS$MEAN_DIFF_ISOT[ii,6] <- MEAN_ANALYSIS$MEAN_DIFF_ISOT[ii,3]-MEAN_ANALYSIS$MEAN_DIFF_ISOT[ii,4]
  
  MEAN_ANALYSIS$MEAN_DIFF_ITPC[ii,1] <- DATA_past1000$CAVES$entity_info$entity_id[ii]
  MEAN_ANALYSIS$MEAN_DIFF_ITPC[ii,2] <- mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ITPC, na.rm = T)
  MEAN_ANALYSIS$MEAN_DIFF_ITPC[ii,3] <- mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$d18O_pw, na.rm = T)
  MEAN_ANALYSIS$MEAN_DIFF_ITPC[ii,4] <- mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq, na.rm = T)
  
  
  MEAN_ANALYSIS$MEAN_DIFF_ITPC[ii,5] <- MEAN_ANALYSIS$MEAN_DIFF_ISOT[ii,2]-MEAN_ANALYSIS$MEAN_DIFF_ITPC[ii,4]
  MEAN_ANALYSIS$MEAN_DIFF_ITPC[ii,6] <- MEAN_ANALYSIS$MEAN_DIFF_ISOT[ii,3]-MEAN_ANALYSIS$MEAN_DIFF_ITPC[ii,4]
}

## Historgramms to whether difference is significantly different from 0

pdf(file = "Plots/Mean/Hist_mean_diff_full_sim.pdf")
hist(MEAN_ANALYSIS$MEAN_DIFF_ISOT$diff_full*mask_mean, 
     breaks = 9,
     border = "black", 
     prob = TRUE,
     ylim = c(0,0.5),
     xlab = TeX("Diff in $\\delta^{18}O$ level (sim. - record)"),
     main = "")
lines(density(MEAN_ANALYSIS$MEAN_DIFF_ISOT$diff_full*mask_mean, na.rm = T),
      lwd = 2, 
      col = "black")
dev.off()

pdf(file = "Plots/Mean/Hist_mean_diff_down_sim.pdf")
hist(MEAN_ANALYSIS$MEAN_DIFF_ISOT$diff_down*mask_mean, 
     breaks = 9,
     border = "black", 
     prob = TRUE,
     ylim = c(0,0.5),
     xlab = TeX("Diff in $\\delta^{18}O$ level (sim. - record)"),
     main = "")
lines(density(MEAN_ANALYSIS$MEAN_DIFF_ISOT$diff_down*mask_mean, na.rm = T),
      lwd = 2, 
      col = "black")
dev.off()

pdf(file = "Plots/Mean/Hist_mean_diff_full_sim_pc-weighted.pdf")
hist(MEAN_ANALYSIS$MEAN_DIFF_ITPC$diff_full*mask_mean, 
     breaks = 9,
     border = "black", 
     prob = TRUE,
     ylim = c(0,0.5),
     xlab = TeX("Diff in $\\delta^{18}O$ level (sim. - record)"),
     main = "")
lines(density(MEAN_ANALYSIS$MEAN_DIFF_ITPC$diff_full*mask_mean, na.rm = T),
      lwd = 2, 
      col = "black")
dev.off()

pdf(file = "Plots/Mean/Hist_mean_diff_down_sim_pc-weighted.pdf")
hist(MEAN_ANALYSIS$MEAN_DIFF_ITPC$diff_down*mask_mean, 
     breaks = 9,
     border = "black", 
     prob = TRUE,
     ylim = c(0,0.5),
     xlab = TeX("Diff in $\\delta^{18}O$ level (sim. - record)"),
     main = "")
lines(density(MEAN_ANALYSIS$MEAN_DIFF_ITPC$diff_down*mask_mean, na.rm = T),
      lwd = 2, 
      col = "black")
dev.off()


###################################################################################################
## Mean Analysis in Latitudes #####################################################################


#Das muss hier anders gehen! sonst kann man keine Fehlerbalken berechnen...

# matrix <- as.tibble(array(dim = c(18, 6)))
# matrix[,1] <- c(90,80,70,60,50,40,30,20,10,0,-10,-20,-30,-40,-50,-60,-70,-80)
# matrix[,2:4] = 0
# colnames(matrix) <- c("until_lat", "n", "mean", "mean_sd", "mean_pw", "mean_pw_sd")
# 
# USED_DATA <- MEAN_ANALYSIS$MEAN_DIFF_ISOT[mask_mean,]
# USED_DATA_pw <- MEAN_ANALYSIS$MEAN_DIFF_ITPC[mask_mean,]
# 
# counter = 1
# for(lat in c(90,80,70,60,50,40,30,20,10,0,-10,-20,-30,-40,-50,-60,-70,-80)){
#   mean_diffs = list()
#   mean_pw_diffs = list()
#   for(entity in entity_lat$entity_id[entity_lat$lat_band == lat]){
#     mean_diffs = c(mean_diffs, USED_DATA$diff_full[USED_DATA$entity_id == entity])
#     mean_pw_diffs = c(mean_pw_diffs, USED_DATA$diff_full[USED_DATA$entity_id == entity])
#   }
#   if(length(mean_diffs)>0){
#     matrix$n[counter]          = length(mean_diffs)
#     matrix$mean[counter]       = mean(as.numeric(mean_diffs), na.rm = T)
#     matrix$mean_sd[counter]    = sd(as.numeric(mean_diffs), na.rm = T)
#     matrix$mean_pw[counter]    = mean(as.numeric(mean_pw_diffs), na.rm = T)
#     matrix$mean_pw_sd[counter] = sd(as.numeric(mean_pw_diffs), na.rm = T)
#   }else{
#     matrix$n[counter] = matrix$mean[counter] = matrix$mean_sd[counter] = matrix$mean_pw[counter] = matrix$mean_pw_sd[counter] = NA
#   }
#   counter = counter +1
# }
# 
# MEAN_ANALYSIS$MEAN_LATITUDES <- matrix
# 
# # for(ii in 1:dim(USED_DATA)[1]){
# #   site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == USED_DATA$entity_id[ii]]
# #   lat = DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site]
# #   for(jj in length(MEAN_ANALYSIS$MEAN_LATITUDES$until_lat):1){
# #     if (lat <= MEAN_ANALYSIS$MEAN_LATITUDES$until_lat[jj]){
# #       MEAN_ANALYSIS$MEAN_LATITUDES$n[jj] = MEAN_ANALYSIS$MEAN_LATITUDES$n[jj] +1
# #       MEAN_ANALYSIS$MEAN_LATITUDES$mean[jj] = MEAN_ANALYSIS$MEAN_LATITUDES$mean[jj] + USED_DATA$diff_full[ii]
# #       MEAN_ANALYSIS$MEAN_LATITUDES$mean_pw[jj] = MEAN_ANALYSIS$MEAN_LATITUDES$mean_pw[jj] + USED_DATA_pw$diff_full[ii]
# #       break
# #     }
# #   }
# # }
# # 
# # for(ii in 1:length(MEAN_ANALYSIS$MEAN_LATITUDES$until_lat)){
# #   if(MEAN_ANALYSIS$MEAN_LATITUDES$n[ii] > 0){
# #     MEAN_ANALYSIS$MEAN_LATITUDES$mean[ii] = MEAN_ANALYSIS$MEAN_LATITUDES$mean[ii]/MEAN_ANALYSIS$MEAN_LATITUDES$n[ii]
# #     MEAN_ANALYSIS$MEAN_LATITUDES$mean_pw[ii] = MEAN_ANALYSIS$MEAN_LATITUDES$mean_pw[ii]/MEAN_ANALYSIS$MEAN_LATITUDES$n[ii]
# #   }
# # }
# 
# pdf(file = "Plots/Mean/Bar_LAT_mean_diff_sim.pdf", width = 8, height = 5)
# barplot(MEAN_ANALYSIS$MEAN_LATITUDES$mean[3:17],
#         names.arg = MEAN_ANALYSIS$MEAN_LATITUDES$until_lat[3:17],
#         #sub = MEAN_ANALYSIS$MEAN_LATITUDES$mean_sd[3:17],
#         xlab = "Lat.",
#         ylab = "Sim-Rec d18O",
#         main = "Latitudal distribution of mean difference",
#         horiz = F,
#         border = NA, 
#         axes = T, axisnames = T)
# dev.off()
# 
# 
# p<-ggplot(data=as.data.frame(MEAN_ANALYSIS$MEAN_LATITUDES), aes(x=until_lat, y=mean, ymin=mean-mean_sd, ymax=mean+mean_sd)) +
#   geom_bar(stat="identity")+
#   geom_errorbar(stat = "identity", width=2.5)+ 
#   geom_text(aes(label=n), vjust=1.6, color="white", size=3.5)+
#   geom_text(aes(label=n), vjust=-1.6, color="white", size=3.5)+
#   xlim(-70,70)+
#   xlab("latitudes")+
#   ylab(paste0("Mean difference (sim-record) ", expression('\u2030')))+
#   theme_minimal()
# p
# p  %>% ggsave(filename = paste('Bar_LAT_mean_diff_sim_2', 'png', sep = '.'), plot = ., path = 'Plots/Mean', 
#               width = 15, height = 10, units = 'cm', dpi = 'print', device = "png")


###################################################################################################
## Mean Diff over time ############################################################################
###################################################################################################

## 1) Make Mean Diff over for mean_mask for each individual 

MEAN_ANALYSIS$MEAN_DIFF_TIME_d18O <- list(
  interp_age = list(),
  diff_d18O = list(),
  entity_id = list()
)

MEAN_ANALYSIS$MEAN_DIFF_TIME_d18O_pw <- list(
  interp_age = list(),
  diff_d18O = list(),
  entity_id = list()
)

TRY <- list()

for(entity in DATA_past1000$CAVES$entity_info$entity_id[mask_mean]){
  cave_data = DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
  sim_data_ds = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]
  
  MEAN_ANALYSIS$MEAN_DIFF_TIME_d18O$interp_age = c(as.numeric(MEAN_ANALYSIS$MEAN_DIFF_TIME_d18O$interp_age), cave_data$interp_age)
  MEAN_ANALYSIS$MEAN_DIFF_TIME_d18O_pw$interp_age = c(as.numeric(MEAN_ANALYSIS$MEAN_DIFF_TIME_d18O_pw$interp_age), cave_data$interp_age)
  
  MEAN_ANALYSIS$MEAN_DIFF_TIME_d18O$entity_id = c(as.numeric(MEAN_ANALYSIS$MEAN_DIFF_TIME_d18O$entity_id), rep(entity, length(cave_data$interp_age)))
  MEAN_ANALYSIS$MEAN_DIFF_TIME_d18O_pw$entity_id = c(as.numeric(MEAN_ANALYSIS$MEAN_DIFF_TIME_d18O_pw$entity_id), rep(entity, length(cave_data$interp_age)))
  
  MEAN_ANALYSIS$MEAN_DIFF_TIME_d18O$diff_d18O = c(as.numeric(MEAN_ANALYSIS$MEAN_DIFF_TIME_d18O$diff_d18O), 
                                                  sim_data_ds$d18O - cave_data$d18O_dw_eq)
  MEAN_ANALYSIS$MEAN_DIFF_TIME_d18O_pw$diff_d18O = c(as.numeric(MEAN_ANALYSIS$MEAN_DIFF_TIME_d18O_pw$diff_d18O), 
                                                     sim_data_ds$d18O_pw - cave_data$d18O_dw_eq)
  
  data <- data.frame(
    interp_age = cave_data$interp_age,
    diff_d18O  = sim_data_ds$d18O - cave_data$d18O_measurement,
    diff_d18O_pw = sim_data_ds$d18O_pw - cave_data$d18O_measurement
  )

  TRY[[paste0("ENTITY", entity)]] <- data
  #MEAN_ANALYSIS$MEAN_DIFF_TIME[[paste0("ENTITY", entity)]] <- data
}

# plot(0, type = "n", xlim = c(-49,1100), ylim = c(-10,10))
# for(entity in DATA_past1000$CAVES$entity_info$entity_id[mask_mean][50:60]){
#   lines(TRY[[paste0("ENTITY", entity)]]$interp_age, TRY[[paste0("ENTITY", entity)]]$diff_d18O, col = "black")
# }

source("Functions/meandiff_SISAL_SIM_segmentplot.R")

plot_data <- as.data.frame(MEAN_ANALYSIS$MEAN_DIFF_TIME_d18O) %>% filter(entity_id != 14 & entity_id != 49 & entity_id != 60) %>% filter(diff_d18O > -10)

plot <- plotfct_time(plot_data)

plot

plot %>% ggsave(filename = 'Plots/Mean/Time_Diff.pdf', width = 10, height = 20, units = 'cm')

## cool wäre vielleicht noch, das hier nach latitudes zu ordnen


# to get everything in propper bins, we need to smuggle in more information and make them all equidistant with dt=1y

time_caves <- seq(from = -49, to = 1100, by = 1)
matrix_list <- matrix()

for(entity in DATA_past1000$CAVES$entity_info$entity_id[mask_mean]){
  print(entity)
  matrix <- PaleoSpec::MakeEquidistant(TRY[[paste0("ENTITY", entity)]]$interp_age, TRY[[paste0("ENTITY", entity)]]$diff_d18O, dt = 1, time_caves)
  if(entity == 14){
    matrix_list <- matrix
    next
  }
  matrix_list = rbind(matrix_list, matrix)
}

matrix <- colMeans(matrix_list, na.rm = T)


pdf(file = "Plots/Mean/Mean_Time_Diff.pdf", width = 8, height = 5)
plot(time_caves, matrix, type = "l", 
     xlim = c(-49, 1050), xlab = "years BP",# ylim = c(-1.5, 0.5), ylab = TeX("$\\delta 18 O_{Sim}-\\delta 18 O_{Rec}$"),
     main = "Average difference in Sim-Record d18O over time")
abline(h = 0)
dev.off()


#count how many are in one lat box
lat_count <- as.data.frame(DATA_past1000$CAVES$entity_lat) %>% filter(entity_id %in% DATA_past1000$CAVES$entity_info$entity_id[mask_mean]) %>% group_by(lat_band) %>% count()

MEAN_ANALYSIS$LAT_TIME_MEANS <- list()

# Jetzt noch das ganze auf latitudinal bands: 
for(jj in 1:dim(lat_count)[1]){
  lat = lat_count$lat_band[jj]
  data <- array(dim = c(lat_count[lat_count$lat_band == lat,2], 1150))
  entities <- as.data.frame(DATA_past1000$CAVES$entity_lat) %>% filter(entity_id %in% DATA_past1000$CAVES$entity_info$entity_id[mask_mean]) %>% filter(lat_band == lat) %>% select(entity_id)
  entities <- entities$entity_id
  
  for(ii in 1:length(entities)){
    entity = entities[ii]
    sub_data <- PaleoSpec::MakeEquidistant(TRY[[paste0("ENTITY", entity)]]$interp_age, TRY[[paste0("ENTITY", entity)]]$diff_d18O, dt = 1, time_caves)
    data[ii,] <- sub_data
  }
  name = paste0("LAT", lat)
  MEAN_ANALYSIS$LAT_TIME_MEANS[[name]] <- colMeans(data, na.rm = T)
}

MEAN_ANALYSIS$LAT_TIME_MEANS$total <- matrix 

colors = c(rev(RColorBrewer::brewer.pal(7, "Reds")[2:6]), RColorBrewer::brewer.pal(9, "Blues")[2:8])

pdf(file = "Plots/Mean/Mean_Time_Diff_LatBands.pdf", width = 8, height = 5)
plot(time_caves, MEAN_ANALYSIS$LAT_TIME_MEANS$`LAT-40`, col = colors[1], type = "l", ylim = c(-5,5), xlim = c(-10,1030), 
     main = "Mean over latitudinal bands with total",
     ylab = TeX("$\\delta 18 O _{sim} - \\delta 18 O _{rec}"), 
     xlab = "years BP")
lines(time_caves, MEAN_ANALYSIS$LAT_TIME_MEANS$`LAT-40`, col = colors[1], lw = 2)
lines(time_caves, MEAN_ANALYSIS$LAT_TIME_MEANS$`LAT-30`, col = colors[2], lw = 2)
lines(time_caves, MEAN_ANALYSIS$LAT_TIME_MEANS$`LAT-20`, col = colors[3], lw = 2)
lines(time_caves, MEAN_ANALYSIS$LAT_TIME_MEANS$`LAT-10`, col = colors[4], lw = 2)
lines(time_caves, MEAN_ANALYSIS$LAT_TIME_MEANS$`LAT0`, col = colors[5], lw = 2)
lines(time_caves, MEAN_ANALYSIS$LAT_TIME_MEANS$`LAT10`, col = colors[6], lw = 2)
lines(time_caves, MEAN_ANALYSIS$LAT_TIME_MEANS$`LAT20`, col = colors[7], lw = 2)
lines(time_caves, MEAN_ANALYSIS$LAT_TIME_MEANS$`LAT30`, col = colors[8], lw = 2)
lines(time_caves, MEAN_ANALYSIS$LAT_TIME_MEANS$`LAT40`, col = colors[9], lw = 2)
lines(time_caves, MEAN_ANALYSIS$LAT_TIME_MEANS$`LAT50`, col = colors[10], lw = 2)
lines(time_caves, MEAN_ANALYSIS$LAT_TIME_MEANS$`LAT60`, col = colors[11], lw = 2)
lines(time_caves, MEAN_ANALYSIS$LAT_TIME_MEANS$`LAT70`, col = colors[12], lw = 2)
lines(time_caves, MEAN_ANALYSIS$LAT_TIME_MEANS$total, col = "black", lw = 3)
abline(h = 0)

legend("top", legend = c("-40°", "0°", "70°", "total"), col = c(colors[1],colors[5],colors[12], "black"), 
       lw = c(2,2,2,2,2,2,2,2,2,2,2,2,3), horiz = T)

dev.off()


#abline(v = 100)
# 
# fields::image.plot(matrix_list, col = RColorBrewer::brewer.pal(11, 'RdBu')) # HIER MUSS AUCH NOCH 0 WEISS WERDEN. NAs? = SCHWARZ?
# 
# trytry <- data.frame(
#   interp_age = c(TRY$ENTITY14$interp_age, TRY$ENTITY20$interp_age, TRY$ENTITY21$interp_age, TRY$ENTITY33$interp_age),
#   diff_d18O = c(TRY$ENTITY14$diff_d18O, TRY$ENTITY20$diff_d18O, TRY$ENTITY21$diff_d18O, TRY$ENTITY33$diff_d18O),
#   entity_id = c(rep(14,length(TRY$ENTITY14$interp_age)),rep(20,length(TRY$ENTITY20$interp_age)),rep(21,length(TRY$ENTITY21$interp_age)),rep(33,length(TRY$ENTITY33$interp_age)))
# )

remove(colors, counter, entities, entity, ii, jj, lat, length_cave, matrix, site, sub_data, 
       data, fld, mean_diffs, mean_pw_diffs, p, plot, plot_data, plot_isot, plot_itpc, sim_data_ds, test, USED_DATA, USED_DATA_pw)


