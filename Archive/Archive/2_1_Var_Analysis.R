###################################################################################################
## ANALYSIS OF VAR VALUES --> VAR BIAS ############################################################
###################################################################################################

# TODO:
# [ ] soll es einen Hintergrund für die Karte geben und wenn ja welchen
# [ ] wenn log dargestellt wird, dann wäre es toll, wenn "0" als weiß in der Mitte wäre
#       --> dazu muss das centercolor funktionieren und eingebaut werden. Vielleicht muss man dafür nochmal STACYmaps von Mo neu runterladen
#       --> ansonsten muss das als Punktekarte drüber gelegt werden mit einem eigenen scale_fill_color wo die Farbskala definiert ist.
# [ ] wie norme ich die Varianz für Temperatur und Niederschlag, damit ich sie mit der Varianz von den Records vergleichen kann
# [ ] berechne die Fehlerbalken für die latitudinal bands


library(dplyr)
library(latex2exp)
source("Functions/STACYmap_3.R")

VAR_ANALYSIS <- list()

## 1) Create map without any Power Test --> just include all ##
##
## 1.1) Use full simulation over 1000 years

## Create var Matrix:

# Hier ist es vielleicht die Frage, ob es um Temperaturvariabilität, Niederschlagsvariabilität oder Isotopievariabilität geht!!!

VAR_ANALYSIS$FIELDS$ISOTlyr <- array(dim = c(96,73))
VAR_ANALYSIS$FIELDS$ITPClyr <- array(dim = c(96,73))
VAR_ANALYSIS$FIELDS$TEMPlyr <- array(dim = c(96,73))
VAR_ANALYSIS$FIELDS$PREClyr <- array(dim = c(96,73))

for(lon in 1:96){
  for(lat in 1:73){
    VAR_ANALYSIS$FIELDS$ISOTlyr[lon,lat] = var(DATA_past1000$SIM_yearly$ISOT[lon,lat,])
    VAR_ANALYSIS$FIELDS$ITPClyr[lon,lat] = var(DATA_past1000$SIM_yearly$ITPC[lon,lat,])
    VAR_ANALYSIS$FIELDS$TEMPlyr[lon,lat] = var(DATA_past1000$SIM_yearly$TEMP[lon,lat,])
    VAR_ANALYSIS$FIELDS$PREClyr[lon,lat] = var(DATA_past1000$SIM_yearly$PREC[lon,lat,])
  }
}



VAR_ANALYSIS$FIELDS$ISOTlyr <- rbind(VAR_ANALYSIS$FIELDS$ISOT[49:96,1:73],
                                     VAR_ANALYSIS$FIELDS$ISOT[1:48, 1:73])

VAR_ANALYSIS$FIELDS$ITPClyr <- rbind(VAR_ANALYSIS$FIELDS$ITPC[49:96,1:73],
                                     VAR_ANALYSIS$FIELDS$ITPC[1:48, 1:73])

VAR_ANALYSIS$FIELDS$TEMPlyr <- rbind(VAR_ANALYSIS$FIELDS$TEMP[49:96,1:73],
                                     VAR_ANALYSIS$FIELDS$TEMP[1:48, 1:73])

VAR_ANALYSIS$FIELDS$PREClyr <- rbind(VAR_ANALYSIS$FIELDS$PREC[49:96,1:73],
                                     VAR_ANALYSIS$FIELDS$PREC[1:48, 1:73])

########################################################################################################################################

VAR_ANALYSIS$POINTS$CAVElyr <- data.frame(
  lon = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  lat = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)), 
  value_record = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  value_VR_temp = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  value_VR_temp_ds = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  value_VR_prec = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  value_VR_prec_ds = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  value_VR_isot = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  value_VR_isot_ds = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  value_VR_itpc = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  value_VR_itpc_ds = numeric(length(DATA_past1000$CAVES$entity_info$entity_id))
)

plotting_limit = -10

Var_temp_norm <- simpleawmean(VAR_ANALYSIS$FIELDS$TEMPlyr, seq(from = -90, to = 90, length.out = 73))
Var_prec_norm <- simpleawmean(VAR_ANALYSIS$FIELDS$PREClyr, seq(from = -90, to = 90, length.out = 73))
Var_isot_norm <- simpleawmean(VAR_ANALYSIS$FIELDS$ISOTlyr, seq(from = -90, to = 90, length.out = 73))
Var_itpc_norm <- simpleawmean(VAR_ANALYSIS$FIELDS$ITPClyr, seq(from = -90, to = 90, length.out = 73))

for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
  site = DATA_past1000$CAVES$entity_info$site_id[ii]
  entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
  VAR_ANALYSIS$POINTS$CAVElyr$lon[ii]   = DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site]
  VAR_ANALYSIS$POINTS$CAVElyr$lat[ii]   = DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site]
  VAR_ANALYSIS$POINTS$CAVElyr$value[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement, na.rm = T)
  #Temp
  # Take normalized Variances
  VAR_ANALYSIS$POINTS$CAVElyr$value_VR_temp[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement, na.rm = T)/
    var(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$TEMP, na.rm = T)*Var_temp_norm/Var_isot_norm
  VAR_ANALYSIS$POINTS$CAVElyr$value_VR_temp_ds[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement, na.rm = T)/
    var(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$Temp, na.rm = T) * Var_temp_norm/Var_isot_norm
  #Prec
  #Take normalized variances
  VAR_ANALYSIS$POINTS$CAVElyr$value_VR_prec[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement, na.rm = T)/
    var(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$PREC, na.rm = T) * Var_prec_norm/Var_isot_norm
  VAR_ANALYSIS$POINTS$CAVElyr$value_VR_prec_ds[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement, na.rm = T)/
    var(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$Prec, na.rm = T) * Var_prec_norm/Var_isot_norm
  #Isot "d18O"
  VAR_ANALYSIS$POINTS$CAVElyr$value_VR_isot[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement, na.rm = T)/
    var(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ISOT, na.rm = T)
  VAR_ANALYSIS$POINTS$CAVElyr$value_VR_isot_ds[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement, na.rm = T)/
    var(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$d18O, na.rm = T)
  #ITPC  "d18O_pw"
  VAR_ANALYSIS$POINTS$CAVElyr$value_VR_itpc[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement, na.rm = T)/
    var(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ITPC, na.rm = T)
  VAR_ANALYSIS$POINTS$CAVElyr$value_VR_itpc_ds[ii] = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement, na.rm = T)/
    var(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$d18O_pw, na.rm = T)
}





Point_Lyr <- data.frame(
  lon = VAR_ANALYSIS$POINTS$CAVElyr$lon,
  lat = VAR_ANALYSIS$POINTS$CAVElyr$lat,
  value = log(VAR_ANALYSIS$POINTS$CAVElyr$value_VR_temp)
)

Plot_lyr <- VAR_ANALYSIS$FIELDS$TEMPlyr
Plot_lyr[Plot_lyr<plotting_limit] = NA


GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE = 1
GLOBAL_STACY_OPTIONS$GLOBAL_FONT_SIZE = 4.5
GLOBAL_STACY_OPTIONS$GLOBAL_FONT_FACE_TITLE = "plain"



plot <- STACYmap(ptlyr = Point_Lyr,
                 zoom = c(-180, -60, 180, 73),
                 legend_names = list(grid  = "", pt = "log(Var_Re/Var_Sim)"),
                 colorscheme = "temp",
                 centercolor = 0) +
  ggtitle("Variance Ratio: past 1000y Temperature")+
  theme(plot.title = element_text(h = 0.5))

plot

###################################################################################################

GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE = 1
GLOBAL_STACY_OPTIONS$GLOBAL_FONT_SIZE = 4.5
GLOBAL_STACY_OPTIONS$GLOBAL_FONT_FACE_TITLE = "plain"


Point_Lyr <- data.frame(
  lon = VAR_ANALYSIS$POINTS$CAVElyr$lon[mask_var],
  lat = VAR_ANALYSIS$POINTS$CAVElyr$lat[mask_var],
  #value = log(VAR_ANALYSIS$POINTS$CAVElyr$value_VR_temp[mask_var])
  #value = log(VAR_ANALYSIS$POINTS$CAVElyr$value_VR_temp_ds[mask_var])
  #value = log(VAR_ANALYSIS$POINTS$CAVElyr$value_VR_prec[mask_var])
  #value = log(VAR_ANALYSIS$POINTS$CAVElyr$value_VR_prec_ds[mask_var])
  #value = log(VAR_ANALYSIS$POINTS$CAVElyr$value_VR_isot[mask_var])
  #value = log(VAR_ANALYSIS$POINTS$CAVElyr$value_VR_isot_ds[mask_var])
  #value = log(VAR_ANALYSIS$POINTS$CAVElyr$value_VR_itpc[mask_var])
  value = log(VAR_ANALYSIS$POINTS$CAVElyr$value_VR_itpc_ds[mask_var])
)

plot <- STACYmap(ptlyr = Point_Lyr,
                 zoom = c(-180, -60, 180, 73),
                 legend_names = list(grid  = "", pt = TeX("log($$Var_{Rec}/Var_{Sim}$$)")),
                 colorscheme = "temp",
                 centercolor = 0) +
  #ggtitle("Variance Ratio: past 1000y Temp")+
  #ggtitle("Variance Ratio: past 1000y Temp downsampled")+
  #ggtitle("Variance Ratio: past 1000y Prec")+
  #ggtitle("Variance Ratio: past 1000y Prec downsampled")+
  #ggtitle("Variance Ratio: past 1000y d18O")+
  #ggtitle("Variance Ratio: past 1000y d18O downsampled")+
  #ggtitle("Variance Ratio: past 1000y d18O prec-weighted")+
  ggtitle("Variance Ratio: past 1000y d18O prec-weighted downsampled")+
  theme(plot.title = element_text(h = 0.5))

plot

#plot %>% ggsave(filename = paste('Map_VarRatio_temp_full_filteredcaves', 'pdf', sep = '.'),
#                plot = ., path = 'Plots/Variance', width = 9, height = 5.5, units = 'cm', dpi = 'print', device = "pdf")

#plot %>% ggsave(filename = paste('Map_VarRatio_temp_down_filteredcaves', 'pdf', sep = '.'),
#                 plot = ., path = 'Plots/Variance', width = 9, height = 5.5, units = 'cm', dpi = 'print', device = "pdf")
# plot %>% ggsave(filename = paste('Map_VarRatio_prec_full_filteredcaves', 'pdf', sep = '.'),
#                  plot = ., path = 'Plots/Variance', width = 9, height = 5.5, units = 'cm', dpi = 'print', device = "pdf")
# plot %>% ggsave(filename = paste('Map_VarRatio_prec_down_filteredcaves', 'pdf', sep = '.'),
#                 plot = ., path = 'Plots/Variance', width = 9, height = 5.5, units = 'cm', dpi = 'print', device = "pdf")
# plot %>% ggsave(filename = paste('Map_VarRatio_isot_full_filteredcaves', 'pdf', sep = '.'),
#                 plot = ., path = 'Plots/Variance', width = 9, height = 5.5, units = 'cm', dpi = 'print', device = "pdf")
# plot %>% ggsave(filename = paste('Map_VarRatio_isot_down_filteredcaves', 'pdf', sep = '.'),
#                 plot = ., path = 'Plots/Variance', width = 9, height = 5.5, units = 'cm', dpi = 'print', device = "pdf")
# plot %>% ggsave(filename = paste('Map_VarRatio_itpc_full_filteredcaves', 'pdf', sep = '.'),
#                 plot = ., path = 'Plots/Variance', width = 9, height = 5.5, units = 'cm', dpi = 'print', device = "pdf")
plot %>% ggsave(filename = paste('Map_VarRatio_itpc_down_filteredcaves', 'pdf', sep = '.'),
                plot = ., path = 'Plots/Variance', width = 9, height = 5.5, units = 'cm', dpi = 'print', device = "pdf")

###################################################################################################
## Now TRY to bring some structure into it --> make histograms ####################################
###################################################################################################
## Now Calculate Mean and put it in histogram and in latitudinal boxes

Var_temp_norm <- simpleawmean(VAR_ANALYSIS$FIELDS$TEMPlyr, seq(from = -90, to = 90, length.out = 73))
Var_prec_norm <- simpleawmean(VAR_ANALYSIS$FIELDS$PREClyr, seq(from = -90, to = 90, length.out = 73))
Var_isot_norm <- simpleawmean(VAR_ANALYSIS$FIELDS$ISOTlyr, seq(from = -90, to = 90, length.out = 73))
Var_itpc_norm <- simpleawmean(VAR_ANALYSIS$FIELDS$ITPClyr, seq(from = -90, to = 90, length.out = 73))


matrix <- as.tibble(array(dim = c(length(DATA_past1000$CAVES$entity_info$entity_id), 9)))
colnames(matrix) <- c("entity_id", 
                      "var_ratio_sim_temp", "var_ratio_sim_temp_ds",
                      "var_ratio_sim_prec", "var_ratio_sim_prec_ds",
                      "var_ratio_sim_isot", "var_ratio_sim_isot_ds",
                      "var_ratio_sim_itpc", "var_ratio_sim_itpc_ds")
VAR_ANALYSIS$VAR_RATIOS <- matrix

for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
  entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
  site = DATA_past1000$CAVES$entity_info$site_id[ii]
  var_rec = var(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_measurement)
  VAR_ANALYSIS$VAR_RATIOS[ii,1] <- DATA_past1000$CAVES$entity_info$entity_id[ii]
  VAR_ANALYSIS$VAR_RATIOS[ii,2] <- var_rec/var(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$TEMP, na.rm = T) *Var_temp_norm/Var_isot_norm
  VAR_ANALYSIS$VAR_RATIOS[ii,3] <- var_rec/var(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$Temp, na.rm = T) *Var_temp_norm/Var_isot_norm
  VAR_ANALYSIS$VAR_RATIOS[ii,4] <- var_rec/var(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$PREC, na.rm = T)*Var_prec_norm/Var_isot_norm
  VAR_ANALYSIS$VAR_RATIOS[ii,5] <- var_rec/var(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$Prec, na.rm = T)*Var_prec_norm/Var_isot_norm
  VAR_ANALYSIS$VAR_RATIOS[ii,6] <- var_rec/var(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ISOT, na.rm = T)
  VAR_ANALYSIS$VAR_RATIOS[ii,7] <- var_rec/var(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$d18O, na.rm = T)
  VAR_ANALYSIS$VAR_RATIOS[ii,8] <- var_rec/var(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE",site)]]$ITPC, na.rm = T) *Var_itpc_norm/Var_isot_norm
  VAR_ANALYSIS$VAR_RATIOS[ii,9] <- var_rec/var(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$d18O_pw, na.rm = T)*Var_itpc_norm/Var_isot_norm
  
}

remove(ii, entity, site, var_rec, Var_isot_norm, Var_itpc_norm, Var_prec_norm, Var_temp_norm)

## Historgramms to whether variance over all larger in sim or in rec

pdf(file = "Plots/Variance/Hist_var_temp_prec.pdf")
par(mfrow=c(2,2))#, mai = c(rep(spacing, 4)), mar = c(3,3,2,0.5))

hist(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_temp*mask_mean), 
     breaks = 9, border = "black", prob = TRUE, 
     ylim = c(0,0.5), 
     xlab = TeX("log($$Var_{Rec}/Var_{Sim}$$)"),
     main = "Temp - full")
lines(density(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_temp*mask_mean), na.rm = T),
      lwd = 2, col = "black")

hist(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_temp_ds*mask_mean), 
     breaks = 9, border = "black", prob = TRUE, 
     ylim = c(0,0.5), 
     xlab = TeX("log($$Var_{Rec}/Var_{Sim}$$)"),
     main = "Temp - down")
lines(density(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_temp_ds*mask_mean), na.rm = T),
      lwd = 2, col = "black")

hist(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_prec*mask_mean), 
     breaks = 9, border = "black", prob = TRUE, 
     ylim = c(0,0.5), 
     xlab = TeX("log($$Var_{Rec}/Var_{Sim}$$)"),
     main = "Prec - full")
lines(density(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_prec*mask_mean), na.rm = T),
      lwd = 2, col = "black")

hist(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_prec_ds*mask_mean), 
     breaks = 9, border = "black", prob = TRUE, 
     ylim = c(0,0.5), 
     xlab = TeX("log($$Var_{Rec}/Var_{Sim}$$)"),
     main = "Prec - down")
lines(density(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_prec_ds*mask_mean), na.rm = T),
      lwd = 2, col = "black")

dev.off()

pdf(file = "Plots/Variance/Hist_var_isot_itpc.pdf")
par(mfrow=c(2,2))#, mai = c(rep(spacing, 4)), mar = c(3,3,2,0.5))

hist(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_isot*mask_mean), 
     breaks = 9, border = "black", prob = TRUE, 
     ylim = c(0,0.5), 
     xlab = TeX("log($$Var_{Rec}/Var_{Sim}$$)"),
     main = "ISOT - full")
lines(density(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_isot*mask_mean), na.rm = T),
      lwd = 2, col = "black")

hist(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_isot_ds*mask_mean), 
     breaks = 9, border = "black", prob = TRUE, 
     ylim = c(0,0.5), 
     xlab = TeX("log($$Var_{Rec}/Var_{Sim}$$)"),
     main = "ISOT - down")
lines(density(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_isot_ds*mask_mean), na.rm = T),
      lwd = 2, col = "black")

hist(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_itpc*mask_mean), 
     breaks = 9, border = "black", prob = TRUE, 
     ylim = c(0,0.5), 
     xlab = TeX("log($$Var_{Rec}/Var_{Sim}$$)"),
     main = "ITPC  - full")
lines(density(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_itpc*mask_mean), na.rm = T),
      lwd = 2, col = "black")

hist(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_itpc_ds*mask_mean), 
     breaks = 9, border = "black", prob = TRUE, 
     ylim = c(0,0.5), 
     xlab = TeX("log($$Var_{Rec}/Var_{Sim}$$)"),
     main = "ITPC - down")
lines(density(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_itpc_ds*mask_mean), na.rm = T),
      lwd = 2, col = "black")

dev.off()

###################################################################################################
## Now start with latitudinal distribution ######
#################################################

matrix <- as.tibble(array(dim = c(18, 10)))
matrix[,1] <- c(90,80,70,60,50,40,30,20,10,0,-10,-20,-30,-40,-50,-60,-70,-80)
matrix[,2:10] = 0
colnames(matrix) <- c("until_lat", "n", 
                      "ratio_isot", "ratio_isot_ds",
                      "ratio_itpc", "ratio_itpc_ds",
                      "ratio_temp", "ratio_temp_ds",
                      "ratio_prec", "ratio_prec_ds")

VAR_ANALYSIS$VAR_LATITUDES <- matrix 

remove(matrix)

USED_DATA <- VAR_ANALYSIS$VAR_RATIOS[mask_var,]

for(ii in 1:dim(USED_DATA)[1]){
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == USED_DATA$entity_id[ii]]
  lat = DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site]
  for(jj in length(VAR_ANALYSIS$VAR_LATITUDES$until_lat):1){
    if (lat <= VAR_ANALYSIS$VAR_LATITUDES$until_lat[jj]){
      if(any(is.na(USED_DATA[ii,]))){
        break
      }
      VAR_ANALYSIS$VAR_LATITUDES$n[jj] = VAR_ANALYSIS$VAR_LATITUDES$n[jj] +1
      VAR_ANALYSIS$VAR_LATITUDES$ratio_isot[jj]    = VAR_ANALYSIS$VAR_LATITUDES$ratio_isot[jj]    + USED_DATA$var_ratio_sim_isot[ii]
      VAR_ANALYSIS$VAR_LATITUDES$ratio_isot_ds[jj] = VAR_ANALYSIS$VAR_LATITUDES$ratio_isot_ds[jj] + USED_DATA$var_ratio_sim_isot_ds[ii]
      VAR_ANALYSIS$VAR_LATITUDES$ratio_itpc[jj]    = VAR_ANALYSIS$VAR_LATITUDES$ratio_itpc[jj]    + USED_DATA$var_ratio_sim_itpc[ii]
      VAR_ANALYSIS$VAR_LATITUDES$ratio_itpc_ds[jj] = VAR_ANALYSIS$VAR_LATITUDES$ratio_itpc_ds[jj] + USED_DATA$var_ratio_sim_itpc_ds[ii]
      VAR_ANALYSIS$VAR_LATITUDES$ratio_temp[jj]    = VAR_ANALYSIS$VAR_LATITUDES$ratio_temp[jj]    + USED_DATA$var_ratio_sim_temp[ii]
      VAR_ANALYSIS$VAR_LATITUDES$ratio_temp_ds[jj] = VAR_ANALYSIS$VAR_LATITUDES$ratio_temp_ds[jj] + USED_DATA$var_ratio_sim_temp_ds[ii]
      VAR_ANALYSIS$VAR_LATITUDES$ratio_prec[jj]    = VAR_ANALYSIS$VAR_LATITUDES$ratio_prec[jj]    + USED_DATA$var_ratio_sim_prec[ii]
      VAR_ANALYSIS$VAR_LATITUDES$ratio_prec_ds[jj] = VAR_ANALYSIS$VAR_LATITUDES$ratio_prec_ds[jj] + USED_DATA$var_ratio_sim_prec_ds[ii]

      break
    }
  }
}

for(ii in 1:length(VAR_ANALYSIS$VAR_LATITUDES$until_lat)){
  if(VAR_ANALYSIS$VAR_LATITUDES$n[ii] > 0){
    VAR_ANALYSIS$VAR_LATITUDES$ratio_isot[ii]    = VAR_ANALYSIS$VAR_LATITUDES$ratio_isot[ii]/   VAR_ANALYSIS$VAR_LATITUDES$n[ii]
    VAR_ANALYSIS$VAR_LATITUDES$ratio_isot_ds[ii] = VAR_ANALYSIS$VAR_LATITUDES$ratio_isot_ds[ii]/VAR_ANALYSIS$VAR_LATITUDES$n[ii]
    VAR_ANALYSIS$VAR_LATITUDES$ratio_itpc[ii]    = VAR_ANALYSIS$VAR_LATITUDES$ratio_itpc[ii]/   VAR_ANALYSIS$VAR_LATITUDES$n[ii]
    VAR_ANALYSIS$VAR_LATITUDES$ratio_itpc_ds[ii] = VAR_ANALYSIS$VAR_LATITUDES$ratio_itpc_ds[ii]/VAR_ANALYSIS$VAR_LATITUDES$n[ii]
    VAR_ANALYSIS$VAR_LATITUDES$ratio_temp[ii]    = VAR_ANALYSIS$VAR_LATITUDES$ratio_temp[ii]/   VAR_ANALYSIS$VAR_LATITUDES$n[ii]
    VAR_ANALYSIS$VAR_LATITUDES$ratio_temp_ds[ii] = VAR_ANALYSIS$VAR_LATITUDES$ratio_temp_ds[ii]/VAR_ANALYSIS$VAR_LATITUDES$n[ii]
    VAR_ANALYSIS$VAR_LATITUDES$ratio_prec[ii]    = VAR_ANALYSIS$VAR_LATITUDES$ratio_prec[ii]/   VAR_ANALYSIS$VAR_LATITUDES$n[ii]
    VAR_ANALYSIS$VAR_LATITUDES$ratio_prec_ds[ii] = VAR_ANALYSIS$VAR_LATITUDES$ratio_prec_ds[ii]/VAR_ANALYSIS$VAR_LATITUDES$n[ii]

  }
}

pdf(file = "Plots/Variance/Bar_LAT_ratio_ISO.pdf", width = 8, height = 5)
par(mfrow=c(2,2))#, mai = c(rep(spacing, 4)), mar = c(3,3,2,0.5))
barplot(log(VAR_ANALYSIS$VAR_LATITUDES$ratio_isot[3:17]),
        names.arg = VAR_ANALYSIS$VAR_LATITUDES$until_lat[3:17],
        xlab = "Lat.",
        ylim = c(-3,3), ylab = TeX("log($$Var_{Rec}/Var_{Sim}$$)"),
        main = "isot full",
        horiz = F,border = NA, axes = T, axisnames = T)
barplot(log(VAR_ANALYSIS$VAR_LATITUDES$ratio_isot_ds[3:17]),
        names.arg = VAR_ANALYSIS$VAR_LATITUDES$until_lat[3:17],
        xlab = "Lat.",
        ylim = c(-3,3), ylab = TeX("log($$Var_{Rec}/Var_{Sim}$$)"),
        main = "isot down",
        horiz = F,border = NA, axes = T, axisnames = T)
barplot(log(VAR_ANALYSIS$VAR_LATITUDES$ratio_itpc[3:17]),
        names.arg = VAR_ANALYSIS$VAR_LATITUDES$until_lat[3:17],
        xlab = "Lat.",
        ylim = c(-3,3), ylab = TeX("log($$Var_{Rec}/Var_{Sim}$$)"),
        main = "itpc full",
        horiz = F,border = NA, axes = T, axisnames = T)
barplot(log(VAR_ANALYSIS$VAR_LATITUDES$ratio_itpc_ds[3:17]),
        names.arg = VAR_ANALYSIS$VAR_LATITUDES$until_lat[3:17],
        xlab = "Lat.",
        ylim = c(-3,3), ylab = TeX("log($$Var_{Rec}/Var_{Sim}$$)"),
        main = "itpc down",
        horiz = F,border = NA, axes = T, axisnames = T)
dev.off()

pdf(file = "Plots/Variance/Bar_LAT_ratio_TP.pdf", width = 8, height = 5)
par(mfrow=c(2,2))#, mai = c(rep(spacing, 4)), mar = c(3,3,2,0.5))
barplot(log(VAR_ANALYSIS$VAR_LATITUDES$ratio_temp[3:17]),
        names.arg = VAR_ANALYSIS$VAR_LATITUDES$until_lat[3:17],
        xlab = "Lat.",
        ylim = c(-3,3), ylab = TeX("log($$Var_{Rec}/Var_{Sim}$$)"),
        main = "temp full",
        horiz = F,border = NA, axes = T, axisnames = T)
barplot(log(VAR_ANALYSIS$VAR_LATITUDES$ratio_temp_ds[3:17]),
        names.arg = VAR_ANALYSIS$VAR_LATITUDES$until_lat[3:17],
        xlab = "Lat.",
        ylim = c(-3,3), ylab = TeX("log($$Var_{Rec}/Var_{Sim}$$)"),
        main = "temp down",
        horiz = F,border = NA, axes = T, axisnames = T)
barplot(log(VAR_ANALYSIS$VAR_LATITUDES$ratio_prec[3:17]),
        names.arg = VAR_ANALYSIS$VAR_LATITUDES$until_lat[3:17],
        xlab = "Lat.",
        ylim = c(-30,30), ylab = TeX("log($$Var_{Rec}/Var_{Sim}$$)"),
        main = "prec full",
        horiz = F,border = NA, axes = T, axisnames = T)
barplot(log(VAR_ANALYSIS$VAR_LATITUDES$ratio_prec_ds[3:17]),
        names.arg = VAR_ANALYSIS$VAR_LATITUDES$until_lat[3:17],
        xlab = "Lat.",
        ylim = c(-30,30), ylab = TeX("log($$Var_{Rec}/Var_{Sim}$$)"),
        main = "prec down",
        horiz = F,border = NA, axes = T, axisnames = T)
dev.off()
