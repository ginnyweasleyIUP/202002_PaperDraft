###################################################################################################
## ANALYSIS OF MEAN VALUES --> MEAN BIAS ##########################################################
###################################################################################################
library(dplyr)
library(latex2exp)
source("Functions/STACYmap_5.R")

# ToDo
# [ ] 



## In Simulation find Time of peak high and low (over 30 year average)

LIA_MCA <- list()

LIA_MCA$SIM_full <- list(
  SIM_high_30 = array(dim = c(96,73)),
  SIM_low_30  = array(dim = c(96,73)),
  SIM_climate_30 = array(dim = c(96,73,1150)),
  SIM_high_10 = array(dim = c(96,73)),
  SIM_low_10  = array(dim = c(96,73)),
  SIM_climate_10 = array(dim = c(96,73,1150))
)

for(lon in 1:96){
  for(lat in 1:73){
    print(paste(lon,lat))
    for(time in 1:1150){
      start_time = time -15
      if(start_time < 0){start_time = 0}
      stop_time = time + 15
      if(stop_time > 1150){stop_time = 1150}
      
      LIA_MCA$SIM_full$SIM_climate_30[lon, lat, time] = mean(DATA_past1000$SIM_yearly$TEMP[lon,lat,start_time:stop_time])
      
      start_time = time -5
      if(start_time < 0){start_time = 0}
      stop_time = time + 5
      if(stop_time > 1150){stop_time = 1150}
      
      LIA_MCA$SIM_full$SIM_climate_10[lon, lat, time] = mean(DATA_past1000$SIM_yearly$TEMP[lon,lat,start_time:stop_time])
      
    }
    
  LIA_MCA$SIM_full$SIM_high_30[lon,lat] = 1100 - which.max(LIA_MCA$SIM_full$SIM_climate_30[lon,lat,1:950])
  LIA_MCA$SIM_full$SIM_high_10[lon,lat] = 1100 - which.max(LIA_MCA$SIM_full$SIM_climate_10[lon,lat,1:950])
  LIA_MCA$SIM_full$SIM_low_30[lon,lat]  = 1100 - which.min(LIA_MCA$SIM_full$SIM_climate_30[lon,lat,])
  LIA_MCA$SIM_full$SIM_low_10[lon,lat]  = 1100 - which.min(LIA_MCA$SIM_full$SIM_climate_10[lon,lat,])
    
  }
}


########################################################
## How with caves? do we want high-high or high-low? ###
########################################################

mask_mean_lia = logical(length = 134)

for(entity in 1:134){
  if(DATA_past1000$CAVES$entity_info$n[entity] > 50 & DATA_past1000$CAVES$entity_info$period[entity] > 600){mask_mean_lia[entity] = T}
}


GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE = 1
GLOBAL_STACY_OPTIONS$GLOBAL_FONT_SIZE = 4.5
GLOBAL_STACY_OPTIONS$GLOBAL_FONT_FACE_TITLE = "plain"

Plot_lyr <- rbind(LIA_MCA$SIM_full$SIM_low_30[49:96,1:73]*PLOTTING_VARIABLES$ls_mask[49:96,1:73],
                 LIA_MCA$SIM_full$SIM_low_30[1:48,1:73]*PLOTTING_VARIABLES$ls_mask[1:48,1:73])
#Plot_lyr <- rbind(LIA_MCA$SIM_full$SIM_low_30[49:96,1:73],
#                  LIA_MCA$SIM_full$SIM_low_30[1:48,1:73])

plot <- STACYmap(gridlyr = Plot_lyr,
                 colorscheme = RColorBrewer::brewer.pal(9, "Reds"),
                 zoom = c(-180,-60,180,73),
                 legend_names = list(grid = "years BP"),
                 legend_cb = FALSE,
                 legend_num_breaks = 8) + 
  ggtitle("Year of coldest 30y average climate past millenium") +
  theme(plot.title = element_text(h = 0.5))

plot

plot %>% ggsave(filename = paste('Map_little_ice_age', 'pdf', sep = '.'), plot = ., path = 'Plots/Mean', 
                width = 9, height = 5.5, units = 'cm', dpi = 'print', device = "pdf")


Plot_lyr <- rbind(LIA_MCA$SIM_full$SIM_high_30[49:96,1:73]*PLOTTING_VARIABLES$ls_mask[49:96,1:73],
                  LIA_MCA$SIM_full$SIM_high_30[1:48,1:73]*PLOTTING_VARIABLES$ls_mask[1:48,1:73])
#Plot_lyr <- rbind(LIA_MCA$SIM_full$SIM_high_30[49:96,1:73],
#                  LIA_MCA$SIM_full$SIM_high_30[1:48,1:73])

plot <- STACYmap(gridlyr = Plot_lyr,
                 colorscheme = RColorBrewer::brewer.pal(9, "Reds"),
                 zoom = c(-180,-60,180,73),
                 legend_names = list(grid = "years BP"),
                 legend_cb = FALSE,
                 legend_num_breaks = 7) + 
  ggtitle("Year of warmest 30y average climate past millenium") +
  theme(plot.title = element_text(h = 0.5))


plot

plot %>% ggsave(filename = paste('Map_Medieval_climate_anomaly', 'pdf', sep = '.'), plot = ., path = 'Plots/Mean', 
                width = 9, height = 5.5, units = 'cm', dpi = 'print', device = "pdf")

source("Functions/aw_mean.R")
means <- list()
means$LIA_30 <- simpleawmean(LIA_MCA$SIM_full$SIM_low_30, seq(from = -90, to = 90, length.out = 73))
means$LIA_10 <- simpleawmean(LIA_MCA$SIM_full$SIM_low_10, seq(from = -90, to = 90, length.out = 73))
means$MCA_30 <- simpleawmean(LIA_MCA$SIM_full$SIM_high_30, seq(from = -90, to = 90, length.out = 73))
means$MCA_10 <- simpleawmean(LIA_MCA$SIM_full$SIM_high_10, seq(from = -90, to = 90, length.out = 73))

means$LIA_30_land <- simpleawmean(LIA_MCA$SIM_full$SIM_low_30*PLOTTING_VARIABLES$ls_mask, seq(from = -90, to = 90, length.out = 73))
means$LIA_10_land <- simpleawmean(LIA_MCA$SIM_full$SIM_low_10*PLOTTING_VARIABLES$ls_mask, seq(from = -90, to = 90, length.out = 73))
means$MCA_30_land <- simpleawmean(LIA_MCA$SIM_full$SIM_high_30*PLOTTING_VARIABLES$ls_mask, seq(from = -90, to = 90, length.out = 73))
means$MCA_10_land <- simpleawmean(LIA_MCA$SIM_full$SIM_high_10*PLOTTING_VARIABLES$ls_mask, seq(from = -90, to = 90, length.out = 73))

LIA_MCA$means <- means

remove(means)


#################################################
## 101y mean ####################################
#################################################

LIA_MCA$SIM_full$SIM_high_101_T = array(dim = c(96,73))
LIA_MCA$SIM_full$SIM_low_101_T = array(dim = c(96,73))
LIA_MCA$SIM_full$SIM_climate_101_T = array(dim = c(96,73,1150))
LIA_MCA$SIM_full$SIM_high_101_P = array(dim = c(96,73))
LIA_MCA$SIM_full$SIM_low_101_P = array(dim = c(96,73))
LIA_MCA$SIM_full$SIM_climate_101_P = array(dim = c(96,73,1150))

for(lon in 1:96){
  for(lat in 1:73){
    print(paste(lon,lat))
    for(time in 1:1150){
      start_time = time -50
      if(start_time < 0){start_time = 0}
      stop_time = time + 50
      if(stop_time > 1150){stop_time = 1150}
      
      LIA_MCA$SIM_full$SIM_climate_101_T[lon, lat, time] = mean(DATA_past1000$SIM_yearly$TEMP[lon,lat,start_time:stop_time])
      LIA_MCA$SIM_full$SIM_climate_101_P[lon, lat, time] = mean(DATA_past1000$SIM_yearly$PREC[lon,lat,start_time:stop_time])
    }
    
    LIA_MCA$SIM_full$SIM_high_101_T[lon,lat] = 1100 - which.max(LIA_MCA$SIM_full$SIM_climate_101_T[lon,lat,1:950])
    LIA_MCA$SIM_full$SIM_low_101_T[lon,lat]  = 1100 - which.min(LIA_MCA$SIM_full$SIM_climate_101_T[lon,lat,])
    LIA_MCA$SIM_full$SIM_high_101_P[lon,lat] = 1100 - which.max(LIA_MCA$SIM_full$SIM_climate_101_P[lon,lat,1:950])
    LIA_MCA$SIM_full$SIM_low_101_P[lon,lat]  = 1100 - which.min(LIA_MCA$SIM_full$SIM_climate_101_P[lon,lat,1:950])
    
  }
}

GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE = 1
GLOBAL_STACY_OPTIONS$GLOBAL_FONT_SIZE = 4.5
GLOBAL_STACY_OPTIONS$GLOBAL_FONT_FACE_TITLE = "plain"

Plot_lyr <- rbind(LIA_MCA$SIM_full$SIM_low_101_T[49:96,1:73]*PLOTTING_VARIABLES$ls_mask[49:96,1:73],
                  LIA_MCA$SIM_full$SIM_low_101_T[1:48,1:73]*PLOTTING_VARIABLES$ls_mask[1:48,1:73])
#Plot_lyr <- rbind(LIA_MCA$SIM_full$SIM_low_30[49:96,1:73],
#                  LIA_MCA$SIM_full$SIM_low_30[1:48,1:73])

plot <- STACYmap(gridlyr = Plot_lyr,
                 colorscheme = RColorBrewer::brewer.pal(9, "Reds"),
                 zoom = c(-180,-60,180,73),
                 legend_names = list(grid = "years BP")) + 
  ggtitle("Year of coldest 30y average climate past millenium") +
  theme(plot.title = element_text(h = 0.5),
        panel.border = element_blank())

source("Functions/aw_mean.R")
LIA_MCA$means$LIA_101_T <- simpleawmean(LIA_MCA$SIM_full$SIM_low_101_T, seq(from = -90, to = 90, length.out = 73))
LIA_MCA$means$MCA_101_T <- simpleawmean(LIA_MCA$SIM_full$SIM_high_101_T, seq(from = -90, to = 90, length.out = 73))
LIA_MCA$means$LIA_101_P <- simpleawmean(LIA_MCA$SIM_full$SIM_low_101_P, seq(from = -90, to = 90, length.out = 73))
LIA_MCA$means$MCA_101_P <- simpleawmean(LIA_MCA$SIM_full$SIM_high_101_P, seq(from = -90, to = 90, length.out = 73))

LIA_MCA$means$LIA_101_land_T <- simpleawmean(LIA_MCA$SIM_full$SIM_low_101_T*PLOTTING_VARIABLES$ls_mask, seq(from = -90, to = 90, length.out = 73))
LIA_MCA$means$MCA_101_land_T <- simpleawmean(LIA_MCA$SIM_full$SIM_high_101_T*PLOTTING_VARIABLES$ls_mask, seq(from = -90, to = 90, length.out = 73))
LIA_MCA$means$LIA_101_land_P <- simpleawmean(LIA_MCA$SIM_full$SIM_low_101_P*PLOTTING_VARIABLES$ls_mask, seq(from = -90, to = 90, length.out = 73))
LIA_MCA$means$MCA_101_land_P <- simpleawmean(LIA_MCA$SIM_full$SIM_high_101_P*PLOTTING_VARIABLES$ls_mask, seq(from = -90, to = 90, length.out = 73))

# dann fehlt noch sd!!! --> Function schon geschrieben!

sd <- list()
sd$LIA_30 <- simpleawsd(LIA_MCA$SIM_full$SIM_low_30, seq(from = -90, to = 90, length.out = 73))
sd$LIA_10 <- simpleawsd(LIA_MCA$SIM_full$SIM_low_10, seq(from = -90, to = 90, length.out = 73))
sd$MCA_30 <- simpleawsd(LIA_MCA$SIM_full$SIM_high_30, seq(from = -90, to = 90, length.out = 73))
sd$MCA_10 <- simpleawsd(LIA_MCA$SIM_full$SIM_high_10, seq(from = -90, to = 90, length.out = 73))
sd$LIA_101_T <- simpleawsd(LIA_MCA$SIM_full$SIM_low_101_T, seq(from = -90, to = 90, length.out = 73))
sd$LIA_101_P <- simpleawsd(LIA_MCA$SIM_full$SIM_low_101_P, seq(from = -90, to = 90, length.out = 73))
sd$MCA_101_T <- simpleawsd(LIA_MCA$SIM_full$SIM_high_101_T, seq(from = -90, to = 90, length.out = 73))
sd$MCA_101_P <- simpleawsd(LIA_MCA$SIM_full$SIM_high_101_P, seq(from = -90, to = 90, length.out = 73))

sd$LIA_30_land <- simpleawsd(LIA_MCA$SIM_full$SIM_low_30*PLOTTING_VARIABLES$ls_mask, seq(from = -90, to = 90, length.out = 73))
sd$LIA_10_land <- simpleawsd(LIA_MCA$SIM_full$SIM_low_10*PLOTTING_VARIABLES$ls_mask, seq(from = -90, to = 90, length.out = 73))
sd$MCA_30_land <- simpleawsd(LIA_MCA$SIM_full$SIM_high_30*PLOTTING_VARIABLES$ls_mask, seq(from = -90, to = 90, length.out = 73))
sd$MCA_10_land <- simpleawsd(LIA_MCA$SIM_full$SIM_high_10*PLOTTING_VARIABLES$ls_mask, seq(from = -90, to = 90, length.out = 73))
sd$LIA_101_T_land <- simpleawsd(LIA_MCA$SIM_full$SIM_low_101_T*PLOTTING_VARIABLES$ls_mask, seq(from = -90, to = 90, length.out = 73))
sd$LIA_101_P_land <- simpleawsd(LIA_MCA$SIM_full$SIM_low_101_P*PLOTTING_VARIABLES$ls_mask, seq(from = -90, to = 90, length.out = 73))
sd$MCA_101_T_land <- simpleawsd(LIA_MCA$SIM_full$SIM_high_101_T*PLOTTING_VARIABLES$ls_mask, seq(from = -90, to = 90, length.out = 73))
sd$MCA_101_P_land <- simpleawsd(LIA_MCA$SIM_full$SIM_high_101_P*PLOTTING_VARIABLES$ls_mask, seq(from = -90, to = 90, length.out = 73))

LIA_MCA$sd <- sd

remove(sd)
