#################################################
## ANALYSIS SEASONS #############################
#################################################

library(plyr)
library(dplyr)
library(tidyverse)

#################################################

ANALYSIS$SEASONS <- list()

## We need Correlation of rec with winter, sommer, spring, autumn temp, prec, isot, itpc?

entity_list <- DATA_past1000$CAVES$entity_info$entity_id[mask_var]

ANALYSIS$SEASONS$Plot_Lyr_TEMP <- data.frame(
  lon = numeric(length(entity_list)),
  lat = numeric(length(entity_list)),
  value = numeric(length(entity_list)),
  season = numeric(length(entity_list))
)

ANALYSIS$SEASONS$Plot_Lyr_PREC <- data.frame(
  lon = numeric(length(entity_list)),
  lat = numeric(length(entity_list)),
  value = numeric(length(entity_list)),
  season = numeric(length(entity_list))
)

ANALYSIS$SEASONS$Plot_Lyr_ISOT <- data.frame(
  lon = numeric(length(entity_list)),
  lat = numeric(length(entity_list)),
  value = numeric(length(entity_list)),
  season = numeric(length(entity_list))
)

counter = 1
for(entity in entity_list){
  print(entity)
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  ANALYSIS$SEASONS$Plot_Lyr_PREC$lon[counter] = DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site]
  ANALYSIS$SEASONS$Plot_Lyr_PREC$lat[counter] = DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site]
  
  ANALYSIS$SEASONS$Plot_Lyr_TEMP$lon[counter] = DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site]
  ANALYSIS$SEASONS$Plot_Lyr_TEMP$lat[counter] = DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site]
  
  ANALYSIS$SEASONS$Plot_Lyr_ISOT$lon[counter] = DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site]
  ANALYSIS$SEASONS$Plot_Lyr_ISOT$lat[counter] = DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site]
  
  corr_temp = numeric(5)
  corr_temp_p = numeric(5)
  corr_prec = numeric(5)
  corr_prec_p = numeric(5)
  corr_isot = numeric(5)
  corr_isot_p = numeric(5)
  
  TS <- list()
  s <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
  # zoo cannot handle objects where order.by has two elements which is why they are sorted out here (no better option found)
  double_time <- s %>% group_by(interp_age) %>% count() %>% filter(n>1)
  s <- s %>% filter(!interp_age %in% double_time$interp_age) %>% filter(!is.na(d18O_measurement))
  
  TS[["Record"]] <- zoo(x = s$d18O_dw_eq, order.by = s$interp_age)
  TS[["SIM"]] <- zoo( x= SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_seasonal[[paste0("CAVE", site)]]$WINTER$temp_mean), 
                                                               start = -50,
                                                               end = 1100),
                                                            s$interp_age), order.by = s$interp_age)
  
  corr <- cor.test(TS$Record, TS$SIM, conflevel = 0.1)
  corr_temp[1] = corr$estimate[[1]]
  corr_temp_p[1] = corr$p.value
  
  TS[["SIM"]] <- zoo( x= SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_seasonal[[paste0("CAVE", site)]]$SPRING$temp_mean), 
                                                               start = -50,
                                                               end = 1100),
                                                            s$interp_age), order.by = s$interp_age)
  
  corr <- cor.test(TS$Record, TS$SIM, conflevel = 0.1)
  corr_temp[2] = corr$estimate[[1]]
  corr_temp_p[2] = corr$p.value
  
  TS[["SIM"]] <- zoo( x= SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_seasonal[[paste0("CAVE", site)]]$SUMMER$temp_mean), 
                                                               start = -50,
                                                               end = 1100),
                                                            s$interp_age), order.by = s$interp_age)
  
  corr <- cor.test(TS$Record, TS$SIM, conflevel = 0.1)
  corr_temp[3] = corr$estimate[[1]]
  corr_temp_p[3] = corr$p.value
  
  TS[["SIM"]] <- zoo( x= SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_seasonal[[paste0("CAVE", site)]]$AUTUMN$temp_mean), 
                                                               start = -50,
                                                               end = 1100),
                                                            s$interp_age), order.by = s$interp_age)
  
  corr <- cor.test(TS$Record, TS$SIM, conflevel = 0.1)
  corr_temp[4] = corr$estimate[[1]]
  corr_temp_p[4] = corr$p.value
  
  TS[["SIM"]] <- zoo( x= SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$TEMP), 
                                                               start = -50,
                                                               end = 1100),
                                                            s$interp_age), order.by = s$interp_age)
  
  corr <- cor.test(TS$Record, TS$SIM, conflevel = 0.1)
  corr_temp[5] = corr$estimate[[1]]
  corr_temp_p[5] = corr$p.value
  
  
  ANALYSIS$SEASONS$Plot_Lyr_TEMP$value[counter] = max(abs(corr_temp))
  if(all(is.na(corr_temp))){
    ANALYSIS$SEASONS$Plot_Lyr_TEMP$season[counter] = NA
  } else {
    ANALYSIS$SEASONS$Plot_Lyr_TEMP$season[counter] = which.max(abs(corr_temp))
  }
  
  #PREC
  TS[["SIM"]] <- zoo( x= SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_seasonal[[paste0("CAVE", site)]]$WINTER$prec_mean), 
                                                               start = -50,
                                                               end = 1100),
                                                            s$interp_age), order.by = s$interp_age)
  
  corr <- cor.test(TS$Record, TS$SIM, conflevel = 0.1)
  corr_prec[1] = corr$estimate[[1]]
  corr_prec_p[1] = corr$p.value
  
  TS[["SIM"]] <- zoo( x= SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_seasonal[[paste0("CAVE", site)]]$SPRING$prec_mean), 
                                                               start = -50,
                                                               end = 1100),
                                                            s$interp_age), order.by = s$interp_age)
  
  corr <- cor.test(TS$Record, TS$SIM, conflevel = 0.1)
  corr_prec[2] = corr$estimate[[1]]
  corr_prec_p[2] = corr$p.value
  
  TS[["SIM"]] <- zoo( x= SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_seasonal[[paste0("CAVE", site)]]$SUMMER$prec_mean), 
                                                               start = -50,
                                                               end = 1100),
                                                            s$interp_age), order.by = s$interp_age)
  
  corr <- cor.test(TS$Record, TS$SIM, conflevel = 0.1)
  corr_prec[3] = corr$estimate[[1]]
  corr_prec_p[3] = corr$p.value
  
  TS[["SIM"]] <- zoo( x= SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_seasonal[[paste0("CAVE", site)]]$AUTUMN$prec_mean), 
                                                               start = -50,
                                                               end = 1100),
                                                            s$interp_age), order.by = s$interp_age)
  
  corr <- cor.test(TS$Record, TS$SIM, conflevel = 0.1)
  corr_prec[4] = corr$estimate[[1]]
  corr_prec_p[4] = corr$p.value
  
  TS[["SIM"]] <- zoo( x= SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$PREC), 
                                                               start = -50,
                                                               end = 1100),
                                                            s$interp_age), order.by = s$interp_age)
  
  corr <- cor.test(TS$Record, TS$SIM, conflevel = 0.1)
  corr_prec[5] = corr$estimate[[1]]
  corr_prec_p[5] = corr$p.value
  
  
  ANALYSIS$SEASONS$Plot_Lyr_PREC$value[counter] = max(abs(corr_prec))
  if(all(is.na(corr_prec))){
    ANALYSIS$SEASONS$Plot_Lyr_PREC$season[counter] = NA
  } else {
    ANALYSIS$SEASONS$Plot_Lyr_PREC$season[counter] = which.max(abs(corr_prec))
  }
  
  #ISOT
  TS[["SIM"]] <- zoo( x= SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_seasonal[[paste0("CAVE", site)]]$WINTER$isot_mean), 
                                                               start = -50,
                                                               end = 1100),
                                                            s$interp_age), order.by = s$interp_age)
  
  corr <- cor.test(TS$Record, TS$SIM, conflevel = 0.1)
  corr_isot[1] = corr$estimate[[1]]
  corr_isot_p[1] = corr$p.value
  
  TS[["SIM"]] <- zoo( x= SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_seasonal[[paste0("CAVE", site)]]$SPRING$isot_mean), 
                                                               start = -50,
                                                               end = 1100),
                                                            s$interp_age), order.by = s$interp_age)
  
  corr <- cor.test(TS$Record, TS$SIM, conflevel = 0.1)
  corr_isot[2] = corr$estimate[[1]]
  corr_isot_p[2] = corr$p.value
  
  TS[["SIM"]] <- zoo( x= SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_seasonal[[paste0("CAVE", site)]]$SUMMER$isot_mean), 
                                                               start = -50,
                                                               end = 1100),
                                                            s$interp_age), order.by = s$interp_age)
  
  corr <- cor.test(TS$Record, TS$SIM, conflevel = 0.1)
  corr_isot[3] = corr$estimate[[1]]
  corr_isot_p[3] = corr$p.value
  
  TS[["SIM"]] <- zoo( x= SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_seasonal[[paste0("CAVE", site)]]$AUTUMN$isot_mean), 
                                                               start = -50,
                                                               end = 1100),
                                                            s$interp_age), order.by = s$interp_age)
  
  corr <- cor.test(TS$Record, TS$SIM, conflevel = 0.1)
  corr_isot[4] = corr$estimate[[1]]
  corr_isot_p[4] = corr$p.value
  
  TS[["SIM"]] <- zoo( x= SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ISOT), 
                                                               start = -50,
                                                               end = 1100),
                                                            s$interp_age), order.by = s$interp_age)
  
  corr <- cor.test(TS$Record, TS$SIM, conflevel = 0.1)
  corr_isot[5] = corr$estimate[[1]]
  corr_isot_p[5] = corr$p.value
  
  
  ANALYSIS$SEASONS$Plot_Lyr_ISOT$value[counter] = max(abs(corr_isot))
  if(all(is.na(corr_isot))){
    ANALYSIS$SEASONS$Plot_Lyr_ISOT$season[counter] = NA
  } else {
    ANALYSIS$SEASONS$Plot_Lyr_ISOT$season[counter] = which.max(abs(corr_isot))
  }
  counter = counter + 1
}

try_prec <- array(dim = c(96,73))


for(lon in 1:96){
  for(lat in 2:72){
    winter_corr = cor.test(DATA_past1000$SIM_seasonal$WINTER$PREC[lon,lat,], DATA_past1000$SIM_yearly$ISOT[lon,lat,])
    spring_corr = cor.test(DATA_past1000$SIM_seasonal$SPRING$PREC[lon,lat,], DATA_past1000$SIM_yearly$ISOT[lon,lat,])
    summer_corr = cor.test(DATA_past1000$SIM_seasonal$SUMMER$PREC[lon,lat,], DATA_past1000$SIM_yearly$ISOT[lon,lat,])
    autumn_corr = cor.test(DATA_past1000$SIM_seasonal$AUTUMN$PREC[lon,lat,], DATA_past1000$SIM_yearly$ISOT[lon,lat,])
    yearly_corr = cor.test(DATA_past1000$SIM_yearly$PREC[lon,lat,], DATA_past1000$SIM_yearly$ISOT[lon,lat,])
    
    corr = c(winter_corr$estimate[[1]], spring_corr$estimate[[1]], summer_corr$estimate[[1]], autumn_corr$estimate[[1]], yearly_corr$estimate[[1]])
    p.value = c(winter_corr$p.value, spring_corr$p.value, summer_corr$p.value, autumn_corr$p.value, yearly_corr$p.value)
    
    corr[p.value>0.2] = NA
    
    if(all(is.na(corr))){
      try_prec[lon,lat] = NA
    }else{
      try_prec[lon,lat] = which.max(abs(corr)) 
    }
    
  }
}


## PLOTTING

source("Functions/Plotting/STACYmap_5.R")
source("Functions/projection_ptlyr.R")

## PREC

Plot_Lyr <- data.frame(lon = ANALYSIS$SEASONS$Plot_Lyr_PREC$lon, 
                       lat = ANALYSIS$SEASONS$Plot_Lyr_PREC$lat, 
                       layer = ANALYSIS$SEASONS$Plot_Lyr_PREC$season, 
                       value = ANALYSIS$SEASONS$Plot_Lyr_PREC$value)

mask_china <- logical(length(Plot_Lyr$lon))

for(ii in 1:length(Plot_Lyr$lon)){
  if(is.na(Plot_Lyr$lon[ii])){next}
  if(Plot_Lyr$lon[ii] > 100 & Plot_Lyr$lon[ii] < 120){
    if(Plot_Lyr$lat[ii] < 35 & Plot_Lyr$lat[ii] > 22){
      mask_china[ii] = T}
  }
}

ptlyr_china <- data.frame(
  lon = Plot_Lyr$lon[mask_china],
  lat = Plot_Lyr$lat[mask_china],
  layer = Plot_Lyr$layer[mask_china],
  value = Plot_Lyr$value[mask_china]
)

ptlyr_rest <- data.frame(
  lon = Plot_Lyr$lon[!mask_china],
  lat = Plot_Lyr$lat[!mask_china],
  layer = Plot_Lyr$layer[!mask_china],
  value = Plot_Lyr$value[!mask_china]
)


ptlyr_china_p <- projection_ptlyr(ptlyr_china, as.character('+proj=robin +datum=WGS84'))#
ptlyr_rest_p <- projection_ptlyr(ptlyr_rest, as.character('+proj=robin +datum=WGS84'))
ptlyr_china_p$value <- ptlyr_china$value
ptlyr_rest_p$value <- ptlyr_rest$value

ptlyr_china_p$layer <- factor(ptlyr_china$layer)
ptlyr_rest_p$layer <- factor(ptlyr_rest$layer)


# plot <- STACYmap(gridlyr = rbind(try_prec[49:96,1:73], try_prec[1:48,1:73]), coastline = T, 
#                  colorscheme = c("#0061FF", "#30B700", "#9B1D20", "#FFB600", "#483D3F")) +
#         new_scale_fill() +
plot <- STACYmap(coastline = T) +
  geom_point(data = ptlyr_china_p, aes(x = long, y = lat, color = layer, size = value), shape = 19, 
             show.legend = c(color =TRUE, size = TRUE), position = position_jitter(width = 1000000, height = 500000)) +
  geom_point(data = ptlyr_rest_p, aes(x = long, y = lat, color = layer, size = value), shape = 19, 
             show.legend = c(color =TRUE, size = TRUE)) +
  scale_color_manual(name = "Seasons", labels = c("DJF", "MAM", "JJA", "SON", "year"), 
                     values = c("#0061FF", "#30B700", "#9B1D20", "#FFB600", "#483D3F")) +
  scale_size(name = "corr.") +
  ggtitle("Precipitation") +
  theme(panel.border = element_blank(),
        legend.background = element_blank(),
        axis.text = element_blank(),
        legend.direction = "horizontal",
        text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(nrow=2,byrow=TRUE), size = guide_legend(nrow = 2, byrow = T))

plot

plot  %>% ggsave(filename = paste('Paper_Plot_A2_Seasons_prec', 'pdf', sep = '.'), plot = ., path = 'Plots', 
                          width = 2*8.3, height =2*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")

## TEMP

Plot_Lyr <- data.frame(lon = ANALYSIS$SEASONS$Plot_Lyr_TEMP$lon, 
                       lat = ANALYSIS$SEASONS$Plot_Lyr_TEMP$lat, 
                       layer = ANALYSIS$SEASONS$Plot_Lyr_TEMP$season, 
                       value = ANALYSIS$SEASONS$Plot_Lyr_TEMP$value)

mask_china <- logical(length(Plot_Lyr$lon))

for(ii in 1:length(Plot_Lyr$lon)){
  if(is.na(Plot_Lyr$lon[ii])){next}
  if(Plot_Lyr$lon[ii] > 100 & Plot_Lyr$lon[ii] < 120){
    if(Plot_Lyr$lat[ii] < 35 & Plot_Lyr$lat[ii] > 22){
      mask_china[ii] = T}
  }
}

ptlyr_china <- data.frame(
  lon = Plot_Lyr$lon[mask_china],
  lat = Plot_Lyr$lat[mask_china],
  layer = Plot_Lyr$layer[mask_china],
  value = Plot_Lyr$value[mask_china]
)

ptlyr_rest <- data.frame(
  lon = Plot_Lyr$lon[!mask_china],
  lat = Plot_Lyr$lat[!mask_china],
  layer = Plot_Lyr$layer[!mask_china],
  value = Plot_Lyr$value[!mask_china]
)


ptlyr_china_p <- projection_ptlyr(ptlyr_china, as.character('+proj=robin +datum=WGS84'))#
ptlyr_rest_p <- projection_ptlyr(ptlyr_rest, as.character('+proj=robin +datum=WGS84'))
ptlyr_china_p$value <- ptlyr_china$value
ptlyr_rest_p$value <- ptlyr_rest$value

ptlyr_china_p$layer <- factor(ptlyr_china$layer)
ptlyr_rest_p$layer <- factor(ptlyr_rest$layer)


# plot <- STACYmap(gridlyr = rbind(try_prec[49:96,1:73], try_prec[1:48,1:73]), coastline = T, 
#                  colorscheme = c("#0061FF", "#30B700", "#9B1D20", "#FFB600", "#483D3F")) +
#         new_scale_fill() +
plot <- STACYmap(coastline = T) +
  geom_point(data = ptlyr_china_p, aes(x = long, y = lat, color = layer, size = value), shape = 19, 
             show.legend = c(color =TRUE, size = TRUE), position = position_jitter(width = 1000000, height = 500000)) +
  geom_point(data = ptlyr_rest_p, aes(x = long, y = lat, color = layer, size = value), shape = 19, 
             show.legend = c(color =TRUE, size = TRUE)) +
  scale_color_manual(name = "Seasons", labels = c("DJF", "MAM", "JJA", "SON", "year"), 
                     values = c("#0061FF", "#30B700", "#9B1D20", "#FFB600", "#483D3F")) +
  scale_size(name = "corr.") +
  ggtitle("Temperature") +
  theme(panel.border = element_blank(),
        legend.background = element_blank(),
        axis.text = element_blank(),
        legend.direction = "horizontal",
        text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(nrow=2,byrow=TRUE), size = guide_legend(nrow = 2, byrow = T))

plot

plot  %>% ggsave(filename = paste('Paper_Plot_A2_Seasons_temp', 'pdf', sep = '.'), plot = ., path = 'Plots', 
                 width = 2*8.3, height =2*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")

## ISOT

Plot_Lyr <- data.frame(lon = ANALYSIS$SEASONS$Plot_Lyr_ISOT$lon, 
                       lat = ANALYSIS$SEASONS$Plot_Lyr_ISOT$lat, 
                       layer = ANALYSIS$SEASONS$Plot_Lyr_ISOT$season, 
                       value = ANALYSIS$SEASONS$Plot_Lyr_ISOT$value)

mask_china <- logical(length(Plot_Lyr$lon))

for(ii in 1:length(Plot_Lyr$lon)){
  if(is.na(Plot_Lyr$lon[ii])){next}
  if(Plot_Lyr$lon[ii] > 100 & Plot_Lyr$lon[ii] < 120){
    if(Plot_Lyr$lat[ii] < 35 & Plot_Lyr$lat[ii] > 22){
      mask_china[ii] = T}
  }
}

ptlyr_china <- data.frame(
  lon = Plot_Lyr$lon[mask_china],
  lat = Plot_Lyr$lat[mask_china],
  layer = Plot_Lyr$layer[mask_china],
  value = Plot_Lyr$value[mask_china]
)

ptlyr_rest <- data.frame(
  lon = Plot_Lyr$lon[!mask_china],
  lat = Plot_Lyr$lat[!mask_china],
  layer = Plot_Lyr$layer[!mask_china],
  value = Plot_Lyr$value[!mask_china]
)


ptlyr_china_p <- projection_ptlyr(ptlyr_china, as.character('+proj=robin +datum=WGS84'))#
ptlyr_rest_p <- projection_ptlyr(ptlyr_rest, as.character('+proj=robin +datum=WGS84'))
ptlyr_china_p$value <- ptlyr_china$value
ptlyr_rest_p$value <- ptlyr_rest$value

ptlyr_china_p$layer <- factor(ptlyr_china$layer)
ptlyr_rest_p$layer <- factor(ptlyr_rest$layer)


# plot <- STACYmap(gridlyr = rbind(try_prec[49:96,1:73], try_prec[1:48,1:73]), coastline = T, 
#                  colorscheme = c("#0061FF", "#30B700", "#9B1D20", "#FFB600", "#483D3F")) +
#         new_scale_fill() +
plot <- STACYmap(coastline = T) +
  geom_point(data = ptlyr_china_p, aes(x = long, y = lat, color = layer, size = value), shape = 19, 
             show.legend = c(color =TRUE, size = TRUE), position = position_jitter(width = 1000000, height = 500000)) +
  geom_point(data = ptlyr_rest_p, aes(x = long, y = lat, color = layer, size = value), shape = 19, 
             show.legend = c(color =TRUE, size = TRUE)) +
  scale_color_manual(name = "Seasons", labels = c("DJF", "MAM", "JJA", "SON", "year"), 
                     values = c("#0061FF", "#30B700", "#9B1D20", "#FFB600", "#483D3F")) +
  scale_size(name = "corr.") +
  ggtitle("Isotopic composition") +
  theme(panel.border = element_blank(),
        legend.background = element_blank(),
        axis.text = element_blank(),
        legend.direction = "horizontal",
        text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(nrow=2,byrow=TRUE), size = guide_legend(nrow = 2, byrow = T))

plot

plot  %>% ggsave(filename = paste('Paper_Plot_A2_Seasons_isot', 'pdf', sep = '.'), plot = ., path = 'Plots', 
                 width = 2*8.3, height =2*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")
  

# plot <- STACYmap(ptlyr = data.frame(lon = ANALYSIS$SEASONS$Plot_Lyr_TEMP$lon, 
#                                     lat = ANALYSIS$SEASONS$Plot_Lyr_TEMP$lat, 
#                                     layer = ANALYSIS$SEASONS$Plot_Lyr_TEMP$season),
#                  colorscheme = c("#0061FF", "#30B700", "#9B1D20", "#FFB600", "#483D3F"), graticules = T) + 
#   theme(panel.border = element_blank(),
#         legend.background = element_blank(),
#         axis.text = element_blank(),
#         text = element_text(size = 12),
#         legend.title = element_text(size = 12)) 
# 
# plot  %>% ggsave(filename = paste('Paper_Plot_A2_Seasons_temp', 'pdf', sep = '.'), plot = ., path = 'Plots/Appendix', 
#                  width = 2*8.3, height = 12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")
# 
# 
# plot <- STACYmap(ptlyr = data.frame(lon = ANALYSIS$SEASONS$Plot_Lyr_ISOT$lon, 
#                                     lat = ANALYSIS$SEASONS$Plot_Lyr_ISOT$lat, 
#                                     layer = ANALYSIS$SEASONS$Plot_Lyr_ISOT$season),
#                  colorscheme = c("#0061FF", "#30B700", "#9B1D20", "#FFB600", "#483D3F"), graticules = T) + 
#   theme(panel.border = element_blank(),
#         legend.background = element_blank(),
#         axis.text = element_blank(),
#         text = element_text(size = 12),
#         legend.title = element_text(size = 12)) 
# 
# plot  %>% ggsave(filename = paste('Paper_Plot_A2_Seasons_isot', 'pdf', sep = '.'), plot = ., path = 'Plots/Appendix', 
#                  width = 2*8.3, height = 12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")
