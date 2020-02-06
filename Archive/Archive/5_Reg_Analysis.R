###################################################################################################
## REGRESSION ANALYSIS ###########################################################################
###################################################################################################

# ToDo:
# [ ] Regression Maps von Temp-Iso, Prec - Iso 
#       [ ] Hierauf kann man dann Punkte legen, wie stark die Records an den einzelnen Punkten mit der Temp oder Prec korrelieren?
# [ ] Latitudinal Plot, wie stark die Iso Werte in der Simulation und die Iso-Werte in den Records übereinstimmen
# [X] Maybe Temperature Stuff is not working because it's in °C??? --> das wars jedenfalls nicht...
# [ ]
# [ ]


library(dplyr)
library(latex2exp)
source("Functions/STACYmap_5.R")
library(PaleoSpec)

REG_ANALYSIS <- list()

## 1) Get matrix analysis from simulation #######

REG_ANALYSIS$GLOBAL_REGRESSION <- list(
  REG_TEMP_PREC_a = array(dim = c(96,73)),
  REG_TEMP_PREC_b = array(dim = c(96,73)),
  REG_TEMP_PREC_r2 = array(dim = c(96,73)),
  REG_TEMP_ISOT_a = array(dim = c(96,73)),
  REG_TEMP_ISOT_b = array(dim = c(96,73)),
  REG_TEMP_ISOT_r2 = array(dim = c(96,73)),
  REG_PREC_ISOT_a = array(dim = c(96,73)),
  REG_PREC_ISOT_b = array(dim = c(96,73)),
  REG_PREC_ISOT_r2 = array(dim = c(96,73))
)

for (lon in 1:96){
  for (lat in 1:73){
    # TEMP-PREC
    temp_data <- DATA_past1000$SIM_yearly$TEMP[lon,lat,] + 273.15
    REG = lm(DATA_past1000$SIM_yearly$PREC[lon,lat,] ~ temp_data)
    REG_ANALYSIS$GLOBAL_REGRESSION$REG_TEMP_PREC_a[lon,lat] = REG$coefficients[[2]]
    REG_ANALYSIS$GLOBAL_REGRESSION$REG_TEMP_PREC_b[lon,lat] = REG$coefficients[[1]]
    REG_ANALYSIS$GLOBAL_REGRESSION$REG_TEMP_PREC_r2[lon,lat] = summary(REG)$r.squared
    
    if(!any(is.na(DATA_past1000$SIM_yearly$ISOT[lon,lat,]))){
      
      REG = lm(DATA_past1000$SIM_yearly$ISOT[lon,lat,] ~ temp_data, na.action = na.omit)
      REG_ANALYSIS$GLOBAL_REGRESSION$REG_TEMP_ISOT_a[lon,lat] = REG$coefficients[[2]]
      REG_ANALYSIS$GLOBAL_REGRESSION$REG_TEMP_ISOT_a[lon,lat]  = REG$coefficients[[2]]
      REG_ANALYSIS$GLOBAL_REGRESSION$REG_TEMP_ISOT_b[lon,lat]  = REG$coefficients[[1]]
      REG_ANALYSIS$GLOBAL_REGRESSION$REG_TEMP_ISOT_r2[lon,lat] = summary(REG)$r.squared

      REG = lm(DATA_past1000$SIM_yearly$ISOT[lon,lat,] ~ DATA_past1000$SIM_yearly$PREC[lon,lat,], na.action = na.omit)
      REG_ANALYSIS$GLOBAL_REGRESSION$REG_PREC_ISOT_a[lon,lat]  = REG$coefficients[[2]]
      REG_ANALYSIS$GLOBAL_REGRESSION$REG_PREC_ISOT_b[lon,lat]  = REG$coefficients[[1]]
      REG_ANALYSIS$GLOBAL_REGRESSION$REG_PREC_ISOT_r2[lon,lat] = summary(REG)$r.squared
    } else {
      REG_ANALYSIS$GLOBAL_REGRESSION$REG_TEMP_ISOT_a[lon,lat]  = NA
      REG_ANALYSIS$GLOBAL_REGRESSION$REG_TEMP_ISOT_b[lon,lat]  = NA
      REG_ANALYSIS$GLOBAL_REGRESSION$REG_TEMP_ISOT_r2[lon,lat] = NA
      REG_ANALYSIS$GLOBAL_REGRESSION$REG_PREC_ISOT_a[lon,lat]  = NA
      REG_ANALYSIS$GLOBAL_REGRESSION$REG_PREC_ISOT_b[lon,lat]  = NA
      REG_ANALYSIS$GLOBAL_REGRESSION$REG_PREC_ISOT_r2[lon,lat] = NA
    }
  }
}

remove(lon,lat, REG, temp_data)

## 2) Regression of CAVE with Simualtion

#Consider all caves of last 1000 years. Caves with too little data points will fall out of the REGelation because of missing significance. 

length_cave = length(DATA_past1000$CAVES$entity_info$entity_id)

REG_ANALYSIS$CAVE_REGRESSION <- data.frame(
  entity_id = numeric(length_cave),
  REG_TEMP_a  = numeric(length_cave),
  REG_TEMP_b  = numeric(length_cave),
  REG_TEMP_r2 = numeric(length_cave),
  REG_PREC_a  = numeric(length_cave),
  REG_PREC_b  = numeric(length_cave),
  REG_PREC_r2 = numeric(length_cave),
  REG_ISOT_a  = numeric(length_cave),
  REG_ISOT_b  = numeric(length_cave),
  REG_ISOT_r2 = numeric(length_cave),
  REG_ITPC_a  = numeric(length_cave),
  REG_ITPC_b  = numeric(length_cave),
  REG_ITPC_r2 = numeric(length_cave)
)

#REG_ANALYSIS$SITE_REGRESSION <- data.frame(
#  entity_id = numeric(length_cave),
#  REG_TI = numeric(length_cave),
#  PVALUE_TI = numeric(length_cave),
#  REG_PI = numeric(length_cave),
#  PVALUE_PI = numeric(length_cave),
#  REG_TI_pw = numeric(length_cave),
#  PVALUE_TI_pw = numeric(length_cave),
#  REG_PI_pw = numeric(length_cave),
#  PVALUE_PI_pw = numeric(length_cave)
#)


for(ii in 1:length_cave){
  print(ii)
  entity <- DATA_past1000$CAVES$entity_info$entity_id[ii]
  site <- DATA_past1000$CAVES$entity_info$site_id[ii]
  REG_ANALYSIS$CAVE_REGRESSION$entity_id[ii] <- entity
  #REG_ANALYSIS$SITE_REGRESSION$entity_id[ii] <- entity
  # CAREFULL --> REGRESSION ONLY WORKS FOR EQUIDISTANT DATAPOINTS
  if(length(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age)>4){
    #### SIM WITH REREGD
    temp_data <- (DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$Temp + 273.15)
    REG <- lm(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement ~ temp_data)
    REG_ANALYSIS$CAVE_REGRESSION$REG_TEMP_a[ii] = REG$coefficients[[2]]
    REG_ANALYSIS$CAVE_REGRESSION$REG_TEMP_b[ii] = REG$coefficients[[1]]
    REG_ANALYSIS$CAVE_REGRESSION$REG_TEMP_r2[ii] = summary(REG)$r.squared
    
    REG <- lm(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement ~ DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$Prec)
    REG_ANALYSIS$CAVE_REGRESSION$REG_PREC_a[ii] = REG$coefficients[[2]]
    REG_ANALYSIS$CAVE_REGRESSION$REG_PREC_b[ii] = REG$coefficients[[1]]
    REG_ANALYSIS$CAVE_REGRESSION$REG_PREC_r2[ii] = summary(REG)$r.squared
    
    REG <- lm(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement ~ DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$d18O)
    REG_ANALYSIS$CAVE_REGRESSION$REG_ISOT_a[ii] = REG$coefficients[[2]]
    REG_ANALYSIS$CAVE_REGRESSION$REG_ISOT_b[ii] = REG$coefficients[[1]]
    REG_ANALYSIS$CAVE_REGRESSION$REG_ISOT_r2[ii] = summary(REG)$r.squared
    
    REG <- lm(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement ~ DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$d18O_pw)
    REG_ANALYSIS$CAVE_REGRESSION$REG_ITPC_a[ii] = REG$coefficients[[2]]
    REG_ANALYSIS$CAVE_REGRESSION$REG_ITPC_b[ii] = REG$coefficients[[1]]
    REG_ANALYSIS$CAVE_REGRESSION$REG_ITPC_r2[ii] = summary(REG)$r.squared
    
    
    
  }else{
    REG_ANALYSIS$CAVE_REGRESSION$REG_TEMP_a[ii]  = NA
    REG_ANALYSIS$CAVE_REGRESSION$REG_TEMP_b[ii]  = NA
    REG_ANALYSIS$CAVE_REGRESSION$REG_TEMP_r2[ii] = NA
    
    REG_ANALYSIS$CAVE_REGRESSION$REG_PREC_a[ii]  = NA
    REG_ANALYSIS$CAVE_REGRESSION$REG_PREC_b[ii]  = NA
    REG_ANALYSIS$CAVE_REGRESSION$REG_PREC_r2[ii] = NA
    
    REG_ANALYSIS$CAVE_REGRESSION$REG_ISOT_a[ii]  = NA
    REG_ANALYSIS$CAVE_REGRESSION$REG_ISOT_b[ii]  = NA
    REG_ANALYSIS$CAVE_REGRESSION$REG_ISOT_r2[ii] = NA
    
    REG_ANALYSIS$CAVE_REGRESSION$REG_ITPC_a[ii]  = NA
    REG_ANALYSIS$CAVE_REGRESSION$REG_ITPC_b[ii]  = NA
    REG_ANALYSIS$CAVE_REGRESSION$REG_ITPC_r2[ii] = NA
    
  }
}

remove(ii, entity, site, REG)

#################################################
## Plot Maps: ###################################
#################################################

source("Functions/projection_ptlyr.R")

Plot_lyr <- REG_ANALYSIS$GLOBAL_REGRESSION$REG_TEMP_ISOT_a
Plot_lyr_r2 <- REG_ANALYSIS$GLOBAL_REGRESSION$REG_TEMP_ISOT_r2
Plot_lyr[Plot_lyr_r2 < 0.3] <- NA

Plot_lyr <- rbind(Plot_lyr[49:96,1:73],
                  Plot_lyr[1:48,1:73])

Point_Lyr <- list(lon = list(), lat = list(), value = list())
Point2_Lyr <- list(lon = list(), lat = list(), value = list())

for(ii in 1:length_cave){
  site <- DATA_past1000$CAVES$entity_info$site_id[ii]
  if(is.na(REG_ANALYSIS$CAVE_REGRESSION$REG_TEMP_r2[ii])){next}
  if(REG_ANALYSIS$CAVE_REGRESSION$REG_TEMP_r2[ii] < 0.3){
    Point2_Lyr$lon = c(Point2_Lyr$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point2_Lyr$lat = c(Point2_Lyr$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point2_Lyr$value = c(Point2_Lyr$value, REG_ANALYSIS$CAVE_REGRESSION$REG_TEMP_a[ii])
  }else{
    Point_Lyr$lon = c(Point_Lyr$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point_Lyr$lat = c(Point_Lyr$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point_Lyr$value = c(Point_Lyr$value, REG_ANALYSIS$CAVE_REGRESSION$REG_TEMP_a[ii])
  }
}

Point_Lyr$lon = as.numeric(Point_Lyr$lon)
Point_Lyr$lat = as.numeric(Point_Lyr$lat)
Point_Lyr$value = as.numeric(Point_Lyr$value)

Point2_Lyr$lon = as.numeric(Point2_Lyr$lon)
Point2_Lyr$lat = as.numeric(Point2_Lyr$lat)
Point2_Lyr$value = as.numeric(Point2_Lyr$value)


GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE <- 3

Point_Lyr <-  projection_ptlyr(as.data.frame(Point_Lyr), as.character('+proj=robin +datum=WGS84'))
Point2_Lyr <-  projection_ptlyr(as.data.frame(Point2_Lyr), as.character('+proj=robin +datum=WGS84'))


plot <- STACYmap(gridlyr = Plot_lyr,
                 centercolor = 0,
                 graticules = T,
                 legend_names = list(grid = "REG", pt = "")) +
  geom_point(data = Point2_Lyr, aes(x = long, y = lat), fill = 'gray', shape = 20,
             size = (GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE-1), show.legend = c(shape =TRUE))+
  geom_point(data = Point_Lyr, aes(x = long, y = lat, fill = layer), shape = 21,
             size = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE, show.legend = c(color = TRUE))+
  ggtitle(label = "Regression (r2 > 0.3) for simulation T-d18O") +
  theme(plot.title = element_text(h = 0.5))

plot

plot %>% ggsave(filename = paste('Map_REG_TI_Records', 'pdf', sep = '.'), plot = ., path = 'Plots/Regression', 
                width = 15, height = 10, units = 'cm', dpi = 'print', device = "pdf")

## Prec-ISOT
Plot_lyr <- REG_ANALYSIS$GLOBAL_REGRESSION$REG_PREC_ISOT_a
Plot_lyr_r2 <- REG_ANALYSIS$GLOBAL_REGRESSION$REG_Prec_ISOT_r2
Plot_lyr[Plot_lyr_r2 < 0.5] <- NA

Plot_lyr <- rbind(Plot_lyr[49:96,1:73],
                  Plot_lyr[1:48,1:73])

Point_Lyr <- list(lon = list(), lat = list(), value = list())
Point2_Lyr <- list(lon = list(), lat = list(), value = list())

for(ii in 1:length_cave){
  site <- DATA_past1000$CAVES$entity_info$site_id[ii]
  if(is.na(REG_ANALYSIS$CAVE_REGRESSION$REG_PREC_r2[ii])){next}
  if(REG_ANALYSIS$CAVE_REGRESSION$REG_TEMP_r2[ii] < 0.5){
    Point2_Lyr$lon = c(Point2_Lyr$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point2_Lyr$lat = c(Point2_Lyr$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point2_Lyr$value = c(Point2_Lyr$value, REG_ANALYSIS$CAVE_REGRESSION$REG_PREC_a[ii])
  }else{
    Point_Lyr$lon = c(Point_Lyr$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point_Lyr$lat = c(Point_Lyr$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point_Lyr$value = c(Point_Lyr$value, REG_ANALYSIS$CAVE_REGRESSION$REG_PREC_a[ii])
  }
}

Point_Lyr$lon = as.numeric(Point_Lyr$lon)
Point_Lyr$lat = as.numeric(Point_Lyr$lat)
Point_Lyr$value = as.numeric(Point_Lyr$value)

Point2_Lyr$lon = as.numeric(Point2_Lyr$lon)
Point2_Lyr$lat = as.numeric(Point2_Lyr$lat)
Point2_Lyr$value = as.numeric(Point2_Lyr$value)


GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE <- 3

Point_Lyr <-  projection_ptlyr(as.data.frame(Point_Lyr), as.character('+proj=robin +datum=WGS84'))
Point2_Lyr <-  projection_ptlyr(as.data.frame(Point2_Lyr), as.character('+proj=robin +datum=WGS84'))


plot <- STACYmap(gridlyr = Plot_lyr/1e5,
                 graticules = T,
                 colorscheme = 
                 legend_names = list(grid = "REG", pt = "")) +
  geom_point(data = Point2_Lyr, aes(x = long, y = lat), fill = 'gray', shape = 20,
             size = (GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE-1), show.legend = c(shape =TRUE))+
  geom_point(data = Point_Lyr, aes(x = long, y = lat, fill = layer/1e5), shape = 21,
             size = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE, show.legend = c(color = TRUE))+
  ggtitle(label = "Regression (r2 > 0.5) for simulation P-d18O") +
  theme(plot.title = element_text(h = 0.5))

plot


plot %>% ggsave(filename = paste('Map_REG_PI_Records', 'pdf', sep = '.'), plot = ., path = 'Plots/Regression', 
                width = 15, height = 10, units = 'cm', dpi = 'print', device = "pdf")

######################################################################################

#Jetzt wäre vielleicht noch eine Karte interessant, mit der Regression zw Records und Sim ISOT

plot(REG_ANALYSIS$CAVE_REGRESSION$REG_ITPC_a[REG_ANALYSIS$CAVE_REGRESSION$REG_ISOT_r2>0.1], ylab = "", main = "ISOT (record) ~ ISOT(sim_ds) mit r2>0.1 (aus 134)")
