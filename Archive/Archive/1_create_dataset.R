#################################################
## CREATE DATASET ###############################
#################################################

## 0.1) Set Data-Structure? #####################
## 1.0) read in simulation data from xnapa ######
## 2.0) SISAL TIME SERIES #######################
## 3.0) Extract data from Caves #################
## 4.0) SIMULATION YEARLY #######################
## 5.0) Seasonal Data ###########################
## 6.0) CLIMATE INFO ############################
## 7.0) DOWNSAMPELING ###########################
## 8.0) SIMULATION SEASONS ######################
## 9.0)SIMULATION MEAN ##########################
## 10 )ENTITY in LAT band #######################
## 11 )ARAGONITE, CALCITE, SMOW #################

#################################################
#################################################

#################################################
##0.1) Set Data-Structure?#######################
#################################################

DATA_past1000 <- list()

DATA_past1000$SIM_yearly <- list(
  TEMP = list(),
  PREC = list(), 
  ISOT = list(),
  WIND = list(
    STRENGTH = list(),
    ANGLE = list()
  ),
  SLPR = list()
)

DATA_past1000$CAVES <- list()
DATA_past1000$CAVES$site_info <- read.csv("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/02_SISAL/SISAL_v2_CARLA/site_countries.csv")
DATA_past1000$CAVES$entity_info <- list()
DATA_past1000$CAVES$record_data <- list()
DATA_past1000$CAVES$sim_data_yearly <- list()
DATA_past1000$CAVES$sim_data_downsampled <- list()

# this will be an extra list, so it can be removed easily later
DATA_past1000_SIM_RAW <- list(
  TEMP = list(),
  PREC = list(), 
  ISOT = list(),
  WIND = list(
    WESTERLY = list(),
    SOUTHERLY = list()
  ),
  SLPR = list()
)

#################################################
##1) read in simulation data from xnapa #########
#################################################

source("Functions/clear_data_matrix.R")

ncf<-ncdf4::nc_open("~/Dokumente/01_Promotion/06_Daten/05_HadCM3/xnap/xnapa_surface_temperature.nc")
DATA_past1000_SIM_RAW$TEMP <- clear_data_matrix(ncdf4::ncvar_get(ncf),1)
ncdf4::nc_close(ncf)

ncf<-ncdf4::nc_open("~/Dokumente/01_Promotion/06_Daten/05_HadCM3/xnap/xnapa_precipitation.nc")
DATA_past1000_SIM_RAW$PREC <- clear_data_matrix(ncdf4::ncvar_get(ncf),2)
DATA_past1000_SIM_RAW$lon <- ncf$dim$longitude$vals
DATA_past1000_SIM_RAW$lat <- ncf$dim$latitude$vals
ncdf4::nc_close(ncf)

ncf<-ncdf4::nc_open("~/Dokumente/01_Promotion/06_Daten/05_HadCM3/xnap/xnapa_isotopes.nc")
DATA_past1000_SIM_RAW$ISOT <- clear_data_matrix(ncdf4::ncvar_get(ncf, ncf$var[[3]]),3)
ncdf4::nc_close(ncf)

ncf<-ncdf4::nc_open("~/Dokumente/01_Promotion/06_Daten/05_HadCM3/xnap/xnapa_westerly_u.nc")
DATA_past1000_SIM_RAW$WIND$WESTERLY <- clear_data_matrix(ncdf4::ncvar_get(ncf),4)
ncdf4::nc_close(ncf)
ncf<-ncdf4::nc_open("~/Dokumente/01_Promotion/06_Daten/05_HadCM3/xnap/xnapa_southerly_v.nc")
DATA_past1000_SIM_RAW$WIND$SOUTHERLY <- clear_data_matrix(ncdf4::ncvar_get(ncf),4)
ncdf4::nc_close(ncf)

ncf<-ncdf4::nc_open("~/Dokumente/01_Promotion/06_Daten/05_HadCM3/xnap/xnapa_sea_level_pressure.nc")
DATA_past1000_SIM_RAW$SLPR <- ncdf4::ncvar_get(ncf)
ncdf4::nc_close(ncf)

ncf<-ncdf4::nc_open("~/Dokumente/01_Promotion/06_Daten/05_HadCM3/hadcm_oro_001cyr.nc")
DATA_past1000$SIM_mean$OROG <- ncdf4::ncvar_get(ncf)
ncdf4::nc_close(ncf)

remove(ncf)

#################################################
##2) SISAL TIME SERIES ##########################
#################################################

# needs to be imported first, as then only the relevant cave sites will be extracted and calculated further

source("Functions/fun_with_SISALv2_Janica.R")

data <- load_sisal_data_janica()

DATA_past1000$CAVES$entity_info <- data[[1]]
DATA_past1000$CAVES$entity_dating <- data[[3]]
DATA_past1000$CAVES$site_to_entity <- data[[4]]

#Schmeißt alle Höhlen raus, die nicht gebraucht werden
DATA_past1000$CAVES$site_info <- DATA_past1000$CAVES$site_info %>% filter(site_id %in% DATA_past1000$CAVES$entity_info$site_id)

for (ii in DATA_past1000$CAVES$entity_info$entity_id){
  name = paste0("ENTITY", ii)
  if(ii%%10 == 0){
    print(name)
  }
  
  site <- DATA_past1000$CAVES$entity_info %>% filter(entity_id == ii) %>% distinct(site_id)
  DATA_past1000$CAVES$record_data[[name]] <- data[[2]] %>% filter(entity_id == ii) %>% distinct(entity_id, mineralogy, arag_corr, interp_age, d18O_measurement) %>%
    mutate(site_id = (site$site_id))
}

remove(data, site, ii, name)

#################################################
## 3) Extract data from Caves such that they are in a grid box that is the average of all surrounding

source("Functions/extract_gridboxes.R")

for (ii in 1:(dim(DATA_past1000$CAVES$site_info)[1])){
  lon_cave = DATA_past1000$CAVES$site_info$longitude[ii]
  if(lon_cave<0){lon_cave = 360+lon_cave}
  lat_cave = DATA_past1000$CAVES$site_info$latitude[ii]
  site_id = DATA_past1000$CAVES$site_info$site_id[ii]
  
  ratios <- extract_gridboxes(lon_cave, lat_cave)
  
  name <- paste0("CAVE",site_id)
  DATA_past1000$CAVES$sim_data_raw[[name]]$TEMP <- ratios$E1*DATA_past1000_SIM_RAW$TEMP[ratios$E1_lon_pos, ratios$E1_lat_pos,] +
    ratios$E2*DATA_past1000_SIM_RAW$TEMP[ratios$E2_lon_pos, ratios$E2_lat_pos,] +
    ratios$E3*DATA_past1000_SIM_RAW$TEMP[ratios$E3_lon_pos, ratios$E3_lat_pos,] +
    ratios$E4*DATA_past1000_SIM_RAW$TEMP[ratios$E4_lon_pos, ratios$E4_lat_pos,]

  DATA_past1000$CAVES$sim_data_raw[[name]]$PREC <- ratios$E1*DATA_past1000_SIM_RAW$PREC[ratios$E1_lon_pos, ratios$E1_lat_pos,] +
    ratios$E2*DATA_past1000_SIM_RAW$PREC[ratios$E2_lon_pos, ratios$E2_lat_pos,] +
    ratios$E3*DATA_past1000_SIM_RAW$PREC[ratios$E3_lon_pos, ratios$E3_lat_pos,] +
    ratios$E4*DATA_past1000_SIM_RAW$PREC[ratios$E4_lon_pos, ratios$E4_lat_pos,]

  DATA_past1000$CAVES$sim_data_raw[[name]]$ISOT <- ratios$E1*DATA_past1000_SIM_RAW$ISOT[ratios$E1_lon_pos, ratios$E1_lat_pos,] +
    ratios$E2*DATA_past1000_SIM_RAW$ISOT[ratios$E2_lon_pos, ratios$E2_lat_pos,] +
    ratios$E3*DATA_past1000_SIM_RAW$ISOT[ratios$E3_lon_pos, ratios$E3_lat_pos,] +
    ratios$E4*DATA_past1000_SIM_RAW$ISOT[ratios$E4_lon_pos, ratios$E4_lat_pos,]

  DATA_past1000$CAVES$sim_data_raw[[name]]$SLPR <- ratios$E1*DATA_past1000_SIM_RAW$SLPR[ratios$E1_lon_pos, ratios$E1_lat_pos,] +
    ratios$E2*DATA_past1000_SIM_RAW$SLPR[ratios$E2_lon_pos, ratios$E2_lat_pos,] +
    ratios$E3*DATA_past1000_SIM_RAW$SLPR[ratios$E3_lon_pos, ratios$E3_lat_pos,] +
    ratios$E4*DATA_past1000_SIM_RAW$SLPR[ratios$E4_lon_pos, ratios$E4_lat_pos,]
  
  
}

elevation_cave_sim <- array(dim = c(length(DATA_past1000$CAVES$entity_info$entity_id), 5))
colnames(elevation_cave_sim) <- c("site_id", "entity_id", "elevation_sim", "elevation_cave", "sim-cave")

for (ii in 1:(length(DATA_past1000$CAVES$entity_info$entity_id))){
  site_id = DATA_past1000$CAVES$entity_info$site_id[ii]
  lon_cave = DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site_id]
  if(lon_cave<0){lon_cave = 360+lon_cave}
  lat_cave = DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site_id]

  ratios <- extract_gridboxes(lon_cave, lat_cave)
  
  elevation_cave_sim[ii,1] <- site_id
  elevation_cave_sim[ii,2] <- DATA_past1000$CAVES$entity_info$entity_id[ii]
  elevation_cave_sim[ii,3] <- ratios$E1*DATA_past1000$SIM_mean$OROG[ratios$E1_lon_pos, ratios$E1_lat_pos] +
    ratios$E2*DATA_past1000$SIM_mean$OROG[ratios$E2_lon_pos, ratios$E2_lat_pos] +
    ratios$E3*DATA_past1000$SIM_mean$OROG[ratios$E3_lon_pos, ratios$E3_lat_pos] +
    ratios$E4*DATA_past1000$SIM_mean$OROG[ratios$E4_lon_pos, ratios$E4_lat_pos]
  
  elevation_cave_sim[ii,4] <- DATA_past1000$CAVES$site_info$elevation[DATA_past1000$CAVES$site_info$site_id == site_id]
  
  elevation_cave_sim[ii,5] <- elevation_cave_sim[ii,3] - elevation_cave_sim[ii,4]

  
}

DATA_past1000$CAVES$elevation_cave_sim <- as.data.frame(elevation_cave_sim)

remove(ratios, ii, lat_cave, lon_cave, name, site_id, elevation_cave_sim)

#################################################
## 4.0) SIMULATION YEARLY #######################
#################################################

DATA_past1000$SIM_yearly$TEMP <- array(dim = c(96,73,1150))
DATA_past1000$SIM_yearly$PREC <- array(dim = c(96,73,1150))
DATA_past1000$SIM_yearly$ISOT <- array(dim = c(96,73,1150))
DATA_past1000$SIM_yearly$ITPC <- array(dim = c(96,73,1150)) # <- prec weighted mean
DATA_past1000$SIM_yearly$WIND$STRENGTH <- array(dim = c(96,72,1149))
DATA_past1000$SIM_yearly$WIND$ANGLE <- array(dim = c(96,72,1149))
DATA_past1000$SIM_yearly$SLPR <- array(dim = c(96,73,1150))

for (lon in (1:96)){
  for (lat in 1:73){
    for(year in 1:1150){
      pos_start = 12*(year-1)+1
      pos_stop  = 12*(year-1)+12
      DATA_past1000$SIM_yearly$TEMP[lon,lat,year] = mean(DATA_past1000_SIM_RAW$TEMP[lon,lat, pos_start:pos_stop], na.rm = T)
      DATA_past1000$SIM_yearly$PREC[lon,lat,year] = mean(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop], na.rm = T)
      DATA_past1000$SIM_yearly$ISOT[lon,lat,year] = mean(DATA_past1000_SIM_RAW$ISOT[lon,lat, pos_start:pos_stop], na.rm = T)
      if(year<1150){DATA_past1000$SIM_yearly$SLPR[lon,lat,year] = mean(DATA_past1000_SIM_RAW$SLPR[lon,lat, pos_start:pos_stop], na.rm = T)}
      else{DATA_past1000$SIM_yearly$SLPR[lon,lat,year] = DATA_past1000$SIM_yearly$SLPR[lon,lat,year-1]}
      DATA_past1000$SIM_yearly$ITPC[lon,lat,year] = sum(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*DATA_past1000_SIM_RAW$ISOT[lon,lat, pos_start:pos_stop],
                                                       na.rm = T)/sum(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop], na.rm = T)
      if(lat < 73 & year < 1150){
        DATA_past1000$SIM_yearly$WIND$STRENGTH[lon,lat,year] = mean(sqrt(DATA_past1000_SIM_RAW$WIND$WESTERLY[lon,lat, pos_start:pos_stop]^2 +
                                                                         DATA_past1000_SIM_RAW$WIND$SOUTHERLY[lon,lat, pos_start:pos_stop]^2), na.rm = T)
        DATA_past1000$SIM_yearly$WIND$ANGLE[lon,lat,year] = mean(atan2(DATA_past1000_SIM_RAW$WIND$SOUTHERLY[lon,lat,pos_start:pos_stop],
                                                                       DATA_past1000_SIM_RAW$WIND$WESTERLY[lon,lat,pos_start:pos_stop]), na.rm = T)
      }
    }
  }
}

remove(lon,lat,year, pos_start, pos_stop, site_id)

for (ii in 1:(dim(DATA_past1000$CAVES$site_info)[1])){
  site_id = DATA_past1000$CAVES$site_info$site_id[ii]
  name <- paste0("CAVE", site_id)
  DATA_past1000$CAVES$sim_data_yearly[[name]]$TEMP <- numeric(1150)
  DATA_past1000$CAVES$sim_data_yearly[[name]]$PREC <- numeric(1150)
  DATA_past1000$CAVES$sim_data_yearly[[name]]$ISOT <- numeric(1150)
  DATA_past1000$CAVES$sim_data_yearly[[name]]$ITPC <- numeric(1150)
  DATA_past1000$CAVES$sim_data_yearly[[name]]$SLPR <- numeric(1150)
  for(year in 1:1150){
    pos_start = 12*(year-1)+1
    pos_stop  = 12*(year-1)+12
    DATA_past1000$CAVES$sim_data_yearly[[name]]$TEMP[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]]$TEMP[pos_start:pos_stop], na.rm = T)
    DATA_past1000$CAVES$sim_data_yearly[[name]]$PREC[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]]$PREC[pos_start:pos_stop], na.rm = T)
    DATA_past1000$CAVES$sim_data_yearly[[name]]$ISOT[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]]$ISOT[pos_start:pos_stop], na.rm = T)
    DATA_past1000$CAVES$sim_data_yearly[[name]]$SLPR[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]]$SLPR[pos_start:pos_stop], na.rm = T)
    DATA_past1000$CAVES$sim_data_yearly[[name]]$ITPC[year] <- sum(DATA_past1000$CAVES$sim_data_raw[[name]]$ISOT[pos_start:pos_stop]*
                                                                    DATA_past1000$CAVES$sim_data_raw[[name]]$PREC[pos_start:pos_stop], na.rm = T)/
      sum(DATA_past1000$CAVES$sim_data_raw[[name]]$PREC[pos_start:pos_stop], na.rm = T)
  }
}

remove(site_id, year, pos_start, pos_stop, name, ii)

#################################################
## 5.0) Seasonal Data ###########################
#################################################

winter_mask = c(1,1,1,NA,NA,NA,NA,NA,NA,NA,NA,NA)
spring_mask = c(NA,NA,NA,1,1,1,NA,NA,NA,NA,NA,NA)
summer_mask = c(NA,NA,NA,NA,NA,NA,1,1,1,NA,NA,NA)
autumn_mask = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,1,1,1)

#DATA_past1000$CAVES$sim_data_seasonal <- vector(mode = "list")


for (ii in 1:(dim(DATA_past1000$CAVES$site_info)[1])){
  site_id = DATA_past1000$CAVES$site_info$site_id[ii]
  name = paste0("CAVE", site_id)
  
  DATA_past1000$CAVES$sim_data_seasonal[[name]] = list(
    SUMMER = list(temp_mean = numeric(1150),
                  prec_mean = numeric(1150),
                  isot_mean = numeric(1150),
                  slpr_mean = numeric(1150)),
    AUTUMN = list(temp_mean = numeric(1150),
                  prec_mean = numeric(1150),
                  isot_mean = numeric(1150),
                  slpr_mean = numeric(1150)),
    WINTER = list(temp_mean = numeric(1150),
                  prec_mean = numeric(1150),
                  isot_mean = numeric(1150),
                  slpr_mean = numeric(1150)),
    SPRING = list(temp_mean = numeric(1150),
                  prec_mean = numeric(1150),
                  isot_mean = numeric(1150),
                  slpr_mean = numeric(1150))
  )


  for(year in 1:1150){
    pos_start = 12*(year-1)+1
    pos_stop  = 12*(year-1)+12
    
    DATA_past1000$CAVES$sim_data_seasonal[[name]]$SUMMER$temp_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]]$TEMP[pos_start:pos_stop]*summer_mask, na.rm = T)
    DATA_past1000$CAVES$sim_data_seasonal[[name]]$SUMMER$prec_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]]$PREC[pos_start:pos_stop]*summer_mask, na.rm = T)
    DATA_past1000$CAVES$sim_data_seasonal[[name]]$SUMMER$isot_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]]$ISOT[pos_start:pos_stop]*summer_mask, na.rm = T)

    DATA_past1000$CAVES$sim_data_seasonal[[name]]$AUTUMN$temp_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]]$TEMP[pos_start:pos_stop]*autumn_mask, na.rm = T)
    DATA_past1000$CAVES$sim_data_seasonal[[name]]$AUTUMN$prec_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]]$PREC[pos_start:pos_stop]*autumn_mask, na.rm = T)
    DATA_past1000$CAVES$sim_data_seasonal[[name]]$AUTUMN$isot_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]]$ISOT[pos_start:pos_stop]*autumn_mask, na.rm = T)

    DATA_past1000$CAVES$sim_data_seasonal[[name]]$WINTER$temp_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]]$TEMP[pos_start:pos_stop]*winter_mask, na.rm = T)
    DATA_past1000$CAVES$sim_data_seasonal[[name]]$WINTER$prec_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]]$PREC[pos_start:pos_stop]*winter_mask, na.rm = T)
    DATA_past1000$CAVES$sim_data_seasonal[[name]]$WINTER$isot_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]]$ISOT[pos_start:pos_stop]*winter_mask, na.rm = T)

    DATA_past1000$CAVES$sim_data_seasonal[[name]]$SPRING$temp_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]]$TEMP[pos_start:pos_stop]*spring_mask, na.rm = T)
    DATA_past1000$CAVES$sim_data_seasonal[[name]]$SPRING$prec_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]]$PREC[pos_start:pos_stop]*spring_mask, na.rm = T)
    DATA_past1000$CAVES$sim_data_seasonal[[name]]$SPRING$isot_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]]$ISOT[pos_start:pos_stop]*spring_mask, na.rm = T)
    
    if(year<1150){
      DATA_past1000$CAVES$sim_data_seasonal[[name]]$WINTER$slpr_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]]$SLPR[pos_start:pos_stop]*winter_mask, na.rm = T)
      DATA_past1000$CAVES$sim_data_seasonal[[name]]$AUTUMN$slpr_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]]$SLPR[pos_start:pos_stop]*autumn_mask, na.rm = T)
      DATA_past1000$CAVES$sim_data_seasonal[[name]]$SUMMER$slpr_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]]$SLPR[pos_start:pos_stop]*summer_mask, na.rm = T)
      DATA_past1000$CAVES$sim_data_seasonal[[name]]$SPRING$slpr_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]]$SLPR[pos_start:pos_stop]*spring_mask, na.rm = T)
    } else {
      DATA_past1000$CAVES$sim_data_seasonal[[name]]$WINTER$slpr_mean[year] <- DATA_past1000$CAVES$sim_data_seasonal[[name]]$WINTER$slpr_mean[year-1]
      DATA_past1000$CAVES$sim_data_seasonal[[name]]$AUTUMN$slpr_mean[year] <- DATA_past1000$CAVES$sim_data_seasonal[[name]]$AUTUMN$slpr_mean[year-1]
      DATA_past1000$CAVES$sim_data_seasonal[[name]]$SUMMER$slpr_mean[year] <- DATA_past1000$CAVES$sim_data_seasonal[[name]]$SUMMER$slpr_mean[year-1]
      DATA_past1000$CAVES$sim_data_seasonal[[name]]$SPRING$slpr_mean[year] <- DATA_past1000$CAVES$sim_data_seasonal[[name]]$SPRING$slpr_mean[year-1]
    }
  }
}

remove(ii, name, pos_start, pos_stop, year, site_id)

#################################################
## 6.0) CLIMATE INFO ############################
#################################################

DATA_past1000$CAVES$climate_info <- list(
  mean_temp = numeric((dim(DATA_past1000$CAVES$site_info)[1])), max_temp = numeric((dim(DATA_past1000$CAVES$site_info)[1])), min_temp = numeric((dim(DATA_past1000$CAVES$site_info)[1])),
  mean_prec = numeric((dim(DATA_past1000$CAVES$site_info)[1])), max_prec = numeric((dim(DATA_past1000$CAVES$site_info)[1])), min_prec = numeric((dim(DATA_past1000$CAVES$site_info)[1])),
  mean_isot = numeric((dim(DATA_past1000$CAVES$site_info)[1])), max_isot = numeric((dim(DATA_past1000$CAVES$site_info)[1])), min_isot = numeric((dim(DATA_past1000$CAVES$site_info)[1])),
  mean_slpr = numeric((dim(DATA_past1000$CAVES$site_info)[1])), max_slpr = numeric((dim(DATA_past1000$CAVES$site_info)[1])), min_slpr = numeric((dim(DATA_past1000$CAVES$site_info)[1]))
)

for (ii in 1:(dim(DATA_past1000$CAVES$site_info)[1])){
  site_id = DATA_past1000$CAVES$site_info$site_id[ii]
  name = paste0("CAVE", site_id)
  DATA_past1000$CAVES$climate_info$mean_temp[ii] = mean(DATA_past1000$CAVES$sim_data_raw[[name]]$TEMP, na.rm = T)
  DATA_past1000$CAVES$climate_info$max_temp[ii] = max(DATA_past1000$CAVES$sim_data_raw[[name]]$TEMP, na.rm = T)
  DATA_past1000$CAVES$climate_info$min_temp[ii] = min(DATA_past1000$CAVES$sim_data_raw[[name]]$TEMP, na.rm = T)
  
  DATA_past1000$CAVES$climate_info$mean_prec[ii] = mean(DATA_past1000$CAVES$sim_data_raw[[name]]$PREC, na.rm = T)
  DATA_past1000$CAVES$climate_info$max_prec[ii] = max(DATA_past1000$CAVES$sim_data_raw[[name]]$PREC, na.rm = T)
  DATA_past1000$CAVES$climate_info$min_prec[ii] = min(DATA_past1000$CAVES$sim_data_raw[[name]]$PREC, na.rm = T)
  
  DATA_past1000$CAVES$climate_info$mean_isot[ii] = mean(DATA_past1000$CAVES$sim_data_raw[[name]]$ISOT, na.rm = T)
  DATA_past1000$CAVES$climate_info$max_isot[ii] = max(DATA_past1000$CAVES$sim_data_raw[[name]]$ISOT, na.rm = T)
  DATA_past1000$CAVES$climate_info$min_isot[ii] = min(DATA_past1000$CAVES$sim_data_raw[[name]]$ISOT, na.rm = T)
  
  DATA_past1000$CAVES$climate_info$mean_slpr[ii] = mean(DATA_past1000$CAVES$sim_data_raw[[name]]$SLPR, na.rm = T)
  DATA_past1000$CAVES$climate_info$max_slpr[ii] = max(DATA_past1000$CAVES$sim_data_raw[[name]]$SLPR, na.rm = T)
  DATA_past1000$CAVES$climate_info$min_slpr[ii] = min(DATA_past1000$CAVES$sim_data_raw[[name]]$SLPR, na.rm = T)
}

remove(ii, name, site_id)

#################################################
## 7.0) DOWNSAMPELING ###########################
#################################################

source("Functions/SubsampleTimeseriesBlock_highresNA.R")

for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
  nameE = paste0("ENTITY", DATA_past1000$CAVES$entity_info$entity_id[ii])
  nameC = paste0("CAVE", DATA_past1000$CAVES$entity_info$site_id[ii])

  data_temp <- SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_yearly[[nameC]]$TEMP), 
                                                     start = -50, 
                                                     end = 1100),
                                                  DATA_past1000$CAVES$record_data[[nameE]]$interp_age)
  data_prec <- SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_yearly[[nameC]]$PREC), 
                                                     start = -50, 
                                                     end = 1100),
                                                  DATA_past1000$CAVES$record_data[[nameE]]$interp_age)
    
  data_isot <- SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_yearly[[nameC]]$ISOT), 
                                                start = -50, 
                                                end = 1100),
                                             DATA_past1000$CAVES$record_data[[nameE]]$interp_age)
  data_itpc <- SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_yearly[[nameC]]$ITPC), 
                                                     start = -50, 
                                                     end = 1100),
                                                  DATA_past1000$CAVES$record_data[[nameE]]$interp_age)
  data_slpr <- SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_yearly[[nameC]]$SLPR), 
                                                     start = -50, 
                                                     end = 1100),
                                                  DATA_past1000$CAVES$record_data[[nameE]]$interp_age)
  
  data <- matrix(c(DATA_past1000$CAVES$record_data[[nameE]]$interp_age, data_temp, data_prec, data_isot, data_itpc, data_slpr), ncol= 6)
  colnames(data) = c("interp_age", "Temp", "Prec","d18O", "d18O_pw", "SLPR")
  
  DATA_past1000$CAVES$sim_data_downsampled[[nameE]] <- as.tibble(data)

}

remove(nameE, nameC, ii, data_isot, data_itpc, data, data_slpr)

#################################################
## 8.0) SIMULATION SEASONS ######################
#################################################

DATA_past1000$SIM_seasonal <- list(
  SUMMER = list(
    TEMP = array(dim = c(96,73,1150)),
    PREC = array(dim = c(96,73,1150)),
    ISOT = array(dim = c(96,73,1150)),
    ITPC = array(dim = c(96,73,1150)),
    SLPR = array(dim = c(96,73,1150)),
    WIND = list(
      STRENGTH = array(dim = c(96,72,1149)),
      ANGLE = array(dim = c(96,72,1149))
    )
  ),
  AUTUMN = list(
    TEMP = array(dim = c(96,73,1150)),
    PREC = array(dim = c(96,73,1150)),
    ISOT = array(dim = c(96,73,1150)),
    ITPC = array(dim = c(96,73,1150)),
    SLPR = array(dim = c(96,73,1150)),
    WIND = list(
      STRENGTH = array(dim = c(96,72,1149)),
      ANGLE = array(dim = c(96,72,1149))
    )
  ),
  WINTER = list(
    TEMP = array(dim = c(96,73,1150)),
    PREC = array(dim = c(96,73,1150)),
    ISOT = array(dim = c(96,73,1150)),
    ITPC = array(dim = c(96,73,1150)),
    SLPR = array(dim = c(96,73,1150)),
    WIND = list(
      STRENGTH = array(dim = c(96,72,1149)),
      ANGLE = array(dim = c(96,72,1149))
    )
  ),
  SPRING = list(
    TEMP = array(dim = c(96,73,1150)),
    PREC = array(dim = c(96,73,1150)),
    ISOT = array(dim = c(96,73,1150)),
    ITPC = array(dim = c(96,73,1150)),
    SLPR = array(dim = c(96,73,1150)),
    WIND = list(
      STRENGTH = array(dim = c(96,72,1149)),
      ANGLE = array(dim = c(96,72,1149))
    )
  )
)

#################################################
##Calc

for(lon in 1:96){
  for(lat in 1:73){
    print(paste(lon,lat,sep = " "))
    for(year in 1:1150){
      #SUMMER
      pos_start = 12*(year-1)+1
      pos_stop  = 12*(year-1)+12
      
      DATA_past1000$SIM_seasonal$SUMMER$TEMP[lon,lat,year] = mean(DATA_past1000_SIM_RAW$TEMP[lon,lat, pos_start:pos_stop]*summer_mask, na.rm = T)
      DATA_past1000$SIM_seasonal$SUMMER$PREC[lon,lat,year] = mean(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*summer_mask, na.rm = T)
      DATA_past1000$SIM_seasonal$SUMMER$ISOT[lon,lat,year] = mean(DATA_past1000_SIM_RAW$ISOT[lon,lat, pos_start:pos_stop]*summer_mask, na.rm = T)
      DATA_past1000$SIM_seasonal$SUMMER$ITPC[lon,lat,year] = sum(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*DATA_past1000_SIM_RAW$ISOT[lon,lat, pos_start:pos_stop]*summer_mask,
                                                        na.rm = T)/sum(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*summer_mask, na.rm = T)
      if(lat < 73 & year < 1150){
        DATA_past1000$SIM_seasonal$SUMMER$WIND$STRENGTH[lon,lat,year] = mean(sqrt((DATA_past1000_SIM_RAW$WIND$WESTERLY[lon,lat, pos_start:pos_stop]*summer_mask)^2 +
                                                                           (DATA_past1000_SIM_RAW$WIND$SOUTHERLY[lon,lat, pos_start:pos_stop]*summer_mask)^2), na.rm = T)
        DATA_past1000$SIM_seasonal$SUMMER$WIND$ANGLE[lon,lat,year] = mean(atan2(DATA_past1000_SIM_RAW$WIND$SOUTHERLY[lon,lat,pos_start:pos_stop]*summer_mask,
                                                                       DATA_past1000_SIM_RAW$WIND$WESTERLY[lon,lat,pos_start:pos_stop]*summer_mask), na.rm = T)
      }
      if(year<1150){DATA_past1000$SIM_seasonal$SUMMER$SLPR[lon,lat,year] = mean(DATA_past1000_SIM_RAW$SLPR[lon,lat, pos_start:pos_stop]*summer_mask, na.rm = T)}
      else{DATA_past1000$SIM_seasonal$SUMMER$SLPR[lon,lat,year] = DATA_past1000$SIM_seasonal$SUMMER$SLPR[lon,lat,year-1]}
      
      #AUTUMN

      DATA_past1000$SIM_seasonal$AUTUMN$TEMP[lon,lat,year] = mean(DATA_past1000_SIM_RAW$TEMP[lon,lat, pos_start:pos_stop]*autumn_mask, na.rm = T)
      DATA_past1000$SIM_seasonal$AUTUMN$PREC[lon,lat,year] = mean(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*autumn_mask, na.rm = T)
      DATA_past1000$SIM_seasonal$AUTUMN$ISOT[lon,lat,year] = mean(DATA_past1000_SIM_RAW$ISOT[lon,lat, pos_start:pos_stop]*autumn_mask, na.rm = T)
      DATA_past1000$SIM_seasonal$AUTUMN$ITPC[lon,lat,year] = sum(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*DATA_past1000_SIM_RAW$ISOT[lon,lat, pos_start:pos_stop]*autumn_mask,
                                                                 na.rm = T)/sum(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*autumn_mask, na.rm = T)
      if(lat < 73 & year < 1150){
        DATA_past1000$SIM_seasonal$AUTUMN$WIND$STRENGTH[lon,lat,year] = mean(sqrt((DATA_past1000_SIM_RAW$WIND$WESTERLY[lon,lat, pos_start:pos_stop]*autumn_mask)^2 +
                                                                                    (DATA_past1000_SIM_RAW$WIND$SOUTHERLY[lon,lat, pos_start:pos_stop]*autumn_mask)^2), na.rm = T)
        DATA_past1000$SIM_seasonal$AUTUMN$WIND$ANGLE[lon,lat,year] = mean(atan2(DATA_past1000_SIM_RAW$WIND$SOUTHERLY[lon,lat,pos_start:pos_stop]*autumn_mask,
                                                                                DATA_past1000_SIM_RAW$WIND$WESTERLY[lon,lat,pos_start:pos_stop]*autumn_mask), na.rm = T)
      }
      if(year<1150){DATA_past1000$SIM_seasonal$AUTUMN$SLPR[lon,lat,year] = mean(DATA_past1000_SIM_RAW$SLPR[lon,lat, pos_start:pos_stop]*autumn_mask, na.rm = T)}
      else{DATA_past1000$SIM_seasonal$AUTUMN$SLPR[lon,lat,year] = DATA_past1000$SIM_seasonal$AUTUMN$SLPR[lon,lat,year-1]}
      
      #WINTER

      DATA_past1000$SIM_seasonal$WINTER$TEMP[lon,lat,year] = mean(DATA_past1000_SIM_RAW$TEMP[lon,lat, pos_start:pos_stop]*winter_mask, na.rm = T)
      DATA_past1000$SIM_seasonal$WINTER$PREC[lon,lat,year] = mean(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*winter_mask, na.rm = T)
      DATA_past1000$SIM_seasonal$WINTER$ISOT[lon,lat,year] = mean(DATA_past1000_SIM_RAW$ISOT[lon,lat, pos_start:pos_stop]*winter_mask, na.rm = T)
      DATA_past1000$SIM_seasonal$WINTER$ITPC[lon,lat,year] = sum(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*DATA_past1000_SIM_RAW$ISOT[lon,lat, pos_start:pos_stop]*winter_mask,
                                                                 na.rm = T)/sum(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*winter_mask, na.rm = T)
      if(lat < 73 & year < 1150){
        DATA_past1000$SIM_seasonal$WINTER$WIND$STRENGTH[lon,lat,year] = mean(sqrt((DATA_past1000_SIM_RAW$WIND$WESTERLY[lon,lat, pos_start:pos_stop]*winter_mask)^2 +
                                                                                    (DATA_past1000_SIM_RAW$WIND$SOUTHERLY[lon,lat, pos_start:pos_stop]*winter_mask)^2), na.rm = T)
        DATA_past1000$SIM_seasonal$WINTER$WIND$ANGLE[lon,lat,year] = mean(atan2(DATA_past1000_SIM_RAW$WIND$SOUTHERLY[lon,lat,pos_start:pos_stop]*winter_mask,
                                                                                DATA_past1000_SIM_RAW$WIND$WESTERLY[lon,lat,pos_start:pos_stop]*winter_mask), na.rm = T)
      }
      
      if(year<1150){DATA_past1000$SIM_seasonal$WINTER$SLPR[lon,lat,year] = mean(DATA_past1000_SIM_RAW$SLPR[lon,lat, pos_start:pos_stop]*winter_mask, na.rm = T)}
      else{DATA_past1000$SIM_seasonal$WINTER$SLPR[lon,lat,year] = DATA_past1000$SIM_seasonal$WINTER$SLPR[lon,lat,year-1]}
      
      #SPRING

      DATA_past1000$SIM_seasonal$SPRING$TEMP[lon,lat,year] = mean(DATA_past1000_SIM_RAW$TEMP[lon,lat, pos_start:pos_stop]*spring_mask, na.rm = T)
      DATA_past1000$SIM_seasonal$SPRING$PREC[lon,lat,year] = mean(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*spring_mask, na.rm = T)
      DATA_past1000$SIM_seasonal$SPRING$ISOT[lon,lat,year] = mean(DATA_past1000_SIM_RAW$ISOT[lon,lat, pos_start:pos_stop]*spring_mask, na.rm = T)
      DATA_past1000$SIM_seasonal$SPRING$SLPR[lon,lat,year] = mean(DATA_past1000_SIM_RAW$SLPR[lon,lat, pos_start:pos_stop]*spring_mask, na.rm = T)
      DATA_past1000$SIM_seasonal$SPRING$ITPC[lon,lat,year] = sum(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*DATA_past1000_SIM_RAW$ISOT[lon,lat, pos_start:pos_stop]*spring_mask,
                                                                 na.rm = T)/sum(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*spring_mask, na.rm = T)
      if(lat < 73 & year < 1150){
        DATA_past1000$SIM_seasonal$SPRING$WIND$STRENGTH[lon,lat,year] = mean(sqrt((DATA_past1000_SIM_RAW$WIND$WESTERLY[lon,lat, pos_start:pos_stop]*spring_mask)^2 +
                                                                                    (DATA_past1000_SIM_RAW$WIND$SOUTHERLY[lon,lat, pos_start:pos_stop]*spring_mask)^2), na.rm = T)
        DATA_past1000$SIM_seasonal$SPRING$WIND$ANGLE[lon,lat,year] = mean(atan2(DATA_past1000_SIM_RAW$WIND$SOUTHERLY[lon,lat,pos_start:pos_stop]*spring_mask,
                                                                                DATA_past1000_SIM_RAW$WIND$WESTERLY[lon,lat,pos_start:pos_stop]*spring_mask), na.rm = T)
      }
      
      if(year<1150){DATA_past1000$SIM_seasonal$SPRING$SLPR[lon,lat,year] = mean(DATA_past1000_SIM_RAW$SLPR[lon,lat, pos_start:pos_stop]*spring_mask, na.rm = T)}
      else{DATA_past1000$SIM_seasonal$SPRING$SLPR[lon,lat,year] = DATA_past1000$SIM_seasonal$SPRING$SLPR[lon,lat,year-1]}
    }
  }
}

#################################################
## SIMULATION MEAN ##############################
#################################################

DATA_past1000$SIM_mean <- list(
  TEMP = array(dim = c(96,73)),
  PREC = array(dim = c(96,73)),
  ISOT = array(dim = c(96,73)),
  ITPC = array(dim = c(96,73)),
  SLPR = array(dim = c(96,73)),
  WIND = list(
    STRENGTH = array(dim = c(96,73)),
    ANGLE = array(dim = c(96,73))
  )
)

for(lon in 1:96){
  for(lat in 1:73){
    DATA_past1000$SIM_mean$TEMP[lon,lat] = mean(DATA_past1000$SIM_yearly$TEMP[lon,lat,], na.rm = T)
    DATA_past1000$SIM_mean$PREC[lon,lat] = mean(DATA_past1000$SIM_yearly$PREC[lon,lat,], na.rm = T)
    DATA_past1000$SIM_mean$ISOT[lon,lat] = mean(DATA_past1000$SIM_yearly$ISOT[lon,lat,], na.rm = T)
    DATA_past1000$SIM_mean$ITPC[lon,lat] = mean(DATA_past1000$SIM_yearly$ITPC[lon,lat,], na.rm = T)
    DATA_past1000$SIM_mean$SLPR[lon,lat] = mean(DATA_past1000$SIM_yearly$SLPR[lon,lat,], na.rm = T)
    #DATA_past1000$SIM_mean$EVAP[lon,lat] = mean(DATA_past1000$SIM_yearly$EVAP[lon,lat,], na.rm = T)
    
    if(lat < 73){
      DATA_past1000$SIM_mean$WIND$STRENGTH[lon,lat] = mean(DATA_past1000$SIM_yearly$WIND$STRENGTH[lon,lat,], na.rm = T)
      DATA_past1000$SIM_mean$WIND$ANGLE[lon,lat] = mean(DATA_past1000$SIM_yearly$WIND$ANGLE[lon,lat,], na.rm = T)
    }
  }
}

remove(lon,lat, year, pos_start, pos_stop)

remove(DATA_past1000_SIM_RAW, autumn_mask, spring_mask, summer_mask, winter_mask)


#################################################
## ENTITY in LAT band ###########################
#################################################

lat_band <- as.numeric(c(90,80,70,60,50,40,30,20,10,0,-10,-20,-30,-40,-50,-60,-70,-80))
entity_lat <- array(dim = c(length(DATA_past1000$CAVES$entity_info$entity_id),2))
entity_lat[,1] <- DATA_past1000$CAVES$entity_info$entity_id

for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
  entity = entity_lat[ii,1]
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  lat = DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site]
  for(jj in length(lat_band):1){
    if(lat < lat_band[jj]){
      entity_lat[ii,2] = lat_band[jj]
      break
    }
  }
}

colnames(entity_lat) <- c("entity_id", "lat_band")

DATA_past1000$CAVES$entity_lat <- entity_lat

#################################################
## ARAGONITE, CALCITE, SMOW #####################
#################################################

## Aragonite and Calcite d18O have to be converted to drip water equivalents to be comparable
## Also in ORder to be comparable to SMOW which is the standard of the simulation we need to convert the dripwater d18O from VPDB to SMOW

for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
  entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
  site = DATA_past1000$CAVES$entity_info$site_id[ii]
  
  print(entity)
  
  data_rec <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
  data_sim <- DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]
  if(DATA_past1000$CAVES$entity_info$mineralogy[ii] == "calcite"){
    dw_eq <- 1.03092 * (data_rec$d18O_measurement - ((16.1*1000)/(data_sim$Temp+273.15)-24.6)) + 30.92
  }else if(DATA_past1000$CAVES$entity_info$mineralogy[ii] == "aragonite"){
    dw_eq <- 1.03092 * (data_rec$d18O_measurement - ((18.34*1000)/(data_sim$Temp+273.15)-31.954)) + 30.92
  }else{
    dw_eq <- numeric(length(data))+NA
  }
  
  DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]] <- data.frame(
    site_id = data_rec$site_id,
    entity_id = data_rec$entity_id,
    interp_age = data_rec$interp_age,
    d18O_measurement = data_rec$d18O_measurement,
    d18O_dw_eq = dw_eq
  )
  
}

remove(entity, dw_eq, ii, lats, site, data_rec, data_sim)

#################################################
## UMWANDLUNGEN #################################
#################################################


DATA_past1000$CAVES$site_info <- DATA_past1000$CAVES$site_info %>% mutate(elevation = as.numeric(as.character(elevation)))

remove(DATA_past1000_SIM_RAW)
