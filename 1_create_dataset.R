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

library(plyr)
library(dplyr)
library(stacy.hadcm.tools)
library(PaleoSpec)
library(nest)

#################################################
##0.1) Set Data-Structure?#######################
#################################################

DATA_past1000 <- list()
for(run in c("a","b","c")){
  DATA_past1000[[paste0("SIM_yearly_",run)]] <- list(
    TEMP = list(),
    PREC = list(), 
    ISOT = list(),
    SLPR = list()
  )
}

DATA_past1000$CAVES <- list()
DATA_past1000$CAVES$site_info <- read.csv("/stacywork/ginnyweasley/02_SISAL/SISAL_v2/site_countries.csv")
DATA_past1000$CAVES$entity_info <- list()
DATA_past1000$CAVES$record_data <- list()
DATA_past1000$CAVES$sim_data_yearly <- list()
DATA_past1000$CAVES$sim_data_downsampled <- list()


DATA_past1000_SIM_RAW <- list()
for(run in c("a","b","c")){
  DATA_past1000_SIM_RAW[[paste0(run)]] <- list(
    TEMP = list(),
    PREC = list(), 
    ISOT = list(),
    SLPR = list()
  )
}

#################################################
## 1) Read in DATA ##############################
#################################################

source("Functions/clear_data_matrix.R")

for(run in c("a", "b", "c")){
  print(run)
  print("TEMP")
  ncf <- (ncdf4::nc_open(paste0("/stacywork/hadcm3/surface_temperature/monthly_fixed/xnap", run, ".nc")))
  DATA_past1000_SIM_RAW[[run]]$TEMP <- clear_data_matrix(ncdf4::ncvar_get(ncf),1)
  DATA_past1000_SIM_RAW[[run]]$TEMP_t <- ncf$dim$t$vals
  ncdf4::nc_close(ncf)

  ncf<-ncdf4::nc_open(paste0("/stacywork/hadcm3/precipitation/monthly_fixed/xnap",run,".nc"))
  DATA_past1000_SIM_RAW[[run]]$PREC <- clear_data_matrix(ncdf4::ncvar_get(ncf),2)
  DATA_past1000_SIM_RAW[[run]]$PREC_t <- ncf$dim$t$vals
  DATA_past1000_SIM_RAW$lon <- ncf$dim$longitude$vals
  DATA_past1000_SIM_RAW$lat <- ncf$dim$latitude$vals
  ncdf4::nc_close(ncf)
  print("ISOT")
  ncf<-ncdf4::nc_open(paste0("/stacywork/hadcm3/isotopes/monthly_fixed/xnap",run,".nc"))
  DATA_past1000_SIM_RAW[[run]]$ISOT <- clear_data_matrix(ncdf4::ncvar_get(ncf, 'dO18'),3)
  DATA_past1000_SIM_RAW[[run]]$ISOT_t <- ncf$dim$t$vals
  ncdf4::nc_close(ncf)
  
  ncf<-ncdf4::nc_open(paste0("/stacywork/hadcm3/sea_level_pressure/monthly_fixed/xnap",run,".nc"))
  DATA_past1000_SIM_RAW[[run]]$SLPR <- ncdf4::ncvar_get(ncf)
  DATA_past1000_SIM_RAW[[run]]$SLPR_t <- ncf$dim$t$vals
  ncdf4::nc_close(ncf)
  
}

DATA_past1000$SIM_mean <- list()

ncf<-ncdf4::nc_open("/stacywork/hadcm3/hadcm_oro_001cyr.nc")
DATA_past1000$SIM_mean$OROG <- ncdf4::ncvar_get(ncf)
ncdf4::nc_close(ncf)

remove(ncf)

remove(run)
remove(clear_data_matrix)

#################################################
## Align timeseries #############################
#################################################

year_start = 810
year_stop = 1920

for(run in c("a","b","c")){
  for(var in c("TEMP", "PREC", "ISOT", "SLPR")){
    #we want to start in December, which is why we start one position earlier!
    pos_start = length(DATA_past1000_SIM_RAW[[run]][[paste0(var, "_t")]][DaysSinceToAD(DATA_past1000_SIM_RAW[[run]][[paste0(var, "_t")]])<year_start]) -1
    pos_stop = length(DATA_past1000_SIM_RAW[[run]][[paste0(var, "_t")]]) - length(DATA_past1000_SIM_RAW[[run]][[paste0(var, "_t")]][DaysSinceToAD(DATA_past1000_SIM_RAW[[run]][[paste0(var, "_t")]])>year_stop]) -1
    
    DATA_past1000_SIM_RAW[[run]][[var]] <- DATA_past1000_SIM_RAW[[run]][[var]][,,pos_start:pos_stop]
    DATA_past1000_SIM_RAW[[run]][[paste0(var,"_t")]] <- DATA_past1000_SIM_RAW[[run]][[paste0(var,"_t")]][pos_start:pos_stop]
  }
  
}

rm(run,var, pos_start, pos_stop)

DATA_past1000$time <- c(year_start, year_stop)
remove(year_start, year_stop)

#################################################
##2) SISAL TIME SERIES ##########################
#################################################

# needs to be imported first, as then only the relevant cave sites will be extracted and calculated further

source("Functions/fun_with_SISALv2_Janica.R")

data <- load_sisal_data_janica(year_start = DATA_past1000$time[1], year_stop = DATA_past1000$time[2])
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

remove(data, site, ii, name, load_sisal_data_janica)


#################################################
## 3) Extract data from Caves such that they are in a grid box that is the average of all surrounding
#################################################


source("Functions/extract_gridboxes.R")

for(run in c("a", "b", "c")){
  print(run)
  for (ii in 1:(dim(DATA_past1000$CAVES$site_info)[1])){
    lon_cave = DATA_past1000$CAVES$site_info$longitude[ii]
    
    if(lon_cave<0){lon_cave = 360+lon_cave}
    
    lat_cave = DATA_past1000$CAVES$site_info$latitude[ii]
    site_id = DATA_past1000$CAVES$site_info$site_id[ii]
    
    ratios <- extract_gridboxes(lon_cave, lat_cave)
    
    name <- paste0("CAVE",site_id)
    DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("TEMP_",run)]] <- ratios$E1*DATA_past1000_SIM_RAW[[run]]$TEMP[ratios$E1_lon_pos, ratios$E1_lat_pos,] +
      ratios$E2*DATA_past1000_SIM_RAW[[run]]$TEMP[ratios$E2_lon_pos, ratios$E2_lat_pos,] +
      ratios$E3*DATA_past1000_SIM_RAW[[run]]$TEMP[ratios$E3_lon_pos, ratios$E3_lat_pos,] +
      ratios$E4*DATA_past1000_SIM_RAW[[run]]$TEMP[ratios$E4_lon_pos, ratios$E4_lat_pos,]
    
    DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("PREC_",run)]] <- ratios$E1*DATA_past1000_SIM_RAW$PREC[[run]][ratios$E1_lon_pos, ratios$E1_lat_pos,] +
      ratios$E2*DATA_past1000_SIM_RAW[[run]]$PREC[ratios$E2_lon_pos, ratios$E2_lat_pos,] +
      ratios$E3*DATA_past1000_SIM_RAW[[run]]$PREC[ratios$E3_lon_pos, ratios$E3_lat_pos,] +
      ratios$E4*DATA_past1000_SIM_RAW[[run]]$PREC[ratios$E4_lon_pos, ratios$E4_lat_pos,]
    
    DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("ISOT_",run)]] <- ratios$E1*DATA_past1000_SIM_RAW[[run]]$ISOT[ratios$E1_lon_pos, ratios$E1_lat_pos,] +
      ratios$E2*DATA_past1000_SIM_RAW[[run]]$ISOT[ratios$E2_lon_pos, ratios$E2_lat_pos,] +
      ratios$E3*DATA_past1000_SIM_RAW[[run]]$ISOT[ratios$E3_lon_pos, ratios$E3_lat_pos,] +
      ratios$E4*DATA_past1000_SIM_RAW[[run]]$ISOT[ratios$E4_lon_pos, ratios$E4_lat_pos,]
    
    DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("SLPR_",run)]] <- ratios$E1*DATA_past1000_SIM_RAW[[run]]$SLPR[ratios$E1_lon_pos, ratios$E1_lat_pos,] +
      ratios$E2*DATA_past1000_SIM_RAW[[run]]$SLPR[ratios$E2_lon_pos, ratios$E2_lat_pos,] +
      ratios$E3*DATA_past1000_SIM_RAW[[run]]$SLPR[ratios$E3_lon_pos, ratios$E3_lat_pos,] +
      ratios$E4*DATA_past1000_SIM_RAW[[run]]$SLPR[ratios$E4_lon_pos, ratios$E4_lat_pos,]
    
    
  }
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

remove(ratios, ii, lat_cave, lon_cave, name, site_id, elevation_cave_sim, run, extract_gridboxes)

#################################################
## 4.0) SIMULATION YEARLY #######################
#################################################
for(run in c("a","b","c")){
  DATA_past1000[[paste0("SIM_yearly_",run)]]$TEMP <- array(dim = c(96,73,diff(DATA_past1000$time)))
  DATA_past1000[[paste0("SIM_yearly_",run)]]$PREC <- array(dim = c(96,73,diff(DATA_past1000$time)))
  DATA_past1000[[paste0("SIM_yearly_",run)]]$ISOT <- array(dim = c(96,73,diff(DATA_past1000$time)))
  DATA_past1000[[paste0("SIM_yearly_",run)]]$ITPC <- array(dim = c(96,73,diff(DATA_past1000$time))) # <- prec weighted mean
  DATA_past1000[[paste0("SIM_yearly_",run)]]$SLPR <- array(dim = c(96,73,diff(DATA_past1000$time)))
  
  for (lon in (1:96)){
    for (lat in 1:73){
      print(paste(lon,lat))
      for(year in 1:diff(DATA_past1000$time)){
        pos_start = 12*(year-1)+1
        pos_stop  = 12*(year-1)+12
        DATA_past1000[[paste0("SIM_yearly_",run)]]$TEMP[lon,lat,year] = mean(DATA_past1000_SIM_RAW[[run]]$TEMP[lon,lat, pos_start:pos_stop], na.rm = T)
        DATA_past1000[[paste0("SIM_yearly_",run)]]$PREC[lon,lat,year] = mean(DATA_past1000_SIM_RAW[[run]]$PREC[lon,lat, pos_start:pos_stop], na.rm = T)
        DATA_past1000[[paste0("SIM_yearly_",run)]]$ISOT[lon,lat,year] = mean(DATA_past1000_SIM_RAW[[run]]$ISOT[lon,lat, pos_start:pos_stop], na.rm = T)
        if(year<length(DATA_past1000_SIM_RAW[[run]]$SLPR_t)/12){DATA_past1000[[paste0("SIM_yearly_",run)]]$SLPR[lon,lat,year] = mean(DATA_past1000_SIM_RAW[[run]]$SLPR[lon,lat, pos_start:pos_stop], na.rm = T)}
        else{DATA_past1000[[paste0("SIM_yearly_",run)]]$SLPR[lon,lat,year] = DATA_past1000[[paste0("SIM_yearly_",run)]]$SLPR[lon,lat,year-1]}
        
        DATA_past1000[[paste0("SIM_yearly_",run)]]$ITPC[lon,lat,year] = sum(DATA_past1000_SIM_RAW[[run]]$PREC[lon,lat, pos_start:pos_stop]*DATA_past1000_SIM_RAW[[run]]$ISOT[lon,lat, pos_start:pos_stop],
                                                          na.rm = T)/sum(DATA_past1000_SIM_RAW[[run]]$PREC[lon,lat, pos_start:pos_stop], na.rm = T)
      }
    }
  }
}


remove(lon,lat,year, pos_start, pos_stop, run)

for(run in c("a", "b", "c")){
  print(run)
  for (ii in 1:(dim(DATA_past1000$CAVES$site_info)[1])){
    print(ii)
    site_id = DATA_past1000$CAVES$site_info$site_id[ii]
    name <- paste0("CAVE", site_id)
    DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("TEMP_", run)]] <- numeric(diff(DATA_past1000$time))
    DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("PREC_", run)]] <- numeric(diff(DATA_past1000$time))
    DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("ISOT_", run)]] <- numeric(diff(DATA_past1000$time))
    DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("ITPC_", run)]] <- numeric(diff(DATA_past1000$time))
    DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("SLPR_", run)]] <- numeric(diff(DATA_past1000$time))
    for(year in 1:diff(DATA_past1000$time)){
      pos_start = 12*(year-1)+1
      pos_stop  = 12*(year-1)+12
      DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("TEMP_", run)]][year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("TEMP_", run)]][pos_start:pos_stop], na.rm = T)
      DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("PREC_", run)]][year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("PREC_", run)]][pos_start:pos_stop], na.rm = T)
      DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("ISOT_", run)]][year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("ISOT_", run)]][pos_start:pos_stop], na.rm = T)
      DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("SLPR_", run)]][year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("SLPR_", run)]][pos_start:pos_stop], na.rm = T)
      DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("ITPC_", run)]][year] <- sum(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("ISOT_", run)]][pos_start:pos_stop]*
                                                                      DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("PREC_", run)]][pos_start:pos_stop], na.rm = T)/
        sum(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("PREC_", run)]][pos_start:pos_stop], na.rm = T)
    }
  }
}


remove(site_id, year, pos_start, pos_stop, name, ii, run)
#################################################
## 5.0) Seasonal Data ###########################
#################################################

winter_mask = c( 1, 1, 1,NA,NA,NA,NA,NA,NA,NA,NA,NA)
spring_mask = c(NA,NA,NA, 1, 1, 1,NA,NA,NA,NA,NA,NA)
summer_mask = c(NA,NA,NA,NA,NA,NA, 1, 1, 1,NA,NA,NA)
autumn_mask = c(NA,NA,NA,NA,NA,NA,NA,NA,NA, 1, 1, 1)

#DATA_past1000$CAVES$sim_data_seasonal <- vector(mode = "list")

for(run in c("a", "b", "c")){
  for (ii in 1:(dim(DATA_past1000$CAVES$site_info)[1])){
    print(ii)
    site_id = DATA_past1000$CAVES$site_info$site_id[ii]
    name = paste0("CAVE", site_id)
    
    DATA_past1000$CAVES$sim_data_seasonal[[name]] = list(
      SUMMER = list(temp_a_mean = numeric(diff(DATA_past1000$time)), prec_a_mean = numeric(diff(DATA_past1000$time)), isot_a_mean = numeric(diff(DATA_past1000$time)),
                    temp_b_mean = numeric(diff(DATA_past1000$time)), prec_b_mean = numeric(diff(DATA_past1000$time)), isot_b_mean = numeric(diff(DATA_past1000$time)),
                    temp_c_mean = numeric(diff(DATA_past1000$time)), prec_c_mean = numeric(diff(DATA_past1000$time)), isot_c_mean = numeric(diff(DATA_past1000$time))),
      AUTUMN = list(temp_a_mean = numeric(diff(DATA_past1000$time)), prec_a_mean = numeric(diff(DATA_past1000$time)), isot_a_mean = numeric(diff(DATA_past1000$time)),
                    temp_b_mean = numeric(diff(DATA_past1000$time)), prec_b_mean = numeric(diff(DATA_past1000$time)), isot_b_mean = numeric(diff(DATA_past1000$time)),
                    temp_c_mean = numeric(diff(DATA_past1000$time)), prec_c_mean = numeric(diff(DATA_past1000$time)), isot_c_mean = numeric(diff(DATA_past1000$time))),
      WINTER = list(temp_a_mean = numeric(diff(DATA_past1000$time)), prec_a_mean = numeric(diff(DATA_past1000$time)), isot_a_mean = numeric(diff(DATA_past1000$time)),
                    temp_b_mean = numeric(diff(DATA_past1000$time)), prec_b_mean = numeric(diff(DATA_past1000$time)), isot_b_mean = numeric(diff(DATA_past1000$time)),
                    temp_c_mean = numeric(diff(DATA_past1000$time)), prec_c_mean = numeric(diff(DATA_past1000$time)), isot_c_mean = numeric(diff(DATA_past1000$time))),
      SPRING = list(temp_a_mean = numeric(diff(DATA_past1000$time)), prec_a_mean = numeric(diff(DATA_past1000$time)), isot_a_mean = numeric(diff(DATA_past1000$time)),
                    temp_b_mean = numeric(diff(DATA_past1000$time)), prec_b_mean = numeric(diff(DATA_past1000$time)), isot_b_mean = numeric(diff(DATA_past1000$time)),
                    temp_c_mean = numeric(diff(DATA_past1000$time)), prec_c_mean = numeric(diff(DATA_past1000$time)), isot_c_mean = numeric(diff(DATA_past1000$time)))
    )
    
    
    for(year in 1:diff(DATA_past1000$time)){
      pos_start = 12*(year-1)+1
      pos_stop  = 12*(year-1)+12
      DATA_past1000$CAVES$sim_data_seasonal[[name]]$SUMMER[[paste0("temp_", run, "_mean")]][year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("TEMP_", run)]][pos_start:pos_stop]*summer_mask, na.rm = T)
      DATA_past1000$CAVES$sim_data_seasonal[[name]]$SUMMER[[paste0("prec_", run, "_mean")]][year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("PREC_", run)]][pos_start:pos_stop]*summer_mask, na.rm = T)
      DATA_past1000$CAVES$sim_data_seasonal[[name]]$SUMMER[[paste0("isot_", run, "_mean")]][year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("ISOT_", run)]][pos_start:pos_stop]*summer_mask, na.rm = T)
      
      DATA_past1000$CAVES$sim_data_seasonal[[name]]$AUTUMN[[paste0("temp_", run, "_mean")]][year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("TEMP_", run)]][pos_start:pos_stop]*autumn_mask, na.rm = T)
      DATA_past1000$CAVES$sim_data_seasonal[[name]]$AUTUMN[[paste0("prec_", run, "_mean")]][year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("PREC_", run)]][pos_start:pos_stop]*autumn_mask, na.rm = T)
      DATA_past1000$CAVES$sim_data_seasonal[[name]]$AUTUMN[[paste0("isot_", run, "_mean")]][year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("ISOT_", run)]][pos_start:pos_stop]*autumn_mask, na.rm = T)
      
      DATA_past1000$CAVES$sim_data_seasonal[[name]]$WINTER[[paste0("temp_", run, "_mean")]][year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("TEMP_", run)]][pos_start:pos_stop]*winter_mask, na.rm = T)
      DATA_past1000$CAVES$sim_data_seasonal[[name]]$WINTER[[paste0("prec_", run, "_mean")]][year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("PREC_", run)]][pos_start:pos_stop]*winter_mask, na.rm = T)
      DATA_past1000$CAVES$sim_data_seasonal[[name]]$WINTER[[paste0("isot_", run, "_mean")]][year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("ISOT_", run)]][pos_start:pos_stop]*winter_mask, na.rm = T)
      
      DATA_past1000$CAVES$sim_data_seasonal[[name]]$SPRING[[paste0("temp_", run, "_mean")]][year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("TEMP_", run)]][pos_start:pos_stop]*spring_mask, na.rm = T)
      DATA_past1000$CAVES$sim_data_seasonal[[name]]$SPRING[[paste0("prec_", run, "_mean")]][year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("PREC_", run)]][pos_start:pos_stop]*spring_mask, na.rm = T)
      DATA_past1000$CAVES$sim_data_seasonal[[name]]$SPRING[[paste0("isot_", run, "_mean")]][year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("ISOT_", run)]][pos_start:pos_stop]*spring_mask, na.rm = T)
      
    }
  }
}



remove(ii, name, pos_start, pos_stop, year, site_id, run)

# #################################################
# ## 6.0) CLIMATE INFO ############################
# #################################################
# 
# source("Functions/aw_mean.R")
# 
# DATA_past1000$CAVES$climate_info <- list(
#   mean_temp = numeric((dim(DATA_past1000$CAVES$site_info)[1])), max_temp = numeric((dim(DATA_past1000$CAVES$site_info)[1])), min_temp = numeric((dim(DATA_past1000$CAVES$site_info)[1])),
#   mean_prec = numeric((dim(DATA_past1000$CAVES$site_info)[1])), max_prec = numeric((dim(DATA_past1000$CAVES$site_info)[1])), min_prec = numeric((dim(DATA_past1000$CAVES$site_info)[1])),
#   mean_isot = numeric((dim(DATA_past1000$CAVES$site_info)[1])), max_isot = numeric((dim(DATA_past1000$CAVES$site_info)[1])), min_isot = numeric((dim(DATA_past1000$CAVES$site_info)[1])),
#   mean_slpr = numeric((dim(DATA_past1000$CAVES$site_info)[1])), max_slpr = numeric((dim(DATA_past1000$CAVES$site_info)[1])), min_slpr = numeric((dim(DATA_past1000$CAVES$site_info)[1]))
# )
# 
# for (ii in 1:(dim(DATA_past1000$CAVES$site_info)[1])){
#   site_id = DATA_past1000$CAVES$site_info$site_id[ii]
#   name = paste0("CAVE", site_id)
#   DATA_past1000$CAVES$climate_info$mean_temp[ii] = mean(DATA_past1000$CAVES$sim_data_raw[[name]]$TEMP, na.rm = T)
#   DATA_past1000$CAVES$climate_info$max_temp[ii] = max(DATA_past1000$CAVES$sim_data_raw[[name]]$TEMP, na.rm = T)
#   DATA_past1000$CAVES$climate_info$min_temp[ii] = min(DATA_past1000$CAVES$sim_data_raw[[name]]$TEMP, na.rm = T)
# 
#   DATA_past1000$CAVES$climate_info$mean_prec[ii] = mean(DATA_past1000$CAVES$sim_data_raw[[name]]$PREC, na.rm = T)
#   DATA_past1000$CAVES$climate_info$max_prec[ii] = max(DATA_past1000$CAVES$sim_data_raw[[name]]$PREC, na.rm = T)
#   DATA_past1000$CAVES$climate_info$min_prec[ii] = min(DATA_past1000$CAVES$sim_data_raw[[name]]$PREC, na.rm = T)
#   
#   DATA_past1000$CAVES$climate_info$mean_isot[ii] = mean(DATA_past1000$CAVES$sim_data_raw[[name]]$ISOT, na.rm = T)
#   DATA_past1000$CAVES$climate_info$max_isot[ii] = max(DATA_past1000$CAVES$sim_data_raw[[name]]$ISOT, na.rm = T)
#   DATA_past1000$CAVES$climate_info$min_isot[ii] = min(DATA_past1000$CAVES$sim_data_raw[[name]]$ISOT, na.rm = T)
#   
#   DATA_past1000$CAVES$climate_info$mean_slpr[ii] = mean(DATA_past1000$CAVES$sim_data_raw[[name]]$SLPR, na.rm = T)
#   DATA_past1000$CAVES$climate_info$max_slpr[ii] = max(DATA_past1000$CAVES$sim_data_raw[[name]]$SLPR, na.rm = T)
#   DATA_past1000$CAVES$climate_info$min_slpr[ii] = min(DATA_past1000$CAVES$sim_data_raw[[name]]$SLPR, na.rm = T)
# }
# 
#   remove(ii, name, site_id)
#   

#################################################
## 7.0) DOWNSAMPELING ###########################
#################################################

source("Functions/SubsampleTimeseriesBlock_highresNA.R")

for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
  print(ii)
  nameE = paste0("ENTITY", DATA_past1000$CAVES$entity_info$entity_id[ii])
  nameC = paste0("CAVE", DATA_past1000$CAVES$entity_info$site_id[ii])
  for(run in c("a", "b", "c")){
    assign(paste0("data_temp_", run), SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_yearly[[nameC]][[paste0("TEMP_",run)]]),
                                                                            start = 1950-DATA_past1000$time[2],
                                                                            end = 1950-DATA_past1000$time[1]),
                                                                         DATA_past1000$CAVES$record_data[[nameE]]$interp_age))
    
    assign(paste0("data_prec_", run), SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_yearly[[nameC]][[paste0("PREC_",run)]]),
                                                       start = 1950-DATA_past1000$time[2],
                                                       end = 1950-DATA_past1000$time[1]),
                                                    DATA_past1000$CAVES$record_data[[nameE]]$interp_age))

    assign(paste0("data_isot_", run), SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_yearly[[nameC]][[paste0("ISOT_",run)]]),
                                                                            start = 1950-DATA_past1000$time[2],
                                                                            end = 1950-DATA_past1000$time[1]),
                                                                         DATA_past1000$CAVES$record_data[[nameE]]$interp_age))
    
    assign(paste0("data_itpc_", run), SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_yearly[[nameC]][[paste0("ITPC_",run)]]),
                                                                            start = 1950-DATA_past1000$time[2],
                                                                            end = 1950-DATA_past1000$time[1]),
                                                                         DATA_past1000$CAVES$record_data[[nameE]]$interp_age))
  }
  data <- matrix(c(DATA_past1000$CAVES$record_data[[nameE]]$interp_age, 
                   data_temp_a, data_prec_a, data_isot_a, data_itpc_a,
                   data_temp_b, data_prec_b, data_isot_b, data_itpc_b,
                   data_temp_c, data_prec_c, data_isot_c, data_itpc_c), ncol= 13)
  colnames(data) = c("interp_age", 
                     "TEMP_a", "PREC_a","ISOT_a", "ITPC_a", 
                     "TEMP_b", "PREC_b","ISOT_b", "ITPC_b", 
                     "TEMP_c", "PREC_c","ISOT_c", "ITPC_c")
  
  DATA_past1000$CAVES$sim_data_downsampled[[nameE]] <- as.tibble(data)
}

remove(nameE, nameC, ii, data,
       data_isot_a, data_itpc_a, data_prec_a, data_temp_a,
       data_isot_b, data_itpc_b, data_prec_b, data_temp_b,
       data_isot_c, data_itpc_c, data_prec_c, data_temp_c, run)
remove(SubsampleTimeseriesBlock_highresNA)

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
    dw_eq_a <- 1.03092 * (data_rec$d18O_measurement - ((16.1*1000)/(data_sim$TEMP_a+273.15)-24.6)) + 30.92
    dw_eq_b <- 1.03092 * (data_rec$d18O_measurement - ((16.1*1000)/(data_sim$TEMP_b+273.15)-24.6)) + 30.92
    dw_eq_c <- 1.03092 * (data_rec$d18O_measurement - ((16.1*1000)/(data_sim$TEMP_c+273.15)-24.6)) + 30.92
  }else if(DATA_past1000$CAVES$entity_info$mineralogy[ii] == "aragonite"){
    dw_eq_a <- 1.03092 * (data_rec$d18O_measurement - ((18.34*1000)/(data_sim$TEMP_a+273.15)-31.954)) + 30.92
    dw_eq_b <- 1.03092 * (data_rec$d18O_measurement - ((18.34*1000)/(data_sim$TEMP_b+273.15)-31.954)) + 30.92
    dw_eq_c <- 1.03092 * (data_rec$d18O_measurement - ((18.34*1000)/(data_sim$TEMP_c+273.15)-31.954)) + 30.92
  }else{
    dw_eq_a <- numeric(length(data))+NA
    dw_eq_b <- numeric(length(data))+NA
    dw_eq_c <- numeric(length(data))+NA
  }

  DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]] <- data.frame(
    site_id = data_rec$site_id,
    entity_id = data_rec$entity_id,
    interp_age = data_rec$interp_age,
    d18O_measurement = data_rec$d18O_measurement,
    d18O_dw_eq_a = dw_eq_a,
    d18O_dw_eq_b = dw_eq_b,
    d18O_dw_eq_c = dw_eq_c
  )
}

remove(entity, dw_eq_a, dw_eq_b, dw_eq_c, ii, site, data_rec, data_sim)

#################################################
## UMWANDLUNGEN #################################
#################################################

DATA_past1000$CAVES$site_info <- DATA_past1000$CAVES$site_info %>% mutate(elevation = as.numeric(as.character(elevation)))
remove(DATA_past1000_SIM_RAW)

#################################################
## MASKS ########################################
#################################################

mask_mean = logical(length = length(DATA_past1000$CAVES$entity_info$entity_id))
mask_var  = logical(length = length(DATA_past1000$CAVES$entity_info$entity_id))
mask_spec = logical(length = length(DATA_past1000$CAVES$entity_info$entity_id))

for(entity in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
  if(DATA_past1000$CAVES$entity_info$n[entity] > 10 & DATA_past1000$CAVES$entity_info$period[entity] > 600){mask_mean[entity] = T}
  if(DATA_past1000$CAVES$entity_info$n[entity] > 20 & DATA_past1000$CAVES$entity_info$period[entity] > 600){mask_var[entity] = T}
  if(DATA_past1000$CAVES$entity_info$n[entity] > 30 & DATA_past1000$CAVES$entity_info$period[entity] > 600){mask_spec[entity] = T}
}

#   #################################################
#   ## 8.0) SIMULATION SEASONS ######################
#   #################################################
#   
#   DATA_past1000$SIM_seasonal <- list(
#     SUMMER = list(
#       TEMP = array(dim = c(96,73,diff(DATA_past1000$time))),
#       PREC = array(dim = c(96,73,diff(DATA_past1000$time))),
#       ISOT = array(dim = c(96,73,diff(DATA_past1000$time))),
#       ITPC = array(dim = c(96,73,diff(DATA_past1000$time))),
#       SLPR = array(dim = c(96,73,diff(DATA_past1000$time)))
#     ),
#     AUTUMN = list(
#       TEMP = array(dim = c(96,73,diff(DATA_past1000$time))),
#       PREC = array(dim = c(96,73,diff(DATA_past1000$time))),
#       ISOT = array(dim = c(96,73,diff(DATA_past1000$time))),
#       ITPC = array(dim = c(96,73,diff(DATA_past1000$time))),
#       SLPR = array(dim = c(96,73,diff(DATA_past1000$time)))
#     ),
#     WINTER = list(
#       TEMP = array(dim = c(96,73,diff(DATA_past1000$time))),
#       PREC = array(dim = c(96,73,diff(DATA_past1000$time))),
#       ISOT = array(dim = c(96,73,diff(DATA_past1000$time))),
#       ITPC = array(dim = c(96,73,diff(DATA_past1000$time))),
#       SLPR = array(dim = c(96,73,diff(DATA_past1000$time)))
#     ),
#     SPRING = list(
#       TEMP = array(dim = c(96,73,diff(DATA_past1000$time))),
#       PREC = array(dim = c(96,73,diff(DATA_past1000$time))),
#       ISOT = array(dim = c(96,73,diff(DATA_past1000$time))),
#       ITPC = array(dim = c(96,73,diff(DATA_past1000$time))),
#       SLPR = array(dim = c(96,73,diff(DATA_past1000$time)))
#     )
#   )
#   
#   #################################################
#   ##Calc
#   
#   for(lon in 1:96){
#     for(lat in 1:73){
#       print(paste(lon,lat,sep = " "))
#         for(year in 1:diff(DATA_past1000$time)){
#         #SUMMER
#         pos_start = 12*(year-1)+1
#         pos_stop  = 12*(year-1)+12
#         
#         DATA_past1000$SIM_seasonal$SUMMER$TEMP[lon,lat,year] = mean(DATA_past1000_SIM_RAW$TEMP[lon,lat, pos_start:pos_stop]*summer_mask, na.rm = T)
#         DATA_past1000$SIM_seasonal$SUMMER$PREC[lon,lat,year] = mean(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*summer_mask, na.rm = T)
#         DATA_past1000$SIM_seasonal$SUMMER$ISOT[lon,lat,year] = mean(DATA_past1000_SIM_RAW$ISOT[lon,lat, pos_start:pos_stop]*summer_mask, na.rm = T)
#         DATA_past1000$SIM_seasonal$SUMMER$ITPC[lon,lat,year] = sum(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*DATA_past1000_SIM_RAW$ISOT[lon,lat, pos_start:pos_stop]*summer_mask,
#                                                           na.rm = T)/sum(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*summer_mask, na.rm = T)
#         if(year<diff(DATA_past1000$time)){DATA_past1000$SIM_seasonal$SUMMER$SLPR[lon,lat,year] = mean(DATA_past1000_SIM_RAW$SLPR[lon,lat, pos_start:pos_stop]*summer_mask, na.rm = T)}
#         else{DATA_past1000$SIM_seasonal$SUMMER$SLPR[lon,lat,year] = DATA_past1000$SIM_seasonal$SUMMER$SLPR[lon,lat,year-1]}
#         
#         #AUTUMN
#         DATA_past1000$SIM_seasonal$AUTUMN$TEMP[lon,lat,year] = mean(DATA_past1000_SIM_RAW$TEMP[lon,lat, pos_start:pos_stop]*autumn_mask, na.rm = T)
#         DATA_past1000$SIM_seasonal$AUTUMN$PREC[lon,lat,year] = mean(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*autumn_mask, na.rm = T)
#         DATA_past1000$SIM_seasonal$AUTUMN$ISOT[lon,lat,year] = mean(DATA_past1000_SIM_RAW$ISOT[lon,lat, pos_start:pos_stop]*autumn_mask, na.rm = T)
#         DATA_past1000$SIM_seasonal$AUTUMN$ITPC[lon,lat,year] = sum(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*DATA_past1000_SIM_RAW$ISOT[lon,lat, pos_start:pos_stop]*autumn_mask,
#                                                                   na.rm = T)/sum(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*autumn_mask, na.rm = T)
#         if(year<diff(DATA_past1000$time)){DATA_past1000$SIM_seasonal$AUTUMN$SLPR[lon,lat,year] = mean(DATA_past1000_SIM_RAW$SLPR[lon,lat, pos_start:pos_stop]*autumn_mask, na.rm = T)}
#         else{DATA_past1000$SIM_seasonal$AUTUMN$SLPR[lon,lat,year] = DATA_past1000$SIM_seasonal$AUTUMN$SLPR[lon,lat,year-1]}
#         
#         #WINTER
#         DATA_past1000$SIM_seasonal$WINTER$TEMP[lon,lat,year] = mean(DATA_past1000_SIM_RAW$TEMP[lon,lat, pos_start:pos_stop]*winter_mask, na.rm = T)
#         DATA_past1000$SIM_seasonal$WINTER$PREC[lon,lat,year] = mean(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*winter_mask, na.rm = T)
#         DATA_past1000$SIM_seasonal$WINTER$ISOT[lon,lat,year] = mean(DATA_past1000_SIM_RAW$ISOT[lon,lat, pos_start:pos_stop]*winter_mask, na.rm = T)
#         DATA_past1000$SIM_seasonal$WINTER$ITPC[lon,lat,year] = sum(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*DATA_past1000_SIM_RAW$ISOT[lon,lat, pos_start:pos_stop]*winter_mask,
#                                                                   na.rm = T)/sum(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*winter_mask, na.rm = T)
#         
#         if(year<diff(DATA_past1000$time)){DATA_past1000$SIM_seasonal$WINTER$SLPR[lon,lat,year] = mean(DATA_past1000_SIM_RAW$SLPR[lon,lat, pos_start:pos_stop]*winter_mask, na.rm = T)}
#         else{DATA_past1000$SIM_seasonal$WINTER$SLPR[lon,lat,year] = DATA_past1000$SIM_seasonal$WINTER$SLPR[lon,lat,year-1]}
#         
#         #SPRING
#         DATA_past1000$SIM_seasonal$SPRING$TEMP[lon,lat,year] = mean(DATA_past1000_SIM_RAW$TEMP[lon,lat, pos_start:pos_stop]*spring_mask, na.rm = T)
#         DATA_past1000$SIM_seasonal$SPRING$PREC[lon,lat,year] = mean(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*spring_mask, na.rm = T)
#         DATA_past1000$SIM_seasonal$SPRING$ISOT[lon,lat,year] = mean(DATA_past1000_SIM_RAW$ISOT[lon,lat, pos_start:pos_stop]*spring_mask, na.rm = T)
#         DATA_past1000$SIM_seasonal$SPRING$ITPC[lon,lat,year] = sum(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*DATA_past1000_SIM_RAW$ISOT[lon,lat, pos_start:pos_stop]*spring_mask,
#                                                                   na.rm = T)/sum(DATA_past1000_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*spring_mask, na.rm = T)
#         
#         if(year<diff(DATA_past1000$time)){DATA_past1000$SIM_seasonal$SPRING$SLPR[lon,lat,year] = mean(DATA_past1000_SIM_RAW$SLPR[lon,lat, pos_start:pos_stop]*spring_mask, na.rm = T)}
#         else{DATA_past1000$SIM_seasonal$SPRING$SLPR[lon,lat,year] = DATA_past1000$SIM_seasonal$SPRING$SLPR[lon,lat,year-1]}
#       }
#     }
#   }
#   
#   remove(lat,lon, pos_start, pos_stop, year)
#   
#   #################################################
#   ## SIMULATION MEAN ##############################
#   #################################################
#   
#   DATA_past1000$SIM_mean <- list(
#     TEMP = array(dim = c(96,73)),
#     PREC = array(dim = c(96,73)),
#     ISOT = array(dim = c(96,73)),
#     ITPC = array(dim = c(96,73)),
#     SLPR = array(dim = c(96,73))
#   )
#   
#   for(lon in 1:96){
#     for(lat in 1:73){
#       DATA_past1000$SIM_mean$TEMP[lon,lat] = mean(DATA_past1000$SIM_yearly$TEMP[lon,lat,], na.rm = T)
#       DATA_past1000$SIM_mean$PREC[lon,lat] = mean(DATA_past1000$SIM_yearly$PREC[lon,lat,], na.rm = T)
#       DATA_past1000$SIM_mean$ISOT[lon,lat] = mean(DATA_past1000$SIM_yearly$ISOT[lon,lat,], na.rm = T)
#       DATA_past1000$SIM_mean$ITPC[lon,lat] = mean(DATA_past1000$SIM_yearly$ITPC[lon,lat,], na.rm = T)
#       DATA_past1000$SIM_mean$SLPR[lon,lat] = mean(DATA_past1000$SIM_yearly$SLPR[lon,lat,], na.rm = T)
#       #DATA_past1000$SIM_mean$EVAP[lon,lat] = mean(DATA_past1000$SIM_yearly$EVAP[lon,lat,], na.rm = T)
#     }
#   }
#   
#   remove(lon,lat, year, pos_start, pos_stop)
#   
#   remove(DATA_past1000_SIM_RAW, autumn_mask, spring_mask, summer_mask, winter_mask)
#   
#   
#   #################################################
#   ## ENTITY in LAT band ###########################
#   #################################################
#   
# lat_band <- as.numeric(c(90,80,70,60,50,40,30,20,10,0,-10,-20,-30,-40,-50,-60,-70,-80))
#   entity_lat <- array(dim = c(length(DATA_past1000$CAVES$entity_info$entity_id),2))
#   entity_lat[,1] <- DATA_past1000$CAVES$entity_info$entity_id
#   
#   for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
#     entity = entity_lat[ii,1]
#     site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
#     lat = DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site]
#     for(jj in length(lat_band):1){
#       if(lat < lat_band[jj]){
#         entity_lat[ii,2] = lat_band[jj]
#         break
#       }
#     }
#   }
#   
#   colnames(entity_lat) <- c("entity_id", "lat_band")
#   
#   DATA_past1000$CAVES$entity_lat <- entity_lat
