########Main#####################################
##
## This Script creates the dataset for the 40 extracted
## caves from Bakers paper.
##
##

##
## Include that caves may not be in the center of a grid box but rather mix different grid boxes together depending on location ratio
##

##0.1) Set Data-Structure?

DATA_SISAL <- list()

DATA_SISAL$SIM_yearly <- list(
  TEMP = list(),
  PREC = list(), 
  ISOT = list()
)



DATA_SISAL$CAVES <- list()
DATA_SISAL$CAVES$site_info <- read.csv("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/02_SISAL/SISAL_v2_CARLA/site_countries.csv")
DATA_SISAL$CAVES$entity_info <- list()
DATA_SISAL$CAVES$record_data <- list()
DATA_SISAL$CAVES$sim_data_yearly <- list()


# this will be an extra list, so it can be removed easily later
DATA_SISAL_SIM_RAW <- list(
  TEMP = list(),
  PREC = list(), 
  ISOT = list()
)

#################################################
##1) read in simulation data from xnapa #########
#################################################

source("Functions/clear_data_matrix.R")

ncf<-ncdf4::nc_open("~/Dokumente/01_Promotion/06_Daten/05_HadCM3/xnapa/xnapa_surface_temperature.nc")
DATA_SISAL_SIM_RAW$TEMP <- clear_data_matrix(ncdf4::ncvar_get(ncf),1)
ncdf4::nc_close(ncf)

ncf<-ncdf4::nc_open("~/Dokumente/01_Promotion/06_Daten/05_HadCM3/xnapa/xnapa_precipitation.nc")
DATA_SISAL_SIM_RAW$PREC <- clear_data_matrix(ncdf4::ncvar_get(ncf),2)
DATA_SISAL_SIM_RAW$lon <- ncf$dim$longitude$vals
DATA_SISAL_SIM_RAW$lat <- ncf$dim$latitude$vals
ncdf4::nc_close(ncf)

ncf<-ncdf4::nc_open("~/Dokumente/01_Promotion/06_Daten/05_HadCM3/xnapa/xnapa_isotopes.nc")
DATA_SISAL_SIM_RAW$ISOT <- clear_data_matrix(ncdf4::ncvar_get(ncf, ncf$var[[3]]),3)
ncdf4::nc_close(ncf)

remove(ncf)

#################################################
##2) SISAL TIME SERIES ############################
#################################################

# needs to be imported first, as then only the relevant cave sites will be extracted and calculated further
library(plyr)
library(dplyr)
source("Functions/fun_with_SISALv2.R")

data <- load_sisal_data()
# afterwards select those that have at least one data point within the last millennium!!! otherwise too much useless data!

site_tb <- data[[1]]
dating_tb <- data[[2]]
sample_tb <- data[[3]]
run <- data[[4]]

entities_mineralogy <- sample_tb %>% filter(interp_age < 1100) %>% distinct(entity_id, mineralogy)
#entities_mineralogy <- sample_tb %>% filter(interp_age < 1100) %>% distinct(entity_id, mineralogy) %>% filter(entity_id != 107, entity_id != 188, 
#                                                                                                              entity_id != 190, entity_id != 197,
#                                                                                                              entity_id != 348, entity_id != 349, 
#                                                                                                              entity_id != 436)

# all entities that have at least one measurement within the past millennium
# get rid of all those entities with multiple mineralogy (id = 107, 188, 190, 197, 348, 349, 436)
entities_SISAL <- sample_tb %>% filter(interp_age < 1100) %>% group_by(entity_id) %>% count() %>% 
  right_join(., entities_mineralogy, by = "entity_id")


sample_min1 <- sample_tb %>% filter(entity_id %in% entities_SISAL$entity_id) %>% distinct(entity_id, interp_age, d18O_measurement) %>% filter(interp_age<1100)
site_min1 <- site_tb %>% filter(entity_id %in% entities_SISAL$entity_id) %>% distinct(site_id, entity_id) %>% right_join(., sample_min1, by = "entity_id")

site_period <- site_min1 %>% group_by(entity_id) %>% 
  summarise(min_corr_age = round(min(interp_age, na.rm = T), digits = 2),
            max_corr_age = round(max(interp_age, na.rm = T), digits = 2)) %>% 
  mutate(period = max_corr_age -min_corr_age)

site_to_entity <- site_tb %>% filter(entity_id %in% entities_SISAL$entity_id) %>%distinct(site_id, entity_id)

# General data for entities with more than 2 measurements
DATA_SISAL$CAVES$entity_info <- entities_SISAL %>% right_join(., site_period, by = "entity_id") %>% right_join(site_to_entity,., by = "entity_id")

remove(entities_mineralogy, entities_SISAL, sample_tb, site_min1, site_period, site_tb, site_to_entity, dating_tb, run)


#Schmeißt alle Höhlen raus, die nicht gebraucht werden
DATA_SISAL$CAVES$site_info <- DATA_SISAL$CAVES$site_info %>% filter(site_id %in% DATA_SISAL$CAVES$entity_info$site_id)

for (ii in DATA_SISAL$CAVES$entity_info$entity_id){
  name = paste0("ENTITY", ii)
  if(ii%%10 == 0){
    print(name)
  }
  
  site <- DATA_SISAL$CAVES$entity_info %>% filter(entity_id == ii) %>% distinct(site_id)
  DATA_SISAL$CAVES$record_data[[name]] <- sample_min1 %>% filter(entity_id == ii) %>% distinct(entity_id, mineralogy, arag_corr, interp_age, d18O_measurement) %>%
    mutate(site_id = (site$site_id))
}

remove(data, sample_min1, site, ii, name)

#################################################
## 3) Extract data from Caves such that they are in a grid box that is the average of all surrounding

source("Functions/extract_gridboxes_SISAL.R")

for (ii in 1:length(DATA_SISAL$CAVES$entity_info$entity_id)){
  entity = DATA_SISAL$CAVES$entity_info$entity_id[ii]
  if(ii%%10 == 0){
    print(entity)
  }
  site   = DATA_SISAL$CAVES$entity_info$site_id[ii]
  lon_cave = DATA_SISAL$CAVES$site_info$longitude[DATA_SISAL$CAVES$site_info$site_id == site]
  if(lon_cave<0){lon_cave = 360+lon_cave}
  lat_cave = DATA_SISAL$CAVES$site_info$latitude[DATA_SISAL$CAVES$site_info$site_id == site]
  
  ratios <- extract_gridboxes_sisal(lon_cave, lat_cave)
  
  name <- paste0("ENTITY", entity)
  DATA_SISAL$CAVES$sim_data_raw[[name]]$TEMP <- ratios$E1*DATA_SISAL_SIM_RAW$TEMP[ratios$E1_lon_pos, ratios$E1_lat_pos,] +
    ratios$E2*DATA_SISAL_SIM_RAW$TEMP[ratios$E2_lon_pos, ratios$E2_lat_pos,] +
    ratios$E3*DATA_SISAL_SIM_RAW$TEMP[ratios$E3_lon_pos, ratios$E3_lat_pos,] +
    ratios$E4*DATA_SISAL_SIM_RAW$TEMP[ratios$E4_lon_pos, ratios$E4_lat_pos,]
  
  DATA_SISAL$CAVES$sim_data_raw[[name]]$PREC <- ratios$E1*DATA_SISAL_SIM_RAW$PREC[ratios$E1_lon_pos, ratios$E1_lat_pos,] +
    ratios$E2*DATA_SISAL_SIM_RAW$PREC[ratios$E2_lon_pos, ratios$E2_lat_pos,] +
    ratios$E3*DATA_SISAL_SIM_RAW$PREC[ratios$E3_lon_pos, ratios$E3_lat_pos,] +
    ratios$E4*DATA_SISAL_SIM_RAW$PREC[ratios$E4_lon_pos, ratios$E4_lat_pos,]
  
  DATA_SISAL$CAVES$sim_data_raw[[name]]$ISOT <- ratios$E1*DATA_SISAL_SIM_RAW$ISOT[ratios$E1_lon_pos, ratios$E1_lat_pos,] +
    ratios$E2*DATA_SISAL_SIM_RAW$ISOT[ratios$E2_lon_pos, ratios$E2_lat_pos,] +
    ratios$E3*DATA_SISAL_SIM_RAW$ISOT[ratios$E3_lon_pos, ratios$E3_lat_pos,] +
    ratios$E4*DATA_SISAL_SIM_RAW$ISOT[ratios$E4_lon_pos, ratios$E4_lat_pos,]
}

remove(ratios, ii, lat_cave, lon_cave, name, site, entity)

#################################################
## SIMULATION YEARLY ############################
#################################################

DATA_SISAL$SIM_yearly$TEMP <- array(dim = c(96,73,1150))
DATA_SISAL$SIM_yearly$PREC <- array(dim = c(96,73,1150))
DATA_SISAL$SIM_yearly$ISOT <- array(dim = c(96,73,1150))
DATA_SISAL$SIM_yearly$ITPC <- array(dim = c(96,73,1150)) # <- prec weighted mean


for (lon in (1:96)){
  for (lat in 1:73){
    for(year in 1:1150){
      pos_start = 12*(year-1)+1
      pos_stop  = 12*(year-1)+12
      DATA_SISAL$SIM_yearly$TEMP[lon,lat,year] = mean(DATA_SISAL_SIM_RAW$TEMP[lon,lat, pos_start:pos_stop], na.rm = T)
      DATA_SISAL$SIM_yearly$PREC[lon,lat,year] = mean(DATA_SISAL_SIM_RAW$PREC[lon,lat, pos_start:pos_stop], na.rm = T)
      DATA_SISAL$SIM_yearly$ISOT[lon,lat,year] = mean(DATA_SISAL_SIM_RAW$ISOT[lon,lat, pos_start:pos_stop], na.rm = T)
      DATA_SISAL$SIM_yearly$ITPC[lon,lat,year] = sum(DATA_SISAL_SIM_RAW$PREC[lon,lat, pos_start:pos_stop]*DATA_SISAL_SIM_RAW$ISOT[lon,lat, pos_start:pos_stop],
                                                        na.rm = T)/sum(DATA_SISAL_SIM_RAW$PREC[lon,lat, pos_start:pos_stop], na.rm = T)
    }
  }
}

remove(lon,lat,year, pos_start, pos_stop, site_id)

for (ii in 1:length(DATA_SISAL$CAVES$entity_info$entity_id)){
  entity = DATA_SISAL$CAVES$entity_info$entity_id[ii]
  name <- paste0("ENTITY", entity)
  DATA_SISAL$CAVES$sim_data_yearly[[name]]$TEMP <- numeric(1150)
  DATA_SISAL$CAVES$sim_data_yearly[[name]]$PREC <- numeric(1150)
  DATA_SISAL$CAVES$sim_data_yearly[[name]]$ISOT <- numeric(1150)
  DATA_SISAL$CAVES$sim_data_yearly[[name]]$ITPC <- numeric(1150)
  for(year in 1:1150){
    pos_start = 12*(year-1)+1
    pos_stop  = 12*(year-1)+12
    DATA_SISAL$CAVES$sim_data_yearly[[name]]$TEMP[year] <- mean(DATA_SISAL$CAVES$sim_data_raw[[name]]$TEMP[pos_start:pos_stop], na.rm = T)
    DATA_SISAL$CAVES$sim_data_yearly[[name]]$PREC[year] <- mean(DATA_SISAL$CAVES$sim_data_raw[[name]]$PREC[pos_start:pos_stop], na.rm = T)
    DATA_SISAL$CAVES$sim_data_yearly[[name]]$ISOT[year] <- mean(DATA_SISAL$CAVES$sim_data_raw[[name]]$ISOT[pos_start:pos_stop], na.rm = T)
    DATA_SISAL$CAVES$sim_data_yearly[[name]]$ITPC[year] <- sum(DATA_SISAL$CAVES$sim_data_raw[[name]]$ISOT[pos_start:pos_stop]*
                                                                    DATA_SISAL$CAVES$sim_data_raw[[name]]$PREC[pos_start:pos_stop], na.rm = T)/
      sum(DATA_SISAL$CAVES$sim_data_raw[[name]]$PREC[pos_start:pos_stop], na.rm = T)
  }
}

remove(site_id, year, pos_start, pos_stop, name, ii)

PAST1000_SISAL <- DATA_SISAL$CAVES
save(PAST1000_SISAL, file = "Past1000_SISAL.RData")
