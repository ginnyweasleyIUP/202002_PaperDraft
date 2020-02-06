###################################################################################################
## CORRELATION ANALYSIS ###########################################################################
###################################################################################################

# ToDo:
# [ ] Correlation Maps von Temp-Iso, Prec - Iso 
#       [ ] Hierauf kann man dann Punkte legen, wie stark die Records an den einzelnen Punkten mit der Temp oder Prec korrelieren?
# [ ] Latitudinal Plot, wie stark die Iso Werte in der Simulation und die Iso-Werte in den Records übereinstimmen
# [ ]
# [ ] Frage: Wie sinnvoll sind die Point-Correlation Plots --> also wie viel Information kann man da raus ziehen?
# [ ]
# [ ] Sign Check!! But HOW???
# [ ]
# [ ]


library(dplyr)
library(latex2exp)
source("Functions/STACYmap_4.R")
library(PaleoSpec)

CORR_ANALYSIS <- list()

## 1) Get matrix analysis from simulation #######

CORR_ANALYSIS$GLOBAL_CORRELATION <- list(
  CORR_TEMP_PREC = array(dim = c(96,73)),
  CORR_TEMP_PREC_P = array(dim = c(96,73)),
  CORR_TEMP_ISOT = array(dim = c(96,73)),
  CORR_TEMP_ISOT_P = array(dim = c(96,73)),
  CORR_PREC_ISOT = array(dim = c(96,73)),
  CORR_PREC_ISOT_P = array(dim = c(96,73))
)

for (lon in 1:96){
  for (lat in 1:73){
    # TEMP-PREC
    COR_TP = cor.test(DATA_past1000$SIM_yearly$TEMP[lon,lat,], DATA_past1000$SIM_yearly$PREC[lon,lat,], na.rm = TRUE)
    CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_TEMP_PREC[lon,lat] <- COR_TP$estimate[[1]]
    CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_TEMP_PREC_P[lon,lat] <- COR_TP$p.value
    
    if(!any(is.na(DATA_past1000$SIM_yearly$ISOT[lon,lat,]))){
      COR_TI = cor.test(DATA_past1000$SIM_yearly$TEMP[lon,lat,], DATA_past1000$SIM_yearly$ISOT[lon,lat,], na.rm = TRUE)
      CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_TEMP_ISOT[lon,lat] = COR_TI$estimate[[1]]
      CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_TEMP_ISOT_P[lon,lat] = COR_TI$p.value
    }else{
      CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_TEMP_ISOT[lon,lat] = NA
      CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_TEMP_ISOT_P[lon,lat] = NA
    }
    
    
    if(!any(is.na(DATA_past1000$SIM_yearly$ISOT[lon,lat,]))){
      COR_PI = cor.test(DATA_past1000$SIM_yearly$PREC[lon,lat,], DATA_past1000$SIM_yearly$ISOT[lon,lat,], na.rm = TRUE)
      CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_PREC_ISOT[lon,lat] = COR_PI$estimate[[1]]
      CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_PREC_ISOT_P[lon,lat] = COR_PI$p.value
    }else{
      CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_PREC_ISOT[lon,lat] = NA
      CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_PREC_ISOT_P[lon,lat] = NA
    }
    
  }
}

remove(lon,lat, COR_TP, COR_TI, COR_PI)

## 2) Correlation of CAVE with Simualtion

#Consider all caves of last 1000 years. Caves with too little data points will fall out of the correlation because of missing significance. 

length_cave = length(DATA_past1000$CAVES$entity_info$entity_id)

CORR_ANALYSIS$CAVE_CORRELATION <- data.frame(
  entity_id = numeric(length_cave),
  CORR = numeric(length_cave),
  PVALUE = numeric(length_cave),
  CORR_TEMP = numeric(length_cave),
  PVALUE_TEMP = numeric(length_cave),
  CORR_PREC = numeric(length_cave),
  PVALUE_PREC = numeric(length_cave),
  CORR_pw = numeric(length_cave),
  PVALUE_pw = numeric(length_cave)
)

CORR_ANALYSIS$SITE_CORRELATION <- data.frame(
  entity_id = numeric(length_cave),
  CORR_TI = numeric(length_cave),
  PVALUE_TI = numeric(length_cave),
  CORR_PI = numeric(length_cave),
  PVALUE_PI = numeric(length_cave),
  CORR_TI_pw = numeric(length_cave),
  PVALUE_TI_pw = numeric(length_cave),
  CORR_PI_pw = numeric(length_cave),
  PVALUE_PI_pw = numeric(length_cave)
)


for(ii in 1:length_cave){
  print(ii)
  entity <- DATA_past1000$CAVES$entity_info$entity_id[ii]
  site <- DATA_past1000$CAVES$entity_info$site_id[ii]
  CORR_ANALYSIS$CAVE_CORRELATION$entity_id[ii] <- entity
  CORR_ANALYSIS$SITE_CORRELATION$entity_id[ii] <- entity
  # CAREFULL --> CORRELATION ONLY WORKS FOR EQUIDISTANT DATAPOINTS
  diff_dt = mean(diff(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age), na.rm = T)
  if(length(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age)>4 & ii != 95 & ii != 53 & ii != 109){
    #### SIM WITH RECORD
    record <- PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age,DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement,
                                         time.target = seq(from = head(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age, n = 1),
                                                           to = tail(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age, n = 1),
                                                           by = diff_dt))
    COR <- cor.test(record, PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$d18O,
                                                       time.target = seq(from = FirstElement(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age),
                                                                         to = LastElement(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age),
                                                                         by = diff_dt)))
    CORR_ANALYSIS$CAVE_CORRELATION$CORR[ii] = COR$estimate[[1]]
    CORR_ANALYSIS$CAVE_CORRELATION$PVALUE[ii] = COR$p.value
    
    COR_T <- cor.test(record, PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$Temp,
                                                         time.target = seq(from = FirstElement(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age),
                                                                           to = LastElement(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age),
                                                                           by = diff_dt)))
    CORR_ANALYSIS$CAVE_CORRELATION$CORR_TEMP[ii] = COR_T$estimate[[1]]
    CORR_ANALYSIS$CAVE_CORRELATION$PVALUE_TEMP[ii] = COR_T$p.value
    
    COR_P <- cor.test(record, PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$Prec,
                                                         time.target = seq(from = FirstElement(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age),
                                                                           to = LastElement(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age),
                                                                           by = diff_dt)))
    CORR_ANALYSIS$CAVE_CORRELATION$CORR_PREC[ii] = COR_P$estimate[[1]]
    CORR_ANALYSIS$CAVE_CORRELATION$PVALUE_PREC[ii] = COR_P$p.value
    
    COR_pw <- cor.test(record, PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$d18O_pw,
                                                          time.target = seq(from = FirstElement(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age),
                                                                            to = LastElement(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age),
                                                                            by = diff_dt)))
    CORR_ANALYSIS$CAVE_CORRELATION$CORR_pw[ii] = COR_pw$estimate[[1]]
    CORR_ANALYSIS$CAVE_CORRELATION$PVALUE_pw[ii] = COR_pw$p.value
    
  }else{
    CORR_ANALYSIS$CAVE_CORRELATION$CORR[ii] = NA
    CORR_ANALYSIS$CAVE_CORRELATION$PVALUE[ii] = NA
    CORR_ANALYSIS$CAVE_CORRELATION$CORR_TEMP[ii] = NA
    CORR_ANALYSIS$CAVE_CORRELATION$PVALUE_TEMP[ii] = NA
    CORR_ANALYSIS$CAVE_CORRELATION$CORR_PREC[ii] = NA
    CORR_ANALYSIS$CAVE_CORRELATION$PVALUE_PREC[ii] = NA
    CORR_ANALYSIS$CAVE_CORRELATION$CORR_pw[ii] = NA
    CORR_ANALYSIS$CAVE_CORRELATION$PVALUE_pw[ii] = NA
    
  }
  
  ###### SIM with SIM at Cave Site
  
  COR <- cor.test(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$TEMP, DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ISOT, na.rm = T)
  CORR_ANALYSIS$SITE_CORRELATION$CORR_TI[ii] <- COR$estimate[[1]]
  CORR_ANALYSIS$SITE_CORRELATION$PVALUE_TI[ii] <- COR$p.value
  
  COR <- cor.test(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$TEMP, DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ITPC, na.rm = T)
  CORR_ANALYSIS$SITE_CORRELATION$CORR_TI_pw[ii] <- COR$estimate[[1]]
  CORR_ANALYSIS$SITE_CORRELATION$PVALUE_TI_pw[ii] <- COR$p.value
  
  COR <- cor.test(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$PREC, DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ISOT, na.rm = T)
  CORR_ANALYSIS$SITE_CORRELATION$CORR_PI[ii] <- COR$estimate[[1]]
  CORR_ANALYSIS$SITE_CORRELATION$PVALUE_PI[ii] <- COR$p.value
  
  COR <- cor.test(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$TEMP, DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ITPC, na.rm = T)
  CORR_ANALYSIS$SITE_CORRELATION$CORR_PI_pw[ii] <- COR$estimate[[1]]
  CORR_ANALYSIS$SITE_CORRELATION$PVALUE_PI_pw[ii] <- COR$p.value

}

remove(ii, entity, site, diff_dt, COR)

#################################################
## Plot Maps: ###################################
#################################################

source("Functions/projection_ptlyr.R")

Plot_lyr <- CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_TEMP_ISOT
Plot_lyr_P <- CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_TEMP_ISOT_P
Plot_lyr[Plot_lyr_P > 0.1] <- NA
Plot_lyr[abs(Plot_lyr) < 0.4] <- NA

Plot_lyr <- rbind(Plot_lyr[49:96,1:73],
                  Plot_lyr[1:48,1:73])

Point_Lyr <- list(lon = list(), lat = list(), value = list())
Point2_Lyr <- list(lon = list(), lat = list(), value = list())

for(ii in 1:length_cave){
  site <- DATA_past1000$CAVES$entity_info$site_id[ii]
  if(is.na(CORR_ANALYSIS$CAVE_CORRELATION$PVALUE_TEMP[ii])){next}
  if(CORR_ANALYSIS$CAVE_CORRELATION$PVALUE_TEMP[ii] > 0.1){
    Point2_Lyr$lon = c(Point2_Lyr$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point2_Lyr$lat = c(Point2_Lyr$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point2_Lyr$value = c(Point2_Lyr$value, CORR_ANALYSIS$CAVE_CORRELATION$CORR_TEMP[ii])
  }else{
    
    Point_Lyr$lon = c(Point_Lyr$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point_Lyr$lat = c(Point_Lyr$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point_Lyr$value = c(Point_Lyr$value, CORR_ANALYSIS$CAVE_CORRELATION$CORR_TEMP[ii])
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
                 legend_names = list(grid = "Corr", pt = "")) +
  geom_point(data = Point2_Lyr, aes(x = long, y = lat), fill = 'gray', shape = 20,
             size = (GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE-1), show.legend = c(shape =TRUE))+
  geom_point(data = Point_Lyr, aes(x = long, y = lat, fill = layer), shape = 21,
             size = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE, show.legend = c(color = TRUE))+
  ggtitle(label = "Correlation ( c>0.4, p<0.1) for simulation T-d18O") +
  theme(plot.title = element_text(h = 0.5))

plot %>% ggsave(filename = paste('Map_Corr_TI_Records', 'pdf', sep = '.'), plot = ., path = 'Plots/Correlation', 
                width = 15, height = 10, units = 'cm', dpi = 'print', device = "pdf")

## Prec-ISOT

Plot_lyr <- CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_PREC_ISOT
Plot_lyr_P <- CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_PREC_ISOT_P
Plot_lyr[Plot_lyr_P > 0.1] <- NA
Plot_lyr[abs(Plot_lyr) < 0.4] <- NA

Point_Lyr <- list(lon = list(), lat = list(), value = list())
Point2_Lyr <- list(lon = list(), lat = list(), value = list())

for(ii in 1:length_cave){
  site <- DATA_past1000$CAVES$entity_info$site_id[ii]
  if(is.na(CORR_ANALYSIS$CAVE_CORRELATION$PVALUE_PREC[ii])){next}
  if(CORR_ANALYSIS$CAVE_CORRELATION$PVALUE_PREC[ii] > 0.1){
    Point2_Lyr$lon = c(Point2_Lyr$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point2_Lyr$lat = c(Point2_Lyr$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point2_Lyr$value = c(Point2_Lyr$value, CORR_ANALYSIS$CAVE_CORRELATION$CORR_PREC[ii])
  }else{
    
    Point_Lyr$lon = c(Point_Lyr$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point_Lyr$lat = c(Point_Lyr$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point_Lyr$value = c(Point_Lyr$value, CORR_ANALYSIS$CAVE_CORRELATION$CORR_PREC[ii])
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


Plot_lyr <- rbind(Plot_lyr[49:96,1:73],
                  Plot_lyr[1:48,1:73])

plot <- STACYmap(gridlyr = Plot_lyr,
                 centercolor = 0,
                 graticules = T,
                 legend_names = list(grid = "Corr", pt = "")) +
  geom_point(data = Point2_Lyr, aes(x = long, y = lat), fill = 'gray', shape = 20,
             size = (GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE-1), show.legend = c(shape =TRUE))+
  geom_point(data = Point_Lyr, aes(x = long, y = lat, fill = layer), shape = 21,
             size = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE, show.legend = c(color = TRUE))+
  ggtitle(label = "Correlation ( c>0.4, p<0.1) for simulation P-d18O") +
  theme(plot.title = element_text(h = 0.5))

plot %>% ggsave(filename = paste('Map_Corr_PI_Records', 'pdf', sep = '.'), plot = ., path = 'Plots/Correlation', 
                width = 15, height = 10, units = 'cm', dpi = 'print', device = "pdf")




###################################################################################################
## CORRELATION SIGN CHECK #########################################################################
###################################################################################################

CHECK_SIGN <- list(
  entity_id = numeric(length_cave),
  same_sign_TI = logical(length_cave),
  same_sign_PI = logical(length_cave),
  same_sign_TI_pw = logical(length_cave),
  same_sign_PI_pw = logical(length_cave)
)

CHECK_SIGN$entity_id <- CORR_ANALYSIS$CAVE_CORRELATION$entity_id

for(ii in 1:length_cave){
  if(!is.na(CORR_ANALYSIS$CAVE_CORRELATION$CORR_TEMP[ii])){
    if(sign(CORR_ANALYSIS$CAVE_CORRELATION$CORR_TEMP[ii]) == sign(CORR_ANALYSIS$SITE_CORRELATION$CORR_TI[ii])){CHECK_SIGN$same_sign_TI[ii] = T}
    if(sign(CORR_ANALYSIS$CAVE_CORRELATION$CORR_TEMP[ii]) == sign(CORR_ANALYSIS$SITE_CORRELATION$CORR_TI_pw[ii])){CHECK_SIGN$same_sign_TI_pw[ii] = T}
    if(sign(CORR_ANALYSIS$CAVE_CORRELATION$CORR_PREC[ii]) == sign(CORR_ANALYSIS$SITE_CORRELATION$CORR_PI[ii])){CHECK_SIGN$same_sign_PI[ii] = T}
    if(sign(CORR_ANALYSIS$CAVE_CORRELATION$CORR_PREC[ii]) == sign(CORR_ANALYSIS$SITE_CORRELATION$CORR_PI_pw[ii])){CHECK_SIGN$same_sign_PI_pw[ii] = T}
    
    if(CORR_ANALYSIS$CAVE_CORRELATION$PVALUE_TEMP[ii] > 0.1 | CORR_ANALYSIS$SITE_CORRELATION$PVALUE_TI[ii] > 0.1){CHECK_SIGN$same_sign_TI[ii] = NA}
    if (CORR_ANALYSIS$CAVE_CORRELATION$PVALUE_TEMP[ii] > 0.1 | CORR_ANALYSIS$SITE_CORRELATION$PVALUE_TI_pw[ii] > 0.1){CHECK_SIGN$same_sign_TI_pw[ii] = NA}
    if(CORR_ANALYSIS$CAVE_CORRELATION$PVALUE_PREC[ii] > 0.1 | CORR_ANALYSIS$SITE_CORRELATION$PVALUE_PI[ii] > 0.1){CHECK_SIGN$same_sign_PI[ii] = NA}
    if(CORR_ANALYSIS$CAVE_CORRELATION$PVALUE_PREC[ii] > 0.1 | CORR_ANALYSIS$SITE_CORRELATION$PVALUE_PI_pw[ii] > 0.1){CHECK_SIGN$same_sign_PI_pw[ii] = NA}
  }else{
      CHECK_SIGN$same_sign_PI[ii] <- NA
      CHECK_SIGN$same_sign_TI[ii] <- NA
      CHECK_SIGN$same_sign_TI_pw[ii] <- NA
      CHECK_SIGN$same_sign_PI_pw[ii] <- NA
    }
}

fields::image.plot(rbind(CHECK_SIGN$same_sign_TI, CHECK_SIGN$same_sign_TI_pw, CHECK_SIGN$same_sign_PI, CHECK_SIGN$same_sign_PI_pw), 
                   horizontal = T, nlevel = 134)



#?????????????????????????# 
# for (ii in entities_used){
#   name = paste0("ENTITY",ii)
#   diff_dt = c(diff_dt,mean(diff(DATA_last1000_mean[[name]]$value_record_age), na.rm = T))
# }
# diff_dt = mean(diff_dt, na.rm = T)
# 
# for (entity in entities_used){
#   site <- sample_final %>% filter(entity_id == entity) %>% count(site_id) %>% pull(site_id)
#   name = paste0("ENTITY", entity)
#   print(name)
#   prepare_corr_matrix_SIM       = c(prepare_corr_matrix_SIM,       CAVES$yearly_data$isot[[site]])
#   prepare_corr_matrix_SIM_ba_ed = c(prepare_corr_matrix_SIM_ba_ed, PaleoSpec::MakeEquidistant(DATA_last1000_mean[[name]]$value_record_age,
#                                                                                               DATA_last1000_mean[[name]]$value_sim_d18O,
#                                                                                               time.target = seq(from = -49, to = 1100, by = diff_dt)))
#   prepare_corr_matrix_REC_ed    = c(prepare_corr_matrix_REC_ed,    PaleoSpec::MakeEquidistant(DATA_last1000_mean[[name]]$value_record_age,
#                                                                                               DATA_last1000_mean[[name]]$value_record_d18O,
#                                                                                               time.target = seq(from = -49, to = 1100, by = diff_dt)))
# }
# 
# prepare_corr_matrix_SIM       <- matrix(prepare_corr_matrix_SIM,       nrow = 1150)
# prepare_corr_matrix_SIM_ba_ed <- matrix(prepare_corr_matrix_SIM_ba_ed, nrow = length(seq(from = -49, to = 1100, by = diff_dt)))
# prepare_corr_matrix_REC_ed    <- matrix(prepare_corr_matrix_REC_ed,    nrow = length(seq(from = -49, to = 1100, by = diff_dt)))
# 
# corr_matrix_SIM       <- Hmisc::rcorr(prepare_corr_matrix_SIM)
# corr_matrix_SIM_ba_ed <- Hmisc::rcorr(prepare_corr_matrix_SIM_ba_ed)
# corr_matrix_REC_ed    <- Hmisc::rcorr(prepare_corr_matrix_REC_ed)
# 
# remove(prepare_corr_matrix_SIM, prepare_corr_matrix_SIM_ba_ed, prepare_corr_matrix_REC_ed)