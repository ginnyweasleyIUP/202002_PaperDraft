#################################################
## Paper Figure 1 ###############################
#################################################

## Here analysis and Plotting

## Einführungsplot (GENERAL --> von Kira kopieren)

library(dplyr)
library(latex2exp)
source("Functions/STACYmap_5.R")
library(PaleoSpec)

ANALYSIS$CORR <- list()

#################################################
## CALCULATION ##################################
#################################################

# 1) FIELD (TEMP-ISOT and PREC-ISOT)

ANALYSIS$CORR$FIELD <- list(
  CORR_TEMP_ISOT = array(dim = c(96,73)),
  CORR_TEMP_ISOT_P = array(dim = c(96,73)),
  CORR_PREC_ISOT = array(dim = c(96,73)),
  CORR_PREC_ISOT_P = array(dim = c(96,73))
)

for (lon in 1:96){
  for (lat in 1:73){
   #TEMP ISOT
    if(!any(is.na(DATA_past1000$SIM_yearly$ISOT[lon,lat,]))){
      COR_TI = cor.test(DATA_past1000$SIM_yearly$TEMP[lon,lat,], DATA_past1000$SIM_yearly$ISOT[lon,lat,], na.rm = TRUE)
      ANALYSIS$CORR$FIELD$CORR_TEMP_ISOT[lon,lat] = COR_TI$estimate[[1]]
      ANALYSIS$CORR$FIELD$CORR_TEMP_ISOT_P[lon,lat] = COR_TI$p.value
    }else{
      ANALYSIS$CORR$FIELD$CORR_TEMP_ISOT[lon,lat] = NA
      ANALYSIS$CORR$FIELD$CORR_TEMP_ISOT_P[lon,lat] = NA
    }
    
    
    if(!any(is.na(DATA_past1000$SIM_yearly$ISOT[lon,lat,]))){
      COR_PI = cor.test(DATA_past1000$SIM_yearly$PREC[lon,lat,], DATA_past1000$SIM_yearly$ISOT[lon,lat,], na.rm = TRUE)
      ANALYSIS$CORR$FIELD$CORR_PREC_ISOT[lon,lat] = COR_PI$estimate[[1]]
      ANALYSIS$CORR$FIELD$CORR_PREC_ISOT_P[lon,lat] = COR_PI$p.value
    }else{
      ANALYSIS$CORR$FIELD$CORR_PREC_ISOT[lon,lat] = NA
      ANALYSIS$CORR$FIELD$CORR_PREC_ISOT_P[lon,lat] = NA
    }
    
  }
}

remove(lon,lat, COR_TP, COR_TI, COR_PI)


# 2) POINT (TEMP-d18O_dw_eq and PREC-d18O_dw_eq)

length_cave = length(DATA_past1000$CAVES$entity_info$entity_id)

ANALYSIS$CORR$POINTS <- data.frame(
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

for(ii in 1:length_cave){
  print(ii)
  entity <- DATA_past1000$CAVES$entity_info$entity_id[ii]
  site <- DATA_past1000$CAVES$entity_info$site_id[ii]
  ANALYSIS$CORR$POINTS$entity_id[ii] <- entity
  ANALYSIS$CORR$POINTS$entity_id[ii] <- entity
  # CAREFULL --> CORRELATION ONLY WORKS FOR EQUIDISTANT DATAPOINTS
  diff_dt = mean(diff(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age), na.rm = T)
  if(length(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age)>4 & ii != 95 & ii != 53 & ii != 109){
    #### SIM WITH RECORD
    record <- PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age,DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_dw_eq,
                                         time.target = seq(from = head(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age, n = 1),
                                                           to = tail(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age, n = 1),
                                                           by = diff_dt))
    COR <- cor.test(record, PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ISOT,
                                                       time.target = seq(from = FirstElement(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age),
                                                                         to = LastElement(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age),
                                                                         by = diff_dt)))
    ANALYSIS$CORR$POINTS$CORR[ii] = COR$estimate[[1]]
    ANALYSIS$CORR$POINTS$PVALUE[ii] = COR$p.value
    
    COR_T <- cor.test(record, PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$TEMP,
                                                         time.target = seq(from = FirstElement(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age),
                                                                           to = LastElement(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age),
                                                                           by = diff_dt)))
    ANALYSIS$CORR$POINTS$CORR_TEMP[ii] = COR_T$estimate[[1]]
    ANALYSIS$CORR$POINTS$PVALUE_TEMP[ii] = COR_T$p.value
    
    COR_P <- cor.test(record, PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$PREC,
                                                         time.target = seq(from = FirstElement(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age),
                                                                           to = LastElement(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age),
                                                                           by = diff_dt)))
    ANALYSIS$CORR$POINTS$CORR_PREC[ii] = COR_P$estimate[[1]]
    ANALYSIS$CORR$POINTS$PVALUE_PREC[ii] = COR_P$p.value
    
    COR_pw <- cor.test(record, PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ITPC,
                                                          time.target = seq(from = FirstElement(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age),
                                                                            to = LastElement(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age),
                                                                            by = diff_dt)))
    ANALYSIS$CORR$POINTS$CORR_pw[ii] = COR_pw$estimate[[1]]
    ANALYSIS$CORR$POINTS$PVALUE_pw[ii] = COR_pw$p.value
    
  }else{
    ANALYSIS$CORR$POINTS$CORR[ii] = NA
    ANALYSIS$CORR$POINTS$PVALUE[ii] = NA
    ANALYSIS$CORR$POINTS$CORR_TEMP[ii] = NA
    ANALYSIS$CORR$POINTS$PVALUE_TEMP[ii] = NA
    ANALYSIS$CORR$POINTS$CORR_PREC[ii] = NA
    ANALYSIS$CORR$POINTS$PVALUE_PREC[ii] = NA
    ANALYSIS$CORR$POINTS$CORR_pw[ii] = NA
    ANALYSIS$CORR$POINTS$PVALUE_pw[ii] = NA
    
  }
}


#################################################
## PLOTS ########################################
#################################################

Plot_lyr_temp <- ANALYSIS$CORR$FIELD$CORR_TEMP_ISOT
Plot_lyr_temp_p <- ANALYSIS$CORR$FIELD$CORR_TEMP_ISOT_P
Plot_lyr_temp[Plot_lyr_temp_p > 0.1] <- NA
Plot_lyr_temp[abs(Plot_lyr_temp) < 0.2] <- NA
Plot_lyr_prec <- ANALYSIS$CORR$FIELD$CORR_PREC_ISOT
Plot_lyr_prec_p <- ANALYSIS$CORR$FIELD$CORR_PREC_ISOT_P
Plot_lyr_prec[Plot_lyr_prec_p > 0.1] <- NA
Plot_lyr_prec[abs(Plot_lyr_prec) < 0.2] <- NA

Plot_lyr_temp <- rbind(Plot_lyr_temp[49:96,1:73],
                  Plot_lyr_temp[1:48,1:73])
Plot_lyr_prec <- rbind(Plot_lyr_prec[49:96,1:73],
                       Plot_lyr_prec[1:48,1:73])

##### Point Layer

Point_Lyr_temp <- list(lon = list(), lat = list(), value = list())
Point_Lyr_prec <- list(lon = list(), lat = list(), value = list())

length_cave = length(DATA_past1000$CAVES$entity_info$site_id)

for(ii in 1:length_cave){
  site <- DATA_past1000$CAVES$entity_info$site_id[ii]
  print(ii)
  if(is.na(ANALYSIS$CORR$POINTS$PVALUE_TEMP[ii])){next}
  # 1) sortiert aus, was nicht signifikant ist
  if(ANALYSIS$CORR$POINTS$PVALUE_TEMP[ii] > 0.1){
    Point_Lyr_temp$lon = c(Point_Lyr_temp$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point_Lyr_temp$lat = c(Point_Lyr_temp$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point_Lyr_temp$value = c(Point_Lyr_temp$value, ANALYSIS$CORR$POINTS$CORR_TEMP[ii])
    # 2) betrachte signifikante Korrelationen:
  }
  if(is.na(ANALYSIS$CORR$POINTS$PVALUE_PREC[ii])){next}
  # 1) sortiert aus, was nicht signifikant ist
  if(ANALYSIS$CORR$POINTS$PVALUE_TEMP[ii] > 0.1){
    Point_Lyr_prec$lon = c(Point_Lyr_prec$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point_Lyr_prec$lat = c(Point_Lyr_prec$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point_Lyr_prec$value = c(Point_Lyr_prec$value, ANALYSIS$CORR$POINTS$CORR_PREC[ii])
    # 2) betrachte signifikante Korrelationen:
  }
}



Point_Lyr_temp$lon = as.numeric(Point_Lyr_temp$lon)
Point_Lyr_temp$lat = as.numeric(Point_Lyr_temp$lat)
Point_Lyr_temp$value = as.numeric(Point_Lyr_temp$value)

Point_Lyr_prec$lon = as.numeric(Point_Lyr_prec$lon)
Point_Lyr_prec$lat = as.numeric(Point_Lyr_prec$lat)
Point_Lyr_prec$value = as.numeric(Point_Lyr_prec$value)


GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE <- 3

plot_temp <- STACYmap(gridlyr = Plot_lyr_temp, centercolor = 0, graticules = T,
                 ptlyr = as.data.frame(Point_Lyr_temp), legend_names = list(grid = 'Temp.-Correlation (p<0.1)')) + 
  theme(panel.border = element_blank(),
        legend.background = element_blank(),
        axis.text = element_blank(),
        text = element_text(size = 12),
        legend.title = element_text(size = 12))

plot_prec <- STACYmap(gridlyr = Plot_lyr_prec, centercolor = 0, graticules = T,
                      ptlyr = as.data.frame(Point_Lyr_prec), legend_names = list(grid = 'Prec.-Correlation (p<0.1)')) + 
  theme(panel.border = element_blank(),
        legend.background = element_blank(),
        axis.text = element_blank(),
        text = element_text(size = 12),
        legend.title = element_text(size = 12))

library(ggpubr)
plot <- ggarrange(plot_temp, plot_prec,
                  labels = c("A", "B"),
                  ncol = 2, nrow = 1)

plot  %>% ggsave(filename = paste('Paper_Plot_5_Correlation', 'pdf', sep = '.'), plot = ., path = 'Plots', 
                 width = 2*12, height = 12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")


#################################################
## Here the all in all Plot #####################
#################################################
# source("Functions/projection_ptlyr.R")
# # Grid Layer for plotting:
# # all areas where d18O correlates better with temperature are marked in red
# # all areas where d18O correlates better with precipitation are marked in blue
# Plot_lyr_temp <- CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_TEMP_ISOT
# Plot_lyr_temp_p <- CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_TEMP_ISOT_P
# Plot_lyr_prec <- CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_PREC_ISOT
# Plot_lyr_prec_p <- CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_PREC_ISOT_P
# Plot_lyr_temp[Plot_lyr_temp_p > 0.1] <- 0
# Plot_lyr_temp[abs(Plot_lyr_temp) < 0.2] <- 0
# Plot_lyr_prec[Plot_lyr_prec_p > 0.1] <- 0
# Plot_lyr_prec[abs(Plot_lyr_prec) < 0.2] <- 0
# 
# Plot_lyr_2 <- Plot_lyr_temp
# Plot_lyr_3 <- Plot_lyr_prec
# 
# Plot_lyr_2[abs(Plot_lyr_prec)>abs(Plot_lyr_temp)] <- 0
# Plot_lyr_3[abs(Plot_lyr_temp)>abs(Plot_lyr_prec)] <- 0
# 
# Plot_lyr <- abs(Plot_lyr_2)- abs(Plot_lyr_3)
# Plot_lyr[Plot_lyr == 0] <- NA
# 
# Plot_lyr <- rbind(Plot_lyr[49:96,1:73],
#                   Plot_lyr[1:48,1:73])
# 
# remove(Plot_lyr_2, Plot_lyr_3, Plot_lyr_prec, Plot_lyr_prec_p, Plot_lyr_temp, Plot_lyr_temp_p)
# 
# ##### Point Layer
# 
# # How should points be colored? Is it so relevant if sign is equal?
# 
# # 0) Check for significance --> if not then, then put in Point_lyr_2
# # 1) Check for what the absolute corellation is stronger
# # 2) make different shapes depending on sign fitting or not
# 
# 
# ### HERE HERE HERE ############################
# ## es muss noch angepasst werden, dass alle Punktlisten mit unterschiedlichem Symbol über eine andere Liste gemacht wird. 
# Point_Lyr_sign <- list(lon = list(), lat = list(), value = list())
# Point_Lyr_notsign <- list(lon = list(), lat = list(), value = list())
# Point2_Lyr <- list(lon = list(), lat = list(), value = list())
# 
# length_cave = length(DATA_past1000$CAVES$entity_info$site_id)
# 
# for(ii in 1:length_cave){
#   site <- DATA_past1000$CAVES$entity_info$site_id[ii]
#   print(ii)
#   if(is.na(ANALYSIS$CORR$POINTS$PVALUE_TEMP[ii])){next}
#   # 1) sortiert aus, was nicht signifikant ist
#   if(ANALYSIS$CORR$POINTS$PVALUE_TEMP[ii] > 0.1 & ANALYSIS$CORR$POINTS$PVALUE_PREC[ii] > 0.1){
#     Point2_Lyr$lon = c(Point2_Lyr$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
#     Point2_Lyr$lat = c(Point2_Lyr$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
#     Point2_Lyr$value = c(Point2_Lyr$value, ANALYSIS$CORR$POINTS$CORR_TEMP[ii])
#     # 2) betrachte signifikante Korrelationen:
#   }else{
#     # 2.1) Nur signifikante Korrelation bei Temp
#     if(ANALYSIS$CORR$POINTS$PVALUE_TEMP[ii] < 0.1 & ANALYSIS$CORR$POINTS$PVALUE_PREC[ii] > 0.1){
#       #Check sign to determine shape
#       if(sign(ANALYSIS$CORR$POINTS$CORR_TEMP[ii]) == sign(CORR_ANALYSIS$SITE_CORRELATION$CORR_TI[ii])){
#         Point_Lyr_sign$lon = c(Point_Lyr_sign$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
#         Point_Lyr_sign$lat = c(Point_Lyr_sign$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
#         Point_Lyr_sign$value = c(Point_Lyr_sign$value, abs(ANALYSIS$CORR$POINTS$CORR_TEMP[ii]))
#       }else{
#         Point_Lyr_notsign$lon = c(Point_Lyr_notsign$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
#         Point_Lyr_notsign$lat = c(Point_Lyr_notsign$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
#         Point_Lyr_notsign$value = c(Point_Lyr_notsign$value, abs(ANALYSIS$CORR$POINTS$CORR_TEMP[ii]))
#       }
#     }
#     
#     # 2.2) Nur signifikante Korrelation bei Prec
#     else if(ANALYSIS$CORR$POINTS$PVALUE_TEMP[ii] > 0.1 & ANALYSIS$CORR$POINTS$PVALUE_PREC[ii] < 0.1){
#       if(sign(ANALYSIS$CORR$POINTS$CORR_PREC[ii]) == sign(CORR_ANALYSIS$SITE_CORRELATION$CORR_PI[ii])){
#         Point_Lyr_sign$lon = c(Point_Lyr_sign$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
#         Point_Lyr_sign$lat = c(Point_Lyr_sign$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
#         Point_Lyr_sign$value = c(Point_Lyr_sign$value, - abs(ANALYSIS$CORR$POINTS$CORR_PREC[ii]))
#       }else{
#         Point_Lyr_notsign$lon = c(Point_Lyr_notsign$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
#         Point_Lyr_notsign$lat = c(Point_Lyr_notsign$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
#         Point_Lyr_notsign$value = c(Point_Lyr_notsign$value, - abs(ANALYSIS$CORR$POINTS$CORR_PREC[ii]))
#       }
#     }
#     
#     # 2.3) Sowohl signifikant für Prec wie für Temp    
#     else{
#       # 2.3.1) absolute CORR größer für Temp als für Prec
#       if(abs(ANALYSIS$CORR$POINTS$CORR_TEMP[ii]) > abs(ANALYSIS$CORR$POINTS$CORR_PREC[ii])){
#         if(sign(ANALYSIS$CORR$POINTS$CORR_TEMP[ii]) == sign(CORR_ANALYSIS$SITE_CORRELATION$CORR_TI[ii])){
#           Point_Lyr_sign$lon = c(Point_Lyr_sign$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
#           Point_Lyr_sign$lat = c(Point_Lyr_sign$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
#           Point_Lyr_sign$value = c(Point_Lyr_sign$value, abs(ANALYSIS$CORR$POINTS$CORR_TEMP[ii]))
#         }else{
#           Point_Lyr_notsign$lon = c(Point_Lyr_notsign$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
#           Point_Lyr_notsign$lat = c(Point_Lyr_notsign$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
#           Point_Lyr_notsign$value = c(Point_Lyr_notsign$value, abs(ANALYSIS$CORR$POINTS$CORR_TEMP[ii]))}
#       }
#       # 2.3.2) absolute CORR größer für Prec als für Temp
#       else{
#         if(sign(ANALYSIS$CORR$POINTS$CORR_PREC[ii]) == sign(CORR_ANALYSIS$SITE_CORRELATION$CORR_PI[ii])){
#           Point_Lyr_sign$lon = c(Point_Lyr_sign$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
#           Point_Lyr_sign$lat = c(Point_Lyr_sign$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
#           Point_Lyr_sign$value = c(Point_Lyr_sign$value, - abs(ANALYSIS$CORR$POINTS$CORR_PREC[ii]))
#         }else{
#           Point_Lyr_notsign$lon = c(Point_Lyr_notsign$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
#           Point_Lyr_notsign$lat = c(Point_Lyr_notsign$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
#           Point_Lyr_notsign$value = c(Point_Lyr_notsign$value, - abs(ANALYSIS$CORR$POINTS$CORR_PREC[ii]))}
#       }
#     }
#   }
# }
# 
# 
# 
# Point_Lyr_sign$lon = as.numeric(Point_Lyr_sign$lon)
# Point_Lyr_sign$lat = as.numeric(Point_Lyr_sign$lat)
# Point_Lyr_sign$value = as.numeric(Point_Lyr_sign$value)
# 
# Point_Lyr_notsign$lon = as.numeric(Point_Lyr_notsign$lon)
# Point_Lyr_notsign$lat = as.numeric(Point_Lyr_notsign$lat)
# Point_Lyr_notsign$value = as.numeric(Point_Lyr_notsign$value)
# 
# Point2_Lyr$lon = as.numeric(Point2_Lyr$lon)
# Point2_Lyr$lat = as.numeric(Point2_Lyr$lat)
# Point2_Lyr$value = as.numeric(Point2_Lyr$value)
# 
# 
# 
# GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE <- 3
# 
# Point_Lyr_sign_p <-  projection_ptlyr(as.data.frame(Point_Lyr_sign), as.character('+proj=robin +datum=WGS84'))
# Point_Lyr_notsign_p <-  projection_ptlyr(as.data.frame(Point_Lyr_notsign), as.character('+proj=robin +datum=WGS84'))
# Point2_Lyr_p <-  projection_ptlyr(as.data.frame(Point2_Lyr), as.character('+proj=robin +datum=WGS84'))
# 
# remove(Point_Lyr_sign, Point_Lyr_notsign, Point2_Lyr)
# 
# # Jetzt existiert ein Plot Layer und 2 Point Layer die man nur noch plotten muss und eine richtige Legende dafür braucht...
# 
# source("Functions/STACYmap_6.R")
# source("Functions/STACYmap_5_2_logscale_corr.R")
# 
# plot <- STACYmap_isot_corr(gridlyr = Plot_lyr, centercolor = 0, graticules = T, 
#                            legend_names = list(grid = "abs(Corr.)"),
#                            breaks_isot = c(-1, -0.5, 0, 0.51, 1),
#                            labels_isot = c(1, "corr prec", "0",  "corr temp", 1)) + 
#   geom_point(data = Point2_Lyr_p, aes(x = long, y = lat, shape = "1"), fill = 'gray', size = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE-1, show.legend = c(shape = T)) + 
#   geom_point(data = Point_Lyr_sign_p, aes(x = long, y = lat, fill = layer, shape = "2"), size = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE, show.legend = c(color = T, shape = T)) + 
#   geom_point(data = Point_Lyr_notsign_p, aes(x = long, y = lat, fill = layer, shape = "3"), size = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE, show.legend = c(color = T, shape = T)) + 
#   scale_shape_manual(name = NULL, labels = c("no corr.-sites", "same sign", "different sign"), 
#                      values = c(20,21,23))+
#   #guides(fill = guide_colorbar(label = F, direction = "horizontal", title = "|Corr.| blue prec, red temp")) + 
#   theme(panel.border = element_blank(),
#         legend.background = element_blank(),
#         axis.text = element_blank(),
#         text = element_text(size = 12),
#         legend.title = element_text(size = 12))
# 
# plot
# 
# plot %>% ggsave(filename = paste('Paper_Plot_5_Correlation', 'pdf', sep = '.'), plot = ., path = 'Plots/Paper', 
#                 width = 2*PLOTTING_VARIABLES$WIDTH, height = 2*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")
# 
