#################################################
## ANALYSIS NETWORKS ############################
#################################################

library(plyr)
library(dplyr)
library(tidyverse)
library(zoo)
library(PaleoSpec)
library(nest)
library(maps)

#time_caves <- seq(from = -49, to = 1100, by = 1)
# blabla

ANALYSIS$NETWORK <- list()

## 1 Old Plot ###################################

lats = c()
longs = c()

#counter = 1
for (entity in DATA_past1000$CAVES$entity_info$entity_id[mask_spec]){
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  #prep_corr_matrix[,counter] = DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ISOT
  #counter = counter + 1
  lats = c(lats, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
  longs = c(longs, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
}

ANALYSIS$NETWORK$entity_meta <- data.frame(
  entity_id = DATA_past1000$CAVES$entity_info$entity_id[mask_spec],
  site_id = DATA_past1000$CAVES$entity_info$site_id[mask_spec],
  lat = lats,
  long = longs,
  gridbox_id = numeric(sum(mask_spec)),
  gridbox_lon = numeric(sum(mask_spec)),
  gridbox_lat = numeric(sum(mask_spec)),
  cluster_id = numeric(sum(mask_spec))
)


dist<-fossil::earth.dist(cbind(lats,longs),dist=TRUE)
dist_matrix <- as.matrix(dist)

ANALYSIS$NETWORK$DIST <- dist_matrix

remove(lats, longs)

## 2 new Plot ###################################

## 2 0 assign gridbox and cluster

# Clustering

hc<-hclust(dist)
plot(hc)

ANALYSIS$hc <- hc

cluster <-cutree(hc,k=8)
#Point_Lyr <- list(
#  lat = ANALYSIS$NETWORK$entity_meta$lat,
#  long = ANALYSIS$NETWORK$entity_meta$long,
#  value =cluster
#)

#map(database='world')
#points(Point_Lyr$long,Point_Lyr$lat,col=Point_Lyr$value+1,pch=3)

#asigning

for(ii in 1:sum(mask_spec)){
  e.long = ANALYSIS$NETWORK$entity_meta$long[ii]
  e.lat = ANALYSIS$NETWORK$entity_meta$lat[ii]
  
  ANALYSIS$NETWORK$entity_meta$gridbox_lon[ii] = ceiling((e.long+180)/3.75)
  ANALYSIS$NETWORK$entity_meta$gridbox_lat[ii] = ceiling((-1*e.lat+90)/180*73)
  ANALYSIS$NETWORK$entity_meta$gridbox_id[ii]  = (ceiling((-1*e.lat+90)/180*73)-1)*96 + ceiling((e.long+180)/3.75)
  ANALYSIS$NETWORK$entity_meta$cluster_id[ii] = cluster[ii]
}

rm(e.lat, e.long, entity, ii, dist, cluster, site)

#################################################
## CAVE CORRELATION #############################
#################################################


#gaussdetr(zoo(x = Timeseries$pages2k$value[850:2000],
#              order.by = Timeseries$pages2k$time[850:2000]), tsc.in = 100)$Xsmooth


#For cross-correlation we use the d18O measurements, as it does not include simulation data

## 2 1 site correlation 

#Create list with sites that have multiple entities and which entities

site_list <- DATA_past1000$CAVES$entity_info %>% select(site_id, entity_id) %>%
  filter(entity_id %in% ANALYSIS$NETWORK$entity_meta$entity_id) %>% group_by(site_id) %>% count() %>% filter(n>1)

# Calculation

ANALYSIS$NETWORK$SITES <- list()
ANALYSIS$NETWORK$SITES$mean <- numeric(length(site_list$site_id))
ANALYSIS$NETWORK$SITES$range <- list()
ANALYSIS$NETWORK$SITES$mean_100gauss <- numeric(length(site_list$site_id))
ANALYSIS$NETWORK$SITES$range_100gauss <- list()
counter = 1

for(site in site_list$site_id){
  
  print(paste("site", site))

  ANALYSIS$NETWORK$SITES[[paste0("SITE", site)]] <- list()
  
  entity_list <- DATA_past1000$CAVES$entity_info %>% select(site_id, entity_id) %>% filter(entity_id %in% ANALYSIS$NETWORK$entity_meta$entity_id) %>%
    filter(site_id == site) %>% select(entity_id)
  
  TS <- list()
  TS_gauss <- list()
  
  for(entity in entity_list$entity_id){
    s <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
    # zoo cannot handle objects where order.by has two elements which is why they are sorted out here (no better option found)
    double_time <- s %>% group_by(interp_age) %>% count() %>% filter(n>1)
    s <- s %>% filter(!interp_age %in% double_time$interp_age) %>% filter(!is.na(d18O_measurement))
    
    TS[[paste0("Entity", entity)]] <- zoo(x = s$d18O_measurement, order.by = s$interp_age)
    TS_gauss[[paste0("Entity", entity)]] <- gaussdetr(zoo(x = s$d18O_measurement, order.by = s$interp_age), tsc.in = 100)$Xsmooth
  }
  
  C<-matrix(NA,nrow=length(TS),ncol=length(TS))
  colnames(C)<-rownames(C)<-names(TS)
  C_gauss <- P_gauss <- P <- C
  
  for (i in 1:(length(TS)-1)){
    for (j in (i+1):length(TS)){
      temp<-nest::nexcf_ci(TS[[i]],TS[[j]],conflevel=0.1)
      temp_gauss<-nest::nexcf_ci(TS_gauss[[i]],TS_gauss[[j]],conflevel=0.1)
      C[i,j]<-temp$rxy
      P[i,j]<-P[j,i]<-temp$pval
      C[j,i]=C[i,j]
      rm(temp)
      C_gauss[i,j]<-temp_gauss$rxy
      P_gauss[i,j]<-P_gauss[j,i]<-temp_gauss$pval
      C_gauss[j,i]=C_gauss[i,j]
      rm(temp_gauss)
    }
  }
  #corrplot::corrplot(C)
  
  ANALYSIS$NETWORK$SITES[[paste0("SITE", site)]]$C <- C
  ANALYSIS$NETWORK$SITES[[paste0("SITE", site)]]$P <- P
  ANALYSIS$NETWORK$SITES$mean[counter] <- mean(C, na.rm = T)
  ANALYSIS$NETWORK$SITES$range[counter] <- range(C, na.rm = T)
  ANALYSIS$NETWORK$SITES[[paste0("SITE", site)]]$C_gauss <- C_gauss
  ANALYSIS$NETWORK$SITES[[paste0("SITE", site)]]$P_gauss <- P_gauss
  ANALYSIS$NETWORK$SITES$mean_100gauss[counter] <- mean(C_gauss, na.rm = T)
  ANALYSIS$NETWORK$SITES$range_100gauss[counter] <- range(C_gauss, na.rm = T)
  counter = counter + 1
  
}

## 2 2 gridbox correlation

gridbox_list <- ANALYSIS$NETWORK$entity_meta %>% select(gridbox_id, entity_id) %>% group_by(gridbox_id) %>% count() %>% filter(n>1)
  
# Calculation

ANALYSIS$NETWORK$GRIDBOX <- list()
ANALYSIS$NETWORK$GRIDBOX$mean <- numeric(length(gridbox_list$gridbox_id))
ANALYSIS$NETWORK$GRIDBOX$range <- list()
ANALYSIS$NETWORK$GRIDBOX$mean_100gauss <- numeric(length(gridbox_list$gridbox_id))
ANALYSIS$NETWORK$GRIDBOX$range_100gauss <- list()
counter = 1

for(gridbox in gridbox_list$gridbox_id){
  
  print(paste("gridbox", gridbox))
  
  ANALYSIS$NETWORK$GRIDBOX[[paste0("gridbox", gridbox)]] <- list()
  
  entity_list <- ANALYSIS$NETWORK$entity_meta %>% select(gridbox_id, entity_id) %>% filter(gridbox_id == gridbox)
  
  TS <- list()
  TS_gauss <- list()
  
  for(entity in entity_list$entity_id){
    s <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
    # zoo cannot handle objects where order.by has two elements which is why they are sorted out here (no better option found)
    double_time <- s %>% group_by(interp_age) %>% count() %>% filter(n>1)
    s <- s %>% filter(!interp_age %in% double_time$interp_age) %>% filter(!is.na(d18O_measurement))
    
    TS[[paste0("Entity", entity)]] <- zoo(x = s$d18O_measurement, order.by = s$interp_age)
    TS_gauss[[paste0("Entity", entity)]] <- gaussdetr(zoo(x = s$d18O_measurement, order.by = s$interp_age), tsc.in = 100)$Xsmooth
  }
  
  C<-matrix(NA,nrow=length(TS),ncol=length(TS))
  colnames(C)<-rownames(C)<-names(TS)
  C_gauss <- P_gauss <- P <- C
  
  for (i in 1:(length(TS)-1)){
    for (j in (i+1):length(TS)){
      temp<-nest::nexcf_ci(TS[[i]],TS[[j]],conflevel=0.1)
      temp_gauss<-nest::nexcf_ci(TS_gauss[[i]],TS_gauss[[j]],conflevel=0.1)
      C[i,j]<-temp$rxy
      P[i,j]<-P[j,i]<-temp$pval
      C[j,i]=C[i,j]
      rm(temp)
      C_gauss[i,j]<-temp_gauss$rxy
      P_gauss[i,j]<-P_gauss[j,i]<-temp_gauss$pval
      C_gauss[j,i]=C_gauss[i,j]
      rm(temp_gauss)
    }
  }

  ANALYSIS$NETWORK$GRIDBOX[[paste0("GRIDBOX", gridbox)]]$C <- C
  ANALYSIS$NETWORK$GRIDBOX[[paste0("GRIDBOX", gridbox)]]$P <- P
  ANALYSIS$NETWORK$GRIDBOX$mean[counter] <- mean(C, na.rm = T)
  ANALYSIS$NETWORK$GRIDBOX$range[counter] <- range(C, na.rm = T)
  ANALYSIS$NETWORK$GRIDBOX[[paste0("GRIDBOX", gridbox)]]$C_gauss <- C_gauss
  ANALYSIS$NETWORK$GRIDBOX[[paste0("GRIDBOX", gridbox)]]$P_gauss <- P_gauss
  ANALYSIS$NETWORK$GRIDBOX$mean_100gauss[counter] <- mean(C_gauss, na.rm = T)
  ANALYSIS$NETWORK$GRIDBOX$range_100gauss[counter] <- range(C_gauss, na.rm = T)
  counter = counter + 1
  
}

## 2 3 cluster correlation

cluster_list <- ANALYSIS$NETWORK$entity_meta %>% select(cluster_id, entity_id) %>% group_by(cluster_id) %>% count() %>% filter(n>1)

# Calculation

ANALYSIS$NETWORK$CLUSTER <- list()
ANALYSIS$NETWORK$CLUSTER$mean <- numeric(length(cluster_list$cluster_id))
ANALYSIS$NETWORK$CLUSTER$range <- list()
ANALYSIS$NETWORK$CLUSTER$mean_100gauss <- numeric(length(cluster_list$cluster_id))
ANALYSIS$NETWORK$CLUSTER$range_100gauss <- list()
counter = 1

for(cluster in cluster_list$cluster_id){
  
  print(paste("cluster", cluster))
  
  ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]] <- list()
  
  entity_list <- ANALYSIS$NETWORK$entity_meta %>% select(cluster_id, entity_id) %>% filter(cluster_id == cluster)
  
  TS <- list()
  TS_gauss <- list()
  
  for(entity in entity_list$entity_id){
    s <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
    # zoo cannot handle objects where order.by has two elements which is why they are sorted out here (no better option found)
    double_time <- s %>% group_by(interp_age) %>% count() %>% filter(n>1)
    s <- s %>% filter(!interp_age %in% double_time$interp_age) %>% filter(!is.na(d18O_measurement))
    
    TS[[paste0("Entity", entity)]] <- zoo(x = s$d18O_measurement, order.by = s$interp_age)
    TS_gauss[[paste0("Entity", entity)]] <- gaussdetr(zoo(x = s$d18O_measurement, order.by = s$interp_age), tsc.in = 100)$Xsmooth
  }
  
  C<-matrix(NA,nrow=length(TS),ncol=length(TS))
  colnames(C)<-rownames(C)<-names(TS)
  C_gauss <- P_gauss <- P <- C
  
  for (i in 1:(length(TS)-1)){
    for (j in (i+1):length(TS)){
      temp<-nest::nexcf_ci(TS[[i]],TS[[j]],conflevel=0.1)
      temp_gauss<-nest::nexcf_ci(TS_gauss[[i]],TS_gauss[[j]],conflevel=0.1)
      C[i,j]<-temp$rxy
      P[i,j]<-P[j,i]<-temp$pval
      C[j,i]=C[i,j]
      rm(temp)
      C_gauss[i,j]<-temp_gauss$rxy
      P_gauss[i,j]<-P_gauss[j,i]<-temp_gauss$pval
      C_gauss[j,i]=C_gauss[i,j]
      rm(temp_gauss)
    }
  }
  
  
  ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]]$C <- C
  ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]]$P <- P
  ANALYSIS$NETWORK$CLUSTER$mean[counter] <- mean(C, na.rm = T)
  ANALYSIS$NETWORK$CLUSTER$range[counter] <- range(C, na.rm = T)
  ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]]$C_gauss <- C_gauss
  ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]]$P_gauss <- P_gauss
  ANALYSIS$NETWORK$CLUSTER$mean_100gauss[counter] <- mean(C_gauss, na.rm = T)
  ANALYSIS$NETWORK$CLUSTER$range_100gauss[counter] <- range(C_gauss, na.rm = T)
  counter = counter + 1
  
}

## 2 4 global correlation

ANALYSIS$NETWORK$GLOBAL <- list()
ANALYSIS$NETWORK$GLOBAL$mean <- list()
ANALYSIS$NETWORK$GLOBAL$range <- list()

entity_list <- ANALYSIS$NETWORK$entity_meta$entity_id

  
TS <- list()
TS_gauss <- list()
  
for(entity in entity_list){
  s <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
  # zoo cannot handle objects where order.by has two elements which is why they are sorted out here (no better option found)
  double_time <- s %>% group_by(interp_age) %>% count() %>% filter(n>1)
  s <- s %>% filter(!interp_age %in% double_time$interp_age) %>% filter(!is.na(d18O_measurement))
  
  TS[[paste0("Entity", entity)]] <- zoo(x = s$d18O_measurement, order.by = s$interp_age)
  TS_gauss[[paste0("Entity", entity)]] <- gaussdetr(zoo(x = s$d18O_measurement, order.by = s$interp_age), tsc.in = 100)$Xsmooth
}
  
C<-matrix(NA,nrow=length(TS),ncol=length(TS))
colnames(C)<-rownames(C)<-names(TS)
C_gauss <- P_gauss <- P <- C

for (i in 1:(length(TS)-1)){
  for (j in (i+1):length(TS)){
    temp<-nest::nexcf_ci(TS[[i]],TS[[j]],conflevel=0.1)
    temp_gauss<-nest::nexcf_ci(TS_gauss[[i]],TS_gauss[[j]],conflevel=0.1)
    C[i,j]<-temp$rxy
    P[i,j]<-P[j,i]<-temp$pval
    C[j,i]=C[i,j]
    rm(temp)
    C_gauss[i,j]<-temp_gauss$rxy
    P_gauss[i,j]<-P_gauss[j,i]<-temp_gauss$pval
    C_gauss[j,i]=C_gauss[i,j]
    rm(temp_gauss)
  }
}

ANALYSIS$NETWORK$GLOBAL$C <- C
ANALYSIS$NETWORK$GLOBAL$P <- P
ANALYSIS$NETWORK$GLOBAL$mean <- mean(C, na.rm = T)
ANALYSIS$NETWORK$GLOBAL$range <- range(C, na.rm = T)
ANALYSIS$NETWORK$GLOBAL$C_gauss <- C_gauss
ANALYSIS$NETWORK$GLOBAL$P_gauss <- P_gauss
ANALYSIS$NETWORK$GLOBAL$mean_100gauss <- mean(C, na.rm = T)
ANALYSIS$NETWORK$GLOBAL$range_100gauss <- range(C, na.rm = T)

remove(cluster, counter, entity, gridbox, i, j, site, TS, site_list, s, P, hc, gridbox_list, entity_list, double_time, cluster_list, C)

#################################################
## SIMULATION ###################################
#################################################

for(run in c("a","b","c")){
  ## 2 3 cluster correlation

  cluster_list <- ANALYSIS$NETWORK$entity_meta %>% select(cluster_id, entity_id) %>% group_by(cluster_id) %>% count() %>% filter(n>1)

  # Calculation

  ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_",run)]] <- list()
  ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_",run)]]$mean <- numeric(length(cluster_list$cluster_id))
  ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_",run)]]$range <- list()
  ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_",run)]]$mean_100gauss <- numeric(length(cluster_list$cluster_id))
  ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_",run)]]$range_100gauss <- list()
  counter = 1

  for(cluster in cluster_list$cluster_id){

    print(paste("cluster", cluster))

    ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_",run)]][[paste0("CLUSTER", cluster)]] <- list()

    entity_list <- ANALYSIS$NETWORK$entity_meta %>% select(cluster_id, entity_id) %>% filter(cluster_id == cluster)

    TS <- list()
    TS_gauss <- list()

    for(entity in entity_list$entity_id){
      site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
      TS[[paste0("Entity", entity)]] <- zoo(x = DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]][[paste0("ISOT_", run)]],
                                            order.by = seq(from = 1950-DATA_past1000$time[1], to = 1950-DATA_past1000$time[2], by = -1))
      TS_gauss[[paste0("Entity", entity)]] <- gaussdetr(zoo(x = DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]][[paste0("ISOT_", run)]],
                                                      order.by = seq(from = 1950-DATA_past1000$time[1], to = 1950-DATA_past1000$time[2], by = -1)), tsc.in = 100)$Xsmooth
    }

    C<-matrix(NA,nrow=length(TS),ncol=length(TS))
    colnames(C)<-rownames(C)<-names(TS)
    C_gauss <- P_gauss <- P <- C

    for (i in 1:(length(TS)-1)){
      for (j in (i+1):length(TS)){
        temp<-nest::nexcf_ci(TS[[i]],TS[[j]],conflevel=0.1)
        temp_gauss<-nest::nexcf_ci(TS_gauss[[i]],TS_gauss[[j]],conflevel=0.1)
        C[i,j]<-temp$rxy
        P[i,j]<-P[j,i]<-temp$pval
        C[j,i]=C[i,j]
        rm(temp)
        C_gauss[i,j]<-temp_gauss$rxy
        P_gauss[i,j]<-P_gauss[j,i]<-temp_gauss$pval
        C_gauss[j,i]=C_gauss[i,j]
        rm(temp_gauss)
      }
    }

    ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_",run)]][[paste0("CLUSTER", cluster)]]$C <- C
    ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_",run)]][[paste0("CLUSTER", cluster)]]$P <- P
    ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_",run)]]$mean[counter] <- mean(C, na.rm = T)
    ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_",run)]]$range[counter] <- range(C, na.rm = T)
    ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_",run)]][[paste0("CLUSTER", cluster)]]$C_gauss <- C_gauss
    ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_",run)]][[paste0("CLUSTER", cluster)]]$P_gauss <- P_gauss
    ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_",run)]]$mean_100gauss[counter] <- mean(C_gauss, na.rm = T)
    ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_",run)]]$range_100gauss[counter] <- range(C_gauss, na.rm = T)
    counter = counter + 1

  }
  
  ## 2 4 global correlation
  
  ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]] <- list()
  ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$mean <- list()
  ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$range <- list()
  ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$mean_100gauss <- list()
  ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$range_100gauss <- list()

  entity_list <- ANALYSIS$NETWORK$entity_meta$entity_id


  TS <- list()
  TS_gauss <- list()

  for(entity in entity_list){
    site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
    TS[[paste0("Entity", entity)]] <- zoo(x = DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]][[paste0("ISOT_", run)]],
                                          order.by = seq(from = 1950-DATA_past1000$time[1], to = 1950-DATA_past1000$time[2], by = -1))
    TS_gauss[[paste0("Entity", entity)]] <- gaussdetr(zoo(x = DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]][[paste0("ISOT_", run)]],
                                                          order.by = seq(from = 1950-DATA_past1000$time[1], to = 1950-DATA_past1000$time[2], by = -1)), tsc.in = 100)$Xsmooth
  }

  C<-matrix(NA,nrow=length(TS),ncol=length(TS))
  colnames(C)<-rownames(C)<-names(TS)
  C_gauss <- P_gauss <- P <- C

  for (i in 1:(length(TS)-1)){
    print(i)
    for (j in (i+1):length(TS)){
      temp<-nest::nexcf_ci(TS[[i]],TS[[j]],conflevel=0.1)
      temp_gauss<-nest::nexcf_ci(TS_gauss[[i]],TS_gauss[[j]],conflevel=0.1)
      C[i,j]<-temp$rxy
      P[i,j]<-P[j,i]<-temp$pval
      C[j,i]=C[i,j]
      rm(temp)
      C_gauss[i,j]<-temp_gauss$rxy
      P_gauss[i,j]<-P_gauss[j,i]<-temp_gauss$pval
      C_gauss[j,i]=C_gauss[i,j]
      rm(temp_gauss)
    }
  }

  ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$C <- C
  ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$P <- P
  ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$mean <- mean(C, na.rm = T)
  ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$range <- range(C, na.rm = T)
  ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$C_gauss <- C_gauss
  ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$P_gauss <- P_gauss
  ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$mean_100gauss <- mean(C_gauss, na.rm = T)
  ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$range_100gauss <- range(C_gauss, na.rm = T)

  remove(cluster, counter, e.lat, e.long, entity, gridbox, i, ii, j, site, TS, site_list, s, Point_Lyr, P, hc, gridbox_list, entity_list, double_time, cluster_list, C, C_sig, a)

}

remove(C_gauss, P_gauss, TS_gauss)

#################################################
## CHRONOLOGY SENSITIVITY #######################
#################################################

ANALYSIS$NETWORK$GLOBAL_CHRONO <- list()

entity_list <- ANALYSIS$NETWORK$entity_meta$entity_id

for(chronology in c("interp_age", "lin_interp_age", "lin_reg_age", "Bchron_age", "Bacon_age", "OxCal_age", "copRa_age", "StalAge_age")){
  print(chronology)
  TS <- list()
  
  for(entity in entity_list){
    if(all(is.na(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$chron[[chronology]]))){
      TS[[paste0("Entity", entity)]] <- NA
      next
      }
    if(entity == 226){next}
    s <- list()
    s$time <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$chron[[chronology]]
    s$value <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$chron$d18O_measurement
    # zoo cannot handle objects where order.by has two elements which is why they are sorted out here (no better option found)
    double_time <- as.tibble(s) %>% group_by(time) %>% count() %>% filter(n>1)
    s <- as.tibble(s) %>% filter(!time %in% double_time$time) %>% filter(!is.na(value)) %>% filter(!is.na(time))
    
    TS[[paste0("Entity", entity)]] <- zoo(x = s$value, order.by = s$time)
    rm(s, double_time)
  }
  
  C<-matrix(NA,nrow=length(TS),ncol=length(TS))
  colnames(C)<-rownames(C)<-names(TS)
  P <- C
  
  for (i in 1:(length(TS)-1)){
    print(i)
    for (j in (i+1):length(TS)){
      if(is.na(TS[[i]]) | is.na(TS[[j]])){
        C[i,j] <- NA
        P[i,j] <- NA
      }
      temp<-nest::nexcf_ci(TS[[i]],TS[[j]],conflevel=0.1)
      C[i,j]<-temp$rxy
      P[i,j]<-P[j,i]<-temp$pval
      C[j,i]=C[i,j]
      rm(temp)
    }
  }
  
  ANALYSIS$NETWORK$GLOBAL_CHRONO[[paste0(chronology,"_C")]] <- C
  ANALYSIS$NETWORK$GLOBAL_CHRONO[[paste0(chronology,"_P")]] <- P
  
}

# Compare and choose highest

C<-matrix(NA,nrow=length(TS),ncol=length(TS))
colnames(C)<-rownames(C)<-names(TS)

mat_comp <- ANALYSIS$NETWORK$GLOBAL_CHRONO
mat_comp$interp_age_C[mat_comp$interp_age_P>0.1] = NA
mat_comp$lin_interp_age_C[mat_comp$lin_interp_age_P>0.1] = NA
mat_comp$lin_reg_age_C[mat_comp$lin_reg_age_P>0.1] = NA
mat_comp$Bchron_age_C[mat_comp$Bchron_age_C>0.1] = NA
mat_comp$Bacon_age_C[mat_comp$Bacon_age_C>0.1] = NA
mat_comp$OxCal_age_C[mat_comp$OxCal_age_C>0.1] = NA
mat_comp$copRa_age_C[mat_comp$copRa_age_P>0.1] = NA
mat_comp$StalAge_age_C[mat_comp$StalAge_age_P>0.1] = NA

for (i in 1:(length(TS)-1)){
  for (j in (i+1):length(TS)){
    if(all(is.na(c(mat_comp$interp_age_C[i,j], mat_comp$lin_interp_age_C[i,j], mat_comp$lin_reg_age_C[i,j], mat_comp$Bchron_age_C[i,j],
                   mat_comp$Bacon_age_C[i,j], mat_comp$OxCal_age_C[i,j], mat_comp$copRa_age_C[i,j], mat_comp$StalAge_age_C[i,j])))){
      C[i,j] = NA
    } else{
      temp <- c(mat_comp$interp_age_C[i,j], mat_comp$lin_interp_age_C[i,j], mat_comp$lin_reg_age_C[i,j], mat_comp$Bchron_age_C[i,j],
                mat_comp$Bacon_age_C[i,j], mat_comp$OxCal_age_C[i,j], mat_comp$copRa_age_C[i,j], mat_comp$StalAge_age_C[i,j])
      index = which.max(abs(temp))
      C[i,j] <- temp[index]
    }
    C[j,i] <- C[i,j]
  }
}

ANALYSIS$NETWORK$GLOBAL_CHRONO$C_max <- C


#################################################
## SUMMARY ######################################
#################################################

table <- array(dim = c(11,14))

colnames(table) <- c("group", "total", 
                     # 3              4               5             6                     7                   8
                     "raw_mean", "raw_95%", "raw_5%", "raw_mean_gauss", "raw_95%_gauss", "raw_5%_gauss", 
                     "sim_mean", "sim_95%", "sim_5%", "sim_mean_gauss", "sim_95%_gauss", "sim_5%_gauss")
rownames(table) <- c("site", "gridbox", 
                     "cluster1-India", "cluster2-SA", "cluster3-Europe", "cluster4-Afica", 
                     "cluster5_Asia", "cluster6-NA", "cluster7-Arabia" , "cluster8-NZ",
                     "global")

table[1,2] <- dim(ANALYSIS$NETWORK$entity_meta %>% select(site_id) %>% group_by(site_id) %>% count())[1]
table[1,1] <- dim(ANALYSIS$NETWORK$entity_meta %>% select(site_id) %>% group_by(site_id) %>% count() %>% filter(n>1))[1]
table[2,2] <- dim(ANALYSIS$NETWORK$entity_meta %>% select(gridbox_id) %>% group_by(gridbox_id) %>% count())[1]
table[2,1] <- dim(ANALYSIS$NETWORK$entity_meta %>% select(gridbox_id) %>% group_by(gridbox_id) %>% count() %>% filter(n>1))[1]
table[3,2] <- (ANALYSIS$NETWORK$entity_meta %>% filter(cluster_id == 1) %>% count())$n
table[4,2] <- (ANALYSIS$NETWORK$entity_meta %>% filter(cluster_id == 2) %>% count())$n
table[5,2] <- (ANALYSIS$NETWORK$entity_meta %>% filter(cluster_id == 3) %>% count())$n
table[6,2] <- (ANALYSIS$NETWORK$entity_meta %>% filter(cluster_id == 4) %>% count())$n
table[7,2] <- (ANALYSIS$NETWORK$entity_meta %>% filter(cluster_id == 5) %>% count())$n
table[8,2] <- (ANALYSIS$NETWORK$entity_meta %>% filter(cluster_id == 6) %>% count())$n
table[9,2] <- (ANALYSIS$NETWORK$entity_meta %>% filter(cluster_id == 7) %>% count())$n
table[10,2] <- (ANALYSIS$NETWORK$entity_meta %>% filter(cluster_id == 8) %>% count())$n
table[11,2] <- length(ANALYSIS$NETWORK$entity_meta$entity_id)

# site ##########################################
site_corr <- list()
site_list <- DATA_past1000$CAVES$entity_info %>% select(site_id, entity_id) %>%
  filter(entity_id %in% ANALYSIS$NETWORK$entity_meta$entity_id) %>% group_by(site_id) %>% count() %>% filter(n>1)
for(site in site_list$site_id){
  site_corr <- c(site_corr, ANALYSIS$NETWORK$SITES[[paste0("SITE",site)]]$C[ANALYSIS$NETWORK$SITES[[paste0("SITE",site)]]$P<0.1])
}
site_corr = as.numeric(site_corr)
table[1,3] <- mean(site_corr, na.rm = T)
table[1,4] <- quantile(site_corr, na.rm = T, probs = seq(0,1,0.05))[20]
table[1,5] <- quantile(site_corr, na.rm = T, probs = seq(0,1,0.05))[2]
site_corr= list()
for(site in site_list$site_id){
  site_corr <- c(site_corr, ANALYSIS$NETWORK$SITES[[paste0("SITE",site)]]$C_gauss[ANALYSIS$NETWORK$SITES[[paste0("SITE",site)]]$P_gauss<0.1])
}
site_corr = as.numeric(site_corr)
table[1,6] <- mean(site_corr, na.rm = T)
table[1,7] <- quantile(site_corr, na.rm = T, probs = seq(0,1,0.05))[20]
table[1,8] <- quantile(site_corr, na.rm = T, probs = seq(0,1,0.05))[2]

# gridbox #######################################
corr <- list()
list <- ANALYSIS$NETWORK$entity_meta %>% select(gridbox_id, entity_id) %>% group_by(gridbox_id) %>% count() %>% filter(n>1)
for(gridbox in list$gridbox_id){
  corr <- c(corr, ANALYSIS$NETWORK$GRIDBOX[[paste0("GRIDBOX",gridbox)]]$C[ANALYSIS$NETWORK$GRIDBOX[[paste0("GRIDBOX",gridbox)]]$P<0.1])
}
corr = as.numeric(corr)

table[2,3] <- mean(corr, na.rm = T)
table[2,4] <- quantile(corr, na.rm = T, probs = seq(0,1,0.05))[20]
table[2,5] <- quantile(corr, na.rm = T, probs = seq(0,1,0.05))[2]

corr = list()
for(gridbox in list$gridbox_id){
  corr <- c(corr, ANALYSIS$NETWORK$GRIDBOX[[paste0("GRIDBOX",gridbox)]]$C_gauss[ANALYSIS$NETWORK$GRIDBOX[[paste0("GRIDBOX",gridbox)]]$P_gauss<0.1])
}
corr = as.numeric(corr)

table[2,6] <- mean(corr, na.rm = T)
table[2,7] <- quantile(corr, na.rm = T, probs = seq(0,1,0.05))[20]
table[2,8] <- quantile(corr, na.rm = T, probs = seq(0,1,0.05))[2]

#cluster
for(cluster in 1:8){
  corr = as.numeric(ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER",cluster)]]$C[ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER",cluster)]]$P<0.1])
  
  table[2+cluster,3] <- mean(corr, na.rm = T)
  table[2+cluster,4] <- quantile(corr, na.rm = T, probs = seq(0,1,0.05))[20]
  table[2+cluster,5] <- quantile(corr, na.rm = T, probs = seq(0,1,0.05))[2]
  
  corr = as.numeric(ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER",cluster)]]$C_gauss[ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER",cluster)]]$P_gauss<0.1])
  table[2+cluster,6] <- mean(corr, na.rm = T)
  table[2+cluster,7] <- quantile(corr, na.rm = T, probs = seq(0,1,0.05))[20]
  table[2+cluster,8] <- quantile(corr, na.rm = T, probs = seq(0,1,0.05))[2]
  
  corr = as.numeric(c(ANALYSIS$NETWORK$CLUSTER_SIM_a[[paste0("CLUSTER",cluster)]]$C[ANALYSIS$NETWORK$CLUSTER_SIM_a[[paste0("CLUSTER",cluster)]]$P<0.1],
                      ANALYSIS$NETWORK$CLUSTER_SIM_b[[paste0("CLUSTER",cluster)]]$C[ANALYSIS$NETWORK$CLUSTER_SIM_b[[paste0("CLUSTER",cluster)]]$P<0.1],
                      ANALYSIS$NETWORK$CLUSTER_SIM_c[[paste0("CLUSTER",cluster)]]$C[ANALYSIS$NETWORK$CLUSTER_SIM_c[[paste0("CLUSTER",cluster)]]$P<0.1]))
  
  table[2+cluster,9] <- mean(corr, na.rm = T)
  table[2+cluster,10] <- quantile(corr, na.rm = T, probs = seq(0,1,0.05))[20]
  table[2+cluster,11] <- quantile(corr, na.rm = T, probs = seq(0,1,0.05))[2]
  
  corr = as.numeric(c(ANALYSIS$NETWORK$CLUSTER_SIM_a[[paste0("CLUSTER",cluster)]]$C_gauss[ANALYSIS$NETWORK$CLUSTER_SIM_a[[paste0("CLUSTER",cluster)]]$P_gauss<0.1],
                      ANALYSIS$NETWORK$CLUSTER_SIM_b[[paste0("CLUSTER",cluster)]]$C_gauss[ANALYSIS$NETWORK$CLUSTER_SIM_b[[paste0("CLUSTER",cluster)]]$P_gauss<0.1],
                      ANALYSIS$NETWORK$CLUSTER_SIM_c[[paste0("CLUSTER",cluster)]]$C_gauss[ANALYSIS$NETWORK$CLUSTER_SIM_c[[paste0("CLUSTER",cluster)]]$P_gauss<0.1]))
  table[2+cluster,12] <- mean(corr, na.rm = T)
  table[2+cluster,13] <- quantile(corr, na.rm = T, probs = seq(0,1,0.05))[20]
  table[2+cluster,14] <- quantile(corr, na.rm = T, probs = seq(0,1,0.05))[2]
}

# Africa Sonderfall weil nur 2 drin

corr = as.numeric(ANALYSIS$NETWORK$GLOBAL$C[ANALYSIS$NETWORK$GLOBAL$P<0.1])
table[11,3] <- mean(corr, na.rm = T)
table[11,4] <- quantile(corr, na.rm = T, probs = seq(0,1,0.05))[20]
table[11,5] <- quantile(corr, na.rm = T, probs = seq(0,1,0.05))[2]
corr = as.numeric(ANALYSIS$NETWORK$GLOBAL$C_gauss[ANALYSIS$NETWORK$GLOBAL$P_gauss<0.1])
table[11,6] <- mean(corr, na.rm = T)
table[11,7] <- quantile(corr, na.rm = T, probs = seq(0,1,0.05))[20]
table[11,8] <- quantile(corr, na.rm = T, probs = seq(0,1,0.05))[2]

corr = as.numeric(c(ANALYSIS$NETWORK$GLOBAL_SIM_a$C[ANALYSIS$NETWORK$GLOBAL_SIM_a$P<0.1], 
                    ANALYSIS$NETWORK$GLOBAL_SIM_b$C[ANALYSIS$NETWORK$GLOBAL_SIM_b$P<0.1],
                    ANALYSIS$NETWORK$GLOBAL_SIM_c$C[ANALYSIS$NETWORK$GLOBAL_SIM_c$P<0.1]))
table[11,9] <- mean(corr, na.rm = T)
table[11,10] <- quantile(corr, na.rm = T, probs = seq(0,1,0.05))[20]
table[11,11] <- quantile(corr, na.rm = T, probs = seq(0,1,0.05))[2]
corr = as.numeric(c(ANALYSIS$NETWORK$GLOBAL_SIM_a$C_gauss[ANALYSIS$NETWORK$GLOBAL_SIM_a$P_gauss<0.1], 
                    ANALYSIS$NETWORK$GLOBAL_SIM_b$C_gauss[ANALYSIS$NETWORK$GLOBAL_SIM_b$P_gauss<0.1],
                    ANALYSIS$NETWORK$GLOBAL_SIM_c$C_gauss[ANALYSIS$NETWORK$GLOBAL_SIM_c$P_gauss<0.1]))
table[11,12] <- mean(corr, na.rm = T)
table[11,13] <- quantile(corr, na.rm = T, probs = seq(0,1,0.05))[20]
table[11,14] <- quantile(corr, na.rm = T, probs = seq(0,1,0.05))[2]



table_round <- round(table, digits = 3)

cairo_pdf(width=17,height=5,file="Plots/Paper_Plot_6_Network_b_table.pdf")
gridExtra::grid.table(table_round)
dev.off()
