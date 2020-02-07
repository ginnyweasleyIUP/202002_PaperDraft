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



#################################################
## SUMMARY ######################################
#################################################

table <- array(dim = c(11,26))

colnames(table) <- c("group", "total", 
                     # 3              4               5             6                     7                   8
                     "raw_mean",   "raw_upper",   "raw_lower",   "raw_mean_gauss",   "raw_upper_gauss",  "raw_lower_gauss", 
                     "sim_mean_a", "sim_upper_a", "sim_lower_a", "sim_mean_a_gauss", "sim_upper_a_gauss", "sim_lower_a_gauss",
                     "sim_mean_b", "sim_upper_b", "sim_lower_b", "sim_mean_b_gauss", "sim_upper_b_gauss", "sim_lower_b_gauss",
                     "sim_mean_c", "sim_upper_c", "sim_lower_c", "sim_mean_c_gauss", "sim_upper_c_gauss", "sim_lower_c_gauss")
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

# site
table[1,3] <- mean(ANALYSIS$NETWORK$SITES$mean)
table[1,4] <- max(as.numeric(ANALYSIS$NETWORK$SITES$range))
table[1,5] <- min(as.numeric(ANALYSIS$NETWORK$SITES$range))

# gridbox
table[2,3] <- mean(ANALYSIS$NETWORK$GRIDBOX$mean)
table[2,4] <- max(as.numeric(ANALYSIS$NETWORK$GRIDBOX$range))
table[2,5] <- min(as.numeric(ANALYSIS$NETWORK$GRIDBOX$range))

#cluster
for(cluster in 1:8){
  table[2+cluster,3] <- mean(ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]]$C[ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]]$P<0.15], na.rm = T)
  table[2+cluster,4] <- range(ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]]$C[ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]]$P<0.15], na.rm = T)[1]
  table[2+cluster,5] <- range(ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]]$C[ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]]$P<0.15], na.rm = T)[2]
  table[2+cluster,6] <- mean(ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]]$C_gauss[ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]]$P_gauss<0.15], na.rm = T)
  table[2+cluster,7] <- range(ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]]$C_gauss[ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]]$P_gauss<0.15], na.rm = T)[1]
  table[2+cluster,8] <- range(ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]]$C_gauss[ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]]$P_gauss<0.15], na.rm = T)[2]
  position <- list("a"<- c(9,10,11), "b"<- c(15,16,17), "c" <- c(21,22,23))
  for(run in c("a", "b", "c")){
    table[2+cluster,position[[run]][1]] <- mean(ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_", run)]][[paste0("CLUSTER", cluster)]]$C[ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_", run)]][[paste0("CLUSTER", cluster)]]$P<0.15], na.rm = T)
    table[2+cluster,position[[run]][2]] <- range(ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_", run)]][[paste0("CLUSTER", cluster)]]$C[ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_", run)]][[paste0("CLUSTER", cluster)]]$P<0.15], na.rm = T)[1]
    table[2+cluster,position[[run]][3]] <- range(ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_", run)]][[paste0("CLUSTER", cluster)]]$C[ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_", run)]][[paste0("CLUSTER", cluster)]]$P<0.15], na.rm = T)[2]  
    table[2+cluster,position[[run]][1]+1] <- mean(ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_", run)]][[paste0("CLUSTER", cluster)]]$C_gauss[ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_", run)]][[paste0("CLUSTER", cluster)]]$P_gauss<0.15], na.rm = T)
    table[2+cluster,position[[run]][2]+1] <- range(ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_", run)]][[paste0("CLUSTER", cluster)]]$C_gauss[ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_", run)]][[paste0("CLUSTER", cluster)]]$P_gauss<0.15], na.rm = T)[1]
    table[2+cluster,position[[run]][3]+1] <- range(ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_", run)]][[paste0("CLUSTER", cluster)]]$C_gauss[ANALYSIS$NETWORK[[paste0("CLUSTER_SIM_", run)]][[paste0("CLUSTER", cluster)]]$P_gauss<0.15], na.rm = T)[2]  
  }

}

# Africa Sonderfall weil nur 2 drin:


table[11,3] <- mean(ANALYSIS$NETWORK$GLOBAL$mean, na.rm = T)
table[11,4] <- max(ANALYSIS$NETWORK$GLOBAL$range, na.rm = T)
table[11,5] <- min(ANALYSIS$NETWORK$GLOBAL$range, na.rm = T)
table[11,6] <- mean(ANALYSIS$NETWORK$GLOBAL$mean_gauss, na.rm = T)
table[11,7] <- max(ANALYSIS$NETWORK$GLOBAL$range_gauss, na.rm = T)
table[11,8] <- min(ANALYSIS$NETWORK$GLOBAL$range_gauss, na.rm = T)

position <- list("a"<- c(9,10,11), "b"<- c(15,16,17), "c" <- c(21,22,23))
for(run in c("a", "b", "c")){
  table[11,position[[run]][1]] <- mean(ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$C[ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$P<0.15], na.rm = T)
  table[11,position[[run]][2]] <- range(ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$C[ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$P<0.15], na.rm = T)[1]
  table[11,position[[run]][3]] <- range(ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$C[ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$P<0.15], na.rm = T)[2]
  table[11,position[[run]][1]+1] <- mean(ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$C_gauss[ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$P_gauss<0.15], na.rm = T)
  table[11,position[[run]][2]+1] <- range(ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$C_gauss[ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$P_gauss<0.15], na.rm = T)[1]
  table[11,position[[run]][3]+1] <- range(ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$C_gauss[ANALYSIS$NETWORK[[paste0("GLOBAL_SIM_",run)]]$P_gauss<0.15], na.rm = T)[2]  
}

table_round <- round(table, digits = 3)

cairo_pdf(width=9,height=5,file="Plots/Paper_Plot_6_Network_b_table.pdf")
gridExtra::grid.table(table_round)
dev.off()
