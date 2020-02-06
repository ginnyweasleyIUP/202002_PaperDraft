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
Point_Lyr <- list(
  lat = lats,
  long = longs,
  value =cluster
)

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


## 2 1 site correlation 

#Create list with sites that have multiple entities and which entities

site_list <- DATA_past1000$CAVES$entity_info %>% select(site_id, entity_id) %>%
  filter(entity_id %in% ANALYSIS$NETWORK$entity_meta$entity_id) %>% group_by(site_id) %>% count() %>% filter(n>1)

# Calculation

ANALYSIS$NETWORK$SITES <- list()
ANALYSIS$NETWORK$SITES$mean <- numeric(length(site_list$site_id))
ANALYSIS$NETWORK$SITES$range <- list()
counter = 1

for(site in site_list$site_id){
  
  print(paste("site", site))

  ANALYSIS$NETWORK$SITES[[paste0("SITE", site)]] <- list()
  
  entity_list <- DATA_past1000$CAVES$entity_info %>% select(site_id, entity_id) %>% filter(entity_id %in% ANALYSIS$NETWORK$entity_meta$entity_id) %>%
    filter(site_id == site) %>% select(entity_id)
  
  TS <- list()
  
  for(entity in entity_list$entity_id){
    s <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
    # zoo cannot handle objects where order.by has two elements which is why they are sorted out here (no better option found)
    double_time <- s %>% group_by(interp_age) %>% count() %>% filter(n>1)
    s <- s %>% filter(!interp_age %in% double_time$interp_age) %>% filter(!is.na(d18O_measurement))
    
    TS[[paste0("Entity", entity)]] <- zoo(x = s$d18O_measurement, order.by = s$interp_age)
  }
  
  C<-matrix(NA,nrow=length(TS),ncol=length(TS))
  colnames(C)<-rownames(C)<-names(TS)
  P <- C
  
  for (i in 1:(length(TS)-1)){
    for (j in (i+1):length(TS)){
      temp<-nest::nexcf_ci(TS[[i]],TS[[j]],conflevel=0.1)
      
      C[i,j]<-temp$rxy
      P[i,j]<-P[j,i]<-temp$pval
      C[j,i]=C[i,j]
      rm(temp)
    }
  }
  #corrplot::corrplot(C)
  
  C_sig <- C
  C_sig[P>0.1] <- NA
  
  ANALYSIS$NETWORK$SITES[[paste0("SITE", site)]]$C <- C
  ANALYSIS$NETWORK$SITES[[paste0("SITE", site)]]$P <- P
  ANALYSIS$NETWORK$SITES$mean[counter] <- mean(C, na.rm = T)
  ANALYSIS$NETWORK$SITES$range[counter] <- range(C, na.rm = T)
  counter = counter + 1
  
}

## 2 2 gridbox correlation

gridbox_list <- ANALYSIS$NETWORK$entity_meta %>% select(gridbox_id, entity_id) %>% group_by(gridbox_id) %>% count() %>% filter(n>1)
  
# Calculation

ANALYSIS$NETWORK$GRIDBOX <- list()
ANALYSIS$NETWORK$GRIDBOX$mean <- numeric(length(gridbox_list$gridbox_id))
ANALYSIS$NETWORK$GRIDBOX$range <- list()
counter = 1

for(gridbox in gridbox_list$gridbox_id){
  
  print(paste("gridbox", gridbox))
  
  ANALYSIS$NETWORK$GRIDBOX[[paste0("gridbox", gridbox)]] <- list()
  
  entity_list <- ANALYSIS$NETWORK$entity_meta %>% select(gridbox_id, entity_id) %>% filter(gridbox_id == gridbox)
  
  TS <- list()
  
  for(entity in entity_list$entity_id){
    s <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
    # zoo cannot handle objects where order.by has two elements which is why they are sorted out here (no better option found)
    double_time <- s %>% group_by(interp_age) %>% count() %>% filter(n>1)
    s <- s %>% filter(!interp_age %in% double_time$interp_age) %>% filter(!is.na(d18O_measurement))
    
    TS[[paste0("Entity", entity)]] <- zoo(x = s$d18O_measurement, order.by = s$interp_age)
  }
  
  C<-matrix(NA,nrow=length(TS),ncol=length(TS))
  colnames(C)<-rownames(C)<-names(TS)
  P <- C
  
  for (i in 1:(length(TS)-1)){
    for (j in (i+1):length(TS)){
      temp<-nest::nexcf_ci(TS[[i]],TS[[j]],conflevel=0.1)
      
      C[i,j]<-temp$rxy
      P[i,j]<-P[j,i]<-temp$pval
      C[j,i]=C[i,j]
      rm(temp)
    }
  }
  #corrplot::corrplot(C)
  
  C_sig <- C
  C_sig[P>0.1] <- NA
  
  ANALYSIS$NETWORK$GRIDBOX[[paste0("GRIDBOX", gridbox)]]$C <- C
  ANALYSIS$NETWORK$GRIDBOX[[paste0("GRIDBOX", gridbox)]]$P <- P
  ANALYSIS$NETWORK$GRIDBOX$mean[counter] <- mean(C, na.rm = T)
  ANALYSIS$NETWORK$GRIDBOX$range[counter] <- range(C, na.rm = T)
  counter = counter + 1
  
}

## 2 3 cluster correlation

cluster_list <- ANALYSIS$NETWORK$entity_meta %>% select(cluster_id, entity_id) %>% group_by(cluster_id) %>% count() %>% filter(n>1)

# Calculation

ANALYSIS$NETWORK$CLUSTER <- list()
ANALYSIS$NETWORK$CLUSTER$mean <- numeric(length(cluster_list$cluster_id))
ANALYSIS$NETWORK$CLUSTER$range <- list()
counter = 1

for(cluster in cluster_list$cluster_id){
  
  print(paste("cluster", cluster))
  
  ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]] <- list()
  
  entity_list <- ANALYSIS$NETWORK$entity_meta %>% select(cluster_id, entity_id) %>% filter(cluster_id == cluster)
  
  TS <- list()
  
  for(entity in entity_list$entity_id){
    s <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
    # zoo cannot handle objects where order.by has two elements which is why they are sorted out here (no better option found)
    double_time <- s %>% group_by(interp_age) %>% count() %>% filter(n>1)
    s <- s %>% filter(!interp_age %in% double_time$interp_age) %>% filter(!is.na(d18O_measurement))
    
    TS[[paste0("Entity", entity)]] <- zoo(x = s$d18O_measurement, order.by = s$interp_age)
  }
  
  C<-matrix(NA,nrow=length(TS),ncol=length(TS))
  colnames(C)<-rownames(C)<-names(TS)
  P <- C
  
  for (i in 1:(length(TS)-1)){
    for (j in (i+1):length(TS)){
      temp<-nest::nexcf_ci(TS[[i]],TS[[j]],conflevel=0.1)
      
      C[i,j]<-temp$rxy
      P[i,j]<-P[j,i]<-temp$pval
      C[j,i]=C[i,j]
      rm(temp)
    }
  }
  #corrplot::corrplot(C)
  
  C_sig <- C
  C_sig[P>0.1] <- NA
  
  ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]]$C <- C
  ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]]$P <- P
  ANALYSIS$NETWORK$CLUSTER$mean[counter] <- mean(C, na.rm = T)
  ANALYSIS$NETWORK$CLUSTER$range[counter] <- range(C, na.rm = T)
  counter = counter + 1
  
}

## 2 4 global correlation

ANALYSIS$NETWORK$GLOBAL <- list()
ANALYSIS$NETWORK$GLOBAL$mean <- list()
ANALYSIS$NETWORK$GLOBAL$range <- list()

entity_list <- ANALYSIS$NETWORK$entity_meta$entity_id

  
TS <- list()
  
for(entity in entity_list){
  s <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
  # zoo cannot handle objects where order.by has two elements which is why they are sorted out here (no better option found)
  double_time <- s %>% group_by(interp_age) %>% count() %>% filter(n>1)
  s <- s %>% filter(!interp_age %in% double_time$interp_age) %>% filter(!is.na(d18O_measurement))
  
  TS[[paste0("Entity", entity)]] <- zoo(x = s$d18O_measurement, order.by = s$interp_age)
}
  
C<-matrix(NA,nrow=length(TS),ncol=length(TS))
colnames(C)<-rownames(C)<-names(TS)
P <- C

for (i in 1:(length(TS)-1)){
  print(i)
  for (j in (i+1):length(TS)){
    temp<-nest::nexcf_ci(TS[[i]],TS[[j]],conflevel=0.1)
    
    C[i,j]<-temp$rxy
    P[i,j]<-P[j,i]<-temp$pval
    C[j,i]=C[i,j]
    rm(temp)
  }
}

corrplot::corrplot(C[1:20, 1:20])

C_sig <- C
C_sig[P>0.1] <- NA

ANALYSIS$NETWORK$GLOBAL$C <- C
ANALYSIS$NETWORK$GLOBAL$P <- P
ANALYSIS$NETWORK$GLOBAL$mean <- mean(C, na.rm = T)
ANALYSIS$NETWORK$GLOBAL$range <- range(C, na.rm = T)

remove(cluster, counter, e.lat, e.long, entity, gridbox, i, ii, j, site, TS, site_list, s, Point_Lyr, P, hc, gridbox_list, entity_list, double_time, cluster_list, C, C_sig, a)

#################################################
## SIMULATION ###################################
#################################################

## 2 3 cluster correlation

cluster_list <- ANALYSIS$NETWORK$entity_meta %>% select(cluster_id, entity_id) %>% group_by(cluster_id) %>% count() %>% filter(n>1)

# Calculation

ANALYSIS$NETWORK$CLUSTER_SIM <- list()
ANALYSIS$NETWORK$CLUSTER_SIM$mean <- numeric(length(cluster_list$cluster_id))
ANALYSIS$NETWORK$CLUSTER_SIM$range <- list()
counter = 1

for(cluster in cluster_list$cluster_id){
  
  print(paste("cluster", cluster))
  
  ANALYSIS$NETWORK$CLUSTER_SIM[[paste0("CLUSTER", cluster)]] <- list()
  
  entity_list <- ANALYSIS$NETWORK$entity_meta %>% select(cluster_id, entity_id) %>% filter(cluster_id == cluster)
  
  TS <- list()
  
  for(entity in entity_list$entity_id){
    site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
    TS[[paste0("Entity", entity)]] <- zoo(x = DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ISOT, order.by = seq(from = 1100, to = -49, by = -1))
  }
  
  C<-matrix(NA,nrow=length(TS),ncol=length(TS))
  colnames(C)<-rownames(C)<-names(TS)
  P <- C
  
  for (i in 1:(length(TS)-1)){
    for (j in (i+1):length(TS)){
      temp<-nest::nexcf_ci(TS[[i]],TS[[j]],conflevel=0.1)
      
      C[i,j]<-temp$rxy
      P[i,j]<-P[j,i]<-temp$pval
      C[j,i]=C[i,j]
      rm(temp)
    }
  }
  #corrplot::corrplot(C)
  
  C_sig <- C
  C_sig[P>0.1] <- NA
  
  ANALYSIS$NETWORK$CLUSTER_SIM[[paste0("CLUSTER", cluster)]]$C <- C
  ANALYSIS$NETWORK$CLUSTER_SIM[[paste0("CLUSTER", cluster)]]$P <- P
  ANALYSIS$NETWORK$CLUSTER_SIM$mean[counter] <- mean(C, na.rm = T)
  ANALYSIS$NETWORK$CLUSTER_SIM$range[counter] <- range(C, na.rm = T)
  counter = counter + 1
  
}

## 2 4 global correlation

ANALYSIS$NETWORK$GLOBAL_SIM <- list()
ANALYSIS$NETWORK$GLOBAL_SIM$mean <- list()
ANALYSIS$NETWORK$GLOBAL_SIM$range <- list()

entity_list <- ANALYSIS$NETWORK$entity_meta$entity_id


TS <- list()

for(entity in entity_list){
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  TS[[paste0("Entity", entity)]] <- DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ISOT
}

C<-matrix(NA,nrow=length(TS),ncol=length(TS))
colnames(C)<-rownames(C)<-names(TS)
P <- C

for (i in 1:(length(TS)-1)){
  for (j in (i+1):length(TS)){
    temp<- cor.test(TS[[i]],TS[[j]])
    #temp$estimate[[1]]
    #j = j+1
    C[i,j]<-temp$estimate[[1]]
    P[i,j]<-P[j,i]<-temp$p.value
    C[j,i]=C[i,j]
    rm(temp)
  }
}

corrplot::corrplot(C[1:20,1:20])

C_sig <- C
C_sig[P>0.1] <- NA

ANALYSIS$NETWORK$GLOBAL_SIM$C <- C
ANALYSIS$NETWORK$GLOBAL_SIM$P <- P
ANALYSIS$NETWORK$GLOBAL_SIM$mean <- mean(C, na.rm = T)
ANALYSIS$NETWORK$GLOBAL_SIM$range <- range(C, na.rm = T)

remove(cluster, counter, e.lat, e.long, entity, gridbox, i, ii, j, site, TS, site_list, s, Point_Lyr, P, hc, gridbox_list, entity_list, double_time, cluster_list, C, C_sig, a)


#################################################
## SUMMARY ######################################
#################################################

table <- array(dim = c(11,8))

colnames(table) <- c("group", "total", "raw_mean", "raw_upper", "raw_lower", "sim_mean", "sim_upper", "sim_lower")
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
  table[2+cluster,6] <- mean(ANALYSIS$NETWORK$CLUSTER_SIM[[paste0("CLUSTER", cluster)]]$C[ANALYSIS$NETWORK$CLUSTER_SIM[[paste0("CLUSTER", cluster)]]$P<0.15], na.rm = T)
  table[2+cluster,7] <- range(ANALYSIS$NETWORK$CLUSTER_SIM[[paste0("CLUSTER", cluster)]]$C[ANALYSIS$NETWORK$CLUSTER_SIM[[paste0("CLUSTER", cluster)]]$P<0.15], na.rm = T)[1]
  table[2+cluster,8] <- range(ANALYSIS$NETWORK$CLUSTER_SIM[[paste0("CLUSTER", cluster)]]$C[ANALYSIS$NETWORK$CLUSTER_SIM[[paste0("CLUSTER", cluster)]]$P<0.15], na.rm = T)[2] 
}

# Africa Sonderfall weil nur 2 drin:


table[11,3] <- mean(ANALYSIS$NETWORK$GLOBAL$mean, na.rm = T)
table[11,4] <- max(ANALYSIS$NETWORK$GLOBAL$range, na.rm = T)
table[11,5] <- min(ANALYSIS$NETWORK$GLOBAL$range, na.rm = T)
table[11,6] <- mean(ANALYSIS$NETWORK$GLOBAL_SIM$mean, na.rm = T)
table[11,7] <- max(ANALYSIS$NETWORK$GLOBAL_SIM$range, na.rm = T)
table[11,8] <- min(ANALYSIS$NETWORK$GLOBAL_SIM$range, na.rm = T)

table_round <- round(table, digits = 3)

cairo_pdf(width=9,height=5,file="Plots/Paper_Plot_6_Network_b_table.pdf")
gridExtra::grid.table(table_round)
dev.off()

# #################################################
# ## TRY Principle Component Analysis for site 2 ##
# #################################################
# 
# site = 2
# entity_list = c(14,620,621,623)
# time_start = -49
# time_stop = 1100
# dt <- list()
# 
# for(ii in 1:length(entity_list)){
#   entity = entity_list[ii]
#   if(PaleoSpec::FirstElement(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$interp_age)>time_start){
#     time_start = PaleoSpec::FirstElement(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$interp_age)
#   }
#   if(PaleoSpec::LastElement(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$interp_age)<time_stop){
#     time_stop = PaleoSpec::LastElement(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$interp_age)
#   }
#   
#   dt = c(dt, mean(diff(DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$interp_age)))
#   
# }
# 
# dt = mean(as.numeric(dt))
# 
# 
# # Jetzt equidistant or blockaveraging?
# 
# TS <- array(dim = c(83,4))
# 
# for(ii in 1:4){
#   entity = entity_list[ii]
#   TS[,ii] = as.numeric(PaleoSpec::MakeEquidistant(t.x = DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$interp_age, 
#                                                              t.y = DATA_past1000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_dw_eq,
#                                                              time.target = time_target))
#   
# }
# 
# PC <- prcomp(TS)
# 
# sum(PC$rotation[,1])/sum(PC$rotation)
# # but what does this tell me?
