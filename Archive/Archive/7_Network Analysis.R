#################################################
## Network Analysis #############################
#################################################


# 1) Network between all caves with mask_spec for Simulation
# 2) Network between all caves wirh mask_spec between Records
#       Here We have to make them Equidistant and make Block Averages, so they fit over whole time period...
#       What happens if they have NA's? Should still work
# 3) Distance plot 

library(plyr)
library(dplyr)
library(tidyverse)
library(zoo)

#################################################

NETWORK_ANALYSIS <- list()


#################################################
# SIMULATION Network: ###########################

prep_corr_matrix <- matrix(nrow = 1150, ncol = sum(mask_spec))

colnames(prep_corr_matrix) <- as.character(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])
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

corr_matrix <- Hmisc::rcorr(prep_corr_matrix)
remove(prep_corr_matrix)

NETWORK_ANALYSIS$CORR_MATRIX_SIM <- corr_matrix
NETWORK_ANALYSIS$longs <- longs
NETWORK_ANALYSIS$lats <- lats

source("Functions/networkmap_simple2.R")

index_p_sim<-corr_matrix$P>0.1 # all that are not-significant, as they will later be that that are NA-ed
C_SIM<- corr_matrix$r
C_SIM_p <- C_SIM
C_SIM_p[index_p_sim] <- NA


## Network-Plot
networkmap_simple2(CMAT = C_SIM_p, 
                   lat = lats, 
                   lon = longs, 
                   title = "Correlation HadCM3 past millenium, sig level = 0.1", 
                   thresh = 0.15)


#################################################
## RECORDS Network ##############################

TS <- list()

for(entity in DATA_past1000$CAVES$entity_info$entity_id[mask_spec]){
  s <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
  TS[[paste0("Entity", entity)]] <- zoo(x = s$d18O_measurement, order.by = s$interp_age)
}

C<-matrix(NA,nrow=length(TS),ncol=length(TS))
colnames(C)<-rownames(C)<-names(TS)
P <- C

for (i in 1:(length(TS)-1)){
  for (j in (i+1):length(TS)){
    if(i == 32 || j == 32 || i == 44 || j == 44){
      next
    }
    temp<-nest::nexcf_ci(TS[[i]],TS[[j]],conflevel=0.1)
    
    C[i,j]<-temp$rxy
    P[i,j]<-P[j,i]<-temp$pval
    C[j,i]=C[i,j]
    rm(temp)
  }
}

C_sig <- C
C_sig[P>0.1] = NA

networkmap_simple2(CMAT = C_sig,
                   lat = lats,
                   lon = longs,
                   title = "Correlation HadCM3 past millenium, sig level = 0.1",
                   thresh = 0.4)

NETWORK_ANALYSIS$CORR_MATRIX_REC <- list(
  r = C,
  p = P
)

remove(TS, C,P, C_sig, i,j,s)

#################################################
## DISTANCE MATRIX ##############################
#################################################

dist_matrix <- matrix(NA,nrow=length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec]),ncol=length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec]))

for (ii in 1:(length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])-1)){
  dist_matrix[ii,ii] = 0
  for(jj in (ii+1):length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])){
    dist_matrix[ii,jj] <- fossil::deg.dist(longs[ii], lats[ii], longs[jj], lats[jj])
    dist_matrix[jj,ii] <- dist_matrix[ii,jj]
  }
}

dist_matrix[length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec]),length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])] = 0

NETWORK_ANALYSIS$DIST <- dist_matrix

plotdist <- dist_matrix
plotdist[lower.tri(plotdist)] <- NA
plot_c <- NETWORK_ANALYSIS$CORR_MATRIX_REC$r
plot_c[lower.tri(NETWORK_ANALYSIS$CORR_MATRIX_REC$r)] <- NA
