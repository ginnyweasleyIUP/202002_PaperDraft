#################################################
## NETWORK ANALYSIS ENSELBLES ###################
#################################################

library(plyr)
library(dplyr)
library(tidyverse)
library(zoo)

# #################################################
# # 1) Get list with chronology-1:1000 - d18O for all
# #################################################
# 
# ENSEMBLE <- list()
# 
# ENSEMBLE$CAVES_RAW <- DATA_past1000$CAVES$record_data
# ENSEMBLE$ensembles <- list()
# ENSEMBLE$names <- read_csv("~/Dokumente/01_Promotion/06_Daten/02_SISAL/SISAL_v2/entity.csv") %>% select(entity_id, entity_name) %>%
#   filter(entity_id %in% DATA_past1000$CAVES$entity_info$entity_id)
# 
# rm_list = list()
# 
# print(".. read in ensembles")
# 
# #for(entity in DATA_past1000$CAVES$entity_info$entity_id){
# for(entity in DATA_past1000$CAVES$entity_info$entity_id[mask_spec]){
#   print(entity)
#   name = ENSEMBLE$names$entity_name[ENSEMBLE$names$entity_id == entity]
#   sample <- read.csv("~/Dokumente/01_Promotion/06_Daten/02_SISAL/SISAL_v2/sample.csv") %>% filter(entity_id == entity)
#   org_chronology <- read.csv("~/Dokumente/01_Promotion/06_Daten/02_SISAL/SISAL_v2/original_chronology.csv") %>% filter(sample_id %in% sample$sample_id) %>%
#     filter(interp_age > 1950-DATA_past1000$time[2] & interp_age < 1950-DATA_past1000$time[1])
#   sample <- sample %>% filter(sample_id %in% org_chronology$sample_id)
#   d18O <- read.csv("~/Dokumente/01_Promotion/06_Daten/02_SISAL/SISAL_v2/d18o.csv") %>% filter(sample_id %in% sample$sample_id)
#   for(chronology in c("Bacon", "Bchron", "copRa", "lininterp", "linReg")){
#     if(file.exists(file = paste0("~/Dokumente/01_Promotion/06_Daten/02_SISAL/SISALv2_ensembles/",entity,"-",name,"-", chronology,".RData"))){
#       rm_list  = c(rm_list, paste0(entity,"-",name,"-", chronology))
#       ens_1000 = get(load(paste0("~/Dokumente/01_Promotion/06_Daten/02_SISAL/SISALv2_ensembles/",entity,"-",name,"-", chronology,".RData"))) %>% 
#         filter(sample_id %in% sample$sample_id)
#       ENSEMBLE$ensembles[[paste0("ENTITY", entity,"_",chronology)]] = ens_1000
#     }
#   }
#   ENSEMBLE$ensembles[[paste0("ENTITY", entity, "_interp_age")]] = org_chronology$interp_age
#   ENSEMBLE$ensembles[[paste0("ENTITY", entity, "_d18O")]] = d18O$d18O_measurement
# }
# 
# rm(list = as.character(rm_list))

#################################################
# 1) Get list with chronology-1:1000 - d18O for all
#################################################

ENSEMBLE <- list()

ENSEMBLE$CAVES_RAW <- DATA_past1000$CAVES$record_data
ENSEMBLE$ensembles <- list()
ENSEMBLE$names <- read_csv("/stacywork/ginnyweasley/02_SISAL/SISAL_v2/entity.csv") %>% select(entity_id, entity_name) %>%
  filter(entity_id %in% DATA_past1000$CAVES$entity_info$entity_id)

rm_list = list()

print(".. read in ensembles")

#for(entity in DATA_past1000$CAVES$entity_info$entity_id){
for(entity in DATA_past1000$CAVES$entity_info$entity_id[mask_spec]){
  print(entity)
  name = ENSEMBLE$names$entity_name[ENSEMBLE$names$entity_id == entity]
  sample <- read.csv("/stacywork/ginnyweasley/02_SISAL/SISAL_v2/sample.csv") %>% filter(entity_id == entity)
  org_chronology <- read.csv("/stacywork/ginnyweasley/02_SISAL/SISAL_v2/original_chronology.csv") %>% filter(sample_id %in% sample$sample_id) %>%
    filter(interp_age > 1950-DATA_past1000$time[2] & interp_age < 1950-DATA_past1000$time[1])
  sample <- sample %>% filter(sample_id %in% org_chronology$sample_id)
  d18O <- read.csv("/stacywork/ginnyweasley/02_SISAL/SISAL_v2/d18o.csv") %>% filter(sample_id %in% sample$sample_id)
  ens_all = list()
  for(chronology in c("Bacon", "Bchron", "copRa", "linInterp", "linReg")){
    if(file.exists(file = paste0("/stacywork/ariana/SISALv2_ensembles/",entity,"-",name,"-", chronology,".RData"))){
      if(entity == 330 && chronology == "linReg"){next}
      if(entity == 514 && chronology == "copRa"){next}
      if(entity == 514 && chronology == "linInterp"){next}
      if(entity == 560 && chronology == "copRa"){next}
      if(entity == 560 && chronology == "linInterp"){next}
      if(entity == 560 && chronology == "linReg"){next}
      if(entity == 613 && chronology == "copRa"){next}
      if(entity == 613 && chronology == "linInterp"){next}
      #if(entity == 613 && chronology == "copRa"){next}
      rm_list  = c(rm_list, paste0(entity,"-",name,"-", chronology))
      ens_1000 = get(load(paste0("/stacywork/ariana/SISALv2_ensembles/",entity,"-",name,"-", chronology,".RData"))) %>% 
        filter(sample_id %in% sample$sample_id) %>% as.data.frame()
      if(is_empty(ens_all)){
        ens_all = ens_1000[3:dim(ens_1000)[2]]
        }else{ens_all = cbind(ens_all, ens_1000[3:dim(ens_1000)[2]])
      }
    }
  }
  ENSEMBLE$ensembles[[paste0("ENTITY", entity, "_age")]] = cbind(ens_all,org_chronology$interp_age)
  ENSEMBLE$ensembles[[paste0("ENTITY", entity, "_d18O")]] = d18O$d18O_measurement
}

rm(list = as.character(rm_list))

#NOTES:
# for eID 85 there are no ensembles

#################################################
## CHRONOLOGY SENSITIVITY #######################
#################################################

print("..start chronology sensitivity")

ENSEMBLE$NETWORK <- list()

entity_list <- DATA_past1000$CAVES$entity_info$entity_id[mask_spec]
#entity_list <- head(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])

C<-matrix(NA,nrow=length(entity_list),ncol=length(entity_list))
colnames(C)<-rownames(C)<-paste0("Entity", entity_list)
P <- C_new <- C

noE_max = 2000

ENSEMBLE$NETWORK$CORRMATRIX <- list()

for(noE in 173:15448){#length eID335
  print(paste0("noE: ", noE))
  name = paste0("corr_",noE)
  ENSEMBLE$NETWORK$CORRMATRIX[[name]] = list()
  ENSEMBLE$NETWORK$CORRMATRIX[[name]]$C = C_new
  ENSEMBLE$NETWORK$CORRMATRIX[[name]]$P = C_new
  for(ecount_1 in 1:(length(entity_list)-1)){
    entity_1 = entity_list[ecount_1]
    #print(paste0("entity1 ", entity_1))
    if(dim(ENSEMBLE$ensembles[[paste0("ENTITY", entity_1, "_age")]])[2]<noE){next}
    TS <- list()
    s <- list()
    s$time <- ENSEMBLE$ensembles[[paste0("ENTITY", entity_1, "_age")]][[noE]]
    s$value <- ENSEMBLE$ensembles[[paste0("ENTITY", entity_1, "_d18O")]]
    if(length(s$time)>length(s$value)){s$time <- s$time[1:length(s$value)]}
    if(length(s$time)<length(s$value)){s$value <- s$value[1:length(s$time)]}
    double_time <- as.tibble(s) %>% group_by(time) %>% count() %>% filter(n>1)
    mask_double = rep(TRUE,length(s$time))
    for(time in double_time$time){
      mask_double[which(time ==s$time)[2:length(which(time ==s$time))]] = F
    }
    rm(time)
    s$time = s$time[mask_double]
    s$value = s$value[mask_double]
    s <- as.tibble(s) %>% filter(!is.na(value)) %>% filter(!is.na(time))
    TS[[1]] <- zoo(x = s$value, order.by = s$time)
    
    for(ecount_2 in (ecount_1+1):length(entity_list)){
      entity_2 = entity_list[ecount_2]
      if(dim(ENSEMBLE$ensembles[[paste0("ENTITY", entity_2, "_age")]])[2]<noE){next}
      s <- list()
      s$time <- ENSEMBLE$ensembles[[paste0("ENTITY", entity_2, "_age")]][[noE]]
      s$value <- ENSEMBLE$ensembles[[paste0("ENTITY", entity_2, "_d18O")]]
      if(length(s$time)>length(s$value)){s$time <- s$time[1:length(s$value)]}
      if(length(s$time)<length(s$value)){s$value <- s$value[1:length(s$time)]}
      double_time <- as.tibble(s) %>% group_by(time) %>% count() %>% filter(n>1)
      mask_double = rep(TRUE,length(s$time))
      for(time in double_time$time){
        mask_double[which(time ==s$time)[2:length(which(time ==s$time))]] = F
      }
      rm(time)
      s$time = s$time[mask_double]
      s$value = s$value[mask_double]
      s <- as.tibble(s) %>% filter(!is.na(value)) %>% filter(!is.na(time))
      TS[[2]] <- zoo(x = s$value, order.by = s$time)
      
      #### CORR
      temp<-nest::nexcf_ci(TS[[1]],TS[[2]],conflevel=0.1)
      
      ENSEMBLE$NETWORK$CORRMATRIX[[name]]$C[ecount_1, ecount_2] = temp$rxy
      ENSEMBLE$NETWORK$CORRMATRIX[[name]]$C[ecount_2, ecount_1] = temp$rxy
      ENSEMBLE$NETWORK$CORRMATRIX[[name]]$P[ecount_1, ecount_2] = temp$pval
      ENSEMBLE$NETWORK$CORRMATRIX[[name]]$P[ecount_2, ecount_1] = temp$pval
      
      ##C max
      if(is.na(temp$pval) | is.na(temp$rxy)){next}
      if(is.na(C[ecount_1,ecount_2]) & temp$pval>0.1){next}
      else if(is.na(C[ecount_1,ecount_2])){
        C[ecount_1,ecount_2] = temp$rxy
        next
      }else if(temp$pval<0.1 & abs(C[ecount_1,ecount_2])<abs(temp$rxy)){
        C[ecount_1,ecount_2] = temp$rxy
      }
      C[ecount_2, ecount_1] = C[ecount_1,ecount_2]
    }
  }
}

ENSEMBLE$NETWORK$C <- C

rm(d18O,double_time, ens_1000, org_chronology, P, rm_list, s, sample, temp, TS, chronology, chronology_2, e_count, e_count2, entity_list, ii,jj, name, noE, entity, entity_2)


# counter = 1
# counter_matrix = 1
# for(e_count in 1:(length(entity_list))-1){
# #for(e_count in 1:(4-1)){  
#   entity = entity_list[e_count]
#   print(paste("Entity:",entity))
#   for(chronology in c("Bacon", "Bchron", "copRa", "lininterp", "linReg")){
#     if(is.null(ENSEMBLE$ensembles[[paste0("ENTITY", entity, "_", chronology)]])){next}
#     noE = dim(ENSEMBLE$ensembles[[paste0("ENTITY", entity,"_", chronology)]])[2]-2
#     #noE = 3
#     for(ii in 1:noE){
#       for(noM in 1:2000){
#       #for(noM in 1:3){
#         for(chron2 in c("Bacon", "Bchron", "copRa", "lininterp", "linReg")){
#           name = paste0("eID_",e_count,"_c1",chronology,"_c2",chron2,"_",noM)
#           ENSEMBLE$NETWORK$CORRMATRIX[[name]] = list()
#           ENSEMBLE$NETWORK$CORRMATRIX[[name]]$C = C_new
#           ENSEMBLE$NETWORK$CORRMATRIX[[name]]$P = C_new
#         }
#       }
#       if(ii%%100 == 0){
#         print(paste("ens Nr.",ii))
#       }
#       
#       
#       for(e_count2 in (e_count+1):length(entity_list)){
#       #for(e_count2 in (e_count+1):4){
#         entity_2 = entity_list[e_count2]
#         
#         print(paste("... entity2:", entity_2))
#         for(chronology_2 in c("Bacon", "Bchron", "copRa", "lininterp", "linReg")){
#           if(is.null(ENSEMBLE$ensembles[[paste0("ENTITY", entity_2, "_", chronology_2)]])){
#             #print(paste(entity_2, chronology_2))
#             next
#             }
#           noE = dim(ENSEMBLE$ensembles[[paste0("ENTITY", entity_2,"_", chronology_2)]])[2]-2
#           #noE = 3
#           for(jj in 1:noE){
#             
#             
#             counter = counter + 1
#             
#             ENSEMBLE$NETWORK$CORRMATRIX[[paste0("eID_",e_count,"_c1",chronology,"_c2",chronology_2,"_",jj)]]$C[e_count,e_count2] = temp$rxy
#             ENSEMBLE$NETWORK$CORRMATRIX[[paste0("eID_",e_count,"_c1",chronology,"_c2",chronology_2,"_",jj)]]$C[e_count2,e_count] = temp$rxy
#             ENSEMBLE$NETWORK$CORRMATRIX[[paste0("eID_",e_count,"_c1",chronology,"_c2",chronology_2,"_",jj)]]$P[e_count,e_count2] = temp$pval
#             ENSEMBLE$NETWORK$CORRMATRIX[[paste0("eID_",e_count,"_c1",chronology,"_c2",chronology_2,"_",jj)]]$P[e_count2,e_count] = temp$pval
#             
#           }
#         }
#       }
#     }
#   }
# }
# 
# ENSEMBLE$NETWORK$C_max <- C

#rm(d18O,double_time, ens_1000, org_chronology, P, rm_list, s, sample, temp, TS, chronology, chronology_2, e_count, e_count2, entity_list, ii,jj, name, noE, entity, entity_2)

# #################################################
# ## CHRONOLOGY SENSITIVITY #######################
# #################################################
# 
# print("..start chronology sensitivity")
# 
# ENSEMBLE$NETWORK <- list()
# 
# entity_list <- DATA_past1000$CAVES$entity_info$entity_id[mask_spec]
# #entity_list <- head(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])
# 
# C<-matrix(NA,nrow=length(entity_list),ncol=length(entity_list))
# colnames(C)<-rownames(C)<-paste0("Entity", entity_list)
# P <- C_new <- C
# 
# noE_max = 2000
# 
# # ENSEMBLE$NETWORK$CORRMATRIX <- list("Bacon" = list(),
# #                                     "Bchron" = list(),
# #                                     "copRa" = list(),
# #                                     "lininterp" = list(), 
# #                                     "linReg" = list())
# 
# ENSEMBLE$NETWORK$CORRMATRIX <- list()
# 
# 
# 
# # Geht alle entities durch, dann alle chronologies und dann alle ensembles. 
# 
# # Nur absolut größte signifikante Correlation wird gespeichert. Rest verworfen. 
# 
# counter = 1
# counter_matrix = 1
# for(e_count in 1:(length(entity_list))-1){
#   #for(e_count in 1:(4-1)){  
#   entity = entity_list[e_count]
#   print(paste("Entity:",entity))
#   for(chronology in c("Bacon", "Bchron", "copRa", "lininterp", "linReg")){
#     if(is.null(ENSEMBLE$ensembles[[paste0("ENTITY", entity, "_", chronology)]])){next}
#     noE = dim(ENSEMBLE$ensembles[[paste0("ENTITY", entity,"_", chronology)]])[2]-2
#     #noE = 3
#     for(ii in 1:noE){
#       for(noM in 1:2000){
#         #for(noM in 1:3){
#         for(chron2 in c("Bacon", "Bchron", "copRa", "lininterp", "linReg")){
#           name = paste0("eID_",e_count,"_c1",chronology,"_c2",chron2,"_",noM)
#           ENSEMBLE$NETWORK$CORRMATRIX[[name]] = list()
#           ENSEMBLE$NETWORK$CORRMATRIX[[name]]$C = C_new
#           ENSEMBLE$NETWORK$CORRMATRIX[[name]]$P = C_new
#         }
#       }
#       if(ii%%100 == 0){
#         print(paste("ens Nr.",ii))
#       }
#       TS <- list()
#       s <- list()
#       s$time <- ENSEMBLE$ensembles[[paste0("ENTITY", entity, "_", chronology)]][[ii+2]]
#       s$value <- ENSEMBLE$ensembles[[paste0("ENTITY", entity, "_d18O")]]
#       double_time <- as.tibble(s) %>% group_by(time) %>% count() %>% filter(n>1)
#       s <- as.tibble(s) %>% filter(!time %in% double_time$time) %>% filter(!is.na(value)) %>% filter(!is.na(time))
#       TS[[1]] <- zoo(x = s$value, order.by = s$time)
#       
#       for(e_count2 in (e_count+1):length(entity_list)){
#         #for(e_count2 in (e_count+1):4){
#         entity_2 = entity_list[e_count2]
#         
#         print(paste("... entity2:", entity_2))
#         for(chronology_2 in c("Bacon", "Bchron", "copRa", "lininterp", "linReg")){
#           if(is.null(ENSEMBLE$ensembles[[paste0("ENTITY", entity_2, "_", chronology_2)]])){
#             #print(paste(entity_2, chronology_2))
#             next
#           }
#           noE = dim(ENSEMBLE$ensembles[[paste0("ENTITY", entity_2,"_", chronology_2)]])[2]-2
#           #noE = 3
#           for(jj in 1:noE){
#             s <- list()
#             s$time <- ENSEMBLE$ensembles[[paste0("ENTITY", entity_2, "_", chronology_2)]][[jj+2]]
#             s$value <- ENSEMBLE$ensembles[[paste0("ENTITY", entity_2, "_d18O")]]
#             if(length(s$time)>length(s$value)){s$time <- s$time[1:length(s$value)]}
#             if(length(s$time)<length(s$value)){s$value <- s$value[1:length(s$time)]}
#             double_time <- as.tibble(s) %>% group_by(time) %>% count() %>% filter(n>1)
#             s <- as.tibble(s) %>% filter(!time %in% double_time$time) %>% filter(!is.na(value)) %>% filter(!is.na(time))
#             TS[[2]] <- zoo(x = s$value, order.by = s$time)
#             
#             counter = counter + 1
#             temp<-nest::nexcf_ci(TS[[1]],TS[[2]],conflevel=0.1)
#             ENSEMBLE$NETWORK$CORRMATRIX[[paste0("eID_",e_count,"_c1",chronology,"_c2",chronology_2,"_",jj)]]$C[e_count,e_count2] = temp$rxy
#             ENSEMBLE$NETWORK$CORRMATRIX[[paste0("eID_",e_count,"_c1",chronology,"_c2",chronology_2,"_",jj)]]$C[e_count2,e_count] = temp$rxy
#             ENSEMBLE$NETWORK$CORRMATRIX[[paste0("eID_",e_count,"_c1",chronology,"_c2",chronology_2,"_",jj)]]$P[e_count,e_count2] = temp$pval
#             ENSEMBLE$NETWORK$CORRMATRIX[[paste0("eID_",e_count,"_c1",chronology,"_c2",chronology_2,"_",jj)]]$P[e_count2,e_count] = temp$pval
#             if(is.na(temp$pval) | is.na(temp$rxy)){next}
#             if(is.na(C[e_count,e_count2]) & temp$pval>0.1){next}
#             else if(is.na(C[e_count,e_count2])){
#               C[e_count,e_count2] = temp$rxy
#               next
#             }else if(temp$pval<0.1 & abs(C[e_count,e_count2])<abs(temp$rxy)){
#               C[e_count,e_count2] = temp$rxy
#             }
#             C[e_count2, e_count] = C[e_count,e_count2]
#           }
#         }
#       }
#     }
#   }
# }
# 
# ENSEMBLE$NETWORK$C_max <- C
# 
# rm(d18O,double_time, ens_1000, org_chronology, P, rm_list, s, sample, temp, TS, chronology, chronology_2, e_count, e_count2, entity_list, ii,jj, name, noE, entity, entity_2)