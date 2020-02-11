#################################################
## Paper Network 2 ##############################
#################################################

## CLUSTER 1 ####################################

cluster = 1

entity_list <- ANALYSIS$NETWORK$entity_meta %>% select(cluster_id, entity_id) %>% filter(cluster_id == cluster)
entity_list <- entity_list$entity_id[c(1,7,8,9,3,5,2,4,6)]

TS_rec <- list()
TS_sim <- list()

for(entity in entity_list){
  s <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
  double_time <- s %>% group_by(interp_age) %>% count() %>% filter(n>1)
  s <- s %>% filter(!interp_age %in% double_time$interp_age) %>% filter(!is.na(d18O_measurement))
  TS_rec[[paste0("Entity", entity)]] <- zoo(x = s$d18O_measurement, order.by = s$interp_age)
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  TS_sim[[paste0("Entity", entity)]] <- zoo(x = DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ISOT, order.by = seq(from = 1100, to = -49, by = -1))
}

C_rec<-matrix(NA,nrow=length(TS_rec),ncol=length(TS_rec))
colnames(C_rec)<-rownames(C_rec)<-entity_list
C_sim <- P_sim <- P_rec <- C_rec

for (i in 1:(length(TS_rec)-1)){
  for (j in (i+1):length(TS_rec)){
    temp<-nest::nexcf_ci(TS_rec[[i]],TS_rec[[j]],conflevel=0.1)
    temp_sim <- nest::nexcf_ci(TS_sim[[i]],TS_sim[[j]],conflevel=0.1)
    C_rec[i,j]<-temp$rxy
    P_rec[i,j]<-P_rec[j,i]<-temp$pval
    C_rec[j,i]=C_rec[i,j]
    C_sim[i,j]<-temp_sim$rxy
    P_sim[i,j]<-P_sim[j,i]<-temp_sim$pval
    C_sim[j,i]=C_sim[i,j]
    rm(temp)
  }
}


point_lyr <- data.frame(
  long = ANALYSIS$NETWORK$entity_meta$long[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster],
  lat = ANALYSIS$NETWORK$entity_meta$lat[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster]
)

network_lyr <- C
network_lyr[P > 0.1] = NA

network_sim <- as.matrix(C_sim)

network_corplot <- C_rec
colnames(network_corplot) <- rownames(network_corplot) <- entity_list
network_corplot[lower.tri(network_corplot)] <- network_sim[lower.tri(network_sim)]

network_p <- P_rec
sim_p <- P_sim
network_p[lower.tri(network_p)] <- sim_p[lower.tri(sim_p)]


noPts <- dim(network_lyr)[1]
rbPal <- colorRampPalette(c("blue", "grey", "red"))
#col = rev(RColorBrewer::brewer.pal(10, 'RdBu'))
COLZ <- array((RColorBrewer::brewer.pal(9, 'RdBu'))[as.numeric(cut(c(-1, 1, c(network_lyr)), 
                                                                   breaks = 9))][-c(1:2)], dim = dim(network_lyr))

cairo_pdf(width=12,height=7,file="Plots/Appendix/Paper_Plot_6_Network_1_India.pdf")
par(mfrow=c(1,2), new = FALSE)

plot(c(60, 100), c(5, 45), type = "n", xlab = "", ylab= "", xaxt='n', yaxt='n')
maps::map("world", add = TRUE, col = "grey", interior = FALSE)
for (i in 1:(noPts - 1)) {
  for (j in (i + 1):(noPts)) {
    if (!is.na(network_lyr[i, j])) {
      lines(c(point_lyr$long[i], point_lyr$long[j]), c(point_lyr$lat[i], point_lyr$lat[j]), col = COLZ[i,j], lwd = 6*network_lyr[i,j])
    }
  }
}
symbols(x = point_lyr$long, y = point_lyr$lat, circles = numeric(length(TS_rec))+1, inches = 1/30, 
        bg = "black" , fg = "black", add = TRUE)
abline(v = seq(from = 0, to = 180, by = 3.75), col = "grey")
abline(h = seq(from = 90, to = 0, by = -2.5), col = "grey")

text(point_lyr$long + c(  0, 0, 0, 0,  0, 0, 0,   0,   0) + 2.5, 
     point_lyr$lat + c(-1.2, 0, 0.6, 1, -0.6, 0, 3, 1.6, 0.2), 
     ANALYSIS$NETWORK$entity_meta$entity_id[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster])

lines(c(78.75,  82.50), c(42.5, 42.5),  col = "black", lwd = "2")
lines(c(78.75,  82.50), c(45,45), col = "black", lwd = "2")
lines(c(78.75, 78.75), c(42.5,  45.00), col = "black", lwd = "2")
lines(c(82.5, 82.5), c(42.5,  45.00), col = "black", lwd = "2")

corrplot::corrplot(network_corplot, p.mat = network_p, 
                   #tl.cex = 0.8, cl.cex = 0.8, 
                   pch.cex = 0.8, number.cex = 0.8,
                   col = rev(RColorBrewer::brewer.pal(10, 'RdBu')))
segments(0.5,9.5,4.5,9.5, lwd=3, col="black")
segments(0.5,5.5,4.5,5.5, lwd=3, col="black")
segments(0.5,5.5,0.5,9.5, lwd=3, col="black")
segments(4.5,5.5,4.5,9.5, lwd=3, col="black")

dev.off()

## CLUSTER 2 ####################################

cluster = 2

entity_list <- ANALYSIS$NETWORK$entity_meta %>% select(cluster_id, entity_id) %>% filter(cluster_id == cluster)
entity_list <- entity_list$entity_id[c(2,3,11,7,1,12,5,6,4,8,10,9)]

TS_rec <- list()
TS_sim <- list()

for(entity in entity_list){
  s <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
  double_time <- s %>% group_by(interp_age) %>% count() %>% filter(n>1)
  s <- s %>% filter(!interp_age %in% double_time$interp_age) %>% filter(!is.na(d18O_measurement))
  TS_rec[[paste0("Entity", entity)]] <- zoo(x = s$d18O_measurement, order.by = s$interp_age)
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  TS_sim[[paste0("Entity", entity)]] <- zoo(x = DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ISOT, order.by = seq(from = 1100, to = -49, by = -1))
}

C_rec<-matrix(NA,nrow=length(TS_rec),ncol=length(TS_rec))
colnames(C_rec)<-rownames(C_rec)<-entity_list
C_sim <- P_sim <- P_rec <- C_rec

for (i in 1:(length(TS_rec)-1)){
  for (j in (i+1):length(TS_rec)){
    temp<-nest::nexcf_ci(TS_rec[[i]],TS_rec[[j]],conflevel=0.1)
    temp_sim <- nest::nexcf_ci(TS_sim[[i]],TS_sim[[j]],conflevel=0.1)
    C_rec[i,j]<-temp$rxy
    P_rec[i,j]<-P_rec[j,i]<-temp$pval
    C_rec[j,i]=C_rec[i,j]
    C_sim[i,j]<-temp_sim$rxy
    P_sim[i,j]<-P_sim[j,i]<-temp_sim$pval
    C_sim[j,i]=C_sim[i,j]
    rm(temp)
  }
}


point_lyr <- data.frame(
  long = ANALYSIS$NETWORK$entity_meta$long[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster],
  lat = ANALYSIS$NETWORK$entity_meta$lat[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster]
)

network_lyr <- C
network_lyr[P > 0.1] = NA

network_sim <- as.matrix(C_sim)

network_corplot <- C_rec
colnames(network_corplot) <- rownames(network_corplot) <- entity_list
network_corplot[lower.tri(network_corplot)] <- network_sim[lower.tri(network_sim)]

network_p <- P_rec
sim_p <- P_sim
network_p[lower.tri(network_p)] <- sim_p[lower.tri(sim_p)]


noPts <- dim(network_lyr)[1]
rbPal <- colorRampPalette(c("blue", "grey", "red"))
#col = rev(RColorBrewer::brewer.pal(10, 'RdBu'))
COLZ <- array((RColorBrewer::brewer.pal(9, 'RdBu'))[as.numeric(cut(c(-1, 1, c(network_lyr)), 
                                                                   breaks = 9))][-c(1:2)], dim = dim(network_lyr))

cairo_pdf(width=12,height=7,file="Plots/Appendix/Paper_Plot_6_Network_2_SouthAmerica.pdf")
par(mfrow=c(1,2), new = FALSE)

plot(c(-80, -40), c(-25, 20), type = "n", xlab = "", ylab= "", xaxt='n', yaxt='n')
maps::map("world", add = TRUE, col = "grey", interior = FALSE)
for (i in 1:(noPts - 1)) {
  for (j in (i + 1):(noPts)) {
    if (!is.na(network_lyr[i, j])) {
      lines(c(point_lyr$long[i], point_lyr$long[j]), c(point_lyr$lat[i], point_lyr$lat[j]), col = COLZ[i,j], lwd = 6*network_lyr[i,j])
    }
  }
}
symbols(x = point_lyr$long, y = point_lyr$lat, circles = numeric(length(TS_rec))+1, inches = 1/30, 
        bg = "black" , fg = "black", add = TRUE)
#abline(v = seq(from = 0, to = 180, by = 3.75), col = "grey")
abline(v = seq(from = 0, to = -180, by = -3.75), col = "grey")
abline(h = seq(from = 90, to = 0, by = -2.5), col = "grey")
abline(h = seq(from = -90, to = 0, by = 2.5), col = "grey")

text(point_lyr$long + c(0,0,0,0,0,0,0,0,0,0,0,0) + 2.5, 
     point_lyr$lat  + c(0,2.8,1.4,0,0,0,0,0,0,0,0,0), 
     ANALYSIS$NETWORK$entity_meta$entity_id[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster])

lines(c(-75.00, -78.75), c(-7.5, -7.5),  col = "black", lwd = "2")
lines(c(-75.00, -78.75), c(-5  ,-5), col = "black", lwd = "2")
lines(c(-78.75, -78.75), c(-7.5,-5.0), col = "black", lwd = "2")
lines(c(-75, -75),       c(-7.5,-5.0), col = "black", lwd = "2")

corrplot::corrplot(network_corplot, p.mat = network_p, 
                   #tl.cex = 0.8, cl.cex = 0.8, 
                   pch.cex = 0.8, number.cex = 0.8,
                   col = rev(RColorBrewer::brewer.pal(10, 'RdBu')))
segments(0.5,12.5,3.5,12.5, lwd=3, col="black")
segments(0.5,9.5,3.5,9.5, lwd=3, col="black")
segments(0.5,9.5,0.5,12.5, lwd=3, col="black")
segments(3.5,9.5,3.5,12.5, lwd=3, col="black")

dev.off()


## CLUSTER 3 ####################################

cluster = 3

entity_list <- ANALYSIS$NETWORK$entity_meta %>% select(cluster_id, entity_id) %>% filter(cluster_id == cluster)
entity_list <- entity_list$entity_id[c(7, 16, 6, 4, 11, 3, 2, 13, 14, 18, 19, 22, 10, 1, 12, 9, 8, 15, 20, 21, 17, 5)]

TS_rec <- list()
TS_sim <- list()

for(entity in entity_list){
  s <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
  double_time <- s %>% group_by(interp_age) %>% count() %>% filter(n>1)
  s <- s %>% filter(!interp_age %in% double_time$interp_age) %>% filter(!is.na(d18O_measurement))
  TS_rec[[paste0("Entity", entity)]] <- zoo(x = s$d18O_measurement, order.by = s$interp_age)
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  TS_sim[[paste0("Entity", entity)]] <- zoo(x = DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ISOT, order.by = seq(from = 1100, to = -49, by = -1))
  }

C_rec<-matrix(NA,nrow=length(TS),ncol=length(TS))
colnames(C_rec)<-rownames(C_rec)<-entity_list
C_sim <- P_sim <- P_rec <- C_rec

for (i in 1:(length(TS_rec)-1)){
  for (j in (i+1):length(TS_rec)){
    temp<-nest::nexcf_ci(TS_rec[[i]],TS_rec[[j]],conflevel=0.1)
    temp_sim <- nest::nexcf_ci(TS_sim[[i]],TS_sim[[j]],conflevel=0.1)
    C_rec[i,j]<-temp$rxy
    P_rec[i,j]<-P_rec[j,i]<-temp$pval
    C_rec[j,i]=C_rec[i,j]
    C_sim[i,j]<-temp_sim$rxy
    P_sim[i,j]<-P_sim[j,i]<-temp_sim$pval
    C_sim[j,i]=C_sim[i,j]
    rm(temp)
  }
}


point_lyr <- data.frame(
  long = ANALYSIS$NETWORK$entity_meta$long[ANALYSIS$NETWORK$entity_meta$cluster_id == 3],
  lat = ANALYSIS$NETWORK$entity_meta$lat[ANALYSIS$NETWORK$entity_meta$cluster_id == 3]
)

network_lyr <- C
network_lyr[P > 0.1] = NA

network_sim <- as.matrix(C_sim)

network_corplot <- C_rec
colnames(network_corplot) <- rownames(network_corplot) <- entity_list
network_corplot[lower.tri(network_corplot)] <- network_sim[lower.tri(network_sim)]

network_p <- P_rec
sim_p <- P_sim
network_p[lower.tri(network_p)] <- sim_p[lower.tri(sim_p)]


noPts <- dim(network_lyr)[1]
rbPal <- colorRampPalette(c("blue", "grey", "red"))
#col = rev(RColorBrewer::brewer.pal(10, 'RdBu'))
COLZ <- array((RColorBrewer::brewer.pal(9, 'RdBu'))[as.numeric(cut(c(-1, 1, c(network_lyr)), 
                                      breaks = 9))][-c(1:2)], dim = dim(network_lyr))

cairo_pdf(width=12,height=7,file="Plots/Appendix/Paper_Plot_6_Network_Europe.pdf")
par(mfrow=c(1,2), new = FALSE)

plot(c(-10, 25), c(30, 70), type = "n", xlab = "", ylab= "", xaxt='n', yaxt='n')
maps::map("world", add = TRUE, col = "grey", interior = FALSE)
for (i in 1:(noPts - 1)) {
  for (j in (i + 1):(noPts)) {
    if (!is.na(network_lyr[i, j])) {
      lines(c(point_lyr$long[i], point_lyr$long[j]), c(point_lyr$lat[i], point_lyr$lat[j]), col = COLZ[i,j], lwd = 6*network_lyr[i,j])
    }
  }
}
symbols(x = point_lyr$long, y = point_lyr$lat, circles = numeric(22)+1, inches = 1/30, 
        bg = "black" , fg = "black", add = TRUE)
abline(v = seq(from = 0, to = 180, by = 3.75), col = "grey")
abline(v = seq(from = 0, to = -180, by = -3.75), col = "grey")
abline(h = seq(from = 90, to = 0, by = -2.5), col = "grey")

text(point_lyr$long + c(0,0,0,0,0,-3.5,0.5,0,0,0,0.4,0,0,0,0.3,3,0,0.5,0.5,0,0,0)+1.7, 
     point_lyr$lat + c(0,0,0,+1,0,0,0,-0.3,0,0,-0.5,0,-0.7,+0.7,0.5,0,0,0,0,1.2,-1.2,0), 
     ANALYSIS$NETWORK$entity_meta$entity_id[ANALYSIS$NETWORK$entity_meta$cluster_id == 3])

lines(c(11.25,11.25), c(67.5, 65), col = "black", lwd = "2")
lines(c(15,15), c(65, 67.25), col = "black", lwd = "2")
lines(c(11.25,15), c(65, 65), col = "black", lwd = "2")
lines(c(11.25,15), c(67.5, 67.5), col = "black", lwd = "2")

corrplot::corrplot(network_corplot, p.mat = network_p, 
                   #tl.cex = 0.8, cl.cex = 0.8, 
                   pch.cex = 0.8, number.cex = 0.8,
                   col = rev(RColorBrewer::brewer.pal(10, 'RdBu')))
segments(0.5,22.5,4.5,22.5, lwd=3, col="black")
segments(0.5,18.5,4.5,18.5, lwd=3, col="black")
segments(0.5,18.5,0.5,22.5, lwd=3, col="black")
segments(4.5,18.5,4.5,22.5, lwd=3, col="black")

dev.off()


## CLUSTER 4 ####################################

cluster = 4

entity_list <- ANALYSIS$NETWORK$entity_meta %>% select(cluster_id, entity_id) %>% filter(cluster_id == cluster)
entity_list <- entity_list$entity_id

TS_rec <- list()
TS_sim <- list()

for(entity in entity_list){
  s <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
  double_time <- s %>% group_by(interp_age) %>% count() %>% filter(n>1)
  s <- s %>% filter(!interp_age %in% double_time$interp_age) %>% filter(!is.na(d18O_measurement))
  TS_rec[[paste0("Entity", entity)]] <- zoo(x = s$d18O_measurement, order.by = s$interp_age)
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  TS_sim[[paste0("Entity", entity)]] <- zoo(x = DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ISOT, order.by = seq(from = 1100, to = -49, by = -1))
}

C_rec<-matrix(NA,nrow=length(TS_rec),ncol=length(TS_rec))
colnames(C_rec)<-rownames(C_rec)<-entity_list
C_sim <- P_sim <- P_rec <- C_rec

for (i in 1:(length(TS_rec)-1)){
  for (j in (i+1):length(TS_rec)){
    temp<-nest::nexcf_ci(TS_rec[[i]],TS_rec[[j]],conflevel=0.1)
    temp_sim <- nest::nexcf_ci(TS_sim[[i]],TS_sim[[j]],conflevel=0.1)
    C_rec[i,j]<-temp$rxy
    P_rec[i,j]<-P_rec[j,i]<-temp$pval
    C_rec[j,i]=C_rec[i,j]
    C_sim[i,j]<-temp_sim$rxy
    P_sim[i,j]<-P_sim[j,i]<-temp_sim$pval
    C_sim[j,i]=C_sim[i,j]
    rm(temp)
  }
}


point_lyr <- data.frame(
  long = ANALYSIS$NETWORK$entity_meta$long[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster],
  lat = ANALYSIS$NETWORK$entity_meta$lat[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster]
)

network_lyr <- C
network_lyr[P > 0.1] = NA

network_sim <- as.matrix(C_sim)

network_corplot <- C_rec
colnames(network_corplot) <- rownames(network_corplot) <- entity_list
network_corplot[lower.tri(network_corplot)] <- network_sim[lower.tri(network_sim)]

network_p <- P_rec
sim_p <- P_sim
network_p[lower.tri(network_p)] <- sim_p[lower.tri(sim_p)]


noPts <- dim(network_lyr)[1]
rbPal <- colorRampPalette(c("blue", "grey", "red"))
#col = rev(RColorBrewer::brewer.pal(10, 'RdBu'))
COLZ <- array((RColorBrewer::brewer.pal(9, 'RdBu'))[as.numeric(cut(c(-1, 1, c(network_lyr)), 
                                                                   breaks = 9))][-c(1:2)], dim = dim(network_lyr))

cairo_pdf(width=12,height=7,file="Plots/Appendix/Paper_Plot_6_Network_4_Africa.pdf")
par(mfrow=c(1,2), new = FALSE)

plot(c(15,50), c(-40, -10), type = "n", xlab = "", ylab= "", xaxt='n', yaxt='n')
maps::map("world", add = TRUE, col = "grey", interior = FALSE)
for (i in 1:(noPts - 1)) {
  for (j in (i + 1):(noPts)) {
    if (!is.na(network_lyr[i, j])) {
      lines(c(point_lyr$long[i], point_lyr$long[j]), c(point_lyr$lat[i], point_lyr$lat[j]), col = COLZ[i,j], lwd = 6*network_lyr[i,j])
    }
  }
}
symbols(x = point_lyr$long, y = point_lyr$lat, circles = numeric(length(TS_rec))+1, inches = 1/30, 
        bg = "black" , fg = "black", add = TRUE)
abline(v = seq(from = 0, to = 180, by = 3.75), col = "grey")
abline(v = seq(from = 0, to = -180, by = -3.75), col = "grey")
abline(h = seq(from = 90, to = 0, by = -2.5), col = "grey")
abline(h = seq(from = -90, to = 0, by = 2.5), col = "grey")

text(point_lyr$long + c(0,0) + 2.5, 
     point_lyr$lat  + c(0,0), 
     ANALYSIS$NETWORK$entity_meta$entity_id[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster])

corrplot::corrplot(network_corplot, p.mat = network_p, 
                   #tl.cex = 0.8, cl.cex = 0.8, 
                   pch.cex = 0.8, number.cex = 0.8,
                   col = rev(RColorBrewer::brewer.pal(10, 'RdBu')))

dev.off()

## CLUSTER 5 ####################################

cluster = 5
View(ANALYSIS$NETWORK$entity_meta %>% filter(cluster_id == cluster))

entity_list <- ANALYSIS$NETWORK$entity_meta %>% select(cluster_id, entity_id) %>% filter(cluster_id == cluster)
entity_list <- entity_list$entity_id[c(5,8,13,16,21,17,3,1,10,11,2,9,4,18,19,15,20,14,22,7,12,6)]

TS_rec <- list()
TS_sim <- list()

for(entity in entity_list){
  print(entity)
  s <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
  double_time <- s %>% group_by(interp_age) %>% count() %>% filter(n>1)
  s <- s %>% filter(!interp_age %in% double_time$interp_age) %>% filter(!is.na(d18O_measurement))
  TS_rec[[paste0("Entity", entity)]] <- zoo(x = s$d18O_measurement, order.by = s$interp_age)
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  TS_sim[[paste0("Entity", entity)]] <- zoo(x = DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ISOT, order.by = seq(from = 1100, to = -49, by = -1))
}

C_rec<-matrix(NA,nrow=length(TS_rec),ncol=length(TS_rec))
colnames(C_rec)<-rownames(C_rec)<-entity_list
C_sim <- P_sim <- P_rec <- C_rec

for (i in 1:(length(TS_rec)-1)){
  for (j in (i+1):length(TS_rec)){
    temp<-nest::nexcf_ci(TS_rec[[i]],TS_rec[[j]],conflevel=0.1)
    temp_sim <- nest::nexcf_ci(TS_sim[[i]],TS_sim[[j]],conflevel=0.1)
    C_rec[i,j]<-temp$rxy
    P_rec[i,j]<-P_rec[j,i]<-temp$pval
    C_rec[j,i]=C_rec[i,j]
    C_sim[i,j]<-temp_sim$rxy
    P_sim[i,j]<-P_sim[j,i]<-temp_sim$pval
    C_sim[j,i]=C_sim[i,j]
    rm(temp)
  }
}


point_lyr <- data.frame(
  long = ANALYSIS$NETWORK$entity_meta$long[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster],
  lat = ANALYSIS$NETWORK$entity_meta$lat[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster]
)

network_lyr <- C
network_lyr[P > 0.1] = NA

network_sim <- as.matrix(C_sim)

network_corplot <- C_rec
colnames(network_corplot) <- rownames(network_corplot) <- entity_list
network_corplot[lower.tri(network_corplot)] <- network_sim[lower.tri(network_sim)]

network_p <- P_rec
sim_p <- P_sim
network_p[lower.tri(network_p)] <- sim_p[lower.tri(sim_p)]


noPts <- dim(network_lyr)[1]
rbPal <- colorRampPalette(c("blue", "grey", "red"))
#col = rev(RColorBrewer::brewer.pal(10, 'RdBu'))
COLZ <- array((RColorBrewer::brewer.pal(9, 'RdBu'))[as.numeric(cut(c(-1, 1, c(network_lyr)), 
                                                                   breaks = 9))][-c(1:2)], dim = dim(network_lyr))

cairo_pdf(width=12,height=7,file="Plots/Appendix/Paper_Plot_6_Network_5_Asia.pdf")
par(mfrow=c(1,2), new = FALSE)

plot(c(90, 140), c(-10, 40), type = "n", xlab = "", ylab= "", xaxt='n', yaxt='n')
maps::map("world", add = TRUE, col = "grey", interior = FALSE)
for (i in 1:(noPts - 1)) {
  for (j in (i + 1):(noPts)) {
    if (!is.na(network_lyr[i, j])) {
      lines(c(point_lyr$long[i], point_lyr$long[j]), c(point_lyr$lat[i], point_lyr$lat[j]), col = COLZ[i,j], lwd = 6*network_lyr[i,j])
    }
  }
}
symbols(x = point_lyr$long, y = point_lyr$lat, circles = numeric(length(TS_rec))+1, inches = 1/30, 
        bg = "black" , fg = "black", add = TRUE)
abline(v = seq(from = -1.65, to = 180, by = 3.75), col = "grey")
#abline(v = seq(from = 0, to = -180, by = -3.75), col = "grey")
abline(h = seq(from = 91.25, to = -90, by = -2.5), col = "grey")

text(point_lyr$long + c(  0,-1.5,   0,-2.5,   0,   0,   0, 5.6, 5.6,   0,   0,   0,   5, 2.5, 2.5,-2.5,   0,   0,-3.5,   5,   0,   0) - 2.5, 
     point_lyr$lat  + c(  2, 0.5,   0,+2.5,  +1,   0,   0,-0.5,   0,   4,   2,   0,-1.1,-1.6,-1.5,   0,   0,   0,   0,  -1, 1.6,   0), 
#                     c( 76, 111, 117, 172, 222, 226, 238, 253, 296, 329, 330, 399, 420, 442, 461, 496, 528, 538, 539, 541, 544, 672)
     ANALYSIS$NETWORK$entity_meta$entity_id[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster])

lines(c((77-0.5)*3.75-180, 77.5*3.75-180), c(90-23.5*2.5, 90-23.5*2.5),  col = "black", lwd = "2")
lines(c(76.5*3.75-180, 77.5*3.75-180), c(90-24.5*2.5, 90-24.5*2.5), col = "black", lwd = "2")
lines(c(76.5*3.75-180, 76.5*3.75-180), c(90-23.5*2.5, 90-24.5*2.5), col = "black", lwd = "2")
lines(c(77.5*3.75-180, 77.5*3.75-180), c(90-23.5*2.5, 90-24.5*2.5), col = "black", lwd = "2")

corrplot::corrplot(network_corplot, p.mat = network_p, 
                   #tl.cex = 0.8, cl.cex = 0.8, 
                   pch.cex = 0.8, number.cex = 0.8,
                   col = rev(RColorBrewer::brewer.pal(10, 'RdBu')))
segments(0.5,22.5,4.5,22.5, lwd=3, col="black")
segments(0.5,18.5,4.5,18.5, lwd=3, col="black")
segments(0.5,18.5,0.5,22.5, lwd=3, col="black")
segments(4.5,18.5,4.5,22.5, lwd=3, col="black")

dev.off()

## CLUSTER 6 ####################################

cluster = 6

entity_list <- ANALYSIS$NETWORK$entity_meta %>% select(cluster_id, entity_id) %>% filter(cluster_id == cluster)
entity_list <- entity_list$entity_id[c(2,3,4,5,1,7,9,10,12,11,8,6)]

TS_rec <- list()
TS_sim <- list()

for(entity in entity_list){
  s <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
  double_time <- s %>% group_by(interp_age) %>% count() %>% filter(n>1)
  s <- s %>% filter(!interp_age %in% double_time$interp_age) %>% filter(!is.na(d18O_measurement))
  TS_rec[[paste0("Entity", entity)]] <- zoo(x = s$d18O_measurement, order.by = s$interp_age)
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  TS_sim[[paste0("Entity", entity)]] <- zoo(x = DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ISOT, order.by = seq(from = 1100, to = -49, by = -1))
}

C_rec<-matrix(NA,nrow=length(TS_rec),ncol=length(TS_rec))
colnames(C_rec)<-rownames(C_rec)<-entity_list
C_sim <- P_sim <- P_rec <- C_rec

for (i in 1:(length(TS_rec)-1)){
  for (j in (i+1):length(TS_rec)){
    temp<-nest::nexcf_ci(TS_rec[[i]],TS_rec[[j]],conflevel=0.1)
    temp_sim <- nest::nexcf_ci(TS_sim[[i]],TS_sim[[j]],conflevel=0.1)
    C_rec[i,j]<-temp$rxy
    P_rec[i,j]<-P_rec[j,i]<-temp$pval
    C_rec[j,i]=C_rec[i,j]
    C_sim[i,j]<-temp_sim$rxy
    P_sim[i,j]<-P_sim[j,i]<-temp_sim$pval
    C_sim[j,i]=C_sim[i,j]
    rm(temp)
  }
}


point_lyr <- data.frame(
  long = ANALYSIS$NETWORK$entity_meta$long[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster][c(2,3,4,5,1,7,9,10,12,11,8,6)],
  lat = ANALYSIS$NETWORK$entity_meta$lat[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster][c(2,3,4,5,1,7,9,10,12,11,8,6)]
)

network_lyr <- C
network_lyr[P > 0.1] = NA

network_sim <- as.matrix(C_sim)

network_corplot <- C_rec
colnames(network_corplot) <- rownames(network_corplot) <- entity_list
network_corplot[lower.tri(network_corplot)] <- network_sim[lower.tri(network_sim)]

network_p <- P_rec
sim_p <- P_sim
network_p[lower.tri(network_p)] <- sim_p[lower.tri(sim_p)]


noPts <- dim(network_lyr)[1]
rbPal <- colorRampPalette(c("blue", "grey", "red"))
#col = rev(RColorBrewer::brewer.pal(10, 'RdBu'))
COLZ <- array((RColorBrewer::brewer.pal(9, 'RdBu'))[as.numeric(cut(c(-1, 1, c(network_lyr)), 
                                                                   breaks = 9))][-c(1:2)], dim = dim(network_lyr))

cairo_pdf(width=12,height=7,file="Plots/Appendix/Paper_Plot_6_Network_6_NorthAmerica.pdf")
par(mfrow=c(1,2), new = FALSE)

plot(c(-130, -80), c(15, 45), type = "n", xlab = "", ylab= "", xaxt='n', yaxt='n')
maps::map("world", add = TRUE, col = "grey", interior = FALSE)
for (i in 1:(noPts - 1)) {
  for (j in (i + 1):(noPts)) {
    if (!is.na(network_lyr[i, j])) {
      lines(c(point_lyr$long[i], point_lyr$long[j]), c(point_lyr$lat[i], point_lyr$lat[j]), col = COLZ[i,j], lwd = 6*network_lyr[i,j])
    }
  }
}
symbols(x = point_lyr$long, y = point_lyr$lat, circles = numeric(length(TS_rec))+1, inches = 1/30, 
        bg = "black" , fg = "black", add = TRUE)
#abline(v = seq(from = 0, to = 180, by = 3.75), col = "grey")
abline(v = seq(from = 0, to = -180, by = -3.75), col = "grey")
abline(h = seq(from = 90, to = 0, by = -2.5), col = "grey")
#abline(h = seq(from = -90, to = 0, by = 2.5), col = "grey")

text(point_lyr$long + c(  2,   2,   0,   0,   0,   0,   0,   0,-2.5,   0,   0,   0) + 2.5, 
     point_lyr$lat  + c(0.5,   0,   0,  -1,   0, 1.1,   0,   0,  -1,   0,   0,   0), 
     #                  178, 209, 286, 289, 147, 395, 443, 506, 613, 577, 422, 294
     ANALYSIS$NETWORK$entity_meta$entity_id[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster][c(2,3,4,5,1,7,9,10,12,11,8,6)])

lines(c(-86.25,  -90.00), c(15.0, 15.0),  col = "black", lwd = "2")
lines(c(-86.25,  -90.00), c(17.5, 17.5), col = "black", lwd = "2")
lines(c(-86.25,  -86.25), c(17.5, 15.0), col = "black", lwd = "2")
lines(c(-90.00,  -90.00), c(17.5, 15.0), col = "black", lwd = "2")

corrplot::corrplot(network_corplot, p.mat = network_p, 
                   #tl.cex = 0.8, cl.cex = 0.8, 
                   pch.cex = 0.8, number.cex = 0.8,
                   col = rev(RColorBrewer::brewer.pal(10, 'RdBu')))
segments(0.5,12.5,2.5,12.5, lwd=3, col="black")
segments(0.5,10.5,2.5,10.5, lwd=3, col="black")
segments(0.5,10.5,0.5,12.5, lwd=3, col="black")
segments(2.5,10.5,2.5,12.5, lwd=3, col="black")

dev.off()

## CLUSTER 7 ####################################

cluster = 7

entity_list <- ANALYSIS$NETWORK$entity_meta %>% select(cluster_id, entity_id) %>% filter(cluster_id == cluster)
entity_list <- entity_list$entity_id[c(3,4,5,2,6,1)]

TS_rec <- list()
TS_sim <- list()

for(entity in entity_list){
  s <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
  double_time <- s %>% group_by(interp_age) %>% count() %>% filter(n>1)
  s <- s %>% filter(!interp_age %in% double_time$interp_age) %>% filter(!is.na(d18O_measurement))
  TS_rec[[paste0("Entity", entity)]] <- zoo(x = s$d18O_measurement, order.by = s$interp_age)
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  TS_sim[[paste0("Entity", entity)]] <- zoo(x = DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ISOT, order.by = seq(from = 1100, to = -49, by = -1))
}

C_rec<-matrix(NA,nrow=length(TS_rec),ncol=length(TS_rec))
colnames(C_rec)<-rownames(C_rec)<-entity_list
C_sim <- P_sim <- P_rec <- C_rec

for (i in 1:(length(TS_rec)-1)){
  for (j in (i+1):length(TS_rec)){
    temp<-nest::nexcf_ci(TS_rec[[i]],TS_rec[[j]],conflevel=0.1)
    temp_sim <- nest::nexcf_ci(TS_sim[[i]],TS_sim[[j]],conflevel=0.1)
    C_rec[i,j]<-temp$rxy
    P_rec[i,j]<-P_rec[j,i]<-temp$pval
    C_rec[j,i]=C_rec[i,j]
    C_sim[i,j]<-temp_sim$rxy
    P_sim[i,j]<-P_sim[j,i]<-temp_sim$pval
    C_sim[j,i]=C_sim[i,j]
    rm(temp)
  }
}


point_lyr <- data.frame(
  long = ANALYSIS$NETWORK$entity_meta$long[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster][c(3,4,5,2,6,1)],
  lat = ANALYSIS$NETWORK$entity_meta$lat[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster][c(3,4,5,2,6,1)]
)

network_lyr <- C
network_lyr[P > 0.1] = NA

network_sim <- as.matrix(C_sim)

network_corplot <- C_rec
colnames(network_corplot) <- rownames(network_corplot) <- entity_list
network_corplot[lower.tri(network_corplot)] <- network_sim[lower.tri(network_sim)]

network_p <- P_rec
sim_p <- P_sim
network_p[lower.tri(network_p)] <- sim_p[lower.tri(sim_p)]


noPts <- dim(network_lyr)[1]
rbPal <- colorRampPalette(c("blue", "grey", "red"))
#col = rev(RColorBrewer::brewer.pal(10, 'RdBu'))
COLZ <- array((RColorBrewer::brewer.pal(9, 'RdBu'))[as.numeric(cut(c(-1, 1, c(network_lyr)), 
                                                                   breaks = 9))][-c(1:2)], dim = dim(network_lyr))

cairo_pdf(width=12,height=7,file="Plots/Appendix/Paper_Plot_6_Network_2_Arabia.pdf")
par(mfrow=c(1,2), new = FALSE)

plot(c(30, 60), c(10, 45), type = "n", xlab = "", ylab= "", xaxt='n', yaxt='n')
maps::map("world", add = TRUE, col = "grey", interior = FALSE)
for (i in 1:(noPts - 1)) {
  for (j in (i + 1):(noPts)) {
    if (!is.na(network_lyr[i, j])) {
      lines(c(point_lyr$long[i], point_lyr$long[j]), c(point_lyr$lat[i], point_lyr$lat[j]), col = COLZ[i,j], lwd = 6*network_lyr[i,j])
    }
  }
}
symbols(x = point_lyr$long, y = point_lyr$lat, circles = numeric(length(TS_rec))+1, inches = 1/30, 
        bg = "black" , fg = "black", add = TRUE)
abline(v = seq(from = 0, to = 180, by = 3.75), col = "grey")
abline(v = seq(from = 0, to = -180, by = -3.75), col = "grey")
abline(h = seq(from = 90, to = 0, by = -2.5), col = "grey")
abline(h = seq(from = -90, to = 0, by = 2.5), col = "grey")

text(point_lyr$long + c(1.5, 1.5, 1.5,   0,   0,   0) + 2.5, 
     point_lyr$lat  + c( +1,   0,  -1,   0,   0,   0), 
     #                  546, 547, 548, 351, 573, 305
     ANALYSIS$NETWORK$entity_meta$entity_id[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster][c(3,4,5,2,6,1)])

lines(c(52.50,  56.25), c(12.5, 12.5),  col = "black", lwd = "2")
lines(c(52.50,  56.25), c(15,  15), col = "black", lwd = "2")
lines(c( 52.5,  52.5), c(12.5, 15.0), col = "black", lwd = "2")
lines(c(56.25,  56.25),c(12.5, 15.0), col = "black", lwd = "2")

corrplot::corrplot(network_corplot, p.mat = network_p, 
                   #tl.cex = 0.8, cl.cex = 0.8, 
                   pch.cex = 0.8, number.cex = 0.8,
                   col = rev(RColorBrewer::brewer.pal(10, 'RdBu')))
segments(0.5,6.5,3.5,6.5, lwd=3, col="black")
segments(0.5,3.5,3.5,3.5, lwd=3, col="black")
segments(0.5,3.5,0.5,6.5, lwd=3, col="black")
segments(3.5,3.5,3.5,6.5, lwd=3, col="black")

dev.off()

## CLUSTER 8 ####################################

cluster = 8

entity_list <- ANALYSIS$NETWORK$entity_meta %>% select(cluster_id, entity_id) %>% filter(cluster_id == cluster)
entity_list <- entity_list$entity_id[c(1,4,3,2)]

TS_rec <- list()
TS_sim <- list()

for(entity in entity_list){
  s <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
  double_time <- s %>% group_by(interp_age) %>% count() %>% filter(n>1)
  s <- s %>% filter(!interp_age %in% double_time$interp_age) %>% filter(!is.na(d18O_measurement))
  TS_rec[[paste0("Entity", entity)]] <- zoo(x = s$d18O_measurement, order.by = s$interp_age)
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  TS_sim[[paste0("Entity", entity)]] <- zoo(x = DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ISOT, order.by = seq(from = 1100, to = -49, by = -1))
}

C_rec<-matrix(NA,nrow=length(TS_rec),ncol=length(TS_rec))
colnames(C_rec)<-rownames(C_rec)<-entity_list
C_sim <- P_sim <- P_rec <- C_rec

for (i in 1:(length(TS_rec)-1)){
  for (j in (i+1):length(TS_rec)){
    temp<-nest::nexcf_ci(TS_rec[[i]],TS_rec[[j]],conflevel=0.1)
    temp_sim <- nest::nexcf_ci(TS_sim[[i]],TS_sim[[j]],conflevel=0.1)
    C_rec[i,j]<-temp$rxy
    P_rec[i,j]<-P_rec[j,i]<-temp$pval
    C_rec[j,i]=C_rec[i,j]
    C_sim[i,j]<-temp_sim$rxy
    P_sim[i,j]<-P_sim[j,i]<-temp_sim$pval
    C_sim[j,i]=C_sim[i,j]
    rm(temp)
  }
}


point_lyr <- data.frame(
  long = ANALYSIS$NETWORK$entity_meta$long[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster][c(1,4,3,2)],
  lat = ANALYSIS$NETWORK$entity_meta$lat[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster][c(1,4,3,2)]
)

network_lyr <- C
network_lyr[P > 0.1] = NA

network_sim <- as.matrix(C_sim)

network_corplot <- C_rec
colnames(network_corplot) <- rownames(network_corplot) <- entity_list
network_corplot[lower.tri(network_corplot)] <- network_sim[lower.tri(network_sim)]

network_p <- P_rec
sim_p <- P_sim
network_p[lower.tri(network_p)] <- sim_p[lower.tri(sim_p)]


noPts <- dim(network_lyr)[1]
rbPal <- colorRampPalette(c("blue", "grey", "red"))
#col = rev(RColorBrewer::brewer.pal(10, 'RdBu'))
COLZ <- array((RColorBrewer::brewer.pal(9, 'RdBu'))[as.numeric(cut(c(-1, 1, c(network_lyr)), 
                                                                   breaks = 9))][-c(1:2)], dim = dim(network_lyr))

cairo_pdf(width=12,height=7,file="Plots/Appendix/Paper_Plot_6_Network_8_NZ.pdf")
par(mfrow=c(1,2), new = FALSE)

plot(c(165, 180), c(-50, -32), type = "n", xlab = "", ylab= "", xaxt='n', yaxt='n')
maps::map("world", add = TRUE, col = "grey", interior = FALSE)
for (i in 1:(noPts - 1)) {
  for (j in (i + 1):(noPts)) {
    if (!is.na(network_lyr[i, j])) {
      lines(c(point_lyr$long[i], point_lyr$long[j]), c(point_lyr$lat[i], point_lyr$lat[j]), col = COLZ[i,j], lwd = 6*network_lyr[i,j])
    }
  }
}
symbols(x = point_lyr$long, y = point_lyr$lat, circles = numeric(length(TS_rec))+1, inches = 1/30, 
        bg = "black" , fg = "black", add = TRUE)
abline(v = seq(from = 0, to = 180, by = 3.75), col = "grey")
abline(v = seq(from = 0, to = -180, by = -3.75), col = "grey")
abline(h = seq(from = 90, to = 0, by = -2.5), col = "grey")
abline(h = seq(from = -90, to = 0, by = 2.5), col = "grey")

text(point_lyr$long + c(0,-1,0,0) + 1, 
     point_lyr$lat  + c(-0.5,0.5,0,0), 
     ANALYSIS$NETWORK$entity_meta$entity_id[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster][c(1,4,3,2)])

lines(c(168.75, 172.50), c(-40.0, -40.0),  col = "black", lwd = "2")
lines(c(168.75, 172.50), c(-42.5, -42.5), col = "black", lwd = "2")
lines(c(168.75, 168.75), c(-42.5, -40.0), col = "black", lwd = "2")
lines(c(172.50, 172.50), c(-42.5, -40.0), col = "black", lwd = "2")

corrplot::corrplot(network_corplot, p.mat = network_p, 
                   #tl.cex = 0.8, cl.cex = 0.8, 
                   pch.cex = 0.8, number.cex = 0.8,
                   col = rev(RColorBrewer::brewer.pal(10, 'RdBu')))
segments(0.5,4.5,2.5,4.5, lwd=3, col="black")
segments(0.5,2.5,2.5,2.5, lwd=3, col="black")
segments(0.5,2.5,0.5,4.5, lwd=3, col="black")
segments(2.5,2.5,2.5,4.5, lwd=3, col="black")

dev.off()
