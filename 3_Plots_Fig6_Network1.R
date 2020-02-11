#################################################
## Paper Figure 6 ###############################
#################################################

## Here analysis and Plotting

## Network Plot 1

#################################################

library(plyr)
library(dplyr)
library(tidyverse)
library(zoo)

## USE nest::network_links

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


source("Functions/Plotting/networkmap_simple3.R")

networkmap_simple2(CMAT = ANALYSIS$NETWORK$GLOBAL$C, 
                   lat = lats, 
                   lon = longs, 
                   title = "Correlation HadCM3 past millenium, sig level = 0.1", 
                   thresh = 0.6)



###################################################################################################

source("Functions/networkmap_simple3.R")

# First Points need to be orderd (lowess doesn't work otherwise)

C_SIM_p <- ANALYSIS$NETWORK$GLOBAL_SIM$C
C_SIM_p[ANALYSIS$NETWORK$GLOBAL_SIM$P>0.1] = NA

C_REC_p <- ANALYSIS$NETWORK$GLOBAL$C
C_REC_p[ANALYSIS$NETWORK$GLOBAL$P>0.1] = NA

plot_dist <- ANALYSIS$NETWORK$DIST
plot_dist[lower.tri(ANALYSIS$NETWORK$DIST)] = NA
lowess_dist <- as.vector(ANALYSIS$NETWORK$DIST[upper.tri(ANALYSIS$NETWORK$DIST)])
o <- order(lowess_dist)
lowess_dist_sorted <- lowess_dist[o]

plot_c_sim <- ANALYSIS$NETWORK$GLOBAL_SIM$C
plot_c_sim[lower.tri(ANALYSIS$NETWORK$GLOBAL_SIM$C, diag = FALSE)] = NA
lowess_c_sim <- as.vector(ANALYSIS$NETWORK$GLOBAL_SIM$C[upper.tri(ANALYSIS$NETWORK$GLOBAL_SIM$C)])
lowess_c_sim_sorted <- lowess_c_sim[o]
plot_c_rec <- ANALYSIS$NETWORK$GLOBAL$C
plot_c_rec[lower.tri(ANALYSIS$NETWORK$GLOBAL$C)] = NA
lowess_c_rec <- as.vector(ANALYSIS$NETWORK$GLOBAL$C[upper.tri(ANALYSIS$NETWORK$GLOBAL$C)])
lowess_c_rec_sorted <- lowess_c_rec[o]

lo <- loess(lowess_c_rec_sorted ~ lowess_dist_sorted)

boxes_sim <- list()
boxes_rec <- list()

for(ii in 1:10){
  boxes_sim[[paste0(ii*2000)]] <- na.omit(as.numeric(plot_c_sim[plot_dist<ii*2000 & plot_dist>(ii-1)*2000]))
  boxes_rec[[paste0(ii*2000)]] <- na.omit(as.numeric(plot_c_rec[plot_dist<ii*2000 & plot_dist>(ii-1)*2000]))
}


# dist_vec_sisal <- as.vector(dist_matrix_sisal[upper.tri(dist_matrix_sisal)])
# o <- order(dist_vec_sisal)
# dist_vec_sisal_sorted <- dist_vec_sisal[o]
# corr_vec_sisal <- as.vector(C_SISAL[upper.tri(C_SISAL)])
# corr_vec_sisal_sorted <- corr_vec_sisal[o]

# x <- 1:10
# y <- c(2,4,6,8,7,12,14,16,18,20)
# lo <- loess(y~x)
# plot(x,y)
# lines(predict(lo), col='red', lwd=2)

# #fitting smoothing splines using smooth.spline(X,Y,df=...)
# fit1<-smooth.spline(age,wage,df=16) #16 degrees of freedom
# #Plotting both cubic and Smoothing Splines 
# plot(age,wage,col="grey",xlab="Age",ylab="Wages")
# points(age.grid,predict(fit,newdata = list(age=age.grid)),col="darkgreen",lwd=2,type="l")
# #adding cutpoints
# abline(v=c(25,40,60),lty=2,col="darkgreen")
# lines(fit1,col="red",lwd=2)
# legend("topright",c("Smoothing Spline with 16 df","Cubic Spline"),col=c("red","darkgreen"),lwd=2)

scaling = 1.5
spacing = 0.7
namcex = 1

pdf(file = "Plots/Paper_Plot_6_Network_a.pdf", height= PLOTTING_VARIABLES$HEIGHT, width = PLOTTING_VARIABLES$WIDTH)
par(mfrow=c(2,2), mai = c(rep(spacing, 4)), mar = c(3,3,2,0.5))
#SIM MAP
networkmap_simple3(CMAT = C_SIM_p, 
                   lat = lats, 
                   lon = longs, 
                   title = "",
                   thresh = 0.15)
mtext("Corr HadCM3 past millenium, p<0.1, c>0.15", side = 3, cex = namcex)
mtext("A", side = 3, adj = 0, cex = namcex)
#SIM Cor-Dist
plot(plot_dist, plot_c_sim, 
     ylim = c(-1,1),
     xlim = c(0,20000),
     ylab = "",
     xlab = "", 
     cex = 1, 
     lwd = 0.5, 
     panel.first = grid(), col = "grey", type = "n")
boxplot(boxes_sim, add = TRUE, at = c(1000,3000,5000,7000,9000,11000,13000,15000,17000,19000), boxwex = 1000, names = "n")
lo <- loess(lowess_c_sim_sorted ~ lowess_dist_sorted, span = 0.2)
lines(lo$x, lo$fitted, lwd = 4, col = "#B2182B")
#lines(lowess(lowess_dist_sorted,lowess_c_sim_sorted, f=0.1), lwd = 4, col = "#B2182B")
mtext("Distance between pairs (km)", side= 1, line = 2)
mtext("B", side = 3, adj = 0, cex = namcex)

#SISAL MAP
networkmap_simple3(CMAT = C_REC_p, 
                   lat = lats, 
                   lon = longs,
                   title = "",
                   thresh = 0.6)
mtext("Corr SISAL past millenium, p<0.1, c>0.6", side = 3, cex = namcex)
mtext("C", side = 3, adj = 0, cex = namcex)

#SISAL Cor-Dist
plot(plot_dist, plot_c_rec, 
     ylim = c(-1,1),
     xlim = c(0,20000),
     ylab = "",
     xlab = "", 
     cex = 1, 
     lwd = 1, 
     panel.first = grid(), col = "grey", type = "n")
boxplot(boxes_rec, add = TRUE, at = c(1000,3000,5000,7000,9000,11000,13000,15000,17000,19000), boxwex = 1000, names = "n")
lo <- loess(lowess_c_rec_sorted ~ lowess_dist_sorted, span = 0.2)
lines(lo$x, lo$fitted, lwd = 4, col = "#B2182B")
mtext("Distance between pairs (km)", side= 1, line = 2)
mtext("D", side = 3, adj = 0, cex = namcex)


dev.off()

#################################################
## cluster matrixes #############################
#################################################

for(cluster in 1:8){
  C_plot <- ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]]$C
  entity_list <- ANALYSIS$NETWORK$entity_meta %>% select(cluster_id, entity_id) %>% filter(cluster_id == cluster)
  colnames(C_plot) <- rownames(C_plot)<- entity_list$entity_id
  pdf(file = paste0("Plots/Network_Plot/corrmatrix_cluster_",cluster, ".pdf"), height= PLOTTING_VARIABLES$HEIGHT, width = PLOTTING_VARIABLES$WIDTH)
  corrplot::corrplot(C_plot, p.mat = ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER", cluster)]]$P)
  dev.off()
}

plot <- STACYmap(coastline = T)+
  theme(panel.border = element_blank(),
        legend.background = element_blank(),
        axis.text = element_blank(),
        legend.text = element_text(size = 8)) 
plot %>% ggsave(filename = paste('Basemap', 'pdf', sep = '.'), plot = ., path = 'Plots/Network_Plot', 
                width = 2*PLOTTING_VARIABLES$WIDTH, height = 2*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")

#################################################
## Curved Network ###############################
#################################################

# tidygraph object erstellen und dann mit geom_conn_bundle bündeln

# Plot_Lyr <- list(
#   lon_start = list(),
#   lat_start = list(),
#   lon_stop = list(),
#   lat_stop = list(),
#   value = list()
# )
# 
# counter = 1
# for (ii in 1:88){
#   entity_1 = DATA_past1000$CAVES$entity_info$entity_id[mask_spec][ii]
#   site_1 = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity_1]
#   for (jj in ii:89){
#     entity_2 = DATA_past1000$CAVES$entity_info$entity_id[mask_spec][jj]
#     site_2 = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity_2]
#     if(site_1 == site_2){next}
#     
#     Plot_Lyr$lon_start[counter] = DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site_1]
#     Plot_Lyr$lat_start[counter] = DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site_1]
#     Plot_Lyr$lon_stop[counter] = DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site_2]
#     Plot_Lyr$lat_stop[counter] = DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site_2]
#     Plot_Lyr$value[counter] = 2*runif(1)-1
#     counter = counter + 1
#   }
# }
# 
# Plot_Lyr$lon_start = as.numeric(Plot_Lyr$lon_start)
# Plot_Lyr$lat_start = as.numeric(Plot_Lyr$lat_start)
# Plot_Lyr$lon_stop = as.numeric(Plot_Lyr$lon_stop)
# Plot_Lyr$lat_stop = as.numeric(Plot_Lyr$lat_stop)
# Plot_Lyr$value = as.numeric(Plot_Lyr$value)
# 
# Plot_Lyr <- as.data.frame(Plot_Lyr)
# 
# # projection: 
# projection = as.character('+proj=robin +datum=WGS84')
# ptlyr_1 <- project(cbind(Plot_Lyr$lon_start, Plot_Lyr$lat_start), proj = as.character(projection)) %>% as_tibble()
# ptlyr_2 <- project(cbind(Plot_Lyr$lon_stop, Plot_Lyr$lat_stop), proj = as.character(projection)) %>% as_tibble()
# 
# Plot_Lyr_p <- data.frame(
#   lon_start = ptlyr_1$V1,
#   lat_start = ptlyr_1$V2,
#   lon_stop = ptlyr_2$V1,
#   lat_stop = ptlyr_2$V2,
#   value = Plot_Lyr$value
# )
# 
# Plot_Lyr_p_2 <- Plot_Lyr_p %>% filter(lon_start != lon_stop & lat_start != lat_stop)
# 
# test <- data.frame(
#   lon_start = Plot_Lyr_p_2$lon_start[1:1000],
#   lat_start = Plot_Lyr_p_2$lat_start[1:1000],
#   lon_stop = Plot_Lyr_p_2$lon_stop[1:1000],
#   lat_stop = Plot_Lyr_p_2$lat_stop[1:1000],
#   value = Plot_Lyr_p_2$value[1:1000]
# )
# 
# plot <- STACYmap(coastline = T) +
#   geom_point(data = test, aes(x = lon_stop, y= lat_stop), color = "green")+
#   geom_point(data = test, aes(x = lon_start, y = lat_start), color = "red")+
#   geom_curve(data = test, aes(x = lon_start, y = lat_start, xend = lon_stop, yend = lat_stop, color = value) ,curvature = 0.33)
# plot
