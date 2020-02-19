#################################################
## Paper Network 2 ##############################
#################################################

order <- list(
  cluster_1 = c(1,6,7,8,2,4,3,5),
  cluster_2 = c(2,3,12,4,8,9,10,1,5,6,7,11,13),
  cluster_3 = c(4,6,7,16,17,20,21,13,14,1,12,2,3,5,8,9,10,11,15,18,19,22),
  cluster_4 = c(1,2),
  cluster_5 = c(5,8,14,17,10,11,19,20,16,21,1,2,3,4,6,7,9,12,13,15,18,22,23),
  cluster_6 = c(7,9,4,5,2,3,1,6,8,10,11,12),
  cluster_7 = c(3,4,5,1,2,6),
  cluster_8 = c(1,3,2)
)

map_window <- list(
  cluster_1 = c(60,  100,  5, 45),
  cluster_2 = c(-80, -40,-25, 20),
  cluster_3 = c(-10,  25, 30, 70),
  cluster_4 = c(15,50, -40, -10),
  cluster_5 = c(100, 140, -15, 40),
  cluster_6 = c(-130, -80, 15, 45),
  cluster_7 = c(30, 60, 10, 45),
  cluster_8 = c(165, 180, -50, -32)
)

box_map <- list(
  cluster_1 = c( 78.75,  82.5,  42.5, 45.0),
  cluster_2 = c(-75.00, -78.75, -7.5, -5),
  cluster_3 = c(11.25, 15, 65, 67.5),
  cluster_4 = c(NA,NA,NA,NA),
  cluster_5 = c((77-0.5)*3.75-180, 77.5*3.75-180, 90-23.5*2.5, 90-24.5*2.5),
  cluster_6 = c(-82.5,  -86.25, 22.5, 20.0),
  cluster_7 = c(52.50,  56.25, 12.5, 15.0),
  cluster_8 = c(168.75, 172.50, -42.5, -40.0)
)

box_corr <- list(
  cluster_1 = c(0.5, 4.5, 4.5,  8.5),
  cluster_2 = c(0.5, 3.5, 10.5, 13.5),
  cluster_3 = c(0.5, 4.5, 18.5, 22.5),
  cluster_4 = c(NA,NA,NA,NA),
  cluster_5 = c(0.5, 4.5, 19.5, 23.5),
  cluster_6 = c(0.5, 2.5, 10.5, 12.5),
  cluster_7 = c(0.5, 3.5, 3.5, 6.5),
  cluster_8 = c(0.5, 2.5, 1.5, 3.5)
)

line_x = 2
line_y = 1.5

names_move <- list(
  cluster_1 = list(x = c(0,0,0,0,0,0,0,0)*line_x, y = c(2,1,0,-1,0,-1,0,0)*line_y),
  cluster_2 = list(x = c(1,1,1,0,0,0,0,0,0,0,0,0,0)*line_x, y = c(1,0,-1,0,0,1,0,0,0,0,0,1,0)*line_y),
  cluster_3 = list(x = c(0,0,1.2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)*line_x, y = c(0.8,0.5,0,0,1,0,-1,1,0,0,0,0,0,0,-0.5,0,0,0,0,0,0,0)*line_y),
  cluster_4 = list(x = c(0,0)*line_x, y = c(0,0)*line_y),
  cluster_5 = list(x = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)*line_x, y = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)*line_y),
  cluster_6 = list(x = c(0,0,0,0,0,0,0,0,0,0,0,0)*line_x, y = c(0.5,0,0.5,0,0.5,0,0,0,0,0.5,0,-0.5)*line_y),
  cluster_7 = list(x = c(0.3,0.3,0.3,0,0,0)*line_x, y = c(1,0,-1,0,0,0)*line_y),
  cluster_8 = list(x = c(-0.3,-0.3,-0.3)*line_x, y = c(0,0,0)*line_y)
)




for(run in c("a", "b","c")){
  for(cluster in 1:8){
    
    plot_data <- list(
      order = order[[paste0("cluster_",cluster)]],
      map_window = map_window[[paste0("cluster_",cluster)]],
      box_map = box_map[[paste0("cluster_",cluster)]],
      box_corr = box_corr[[paste0("cluster_",cluster)]],
      names_move_x = names_move[[paste0("cluster_",cluster)]]$x,
      names_move_y = names_move[[paste0("cluster_",cluster)]]$y
    )

    entity_list <- ANALYSIS$NETWORK$entity_meta %>% select(cluster_id, entity_id, site_id, gridbox_id) %>% filter(cluster_id == cluster)
    entity_list <- entity_list$entity_id[plot_data$order]
    
    
    ##CALC########################
    TS_rec <- list()
    TS_sim <- list()
    
    for(entity in entity_list){
      s <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
      ds <- DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]
      double_time <- s %>% group_by(interp_age) %>% count() %>% filter(n>1)
      double_time <- ds %>% group_by(interp_age) %>% count() %>% filter(n>1)
      s <- s %>% filter(!interp_age %in% double_time$interp_age) %>% filter(!is.na(d18O_measurement))
      ds <- ds %>% filter(!interp_age %in% double_time$interp_age)
      TS_rec[[paste0("Entity", entity)]] <- zoo(x = s$d18O_measurement, order.by = s$interp_age)
      site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
      TS_sim[[paste0("Entity", entity)]] <- zoo(x = ds[[paste0("ISOT_",run)]], order.by = ds$interp_age)
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
      long = ANALYSIS$NETWORK$entity_meta$long[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster][plot_data$order],
      lat = ANALYSIS$NETWORK$entity_meta$lat[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster][plot_data$order]
    )
    
    network_lyr_sim <- C_sim
    network_lyr_sim[P_sim > 0.1] = NA
    network_lyr_rec <- C_rec
    network_lyr_rec[P_rec> 0.1] = NA
    
    network_sim <- as.matrix(C_sim)
    network_corplot <- C_rec
    colnames(network_corplot) <- rownames(network_corplot) <- entity_list
    network_corplot[lower.tri(network_corplot)] <- network_sim[lower.tri(network_sim)]
    
    network_p <- P_rec
    sim_p <- P_sim
    network_p[lower.tri(network_p)] <- sim_p[lower.tri(sim_p)]
    
    
    noPts <- dim(network_lyr_rec)[1]
    col = rev(RColorBrewer::brewer.pal(10, 'RdBu'))
    #col = rev(RColorBrewer::brewer.pal(10, 'RdBu'))
    COLZ_rec <- matrix(col[round((network_lyr_rec+1)*5)+1], ncol = dim(network_lyr_rec)[1])
    COLZ_sim <- matrix(col[round((network_lyr_sim+1)*5)+1], ncol = dim(network_lyr_rec)[1])
    ##PLOT #################################
    for(plot in 1:1){
      cairo_pdf(width=12,height=5,file=paste0("Plots/Appendix/Network_ds_",run,"/Paper_Plot_6_Network_",cluster,".pdf"))
      par(mfrow=c(1,3), new = FALSE, oma = c(1,1,1,1), mar = c(3,3,1,1))
      
      plot(c(plot_data$map_window[1], plot_data$map_window[2]),c(plot_data$map_window[3], plot_data$map_window[4]), 
           type = "n", xlab = "", ylab= "")
      maps::map("world", add = TRUE, col = "grey", interior = FALSE)
      for (i in 1:(noPts - 1)) {
        for (j in (i + 1):(noPts)) {
          if (!is.na(network_lyr_sim[i, j])) {
            lines(c(point_lyr$long[i], point_lyr$long[j]), c(point_lyr$lat[i], point_lyr$lat[j]), col = COLZ_sim[i,j], lwd = 8*abs(network_lyr_sim[i,j]))
          }
        }
      }
      points(x = point_lyr$long, y = point_lyr$lat, pch = 21, cex = 1.5)
      abline(v = seq(from = 0, to = 180, by = 3.75), col = "grey")
      abline(v = seq(from = 0, to = -180, by = -3.75), col = "grey")
      abline(h = seq(from = 90, to = 0, by = -2.5), col = "grey")
      abline(h = seq(from = -90, to = 0, by = 2.5), col = "grey")
      
      text(point_lyr$long + names_move[[paste0("cluster_",cluster)]]$x + 2.5, 
           point_lyr$lat + names_move[[paste0("cluster_",cluster)]]$y, 
           ANALYSIS$NETWORK$entity_meta$entity_id[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster][plot_data$order])
      
      lines(c(plot_data$box_map[1], plot_data$box_map[2]), c(plot_data$box_map[3], plot_data$box_map[3]), col = "black", lwd = "2")
      lines(c(plot_data$box_map[1], plot_data$box_map[2]), c(plot_data$box_map[4], plot_data$box_map[4]), col = "black", lwd = "2")
      lines(c(plot_data$box_map[1], plot_data$box_map[1]), c(plot_data$box_map[3],  plot_data$box_map[4]), col = "black", lwd = "2")
      lines(c(plot_data$box_map[2], plot_data$box_map[2]), c(plot_data$box_map[3],  plot_data$box_map[4]), col = "black", lwd = "2")
      
      mtext(" (a)", side = 3, adj = 0, cex = 1.5, col = "black", line = -2)
      mtext("lat (° N)", side = 2, cex = 1, col = "black", line = 2)
      mtext("long (° E)", side = 1, cex = 1, col = "black", line = 2)
      
      corrplot::corrplot(network_corplot, p.mat = network_p, 
                         tl.cex = 1, cl.cex = 1, 
                         pch.cex = 1, number.cex = 1, tl.col = "black",
                         col = rev(RColorBrewer::brewer.pal(10, 'RdBu')), na.label = " ")
      segments(plot_data$box_corr[1],plot_data$box_corr[4],plot_data$box_corr[2],plot_data$box_corr[4], lwd=3, col="black")
      segments(plot_data$box_corr[1],plot_data$box_corr[3],plot_data$box_corr[2],plot_data$box_corr[3], lwd=3, col="black")
      segments(plot_data$box_corr[1],plot_data$box_corr[3],plot_data$box_corr[1],plot_data$box_corr[4], lwd=3, col="black")
      segments(plot_data$box_corr[2],plot_data$box_corr[3],plot_data$box_corr[2],plot_data$box_corr[4], lwd=3, col="black")
      segments(length(entity_list)+0.5,0.5, 0.5, length(entity_list)+0.5, lwd=5, col="black")
      mtext(" (b)", side = 3, adj = 0, cex = 1.5, col = "black", line = -2)
      
      plot(c(plot_data$map_window[1], plot_data$map_window[2]),c(plot_data$map_window[3], plot_data$map_window[4]),
           type = "n", xlab = "", ylab= "")
      maps::map("world", add = TRUE, col = "grey", interior = FALSE)
      for (i in 1:(noPts - 1)) {
        for (j in (i + 1):(noPts)) {
          if (!is.na(network_lyr_rec[i, j])) {
            lines(c(point_lyr$long[i], point_lyr$long[j]), c(point_lyr$lat[i], point_lyr$lat[j]), col = COLZ_rec[i,j] , lwd = 8*abs(network_lyr_rec[i,j]))
          }
        }
      }
      points(x = point_lyr$long, y = point_lyr$lat, pch = 21, cex = 1.5)
      abline(v = seq(from = 0, to = 180, by = 3.75), col = "grey")
      abline(v = seq(from = 0, to = -180, by = -3.75), col = "grey")
      abline(h = seq(from = 90, to = 0, by = -2.5), col = "grey")
      abline(h = seq(from = -90, to = 0, by = 2.5), col = "grey")
      
      text(point_lyr$long + names_move[[paste0("cluster_",cluster)]]$x + 2.5, 
           point_lyr$lat + names_move[[paste0("cluster_",cluster)]]$y, 
           ANALYSIS$NETWORK$entity_meta$entity_id[ANALYSIS$NETWORK$entity_meta$cluster_id == cluster][plot_data$order])
      lines(c(plot_data$box_map[1], plot_data$box_map[2]), c(plot_data$box_map[3], plot_data$box_map[3]), col = "black", lwd = "2")
      lines(c(plot_data$box_map[1], plot_data$box_map[2]), c(plot_data$box_map[4], plot_data$box_map[4]), col = "black", lwd = "2")
      lines(c(plot_data$box_map[1], plot_data$box_map[1]), c(plot_data$box_map[3],  plot_data$box_map[4]), col = "black", lwd = "2")
      lines(c(plot_data$box_map[2], plot_data$box_map[2]), c(plot_data$box_map[3],  plot_data$box_map[4]), col = "black", lwd = "2")
      mtext(" (c)", side = 3, adj = 0, cex = 1.5, col = "black", line = -2)
      mtext("lat (° N)", side = 2, cex = 1, col = "black", line = 2)
      mtext("long (° E)", side = 1, cex = 1, col = "black", line = 2)
      
      dev.off()
    }
  }
}

rm(C_rec, C_sim, COLZ_rec, COLZ_sim, double_time, ds, names_move, network_corplot, network_lyr_rec, network_lyr_sim, network_p, network_sim, plot_data, 
   point_lyr, s, sim_p, temp_sim, TS_rec, TS_sim, cluster, entity, entity_list, i, j, line_x, line_y, noPts, plot, run, site, var, P_rec, P_sim, order)
