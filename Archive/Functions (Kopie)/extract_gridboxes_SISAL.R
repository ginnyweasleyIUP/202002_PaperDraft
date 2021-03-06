##############################################################
## extract Data for Cave Site from surrounding grid boxes ####
##############################################################

extract_gridboxes_sisal <- function(lon_cave, lat_cave){
  result_list <- list()
  
  edges <- list(
    E1_lon = lon_cave+3.75/2, E1_lat = lat_cave+2.5/2,
    E2_lon = lon_cave-3.75/2, E2_lat = lat_cave+2.5/2,
    E3_lon = lon_cave+3.75/2, E3_lat = lat_cave-2.5/2,
    E4_lon = lon_cave-3.75/2, E4_lat = lat_cave-2.5/2
  )
  
  if(edges$E1_lon<0){edges$E1_lon = 360 + edges$E1_lon}
  if(edges$E2_lon<0){edges$E2_lon = 360 + edges$E2_lon}
  if(edges$E3_lon<0){edges$E3_lon = 360 + edges$E3_lon}
  if(edges$E4_lon<0){edges$E4_lon = 360 + edges$E4_lon}
  
  #we desire that lon_real>long_grid and lat_real<lat_grid
  edges$E1_lon_pos <- which.min(abs(DATA_SISAL_SIM_RAW$lon-edges$E1_lon))
  if(DATA_SISAL_SIM_RAW$lon[edges$E1_lon_pos]>edges$E1_lon){edges$E1_lon_pos <- edges$E1_lon_pos -1}
  edges$E1_lat_pos <- which.min(abs(DATA_SISAL_SIM_RAW$lat-edges$E1_lat))
  if(DATA_SISAL_SIM_RAW$lon[edges$E1_lat_pos]<edges$E1_lat){edges$E1_lat_pos <- edges$E1_lat_pos +1}
  
  edges$E2_lon_pos <- which.min(abs(DATA_SISAL_SIM_RAW$lon-edges$E2_lon))
  if(DATA_SISAL_SIM_RAW$lon[edges$E2_lon_pos]>edges$E2_lon){edges$E2_lon_pos <- edges$E2_lon_pos -1}
  edges$E2_lat_pos <- which.min(abs(DATA_SISAL_SIM_RAW$lat-edges$E2_lat))
  if(DATA_SISAL_SIM_RAW$lon[edges$E2_lat_pos]<edges$E2_lat){edges$E2_lat_pos <- edges$E2_lat_pos +1}
  
  edges$E3_lon_pos <- which.min(abs(DATA_SISAL_SIM_RAW$lon-edges$E3_lon))
  if(DATA_SISAL_SIM_RAW$lon[edges$E3_lon_pos]>edges$E3_lon){edges$E3_lon_pos <- edges$E3_lon_pos -1}
  edges$E3_lat_pos <- which.min(abs(DATA_SISAL_SIM_RAW$lat-edges$E3_lat))
  if(DATA_SISAL_SIM_RAW$lon[edges$E3_lat_pos]<edges$E3_lat){edges$E3_lat_pos <- edges$E3_lat_pos +1}
  
  edges$E4_lon_pos <- which.min(abs(DATA_SISAL_SIM_RAW$lon-edges$E4_lon))
  if(DATA_SISAL_SIM_RAW$lon[edges$E4_lon_pos]>edges$E4_lon){edges$E4_lon_pos <- edges$E4_lon_pos -1}
  edges$E4_lat_pos <- which.min(abs(DATA_SISAL_SIM_RAW$lat-edges$E4_lat))
  if(DATA_SISAL_SIM_RAW$lon[edges$E4_lat_pos]<edges$E4_lat){edges$E4_lat_pos <- edges$E4_lat_pos +1}
  
  ratio <- list()
  
  ratio$E1 <- (edges$E1_lon - DATA_SISAL_SIM_RAW$lon[edges$E1_lon_pos])*(edges$E1_lat - DATA_SISAL_SIM_RAW$lat[edges$E1_lat_pos] + 2.5)/(3.75*2.5)
  ratio$E2 <- (3.75-(edges$E2_lon - DATA_SISAL_SIM_RAW$lon[edges$E2_lon_pos]))*(edges$E2_lat - DATA_SISAL_SIM_RAW$lat[edges$E2_lat_pos] + 2.5)/(3.75*2.5)
  ratio$E3 <- (3.75-(edges$E3_lon - DATA_SISAL_SIM_RAW$lon[edges$E3_lon_pos]))*((DATA_SISAL_SIM_RAW$lat[edges$E3_lat_pos]- edges$E3_lat))/(3.75*2.5)
  ratio$E4 <- (edges$E4_lon - DATA_SISAL_SIM_RAW$lon[edges$E4_lon_pos])*((DATA_SISAL_SIM_RAW$lat[edges$E4_lat_pos] - edges$E4_lat))/(3.75*2.5)
  
  result_list <- c(edges, ratio)
  
  return(result_list)
}