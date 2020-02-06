library(plyr)
library(dplyr)
library(rgdal)
library(latex2exp)
library(sp)

source("Functions/STACYmap_5.R")
source("Functions/projection_ptlyr.R")
source("Functions/projection_gridlyr.R")

coastline_map <- rgdal::readOGR(dsn ="/home/ginnyweasley/Dokumente/01_Promotion/07_R_Code/naturalearth_10m_physical/ne_10m_coastline.shp", verbose = FALSE)

isot_log_map_plot <- function(Point_lyr,
                              Plot_lyr, 
                              projection = as.character('+proj=robin +datum=WGS84'),
                              pt_size){
  
  coastline_map2 <- coastline_map %>% sp::spTransform(., CRSobj = CRS(projection)) %>% fortify(.)
  
  Point_lyr_p <- projection_ptlyr(Point_lyr, projection)
  Plot_lyr_p <- projection_gridlyr(Plot_lyr3, projection)
  
  
  # map_plot <- map_plot + 
  #   geom_sf(data = gridlyr,
  #           mapping = aes(color = layer, fill = layer), show.legend = T) + 
  #   scale_fill_gradientn(colors = {if (splcol) {colorscheme$grid} else {colorscheme}},
  #                        limits = lims$grid,
  #                        guide = FALSE) +
  #   scale_color_gradientn(colors = {if (splcol) {colorscheme$grid} else {colorscheme}}, 
  #                         limits = lims$grid, guide = FALSE)
  
  plot <- ggplot() +
    geom_sf(data = Plot_lyr_p, aes(color = layer, fill = layer)) +
    scale_colour_gradientn(colors = RColorBrewer::brewer.pal(9, 'BrBG'),
                           limits = c(-allmax, allmax), 
                           breaks = c(-allmax, -log(2), 0, log(2), allmax),
                           labels = c(-(exp(allmax)+1), -1, 0, 1, (exp(allmax)-1)))
  plot
  
    
  geom_point(data = NOT_SITES_p, aes(x = long, y = lat, shape = "4", color = '4'),# shape = 20,
               size = (pt_size), alpha = 0.7, show.legend = c(shape =TRUE)) +
    geom_point(data = USED_SITES_spec_p, aes(x = long, y = lat, shape = "1", color = '1'),# shape = 17,
               size = (pt_size), alpha = 0.7, show.legend = c(shape =TRUE)) +
    geom_point(data = USED_SITES_var_p, aes(x = long, y = lat, shape = "2", color = '2'),# shape = 17,
               size = (pt_size), alpha = 0.7, show.legend = c(shape =TRUE)) +
    geom_point(data = USED_SITES_mean_p, aes(x = long, y = lat, shape = "3", color = '3'),# shape = 17,
               size = (pt_size), alpha = 0.7, show.legend = c(shape =TRUE)) +
    scale_color_manual(name = "SISAL v2 sites", labels = c("600y last mill, d18O res >30", "600y last mill, d18O res >20", "600y last mill, d18O res >10", "other"), 
                        values = c('#DD2525', '#F98D11', '#8D0DC4', '#000000')) +
    scale_shape_manual(name = "SISAL v2 sites", labels = c("600y last mill, d18O res >30", "600y last mill, d18O res >20", "600y last mill, d18O res >10", "other"), 
                       values = c(17, 17, 17, 20)) +
    #geom_polygon(data = coastline_map2, aes(x=long, y=lat, group = group), color = 'black',  size = 0.2, fill = NA, alpha = 0.8) +
    theme(legend.position = c(0.02, 0.02),
          legend.justification = c(0, 0),
          legend.box = 'vertical',
          legend.box.background = element_blank(),
          legend.background = element_rect(colour = 'black'),
          panel.border = element_blank())
  
  return(plot)
}

