library(plyr)
library(dplyr)
library(rgdal)
library(latex2exp)

source("Functions/STACYmap_5.R")
source("Functions/projection_ptlyr.R")


var_map_plot <- function(Point_Lyr, 
                         projection = as.character('+proj=robin +datum=WGS84'),
                         pt_size,
                         txt_size){
  
  Point_Lyr_p <- projection_ptlyr(Point_Lyr, projection)
  
  allmax = max(abs(Point_Lyr$value))
  
  plot <- STACYmap(coastline = TRUE) +
    geom_point(data = Point_Lyr_p, aes(x = long, y = lat, fill = layer), color = "black", shape = 21,
               size = (pt_size), alpha = 1.0, show.legend = c(fill =TRUE)) +
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, 'RdBu')), 
                         limits = c(log(0.025),log(40)),
                         breaks = c(log(0.025), log(0.05), log(0.25), log(0.5), 0, log(2), log(4), log(20), log(40)),
                         labels = c(0.025, "", 0.25, "", 1, "", 4, "", 40),
                         name = TeX("$Var_{Rec}/Var_{Sim}$")) +
    theme(legend.direction = "horizontal", 
          panel.border = element_blank(),
          legend.background = element_blank(),
          axis.text = element_blank(),
          text = element_text(size = txt_size),
          legend.title = element_text(size = txt_size))
  
  return(plot)
}

