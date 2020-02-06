#################################################
## Plot Karst Map with SISAL sites ############## ------------------------------------
#################################################


coastline_map <- rgdal::readOGR(dsn ="/home/ginnyweasley/Dokumente/01_Promotion/07_R_Code/naturalearth_10m_physical/ne_10m_coastline.shp", verbose = FALSE)
karst_map <- rgdal::readOGR(dsn ="/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/07_Karst/karst_wgs.shp", verbose = FALSE)

ALL_SITES <- read.csv("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/02_SISAL/SISAL_v2_CARLA/site_countries.csv")
USED_SITES <- ALL_SITES %>% filter(site_id %in% DATA_past1000$CAVES$site_info$site_id) %>% distinct(site_id, longitude, latitude)
USED_SITES <- data.frame(
  lon = USED_SITES$longitude,
  lat = USED_SITES$latitude,
  value = USED_SITES$site_id
)
NOT_SITES <- ALL_SITES %>% filter(!site_id %in% DATA_past1000$CAVES$site_info$site_id) %>% distinct(site_id, longitude, latitude)
NOT_SITES <- data.frame(
  lon = NOT_SITES$longitude,
  lat = NOT_SITES$latitude,
  value = NOT_SITES$site_id
)


#Make Projection for points and karst_lyr
projection = as.character('+proj=robin +datum=WGS84')
karst_map2 <- karst_map %>% spTransform(., CRSobj = CRS(projection)) %>% fortify(.)
coastline_map2 <- coastline_map %>% spTransform(., CRSobj = CRS(projection)) %>% fortify(.)
source("Functions/projection_ptlyr.R")

USED_SITES_p <- projection_ptlyr(USED_SITES, projection)
NOT_SITES_p <- projection_ptlyr(NOT_SITES, projection)

remove(ALL_SITES, USED_SITES, NOT_SITES)

plot <- STACYmap(coastline = TRUE, filledbg = TRUE) +
  geom_polygon(data = karst_map2, aes(x=long, y = lat, group = group), fill = '#A57F7C', color = NA) +
  new_scale_color() +
  geom_point(data = NOT_SITES_p, aes(x = long, y = lat, shape = "20", colour = '#3C46B2'),# shape = 20,
             size = (GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE-1), show.legend = c(shape =TRUE)) +
  geom_point(data = USED_SITES_p, aes(x = long, y = lat, shape = "17", color = '#B2463C'),# shape = 17,
             size = (GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE), show.legend = c(shape =TRUE)) +
  scale_colour_manual(name = "SISALv2 available sites", labels = c("unused sites", "sites with past 1000y records"), values = c('#3C46B2','#B2463C')) +
  scale_shape_manual(name = "SISALv2 available sites",labels = c("unused sites", "sites with past 1000y records"), values = c(20,17)) +
  #geom_polygon(data = coastline_map2, aes(x=long, y=lat, group = group), color = 'black',  size = 0.2, fill = NA, alpha = 0.8) +
  theme(legend.position = c(0.02, 0.02),
        legend.justification = c(0, 0),
        legend.box = 'vertical',
        legend.box.background = element_blank(),
        legend.background = element_rect(colour = 'black'))

plot %>% ggsave(filename = paste('map_sisal_used_sites', 'pdf', sep = '.'), plot = ., path = 'Plots', 
                width =20, height = 12.5, units = 'cm', dpi = 'print')

