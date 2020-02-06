library(ggplot2)
library(ggnewscale)
library(maps)
library(marmap)
library(rgdal)
library(sp)
library(RColorBrewer)
library(magrittr)
library(knitr)

# plotting options
GLOBAL_FONT_FACE_TITLE <- 'bold'
GLOBAL_FONT_FACE_TEXT <- 'plain'
GLOBAL_FONT_SIZE <- 17
GLOBAL_FONT_FAMILY <- 'sans'

GLOBAL_BLUE_LIGHT <- viridis(1, begin = 0.6)
GLOBAL_BLUE_LIGHT_ALPHA_LOW <- viridis(1, alpha = 0.45, begin = 0.6)
GLOBAL_BLUE_DARK <- viridis(1, begin = 0.4)
GLOBAL_GREEN_LIGHT <- viridis(1, begin = 0.9)
GLOBAL_GREEN_LIGHT_ALPHA_LOW <- viridis(1, alpha = 0.45, begin = 0.9)
GLOBAL_GREEN_DARK <- viridis(1, begin = 0.7)
GLOBAL_GREY_LIGHT <- grey.colors(1, start = 0.97, end = 0.99)
GLOBAL_GREY_LIGHT_ALPHA_LOW <- grey.colors(1, start = 0.97, end = 0.99, alpha = 0.45)
GLOBAL_GREY_DARK <- grey.colors(1, start = 0.6)
GLOBAL_GREY_MEDIUM <- grey.colors(1, start = 0.88)
GLOBAL_GREY_MEDIUM_LIGHT <- grey.colors(1, start = 0.95)
GLOBAL_RED_DARK <- viridis(1, begin = 0.6, option = 'B')
GLOBAL_RED_LIGHT <- viridis(1, begin = 0.7, option = 'B')
GLOBAL_BLUE_DDARK <- viridis(1, begin = 0.2, option = 'E')
GLOBAL_BLUE_DLIGHT <- viridis(1, begin = 0.4, option = 'E')
GLOBAL_ORANGE_DARK <- viridis(1, begin = 0.75, option = 'A')
GLOBAL_ORANGE_LIGHT <- viridis(1, begin = 0.85, option = 'A')
GLOBAL_YELLOW_DARK <- viridis(1, begin = 0.9, option = 'E')
GLOBAL_YELLOW_LIGHT <- viridis(1, begin = 0.98, option = 'E')
GLOBAL_VIOLET_DARK <- viridis(1, begin = 0.4, option = 'C')
GLOBAL_VIOLET_LIGHT <- viridis(1, begin = 0.5, option = 'C')

LABELED_SITES_FOR_LAT_ORDERED_PLOT <- c(10, 7, 19, 90, 29, 82, 17, 38, 28, 15)
SAMPLING_RATE_INTERVALS_SIMPLE <- list(list(start = 0, stop = 8000), list(start = 15000, stop = 73000))
SAMPLING_RATE_INTERVALS_DEFAULT <- list(list(start = 0, stop = 8000), list(start = 8001, stop = 18999),
                                        list(start = 19000, stop = 27000))
SAMPLING_RATE_INTERVALS_10k_WINDOW <- list(list(start = 0, stop = 10000), list(start = 10000, stop = 20000),
                                           list(start = 20000, stop = 30000), list(start = 30000, stop = 40000), 
                                           list(start = 40000, stop = 50000), list(start = 50000, stop = 60000), 
                                           list(start = 60000, stop = 70000))

# for further projections see https://proj4.org/operations/projections/index.html
GLOBAL_CRS <-list(fixed = 'fixed', robinson = '+proj=robin', wintri = '+proj=wintri', azequidistant = '+proj=aeqd', area = '+proj=aea +lat_1=29.5 +lat_2=42.5')

plot_CAVES_sites_on_map <- function(sites, sample_dating,
                                   bathy = TRUE, 
                                   topo = FALSE,
                                   pnv = FALSE, 
                                   lgm_ice = FALSE, 
                                   names = FALSE,
                                   site_id = FALSE, 
                                   ISR = FALSE, ISR_stat = 'med',
                                   save_plot = NULL, 
                                   legend_inside = FALSE, 
                                   force_shapefiles = FALSE, 
                                   projection = 'fixed',
                                   zoom = NULL, 
                                   ...) {
  if (ISR) {if (!(ISR_stat %in% c('med', 'mean'))) {stop('unknown ISR_stat')}}
  if (!(projection %in% names(GLOBAL_CRS))) {
    stop('unknown projection')
  } else if (!(projection == 'fixed')) {
    sites <- arrange(sites, site_id)
    sites_transformed <- as_tibble(project(cbind(sites$long, sites$lat), proj = GLOBAL_CRS[[projection]])) %>% 
      rename(long_new = V1, lat_new = V2) %>% 
      bind_cols(select(sites, site_id))
    
    sites <- sites %>% 
      full_join(sites_transformed,
                by = 'site_id') %>% 
      select(-lat, -long) %>% 
      rename(long = long_new, lat = lat_new)
    if (!is.null(zoom)) {
      yaxislclip <- zoom[2]; xaxislclip <- zoom[1]; xaxisrclip <- zoom[3]
      zoom <- c(zoom[1], 0, 0, zoom[2], zoom[3], 0, 0, zoom[4])
      zoom <- project(matrix(c(unlist(zoom)), nrow = 4, byrow = T), proj = GLOBAL_CRS[[projection]]) #;zoom <- project(matrix(c(unlist(zoom)), nrow = 2, byrow = T), proj = GLOBAL_CRS[[projection]])
      zoom <- c(zoom[1,1], zoom[2,2], zoom[3,1], zoom[4,2])
      zoom <- matrix(zoom, nrow = 2, byrow = T)
    } else {
      yaxislclip <- -90; xaxislclip <- -180; xaxisrclip <- 180
    }
    
    cat(paste('transformed locations to projection:', projection, '\n'))
  }
  #data <- get_core_existance(sites, core_types = list('pollen' = pollen_data, 'charcoal' = charcoal_data, 'biome_percentages' = biome_percentages)) %>% 
  data <- inner_join(sites, compute_age_interval_stats(samples = sample_dating,
                                                       ranges = list(list(start = 0, stop = 65000))) %>% 
                       categorize_age_interval_stats(., ...) %>% 
                       select(-oldest_dp_in_interv, -youngest_dp_in_interv, -status_smp_per_10k, -status_ev_res_in_interv),
                     by = 'site_id')
  #data$interv <- ordered(data$interv, levels = c('0-8000', '19000-27000', '0-73000'))
  
  if (pnv | force_shapefiles | projection != 'fixed') {
    map_plot <- shapefile_map_plot_(ggobj = ggplot(), projection = projection, 
                                    graticules = 30, topo = FALSE, bathy = bathy,
                                    pnv = pnv, lgm_ice = lgm_ice, cb_friendly = TRUE, 
                                    yaxislclip = yaxislclip, xaxisrclip = xaxisrclip, xaxislclip = xaxislclip)
  }else {
    map_plot <- bathymetric_map_plot_(ggobj = ggplot(),
                                      plot_bathy = bathy, resolution = 15, 
                                      use_contours = FALSE, plot_topo = topo)
  }
  
  if(names){
    map_plot <- map_plot + geom_label_repel(data = data,
                                            mapping = aes(x=long, y=lat, label = paste(site_id, ' - ', site_name, sep = '')), nudge_y = 2.5, show.legend = FALSE,
                                            size = GLOBAL_FONT_SIZE - 10, fill = 'black', colour='white', label.padding = 0.15, label.r = 0,
                                            min.segment.length = 0, segment.color = 'black')
  }
  if(site_id){
    if(names){
      print('name labels overwrite site id labels')
    }
    else{
      map_plot <- map_plot + geom_label_repel(data = data,
                                              mapping = aes(x=long, y=lat, label = site_id), nudge_y = 2.5, show.legend = FALSE,
                                              size = GLOBAL_FONT_SIZE - 10, fill = 'black', colour='white', label.padding = 0.15, label.r = 0,
                                              min.segment.length = 0, segment.color = 'black')
    }
  }
  
  if (ISR) {
    if (pnv) {
      map_plot <- map_plot + 
        geom_point(data = data,
                   mapping = aes(x = long, y = lat, shape = {if (ISR_stat == 'med') {status_med_smp_res} else {status_mean_smp_res}}), size = 4.5,
                   alpha = I(0.7),
                   fill = grey.colors(1, start = 0.94, gamma = 0.1)) + 
        scale_shape_manual(guide = global_legend(title = {if (ISR_stat == 'med') {'median ISR [y]\n0-65 kyr BP'} else {'mean ISR [y]\n0-65 kyr BP'}}, order = 2), 
                           values = c(24, 23, 21))
    } else {
      map_plot <- map_plot + 
        geom_point(data = data,
                   mapping = aes(x = long, y = lat,
                                 shape = {if (ISR_stat == 'med') {status_med_smp_res} else {status_mean_smp_res}}, 
                                 fill = {if (ISR_stat == 'med') {status_med_smp_res} else {status_mean_smp_res}}), 
                   size = 5,
                   alpha = I(0.6)) + 
        scale_shape_manual(guide = global_legend(title = {if (ISR_stat == 'med') {'median ISR [y]\n0-65 kyr BP'} else {'mean ISR [y]\n0-65 kyr BP'}}, order = 2), 
                           values = c(24, 23, 21)) + 
        scale_fill_manual(guide = global_legend(title = {if (ISR_stat == 'med') {'median ISR [y]\n0-65 kyr BP'} else {'mean ISR [y]\n0-65 kyr BP'}}, order = 2), 
                          values = c(GLOBAL_GREEN_DARK, GLOBAL_YELLOW_DARK, GLOBAL_BLUE_DDARK))
    }
  } else {
    map_plot <- map_plot + 
      geom_point(data = data,
                 mapping = aes(x = long, y = lat), size = 3.5, alpha=I(0.6)) 
  }
  if (!is.null(zoom)) {
    map_plot <- map_plot + coord_equal(xlim = zoom[,1], ylim = zoom[,2])
  }
  if(legend_inside){map_plot <- map_plot + 
    theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0), legend.box = 'vertical')
  }
  if(!is.null(save_plot)){
    global_save_plot(save_plot = save_plot, plot = map_plot)
  }
  return(map_plot)
}

###############################################################################################################
shapefile_map_plot_ <- function(ggobj, 
                                projection = 'robinson', 
                                graticules = 30, 
                                topo = FALSE, 
                                bathy = FALSE,
                                pnv = FALSE, 
                                lgm_ice = FALSE, 
                                black = FALSE, 
                                cb_friendly = TRUE, 
                                bg = NULL,
                                yaxislclip = -90,
                                xaxisrclip = 180, 
                                xaxislclip = -180) {
  if (pnv & projection != 'fixed') {warning('using fixed projection for pnv = TRUE'); projection <- 'fixed'; topo <- FALSE}
  if (!exists(x = 'shapefile_map_plot_data', envir = .GlobalEnv)) {load <- TRUE; exists <- FALSE}
  else if(exists(x = 'shapefile_map_plot_data', envir = .GlobalEnv) & !(projection %in% names(shapefile_map_plot_data))) {load <- TRUE; exists <- TRUE}
  else {load <- FALSE}
  
  
  #if (topo & !exists('shapefile_topo_data', envir = .GlobalEnv)) {load_topo <- TRUE; exists_topo <- FALSE}
  #else if (topo & exists(x = 'shapefile_topo_data', envir = .GlobalEnv) & !(projection %in% names(shapefile_map_plot_data))) {load_topo <- TRUE; exists_topo <- TRUE}
  #else {load_topo <- FALSE}
  load_topo <- FALSE
  
  if (load) {
    cat(cr('caching shapefile_map_plot_data to .GlobalEnv for projection', projection, '\n'))
    map_data <- load_natural_earth_data(file = 'ne_10m_land.shp') %>% 
      transform_shapefile(data = ., transform = GLOBAL_CRS[[projection]]) %>% 
      fortify(.)
    ocean_data <- load_natural_earth_data(file = 'ne_10m_ocean.shp') %>% 
      transform_shapefile(data = ., transform = GLOBAL_CRS[[projection]]) %>% 
      fortify(.)
    dupl_id <- {if (is.logical(graticules)) {180/30-1} else if (90 %% graticules == 0) {180/graticules-1} else {180/graticules}}
    grat <- load_natural_earth_data(file = paste('ne_10m_graticules_', {if (is.logical(graticules)) {30}else {graticules}}, '.shp', sep = '')) %>% 
      transform_shapefile(data = ., transform = GLOBAL_CRS[[projection]]) %>% 
      fortify(.) %>% 
      mutate_at(vars(id), as.integer)
    grat <- grat %>% mutate_at(vars(group), as.character) %>% 
      bind_rows(grat %>% filter(id == dupl_id) %>%
                  mutate(long=-long, id = max(grat$id) + 1, group = as.character(max(grat$group %>% levels() %>% as.numeric()) + 1)))
    bbox <- load_natural_earth_data(file = 'ne_10m_wgs84_bounding_box.shp') %>% 
      transform_shapefile(data = ., transform = GLOBAL_CRS[[projection]]) %>% 
      fortify(.)
    bathy_data <- lapply(sort(c(200, seq(0, 10000, 1000))), function (x, L) {
      load_natural_earth_data(file = paste0('ne_10m_bathymetry_', L[[as.character(x)]], '_', x, '.shp')) %>% 
        transform_shapefile(data = ., transform = GLOBAL_CRS[[projection]]) %>% 
        fortify(.)}, L = LETTERS[seq(1, 12, 1) %>% rev()] %>% setNames(., c(200, seq(0, 10000, 1000)) %>% sort())) %>% setNames(., sort(c(200, seq(0, 10000, 1000))))
    lgm_ice_data <- load_natural_earth_data(dir = '../data_supplementary/lgm_18k_ice_extent', file = 'lgm.shp') %>% 
      transform_shapefile(data = ., transform = GLOBAL_CRS[[projection]]) %>% 
      fortify(.)
    sfile_map_plot_data <- list(); sfile_map_plot_data[[projection]] <- list(map_data = map_data, ocean_data = ocean_data, bathy_data = bathy_data, grat = grat, bbox = bbox, lgm_ice_data = lgm_ice_data)
    if(!exists) {shapefile_map_plot_data <<- sfile_map_plot_data}else {shapefile_map_plot_data <<- merge.list(shapefile_map_plot_data, sfile_map_plot_data)}
  }
  else {
    map_data <- shapefile_map_plot_data[[projection]]$map_data; ocean_data <- shapefile_map_plot_data[[projection]]$ocean_data
    grat <- shapefile_map_plot_data[[projection]]$grat; bbox <- shapefile_map_plot_data[[projection]]$bbox
    bathy_data <- lapply(as.vector(names(shapefile_map_plot_data[[projection]]$bathy_data)), function(d, x){as_tibble(x[[d]]) %>% mutate(depth = d)},
                         x = shapefile_map_plot_data[[projection]]$bathy_data) %>% setNames(., sort(c(200, seq(0, 10000, 1000))))
    if (lgm_ice) {lgm_ice_data <- shapefile_map_plot_data[[projection]]$lgm_ice_data}
  }
  
  if (load_topo) {
    stop('topo not allowed at the moment')
    cat('caching shapefile_topo to .GlobalEnv for projection', projection, '\n')
    topo_data <- raster::stack(file.path(NAT_EARTH_DATA_PATH, 'HYP_LR', 'HYP_LR.tif')) %>% as(., 'SpatialPixelsDataFrame') %>% as.data.frame(); colnames(topo_data) <- c('z', 'x', 'y')
    #topo_data <- data.frame(raster::xyFromCell(topo_data, 1:raster::ncell(topo_data)), raster::getValues(topo_data/255))
    #topo_data$rgb <- with(topo_data, rgb(NE2_HR_LC_SR_W_DR.1, NE2_HR_LC_SR_W_DR.2, NE2_HR_LC_SR_W_DR.3, 1))
  }
  
  gobj <- ggobj + 
    geom_polygon(data = bbox, mapping = aes(x = long, y = lat, group = group), fill = GLOBAL_GREY_DARK)
  if (topo) {
    ggobj <- ggobj + 
      geom_tile(data = topo_data, mapping = aes(x = x, y = y, fill = z))
  }
  else if (bathy) {
    gobj <- gobj + 
      geom_polygon(data = bathy_data[['0']], mapping = aes(x = long, y = lat, group = group), fill = '#d9ebf9', na.rm = T) + 
      geom_polygon(data = bathy_data[['200']], mapping = aes(x = long, y = lat, group = group), fill = '#cae1f4', na.rm = T) + 
      geom_polygon(data = bathy_data[['1000']], mapping = aes(x = long, y = lat, group = group), fill = '#afd3ef', na.rm = T) + 
      geom_polygon(data = bathy_data[['2000']], mapping = aes(x = long, y = lat, group = group), fill = '#aacde9', na.rm = T) + 
      geom_polygon(data = bathy_data[['3000']], mapping = aes(x = long, y = lat, group = group), fill = '#96c1e3', na.rm = T) + 
      geom_polygon(data = bathy_data[['4000']], mapping = aes(x = long, y = lat, group = group), fill = '#83b9df', na.rm = T) + 
      geom_polygon(data = bathy_data[['5000']], mapping = aes(x = long, y = lat, group = group), fill = '#6fadd6', na.rm = T) + 
      geom_polygon(data = bathy_data[['6000']], mapping = aes(x = long, y = lat, group = group), fill = '#5ba2d0', na.rm = T) + 
      geom_polygon(data = bathy_data[['7000']], mapping = aes(x = long, y = lat, group = group), fill = '#589cc9', na.rm = T) + 
      geom_polygon(data = bathy_data[['8000']], mapping = aes(x = long, y = lat, group = group), fill = '#337fb2', na.rm = T) + 
      geom_polygon(data = bathy_data[['9000']], mapping = aes(x = long, y = lat, group = group), fill = '#2a77ac', na.rm = T) + 
      geom_polygon(data = bathy_data[['10000']], mapping = aes(x = long, y = lat, group = group), fill = '#2371a6', na.rm = T)
    if (!pnv) {
      gobj <- gobj + 
        geom_polygon(data = map_data, mapping = aes(x = long, y = lat, group = group, fill = hole)) + 
        scale_fill_manual(values = c('#f0e6c2', NA), guide = FALSE)
    }
  }else if(!pnv & is.null(bg)) {
    #map_data <- mutate_at(map_data, vars(group, id), as.numeric)
    #bathymap <- bathy_data[['0']] %>% mutate_at(vars(group, id), as.numeric) %>% 
    #  mutate(id = max(map_data$id) + 1, group = max(map_data$group) + 1) %>% 
    #  bind_rows(map_data)
    #bathy$group <- bathy$group
    gobj <- gobj + 
      #geom_polygon(data = bathymap, mapping = aes(x = long, y = lat, group = group, fill = hole), na.rm = T) + 
      geom_polygon(data = bathy_data[['0']], mapping = aes(x = long, y = lat, group = group,  fill = TRUE, color = TRUE)) + 
      geom_polygon(data = map_data, mapping = aes(x = long, y = lat, group = group, fill = hole, color = hole)) + 
      scale_fill_manual(values = c('#f0e6c2', '#83b9df'), guide = FALSE) + 
      scale_color_manual(values = c('#f0e6c2', '#83b9df'), guide = FALSE) + 
      new_scale_fill() + new_scale_color() # grey scale: scale_fill_manual(values = c(GLOBAL_GREY_LIGHT, GLOBAL_GREY_DARK), guide = FALSE) 
  }
  
  if (pnv & is.null(bg)) {
    pnv <- load_pnv_nc_file(NAME_POT_VEG_DISTR)
    gobj <- gobj + 
      geom_polygon(data = filter(map_data, lat < -60), mapping = aes(x = long, y = lat, group = group), fill = grey.colors(1, start = 0.94, gamma = 0.1)) +  
      geom_tile(data = pnv,
                mapping = aes(x=lat, y=lon, fill=megabiome, color=megabiome),
                show.legend = c(color = FALSE, fill = TRUE))
    if (cb_friendly) {
      gobj <- gobj + 
        scale_fill_viridis_d(begin = 0.2, guide = global_legend_horiz('megabiome', ncol = 5)) + 
        scale_color_viridis_d(begin = 0.1, guide = global_legend_horiz('megabiome', ncol = 5))
    }else {
      gobj <- gobj + 
        scale_fill_manual(values = PNV_COLORING) + 
        scale_color_manual(values = PNV_COLORING)
    }
    gobj <- gobj + 
      new_scale('fill') + new_scale('colour')
  }
  
  if (lgm_ice & is.null(bg)) {
    gobj <- gobj + 
      geom_polygon(data = lgm_ice_data, mapping = aes(x = long, y = lat, group = group), fill = NA, #fill = grey.colors(1, 0.96, gamma = 0.2),
                   alpha = 0.15, color = grey.colors(1, 0.2, gamma = 0.1), size = 0.45)
  }
  
  if (!is.null(bg)) {
    if (projection == 'fixed') {
      bg <- bg %>% as_tibble(); colnames(bg) <- famous_lat
      bg <- tibble(long = famous_lon) %>% bind_cols(bg) %>% gather(key = 'lat', value = 'slope', -long) %>% mutate_at(vars(lat), as.numeric)
      colnames(bg) <- c(colnames(bg)[1:2], 'z')
    } else if (projection != 'fixed') {
      print('#')
      bg <- raster::raster(t(bg), crs = "+proj=longlat +datum=WGS84", xmn = -180, xmx = 180, ymn = -90, ymx = 90) %>% raster::flip('y') %>% 
        raster::rasterToPolygons() %>% 
        spTransform(GLOBAL_CRS[[projection]]) %>% as('sf')
      #bg <- fortify(bg)
      #bg <- tibble(long = coordinates(bg)[, 1], lat = coordinates(bg)[, 2], z = bg@data$layer)
    }
    #print(bg)
    #stop()
    print(bg)
    gobj <- ggplot() + 
      geom_sf(data = bg, mapping = aes(color = layer, fill = layer), ) + coord_sf() #+ 
    #geom_tile(data = bg,
    #            mapping = aes(x=long, y=lat, fill=z, colour = z), width=600000, height=300000,
    #            show.legend = c(fill = FALSE)) +
    #  scale_fill_viridis_c(begin = 0.2) + 
    geom_polygon(data = map_data, mapping = aes(x = long, y = lat, group = group), fill = NA, colour = 'black', size = 0.25)
  }
  
  if (is.logical(graticules)) {
    xgrats <- c(-90, 0, 90); ygrats <- c(-45, 0, 45); xgrats_names <- c(-90, 0, 90); ygrats_names <- c(-45, 0, 45)
  } else {
    xgrats <- seq(-180, 180, graticules) %>% ifelse(. >= xaxisrclip | . <= xaxislclip, NA, .) %>% .[!is.na(.)]
    ygrats <- seq(-90, 90, graticules)
    xgrats_names <- xgrats
    ygrats_names <- ygrats
  }
  if (!projection == 'fixed') {
    graticules <- tibble(x = c(xgrats, rep(0, length(ygrats))), y = c(rep(yaxislclip, length(xgrats)), ygrats))
    graticules <- as_tibble(project(cbind(graticules$x, graticules$y), proj = GLOBAL_CRS[[projection]])) %>% 
      rename(x = V1, y = V2)
    xgrats <- graticules$x[1:length(xgrats)]; ygrats <- graticules$y[(length(xgrats) + 1):(length(ygrats) + length(xgrats))]
  }
  
  if (!is.logical(graticules)) {plot_grats <- TRUE} else {if (graticules) {plot_grats <- FALSE} else {plot_grats <- FALSE}}
  
  gobj <- gobj + 
    #geom_polygon(data = ocean_data, mapping = aes(x = long, y = lat, group = group), fill = GLOBAL_GREY_DARK) + 
  {if (!plot_grats){} else {geom_path(data = grat, aes(long, lat, group=group, fill = NULL), linetype="dashed", color="grey50")}} + 
    scale_y_continuous(breaks = ygrats, labels = ygrats_names, expand = expand_scale(add = 0.01)) + 
    scale_x_continuous(breaks = xgrats, labels = xgrats_names, expand = expand_scale(add = 0.01)) + 
    coord_equal() + 
    global_title_and_axis() + 
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(), 
          legend.position = 'bottom', legend.box = 'horizontal') + 
          {if (black) {theme(rect = element_rect(fill = 'black'), text = element_text(colour = 'white'),
                             panel.background = element_rect('black'), panel.border = element_rect(color = 'white'),
                             legend.background = element_rect(color = 'white'), strip.text = element_text(colour = 'white'))}}
  return(gobj)
}