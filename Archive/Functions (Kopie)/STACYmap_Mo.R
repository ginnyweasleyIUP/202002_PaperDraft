library(tidyverse)
library(ggnewscale)
library(rgdal)
library(RColorBrewer)
library(sp)
library(sf)
library(viridisLite)

# environment for data caching
.STACYmap <- environment()

# helpers
title_and_axis <- function(){
  #styles
  
  GLOBAL_FONT_FACE_TITLE <- 'bold'
  GLOBAL_FONT_FACE_TEXT <- 'plain'
  GLOBAL_FONT_SIZE <- 17
  GLOBAL_FONT_FAMILY <- 'sans'
  
  ax <- theme_bw() + 
    theme(title = element_text(face = GLOBAL_FONT_FACE_TITLE, size = GLOBAL_FONT_SIZE + 1, family = GLOBAL_FONT_FAMILY),
          text = element_text(face = GLOBAL_FONT_FACE_TEXT, size = GLOBAL_FONT_SIZE, family = GLOBAL_FONT_FAMILY), 
          #line = element_blank(),
          legend.direction = 'vertical', legend.box = 'vertical', legend.box.background = element_blank(), legend.background = element_rect(colour = 'black'), 
          strip.background.x = element_blank(), strip.text.x = element_text(face = GLOBAL_FONT_FACE_TITLE, hjust = 0), 
          strip.background.y = element_rect(fill = 'white'), strip.text.y = element_text(face = GLOBAL_FONT_FACE_TITLE, vjust = 0))
  return(ax)
}


load_natural_earth_data <- function(file, dir = '/home/ginnyweasley/Dokumente/01_Promotion/07_R_Code/naturalearth_10m_physical', ...) {
  NAT_EARTH_PATH <- '/home/ginnyweasley/Dokumente/01_Promotion/07_R_Code/naturalearth_10m_physical'
  if (dir == NAT_EARTH_PATH & str_detect(file, '.shp')) {
    data <- readOGR(dsn = file.path(NAT_EARTH_PATH, file), verbose = FALSE, ...)
  } else if (str_detect(file, '.shp')) {
    data <- readOGR(dsn = file.path(dir, file), verbose = FALSE, ...)
  } else if (str_detect(file, '.tif')) {
    data <- raster::raster(x = file.path(NAT_EARTH_PATH, file))
  }
  return(data)
}


transform_shapefile <- function(data, transform) {
  if (transform == 'fixed') {
    return(data)
  }
  else {
    return(spTransform(data, CRSobj = CRS(transform)))
  }
}


#' Plot appealing maps of latlon-gridded and point-wise data
#'
#' @param gridlyr latlon-grid to plot of class raster, matrix (both lon columns, lat rows, cell value) or data.frame/tibble (either lon columns, lat rows, cell value or lon, lat, cell columns)
#' @param plgnlyr NOT YET IMPLEMENTED latlon-polygon layer, can be a raster object (will be converted as is to SpatialPolygonsDataFrame) or SpatialPolyGonsDataFrame or sf, will be plotted on top of gridlyr
#' @param ptlyr point layer of class matrix or data.frame (lon, lat, cell columns)
#' @param coastline logical - plot coastline?
#' @param bathymetry logical - plot ocean bathymetry?, option can be used independently of coastlines and will underlay a climate field
#' @param projection object of class CRS or character (is converted to CRS), see package rgdal
#' @param graticules logical, should a projected panel.grid be plotted?
#' @param zoom numeric, c(west long, south lat, east long, north lat) to be cut off
#' @param legend_names list(grid = nm1, pt = nm2)
#' @param colorscheme color palette(s) or character from 'temp', 'prcp_div', 'prcp_grd', 'spec' - pass one color palette to be respected by grid and points or list(grid = pal1, pt = pal2) for independent color schemes
#' @param legend_inside logical - should thelegend be placed inside the map plot?
#' @param revcolorscheme logical - has effect only on built-in schemes
#' @param flipy logical - flip the grid input in y direction, try option if map orientation is wrong
#' @param transpose logical - transpose the grid input, try option if map orientation is wrong
#' @param box logical - should the entire plot be bounded by a box?
#' @param filledbg logical - should land and ocean surface be coloured?, option can be used independently of coastlines and will underlay a climate field
#' @param centercolor numeric - or list(grid = val1, pt = val2) of numeric to center colorscheme to a specific value 
#'
#' @return
#' @export
#'
#' @examples 
#' # plotting a field and points
#' nc <- nc_open(PLASIM_TESTSLICE_PI)
#' testfld <- ncvar_get(nc, varid = "ts")
#' lat <- ncvar_get(nc, varid = "lat")
#' lon <- ncvar_get(nc, varid = "lon")
#' nc_close(nc)
#' testfld <- testfld[,,1]
#' ptlyr <- tibble(long = c(65, 30), lat = c(45, -60), value = c(300, 260))
#' STACYmap(gridlyr = testfld, ptlyr = ptlyr, colorscheme = 'temp',
#'          graticules = T, zoom = c(-100, -45, 100, 45), 
#'          legend_names = list(grid = 'temp', pt = 'temp'))
#' 
#' # plotting points on a world map
#' STACYmap(ptlyr = ptlyr, colorscheme = 'temp', filledbg = T, bathymetry = T)
#' 
#' # using other projections than the default robinson
#' # disclaimer: some combinations of projection and raster plots are not working yet
#' plt <- STACYmap(ptlyr = ptlyr, colorscheme = 'temp', filledbg = T, projection = '+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000')
STACYmap <- function(gridlyr = NULL, 
                     ptlyr = NULL, 
                     coastline = TRUE,
                     filledbg = FALSE,
                     bathymetry = FALSE,
                     projection = '+proj=robin +datum=WGS84', #+proj=longlat +datum=WGS84
                     graticules = TRUE, 
                     zoom = NULL, 
                     legend_names = list(grid = 'grid', pt = 'points'), 
                     legend_inside = FALSE,
                     colorscheme = rev(RColorBrewer::brewer.pal(11, 'RdBu'))[2:10],
                     centercolor = NULL,
                     revcolorscheme = F,
                     flipy = FALSE, 
                     transpose = TRUE, 
                     box = TRUE) {
  
  # style options
  GLOBAL_CRS <- list(equirectangular = 'fixed',
                     robinson = '+proj=robin',
                     wintri = '+proj=wintri', 
                     azequidistant = '+proj=aeqd',
                     area = '+proj=aea +lat_1=29.5 +lat_2=42.5')
  GLOBAL_GREY_LIGHT <- grey.colors(1, start = 0.97, end = 0.99)
  GLOBAL_GREY_LIGHT_ALPHA_LOW <- grey.colors(1, start = 0.97, end = 0.99, alpha = 0.45)
  GLOBAL_GREY_DARK <- grey.colors(1, start = 0.6)
  GLOBAL_GREY_MEDIUM <- grey.colors(1, start = 0.88)
  GLOBAL_GREY_MEDIUM_LIGHT <- grey.colors(1, start = 0.95)
  GLOBAL_FONT_FACE_TITLE <- 'bold'
  GLOBAL_FONT_FACE_TEXT <- 'plain'
  GLOBAL_FONT_SIZE <- 17
  GLOBAL_FONT_FAMILY <- 'sans'
  GLOBAL_POINT_SIZE <- 4
  GLOBAL_BATHYMETRY_COLORS <- c('#d9ebf9', '#cae1f4', '#afd3ef', '#aacde9', '#96c1e3', '#83b9df', '#6fadd6', '#5ba2d0', '#589cc9', '#337fb2', '#2a77ac', '#2371a6')
  GLOBAL_LAND_COLOR <- '#f0e6c2'
  GLOBAL_OCEAN_COLOR <- '#aacde9'
  
  # TODO enable lapply -> option to pass multiple gridlyr, plgnlyr, ptlyr as a list
  # projections
  #if (class(projection) == 'character') {
  #  projection <- CRS(projection)
  #}
  splcol <- FALSE; if (all(class(colorscheme) == 'list')) {splcol <- TRUE}
  if ((class(colorscheme) == 'list' & all(lapply(colorscheme[1:length(colorscheme)], length) == 1)) | (class(colorscheme) == 'character' & length(colorscheme) == 1)) {
    colorschemeo <- colorscheme
    colorscheme <- lapply(colorscheme, 
                          function (x) {
                            if (is.character(x) & length(x) == 1) {
                              x <- list(
                                temp = rev(RColorBrewer::brewer.pal(11, 'RdBu'))[2:10],
                                prcp_grd = RColorBrewer::brewer.pal(9, 'YlGnBu'), 
                                prcp_div = rev(RColorBrewer::brewer.pal(11, 'BrBG')),
                                spec = rev(RColorBrewer::brewer.pal(11, 'RdBu'))
                              )[[x]]
                              
                              if (is.null(x)) {stop('unknwon colorscheme given to STACYmap')}
                              if (revcolorscheme) {x <- rev(x)}
                              return(x)
                            } else {
                              return(x)
                            }
                          })
    if (class(colorschemeo) == 'character' & length(colorschemeo) == 1) {
      colorscheme <- unlist(colorscheme)
    }
  } else if (any(lapply(colorscheme[1:length(colorscheme)], length) != 1)) {
    stop('mixing pre-defined and built-in colorschemes is not supported, please choose one option for all layers')
  }
  ## checks and data preparation/projection
  #if (length(colorscheme) == 1) {colorscheme <- unlist(colorscheme)}
  
  xlim <- NULL; ylim <- NULL
  if (!is.null(zoom)) {
    zoom <- matrix(c(zoom[1], 0, zoom[3], 0, 0, zoom[2], 0, zoom[4]), nrow = 4)
    zoom <- project(zoom,
                    proj = as.character(projection))
    xlim <- c(zoom[1,1], zoom[3,1])
    ylim <- c(zoom[2,2], zoom[4,2])
  }
  
  if (!is.null(gridlyr)) {
    if (!any(class(gridlyr) %in% c('raster', 'matrix', 'data.frame'))) {
      stop('class(gridlyr) has to be one of \'raster\', \'matrix\', \'data.frame\'')
    } else {
      prepg <- list(raster = function(x) x, 
                    matrix = function(x) 
                      raster::raster({if (transpose) {t(x)} else {x}}, 
                                                        crs = "+proj=longlat +datum=WGS84", xmn = -180, xmx = 180, ymn = -90, ymx = 90) %>% 
                      {if (flipy) {raster::flip(., 'y')} else {.}},
                    data.frame = function(x) as.matrix(x) %>% prepg$raster(.)
                    )
      
      gridlyr <- gridlyr %>% 
        prepg[[class(gridlyr)]](.) %>% 
        raster::rasterToPolygons() %>% 
        spTransform(CRS(projection)) %>%
        as('SpatialPolygonsDataFrame') %>% 
        as('sf')
    }
    
  }
  if (!is.null(ptlyr)) {
    if (!any(class(ptlyr) %in% c('matrix', 'data.frame'))) {
      stop('class(ptlyr) has to be one of \'matrix\', \'data.frame\'')
    } else {
      prep <- list(matrix = function(x) as_tibble(as.data.frame(x)), 
                   data.frame = function(x) as_tibble(x)
                   )
      ptlyr <- ptlyr %>% 
        prep[[class(ptlyr)[length(class(ptlyr))]]](.) %>% 
        rename(long = 1, lat = 2, layer = 3)
      
      ptlyr <- ptlyr %>% 
        rownames_to_column('ptid')
      ptlyr_trf <- project(cbind(ptlyr$long, ptlyr$lat),
                              proj = as.character(projection)) %>% 
        as_tibble() %>% 
        rename(long = V1, lat = V2) %>% 
        bind_cols(select(ptlyr, ptid, layer))
      ptlyr <- ptlyr_trf
      rm(ptlyr_trf)
    }
  }
  
  if (!is.null(gridlyr) & !is.null(ptlyr)) {
    if (splcol) {
      lims <- list(
        grid = c(min(gridlyr), max(gridlyr)),
        pt = c(min(ptlyr$layer), max(ptlyr$layer))
      )
    } else {
      lims <- list (
        grid = c(min(c(min(ptlyr$layer, na.rm = T), min(gridlyr$layer))), max(c(max(ptlyr$layer, na.rm = T), max(gridlyr$layer))))
      )
    }
  } else if (!is.null(gridlyr)) {
    lims <- list(
      grid = c(min(gridlyr$layer), max(gridlyr$layer))
    )
  } else if (!is.null(ptlyr)) {
    lims <- list(
      grid = c(min(ptlyr$layer), max(ptlyr$layer))
    )
  }
  
  ## plotting
  map_plot <- ggplot()
  
  # bathymetry
  load_bathy <- FALSE
  if (bathymetry | filledbg) {
    if (exists('bathy_data', envir = .STACYmap)) {
      if (!is.null(.STACYmap$bathy_data[[as.character(projection)]])) {
        bathy_data <- .STACYmap$bathy_data[[as.character(projection)]]#lapply(as.vector(names(.STACYmap$bathy_data[[projection]])), function(d, x){as_tibble(x[[d]]) %>% mutate(depth = d)},
                      #       x = .STACYmap$bathy_data[[projection]]) %>% setNames(., sort(c(200, seq(0, 10000, 1000))))
      } else {
        load_bathy <- TRUE
      }
    } else {
      load_bathy <- TRUE
    }
  }
  
  if (load_bathy) {
     bathy_data <- lapply(sort(c(200, seq(0, 10000, 1000))), function (x, L) {
      load_natural_earth_data(file = paste0('ne_10m_bathymetry_', L[[as.character(x)]], '_', x, '.shp')) %>% 
        spTransform(., CRSobj = CRS(projection)) %>% 
        fortify(.)}, L = LETTERS[seq(1, 12, 1) %>% rev()] %>% 
        setNames(., c(200, seq(0, 10000, 1000)) %>% sort())) %>% 
      setNames(., sort(c(200, seq(0, 10000, 1000))))
     .STACYmap$bathy_data[[as.character(projection)]] <- bathy_data
  }
  
  if (filledbg) {
    map_data <- load_natural_earth_data(file = 'ne_10m_land.shp') %>% 
      spTransform(., CRSobj = CRS(projection)) %>% 
      fortify(.)
    
    if (bathymetry) {
      deps <- c(0, 200, seq(1e3, 1e4, 1e3)) %>% as.character()
      cols <- GLOBAL_BATHYMETRY_COLORS
      for (i in 1:12) {
        map_plot <- map_plot + 
          geom_polygon(data = bathy_data[[deps[i]]],
                       mapping = aes(x = long, y = lat, group = group),
                       fill = cols[i], 
                       na.rm = T)
      }
      map_plot <- map_plot + 
        geom_polygon(data = map_data, mapping = aes(x = long, y = lat, group = group, fill = hole)) + 
        scale_fill_manual(values = c(GLOBAL_LAND_COLOR, NA), guide = FALSE) #+ 
        #new_scale_fill() + new_scale_color()
    } else {
      # TODO fix Caspian sea hole
      map_plot <- map_plot + 
        geom_polygon(data = bathy_data[['0']], mapping = aes(x = long, y = lat, group = group,  fill = TRUE, color = TRUE)) + 
        geom_polygon(data = map_data, mapping = aes(x = long, y = lat, group = group, fill = hole, color = hole)) + 
        scale_fill_manual(values = c(GLOBAL_LAND_COLOR, GLOBAL_OCEAN_COLOR), guide = FALSE) + 
        scale_color_manual(values = c(GLOBAL_LAND_COLOR, GLOBAL_OCEAN_COLOR), guide = FALSE) #+ 
        #new_scale_fill() + new_scale_color()
    }
    
    if (!is.null(gridlyr) | !is.null(ptlyr)) {
      map_plot <- map_plot + 
        new_scale_fill() + new_scale_color()
    }
  }
  #return(map_plot)
  # gridlyr
  # TODO: option fopr multiple gridlyrs
  if (!is.null(gridlyr)) {
    map_plot <- map_plot + 
      geom_sf(data = gridlyr,
              mapping = aes(color = layer, fill = layer)) + 
      scale_fill_gradientn(colors = {if (splcol) {colorscheme$grid} else {colorscheme}},
                           limits = lims$grid,
                           guide = guide_legend(legend_names$grid,
                                                direction = 'horizontal'),
                           breaks = round(seq(lims$grid[1], lims$grid[2], length.out = 5), digits = 1)) +
      scale_color_gradientn(colors = {if (splcol) {colorscheme$grid} else {colorscheme}}, 
                            limits = lims$grid, guide = FALSE) 
    
    if (splcol) {
      map_plot <- map_plot + 
        new_scale_color() + 
        new_scale_fill()
    }
  }
  
  # coastline contour
  if (coastline) {
    load_map <- FALSE
    if (exists('map_data', envir = .STACYmap)) {
      if (!is.null(.STACYmap$map_data[[as.character(projection)]])) {
        map_data <- .STACYmap$map_data[[as.character(projection)]]
      } else {
        load_map <- TRUE
      }
    } else {
      load_map <- TRUE
    }
    
    if (load_map) {
      map_data <- load_natural_earth_data(file = 'ne_10m_land.shp') %>% 
        spTransform(., CRSobj = CRS(projection)) %>% 
        fortify(.)
      .STACYmap$map_data[[as.character(projection)]] <- map_data
    }
    
    map_data$id <- as.numeric(map_data$id)
    map_data[map_data$id <= 0, 'id'] <- NA
    
    map_plot <- map_plot + 
      geom_polygon(data = map_data,
                   mapping = aes(x = long, y = lat, group = group),
                   fill = NA,
                   colour = 'black',
                   size = 0.25)
  }
  
  # pointlyr
  if (!is.null(ptlyr)) {
  map_plot <- map_plot + 
    geom_point(data = ptlyr,
               mapping = aes(x = long, y = lat, fill = layer),#, size = has_slope), 
               alpha = I(0.7), shape = 21, size = GLOBAL_POINT_SIZE,
               show.legend = c(fill = {if (is.null(gridlyr) | splcol) {TRUE} else {FALSE}}, size = TRUE))
  
    if (splcol | is.null(gridlyr)) {
      map_plot <- map_plot + 
        scale_fill_gradientn(colors = {if (splcol) {colorscheme$pt} else {colorscheme}}, 
                             limits = lims$point, 
                             guide = guide_legend(legend_names$pt,
                                                  direction = 'horizontal')) 
    }
  }
    
  map_plot <- map_plot + 
    title_and_axis() + 
    theme(panel.ontop = T,
          panel.background = element_blank(),
          axis.title = element_blank(), 
          legend.position = 'bottom',
          panel.grid = element_line(colour = "black", size = 0.1), 
          plot.background = element_blank(), 
          axis.ticks = element_blank()) # element_line(linetype = 'af', colour = 'grey50'))
  
  if (!graticules) {
    map_plot <- map_plot + 
      coord_sf(label_graticule = 'SW', label_axes = '--EN',
               crs = CRS(projection),
               datum = NA,
               xlim = {if (!is.null(xlim)) {xlim}},
               ylim = {if (!is.null(ylim)) {ylim}},#c(-0.63*1e7, 0.77*1e7),
               expand = F) + 
      theme(axis.text = element_blank(),
            panel.grid = element_blank())
  } else {
    map_plot <- map_plot  + 
      coord_sf(label_graticule = 'SW', label_axes = '--EN',
               crs = CRS(projection),
               xlim = xlim,#{if (!is.null(xlim)) {xlim}},
               ylim = ylim,#{if (!is.null(ylim)) {ylim}},#c(-0.63*1e7, 0.77*1e7),
               expand = F)
  }
  
  map_plot <- map_plot + 
    scale_x_continuous() +#seq(-90, 90, 90)) + 
    scale_y_continuous() #seq(-60, 60, 30))
  
  if (!box) {
    map_plot <- map + 
      theme(panel.border = element_blank())
  }
  
  if (legend_inside) {map_plot <- map_plot + 
    theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0), legend.box = 'vertical')
  }
  
  return(map_plot)
}
