library(dplyr)
library(tidyverse)
library(rgdal)

projection_gridlyr <- function(gridlyr, projection){
  
  rotate = FALSE
  transpose = TRUE
  flipy = FALSE

  if (!is.null(gridlyr)) {
    if (!any(class(gridlyr) %in% c('raster', 'matrix', 'data.frame'))) {
      stop('class(gridlyr) has to be one of \'raster\', \'matrix\', \'data.frame\'')
    } else {
      prepg <- list(raster = function(x) x, 
                    matrix = function(x) {
                      if (!rotate) {
                        o <- raster::raster({if (transpose) {t(x)} else {x}}, 
                                            crs = "+proj=longlat +datum=WGS84", xmn = -180, xmx = 180, ymn = -90, ymx = 90) %>% 
                                            {if (flipy) {raster::flip(., 'y')} else {.}}
                      } else {
                        o <- raster::raster({if (transpose) {t(x)} else {x}}, 
                                            crs = "+proj=longlat +datum=WGS84", xmn = 0, xmx = 360, ymn = -90, ymx = 90) %>% 
                          raster::rotate() %>% 
                          {if (flipy) {raster::flip(., 'y')} else {.}}
                      }
                      o
                    },
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
  
  return(gridlyr)
}