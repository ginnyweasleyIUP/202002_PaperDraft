library(dplyr)
library(tidyverse)
library(rgdal)

projection_ptlyr <- function(ptlyr, projection){

  if (!is.null(ptlyr)) {
    if (!any(class(ptlyr) %in% c('matrix', 'data.frame'))) {
      stop('class(ptlyr) has to be one of \'matrix\', \'data.frame\'')
    } else {
      prep <- list(matrix = function(x) as_tibble(as.data.frame(x)), 
                   data.frame = function(x) as_tibble(x)
      )
      
      colnames(ptlyr) <- c("long", "lat", "layer")
      ptlyr <- as.tibble(as.data.frame(ptlyr))
      
      ptlyr <- ptlyr %>% 
        rownames_to_column('ptid')
      ptlyr_trf <- project(cbind(ptlyr$long, ptlyr$lat),
                           proj = as.character(projection)) %>% as_tibble() 
      
      colnames(ptlyr_trf) <- c("long", "lat")
      
      ptlyr <- data.frame(
        long = ptlyr_trf$long,
        lat = ptlyr_trf$lat,
        layer = ptlyr$layer
      )
      rm(ptlyr_trf)
    }
  }
  
  return(ptlyr)
}