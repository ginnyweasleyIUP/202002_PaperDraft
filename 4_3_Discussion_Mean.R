#################################################
## DISCUSSION MEAN ##### ########################
#################################################

library(plyr)
library(dplyr)
library(tidyverse)
library(zoo)
library(PaleoSpec)
library(nest)
library(latex2exp)

DISCUSSION$MEAN <- list()

#################################################

# berechne mittleren Messfehler der Simulation 

source("Functions/aw_mean.R")

gmd18O = numeric(diff(DATA_past1000$time))
for(ii in 1:length(gmd18O)){
  gmd18O[ii] <- simpleawmean(DATA_past1000$SIM_yearly_a$ITPC[,,ii])
}
#mean(value$a)*500000/(4419*mean(value$a)-500000) Temp Änderung pro Delta Änderung

gmd18O[c(582,619)] = NA
Timeseries$HadCM3_GMST_a <- zoo(x = value$a-mean(value$a[1080:1110]), order.by = seq(from = -1*(1950- DATA_past1000$time[1]), to = -1*(1950 - DATA_past1000$time[2]), by = 1))

DISCUSSION$MEAN$std.dev <- sd(gmd18O, na.rm = T)

#################################################
## Regression var ratio and offset

reg <- lm(log(DISCUSSION$scatter_data_variance$vr_ds_itpc) ~ DISCUSSION$scatter_data_variance$diff_down)
confint(reg, 'DISCUSSION$scatter_data_variance$diff_down')


#################################################

# Mean resolution of records:

mean_resolution <- c()

for(entity in DATA_past1000$CAVES$entity_info$entity_id[mask_spec]){
  mean_resolution = c(mean_resolution, mean(diff(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age), na.rm = T))
}

DISCUSSION$MEAN$mean_resolution <- mean_resolution

rm(mean_resolution)

mean(DISCUSSION$MEAN$mean_resolution, na.rm = T)
median(DISCUSSION$MEAN$mean_resolution, na.rm = T)
quantile(DISCUSSION$MEAN$mean_resolution, prob = seq(0,1,0.05))[c(2,20)]
