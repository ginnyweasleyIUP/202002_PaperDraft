#################################################
## PLOTTING SKRIPT ##############################
#################################################

# Das Paper soll aus 6-8 Plots bestehen. Diese werden hier geplottet!

# 1 Einführungsgraphik:   mit Daten aus solar Forcing, Modeled, Bunker Höhle, und GNIP
# 2 SISAL Database:       past millenium caves 
# 3 Mean:                 Map mit Balken über latitudes daneben
#                         APPENDIX: prec-weighted, außerdem mean evolution, scatterplots
# 4 Variance:             Map mit 4 Histogramen
# 5 Correlation:          negativ Plot mit Temp oder Prec in Rot oder blau
#                         APPENDIX: regression
# 6 Network:              4er Plot (2 Karten, 2 distance vs correlation plot)
# 7 Spektrum:             be creative
# 8 LIA and MCA:          4er Karte (warmest, coldest, dryest, wettest)

library(plyr)
library(dplyr)

PLOTTING_VARIABLES$WIDTH = 8.3
PLOTTING_VARIABLES$HEIGHT = 5.5

#################################################

## ToDo:

# [ ] einheitliche Schriftgröße zu allem!
# [ ] Kasten um Legende entfernen
# [ ] Area weighing bei den Spektren
# [ ] Ohne Land-Sea mask LIC MCA
# [ ] logarithmische Skala bei Isot.Plot
# [ ] Scatter Plots ausbauen
# [ ] Variance Plot über zwei Spalten ziehen
# [ ] EINLEITUNGSPLOT
# [ ] SISAL Plot legende in Karte rein
# [ ] Variance nicht über log(v1/v2) ausdrücken sondern andere Skala (immer noch logarithmisch, aber nur die Achse)
# [ ] Network Plots coastline nicht grau
# [ ] LIA und MCA Zeiten raussuchen...


#################################################
## 1 EINFÜHRUNGSPLOT ############################
#################################################

#################################################
## 2 SISAL DATABASE #############################
#################################################


source("Functions/Plotting/karst_map_plot.R")
library(tidyverse)

#aus 1_5_create_dataset_SISAL.R --> site_past1000_min1



ALL_SITES <- read.csv("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/02_SISAL/SISAL_v2_CARLA/site_countries.csv")

sites_spec <- DATA_past1000$CAVES$entity_info$site_id[mask_spec]
sites_var  <- DATA_past1000$CAVES$entity_info$site_id[mask_var]
sites_mean <- DATA_past1000$CAVES$entity_info$site_id[mask_mean]

USED_SITES_spec <- ALL_SITES %>% filter(site_id %in% DATA_past1000$CAVES$site_info$site_id) %>% filter(site_id %in% sites_spec) %>% distinct(site_id, longitude, latitude)
USED_SITES_spec <- data.frame(
  lon = USED_SITES_spec$longitude,
  lat = USED_SITES_spec$latitude,
  value = USED_SITES_spec$site_id
)
USED_SITES_var <- ALL_SITES %>% filter(site_id %in% DATA_past1000$CAVES$site_info$site_id) %>% filter(!site_id %in% sites_spec) %>% 
  filter(site_id %in% sites_var) %>% distinct(site_id, longitude, latitude)
USED_SITES_var <- data.frame(
  lon = USED_SITES_var$longitude,
  lat = USED_SITES_var$latitude,
  value = USED_SITES_var$site_id
)
USED_SITES_mean <- ALL_SITES %>% filter(site_id %in% DATA_past1000$CAVES$site_info$site_id) %>% filter(!site_id %in% sites_spec) %>% 
  filter(!site_id %in% sites_var) %>% 
  filter(site_id %in% sites_mean) %>% distinct(site_id, longitude, latitude)
USED_SITES_mean <- data.frame(
  lon = USED_SITES_mean$longitude,
  lat = USED_SITES_mean$latitude,
  value = USED_SITES_mean$site_id
)
NOT_SITES <- ALL_SITES %>% filter(!site_id %in% sites_mean) %>% distinct(site_id, longitude, latitude)
NOT_SITES <- data.frame(
  lon = NOT_SITES$longitude,
  lat = NOT_SITES$latitude,
  value = NOT_SITES$site_id
)

GLOBAL_STACY_OPTIONS$GLOBAL_FONT_SIZE = 10
plot <- karst_map_plot(USED_SITES_spec = USED_SITES_spec,
                       USED_SITES_var = USED_SITES_var,
                       USED_SITES_mean = USED_SITES_mean,
                       NOT_SITES = NOT_SITES, pt_size = 3) + 
  #theme(legend.position = "bottom", legend.direction = "horizontal")#c(0.0, 0.02))
  theme(legend.position = c(-0.01, 0), legend.justification = c(0, 0), legend.box = 'vertical',
        axis.text = element_blank())
#plot
plot <- plot + theme(panel.border = element_blank()) #,legend.text = element_text(size = 8)

#plot
plot %>% ggsave(filename = paste('Paper_Plot_2_SISAL_database', 'pdf', sep = '.'), plot = ., path = 'Plots/Paper', 
                width = 2*PLOTTING_VARIABLES$WIDTH, height = 2*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")

remove(plot, NOT_SITES, USED_SITES_mean, USED_SITES_spec, USED_SITES_var, ALL_SITES, sites_mean, sites_spec, sites_var)

#################################################
## 3 MEAN #######################################
#################################################

library(plyr)
library(dplyr)
library(rgdal)
source("Functions/STACYmap_5.R")
source("Functions/STACYmap_5_1_NAgrid.R")
source("Functions/STACYmap_5_2_logscale.R")

GLOBAL_STACY_OPTIONS$GLOBAL_FONT_SIZE = 10

temp <- simpleawmean(DATA_past1000$SIM_mean$TEMP, seq(from = -90, to = 90, length.out = 73))
temp_sd <- simpleawsd(DATA_past1000$SIM_mean$TEMP, seq(from = -90, to = 90, length.out = 73))
prec <- simpleawmean(DATA_past1000$SIM_mean$PREC, seq(from = -90, to = 90, length.out = 73))*8.6148e4
prec_sd <- simpleawsd(DATA_past1000$SIM_mean$PREC, seq(from = -90, to = 90, length.out = 73))*8.6148e4

isot <- simpleawmean(DATA_past1000$SIM_mean$ISOT, seq(from = -90, to = 90, length.out = 73))
isot_sd <- simpleawsd(DATA_past1000$SIM_mean$ISOT, seq(from = -90, to = 90, length.out = 73))
slpr <- simpleawmean(DATA_past1000$SIM_mean$SLPR, seq(from = -90, to = 90, length.out = 73))
slpr_sd <- simpleawsd(DATA_past1000$SIM_mean$SLPR, seq(from = -90, to = 90, length.out = 73))


plot_temp <- STACYmap(gridlyr = rbind(DATA_past1000$SIM_mean$TEMP[49:96,1:73],DATA_past1000$SIM_mean$TEMP[1:48,1:73]),
                      #zoom = c(-180, -60, 180, 73),
                      legend_names = list(grid = "Mean Temperature (°C)"),
                      graticules = TRUE,
                      colorscheme = "temp", 
                      centercolor = 0) +
  theme(panel.border = element_blank(),
        legend.background = element_blank(),
        axis.text = element_blank(),
        legend.text = element_text(size = 8)) 

plot_prec <- STACYmap(gridlyr = 8.6148e4*rbind(DATA_past1000$SIM_mean$PREC[49:96,1:73],DATA_past1000$SIM_mean$PREC[1:48,1:73]),
                      #zoom = c(-180, -60, 180, 73),
                      legend_names = list(grid = "Mean Precipitation (mm/day)"),
                      graticules = TRUE,
                      colorscheme = "prcp_grd") +
  theme(panel.border = element_blank(),
        legend.background = element_blank(),
        axis.text = element_blank(),
        legend.text = element_text(size = 8)) 

color_slpr = c("#edf7fb", "#ddedf4", "#cfe2ef", "#c1d6e8", "#b3cde2", "#a8bfdb", "#9eb1d4", "#95a3cd",
               "#8c95c6", "#8a85bd", "#8975b5", "#8966ae", "#8756a7", "#85449c", "#811e84", "#800e7c")

PLOTTING_VARIABLES$COLORS$SLPR <- color_slpr
remove(color_slpr)

plot_slpr <- STACYmap(gridlyr = rbind(DATA_past1000$SIM_mean$SLPR[49:96,1:73],DATA_past1000$SIM_mean$SLPR[1:48,1:73]),
                      #zoom = c(-180, -60, 180, 73),
                      legend_names = list(grid = "Mean sea level pressure (mbar)"),
                      graticules = TRUE,
                      colorscheme = PLOTTING_VARIABLES$COLORS$SLPR) +
  theme(panel.border = element_blank(),
        legend.background = element_blank(),
        axis.text = element_blank(),
        legend.text = element_text(size = 8)) 

Plot_lyr1 <-rbind(DATA_past1000$SIM_mean$ISOT[49:96,1:73],DATA_past1000$SIM_mean$ISOT[1:48,1:73])

Plot_lyr3 <- Plot_lyr1
Plot_lyr3[is.na(Plot_lyr3)] = 1000
Plot_lyr3[Plot_lyr3>0] <- Plot_lyr3[Plot_lyr3>0]+1 
Plot_lyr3[Plot_lyr3<0] <- Plot_lyr3[Plot_lyr3<0]-1
Plot_lyr3[Plot_lyr3>0] <- log(Plot_lyr3[Plot_lyr3>0])
Plot_lyr3[Plot_lyr3<0] <- - log(abs(Plot_lyr3[Plot_lyr3<0]))
Plot_lyr3[abs(Plot_lyr3)>5] <- NA
#Plot_lyr3[,1:6] <- NA
#Plot_lyr3[,60:73] <- NA
#Plot_lyr1[,1:6] <- NA
#Plot_lyr1[,60:73] <- NA


Point_lyr <- data.frame(
  lon = MEAN_ANALYSIS$POINTS$CAVElyr_isot_used$lon,
  lat = MEAN_ANALYSIS$POINTS$CAVElyr_isot_used$lat,
  value = - log(abs(MEAN_ANALYSIS$POINTS$CAVElyr_isot_used$value -1))
)

Point_lyr$value[[57]] <- log(MEAN_ANALYSIS$POINTS$CAVElyr_isot_used$value[[57]]+1)

# Plot_lyr1[Plot_lyr1 < -17] = -17 -0.1*Plot_lyr1
# Plot_lyr2 <-rbind(DATA_past1000$SIM_mean$ISOT[49:96,1:73],DATA_past1000$SIM_mean$ISOT[1:48,1:73])
# Plot_lyr2[Plot_lyr2 > -17] = NA
# Plot_lyr2[Plot_lyr2 < -17] = -17

allmax = - min(Plot_lyr3, na.rm = T)
allmax_real = - min(Plot_lyr1, na.rm = T)

GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE = 4

plot_isot <- STACYmap_isot(gridlyr = Plot_lyr3,
                      ptlyr = Point_lyr,
                      #zoom = c(-180, -60, 180, 73),
                      legend_names = list(grid = "Mean d18O (%)"),
                      graticules = TRUE,
                      colorscheme = RColorBrewer::brewer.pal(9, 'BrBG'),
                      centercolor = 0, 
                      breaks_isot = c(-allmax,-log(11), -log(2), 0, log(2), log(11), allmax),
                      labels_isot = c(round(-allmax_real), -10, -1, 0, 1, 10, round(allmax_real))) +
    theme(panel.border = element_blank(),
          legend.background = element_blank(),
          axis.text = element_blank(),
          text = element_text(size = 8)) 

plot_isot

#plot_isot %>% ggsave(filename = "Nadine_hadcm3_pmil_d18Oinprecip.pdf", path = "Plots", width = PLOTTING_VARIABLES$WIDTH, height = PLOTTING_VARIABLES$HEIGHT, units = "cm", dpi = 'print')

library(ggpubr)
plot <- ggarrange(plot_temp, plot_prec, plot_slpr, plot_isot,
                  labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2)

plot  %>% ggsave(filename = paste('Paper_Plot_3_Mean', 'pdf', sep = '.'), plot = ., path = 'Plots/Paper', 
                 width = 2*12, height = 2*12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")

remove(plot_temp, plot_prec, plot_isot, plot_slpr, Plot_lyr1, Plot_lyr2, plot)

#################################################
## 4 VARIANCE ###################################
#################################################

library(latex2exp)
source("Functions/STACYmap_6.R")

## Total of 5 Plots that then have to be arranged
## 1) Variance Map
## 2) ISOT density
## 3) ITPC density
## 4) TEMP density
## 5) PREC density

## Variance Map #################################

source("Functions/Plotting/var_map_plot.R")

# mask_var_2 <- mask_var
# mask_var_2[85] <- FALSE
# mask_var_2[63] <- FALSE
# mask_var_2[64] <- FALSE
# mask_var_2[104] <- FALSE
# mask_var_2[38] <- FALSE


Point_Lyr <- data.frame(
  lon = VAR_ANALYSIS$POINTS$CAVElyr$lon[mask_var],
  lat = VAR_ANALYSIS$POINTS$CAVElyr$lat[mask_var],
  value = log(VAR_ANALYSIS$POINTS$CAVElyr$value_VR_isot[mask_var])
  #value = log(VAR_ANALYSIS$POINTS$CAVElyr$value_VR_isot_ds[mask_var])
)

GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE = 2.5

plot_var <- var_map_plot(Point_Lyr =  Point_Lyr, pt_size =  3, txt_size =  10)
plot_var

plot_var %>% ggsave(filename = paste('Paper_Plot_4_Variance_1_map', 'pdf', sep = '.'), plot = ., path = 'Plots/Paper', 
                width = 2*6, height = 2*PLOTTING_VARIABLES$HEIGHT/1.5, units = 'cm', dpi = 'print', device = "pdf")


#plot_var

## Density Plots ################################

pdf(file = "Plots/Paper/Paper_Plot_4_Variance_2_histo.pdf", width = 2*6, height = 2*PLOTTING_VARIABLES$HEIGHT/1.5)
par(mfrow=c(2,2),oma = c(1,3,0,0) + 0.1,mar = c(3,1,3,1) + 0.1, new = FALSE)

hist(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_temp*mask_mean), 
     breaks = 9, border = "white", prob = TRUE, 
     ylim = c(0,0.25), 
     xlim = c(-6, 6),
     xlab = "",
     xaxt = 'n',
     main = "Temperature", cex.main = 1.5, cex.axis = 1.5)
axis(side = 1, at = c(log(0.01), log(0.1), 0, log(10), log(100)), 
     labels = c(0.01, 0.1, 1, 10, 100), cex.axis = 1.5)
lines(density(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_temp*mask_mean), na.rm = T),
      lwd = 2, col = "black")
lines(density(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_temp_ds*mask_mean), na.rm = T),
      lwd = 2, col = "#B2182B")
abline(v=0, col = "grey60", lty = 3)
#mtext(text = TeX("$$Var_{Rec}/Var_{Sim}$$"),side = 1,line = 2.5, cex = 1.5)
mtext(text = "density",side = 2,line = 2.5, cex = 1.5)
text(2.0, 0.18, "down-sampled", col = "#B2182B", cex = 1.5)
text(-4.5, 0.15, "full", col = "black", cex = 1.5)

hist(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_prec*mask_mean), 
     breaks = 9, border = "white", prob = TRUE, 
     ylim = c(0,0.25),
     xlim = c(-6,6),
     xlab = "",
     xaxt = 'n',
     main = "Precipitation", cex.main = 1.5, cex.axis = 1.5)
axis(side = 1, at = c(log(0.01), log(0.1), 0, log(10), log(100)), 
     labels = c(0.01, 0.1, 1, 10, 100), cex.axis = 1.5)
lines(density(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_prec*mask_mean), na.rm = T),
      lwd = 2, col = "black")
lines(density(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_prec_ds*mask_mean), na.rm = T),
      lwd = 2, col = "#B2182B")
abline(v=0, col = "grey60", lty = 3)
#mtext(text = TeX("$$Var_{Rec}/Var_{Sim}$$"),side = 1,line = 2.5, cex = 1.5)

hist(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_isot*mask_mean), 
     breaks = 9, border = "white", prob = TRUE, 
     ylim = c(0,0.25),
     xlim = c(-6,6),
     xaxt = 'n',
     xlab = "",
     main = "d18O in precipitation", cex.main = 1.5, cex.axis = 1.5)
axis(side = 1, at = c(log(0.01), log(0.1), 0, log(10), log(100)), 
     labels = c(0.01, 0.1, 1, 10, 100), cex.axis = 1.5)
lines(density(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_isot*mask_mean), na.rm = T),
      lwd = 2, col = "black")
lines(density(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_isot_ds*mask_mean), na.rm = T),
      lwd = 2, col = "#B2182B")
abline(v=0, col = "grey60", lty = 3)
mtext(text = TeX("$$Var_{Rec}/Var_{Sim}$$"),side = 1,line = 2.7, cex = 1.2)
mtext(text = "density",side = 2,line = 2.5, cex = 1.5)

hist(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_itpc*mask_mean), 
     breaks = 9, border = "white", prob = TRUE, 
     ylim = c(0,0.25),
     xlim = c(-6,6),
     xlab = "",
     xaxt = 'n',
     main = "precipitation weighted d18O", cex.main = 1.5, cex.axis = 1.5)
axis(side = 1, at = c(log(0.01), log(0.1), 0, log(10), log(100)), 
     labels = c(0.01, 0.1, 1, 10, 100), cex.axis = 1.5)
lines(density(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_itpc*mask_mean), na.rm = T),
      lwd = 2, col = "black")
lines(density(log(VAR_ANALYSIS$VAR_RATIOS$var_ratio_sim_itpc_ds*mask_mean), na.rm = T),
      lwd = 2, col = "#B2182B")
abline(v=0, col = "grey60", lty = 3)
mtext(text = TeX("$$Var_{Rec}/Var_{Sim}$$"),side = 1,line = 2.7, cex = 1.2)

dev.off()




## ITPC density #################################
## TEMP density #################################
## PREC density #################################

#################################################
## 5 CORRELATION ################################
#################################################

source("Functions/projection_ptlyr.R")
# Grid Layer for plotting:
# all areas where d18O correlates better with temperature are marked in red
# all areas where d18O correlates better with precipitation are marked in blue
Plot_lyr_temp <- CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_TEMP_ISOT
Plot_lyr_temp_p <- CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_TEMP_ISOT_P
Plot_lyr_prec <- CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_PREC_ISOT
Plot_lyr_prec_p <- CORR_ANALYSIS$GLOBAL_CORRELATION$CORR_PREC_ISOT_P
Plot_lyr_temp[Plot_lyr_temp_p > 0.1] <- 0
Plot_lyr_temp[abs(Plot_lyr_temp) < 0.2] <- 0
Plot_lyr_prec[Plot_lyr_prec_p > 0.1] <- 0
Plot_lyr_prec[abs(Plot_lyr_prec) < 0.2] <- 0

Plot_lyr_2 <- Plot_lyr_temp
Plot_lyr_3 <- Plot_lyr_prec

Plot_lyr_2[abs(Plot_lyr_prec)>abs(Plot_lyr_temp)] <- 0
Plot_lyr_3[abs(Plot_lyr_temp)>abs(Plot_lyr_prec)] <- 0

Plot_lyr <- abs(Plot_lyr_2)- abs(Plot_lyr_3)
Plot_lyr[Plot_lyr == 0] <- NA

Plot_lyr <- rbind(Plot_lyr[49:96,1:73],
                  Plot_lyr[1:48,1:73])

remove(Plot_lyr_2, Plot_lyr_3, Plot_lyr_prec, Plot_lyr_prec_p, Plot_lyr_temp, Plot_lyr_temp_p)

##### Point Layer

# How should points be colored? Is it so relevant if sign is equal?

# 0) Check for significance --> if not then, then put in Point_lyr_2
# 1) Check for what the absolute corellation is stronger
# 2) make different shapes depending on sign fitting or not


### HERE HERE HERE ############################
## es muss noch angepasst werden, dass alle Punktlisten mit unterschiedlichem Symbol über eine andere Liste gemacht wird. 
Point_Lyr_sign <- list(lon = list(), lat = list(), value = list())
Point_Lyr_notsign <- list(lon = list(), lat = list(), value = list())
Point2_Lyr <- list(lon = list(), lat = list(), value = list())

length_cave = length(DATA_past1000$CAVES$entity_info$site_id)

for(ii in 1:length_cave){
  site <- DATA_past1000$CAVES$entity_info$site_id[ii]
  print(ii)
  if(is.na(CORR_ANALYSIS$CAVE_CORRELATION$PVALUE_TEMP[ii])){next}
  # 1) sortiert aus, was nicht signifikant ist
  if(CORR_ANALYSIS$CAVE_CORRELATION$PVALUE_TEMP[ii] > 0.1 & CORR_ANALYSIS$CAVE_CORRELATION$PVALUE_PREC[ii] > 0.1){
    Point2_Lyr$lon = c(Point2_Lyr$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point2_Lyr$lat = c(Point2_Lyr$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
    Point2_Lyr$value = c(Point2_Lyr$value, CORR_ANALYSIS$CAVE_CORRELATION$CORR_TEMP[ii])
  # 2) betrachte signifikante Korrelationen:
  }else{
    # 2.1) Nur signifikante Korrelation bei Temp
    if(CORR_ANALYSIS$CAVE_CORRELATION$PVALUE_TEMP[ii] < 0.1 & CORR_ANALYSIS$CAVE_CORRELATION$PVALUE_PREC[ii] > 0.1){
      #Check sign to determine shape
      if(sign(CORR_ANALYSIS$CAVE_CORRELATION$CORR_TEMP[ii]) == sign(CORR_ANALYSIS$SITE_CORRELATION$CORR_TI[ii])){
        Point_Lyr_sign$lon = c(Point_Lyr_sign$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
        Point_Lyr_sign$lat = c(Point_Lyr_sign$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
        Point_Lyr_sign$value = c(Point_Lyr_sign$value, abs(CORR_ANALYSIS$CAVE_CORRELATION$CORR_TEMP[ii]))
      }else{
        Point_Lyr_notsign$lon = c(Point_Lyr_notsign$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
        Point_Lyr_notsign$lat = c(Point_Lyr_notsign$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
        Point_Lyr_notsign$value = c(Point_Lyr_notsign$value, abs(CORR_ANALYSIS$CAVE_CORRELATION$CORR_TEMP[ii]))
      }
    }
    
    # 2.2) Nur signifikante Korrelation bei Prec
    else if(CORR_ANALYSIS$CAVE_CORRELATION$PVALUE_TEMP[ii] > 0.1 & CORR_ANALYSIS$CAVE_CORRELATION$PVALUE_PREC[ii] < 0.1){
      if(sign(CORR_ANALYSIS$CAVE_CORRELATION$CORR_PREC[ii]) == sign(CORR_ANALYSIS$SITE_CORRELATION$CORR_PI[ii])){
        Point_Lyr_sign$lon = c(Point_Lyr_sign$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
        Point_Lyr_sign$lat = c(Point_Lyr_sign$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
        Point_Lyr_sign$value = c(Point_Lyr_sign$value, - abs(CORR_ANALYSIS$CAVE_CORRELATION$CORR_PREC[ii]))
      }else{
        Point_Lyr_notsign$lon = c(Point_Lyr_notsign$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
        Point_Lyr_notsign$lat = c(Point_Lyr_notsign$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
        Point_Lyr_notsign$value = c(Point_Lyr_notsign$value, - abs(CORR_ANALYSIS$CAVE_CORRELATION$CORR_PREC[ii]))
        }
    }
    
    # 2.3) Sowohl signifikant für Prec wie für Temp    
    else{
      # 2.3.1) absolute CORR größer für Temp als für Prec
      if(abs(CORR_ANALYSIS$CAVE_CORRELATION$CORR_TEMP[ii]) > abs(CORR_ANALYSIS$CAVE_CORRELATION$CORR_PREC[ii])){
        if(sign(CORR_ANALYSIS$CAVE_CORRELATION$CORR_TEMP[ii]) == sign(CORR_ANALYSIS$SITE_CORRELATION$CORR_TI[ii])){
          Point_Lyr_sign$lon = c(Point_Lyr_sign$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
          Point_Lyr_sign$lat = c(Point_Lyr_sign$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
          Point_Lyr_sign$value = c(Point_Lyr_sign$value, abs(CORR_ANALYSIS$CAVE_CORRELATION$CORR_TEMP[ii]))
        }else{
          Point_Lyr_notsign$lon = c(Point_Lyr_notsign$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
          Point_Lyr_notsign$lat = c(Point_Lyr_notsign$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
          Point_Lyr_notsign$value = c(Point_Lyr_notsign$value, abs(CORR_ANALYSIS$CAVE_CORRELATION$CORR_TEMP[ii]))}
      }
      # 2.3.2) absolute CORR größer für Prec als für Temp
      else{
        if(sign(CORR_ANALYSIS$CAVE_CORRELATION$CORR_PREC[ii]) == sign(CORR_ANALYSIS$SITE_CORRELATION$CORR_PI[ii])){
          Point_Lyr_sign$lon = c(Point_Lyr_sign$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
          Point_Lyr_sign$lat = c(Point_Lyr_sign$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
          Point_Lyr_sign$value = c(Point_Lyr_sign$value, - abs(CORR_ANALYSIS$CAVE_CORRELATION$CORR_PREC[ii]))
        }else{
          Point_Lyr_notsign$lon = c(Point_Lyr_notsign$lon, DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site])
          Point_Lyr_notsign$lat = c(Point_Lyr_notsign$lat, DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site])
          Point_Lyr_notsign$value = c(Point_Lyr_notsign$value, - abs(CORR_ANALYSIS$CAVE_CORRELATION$CORR_PREC[ii]))}
      }
    }
  }
}



Point_Lyr_sign$lon = as.numeric(Point_Lyr_sign$lon)
Point_Lyr_sign$lat = as.numeric(Point_Lyr_sign$lat)
Point_Lyr_sign$value = as.numeric(Point_Lyr_sign$value)

Point_Lyr_notsign$lon = as.numeric(Point_Lyr_notsign$lon)
Point_Lyr_notsign$lat = as.numeric(Point_Lyr_notsign$lat)
Point_Lyr_notsign$value = as.numeric(Point_Lyr_notsign$value)

Point2_Lyr$lon = as.numeric(Point2_Lyr$lon)
Point2_Lyr$lat = as.numeric(Point2_Lyr$lat)
Point2_Lyr$value = as.numeric(Point2_Lyr$value)



GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE <- 3

Point_Lyr_sign_p <-  projection_ptlyr(as.data.frame(Point_Lyr_sign), as.character('+proj=robin +datum=WGS84'))
Point_Lyr_notsign_p <-  projection_ptlyr(as.data.frame(Point_Lyr_notsign), as.character('+proj=robin +datum=WGS84'))
Point2_Lyr_p <-  projection_ptlyr(as.data.frame(Point2_Lyr), as.character('+proj=robin +datum=WGS84'))

remove(Point_Lyr_sign, Point_Lyr_notsign, Point2_Lyr)

# Jetzt existiert ein Plot Layer und 2 Point Layer die man nur noch plotten muss und eine richtige Legende dafür braucht...

source("Functions/STACYmap_6.R")
source("Functions/STACYmap_5_2_logscale_corr.R")

plot <- STACYmap_isot_corr(gridlyr = Plot_lyr, centercolor = 0, graticules = T, 
                      legend_names = list(grid = "abs(Corr.)"),
                      breaks_isot = c(-1, -0.5, 0, 0.51, 1),
                      labels_isot = c(1, "corr prec", "0",  "corr temp", 1)) + 
  geom_point(data = Point2_Lyr_p, aes(x = long, y = lat, shape = "1"), fill = 'gray', size = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE-1, show.legend = c(shape = T)) + 
  geom_point(data = Point_Lyr_sign_p, aes(x = long, y = lat, fill = layer, shape = "2"), size = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE, show.legend = c(color = T, shape = T)) + 
  geom_point(data = Point_Lyr_notsign_p, aes(x = long, y = lat, fill = layer, shape = "3"), size = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE, show.legend = c(color = T, shape = T)) + 
  scale_shape_manual(name = NULL, labels = c("no corr.-sites", "same sign", "different sign"), 
                     values = c(20,21,23))+
  #guides(fill = guide_colorbar(label = F, direction = "horizontal", title = "|Corr.| blue prec, red temp")) + 
  theme(panel.border = element_blank(),
        legend.background = element_blank(),
        axis.text = element_blank(),
        text = element_text(size = 12),
        legend.title = element_text(size = 12))

plot

plot %>% ggsave(filename = paste('Paper_Plot_5_Correlation', 'pdf', sep = '.'), plot = ., path = 'Plots/Paper', 
                width = 2*PLOTTING_VARIABLES$WIDTH, height = 2*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")


#################################################
## 6 NETWORK ####################################
#################################################

source("Functions/networkmap_simple3.R")

# First Points need to be orderd (lowess doesn't work otherwise)

C_SIM_p <- NETWORK_ANALYSIS$CORR_MATRIX_SIM$r
C_SIM_p[NETWORK_ANALYSIS$CORR_MATRIX_SIM$P>0.1] = NA

C_REC_p <- NETWORK_ANALYSIS$CORR_MATRIX_REC$r
C_REC_p[NETWORK_ANALYSIS$CORR_MATRIX_REC$P>0.1] = NA

plot_dist <- NETWORK_ANALYSIS$DIST
plot_dist[lower.tri(NETWORK_ANALYSIS$DIST)] = NA
lowess_dist <- as.vector(NETWORK_ANALYSIS$DIST[upper.tri(NETWORK_ANALYSIS$DIST)])
o <- order(lowess_dist)
lowess_dist_sorted <- lowess_dist[o]

plot_c_sim <- NETWORK_ANALYSIS$CORR_MATRIX_SIM$r
plot_c_sim[lower.tri(NETWORK_ANALYSIS$CORR_MATRIX_SIM$r, diag = FALSE)] = NA
lowess_c_sim <- as.vector(NETWORK_ANALYSIS$CORR_MATRIX_SIM$r[upper.tri(NETWORK_ANALYSIS$CORR_MATRIX_SIM$r)])
lowess_c_sim_sorted <- lowess_c_sim[o]
plot_c_rec <- NETWORK_ANALYSIS$CORR_MATRIX_REC$r
plot_c_rec[lower.tri(NETWORK_ANALYSIS$CORR_MATRIX_REC$r)] = NA
lowess_c_rec <- as.vector(NETWORK_ANALYSIS$CORR_MATRIX_REC$r[upper.tri(NETWORK_ANALYSIS$CORR_MATRIX_REC$r)])
lowess_c_rec_sorted <- lowess_c_rec[o]

lo <- loess(lowess_c_rec_sorted ~ lowess_dist_sorted)

# dist_vec_sisal <- as.vector(dist_matrix_sisal[upper.tri(dist_matrix_sisal)])
# o <- order(dist_vec_sisal)
# dist_vec_sisal_sorted <- dist_vec_sisal[o]
# corr_vec_sisal <- as.vector(C_SISAL[upper.tri(C_SISAL)])
# corr_vec_sisal_sorted <- corr_vec_sisal[o]

# x <- 1:10
# y <- c(2,4,6,8,7,12,14,16,18,20)
# lo <- loess(y~x)
# plot(x,y)
# lines(predict(lo), col='red', lwd=2)

# #fitting smoothing splines using smooth.spline(X,Y,df=...)
# fit1<-smooth.spline(age,wage,df=16) #16 degrees of freedom
# #Plotting both cubic and Smoothing Splines 
# plot(age,wage,col="grey",xlab="Age",ylab="Wages")
# points(age.grid,predict(fit,newdata = list(age=age.grid)),col="darkgreen",lwd=2,type="l")
# #adding cutpoints
# abline(v=c(25,40,60),lty=2,col="darkgreen")
# lines(fit1,col="red",lwd=2)
# legend("topright",c("Smoothing Spline with 16 df","Cubic Spline"),col=c("red","darkgreen"),lwd=2)

scaling = 1.5
spacing = 0.7

pdf(file = "Plots/Paper/Paper_Plot_6_Network.pdf", height= PLOTTING_VARIABLES$HEIGHT/2, width = PLOTTING_VARIABLES$WIDTH/2)
par(mfrow=c(2,2), mai = c(rep(spacing, 4)), mar = c(3,3,2,0.5))
#SIM MAP
networkmap_simple3(CMAT = C_SIM_p, 
                   lat = NETWORK_ANALYSIS$lats, 
                   lon = NETWORK_ANALYSIS$longs, 
                   title = "Corr HadCM3 past millenium, p<0.1", 
                   thresh = 0.15)
#SIM Cor-Dist
plot(plot_dist, plot_c_sim, 
     ylim = c(-1,1),
     xlim = c(0,20000),
     ylab = "",
     xlab = "", 
     cex = 1, 
     lwd = 0.5, 
     panel.first = grid())
lo <- loess(lowess_c_sim_sorted ~ lowess_dist_sorted, span = 0.2)
lines(lo$x, lo$fitted, lwd = 4, col = "#B2182B")
#lines(lowess(lowess_dist_sorted,lowess_c_sim_sorted, f=0.1), lwd = 4, col = "#B2182B")
mtext("Distance between pairs (km)", side= 1, line = 2)

#SISAL MAP
networkmap_simple3(CMAT = C_REC_p, 
                   lat = NETWORK_ANALYSIS$lats, 
                   lon = NETWORK_ANALYSIS$longs,
                   title = "Corr SISAL past millenium, p<0.1",
                   thresh = 0.4)
#SISAL Cor-Dist
plot(plot_dist, plot_c_rec, 
     ylim = c(-1,1),
     xlim = c(0,20000),
     ylab = "",
     xlab = "", 
     cex = 1, 
     lwd = 1, 
     panel.first = grid())
lo <- loess(lowess_c_rec_sorted ~ lowess_dist_sorted, span = 0.2)
lines(lo$x, lo$fitted, lwd = 4, col = "#B2182B")
mtext("Distance between pairs (km)", side= 1, line = 2)
dev.off()

#################################################
## 7 SPECTRA ####################################
#################################################


pdf(file = "Plots/Paper/Paper_Plot_7_Spectra.pdf", width = PLOTTING_VARIABLES$WIDTH/2, height = PLOTTING_VARIABLES$HEIGHT/2*1.2)
LPlot(SPEC_ANALYSIS$MEAN_SPEC$SIM_full$spec, col = "#074893", 
      ylim = c(0.00001,1000), xlim = c(1/500, 0.5),
      xaxt = 'n',
      xlab = "")#,
      #main = TeX("Mean Spectra from cave locations (res>8)"))
mtext("Periode (years)", side = 1, line= 2)
LLines(SPEC_ANALYSIS$MEAN_SPEC$SIM_full_down$spec, col = "#074893", lty = 3 )
LLines(SPEC_ANALYSIS$MEAN_SPEC$SIM_full_rec$spec, col = "#074893", lty = 3)
#text(0.2, 8e2, "HadCM3 yearly res.", col = "#074893")
#text(0.3, 3e2, "5y filter", col = "#074893")
#text(0.3, 1e2, "50y filter", col = "#074893")

axis(side = 1, at = c(0.002, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5), 
     labels = c(1/0.001, 1/0.005, 1/0.01, 1/0.02, 1/0.05, 1/0.2, 1/0.5))

LLines(SPEC_ANALYSIS$MEAN_SPEC$SIM_ds$spec, col = "#91002B")
LLines(SPEC_ANALYSIS$MEAN_SPEC$SIM_down_rec$spec, col = "#91002B", lty = 3)
#text(0.02, 0.01, "HadCM3 down-sampled", col = "#91002B")
#text(0.02, 0.001, "10y filter", col = "#91002B")

LLines(SPEC_ANALYSIS$MEAN_SPEC$Record$spec, col = "black", lw = 2)
#text(0.02, 0.0005, "Records", col = "black")



legend("bottomleft", legend = c("HadCM3 yearly res.", " ... 5y filter", " ... 50y filter", "HadCM3 down-sampled to record res.", " ... 10y filter", "Records"), 
       col = c("#074893","#074893","#074893","#91002B","#91002B","black"), lwd = c(1,1,1,1,1,2), lty = c(1,3,3,1,3,1), bty = "n")
dev.off()



#################################################
## 8 LIA and MCA ################################
#################################################

library(ggpubr)

Plot_lyr <- rbind(LIA_MCA$SIM_full$SIM_low_101_T[49:96,1:73]*PLOTTING_VARIABLES$ls_mask[49:96,1:73],
                  LIA_MCA$SIM_full$SIM_low_101_T[1:48,1:73]*PLOTTING_VARIABLES$ls_mask[1:48,1:73])

plot_low_T <- STACYmap(gridlyr = Plot_lyr, colorscheme = RColorBrewer::brewer.pal(9, "Reds"), zoom = c(-180,-60,180,73), legend_names = list(grid = "years BP")) + 
  ggtitle("Coldest 101y average T") + theme(plot.title = element_text(h = 0.5), panel.border = element_blank(),
                                            legend.background = element_blank(),
                                            axis.text = element_blank())

Plot_lyr <- rbind(LIA_MCA$SIM_full$SIM_high_101_T[49:96,1:73]*PLOTTING_VARIABLES$ls_mask[49:96,1:73],
                  LIA_MCA$SIM_full$SIM_high_101_T[1:48,1:73]*PLOTTING_VARIABLES$ls_mask[1:48,1:73])

plot_high_T <- STACYmap(gridlyr = Plot_lyr, colorscheme = RColorBrewer::brewer.pal(9, "Reds"), zoom = c(-180,-60,180,73), legend_names = list(grid = "years BP")) + 
  ggtitle("Warmest 101y average T") + theme(plot.title = element_text(h = 0.5), panel.border = element_blank(),
                                            legend.background = element_blank(),
                                            axis.text = element_blank())

Plot_lyr <- rbind(LIA_MCA$SIM_full$SIM_low_101_P[49:96,1:73]*PLOTTING_VARIABLES$ls_mask[49:96,1:73],
                  LIA_MCA$SIM_full$SIM_low_101_P[1:48,1:73]*PLOTTING_VARIABLES$ls_mask[1:48,1:73])

plot_low_P <- STACYmap(gridlyr = Plot_lyr, colorscheme = RColorBrewer::brewer.pal(9, "Reds"), zoom = c(-180,-60,180,73), legend_names = list(grid = "years BP")) + 
  ggtitle("Driest 101y average P") + theme(plot.title = element_text(h = 0.5), panel.border = element_blank(),
                                           legend.background = element_blank(),
                                           axis.text = element_blank())

Plot_lyr <- rbind(LIA_MCA$SIM_full$SIM_high_101_P[49:96,1:73]*PLOTTING_VARIABLES$ls_mask[49:96,1:73],
                  LIA_MCA$SIM_full$SIM_high_101_P[1:48,1:73]*PLOTTING_VARIABLES$ls_mask[1:48,1:73])

plot_high_P <- STACYmap(gridlyr = Plot_lyr, colorscheme = RColorBrewer::brewer.pal(9, "Reds"), zoom = c(-180,-60,180,73), legend_names = list(grid = "years BP")) + 
  ggtitle("Wettest 101y average P") + theme(plot.title = element_text(h = 0.5), panel.border = element_blank(),
                                            legend.background = element_blank(),
                                            axis.text = element_blank())


plot <- ggarrange(plot_low_T, plot_high_T, plot_low_P, plot_high_P,
                  labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2)

plot  %>% ggsave(filename = paste('Paper_Plot_8_LIA_MCA', 'pdf', sep = '.'), plot = ., path = 'Plots/Paper', 
                 width = PLOTTING_VARIABLES$WIDTH, height = PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")

remove(plot_low_T, plot_high_T, plot_low_P, plot_high_P, plot)


###################################################################################################
######################################### APPENDIX ################################################
###################################################################################################

#################################################
## MEAN SCATTER PLOTS ###########################
#################################################

## What can we try?

# 1) elevation_lat and mean diff

# TODO: Hier muss noch Gesteinsart rein, damit man zwischen Calcit und Aragonit unterscheiden kann!
# Oder man mach Masken!

scatter_data = array(dim = c(length(MEAN_ANALYSIS$MEAN_DIFF_ISOT$entity_id[mask_mean]),13))
scatter_data[,1] = MEAN_ANALYSIS$MEAN_DIFF_ISOT$entity_id[mask_mean]
scatter_data[,3] = MEAN_ANALYSIS$MEAN_DIFF_ISOT$diff_full[mask_mean]
scatter_data[,4] = MEAN_ANALYSIS$MEAN_DIFF_ISOT$diff_down[mask_mean]
for(ii in 1:length(MEAN_ANALYSIS$MEAN_DIFF_ISOT$entity_id[mask_mean])){
  entity = scatter_data[[ii,1]]
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  index = match(site, DATA_past1000$CAVES$site_info$site_id)
  scatter_data[ii,2] = site
  scatter_data[ii,5] = DATA_past1000$CAVES$site_info$elevation[index]
  scatter_data[ii,6] = DATA_past1000$CAVES$site_info$latitude[index]
  scatter_data[ii,7] = mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$TEMP, na.rm = T)
  scatter_data[ii,8] = mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$PREC, na.rm = T)
  scatter_data[ii,9] = mean(DATA_past1000$CAVES$sim_data_seasonal[[paste0("CAVE", site)]]$WINTER$prec_mean, na.rm = T)
  scatter_data[ii,10] = mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$SLP, na.rm = T)
  scatter_data[ii,11] = as.numeric(DATA_past1000$CAVES$elevation_cave_sim$`sim-cave`[DATA_past1000$CAVES$elevation_cave_sim$entity_id == entity])
  scatter_data[ii,12] = as.numeric(DATA_past1000$CAVES$entity_info$distance_entrance[match(entity, DATA_past1000$CAVES$entity_info$entity_id)])
  scatter_data[ii,13] = as.numeric(DATA_past1000$CAVES$entity_info$geology[match(entity, DATA_past1000$CAVES$entity_info$entity_id)])
}

colnames(scatter_data) = c("entity_id", "site_id", "diff_full", "diff_down", "elevation", "latitude", 
                           "mean_temp", "mean_prec", "winter_mean_prec", "SLPR", "elevation_diff", "dist_entrance", "geology")
scatter_data = as.tibble(scatter_data)

#View(scatter_data)
# mask for aragonite and calcite


mask_mean_calcite = logical(length = length(MEAN_ANALYSIS$MEAN_DIFF_ISOT$entity_id[mask_mean]))
mask_mean_aragonite  = logical(length = length(MEAN_ANALYSIS$MEAN_DIFF_ISOT$entity_id[mask_mean]))

for(ii in 1:length(MEAN_ANALYSIS$MEAN_DIFF_ISOT$entity_id[mask_mean])){
  entity = MEAN_ANALYSIS$MEAN_DIFF_ISOT$entity_id[mask_mean][ii]
  if(DATA_past1000$CAVES$entity_info$mineralogy[DATA_past1000$CAVES$entity_info$entity_id == entity] == "calcite") {mask_mean_calcite[ii] = T}
  else{mask_mean_aragonite[ii] = T}
}


##diff_full

pdf(file = "Plots/Paper/Appendix_01_ScatterMean_diff-full.pdf", width = 8, height = 8)
par(mfrow=c(3,2),oma = c(1,3,0,0) + 0.1,mar = c(3,0,1,1) + 0.1)
# elevation
plot(scatter_data$elevation[mask_mean_calcite], scatter_data$diff_full[mask_mean_calcite], xlab = "", ylab = "", ylim = c(-10,10), xlim = c(0,4000), panel.first = grid())
abline(h=0)
points(scatter_data$elevation[mask_mean_aragonite], scatter_data$diff_full[mask_mean_aragonite], pch = 8, col = "blue")
lines(lowess(scatter_data$elevation, scatter_data$diff_full, f = 2/3, delta = 0.01*diff(range(scatter_data$elevation, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "elevation",side = 1,line = 2)
mtext(text = "d18O-d18Oc",side = 2,line = 2)
text(350, 5.7, "calcite")
text(1150, -7., "aragonite", col = "blue")

plot(scatter_data$elevation_diff[mask_mean_calcite], scatter_data$diff_full[mask_mean_calcite], yaxt = 'n', xlab = "", ylab = "", ylim = c(-10,10), panel.first = grid())
abline(h=0)
points(scatter_data$elevation_diff[mask_mean_aragonite], scatter_data$diff_full[mask_mean_aragonite], pch = 8, col = "blue")
lines(lowess(scatter_data$elevation_diff, scatter_data$diff_full, f = 2/5, delta = 0.01*diff(range(scatter_data$elevation, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "elevation diff. (sim-rec)",side = 1,line = 2)

plot(scatter_data$latitude[mask_mean_calcite], scatter_data$diff_full[mask_mean_calcite], xlab = "", ylab = "", ylim = c(-10,10), panel.first = grid())
points(scatter_data$latitude[mask_mean_aragonite], scatter_data$diff_full[mask_mean_aragonite], pch = 8, col = "blue")
abline(h=0)
lines(lowess(scatter_data$latitude, scatter_data$diff_full, f = 2/3, delta = 0.01*diff(range(scatter_data$latitude, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "latitude",side = 1,line = 2)
mtext(text = "d18O-d18Oc",side = 2,line = 2)

plot(scatter_data$mean_temp[mask_mean_calcite], scatter_data$diff_full[mask_mean_calcite], yaxt = 'n', xlab = "", ylab = "", ylim = c(-10,10), panel.first = grid())
abline(h=0)
points(scatter_data$mean_temp[mask_mean_aragonite], scatter_data$diff_full[mask_mean_aragonite], pch = 8, col = "blue")
lines(lowess(scatter_data$mean_temp, scatter_data$diff_full, f = 2/3, delta = 0.01*diff(range(scatter_data$mean_temp, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "mean temp",side = 1,line = 2)

plot(scatter_data$mean_prec[mask_mean_calcite]*8.6148e4, scatter_data$diff_full[mask_mean_calcite], xlab = "", ylab = "", log = "x", ylim = c(-10,10), xlim = c(0.2,15), panel.first = grid(equilogs = FALSE))
abline(h=0)
points(scatter_data$mean_prec[mask_mean_aragonite]*8.6148e4, scatter_data$diff_full[mask_mean_aragonite], pch = 8, col = "blue")
lines(lowess(scatter_data$mean_prec*8.6148e4, scatter_data$diff_full, f = 2/3, delta = 0.01*diff(range(scatter_data$winter_mean_prec, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "mean prec",side = 1,line = 2)
mtext(text = "d18O-d18Oc",side = 2,line = 2)

plot(scatter_data$winter_mean_prec[mask_mean_calcite]*8.6148e4, scatter_data$diff_full[mask_mean_calcite], yaxt = 'n', xlab = "", ylab = "", xlim = c(0.2,15), log = "x", ylim = c(-10,10), yaxt = "n", panel.first = grid(equilogs = FALSE))
abline(h=0)
points(scatter_data$winter_mean_prec[mask_mean_aragonite]*8.6148e4, scatter_data$diff_full[mask_mean_aragonite], pch = 8, col = "blue")
lines(lowess(scatter_data$winter_mean_prec*8.6148e4, scatter_data$diff_full, f = 2/3, delta = 0.01*diff(range(scatter_data$winter_mean_prec, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "winter prec",side = 1,line = 2)



dev.off()

##down_sampled

pdf(file = "Plots/Paper/Appendix_01_ScatterMean_diff-down.pdf", width = 8, height = PLOTTING_VARIABLES$HEIGHT/2)
par(mfrow=c(2,2),oma = c(1,3,0,0) + 0.1,mar = c(3,0,1,1) + 0.1)
plot(scatter_data$elevation[mask_mean_calcite], scatter_data$diff_down[mask_mean_calcite], xlab = "", ylab = "", ylim = c(-10,10), xlim = c(0,4000), panel.first = grid())
points(scatter_data$elevation[mask_mean_aragonite], scatter_data$diff_down[mask_mean_aragonite], pch = 8, col = "blue")
lines(lowess(scatter_data$elevation, scatter_data$diff_down, f = 2/3, delta = 0.01*diff(range(scatter_data$elevation, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "elevation",side = 1,line = 2)
mtext(text = "diff-d18O",side = 2,line = 2)
text(350, 5.7, "calcite")
text(1150, -7., "aragonite", col = "blue")
plot(scatter_data$latitude[mask_mean_calcite], scatter_data$diff_down[mask_mean_calcite], xlab = "", ylab = "", ylim = c(-10,10), panel.first = grid(), yaxt ="n")
points(scatter_data$latitude[mask_mean_aragonite], scatter_data$diff_down[mask_mean_aragonite], pch = 8, col = "blue")
lines(lowess(scatter_data$latitude, scatter_data$diff_down, f = 2/3, delta = 0.01*diff(range(scatter_data$latitude, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "latitude",side = 1,line = 2)

plot(scatter_data$mean_temp[mask_mean_calcite], scatter_data$diff_down[mask_mean_calcite], xlab = "", ylab = "", ylim = c(-10,10), panel.first = grid())
points(scatter_data$mean_temp[mask_mean_aragonite], scatter_data$diff_down[mask_mean_aragonite], pch = 8, col = "blue")
lines(lowess(scatter_data$mean_temp, scatter_data$diff_down, f = 2/3, delta = 0.01*diff(range(scatter_data$mean_temp, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "mean temp",side = 1,line = 2)
mtext(text = "diff-d18O",side = 2,line = 2)
plot(scatter_data$winter_mean_prec[mask_mean_calcite], scatter_data$diff_down[mask_mean_calcite], xlab = "", ylab = "", log = "x", ylim = c(-10,10), yaxt = "n", panel.first = grid(equilogs = FALSE))
points(scatter_data$winter_mean_prec[mask_mean_aragonite], scatter_data$diff_down[mask_mean_aragonite], pch = 8, col = "blue")
lines(lowess(scatter_data$winter_mean_prec, scatter_data$diff_down, f = 2/3, delta = 0.01*diff(range(scatter_data$winter_mean_prec, na.rm = T))), lwd = 4, col = "#B2182B")
mtext(text = "winter prec",side = 1,line = 2)
dev.off()


