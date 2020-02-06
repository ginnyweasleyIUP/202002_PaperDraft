library(ggplot2)
library(ggrepel)
library(plotly)
library(maps)
library(RColorBrewer)
library(tidyverse)
library(magrittr)
library(hash)


#source('util_file_processing.R')
#source('util_global_plotting_options.R')
#source('util_taxa_homogenization.R')
#source('plot_maps.R')
#source('util_pca_analysis.R')
#source('plot_pca_analysis.R')
#source('constants.R')
#source('plot_dataset_overview.R')
#source('util_dating.R')

DIR_DATASETS <- "~/Dokumente/01_Promotion/07_R_Code/201907_SpatialAnalysis"


# time series class (with information from overview)
CAVES_dataset <- setRefClass('CAVES_dataset', 
                            fields = list(site_info = 'data.frame',
                                          climate_info = 'data.frame',
                                          raw_data = 'list',
                                          yearly_data = 'list',
                                          reg_data = 'data.frame',
                                          corr_data = 'data.frame'),
                            methods = list(
                              initialize = function(){
                                site_info <<- as.tibble(sisal_data)
                                raw_data <<-  list(time = 'data.frame',
                                                   temp = 'data.frame',
                                                   prec = 'data.frame', 
                                                   isot = 'data.frame')
                                yearly_data <<-  list(time = 'data.frame',
                                                   temp = 'data.frame',
                                                   prec = 'data.frame', 
                                                   isot = 'data.frame')
                              }
                            )
                            )


# initialize_raw_data <- function(NoC){
#   my_list = c()
#   for (ii in 1:NoC){
#     nam <- paste('cave', ii, sep='')
#     input <- tibble(temp = get(paste('cave', ii, '_temp', sep = '')), 
#                     prec = get(paste('cave', ii, '_prec', sep = '')), 
#                     isot = get(paste('cave', ii, '_isot', sep = '')))
#     assign(nam, input)
#     new_list = c(my_list, nam)
#     my_list = new_list
#   }
#   return(my_list)
# }

# CAVES_dataset$methods = list(
#   initialize = function(){
#     site_info <- as.tibble(cave_data)
#   }
# )
    
    # raw_data$cave_id <- seq(1:40)
    #  for (ii in 1:40){
    #    raw_data$temp[ii] = paste('cave', ii, '_temp', sep= '')
    #    raw_data$prec[ii] = paste('cave', ii, '_prec', sep= '')
    #    raw_data$isot[ii] = paste('cave', ii, '_isot', sep= '')
    #  }
    
  


# ACER_dataset <- setRefClass('ACER_dataset', 
#                             fields = list(db = 'list', 
#                                           pollen_data_percentized = 'logical',
#                                           PCA_calculated = 'logical', 
#                                           nw = 'list'),
#                             methods = list(
#                               initialize = function(){
#                                 db <<- list('sites' = load_db_flat('site'), 
#                                             'sample_dating' = load_db_flat('sample_ct_corrected'), 
#                                             'pollen_data' = load_db_flat('pollen_data'), 
#                                             'charcoal_data' = load_db_flat('charcoal_data'), 
#                                             'biome_percentages' = load_db_flat('biome_percentages'))
#                                 
#                                 pollen_data_percentized <<- FALSE
#                                 PCA_calculated <<- FALSE
#                               }
#                             )
#)

#load_db_flat <- function(filename, dir = DIR_DATASETS){
#  path <- paste(dir, filename, sep = '/')
#  print(paste('loading flat database file:', filename))
#  return(as.tibble(read.csv(file = path, header = TRUE, sep = ',', encoding = 'UTF-8')))
#}

# load_db_raw <-function(NumoCaves){
#   raw_data$cave_id <- 1:NumoCaves
#   for (ii in 1:NumoCaves){
#     raw_data$temp[ii] = paste('cave', ii, '_temp', sep= '')
#     raw_data$prec[ii] = paste('cave', ii, '_prec', sep= '')
#     raw_data$isot[ii] = paste('cave', ii, '_isot', sep= '')
#   }
# }


# # data processing
# ACER_dataset$methods(pollen_percentize = function(overwrite_original_col = FALSE, in_place = TRUE){
#   if(!pollen_data_percentized){
#     pollen <- inner_join(db$pollen_data, select(db$sample_dating, sample_id, count_type), by = 'sample_id')
#     
#     coun <- pollen[pollen$count_type == 'COUN',] %>% group_by(sample_id) %>% summarise(total_taxon_count = sum(taxon_count)) %>% 
#       full_join(., pollen[pollen$count_type == 'COUN',], by = 'sample_id') %>% 
#       {if(overwrite_original_col){
#         mutate(., taxon_count = taxon_count / total_taxon_count * 100.0)
#       }
#       else{
#         mutate(., taxon_pcnt = taxon_count / total_taxon_count * 100.0)
#       }} %>% select(-total_taxon_count)
#     pcnt <- pollen[pollen$count_type %in% list('DIGI', 'PCNT'),] %>% 
#       {if(overwrite_original_col){.}
#         else{
#           mutate(., taxon_pcnt = taxon_count)
#         }}
#     
#     if(overwrite_original_col){
#       pollen <- full_join(coun, pcnt, by = c('pollen_data_id', 'sample_id', 'site_id', 'taxon', 'taxon_count', 'count_type')) %>%
#         select(pollen_data_id, sample_id, site_id, taxon, taxon_count)
#       if(in_place){
#         db$pollen_data <<- select(pollen, pollen_data_id, sample_id, site_id, taxon, taxon_count)
#         pollen_data_percentized <<- TRUE
#         cat('wrote taxon percentages to db$pollen_data overwriting existing taxon_count\n')
#       }
#       else{
#         cat('calculated taxon percentages to pollen_data\n')
#         return(pollen)
#       }
#     }
#     else{
#       pollen <- full_join(coun, pcnt, by = c('pollen_data_id', 'sample_id', 'site_id', 'taxon', 'taxon_count', 'taxon_pcnt', 'count_type')) %>%
#         select(pollen_data_id, sample_id, site_id, taxon, taxon_count, taxon_pcnt)
#       if(in_place){
#         db$pollen_data <<- select(pollen, pollen_data_id, sample_id, site_id, taxon, taxon_count, taxon_pcnt)
#         pollen_data_percentized <<- TRUE
#         cat('wrote taxon percentages to db$pollen_data\n')
#       }
#       else{
#         cat('calculated taxon percentages to pollen_data\n')
#         return(pollen)
#       }
#     }
#   }
#   else if(in_place){
#     cat('wrote taxon percentages to db$pollen_data\n')
#     return(db$pollen_data)
#   }
# })
# 
# 
# # PCA
# ACER_dataset$methods(pca = function(type_data = 'pollen', rank. = NULL, scale. = FALSE, windowing = NULL, inplace = TRUE){
#   # windowing = list(windows = list(0, 12000, 19000, 27000), labels = list('0-12', '12-19', '19-27), compare_windows = list('0-12', '19-27'))
#   db_files <- list(pollen = 'pollen_data', harm_pollen = 'harm_pollen_data')
#   dir <- paste(DIR_CACHE, paste('pca', type_data, sep = '_'), sep = '/')
#   if(is.null(rank.)){filename <- 'no_limit'}else {filename <- as.character(rank.)}
#   if(scale.) {filename <- paste0(filename, '_scaled')}else {filename <- paste0(filename, '_unscaled')}
#   if(!is.null(windowing)) {filename <- paste0(filename, '_win_', windowing$compare_windows[1], '_', windowing$compare_windows[2])}
#   if(!dir.exists(dir)){dir.create(dir)}
#   if(!(type_data %in% names(db_files))) {stop('unknown type_data provided to ACER_dataset$pca')}
#   
#   if (type_data == 'pollen' | type_data == 'harm_pollen') {
#     pollen_percentize()
#     if(file.exists(paste(paste(dir, filename, sep = '/'), DATAFILES_TYPE, sep = '.'))){
#       cat(paste('loading PCA for', type_data, 'from file', rank., '\n'))
#       pca <- load_db_flat(filename = filename, dir = dir) %>% 
#         select(pca_id, site_id, sample_id, pc, pc_value, variance, var_exp, cum_var_exp)
#     }
#     else{
#       cat('calculating \n')
#       if(!is.null(windowing)) {
#         pollen_data <- inner_join(db[[db_files[[type_data]]]], select(db$sample_dating, site_id, sample_id, mixed_age), by = c('site_id', 'sample_id')) %>% 
#           window_data(data = ., type_data = type_data, transform = NULL, windows = windowing$windows, labels = windowing$labels) %>% 
#           ungroup(.) %>% 
#           filter(window %in% windowing$compare_windows) %>% 
#           select(-mixed_age, -window)
#       }
#       else {pollen_data <- db[[db_files[[type_data]]]]}
#       pca <- pca_pollen(pollen_data = pollen_data, type_data = type_data, rank. = rank., scale. = scale.)
#       cat(paste('writing to file', filename, '\n'))
#       write.csv(x = pca,
#                 file = paste(paste(dir, filename, sep = '/'), DATAFILES_TYPE, sep = '.'),
#                 sep = ',', fileEncoding = 'UTF-8')
#     }
#     PCA_calculated <<- TRUE
#   }
#   
#   if(inplace){
#     db[[paste('pca', type_data, filename, sep = '_')]] <<- pca
#   }
#   else{
#     return(pca)
#   }
# 
# })
# 
# 
# # call only working from outside the Class
# ACER_dataset$methods(get_plot_pca_var = function(pca_data = db$pca_data){
#   plot <- get_plot_pca_var(pca_data = pca_data)
#   return(plot)
# })
# 
# 
# # plot site locations and elevations on map
# ACER_dataset$methods(get_plot_sites_geographic = function(bathymetry = TRUE, 
#                                                           topo = FALSE,
#                                                           pnv = FALSE, lgm_ice = FALSE, 
#                                                           names = FALSE,
#                                                           site_id = FALSE, ISR = FALSE,
#                                                           save_plot = NULL, 
#                                                           legend_inside = TRUE, 
#                                                           labs = FALSE,
#                                                           return = TRUE){
#   map_plot <- plot_ACER_sites_on_map(sites = db$sites, db$sample_dating,
#                                      bathy = bathymetry, topo = topo, pnv = pnv, lgm_ice = lgm_ice,
#                                      names = names, site_id = site_id,
#                                      ISR = ISR, save_plot = save_plot, 
#                                      legend_inside = legend_inside)
#   if (labs) {map_plot <- map_plot + labs(title = 'ACER sampling sites')}
#   if (return) {return(map_plot)}
# })
# 
# 
# # plot sampling & dating info for all sites -> age comparison & overlapping; number of samples, sampling density
# ACER_dataset$methods(get_pollen_sampling_info_plot = function(){
#   pollen <- inner_join(db$sites, db$sample_dating, by = 'site_id')
#   
#   plot <- ggplot(data = pollen, mapping = aes(x = site_id, color = depth * 0.1, alpha = I(0.5), size = I(0.1))) +
#     geom_point(data = filter(pollen, CLAM_best >= 0.0), mapping = aes(y = CLAM_best)) +
#     geom_point(data = filter(pollen, CLAM_best < 0.0), mapping = aes(y = original_age)) + theme_minimal() + scale_y_reverse() + 
#     labs(title = '', y = 'age / yrs B.P.', x = 'site id') +
#     scale_color_gradientn(colors = rev(brewer.pal(11, 'RdYlBu'))) + 
#     guides(colour = global_colorbar(title = 'depth below groundlevel [m]'))
#   return(plot)
# })
# 
# 
# ACER_dataset$methods(get_multi_site_plot = function(type_data = 'pollen', y_vars_subset = 'all',
#                                                     y_vars_top_n = 'all', site_names = 'all',
#                                                     site_ids = 'all',
#                                                     time_restrict = NULL, lat_restrict = NULL, 
#                                                     plot_errors = FALSE, use_bars = FALSE, use_lines = FALSE, use_points = TRUE,
#                                                     mark_events = NULL, mark_time_periods = NULL, 
#                                                     high_res_only = FALSE,
#                                                     plot_ice_cores = list(activate = FALSE),
#                                                     interactive = list('activate' = FALSE),
#                                                     save_plot = NULL){
#   
#   types <- list(pollen = list(data = 'pollen_data', signal_name = 'taxon_pcnt', color_fill = 'taxon', color_fill_title = 'taxon',
#                               title = 'taxa time series',
#                               subtitle = paste('sorting by site latitude', {if(!is.null(time_restrict)) {paste0(' - time restriction ', time_restrict$upper/1000, ' - ', time_restrict$lower/1000, ' kyrs b.p.')}}, 
#                                                {if(y_vars_top_n != 'all') {paste0(' - most common taxa: ', y_vars_top_n)}}, 
#                                                {if(!is.null(lat_restrict)) {paste0(subtitle, ' - latitude restrict: ', lat_restrict$max, ' to ', lat_restrict$min, ' degrees')}}), 
#                               x_title = 'age [ky b.p.]', 
#                               y_title = 'taxon percentage / %'), 
#                 harmonized_pollen = list(data = 'harm_pollen_data', signal_name = 'taxon_pcnt', color_fill = 'taxon_harmonized', color_fill_title = 'taxon harmonized', 
#                               title = 'harmonized taxa time series',
#                               subtitle = paste('sorting by site latitude', {if(!is.null(time_restrict)) {paste0(' - time restriction ', time_restrict$upper/1000, ' - ', time_restrict$lower/1000, ' kyrs b.p.')}}, 
#                                                {if(y_vars_top_n != 'all') {paste0(' - most common taxa: ', y_vars_top_n)}}, 
#                                                {if(!is.null(lat_restrict)) {paste0(subtitle, ' - latitude restrict: ', lat_restrict$max, ' to ', lat_restrict$min, ' degrees')}}), 
#                               x_title = 'age [ky b.p.]', 
#                               y_title = 'harmonized taxon percentage / %'), 
#                 arboreal_pollen = list(data = 'arb_pollen_data', signal_name = 'pcnt_arb_pollen', 
#                                        title = 'arboreal pollen time series',
#                                        subtitle = paste('sorting by site latitude', {if(!is.null(time_restrict)) {paste0(' - time restriction ', time_restrict$upper/1000, ' - ', time_restrict$lower/1000, ' kyrs b.p.')}},  
#                                                         {if(!is.null(lat_restrict)) {paste0(subtitle, ' - latitude restrict: ', lat_restrict$max, ' to ', lat_restrict$min, ' degrees')}}, 
#                                                         {if(high_res_only) {' - only high resolution sites considered'}}), 
#                                        x_title = 'age [ky b.p.]', 
#                                        y_title = 'percentage arboreal pollen / %'), 
#                 logit_arboreal_pollen = list(data = 'logit_arb_pollen_data', signal_name = 'logit_pcnt_arb_pollen', 
#                                        title = 'logit transformed arboreal pollen time series',
#                                        subtitle = paste('sorting by site latitude', {if(!is.null(time_restrict)) {paste0(' - time restriction ', time_restrict$upper/1000, ' - ', time_restrict$lower/1000, ' kyrs b.p.')}},  
#                                                         {if(!is.null(lat_restrict)) {paste0(subtitle, ' - latitude restrict: ', lat_restrict$max, ' to ', lat_restrict$min, ' degrees')}}, 
#                                                         {if(high_res_only) {' - only high resolution sites considered'}}), 
#                                        x_title = 'age [ky b.p.]', 
#                                        y_title = 'transformed percentage arboreal pollen / a.u.'), 
#                 probit_arboreal_pollen = list(data = 'probit_arb_pollen_data', signal_name = 'probit_pcnt_arb_pollen', 
#                                              title = 'probit transformed arboreal pollen time series',
#                                              subtitle = paste('sorting by site latitude', {if(!is.null(time_restrict)) {paste0(' - time restriction ', time_restrict$upper/1000, ' - ', time_restrict$lower/1000, ' kyrs b.p.')}},  
#                                                               {if(!is.null(lat_restrict)) {paste0(subtitle, ' - latitude restrict: ', lat_restrict$max, ' to ', lat_restrict$min, ' degrees')}}, 
#                                                               {if(high_res_only) {' - only high resolution sites considered'}}), 
#                                              x_title = 'age [ky b.p.]', 
#                                              y_title = 'transformed percentage arboreal pollen / a.u.'), 
#                 detrend_logit_arboreal_pollen = list(data = 'detrend_logit_arb_pollen_data', signal_name = 'detrend_logit_pcnt_arb_pollen', 
#                                              title = 'logit transformed and 10-ky-windowed-linear detrended arboreal pollen time series',
#                                              subtitle = paste('sorting by site latitude', {if(!is.null(time_restrict)) {paste0(' - time restriction ', time_restrict$upper/1000, ' - ', time_restrict$lower/1000, ' kyrs b.p.')}},  
#                                                               {if(!is.null(lat_restrict)) {paste0(subtitle, ' - latitude restrict: ', lat_restrict$max, ' to ', lat_restrict$min, ' degrees')}}, 
#                                                               {if(high_res_only) {' - only high resolution sites considered'}}), 
#                                              x_title = 'age [ky b.p.]', 
#                                              y_title = 'transf. & detr. percentage arboreal pollen / a.u.'), 
#                 detrend_probit_arboreal_pollen = list(data = 'detrend_probit_arb_pollen_data', signal_name = 'detrend_probit_pcnt_arb_pollen', 
#                                                      title = 'probit transformed and 10-ky-windowed-linear detrended arboreal pollen time series',
#                                                      subtitle = paste('sorting by site latitude', {if(!is.null(time_restrict)) {paste0(' - time restriction ', time_restrict$upper/1000, ' - ', time_restrict$lower/1000, ' kyrs b.p.')}},  
#                                                                       {if(!is.null(lat_restrict)) {paste0(subtitle, ' - latitude restrict: ', lat_restrict$max, ' to ', lat_restrict$min, ' degrees')}}, 
#                                                                       {if(high_res_only) {' - only high resolution sites considered'}}), 
#                                                      x_title = 'age [ky b.p.]', 
#                                                      y_title = 'transf. & detr. percentage arboreal pollen / a.u.'),
#                 biomes = list(data = 'biome_percentages', signal_name = 'percentages', color_fill = 'biome_long', color_fill_title = 'biome',
#                               title = 'biome percentage time series',
#                               subtitle = paste('sorting by site latitude', {if(!is.null(time_restrict)) {paste0(' - time restriction ', time_restrict$upper/1000, ' - ', time_restrict$lower/1000, ' kyrs b.p.')}}, 
#                                                {if(y_vars_top_n != 'all') {paste0(' - most common taxa: ', y_vars_top_n)}}, 
#                                                {if(!is.null(lat_restrict)) {paste0(subtitle, ' - latitude restrict: ', lat_restrict$max, ' to ', lat_restrict$min, ' degrees')}}), 
#                               x_title = 'age [ky b.p.]', 
#                               y_title = 'biome percentage / %'), 
#                 pca_pollen = list(data = 'pca_data', signal_name = 'pc_value', color_fill = 'pc', color_fill_title = 'principal component',
#                                   title = 'principal components time series',
#                                   subtitle = paste('sorting by site latitude', {if(!is.null(time_restrict)) {paste0(' - time restriction ', time_restrict$upper/1000, ' - ', time_restrict$lower/1000, ' kyrs b.p.')}}, 
#                                                    {if(!is.null(lat_restrict)) {paste0(subtitle, ' - latitude restrict: ', lat_restrict$max, ' to ', lat_restrict$min, ' degrees')}}), 
#                                   x_title = 'age [ky b.p.]', 
#                                   y_title = 'pc value / a.u.'), 
#                 hybrid_ice_core = list(data = 'hic_data', signal_name = 'hic', 
#                                        title = 'hybrid ice core',
#                                        subtitle = paste('sorting by site latitude', {if(!is.null(time_restrict)) {paste0(' - time restriction ', time_restrict$upper/1000, ' - ', time_restrict$lower/1000, ' kyrs b.p.')}},  
#                                                         {if(!is.null(lat_restrict)) {paste0(subtitle, ' - latitude restrict: ', lat_restrict$max, ' to ', lat_restrict$min, ' degrees')}}, 
#                                                         {if(high_res_only) {' - only high resolution sites considered'}}), 
#                                        x_title = 'age [ky b.p.]', 
#                                        y_title = 'hybrid isotope signal'), 
#                 detrend_identity_hybrid_ice_core = list(data = 'hic_data', signal_name = 'detrend_identity_hic', 
#                                        title = 'hybrid ice core variability',
#                                        subtitle = paste('sorting by site latitude', {if(!is.null(time_restrict)) {paste0(' - time restriction ', time_restrict$upper/1000, ' - ', time_restrict$lower/1000, ' kyrs b.p.')}},  
#                                                         {if(!is.null(lat_restrict)) {paste0(subtitle, ' - latitude restrict: ', lat_restrict$max, ' to ', lat_restrict$min, ' degrees')}}, 
#                                                         {if(high_res_only) {' - only high resolution sites considered'}}), 
#                                        x_title = 'age [ky b.p.]', 
#                                        y_title = 'hybrid isotope signal'))
#   
#   if(!is.null(lat_restrict)){
#     site_names <- lat_restrict_sites(sites = db$sites, lat = lat_restrict)
#   }
#   if(str_detect(type_data, 'pca') & str_detect(type_data, 'pollen')) {types[['pca_pollen']][['data']] <- type_data; type_data <- 'pca_pollen'}
#   if(str_detect(type_data, 'hybrid_ice_core')) {if(str_detect(type_data, 'detrend') & str_detect(type_data, 'identity')) {types[['detrend_identity_hybrid_ice_core']][['data']] <- str_replace(type_data, 'hybrid_ice_core', 'hic_data'); type_data <- 'detrend_identity_hybrid_ice_core'}else {types[['hybrid_ice_core']][['data']] <- str_replace(type_data, 'hybrid_ice_core', 'hic_data'); type_data <- 'hybrid_ice_core'}}
#   if(str_detect(type_data, 'detrend') & str_detect(type_data, 'probit_arboreal_pollen')) {types[['detrend_probit_arboreal_pollen']][['data']] <- paste0(str_replace(type_data, 'arboreal', 'arb'), '_data'); type_data <- 'detrend_probit_arboreal_pollen'}
#   if(str_detect(type_data, 'detrend') & str_detect(type_data, 'logit_arboreal_pollen')) {types[['detrend_logit_arboreal_pollen']][['data']] <- paste0(str_replace(type_data, 'arboreal', 'arb'), '_data'); type_data <- 'detrend_logit_arboreal_pollen'}
#   if(!(type_data %in% names(types))) {stop(paste('unknown type_data', type_data, 'passed to ACER_dataset$get_multi_site_plot'))}
#   if(!(types[[type_data]]$data %in% names(db))) {stop('given type not yet in ACER_dataset$db, need to call a ..._to_db fct first')}
#   mix_sample_dating()
#   print(type_data)
#   # type specific data preparation
#   if(type_data == 'pollen' | type_data == 'harmonized_pollen') {
#     pollen_percentize()
#     data <- db[[types[[type_data]]$data]]
#     
#     if(y_vars_top_n != 'all'){
#       tmp <- get_most_common_taxa(pollen_dated = data, top_n = y_vars_top_n, taxon_signal = types[[type_data]]$color_fill)
#       y_vars_subset <- tmp$use_taxon
#       data <- tmp$pollen
#     }
#     data <- filter_pollen_for_taxa_sites(pollen = data %>% inner_join(select(db$sites, site_id, site_name), by = c('site_id')),
#                                          site_names = site_names, taxon_signal = types[[type_data]]$color_fill,
#                                          use_taxon = y_vars_subset, most_common_taxa = y_vars_top_n) %>% 
#       select(-site_name)
#     }
#   else if(type_data == 'biomes') {if(y_vars_top_n != 'all'){warning('variable y_vars_top_n has no effect for type_data = biome')}; data <- db[[types[[type_data]]$data]] %>% inner_join(BIOME_PERCENTAGES_ABBR, by = 'biome')}
#   else if(type_data == 'pca_pollen') {
#     if(!PCA_calculated) {stop('no PCA data found - need to calculate PCA first using ACER$get_pca')}
#     else{
#       data <- db[[types[[type_data]]$data]]
#       
#       if(!is.null(y_vars_top_n)){
#         data <- data %>% 
#           filter(., as.character(.$pc) %in% paste('pc', seq(1, y_vars_top_n), sep = ''))
#       }
#     }
#   }
#   else if(str_detect(type_data, 'arboreal_pollen')) {if(y_vars_top_n != 'all'){warning('variable y_vars_top_n has no effect for type_data = arboreal_pollen')}
#     data <- db[[types[[type_data]]$data]]}
#   else if(str_detect(type_data, 'hybrid_ice_core')) {if(y_vars_top_n != 'all'){warning('variable y_vars_top_n has no effect for type_data = hybrid_ice_core')}
#     data <- db[[types[[type_data]]$data]]}
#   else {stop('plotting not implemented for given type_data')}
#   
#   # general data preparation
#   sites <- db$sites %>% 
#   {if (!(site_names == 'all')) {filter(., as.character(site_name) %in% site_names)}else {.}} %>% 
#   {if (!(site_ids == 'all')) {filter(., site_id %in% site_ids)}else {.}} %>% 
#     inner_join(categorize_age_interval_stats(compute_age_interval_stats(db$sample_dating, ranges = list(list(start = 0, stop = 73000)))), by = 'site_id') %>%
#     {if (high_res_only) {filter(., status_med_smp_res %in% c('0 - 200', '200 - 300'))}else {.}}
#   
#   data <- sites_by_lat(sites) %>% 
#     inner_join(db$sample_dating, by = 'site_id') %>% 
#     inner_join(data, by = c('sample_id', 'site_id')) %>% 
#     time_restrict_data(., time_restrict) %>% 
#     filter(!is.na(UQ(sym(types[[type_data]]$signal_name))))
#   
#   if (length(distinct(data, site_id)$site_id) > 50) {
#     data <- data %>% 
#       mutate(new_name = as.character(site_id)) %>% 
#       select(-site_name) %>% 
#       rename(site_name = new_name) #%>% 
#       sites_by_lat(.)
#   }
#   
#   # make plot
#   plot <- {if('color_fill' %in% names(types[[type_data]])) {ggplot(data = data, mapping = aes(x = mixed_age/1000, y = UQ(sym(types[[type_data]]$signal_name)), color = UQ(sym(types[[type_data]]$color_fill)), alpha = I(0.7)))}
#     else {ggplot(data = data, mapping = aes(x = mixed_age/1000, y = UQ(sym(types[[type_data]]$signal_name))))}} + 
#     guides(colour = global_legend(title = types[[type_data]]$color_fill_title)) + 
#     {if('color_fill' %in% names(types[[type_data]])) {scale_color_viridis_d(alpha = 0.8, begin = 0.5, end = 0.9)}}
#   
#   if(use_bars){
#     plot <- plot +
#       {if('color_fill' %in% names(types[[type_data]])) {geom_bar(stat = 'identity', aes(fill = UQ(sym(types[[type_data]]$color_fill)), alpha = I(0.7)), show.legend = c(color = TRUE, fill = TRUE))}
#         else {geom_bar(stat = 'identity', color = GLOBAL_BLUE_LIGHT, alpha = I(0.7))}} + 
#       {if('color_fill' %in% names(types[[type_data]])) {scale_fill_viridis_d(alpha = 0.8, begin = 0.5, end = 0.9)}} + 
#       guides(fill = global_legend(title = types[[type_data]]$color_fill_title))
#   }
#   else{
#     if(use_points) {
#       plot <- plot +
#         {if('color_fill' %in% names(types[[type_data]])) {geom_point(size = I(0.5))}else {geom_point(color = GLOBAL_BLUE_LIGHT, alpha = I(0.7), size = I(0.5))}}
#     }
#     if(plot_errors){
#       plot <- plot + 
#       {if('color_fill' %in% names(types[[type_data]])) {geom_errorbarh(data = filter(data, CLAM_best >= 0), mapping = aes(xmax = CLAM_max95/1000, xmin = CLAM_min95/1000, size = I(0.5), height = I(1), alpha = I(0.7)))}
#         else {geom_errorbarh(data = filter(data, CLAM_best >= 0), mapping = aes(xmax = CLAM_max95/1000, xmin = CLAM_min95/1000, size = I(0.5), height = I(1), alpha = I(0.7)), color = GLOBAL_BLUE_LIGHT)}}
#     }
#     if(use_lines & !plot_ice_cores$activate){
#       plot <- plot + 
#         {if('color_fill' %in% names(types[[type_data]])) {geom_line()}else {geom_line(color = GLOBAL_BLUE_LIGHT, alpha = I(0.7))}}
#     }
#   }
#   
#   # eventual ice plotting
#   if (plot_ice_cores$activate) {
#     if (use_bars) {warning('plot_ice_cores only if !use_bars')}
#     else {
#       if(!all(c('transform', 'detrend') %in% names(plot_ice_cores$process))) {stop('some processing steps for ice core signal undefined')}
#       if('window' %in% names(data)) {windows <- ungroup(data) %>% distinct(window) %>% .$window %>% as.character(.) %>% sort(.) %>% as.list(.)}
#       else {windows <- list('0-150')}
#       name <- paste(plot_ice_cores$process$transform, 'hic', sep = '_')
#       #if(plot_ice_cores$process$detrend$activate) {name <- paste('detrend', name, sep = '_')}
#       data_oc <- signal.orig.ice.cores() %>% 
#         mutate_at(vars(sample_id), as.numeric) %>% 
#         mutate(hic_id = sample_id) %>% 
#         mutate(sample_id = add(sample_id, 100000)) %>% 
#         correct_data_stage(data = ., orig_data = 'hybrid_ice_core', transform = plot_ice_cores$process$transform, detrend = plot_ice_cores$process$detrend, unnest_data = TRUE, 
#                            windows = window.char.to.num(windows), labels = windows, sd_one = TRUE, signal_id = 'hic_id', signal_name = 'hic', remove_age = FALSE) %>% 
#         mutate(site_name = if_else(site_id == 101, 'NGRIP', 'EPICA'),
#                type = if_else(site_id == 101, 'NGRIP', 'EPICA')) %>% 
#         ungroup() %>% 
#         #{if (!str_detect(type_data, 'hybrid_ice_core')) {rename(., UQ(sym()))}else {.}}
#         {if(plot_ice_cores$process$detrend$activate == FALSE & str_detect(type_data, 'hybrid_ice_core')) {select(., -UQ(sym(types[[type_data]]$signal_name))) %>% rename(UQ(sym(types[[type_data]]$signal_name)) := UQ(sym(name)))}else if(plot_ice_cores$process$detrend$activate == FALSE){rename(., UQ(sym(types[[type_data]]$signal_name)) := identity_hic)} else {.}}
#       data <- mutate(data, type = types[[type_data]]$signal_name, site_id2 = site_id) %>% 
#         group_by(site_id2) %>% 
#         do((full_join(., mutate(data_oc, site_id = .$site_id %>% unique(), site_name = .$site_name %>% unique(), lat = .$lat %>% unique())))) %>% ungroup %>% select(-site_id2)
#       if (!str_detect(type_data, 'hybrid_ice_core')) {
#         data <- group_by(data, site_id, type) %>% 
#           mutate(UQ(sym(types[[type_data]]$signal_name)) := (UQ(sym(types[[type_data]]$signal_name)) - mean(UQ(sym(types[[type_data]]$signal_name)))) / sd(UQ(sym(types[[type_data]]$signal_name)))) %>% 
#           ungroup() %>% 
#           {if (!is.null(time_restrict)) {filter(., mixed_age <= time_restrict$lower & mixed_age >= time_restrict$upper)}else {.}} %>% 
#           filter(type != 'EPICA')
#       }
#       data$type <- ordered(data$type, levels = c(types[[type_data]]$signal_name, 'NGRIP', 'EPICA'))
#       upper_time_lim <- filter(data, !(type %in% c('NGRIP', 'EPICA'))) %>% .$mixed_age %>% max(.)
#       data <- filter(data, mixed_age <= upper_time_lim) %>% select(site_id, UQ(sym(types[[type_data]]$signal_name)), type) %>% 
#         group_by(site_id) %>% 
#         filter(., !(type %in% c('NGRIP', 'EPICA'))) %>%
#         summarise(max = max(UQ(sym(types[[type_data]]$signal_name))), min = min(UQ(sym(types[[type_data]]$signal_name)))) %>% 
#         ungroup() %>% 
#         inner_join(data, by = 'site_id') %>% 
#         filter(UQ(sym(types[[type_data]]$signal_name)) >= min & UQ(sym(types[[type_data]]$signal_name)) <= max)
#         
#       plot <- plot + 
#         geom_line(data = data, aes(color = type, alpha = type)) + 
#         scale_color_manual(values = c(GLOBAL_BLUE_LIGHT, GLOBAL_RED_LIGHT, GLOBAL_RED_DARK), guide = global_legend('signal type')) + 
#         scale_alpha_manual(values = c(1, 0.4, 0.4), guide = FALSE)
#     }
#   }
#   
#   plot <- plot +
#   {if(length(distinct(data, site_id)$site_id) > 50) {facet_wrap(vars(site_name), scales = 'free_y')}
#     else {facet_grid(interaction(site_id, site_name, sep = ' - ') + round(lat, digits = 2) ~ ., scales = 'free_y', switch = 'x')}} + 
#     labs(title = types[[type_data]]$title, subtitle = types[[type_data]]$subtitle, x = types[[type_data]]$x_title, y = types[[type_data]]$y_title) + 
#     global_title_and_axis() + 
#     theme(legend.box = 'horizontal', legend.position = 'bottom', legend.direction = 'horizontal', strip.text.y = element_text(face = GLOBAL_FONT_FACE_TITLE, size = (GLOBAL_FONT_SIZE - 8), vjust = 0))
#   
#   if(!is.null(save_plot)){
#     global_save_plot(save_plot = save_plot, plot = plot)
#   }
#   
#   return(global_interactive_plot(ggplot = plot, interactive = interactive))
#   })



#source('dataset_methods.R')

