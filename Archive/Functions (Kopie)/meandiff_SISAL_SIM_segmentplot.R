library(ggplot2)
library(dplyr)
library(tibble)

plotfct_time <- function(data, ncol = 1, pyr = F) {
  if (pyr) {
    data <- group_by(data, entity_id) %>% mutate(mins = median(interp_age, na.rm = T)) %>% arrange(mins)
  }
  
  data <- data %>% group_by(entity_id) %>% tidyr::nest() %>% rownames_to_column(var = 'rentity_id') %>% mutate_at(vars(rentity_id), as.numeric) %>% tidyr::unnest()
  brks <- (length(unique(data$rentity_id)) / ncol)
  brks <- seq(0, max(data$rentity_id), by = brks)
  data <- filter(data, !is.na(diff_d18O)) %>% mutate(entity_cat = cut(rentity_id, breaks = brks))
  
  labs <- unique(data$entity_id)[seq(1, length(unique(data$entity_id)), 1)]
  labbrks <- unique(data$rentity_id)[seq(1, length(unique(data$rentity_id)), 1)]
  
  plot <- ggplot(data = data,
                 mapping = aes(x = interp_age, y = rentity_id, group = rentity_id)) + 
    geom_line(mapping = aes(color = diff_d18O), position = position_identity(), size = 1.7, stat = 'unique') + 
    scale_color_gradient2(low = "#053061", high = "#67001F", 
                          guide = guide_colorbar(title = 'diff d18O Sim-Rec', barwidth = 10)) +
    scale_y_reverse(breaks = labbrks, labels = labs) + 
    facet_wrap(~ entity_cat, ncol = ncol, scales = 'free') +
    labs(y = 'entity id', x = 'years BP') + 
    theme(strip.text = element_blank(), axis.text.y = element_text(size = 6), panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
          legend.direction = 'horizontal', legend.position = 'bottom')
  return(plot)
}

#load('/home/ariana/Documents/full Run 2/gr_v1b.RData')
#gr <- gr_new %>% select(entity_id, bacon_age, bacon_gr) %>% group_by(entity_id) %>% mutate(bacon_gr = abs(bacon_gr)) %>% filter(bacon_gr < 100)
#plotfct(gr, 5, T) %>% 
#  ggsave(filename = '~/Documents/code/STACYplots/gr_SISAL.pdf', width = 30, height = 20, units = 'cm')






