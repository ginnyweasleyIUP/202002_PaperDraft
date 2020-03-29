#################################################
## Paper Network 3 ##############################
#################################################

# Boxplot

line_names = 1
pdf(file = paste0("Plots/Paper_Plot_6_Network_c.pdf"), height= PLOTTING_VARIABLES$HEIGHT, width = PLOTTING_VARIABLES$WIDTH)
par(mar = c(2,9,2,2))
plot(c(-1,1), c(1,25+3), type = "n", axes = FALSE, xlab = "", ylab = "" )
abline(v=0)
## SITES
site_corr <- list(full = list(), gauss = list())
site_list <- DATA_past1000$CAVES$entity_info %>% select(site_id, entity_id) %>%
  filter(entity_id %in% ANALYSIS$NETWORK$entity_meta$entity_id) %>% group_by(site_id) %>% count() %>% filter(n>1)
for(site in site_list$site_id){
  site_corr$full <- c(site_corr$full, ANALYSIS$NETWORK$SITES[[paste0("SITE",site)]]$C[ANALYSIS$NETWORK$SITES[[paste0("SITE",site)]]$P<0.1])
  site_corr$gauss <- c(site_corr$gauss, ANALYSIS$NETWORK$SITES[[paste0("SITE",site)]]$C_gauss[ANALYSIS$NETWORK$SITES[[paste0("SITE",site)]]$P_gauss<0.1])
}
site_corr$full = as.numeric(site_corr$full)
site_corr$gauss = as.numeric(site_corr$gauss)

boxplot(site_corr$full, add = TRUE, at = 25+3,boxwex = 1, names = "n", horizontal = T, outline = F) 
boxplot(site_corr$gauss, add = TRUE, at = 24.5+3, boxwex = 1, names = "n", horizontal = T, col = "grey", axes = F, outline = F)

#mtext(side=2,"GMST",                cex = unitscex,    line = unitslinno, las = 1, col = "black", at = 1)
mtext(side = 2, "sites (27/12)", cex = 1, line = line_names, las = 1, col = "black", at = 25+3)
abline(h=23.75+3, lty = 3)
##GRIGBOX

corr <- list(full = list(), gauss = list())
list <- ANALYSIS$NETWORK$entity_meta %>% select(gridbox_id, entity_id) %>% group_by(gridbox_id) %>% count() %>% filter(n>1)
for(gridbox in list$gridbox_id){
  corr$full <- c(corr$full, ANALYSIS$NETWORK$GRIDBOX[[paste0("GRIDBOX",gridbox)]]$C[ANALYSIS$NETWORK$GRIDBOX[[paste0("GRIDBOX",gridbox)]]$P<0.1])
  corr$gauss <- c(corr$full, ANALYSIS$NETWORK$GRIDBOX[[paste0("GRIDBOX",gridbox)]]$C_gauss[ANALYSIS$NETWORK$GRIDBOX[[paste0("GRIDBOX",gridbox)]]$P_gauss<0.1])
}
corr$full = as.numeric(corr$full)
corr$gauss = as.numeric(corr$gauss)

boxplot(corr$full, add = TRUE, at = 23+3 ,boxwex = 1, names = "n", horizontal = T, axes = F, outline = F) 
boxplot(corr$gauss, add = TRUE, at = 22.5+3, boxwex = 1, names = "n", horizontal = T, col = "grey", axes = F, outline = F)

#mtext(side=2,"GMST",                cex = unitscex,    line = unitslinno, las = 1, col = "black", at = 1)
mtext(side = 2, "gridbox (45/18)", cex = 1, line = line_names, las = 1, col = "black", at = 23+3)
abline(h=21.75+3, lty = 3)


position <- list(cluster6 = c(17, 16.5, 16, 15.5)+4+3,
                 cluster2 = c(14, 13.5, 13, 12.5)+4+3,
                 cluster3 = c(11, 10.5, 10, 9.5)+4+3,
                 cluster7 = c(8, 7.5, 7, 6.5)+4+3,
                 cluster1 = c(5, 4.5, 4, 3.5)+4+3,
                 cluster5 = c(2, 1.5, 1, 0.5)+4+3,
                 cluster9 = c(2, 1.5, 1, 0.5)+4)
text <- list(cluster1 = "India", cluster2 = "SouthAm", cluster3 = "Europe", cluster4 = "Africa", 
             cluster5 = "China", cluster6 = "NorthAm", cluster7 = "Arabia", cluster8 = "NZ", cluster9 = "SE Asia")
cluster_list <- ANALYSIS$NETWORK$entity_meta %>% select(cluster_id, entity_id) %>% group_by(cluster_id) %>% count() %>% filter(n>1)
cluster_number <- c(6,2,3,4,7,1,5,9,8)

## CLUSTER
#for(cluster in c(1,2,3,5,6,7)){
for(cluster in c(6,2,3,7,1,5,9)){
  corr <- list(full = as.numeric(ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER",cluster)]]$C[ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER",cluster)]]$P<0.1]), 
               gauss = list(as.numeric(ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER",cluster)]]$C_gauss[ANALYSIS$NETWORK$CLUSTER[[paste0("CLUSTER",cluster)]]$P_gauss<0.1])),
               full_sim = list(as.numeric(c(ANALYSIS$NETWORK$CLUSTER_SIM_a[[paste0("CLUSTER",cluster)]]$C[ANALYSIS$NETWORK$CLUSTER_SIM_a[[paste0("CLUSTER",cluster)]]$P<0.1],
                                            ANALYSIS$NETWORK$CLUSTER_SIM_b[[paste0("CLUSTER",cluster)]]$C[ANALYSIS$NETWORK$CLUSTER_SIM_b[[paste0("CLUSTER",cluster)]]$P<0.1],
                                            ANALYSIS$NETWORK$CLUSTER_SIM_c[[paste0("CLUSTER",cluster)]]$C[ANALYSIS$NETWORK$CLUSTER_SIM_c[[paste0("CLUSTER",cluster)]]$P<0.1]))), 
               gauss_sim = list(as.numeric(c(ANALYSIS$NETWORK$CLUSTER_SIM_a[[paste0("CLUSTER",cluster)]]$C_gauss[ANALYSIS$NETWORK$CLUSTER_SIM_a[[paste0("CLUSTER",cluster)]]$P_gauss<0.1],
                                             ANALYSIS$NETWORK$CLUSTER_SIM_b[[paste0("CLUSTER",cluster)]]$C_gauss[ANALYSIS$NETWORK$CLUSTER_SIM_b[[paste0("CLUSTER",cluster)]]$P_gauss<0.1],
                                             ANALYSIS$NETWORK$CLUSTER_SIM_c[[paste0("CLUSTER",cluster)]]$C_gauss[ANALYSIS$NETWORK$CLUSTER_SIM_c[[paste0("CLUSTER",cluster)]]$P_gauss<0.1]))))
  
  boxplot(corr$full, add = TRUE, at = position[[paste0("cluster",cluster)]][1] ,boxwex = 1, names = "n", horizontal = T, axes = F, outline = F) 
  boxplot(corr$gauss, add = TRUE, at = position[[paste0("cluster",cluster)]][2] , boxwex = 1, names = "n", horizontal = T, col = "grey", axes = F, outline = F)
  boxplot(corr$full_sim, add = TRUE, at = position[[paste0("cluster",cluster)]][3]  ,boxwex = 1, names = "n", horizontal = T, col = "dodgerblue3", axes = F, outline = F) 
  boxplot(corr$gauss_sim, add = TRUE, at = position[[paste0("cluster",cluster)]][4] , boxwex = 1, names = "n", horizontal = T, col = adjustcolor("dodgerblue3", alpha.f = 0.5), axes = F, outline = F)
  
  mtext(side = 2, paste0("c",cluster_number[cluster],"(",text[[cluster]],") (",cluster_list$n[cluster],")"), cex = 1, line = line_names, las = 1, col = "black", at = position[[paste0("cluster",cluster)]][1])
  
}

abline(h=3.75, lty =3)

##Global

corr$full = as.numeric(as.numeric(ANALYSIS$NETWORK$GLOBAL$C[ANALYSIS$NETWORK$GLOBAL$P<0.1]))
corr$gauss = as.numeric(as.numeric(ANALYSIS$NETWORK$GLOBAL$C_gauss[ANALYSIS$NETWORK$GLOBAL$P_gauss<0.1]))
corr$full_sim = as.numeric(ANALYSIS$NETWORK$GLOBAL_SIM_a$C[ANALYSIS$NETWORK$GLOBAL_SIM_a$P<0.1], 
                           ANALYSIS$NETWORK$GLOBAL_SIM_b$C[ANALYSIS$NETWORK$GLOBAL_SIM_b$P<0.1],
                           ANALYSIS$NETWORK$GLOBAL_SIM_c$C[ANALYSIS$NETWORK$GLOBAL_SIM_c$P<0.1])
corr$gauss_sim = as.numeric(ANALYSIS$NETWORK$GLOBAL_SIM_a$C_gauss[ANALYSIS$NETWORK$GLOBAL_SIM_a$P_gauss<0.1], 
                            ANALYSIS$NETWORK$GLOBAL_SIM_b$C_gauss[ANALYSIS$NETWORK$GLOBAL_SIM_b$P_gauss<0.1],
                            ANALYSIS$NETWORK$GLOBAL_SIM_c$C_gauss[ANALYSIS$NETWORK$GLOBAL_SIM_c$P_gauss<0.1])

boxplot(corr$full, add = TRUE, at = 3 ,boxwex = 1, names = "n", horizontal = T, axes = F) 
boxplot(corr$gauss, add = TRUE, at = 2.5 , boxwex = 1, names = "n", horizontal = T, col = "grey", axes = F)
boxplot(corr$full_sim, add = TRUE, at = 2,boxwex = 1, names = "n", horizontal = T, col = "dodgerblue3", axes = F, outline = F) 
boxplot(corr$gauss_sim, add = TRUE, at = 1.5 , boxwex = 1, names = "n", horizontal = T, col = adjustcolor("dodgerblue3", alpha.f = 0.5), axes = F, outline = F)

mtext(side = 2, paste0("global (",dim(ANALYSIS$NETWORK$GLOBAL$C)[[1]],")"), cex = 1, line = line_names, las = 1, col = "black", at = 2)

legend(-1.075,28.5+3,xpd = T,inset=-0.2, bty='n', x.intersp=0.5,text.width=c(1,0.25,0.4, 0.37),
       c("record","record 100gauss","xnap(a/b/c)", "xnap(a/b/c) 100gauss"), fill=c("white", "grey", "dodgerblue3", adjustcolor("dodgerblue3", alpha.f = 0.5)), horiz=TRUE, cex=1)

dev.off()