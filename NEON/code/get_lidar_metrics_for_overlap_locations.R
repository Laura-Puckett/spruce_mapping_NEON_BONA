setwd('C:/Users/lkp58/Documents/research/')
source('./analyses/BONA_species_mapping/code/structural_diversity_metrics.R')

options(stringsAsFactors=F)
library(neonUtilities);library(geoNEON);library(sp);library(ggplot2);library(dplyr);
library(ggforce);library(rgdal);library(raster);library(lidR);library(gstat);
library(data.table);library(kableExtra)

# ctg = lidR::readLAScatalog('./datasets/NEON/BONA/DP1.30003.001/2019/FullSite/D19/2019_BONA_3/L1/DiscreteLidar/ClassifiedPointCloud/')
# ctg = saveRDS(ctg, './analyses/BONA_species_mapping/ctg.rds')
ctg = readRDS('./analyses/BONA_species_mapping/ctg.rds')
plot(ctg)
veg = readOGR('./datasets/NEON/BONA/veg_shp/NEON_BONA_veg.shp')
veg_data %>% group_by(indvdID) %>% arrange(date_x, descending = T) %>% slice(1) %>% ungroup()

scale = 5
metrics = data.frame()
for(i in 1:nrow(veg@data)){
  print(i)
  tree = veg[i,]
  buffer_small = rgeos::gBuffer(tree, width = scale/2, byid = T, capStyle = 'SQUARE')
  buffer_large = rgeos::gBuffer(tree, width = scale/2 + 5, byid = T, capStyle = 'SQUARE')
  ctg_intsct = lidR::catalog_intersect(ctg, buffer_large)
  # plot(ctg_intsct)
  # plot(buffer_large, add = T)
  
  las = readLAS(files = ctg_intsct)
  las = lidR::clip_roi(las, buffer_large)
  Z_sd=sd(las@data$Z)
  Z_mean=mean(las@data$Z)
  f = paste("-drop_z_below",(Z_mean-4*Z_sd),"-drop_z_above",(Z_mean+4*Z_sd))
  
  las = readLAS(files = ctg_intsct[1,], filter = f)
  las = lidR::clip_roi(las, buffer_large)
  # plot(las)
  
  # calculate dem and normalize heights to that over entire buffered region
  dtm <- grid_terrain(las, 1, kriging(k = 10L))
  las_norm <- lasnormalize(las, dtm)
  
  
  las_norm_clipped = clip_roi(las_norm, buffer_small)
    
  chm <- grid_canopy(las_norm_clipped, res = 1, dsmtin()) 
  # plot(chm)
  
  metrics_tmp = structural_diversity_metrics(las_norm_clipped, tree@data$adjEstn, tree@data$adjNrth)
  metrics_tmp = cbind(metrics_tmp, "scale" = 5)
  rbind(metrics, metrics_tmp)
}

veg_colnames_to_keep = c("indvdID", "plotID","date_x","grwthFr","plntStt","stmDmtr",
                         "mstmntH","rcrdTyp","scntfcn","adjDcmlLt","adjDcmLn","adjElvt")

veg_data = veg@data %>% select(veg_colnames_to_keep)

all_data = rbind(metrics, veg_data)

