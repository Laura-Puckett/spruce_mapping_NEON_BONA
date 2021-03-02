#*********************************************************************************#
# -----0: Packages
# ********************************************************************************#
list.of.packages <- c("rgeos","dplyr", "rgdal", "sp", "raster")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages, repos='http://cloud.r-project.org')

library(rgeos)
library(raster)
library(rgdal)
library(dplyr)
library(sp)
library(neonUtilities)
library(geoNEON)
library(stringr)
library(ggplot2)

#*********************************************************************************#
# -----1: Inputs
# ********************************************************************************#
setwd('C:/Users/lkp58/Documents/research/')
output_dir='./analyses/BONA_species_mapping/AVIRIS/'

AVIRIS_path = '//nau.froot.nau.edu/cirrus/scratch/lkp58/datasets/AVIRIS/sites/BONA/2_gdalwarp_to_tiff/ang20180723t193341_corr_v2r2_img.tif'
LVIS_path = './analyses/BONA_species_mapping/LVIS/3_selected_LVIS_with_AVIRIS.shp'
# UTM_6N = CRS("+proj=utm +zone=6 +datum=WGS84 +units=m +no_defs")
# veg_path = './datasets/NEON/BONA/veg_shp/NEON_BONA_veg.shp'

# above_crs =  sp::CRS("+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# ABoVE projection: Canada Albers Equal Area Conic
# https://above.nasa.gov/implementation_plan/standard_projection.html
# Canada_Albers_Equal_Area_Conic
# WKID (EPSG): 102001 Authority: ESRI


#################################################################
#### read NEON data
#################################################################
sites = readOGR(LVIS_path)

#################################################################
#### read AVIRIS data
#################################################################
img = brick(AVIRIS_path)

#iterate through each LVIS footprint

summary_df = NULL
for (i in 1:nrow(sites)){
  print(i)
  site = sites[i,]
  overlap = raster::extract(img, site, weights = T, normalizeWeights = T)[[1]]
  values = overlap[,colnames(overlap) != "weight"]
  weight = overlap[,"weight"]
  rm(overlap)
  
  if(all(is.na(values))==F){ # check if overlap has values
    bands = colnames(values)
    bands = stringr::str_remove(bands, ".Nanometers")
    bands = stringr::str_remove(bands, "X")
    bands = as.numeric(bands) %>% round(digits = 2)
    bands
    
    # make a quick viz of the data while it's in an easy format for that
    values_for_plot = cbind("Nanometers" = bands, t(values))
    rownames(values_for_plot) = NULL
    
    values_for_plot = values_for_plot %>%
      as.data.frame() %>%
      gather("pixel", "reflectance", -Nanometers) %>%
      mutate(reflectance = as.numeric(reflectance))
    
    # p1 = ggplot(values_for_plot, aes(x = Nanometers, y = reflectance)) +
    #   geom_point(aes(color = pixel)) +
    #   ggtitle(paste0("unfiltered, plot # ", i))
    # print(p1)
    
    values_for_plot = values_for_plot %>%
      filter((Nanometers >= 1330 & Nanometers <= 1435) == F,
             (Nanometers >= 1814 & Nanometers <= 1960) == F,
             ((Nanometers >= 2480) == F)) %>%
      mutate(reflectance = ifelse(reflectance < 0, NA, reflectance),
             refelctance = ifelse(reflectance > 1, NA, reflectance))
    
    p2 = ggplot(values_for_plot, aes(x = Nanometers, y = reflectance)) +
      geom_point(aes(color = pixel)) +
      ggtitle(paste0("filtered, plot # ", i))

    # multiply values by weights and sum for each band
    weighted_mean = crossprod(values, weight)
    rownames(weighted_mean) = NULL
    weighted_mean = cbind("Nanometers" = bands, "reflectance" = weighted_mean) %>%
      as.data.frame()
    colnames(weighted_mean) = c("Nanometers", "reflectance")
    

    
    weighted_mean = weighted_mean %>%
      filter((Nanometers >= 1330 & Nanometers <= 1435) == F,
             (Nanometers >= 1814 & Nanometers <= 1960) == F,
             ((Nanometers >= 2480) == F)) %>%
      mutate(reflectance = ifelse(reflectance < 0, NA, reflectance),
             reflectance = ifelse(reflectance > 1, NA, reflectance))
    
    # add the weighted mean to the plot from earlier that included spectra of multiple pixels
    p3 = p2 + 
      geom_point(data = weighted_mean, aes(x = Nanometers, y = reflectance))
    
    print(p3)
    # swich rows and columns
    weighted_mean_tidy = weighted_mean %>%
      select(-Nanometers) %>% (t) %>% as.data.frame()
    colnames(weighted_mean_tidy) = paste0(weighted_mean$Nanometers, "_nm")
    rownames(weighted_mean_tidy) = NULL
    
    rm(weighted_mean)

    summary_df_row = cbind(site@data, weighted_mean_tidy)
    summary_df = rbind(summary_df, summary_df_row)
    
    }
}
# # or do the same without a loop
# overlap = raster::extract(img, sites, weights = T, normalizeWeights = T)


#### write output
write.csv(x = summary_df, file = paste0(output_dir, '3a_AVIRIS_for_LVIS_footprints_with_field_data.csv'))


summary_df = read.csv(paste0(output_dir, '3a_AVIRIS_for_LVIS_footprints_with_field_data.csv'))
