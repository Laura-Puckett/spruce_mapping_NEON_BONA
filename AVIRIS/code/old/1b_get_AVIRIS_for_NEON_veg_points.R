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
veg_path = './datasets/NEON/BONA/veg_shp/NEON_BONA_veg.shp'

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
sites = readOGR(veg_path)

#################################################################
#### read AVIRIS data
#################################################################
img = brick(AVIRIS_path)

#iterate through each LVIS footprint

output_df = NULL
for (i in 1:nrow(sites)){
  # print(i)
  site = sites[i,]
  overlap = raster::extract(img, site, cellnumbers = T) %>% as.data.frame()
  cell_num = overlap$cells
  values = overlap[colnames(overlap) != "cells"]
  # print(values)
  # check if the site overlapped with a pixel, only continue if so
  if(is.na(cell_num)==F){ 
    # if the values are all na, don't continue
    if(all(is.na(values))==F){
      
      bands = colnames(values)
      bands = stringr::str_remove(bands, ".Nanometers")
      bands = stringr::str_remove(bands, "X")
      bands = as.numeric(bands) %>% round(digits = 2)
      bands
      
      # make a quick viz of the data while it's in an easy format for that
      values_for_plot = cbind("Nanometers" = bands, t(values))
      rm(values)
      
      rownames(values_for_plot) = NULL
      
      values_for_plot = values_for_plot %>%
        as.data.frame() %>%
        gather("pixel", "reflectance", -Nanometers) %>%
        mutate(reflectance = as.numeric(reflectance)) %>%
        select(-pixel) # there's only one pixel for this comparison. 
      
      p1 = ggplot(values_for_plot, aes(x = Nanometers, y = reflectance)) +
        geom_point() +
        ggtitle(paste0("unfiltered, plot # ", i))
      # print(p1)
      values_for_plot = values_for_plot %>%
        filter((Nanometers >= 1330 & Nanometers <= 1435) == F,
               (Nanometers >= 1814 & Nanometers <= 1960) == F,
               ((Nanometers >= 2480) == F)) 
      
      # # which wavelengths have data outside 0,1?
      # print(values_for_plot[(values_for_plot$reflectance < 0 | values_for_plot$reflectance > 1),])
      
      p2 = ggplot(values_for_plot, aes(x = Nanometers, y = reflectance)) +
        geom_point() +
        ggtitle(paste0("filtered, plot # ", i))
      
      p2
      
      values_for_output = values_for_plot %>%
        spread("Nanometers","reflectance")
      
      colnames(values_for_output) = paste0(colnames(values_for_output), "_nm")
      
      # # swich rows and columns
      # values_tidy = values %>%
      #   as.data.frame() 
      # colnames(values_tidy) = colnames(values) %>%
      #   stringr::str_remove("X") %>%
      #   stringr::str_replace(".Nanometers","_nm")
      # rownames(values_tidy) = NULL
      
      output_df_row = cbind("uid" = site@data$uid_x, "cell_number" = cell_num, values_for_output)
      output_df = rbind(output_df, output_df_row)
    }
  }
}
# # or do the same without a loop
# overlap = raster::extract(img, sites, weights = T, normalizeWeights = T)


#################################################################
#### write output
#################################################################
write.csv(x = output_df, file = paste0(output_dir, '3b_AVIRIS_for_NEON_veg_points.csv'))

#################################################################
#### summarize by pixel
#################################################################
veg_path = './datasets/biomass/NEON/BONA/veg.rds'

AVIRIS_veg_df = read.csv(paste0(output_dir, '3b_AVIRIS_for_NEON_veg_points.csv')) %>%
  select(-X) %>%
  filter(is.na(cell_num) == F)

# this will be for joining to the summary values per pixel later on
AVIRIS_no_veg = AVIRIS_veg_df %>%
  select(-uid) %>%
  distinct()

# relevant veg data from NEON field measurements
veg = readRDS(veg_path) %>%
  select(uid.x, stemDiameter, height, taxonID, plantStatus, date.x) %>%
  rename("uid" = uid.x,
         "date" = date.x)

# summarize veg data by AVIRIS pixel 
summary_df = AVIRIS_veg_df %>%
  select(cell_number, uid) %>%
  filter(is.na(cell_number) == F) %>%
  full_join(veg, by = "uid") %>%
  filter(substr(plantStatus,0,4) == "Live") %>%
  mutate(ind_BA = 3.14*(stemDiameter/200)^2,
         multiplier = ifelse(taxonID == "PIMA", 1, 0), # PIMA: picea mariana (black spruce)
         bs_ind_BA = ind_BA * multiplier) %>%
  group_by(cell_number) %>%
  dplyr::mutate(n_stems = n(),
                   BA = sum(ind_BA)/25 * 10000, # in units of m2/ha
                   bs_BA = sum(bs_ind_BA)/25 * 10000,
                   prop_bs = bs_BA/BA,
                   prop_bs = ifelse(BA == 0, 0, prop_bs) # set proportion of bs to 0 if there's no biomass
                   ) %>%
  ungroup()

summary_df = summary_df %>%
  inner_join(AVIRIS_no_veg, by = "cell_number")

write.csv(summary_df, paste0(output_dir, '4b_AVIRIS_NEON_for_modeling.csv'))
