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

#*********************************************************************************#
# -----1: Inputs
# ********************************************************************************#
source('//nau.froot.nau.edu/cirrus/home/LVIS/LVIS_colnames.R')

sensor = "facility_2019"
output_dir='./analyses/BONA_species_mapping/LVIS/'
LVIS_dir = './datasets/lidar/LVIS/BONA/0_dl/'

UTM_6N = CRS("+proj=utm +zone=6 +datum=WGS84 +units=m +no_defs")
radius=5 # distance for buffering LVIS

# above_crs =  sp::CRS("+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# ABoVE projection: Canada Albers Equal Area Conic
# https://above.nasa.gov/implementation_plan/standard_projection.html
# Canada_Albers_Equal_Area_Conic
# WKID (EPSG): 102001 Authority: ESRI

if(sensor == "facility_2017"){
  LVIS_colnames = facility_2017_colnames 
  
}else if(sensor == "facility_2019"){
  LVIS_colnames = facility_2019_colnames
  
}else if(sensor == "classic_2019"){
  LVIS_colnames = classic_2019_colnames
}

#################################################################
#### read NEON data
#################################################################
NEON_sites = readOGR('./datasets/biomass/NEON/NEON_TOS_Plot_Polygons_view-shp/NEON_TOS_Plot_Polygons_view.shp')
NEON_sites = NEON_sites[NEON_sites@data$siteID == "BONA" & NEON_sites@data$subtype == "basePlot",]
NEON_sites = spTransform(NEON_sites, UTM_6N)

#################################################################
#### Read LVIS TXT data to points
#########################################################################
files = list.files(LVIS_dir, recursive = T, full.names=T, pattern = "TXT$|txt$")
nfiles = length(files)
LVIS_overlap = NULL
# NEON_sites = gBuffer(NEON_sites, width = 5, capStyle = "SQUARE", byid = T)
for(fileNum in 1:nfiles){
  filepath = files[fileNum]
  print(paste0(fileNum/nfiles, " ", fileNum, " ", filepath))
  
  
  LVIS = read.table(paste0(filepath), skip = 0, header = FALSE)
  colnames(LVIS) = LVIS_colnames 
  LVIS = LVIS %>%
    mutate(GLON = ifelse(GLON >0, -1*(360-GLON), GLON))
  
  LVIS_points = sp::SpatialPointsDataFrame(cbind(LVIS$GLON, LVIS$GLAT), LVIS, proj4string = sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  rm(LVIS)
  LVIS_points = spTransform(LVIS_points, UTM_6N)
  
  
  # # a quick n dirty way to buffer the LVIS file extent before clipping field points to that area
  # # The objective here is to deal with the edge case where a field plot area might overlap with 
  # # LVIS data, but the plot center is slightly outside the extent of the LVIS data
  # 
  # LVIS_buffered_extent = extent(LVIS_points)
  # LVIS_buffered_extent@xmin=LVIS_buffered_extent@xmin-radius 
  # LVIS_buffered_extent@xmax=LVIS_buffered_extent@xmax+radius 
  # LVIS_buffered_extent@ymin=LVIS_buffered_extent@ymin-radius 
  # LVIS_buffered_extent@ymax=LVIS_buffered_extent@ymax+radius 
  # sites_cropped = raster::crop(NEON_sites, LVIS_buffered_extent)
  # if(length(field_cropped)>0){
  
  
  #iterate through each plot
  for (i in 1:nrow(NEON_sites)){
    site = NEON_sites[i,]
    LVIS_overlap_tmp = raster::intersect(LVIS_points, site)
    
    if(length(LVIS_overlap_tmp)>0){
      print('intersection')
      plot(site)
      plot(LVIS_overlap_tmp, add = T, col = "green")
      LVIS_overlap_tmp@data = cbind(LVIS_overlap_tmp@data, "filename" = filepath, "sensor" = sensor)
      
      if(is.null(LVIS_overlap)){
        LVIS_overlap = LVIS_overlap_tmp
      }else{
        LVIS_overlap = rbind(LVIS_overlap, LVIS_overlap_tmp)}
    }
  }
}

LVIS_overlap_footprints = rgeos::gBuffer(LVIS_overlap, width = 5, capStyle = "ROUND", byid = T)


#### write output
writeOGR(obj = LVIS_overlap_footprints,
         layer = "LVIS_overlap",
         driver = "ESRI Shapefile",
         dsn = paste0(output_dir, '1_',sensor, '_footprints_near_field.shp'),
         overwrite_layer = T)

writeOGR(obj = LVIS_overlap,
         layer = "LVIS_overlap",
         driver = "ESRI Shapefile",
         dsn = paste0(output_dir, '1_',sensor, '_points_near_field.shp'),
         overwrite_layer = T)

if(!is.null(LVIS_overlap)){
  write.csv(x = LVIS_overlap@data, file = paste0(output_dir, '1_',sensor, '_near_field.csv'))
}
