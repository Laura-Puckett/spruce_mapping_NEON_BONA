library(stringr); library(ggplot2); library(raster); library(rgdal); library(sp)

LVIS_overlap_path = './analyses/BONA_species_mapping/LVIS/output/1_facility_2019_near_field.csv'
output_dir = './/analyses/BONA_species_mapping/LVIS/output/'
subplots_path = './analyses/BONA_species_mapping/NEON/data/subplots/subplots.shp'
LVIS_crs = sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
UTM_6N = sp::CRS("+proj=utm +zone=6 +datum=WGS84 +units=m +no_defs")
LVIS_footprint_radius = 5
res = 1 # resolution for gridding the LVIS footprints for getting area weighted mean
PLOT = FALSE
metrics = c("ZG","RH10","RH15","RH20","RH25","RH30","RH35","RH40","RH45","RH50",
            "RH55","RH60","RH65","RH70","RH75", "RH80","RH85","RH90","RH95",
            "RH96","RH97","RH98","RH99","RH100") # what to grab from LVIS file (everything rn)

#################################################################
### read data
#################################################################
LVIS = read.csv(LVIS_overlap_path)
sites = readOGR(subplots_path)


# # remove waveforms that don't have any trees
# LVIS = LVIS %>%
#   filter(RH100>2)

# ************************************************************#
### Grid and get area-weighted mean of LVIS footprints per site
# ************************************************************#
output_df = data.frame()
for(i in 1:nrow(sites)){
  print(i)
  site = sites[i,]
  site_id = paste0(site@data$plotID, '_', site@data$sbpltID)
  site_data = LVIS %>% filter(plotID == site@data$plotID)
  summary_data = data.frame()
  if(nrow(site_data) > 0){
    LVIS_points = sp::SpatialPointsDataFrame(cbind(site_data$GLON, site_data$GLAT, site_data$ZG), site_data, proj4string = LVIS_crs)
    
    # transform to projected coordinate system and buffer to footprint width
    LVIS_points = spTransform(LVIS_points, UTM_6N)
    LVIS_footprints = rgeos::gBuffer(spgeom = LVIS_points, byid=T, width = LVIS_footprint_radius)
    
    # force extent to be in even multiples of res. This will allow
    # alignment between the output rasters from different sites.
    points_extent = extent(LVIS_footprints)
    new_extent = points_extent
    new_extent[1] = floor(points_extent[1]/res)*res
    new_extent[2] = ceiling(points_extent[2]/res)*res
    new_extent[3] = floor(points_extent[3]/res)*res
    new_extent[4] = ceiling(points_extent[4]/res)*res
    
    # create a template, "r", for the raster
    r <- raster::raster(crs = LVIS_footprints@proj4string)
    raster::extent(r) = new_extent
    res(r) = c(res, res)
    
    for(metric in metrics){
      rast <- raster::rasterize(x = LVIS_footprints,
                                y = r,
                                field = metric,
                                fun=function(x,...){mean(x, na.rm = T)})
      if(metric == "RH98" & PLOT == TRUE){
        plot(site)
        plot(rast, add = T)
        plot(site, add = T)
        writeRaster(rast,
                    filename = paste0(output_dir, '2_LVIS_rasters_by_plot/', site_id, '.tif'),
                    format="GTiff",
                    overwrite=T)
      }
      
      extracted = raster::extract(x = rast, y = site, weights=T)[[1]]
      if(!is.null(extracted)){
        extracted = extracted %>%
          as.data.frame() %>% filter(!is.na(value)) %>%
          mutate(weighted_value = value * weight,
                 weighted_sq_value = value^2 * weight)
        
        # calculate area-weighted average
        coverage = sum(extracted$weight)
        weighted_average = sum(extracted$weighted_value)/coverage
        # weighted_rms = sqrt(sum(extracted$weighted_sq_value)/coverage)
        
        if(metric=="ZG"){
          # only want to record these once per site
          summary_data[1, "plotID"] = site@data$plotID
          summary_data[1, "sbpltID"] = site@data$sbpltID
          summary_data[1, "sensor"] = site_data$sensor[1]
          summary_data[1, "LVIS_coverage"] = coverage
        }else if(metric == "RH98"){
          summary_data[1, "RH_98_plot_max"] = max(extracted$value, na.rm = T)
        }
        summary_data[1, metric] = weighted_average
        # summary_data[1, paste0(metric, '_RMS')] = weighted_rms
        
      }
    }
  }
  output_df = rbind(output_df, summary_data)
}

write.csv(output_df, paste0(output_dir, '2a_LVIS_summarized_by_NEON_subplot.csv'))



# profile = tidyr::gather(subset, key = "RH", value = "value",
#                         RH10, RH15, RH20, RH25, RH30, RH35, RH40,
#                         RH45, RH50, RH55, RH60, RH65, RH70, RH75,
#                         RH80, RH85, RH90, RH95, RH96, RH97, RH98, RH99, RH100) %>%
#   mutate(RH = str_remove(RH, "RH")) %>%
#   mutate(RH = ifelse(RH == "ZG", 0, RH)) %>%
#   mutate(RH = as.numeric(RH))   %>%
#   group_by(SHOTNUMBER) %>%
#   mutate(rel_value = value/max(value)) %>%
#   ungroup()
# 
# profile = filter(profile, value < 1000)
# 
# ggplot(profile,
#        alpha = 0.5) +
#   # facet_grid(cols = vars(BA)) + 
#   geom_line(aes(x = rel_value, y = RH, group = SHOTNUMBER, 
#                 # color = scientificName,
#                 color = plotID,
#                 alpha = BA))


