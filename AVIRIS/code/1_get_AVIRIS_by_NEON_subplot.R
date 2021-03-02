#*********************************************************************************#
# -----0: Packages
# ********************************************************************************#
list.of.packages <- c("rgeos","dplyr", "rgdal", "sp", "raster", "stringr", "tidyr", "ggplot2")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages, repos='http://cloud.r-project.org')

library(rgeos)
library(raster)
library(rgdal)
library(dplyr)
library(sp)
library(stringr)
library(tidyr)
library(ggplot2)
library(plotly)
library(orca)

#*********************************************************************************#
# -----1: Inputs
# ********************************************************************************#
setwd('C:/Users/lkp58/Documents/research/')
output_dir='./analyses/BONA_species_mapping/AVIRIS/'

AVIRIS_path = '//nau.froot.nau.edu/cirrus/scratch/lkp58/datasets/AVIRIS/sites/BONA/2_gdalwarp_to_tiff/ang20180723t193341_corr_v2r2_img.tif'
subplots_path = './analyses/BONA_species_mapping/NEON/data/subplots/subplots.shp'

PLOT = FALSE
#*********************************************************************************#
# -----2: Read Files
# ********************************************************************************#

img = brick(AVIRIS_path)
sites = readOGR(subplots_path)


output_df = NULL
for (i in 1:nrow(sites)){
  # print(i)
  site = sites[i,]
  overlap = raster::extract(img, site, cellnumbers = F, weights = T, normalizeWeights = T)[[1]]
  values = overlap[,colnames(overlap) != "weight"]
  weights = overlap[,"weight"]
  
  # check if the site overlapped with a pixel, only continue if so
  
  if(all(is.na(values))==F){
    
    bands = colnames(values)
    bands = stringr::str_remove(bands, ".Nanometers")
    bands = stringr::str_remove(bands, "X")
    bands = as.numeric(bands) %>% round(digits = 2)
    bands
    
    # multiply values by weights and sum for each band
    weighted_mean = crossprod(values, weights)
    rownames(weighted_mean) = NULL
    weighted_mean = cbind("Nanometers" = bands, "reflectance" = weighted_mean) %>%
      as.data.frame()
    colnames(weighted_mean) = c("Nanometers", "reflectance")
    
    weighted_mean = weighted_mean %>%
      filter((Nanometers >= 1330 & Nanometers <= 1435) == F,
             (Nanometers >= 1814 & Nanometers <= 1960) == F,
             ((Nanometers >= 2480) == F))
    
    # swich rows and columns
    weighted_mean_tidy = weighted_mean %>%
      select(-Nanometers) %>% (t) %>% as.data.frame()
    colnames(weighted_mean_tidy) = paste0("X", weighted_mean$Nanometers, "_nm")
    rownames(weighted_mean_tidy) = NULL
    
    output_df_row = cbind(site@data, weighted_mean_tidy)
    output_df = rbind(output_df, output_df_row)
    
    if(PLOT == T){
      # cropped = crop(img, rgeos::gBuffer(site, width = 5))
      # plot(cropped[[65]])
      # plot(site, add = T)
      
      values_for_plot = cbind("Nanometers" = bands, t(values))
      rownames(values_for_plot) = NULL
      
      values_for_plot = values_for_plot %>%
        as.data.frame() %>%
        gather("pixel", "reflectance", -Nanometers) %>%
        mutate(reflectance = as.numeric(reflectance)) %>%
        group_by(pixel)
      
      values_for_plot = values_for_plot %>%
        filter((Nanometers >= 1330 & Nanometers <= 1435) == F,
               (Nanometers >= 1809 & Nanometers <= 1960) == F,
               ((Nanometers >= 2480) == F))
      
      # ggplot(values_for_plot, aes(x = Nanometers, y = reflectance)) +
      #   geom_point(aes(color = pixel)) +
      #   ggtitle(paste0("AVIRIS data overlapping ", site@data$plotID, " subplot ", site@data$sbpltID)) +
      #   theme_bw()
      
   #### Spectral figure using plotly ####
      
      # fig <- plot_ly() %>% layout(
      #   yaxis=list(title='Reflectance', tickfont = list(size = 15), titlefont = list(size = 20)),
      #   xaxis = list(title = "Wavelength [nm]", tickfont = list(size = 15), titlefont = list(size = 20)),
      #   legend = list(x = 0.6, y = 0.9, font = list(size = 20)),
      #   title = paste0("AVIRIS data overlapping ", site@data$plotID, " subplot ", site@data$sbpltID)) %>%
      # # dummy line so I can splot a single legend entry for "spectra by pixel"
      # add_lines(x = weighted_mean$Nanometers,
      #                  y = weighted_mean$reflectance,
      #                  name = "By pixel",
      #                  line=list(color="crimson", width=2),
      #                  opacity = 0.8) %>%
      # # the actual spectra by pixel data
      #   add_lines(x = values_for_plot$Nanometers,
      #                 y = values_for_plot$reflectance,
      #                 split = values_for_plot$pixel,
      #                 line=list(color="crimson", width=2),
      #                 opacity=.8,
      #                 showlegend=F) %>%
      # # the weight average (will plot over that dummy line)
      #  add_lines(x = weighted_mean$Nanometers,
      #                   y = weighted_mean$reflectance,
      #                   name = "Weighted average",
      #                   line=list(color="black", width=4),
      #                   opacity=1,
      #                   showlegend=TRUE)
      # 
      # orca(p = fig, file = paste0(output_dir, 'figures/1a_AVIRIS_spectra_for_NEON', site@data$plotID, '_', site@data$sbpltID, '.png'))

      fig <- plot_ly() %>% layout(
        yaxis=list(title='Reflectance', tickfont = list(size = 15), titlefont = list(size = 20)),
        xaxis = list(title = "Wavelength [nm]", tickfont = list(size = 15), titlefont = list(size = 20)),
        legend = list(x = 0.6, y = 0.9, font = list(size = 20)),
        title = paste0("AVIRIS spectra by pixel" )) %>%
        add_lines(x = values_for_plot$Nanometers,
                      y = values_for_plot$reflectance,
                      split = values_for_plot$pixel,
                      line=list(color="crimson", width=2),
                      opacity=.8,
                      showlegend=F)
      fig
      
    }
  }
}


write.csv(output_df, paste0(output_dir, 'output/1a_AVIRIS_by_NEON_subplot.csv'))
