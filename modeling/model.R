setwd('C:/Users/lkp58/Documents/research/')

library(rgdal); library(dplyr); library(raster); library(ggplot2); library(tidyr)
library(plotly)

LVIS_path = './analyses/BONA_species_mapping/LVIS/output/2a_LVIS_summarized_by_NEON_subplot.csv'
AVIRIS_path = './analyses/BONA_species_mapping/AVIRIS/output/1a_AVIRIS_by_NEON_subplot.csv'
veg_path = './analyses/BONA_species_mapping/NEON/data/veg_summary.rds'
subplots_path = './analyses/BONA_species_mapping/NEON/data/subplots/subplots.shp'

options(stringsAsFactors = F)
#################################################################
#### Vegetation Data
#################################################################

sites = readOGR(subplots_path)

veg = readRDS(veg_path)
# aggregate 4x100m2 distributed site groups to a single 400m2 site 
# inner join with the subplot data and filter for size will remove any merged 2x400m sites
merged_subplots_veg = veg %>% 
  group_by(plotID) %>%
  summarize(n_stems = sum(n_stems),
            BA = sum(BA),
            BA_PIMA = sum(BA_PIMA),
            BA_PIGL = sum(BA_PIGL),
            prop_PIMA = mean(prop_PIMA),
            prop_PIGL = mean(prop_PIGL),
            prop_spruce = mean(prop_spruce)) %>%
  mutate(subplotID = "merged")

# append merged version to single plot version and filter
veg_for_model = rbind(veg, merged_subplots_veg) %>%
  inner_join(sites@data %>%
               rename(subplotID = sbpltID)) %>%
  filter(sbpltSz == 400) %>%
  select(plotID, subplotID, n_stems, prop_spruce)


#################################################################
#### LVIS Data
#################################################################

LVIS = read.csv(LVIS_path) %>%
  select(-X) %>%
  rename(subplotID = sbpltID) %>% 
  filter(LVIS_coverage > 0.3) %>%
  mutate(d5 = RH10, # delta height, centered on 5%
         d15 = RH20-RH10,
         d25 = RH30-RH20,
         d35 = RH40-RH30,
         d45 = RH50-RH40,
         d55 = RH60-RH50,
         d65 = RH70-RH60,
         d75 = RH80-RH70,
         d85 = RH90-RH80,
         d95 = RH98-RH90,
         dQ1 = RH25,
         dQ2 = RH50-RH25,
         dQ3 = RH75-RH50,
         dQ4 = RH98-RH75,
         rd5 = d5/RH98,
         rd15 = d15/RH98,
         rd25 = d25/RH98,
         rd35 = d35/RH98,
         rd45 = d45/RH98,
         rd55 = d55/RH98,
         rd65 = d65/RH98,
         rd75 = d75/RH98,
         rd85 = d85/RH98,
         rd95 = d95/RH98) %>%
  group_by(plotID, subplotID) %>% # for some reason, dplyr was getting sd for all rows without grouping first
  mutate(sd_d = sd(c(d5, d15, d25, d35, d45, d55, d65, d75, d85, d95))) %>%
  ungroup()

#### plot
LVIS_for_plot = pivot_longer(LVIS %>% select(subplotID, plotID,
                                             RH10, RH15, RH20, RH25, RH30, RH35, RH40, RH45,
                                             RH50, RH55, RH60, RH65, RH70, RH75, RH80, RH85,
                                             RH90, RH95, RH96, RH97, RH98, RH99, RH100),
                             c(RH10, RH15, RH20, RH25, RH30, RH35, RH40, RH45,
                               RH50, RH55, RH60, RH65, RH70, RH75, RH80, RH85,
                               RH90, RH95, RH96, RH97, RH98, RH99, RH100),
             names_to = "RH", values_to = "value") %>%
  mutate(RH = stringr::str_remove(RH, 'RH'),
         id = paste0(plotID, '_', subplotID),
         RH = as.numeric(RH))

# just grab a few lines so the plots isn't as cluttered
LVIS_for_plot = LVIS_for_plot %>% 
  filter(id %in% unique(LVIS_for_plot$id)[1:5])

plot_ly() %>% layout(
  yaxis=list(title='Relative Height (m)', tickfont = list(size = 15), titlefont = list(size = 20)),
  xaxis = list(title = "Cumulative Return Energy (%)", tickfont = list(size = 15), titlefont = list(size = 20)),
  legend = list(x = 0.6, y = 0.9, font = list(size = 20)),
  title = paste0("")) %>%
  add_lines(x = LVIS_for_plot$RH,
            y = LVIS_for_plot$value,
            split = LVIS_for_plot$id,
            color = LVIS_for_plot$id,
            line=list(color="black", width=2),
            opacity=.8,
            showlegend=F)



max_location_of_waveform = pivot_longer(LVIS %>% select(subplotID, plotID, rd5, rd15, rd25, rd35, rd45, rd55, rd65, rd75, rd85, rd95),
                                        c(rd5, rd15, rd25, rd35, rd45, rd55, rd65, rd75, rd85, rd95),
                                        names_to = "key", values_to = "value") %>%
  mutate(key = stringr::str_remove(key, 'rd')) %>%
  group_by(subplotID, plotID) %>%
  mutate(waveform_max = max(value),
         percentile_of_waveform_max = ifelse(value == waveform_max, key, NA )) %>%
  ungroup() %>%
  filter(is.na(percentile_of_waveform_max) ==F) %>%
  select(plotID, subplotID, waveform_max, percentile_of_waveform_max)

LVIS_for_model = inner_join(LVIS, max_location_of_waveform) %>%
  select(plotID,
         subplotID,
         #ZG,
         RH25, RH50, RH75, RH98, RH_98_plot_max,
         dQ1, dQ2, dQ3, dQ4,
         sd_d,
         waveform_max,
         percentile_of_waveform_max)

#################################################################
#### AVIRIS Data
#################################################################

AVIRIS = read.csv(AVIRIS_path) %>%
  select(-X, -siteID, -subtype, -plotTyp, -sbpltSz) %>%
  rename(subplotID = sbpltID)

AVIRIS_indices = AVIRIS %>%
  group_by(plotID, subplotID) %>%
  transmute(BLUE = mean(c(X451.99_nm, X457_nm, X462.01_nm, X467.02_nm, X472.02_nm, X477.03_nm, X482.04_nm, X487.05_nm, 
                          X492.06_nm, X497.07_nm, X502.08_nm, X507.09_nm)),
            RED = mean(c(X642.32_nm, X647.33_nm, X652.34_nm, X657.35_nm, X662.35_nm, X667.36_nm)),
            NIR = mean(c(X852.68_nm, X857.69_nm, X862.7_nm, X867.71_nm, X872.72_nm, X877.73_nm)),
            NDVI = (NIR - RED)/(NIR + RED),
            EVI = 2.5 * ( NIR - RED ) / ( NIR + 6.0 * RED - 7.5 * BLUE+ 1.0 )) %>%
  ungroup()

AVIRIS_bands_only = AVIRIS %>% select(-plotID, -subplotID)
AVIRIS_bands_subset = AVIRIS_bands_only[,seq(1, ncol(AVIRIS_bands_only), 5), ]

# a better method would be to take the average of bands to reduce noise and dimensionality...revisit this later
# AVIRIS_bands = colnames(AVIRIS_bands_only) %>%
# stringr::str_remove("X") %>%
# stringr::str_remove("_nm") %>%
# as.numeric()
# AVIRIS_summarized_bands[,paste0("X",mean(AVIRIS_bands[1:5]),"_nm_mean")] = mean(AVIRIS_bands_only[,1:5])
# %>% pivot_longer(c(-plotID, -subplotID), "Nanometers","Reflectance" )

AVIRIS_for_model = cbind(AVIRIS_indices, AVIRIS_bands_subset)

#################################################################
#### Select sites
#################################################################
data = inner_join(LVIS_for_model, AVIRIS_for_model) %>%
  inner_join(veg_for_model)

ggplot(data) +
  geom_histogram(aes(x = prop_spruce, fill = plotID), bins = 15) +
  theme_bw() +
  xlab('spruce proportion by basal area')

colnames(data)

data_for_model = data %>%
  select(-plotID, -subplotID, -n_stems)

#### Fit model for proportial black spruce BA
#################################################################
library(caret); library(dplyr); library(tidyr); library(stringr);

set.seed(998)
# trainIndex <- createDataPartition(data_for_model$prop_spruce, p = .75, list = FALSE)
# training <- data_for_model[ trainIndex,]
# testing  <- data_for_model[-trainIndex,]

fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  repeats = 30)

set.seed(825)

rfFit1 <- train(prop_spruce ~ ., data = data_for_model, 
                method = "rf", 
                trControl = fitControl,
                verbose = FALSE)
rfFit1
varImp(rfFit1, scale = TRUE)


