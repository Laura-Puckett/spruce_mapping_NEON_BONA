library(stringr); library(ggplot2)

LVIS_overlap_path = './analyses/BONA_species_mapping/LVIS/1_facility_2019_near_field.csv'
output_path = './analyses/BONA_species_mapping/LVIS/2_profile_data.csv'
veg = readRDS('./datasets/biomass/NEON/BONA/veg.rds')


#################################################################
### sort by species
#################################################################
LVIS_overlap = read.csv(LVIS_overlap_path)
plot_species_summary = veg %>% 
  group_by(plotID, scientificName) %>%
  mutate(BA = 3.14*(stemDiameter/2)^2) %>%
  summarize(BA = sum(BA, na.rm=T)) %>%
  ungroup()

plot_summary = veg %>% 
  group_by(plotID) %>%
  mutate(BA = 3.14*(stemDiameter/2)^2) %>%
  summarize(plot_BA = sum(BA, na.rm = T))

sp_summary = full_join(plot_species_summary, plot_summary) %>%
  mutate(sp_prop = BA/plot_BA*100)

sp_summary = sp_summary %>% filter(sp_prop>90)

NEON_subset = NEON_sites[NEON_sites$plotID %in% sp_summary$plotID,]

#####################################################################
LVIS_overlap = read.csv('./lidar/LVIS/LVIS_and_NEON_BONA/1_facility_2019_near_field.csv')

LVIS_overlap =  LVIS_overlap[LVIS_overlap$plotID %in% sp_summary$plotID,]

subset = LVIS_overlap %>% filter(RH100>2) # filter for trees

# # get average per site
# subset = subset %>%
#   group_by(plotID) %>%
#   summarise_if(is.numeric, mean, na.rm = TRUE)

subset = subset %>% left_join(sp_summary, by = "plotID")

profile = tidyr::gather(subset, key = "RH", value = "value",
                        RH10, RH15, RH20, RH25, RH30, RH35, RH40,
                        RH45, RH50, RH55, RH60, RH65, RH70, RH75,
                        RH80, RH85, RH90, RH95, RH96, RH97, RH98, RH99, RH100) %>%
  mutate(RH = str_remove(RH, "RH")) %>%
  mutate(RH = ifelse(RH == "ZG", 0, RH)) %>%
  mutate(RH = as.numeric(RH))   %>%
  group_by(SHOTNUMBER) %>%
  mutate(rel_value = value/max(value)) %>%
  ungroup()

profile = filter(profile, value < 1000)

ggplot(profile,
       alpha = 0.5) +
  # facet_grid(cols = vars(BA)) + 
  geom_line(aes(x = rel_value, y = RH, group = SHOTNUMBER, 
                # color = scientificName,
                color = plotID,
                alpha = BA))

write.csv(profile, output_path)

