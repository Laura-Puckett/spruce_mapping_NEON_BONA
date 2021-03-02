veg_shp_path = './analyses/BONA_species_mapping/NEON/data/veg_shp/NEON_BONA_veg.shp'
subplots_path = './analyses/BONA_species_mapping/NEON/data/subplots/subplots.shp'
veg_shp = readOGR(veg_shp_path)
subplots = readOGR(subplots_path)

site = subplots[55,]
plot(site)

site_veg = veg_shp[veg_shp@data$plotID == site@data$plotID,]
plot(veg_shp, color = veg_shp@data$taxonID, add = T)

veg_coords = site_veg@coords
veg_for_ggplot = cbind(veg_coords, site_veg@data)

ggplot(veg_for_ggplot, aes(x = coords.x1, y = coords.x2)) +
  geom_point(aes(color = scntfcN, size = stmDmtr)) + 
  theme_bw()
