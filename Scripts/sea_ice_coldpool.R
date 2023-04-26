library("sf")
library("dplyr")
library("akgfmaps")
library("ggplot2")

# API
#ice_poly <- st_read("https://nsidc.org/api/mapservices/NSIDC/ows?service=WFS&version=1.0.0&request=GetFeature&typeName=NSIDC:g02135_polyline_n") # &filter=time@filter_from=1982-03
# filter to march of every year from ice_ext_2$timestamp

# Data downloaded by FTP 2023-04-25: ftp://sidads.colorado.edu/pub/DATASETS/NOAA/G02135/
# More information: https://nsidc.org/data/g02135/versions/3
ice_ras <- read.csv("Data/Sea_ice_data/NSIDC.csv")
#TODO: replace with same product for full range of latitudes and redo (these data are cut off halfway through EBS)

# filter data by region, time, sea ice cover (only available from 1988 onward for March in this version and need lower latitudes)
#Latitude <= 73,
#Longitude <= -156,
#Longitude >= -175.3,
ice_df <- ice_ras %>% filter(Month == 3, 
                             Year > 1981) %>%
                      na.omit()

# Minimum latitude by year
# TODO: try Seaice >= 0.1
ice_ext <- ice_df %>% filter(Seaice >= 0.15) %>% group_by(Year) %>% summarise(Extent = min(Latitude))
colnames(ice_ext) <- tolower(colnames(ice_ext))

# compare to coldpool extent
coldpool <- coldpool:::cold_pool_index
colnames(coldpool) <- tolower(colnames(coldpool))

df <- left_join(ice_ext, coldpool) %>% na.omit()

m_lte2 <- lm(log(area_lte2_km2) ~ extent, df)
m_lte1 <- lm(area_lte1_km2 ~ log(extent), df)
m_lte0 <- lm(area_lte0_km2 ~ log(extent), df)
m_lteminus1 <- lm(area_lteminus1_km2 ~ log(extent), df)

print(summary(m_lte2))
print(summary(m_lte1))
print(summary(m_lte0))
print(summary(m_lteminus1))
# r-square degrades with decreasing temperature threshold

#best relationship plotted
plot(df$extent, log(df$area_lte2_km2))
abline(m_lte2)

# TODO: separate by outer/inner/middle domain boxes from stratum shapefile and calculate extent specific to each?
ice <- ice_df %>% filter(Seaice >= 0.15)
ice_sf <- st_as_sf(ice, coords = c("Longitude", "Latitude"))
st_crs(ice_sf) <- 4326

ebs <- akgfmaps::get_base_layers("ebs", "EPSG:4326")

# TODO: use bathymetry to define inner/outer/middle and then do spatial join between areas and ice_sf
# Note that there may be complications as bathymetry includes geometries for islands and no inner boundary for land
ggplot() +
  geom_sf(data = ebs$akland) +
  geom_sf(data = ebs$bathymetry) +
  geom_sf(data = ice_sf, color = "red") +
  geom_sf(data = ebs$graticule, color = "grey70", alpha = 0.5) +
  coord_sf(xlim = ebs$plot.boundary$x,
           ylim = ebs$plot.boundary$y) +
  scale_x_continuous(name = "Longitude",
                     breaks = ebs$lon.breaks) +
  scale_y_continuous(name = "Latitude",
                     breaks = ebs$lat.breaks) +
  theme_bw()

#st_filter(ice_sf, ebs$bathymetry) #could use this if make bathymetry into polygon with 0 depth/land
