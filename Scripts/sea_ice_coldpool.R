library("sf")
library("dplyr")

# API
#ice_poly <- st_read("https://nsidc.org/api/mapservices/NSIDC/ows?service=WFS&version=1.0.0&request=GetFeature&typeName=NSIDC:g02135_polyline_n") # &filter=time@filter_from=1982-03
# filter to march of every year from ice_ext_2$timestamp

# Data downloaded by FTP 2023-04-25: ftp://sidads.colorado.edu/pub/DATASETS/NOAA/G02135/
# More information: https://nsidc.org/data/g02135/versions/3
ice_ras <- read.csv("Data/Sea_ice_data/NSIDC.csv")
#TODO: replace with same product for full range of latitudes and redo

# filter data by region, time, sea ice cover (only available from 1988 onward for March in this version and need lower latitudes)
#Latitude <= 73,
#Longitude <= -156,
#Longitude >= -175.3,
ice_df <- ice_ras %>% filter(Month == 3, 
                             Year > 1981) %>%
                      na.omit()

# Minimum latitude by year
ice_ext <- ice_df %>% filter(Seaice >= 0.15) %>% group_by(Year) %>% summarise(Extent = min(Latitude))
colnames(ice_ext) <- tolower(colnames(ice_ext))

# compare to coldpool extent
coldpool <- coldpool:::cold_pool_index
colnames(coldpool) <- tolower(colnames(coldpool))

df <- left_join(ice_ext, coldpool)

m_lte2 <- lm(log(area_lte2_km2) ~ extent, df)
m_lte1 <- lm(log(area_lte1_km2) ~ extent, df)
m_lte0 <- lm(log(area_lte0_km2) ~ extent, df)
m_lteminus1 <- lm(log(area_lteminus1_km2) ~ extent, df)

print(summary(m_lte2))
print(summary(m_lte1))
print(summary(m_lte0))
print(summary(m_lteminus1))
# r-square degrades with decreasing temperature threshold

#best relationship plotted
plot(df$extent, log(df$area_lte2_km2))
abline(m_lte2)

# TODO: separate by outer/inner/middle domain boxes from stratum shapefile and calculate extent specific to each?