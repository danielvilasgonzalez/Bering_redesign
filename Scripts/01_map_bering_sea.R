####################################################################
####################################################################
##
##    Create a Bering Sea map figure
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu)
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c('rasterVis','rgeos','scales','rnaturalearth','grid','ggplot2')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#install akgfmaps to extract shapefile of Alaska
if (!('akgfmaps' %in% installed.packages())) {
  devtools::install_github('afsc-gap-products/akgfmaps')};library(akgfmaps)

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#set working directory
mydir<-'E:/UW/Adapting Monitoring to a Changing Seascape/Resources/'
setwd(mydir)

# Define plot extent (through trial end error) units km
panel_extent <- data.frame(x = c(-1716559.21, -77636.05), #x = c(-1326559.21, -87636.05),
                           y = c(483099.5, 2194909.7)) #y = c(533099.5, 1894909.7))

#Alaska layer
ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")
ak_sppoly<-as(ebs_layers$akland, 'Spatial')

#bathy from GEBCO (https://download.gebco.net/)
#ak_bathy0<-raster("C:/Users/danie/Desktop/UW/GEBCO_21_Sep_2022_895fcb2e2466/GEBCO_21_Sep_2022_895fcb2e2466/gebco_2022_n74.2676_s50.1416_w175.2539_e180.0.asc")
#ak_bathy1<-raster("C:/Users/danie/Desktop/UW/GEBCO_21_Sep_2022_895fcb2e2466/GEBCO_21_Sep_2022_895fcb2e2466/gebco_2022_n69.8291_s50.0977_w-179.9121_e-167.0.asc")
#ak_bathy2<-raster("C:/Users/danie/Desktop/UW/GEBCO_21_Sep_2022_895fcb2e2466/GEBCO_21_Sep_2022_895fcb2e2466/gebco_2022_n69.8291_s50.0977_w-167.0_e-154.0.asc")
#ak_bathy<-raster::merge(ak_bathy1,ak_bathy2)
#ak_bathy_2<-projectRaster(ak_bathy,
#                          crs='+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
#
#save raster with the new projection
#writeRaster(ak_bathy_2,
#            'E:/UW/Adapting Monitoring to a Changing Seascape/Resources/ak_bathy_NAD83.tiff',overwrite=FALSE)
ak_bathy_2<-raster('./Bathymetry/ak_bathy_NAD83.tiff')
#plot(ak_bathy_2)

#shapefile EBS
ebs_sh<-rgdal::readOGR(dsn='./Shapefiles/',layer = 'EBSshelfThorson')

#shapefile NBS
nbs_sh<-rgdal::readOGR(dsn='./Shapefiles/',layer = 'NBSThorson')

#shapefile NBS
slo_sh<-rgdal::readOGR(dsn='./Shapefiles/',layer = 'EBSslopeThorson')
proj4string(slo_sh) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
slo_sh<-spTransform(slo_sh,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))

#merge shapefiles
bs_sh<-union(ebs_sh,nbs_sh)
bs_sh<-union(bs_sh,slo_sh)

#shapefile EEZ
eez_sh<-rgdal::readOGR(dsn='./EEZ',layer = 'EEZ_Land_v3_202030')

#clip EEZ
bbox = c(latN = 70, latS = 50, lonW = -200, lonE = -150)
pol <- extent(bbox[3],bbox[4], bbox[2],bbox[1])
eez_sh1<-crop(eez_sh, pol)
bbox = c(latN = 70, latS = 50, lonW = 160, lonE = 180)
pol <- extent(bbox[3],bbox[4], bbox[2],bbox[1])
eez_sh2<-crop(eez_sh, pol)

#change CRS projection of EEZ files
proj4string(eez_sh1) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
eez_sh1<-spTransform(eez_sh1,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))
proj4string(eez_sh2) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
eez_sh2<-spTransform(eez_sh2,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))

## crop and mask
ak_bathy_3 <- crop(ak_bathy_2, extent(bs_sh))
ak_bathy_4 <- mask(ak_bathy_3, bs_sh)

#positive values equal to zero and convert negative to positive values
ak_bathy_4[ak_bathy_4>0]<-0 
ak_bathy_4<--ak_bathy_4 

#extract station EBS bottom trawl
st_EBS<-read.csv('./Stations/ebs_nbs_temperature_full_area.csv')

#filter 2019 stations, an example year where EBS and NBS surveys were carried out
st_EBS<-subset(st_EBS,year==2019 ) #& survey_definition_id ==98

#convert to spatial object
coordinates(st_EBS)<-c('longitude','latitude')
proj4string(st_EBS) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 

#change CRS
st_EBS<-spTransform(st_EBS,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))

#convert to dataframe
st_EBS1<-as.data.frame(st_EBS)
st_EBS1$survey_definition_id<-as.factor(st_EBS1$survey_definition_id)

#remove some polygons of EEZ object
eez_sh11 <- eez_sh1[eez_sh1$AREA_KM2 == 5193061,] 
eez_sh22 <- eez_sh2[eez_sh2$AREA_KM2 == 5193061,]  #"5193061"  "24614858" "8521"    

#join both polygons without inner line
eez_sh3<-aggregate(rbind(eez_sh11,eez_sh22),dissolve=T)
eez_sh33<-rgeos::gUnaryUnion(eez_sh3)

#zoomin plot
zoomin<-
  gplot(ak_bathy_4) +
    geom_tile(aes(fill=value))+
    geom_point(data=st_EBS1,aes(x=longitude,y=latitude),shape=4,size=1)+
    geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
    geom_polygon(data=eez_sh33,aes(x=long,y=lat,group=group),fill=NA,color='grey40')+
    scale_x_continuous(expand = c(0,0),position = 'bottom',
                       breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
    geom_polygon(data=nbs_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
    geom_polygon(data=ebs_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
    geom_polygon(data=slo_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
    coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
             xlim = panel_extent$x,
             ylim = panel_extent$y,
             label_axes = "-NE-")+
    scale_fill_gradient2(low = '#c1f2fe',
                        high = '#007c9b',
                       limits=c(0,200),oob = scales::squish,breaks=c(0,50,100,200),
                       labels=c('0','50','100',paste0('200 - ',round(maxValue(ak_bathy_4)))),
                        na.value = 'white',
                        name='depth (m)',
                        guide = guide_colorbar(frame.colour = 'black',ticks.colour = 'black'))+
    theme(aspect.ratio = 1,panel.grid.major = element_line(color = rgb(0, 0, 0,20, maxColorValue = 285), linetype = 'dashed', linewidth =  0.5),
          panel.background = element_rect(fill = NA),panel.ontop = TRUE,text = element_text(size=10),
          legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(25, 'points'),
          legend.key.width= unit(25, 'points'),axis.title = element_blank(),legend.position = c(0.12,0.47),
          panel.border = element_rect(fill = NA, colour = 'black'),legend.key = element_rect(color="black"),
          axis.text = element_text(color='black'),legend.spacing.y = unit(10, 'points'),
          axis.text.y.right = element_text(hjust= 0.1 ,margin = margin(0,7,0,-25, unit = 'points'),color='black'),
          axis.text.x = element_text(vjust = 6, margin = margin(-7,0,7,0, unit = 'points'),color='black'),
          axis.ticks.length = unit(-5,"points"))+
    annotate("text", x = -256559, y = 1354909, label = "Alaska",parse=TRUE,size=7)+
    annotate("text", x = -1376559, y = 2049090, label = "Russia",parse=TRUE,size=7)+
    scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())+
    annotate("text", x = -1376559, y = 744900, label = "italic('Bering Sea')",parse=TRUE,size=9)

#extract countries and lake for the zoomout plot
count<-ne_countries(scale = 'large',
                    country = c('united states of america',"russia","canada",'mexico'),
                    returnclass = c("sp", "sf"))
lakes <- ne_download(category = "physical", type = "lakes", returnclass = "sp", 
                     scale = 'small')

#reproject shapefiles
count83<-spTransform(count,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))
lakes83<-spTransform(lakes,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))

#zoomout plot
zoomout<-
  ggplot() + 
    geom_polygon(data=count83,aes(x=long,y=lat,group=group),color='black',fill='grey80') + 
    geom_polygon(data=lakes83,aes(x=long,y=lat,group=group),color='black',fill='white') + 
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(expand = c(0,0),
                       breaks = c(160,170,-180,-170,-160,-150,-140,-130))+
    coord_sf(crs ='+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +no_defs',
             xlim=c(-3500000, 3500000),ylim = c(-2000000, 3000000))+
    theme(panel.grid.major = element_line(color = rgb(0, 0, 0,20, maxColorValue = 285), linetype = 'dashed', linewidth =  0.5),
          panel.ontop = TRUE,text = element_text(size=10),
          legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(20, 'points'),
          legend.key.width= unit(20, 'points'),axis.title = element_blank(),legend.position = c(0.12,0.155),
          panel.border = element_rect(fill = NA, colour = 'black'),
          axis.text = element_blank(),axis.ticks.length = unit(-3, "points"),
          panel.background = element_rect(fill = NA),plot.margin=grid::unit(c(0,0,0,0), "mm"))+
    geom_rect(aes(xmin = panel_extent$x[1], xmax = panel_extent$x[2], ymin = panel_extent$y[1], ymax = panel_extent$y[2]),
          colour = '#800000', linetype='solid', fill = NA,size=0.7) +
    annotate("text", x = -2300000, y = 2600000, label = "Russia",parse=TRUE,size=4)+
    annotate("text", x = 0, y = -1000000, label = "italic('Pacific Ocean')",parse=TRUE,size=4)+
    annotate("text", x = +2000000, y = 1400000, label = "Canada",parse=TRUE,size=4)+
    annotate("text", x = +3000000, y = 0, label = "USA",parse=TRUE,size=4)+
    annotation_scale()


#save plot
tiff(filename = 'E:/UW/Adapting Monitoring to a Changing Seascape/Figures/map_bering_sea.tif',res = 220,width = 1500,height = 1900)
grid.newpage()
vp_b <- viewport(width = 1, height = 1, x = 0.5, y = 0.5)  # the larger map
vp_a <- viewport(width = 0.4, height = 0.3, x = 0.211, y = 0.78)  # the inset in upper left
print(zoomin , vp = vp_b)
print(zoomout , vp = vp_a)
dev.off()
