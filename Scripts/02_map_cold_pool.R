####################################################################
####################################################################
##
##    Create a cold pool extension map
##    from a) coldpool package and b) ROMS Bering Sea
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu)
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c('rasterVis','scales','rnaturalearth','cowplot','ggplot2','ncdf4','googledrive','rgeos','rgdal')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#install coldpool to extract SBT for the EBS
if (!('coldpool' %in% installed.packages())) {
   devtools::install_github('afsc-gap-products/coldpool')};library(coldpool)

#install akgfmaps to extract shapefile of Alaska
if (!('akgfmaps' %in% installed.packages())) {
  devtools::install_github('afsc-gap-products/akgfmaps')};library(akgfmaps)

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#set working directory
mydir<-'E:/UW/Adapting Monitoring to a Changing Seascape/'
setwd(mydir)

# Define plot extent (through trial end error) units km
panel_extent <- data.frame(x = c(-1856559.21, 55636.05), #x = c(-1326559.21, -87636.05),
                           y = c(353099.5, 1950909.7)) #y = c(533099.5, 1894909.7))

#Alaska layer
ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")
ak_sppoly<-as(ebs_layers$akland, 'Spatial')

#shapefile EBS
ebs_sh<-rgdal::readOGR(dsn='./Resources/Shapefiles/',layer = 'EBSshelfThorson')

####################################################
# SBT from coldpool package
####################################################

#sbt ebs
sbt_brick<-coldpool::ebs_bottom_temperature

#names
years<-sort(as.numeric(gsub(".*?([0-9]+).*", "\\1", names(sbt_brick))))

#list to store values
plot_list<-list()

#loop over years, coldpool only has data between 1982 and 2022
for (y in min(years):max(years)) {

  #y<-1982
  
  #print year to check progress
  cat(paste("    ----- year", y, "-----\n"))  
  
  #2020 empty, so all NAs
  if (y == 2020) {
    r<-sbt_brick[[grep(y-1,names(sbt_brick))]]
    r[r]<-NA
  } else {
    r<-sbt_brick[[grep(y,names(sbt_brick))]]
  }
  
  #plot
  p<-
   gplot(r) +
    geom_tile(aes(fill=value))+
    geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
    scale_x_continuous(expand = c(0,0),position = 'bottom',
                       breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
    geom_polygon(data=ebs_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
    coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
             xlim =c(-1716559.21,57636.05),
             ylim = c(403099.5, 1754909.7),
             label_axes = "-NE-")+
    scale_fill_binned(limits=c(-2,2),na.value = 'white',low = "#007c9b", high = "#c1f2fe",name='SBT (°C)',guide = guide_colorbar(frame.colour = 'black',ticks.colour = NA))+
    theme(aspect.ratio = 1,panel.grid.major = element_line(color = rgb(0, 0, 0,20, maxColorValue = 285), linetype = 'dashed', linewidth = 0.5),
        panel.background = element_rect(fill = NA),panel.ontop = TRUE,text = element_text(size=10),
        legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(10, 'points'),
        legend.key.width= unit(10, 'points'),axis.title = element_blank(),legend.position = c(0.15,0.35),
        panel.border = element_rect(fill = NA, colour = 'black'),legend.key = element_rect(color="black"),
        axis.text = element_text(color='black'),legend.spacing.y = unit(10, 'points'),
        axis.text.y.right = element_text(hjust= 0.1 ,margin = margin(0,0,0,-25, unit = 'points'),color='black'),
        axis.text.x = element_text(vjust = 6, margin = margin(-7,0,0,0, unit = 'points'),color='black'),
        axis.ticks.length = unit(-5,"points"))+
    scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())+
    annotate("text", x = -1406559, y = 1624909, label = y,size=8)

  #remove legend from plots except first (topleft) in the panel
  if (!(y %in% years[c(1,17,33)])) {
    p<-p+theme(legend.position = 'none')}
  
  #store plot on list
  plot_list[[paste0(y)]]<-p
  
}

#save plot 16 first years (4x4)
tiff(paste0('./Figures/coldpool_',years[1],'_',years[16],'.tiff'),height=2500,width = 2500,res = 250)
cowplot::plot_grid(plotlist = plot_list[1:16],nrow = 4,ncol = 4)
dev.off()

#save plot 16 next years (4x4)
tiff(paste0('./Figures/coldpool_',years[17],'_',years[32],'.tiff'),height=2500,width = 2500,res = 250)
cowplot::plot_grid(plotlist = plot_list[17:32],nrow = 4,ncol = 4)
dev.off()

#save plot last years (4x4)
tiff(paste0('./Figures/coldpool_',years[33],'_',years[length(years)],'.tiff'),height=2500,width = 2500,res = 250)
cowplot::plot_grid(plotlist = plot_list[33:length(min(years):max(years))],nrow = 4,ncol = 4)
dev.off()

####################################################
# SBT from Bering 10K ROMS
####################################################

#shapefile EBS
ebs_sh<-rgdal::readOGR(dsn='./Resources/Shapefiles/',layer = 'EBSshelfThorson')

#shapefile NBS
nbs_sh<-rgdal::readOGR(dsn='./Resources/Shapefiles/',layer = 'NBSThorson')

#shapefile NBS
slo_sh<-rgdal::readOGR(dsn='./Resources/Shapefiles/',layer = 'EBSslopeThorson')
proj4string(slo_sh) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
slo_sh<-spTransform(slo_sh,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))

#merge shapefiles
bs_sh<-union(ebs_sh,nbs_sh)
bs_sh<-union(bs_sh,slo_sh)

#get files from google drive and set up
files<-googledrive::drive_find()
1 #for dvilasg@uw.edu

#get id shared folder from google drive
id.roms.folder<-files[which(files$name=='Bering 10K ROMS'),'id']

#list of files and folder
files.1<-googledrive::drive_ls(id.roms.folder$id)
id.data<-files.1[which(files.1$name=='netcdf_historical'),'id']
nc_histfiles<-googledrive::drive_ls(id.data$id)
id.data<-files.1[which(files.1$name=='netcdf_forecast'),'id']
nc_forfiles<-googledrive::drive_ls(id.data$id)
files.hist<-sort(nc_histfiles$name)
files.for<-sort(nc_forfiles$name)

#open downloaded SBT weekly netcdf files from Bering 10K ROMS (https://data.pmel.noaa.gov/aclim/thredds/catalog/files.html)
# nc_histfiles<-list.files('./Resources/ACLIM2/Data/out/netcdf_historical/',full.names = TRUE)
# nc_forfiles<-list.files('./Resources/ACLIM2/Data/out/netcdf_forecast/',full.names = TRUE)

#create folder to store netcdf
dir.create('./Data/Bering 10K ROMS/',showWarnings = FALSE)

#list to store values
plot_list<-list()

#loop over years to incorporate values into the Bering Sea grid
for (y in min(years):max(years)) {
  
  #y<-1982
  
  #print year to check progress
  cat(paste("    ----- year", y, "-----\n"))  
  
  #open netcdf file for each year
  if (y %in% c(1980:1984)) {
    f<-files.hist[1]
    file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
  } else  if (y %in% c(1985:1989)) {
    f<-files.hist[2]
    file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
  } else  if (y %in% c(1990:1994)) {
    f<-files.hist[3]
    file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
  } else  if (y %in% c(1995:1999)) {
    f<-files.hist[4]
    file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
  } else if (y %in% c(2000:2004)) {
    f<-files.hist[5]
    file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
  } else if (y %in% c(2005:2009)) {
    f<-files.hist[6]
    file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
  } else if (y %in% c(2010:2014)) {
    f<-files.hist[7]
    file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
  } else if (y %in% c(2015:2019)) {
    f<-files.hist[8]
    file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
  } else if (y %in% c(2020)) {
    f<-files.hist[9]
    file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
  } else if (y %in% c(2021:2022)) {
    f<-files.for[1]
    file.id<-nc_forfiles[which(nc_forfiles$name==f),'id']} #for year>2020 have to select projection
 
  #download file into temp folder
  googledrive::drive_download(file=file.id$id,
                              path = paste0('./Data/Bering 10K ROMS/',f),
                              overwrite = TRUE)
  
  #open netcdf
  nc<-nc_open(paste0('./Data/Bering 10K ROMS/',f))
  
  #dimensions netcdf files
  #258 rows
  #182 cols
  #46956 cells
  #259 time steps
  
  #get variables
  #names(nc$var)
  
  #get latitude
  lats <- ncvar_get(nc,"lat_rho")
  
  #get longitude
  lons <- ncvar_get(nc,"lon_rho")
  
  #get SBT
  temp<-ncvar_get(nc,'temp')
  
  #get time
  t_axis<-ncvar_get(nc,"ocean_time")
  
  #convert time
  time_axis <- as.POSIXct(t_axis, origin = "1900-01-01", tz = "GMT") 
  
  #get weekly temp slices from specific year y
  nc_y<-ncvar_get(nc, "temp")[,,which(grepl(paste0(y,'-'),time_axis))]
  
  #get mean matrix for this year
  mean_nc<-apply(nc_y,c(1,2),mean,na.rm=TRUE)
  
  #create dataframe with lats, lons and mean year SBT
  df_nc<-data.frame(Lat=as.vector(lats),
                    Lon=as.vector(lons),
                    temp=as.vector(mean_nc))
  
  #longitude are in eastern. get SBT for the western hemisphere (Bering Sea). longitude greater than 180 degrees
  df_nc1<-df_nc[which(df_nc$Lon>=180),]
  
  #convert eastern longitude to western values (higher). longitude should be negative
  df_nc1$Lon<-df_nc1$Lon-180-180
  
  #remove rows with NAs
  df_nc2<-df_nc1[complete.cases(df_nc1),]
  
  #filter values from the grid 
  #df_nc3<-subset(df_nc2,Lat >= 54 & Lat <= 67 & Lon >= -179 & Lon <= -155)
  df_nc3<-df_nc2
  
  #create spatial object from df
  coordinates(df_nc3) <- ~ Lon + Lat
  proj4string(df_nc3) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
  df_nc4<-spTransform(df_nc3,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))
  #proj4string(df_nc3) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 

  # create a template raster
  r1 <- raster(ext=extent(df_nc4),res=c(15800,15800))
  r1<-rasterize(df_nc4, r1 )
  r1<-dropLayer(r1,'ID')
  #plot(r1)
  
  ## crop and mask
  r2 <- crop(r1, extent(bs_sh))
  r3<- mask(r2, bs_sh)
  
  p<-
    gplot(r3) +
    geom_tile(aes(fill=value))+
    geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
    scale_x_continuous(expand = c(0,0),position = 'bottom',
                       breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
    geom_polygon(data=ebs_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
    geom_polygon(data=nbs_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
    geom_polygon(data=slo_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
    coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
             xlim = panel_extent$x,
             ylim = panel_extent$y,
             label_axes = "-NE-")+
    scale_fill_binned(limits=c(-2,2),na.value = 'white',low = "#007c9b", high = "#c1f2fe",name='SBT (°C)',guide = guide_colorbar(frame.colour = 'black',ticks.colour = NA))+
    theme(aspect.ratio = 1,panel.grid.major = element_line(color = rgb(0, 0, 0,20, maxColorValue = 285), linetype = 'dashed', linewidth = 0.5),
          panel.background = element_rect(fill = NA),panel.ontop = TRUE,text = element_text(size=10),
          legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(10, 'points'),
          legend.key.width= unit(10, 'points'),axis.title = element_blank(),legend.position = c(0.15,0.33),
          panel.border = element_rect(fill = NA, colour = 'black'),legend.key = element_rect(color="black"),
          axis.text = element_text(color='black'),legend.spacing.y = unit(10, 'points'),
          axis.text.y.right = element_text(hjust= 0.1 ,margin = margin(0,0,0,-25, unit = 'points'),color='black'),
          axis.text.x = element_text(vjust = 6, margin = margin(-7,0,0,0, unit = 'points'),color='black'),
          axis.ticks.length = unit(-5,"points"))+
    scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())+
    annotate("text", x = -1526559, y = 1804909, label = y,size=8)
  
  #remove legend from plots except first (topleft) in the panel
  if (!(y %in% years[c(1,17,33)])) {
    p<-p+theme(legend.position = 'none')}
  
  #store plot on list
  plot_list[[paste0(y)]]<-p
  
  #close netcdf file
  nc_close(nc)
  
}

#save plot 16 first years (4x4)
tiff(paste0('./Figures/coldpool_ROMS_',years[1],'_',years[16],'.tiff'),height=2500,width = 2500,res = 250)
cowplot::plot_grid(plotlist = plot_list[1:16],nrow = 4,ncol = 4)
dev.off()

#save plot 16 next years (4x4)
tiff(paste0('./Figures/coldpool_ROMS_',years[17],'_',years[32],'.tiff'),height=2500,width = 2500,res = 250)
cowplot::plot_grid(plotlist = plot_list[17:32],nrow = 4,ncol = 4)
dev.off()

#save plot last years (4x4)
tiff(paste0('./Figures/coldpool_ROMS_',years[33],'_',years[length(years)],'.tiff'),height=2500,width = 2500,res = 250)
cowplot::plot_grid(plotlist = plot_list[33:length(min(years):max(years))],nrow = 4,ncol = 4)
dev.off()
