####################################################################
####################################################################
##
##    Create a Bering Sea maps
##    Sampling approach maps
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu)
##
####################################################################
####################################################################

#####################################
# Settings
#####################################

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
out_dir<-'E:/UW/Adapting Monitoring to a Changing Seascape/'
out_dir<-'C:/Users/Daniel.Vilas/Work//Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

# Define plot extent (through trial end error) units km
panel_extent <- data.frame(x = c(-1716559.21, -77636.05), #x = c(-1326559.21, -87636.05),
                           y = c(483099.5, 2194909.7)) #y = c(533099.5, 1894909.7))

#Alaska layer
ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")
ak_sppoly<-as(ebs_layers$akland, 'Spatial')

#get files from google drive and set up
files<-googledrive::drive_find()
2 #for dvilasg@uw.edu

#####################################
# Depth raster (from gebco)
#####################################

#create directory
dir.create('./bathymetry/',showWarnings = FALSE)

#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='Bathymetry'),'id']

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)
id.data<-files.1[which(files.1$name=='ak_bathy_NAD83.tiff'),]

#download file
googledrive::drive_download(file=id.data$id,
                            path = paste0('./bathymetry/',id.data$name),
                            overwrite = TRUE)

#read raster
ak_bathy_2<-raster('./bathymetry/ak_bathy_NAD83.tiff')

#####################################
# Region Shapefiles
#####################################

#create directory
dir.create('./shapefiles/',showWarnings = FALSE)

#name shapefiles 
shfiles<-c('EBSshelfThorson','NBSThorson','EBSslopeThorson')

#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='Shapefiles'),'id']

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)

#loop over shapefiles
for (i in shfiles) {
  
  #i=shfiles[1]
  
  id.data<-files.1[which(grepl(i,files.1$name)),]
  
  for (j in 1:nrow(id.data)) {
    
    #download data
    googledrive::drive_download(file=id.data$id[j],
                                path = paste0('./shapefiles/',id.data$name[j]),
                                overwrite = TRUE)
    
  }
  
  #shapefile EBS
  sh<-rgdal::readOGR(dsn='./shapefiles/',layer = i)
  
  if (i=='EBSslopeThorson') {
    
    #reproject shapefile
    proj4string(sh) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
    sh<-spTransform(sh,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))
    
  }
  
  #shapefile name
  shname<-paste0(gsub('Thorson','',i),'_sh')
  
  #assign shapefiles
  assign(shname,sh)
  
}

#merge shapefiles
bs_sh1<-union(EBSshelf_sh,NBS_sh)
bs_sh<-union(bs_sh1,EBSslope_sh)

#####################################
# Exclusive Economic Zone (EEZ)
#####################################

#create directory
dir.create('./shapefiles/',showWarnings = FALSE)

#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='EEZ'),'id']

#list of files and folder
id.data<-googledrive::drive_ls(id.bering.folder$id)

for (j in 1:nrow(id.data)) {
  
  googledrive::drive_download(file=id.data$id[j],
                              path = paste0('./shapefiles/',id.data$name[j]),
                              overwrite = TRUE)
  
}

#shapefile EEZ
eez_sh<-rgdal::readOGR(dsn='./shapefiles',layer = 'EEZ_Land_v3_202030')

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
ak_bathy_3 <- crop(ak_bathy_2, extent(bs_sh1))
ak_bathy_4 <- mask(ak_bathy_3, bs_sh1)

#positive values equal to zero and convert negative to positive values
ak_bathy_4[ak_bathy_4>0]<-0 
ak_bathy_4<--ak_bathy_4 

#####################################
# Current sampling stations 
#####################################

#create directory
dir.create('./additional/',showWarnings = FALSE)

#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='Additional'),'id']

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)
id.data<-files.1[which(files.1$name=='ebs_nbs_temperature_full_area.csv'),]

#download file
googledrive::drive_download(file=id.data$id,
                            path = paste0('./additional/',id.data$name),
                            overwrite = TRUE)

#extract station EBS bottom trawl
st_EBS<-read.csv('./additional//ebs_nbs_temperature_full_area.csv')

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

#####################################
# Baseline strata
#####################################

#load baseline strata
load('./output/baseline_strata.RData')

baseline_strata$locations[grep('GF|HG|JI|IH|ON|QP|PO',baseline_strata$locations$stationid),]
baseline_strata$locations$corner<-ifelse(grepl('GF|HG|JI|IH|ON|QP|PO',baseline_strata$locations$stationid),'TRUE','FALSE')

baseline_strata$locations<-baseline_strata$locations[order(baseline_strata$locations$cell),]
baseline_strata$locations$difference <- c( NA, diff( baseline_strata$locations$cell ) )
mean(baseline_strata$locations$difference,na.rm=TRUE)/2 #so 50

#baseline_strata
pts<-baseline_strata$locations
#coordinates(pts)<-~longitude + latitude
#proj4string(pts) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
#pts<-spTransform(pts,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))
#pts<-as.data.frame(pts)

#####################################
# EBS+NBS grid 
#####################################

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)
#add col and row number
x1<-grid[,c('Lon','Lat','cell',"Stratum")]
names(x1)<-c('x','y','z','strata')
coordinates(x1)=~x + y
crs(x1)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
x2<-spTransform(x1,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
x3<-data.frame(x2)
x3$x<-as.integer(x3$x)
x3$y<-as.integer(x3$y)
lon<-sort(unique(x3$x),decreasing = FALSE) #1556
lat<-sort(unique(x3$y),decreasing = TRUE) #1507
lons<-data.frame(x=lon,col=1:length(lon))
lats<-data.frame(y=lat,row=1:length(lat))
x4<-merge(x3,lons,by='x',all.x=TRUE)
x5<-merge(x4,lats,by='y',all.x=TRUE)
colnames(x5)<-c('Lat','Lon','cell','strata','optional','col','row')
grid<-x5[,c('Lat','Lon','cell','strata','col','row')]
grid1<-as.data.frame(grid)

#####################################
# Stations example for buffer approach
#####################################

#cell1 to represent
i_cell2<-10000
#row
row<-grid1[which(grid1$cell %in% i_cell2),c('row','col')][1,1]
#col
col<-grid1[which(grid1$cell %in% i_cell2),c('row','col')][1,2]
#get 100 adjacent cells
adj_cells<-expand.grid(row=c(row,row+(1:100),row-(1:100)),
                       col=c(col,col+(1:100),col-(1:100)))
adj_cells1<-merge(adj_cells,grid1,by=c('row', 'col'))#[,'cell']
adj_cells2<-adj_cells1

#select buffer and select 10 cells
dff<-subset(grid,strata==70)

#create vectors to store cells
dropcell<-c()
selcell<-c()

#loop over required samples
for (iii in rep(1,times=10)) {
  
  #iii<-1
  
  #random sample
  cell_i<-sample(dff$cell,iii)
  #get row of selected sample
  row<-dff[which(dff$cell==cell_i),c('row','col')][1,1]
  #get col of selected sample
  col<-dff[which(dff$cell==cell_i),c('row','col')][1,2]
  #get adjacent cells of selected sample
  adj_cells<-expand.grid(row=c(row,row+(1:100),row-(1:100)),
                         col=c(col,col+(1:100),col-(1:100)))
  adj_cells1<-merge(adj_cells,grid,by=c('row', 'col'))[,'cell']
  #drop cells is equal to selected cell and adjacent cells
  dropcell_i<-c(cell_i,adj_cells1)
  
    dff<-subset(dff, !(cell %in% dropcell_i))
    selcell<-c(selcell,cell_i)
    dropcell<-c(dropcell,dropcell_i)
  
}
  pointsb<-subset(grid,cell %in% selcell)#[,c('Lon','Lat','cell','strata')]
  buffer<-subset(grid,cell %in% dropcell)
  #names(pointsb)<-c('Lon','Lat','cell','strata')


#####################################
# Plot approaches
#####################################

#CURRENT
plot_list<-list()

for (y in c(1,2,3,4)) {
  
  x<-
    ggplot()+
    geom_sf(data=ebs_layers$survey.strata,fill = NA)+
    geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
    scale_x_continuous(expand = c(0,0),position = 'bottom',
                       breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
    geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
    geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
    coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
             #xlim = c(panel_extent$x[1]+200000,panel_extent$x[2]),
             xlim = c(panel_extent$x[1]+200000,panel_extent$x[2]+100000),
             ylim = c(panel_extent$y[1]-100000,panel_extent$y[2]-200000),
             label_axes = "-NE-")+
    theme(aspect.ratio = 1,panel.grid.major = element_line(color = rgb(0, 0, 0,20, maxColorValue = 285), linetype = 'dashed', linewidth =  0.5),
          panel.background = element_rect(fill = NA),panel.ontop = TRUE,text = element_text(size=10),
          legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(25, 'points'),
          legend.key.width= unit(25, 'points'),axis.title = element_blank(),legend.position = c(0.12,0.47), #c(0.12,0.47)
          panel.border = element_rect(fill = NA, colour = 'black'),legend.key = element_rect(color="black"),
          axis.text = element_text(color='black'),legend.spacing.y = unit(10, 'points'),
          axis.text.y.right = element_text(hjust= 0.1 ,margin = margin(0,7,0,-25, unit = 'points'),color='black'),
          axis.text.x = element_text(vjust = 6, margin = margin(-7,0,7,0, unit = 'points'),color='black'),
          axis.ticks.length = unit(-5,"points"))+
    scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())#+
  
  x<-
    if (y==1) {
      x + geom_point(data=baseline_strata$locations,aes(x=longitude,y=latitude),size=0.1,alpha=0.1) #53464
    } else if (y==2) {
      x + geom_point(data=subset(baseline_strata$locations,stratum==70),aes(x=longitude,y=latitude),size=0.1,alpha=0.1) #6002 available points
    } else if (y==3) {
      x +
        geom_point(data=subset(baseline_strata$locations,stratum==70),aes(x=longitude,y=latitude),size=0.1,alpha=0.1)+
        geom_point(data=subset(baseline_strata$locations,stratum==70 & cell==14595),aes(x=longitude,y=latitude),size=1,color='white',fill='black',shape=21)
    } else if (y==4){
      x + 
        geom_point(data=subset(baseline_strata$locations,stratum==70),aes(x=longitude,y=latitude),size=0.1,alpha=0.1)+
        geom_point(data=subset(baseline_strata$locations,stratum==70 & cell==14595),aes(x=longitude,y=latitude),size=1,color='white',fill='black',shape=21)+
        geom_point(data=subset(baseline_strata$locations,stratum==70)[sample(nrow(subset(baseline_strata$locations,stratum==70 | cell==14595)),size=19),],aes(x=longitude,y=latitude),size=1,color='white',fill='black',shape=21)
      
    }
  
  plot_list[[y]]<-x
  
}

#save plot
ragg::agg_png(paste0('./figures/current_sampling.png'), width = 12, height = 4, units = "in", res = 300)
gridExtra::grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],nrow=1)
dev.off()

#RANDOM
plot_list<-list()

  for (y in c(1,2,3,4)) {
    
  x<-
    ggplot()+
    geom_sf(data=ebs_layers$survey.strata,fill = NA)+
    geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
    scale_x_continuous(expand = c(0,0),position = 'bottom',
                       breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
    geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
    geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
    coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
             #xlim = c(panel_extent$x[1]+200000,panel_extent$x[2]),
             xlim = c(panel_extent$x[1]+200000,panel_extent$x[2]+100000),
             ylim = c(panel_extent$y[1]-100000,panel_extent$y[2]-200000),
             label_axes = "-NE-")+
    theme(aspect.ratio = 1,panel.grid.major = element_line(color = rgb(0, 0, 0,20, maxColorValue = 285), linetype = 'dashed', linewidth =  0.5),
          panel.background = element_rect(fill = NA),panel.ontop = TRUE,text = element_text(size=10),
          legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(25, 'points'),
          legend.key.width= unit(25, 'points'),axis.title = element_blank(),legend.position = c(0.12,0.47), #c(0.12,0.47)
          panel.border = element_rect(fill = NA, colour = 'black'),legend.key = element_rect(color="black"),
          axis.text = element_text(color='black'),legend.spacing.y = unit(10, 'points'),
          axis.text.y.right = element_text(hjust= 0.1 ,margin = margin(0,7,0,-25, unit = 'points'),color='black'),
          axis.text.x = element_text(vjust = 6, margin = margin(-7,0,7,0, unit = 'points'),color='black'),
          axis.ticks.length = unit(-5,"points"))+
    scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())#+
  
  x<-
    if (y==1) {
      x + geom_point(data=grid,aes(x=Lon,y=Lat),size=0.1,alpha=0.1) #53464
    } else if (y==2) {
      x + geom_point(data=subset(grid,strata==70),aes(x=Lon,y=Lat),size=0.1,alpha=0.1) #6002 available points
    } else if (y==3) {
      x +
        geom_point(data=subset(grid,strata==70),aes(x=Lon,y=Lat),size=0.1,alpha=0.1)+
        geom_point(data=subset(grid,strata==70 & cell==12380),aes(x=Lon,y=Lat),size=1,color='white',fill='black',shape=21)
    } else if (y==4){
      x + 
        geom_point(data=subset(grid,strata==70),aes(x=Lon,y=Lat),size=0.1,alpha=0.1)+
        geom_point(data=subset(grid,strata==70 & cell==12380),aes(x=Lon,y=Lat),size=1,color='white',fill='black',shape=21)+
        geom_point(data=subset(grid,strata==70)[sample(nrow(subset(grid,strata==70 | cell==12380)),size=19),],aes(x=Lon,y=Lat),size=1,color='white',fill='black',shape=21)
      
    }
  
  plot_list[[y]]<-x
  
  }

#save plot
ragg::agg_png(paste0('./figures/random_sampling.png'), width = 12, height = 4, units = "in", res = 300)
gridExtra::grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],nrow=1)
dev.off()

#BUFFER
plot_list<-list()

for (y in c(1,2,3,4)) {
  
  x<-
    ggplot()+
    geom_sf(data=ebs_layers$survey.strata,fill = NA)+
    geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
    scale_x_continuous(expand = c(0,0),position = 'bottom',
                       breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
    geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
    geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
    coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
             #xlim = c(panel_extent$x[1]+200000,panel_extent$x[2]),
             xlim = c(panel_extent$x[1]+200000,panel_extent$x[2]+100000),
             ylim = c(panel_extent$y[1]-100000,panel_extent$y[2]-200000),
             label_axes = "-NE-")+
    theme(aspect.ratio = 1,panel.grid.major = element_line(color = rgb(0, 0, 0,20, maxColorValue = 285), linetype = 'dashed', linewidth =  0.5),
          panel.background = element_rect(fill = NA),panel.ontop = TRUE,text = element_text(size=10),
          legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(25, 'points'),
          legend.key.width= unit(25, 'points'),axis.title = element_blank(),legend.position = c(0.12,0.47), #c(0.12,0.47)
          panel.border = element_rect(fill = NA, colour = 'black'),legend.key = element_rect(color="black"),
          axis.text = element_text(color='black'),legend.spacing.y = unit(10, 'points'),
          axis.text.y.right = element_text(hjust= 0.1 ,margin = margin(0,7,0,-25, unit = 'points'),color='black'),
          axis.text.x = element_text(vjust = 6, margin = margin(-7,0,7,0, unit = 'points'),color='black'),
          axis.ticks.length = unit(-5,"points"))+
    scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())#+
  
  x<-
    if (y==1) {
      x + geom_point(data=grid,aes(x=Lon,y=Lat),size=0.1,alpha=0.1) #53464
    } else if (y==2) {
      x + geom_point(data=subset(grid,strata==70),aes(x=Lon,y=Lat),size=0.1,alpha=0.1) #6002 available points
    } else if (y==3) {
      x +
        geom_point(data=subset(grid,strata==70),aes(x=Lon,y=Lat),size=0.1,alpha=0.1)+
        geom_point(data=adj_cells2,aes(x=Lon,y=Lat),size=0.1,color='white')+
        geom_point(data=subset(grid,strata==70 & cell==i_cell2),aes(x=Lon,y=Lat),size=1,color='white',fill='black',shape=21)
    } else if (y==4){
      x + 
        geom_point(data=subset(grid,strata==70),aes(x=Lon,y=Lat),size=0.1,alpha=0.1)+
        geom_point(data=buffer,aes(x=Lon,y=Lat),size=0.1,color='white')+
        geom_point(data=pointsb,aes(x=Lon,y=Lat),size=1,color='white',fill='black',shape=21)
      
    }
  
  plot_list[[y]]<-x
  
}

#save plot
ragg::agg_png(paste0('./figures/buffer_sampling.png'), width = 12, height = 4, units = "in", res = 300)
gridExtra::grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],nrow=1)
dev.off()

###########################
# Plot study location map
###########################

#zoomin plot
zoomin<-
  gplot(ak_bathy_4) +
  #x<-ggplot()+
    geom_tile(aes(fill=value))+
    #geom_sf(data=ebs_layers$survey.strata,fill = NA)+
    #geom_point(data=subset(pts,corner=='FALSE'),aes(x=longitude,y=latitude),size=1)+
    #geom_point(data=subset(pts,rm25==0),aes(x=longitude,y=latitude),shape=4,size=1,col='red')+
    geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
    geom_polygon(data=eez_sh33,aes(x=long,y=lat,group=group),fill=NA,color='grey40')+
    scale_x_continuous(expand = c(0,0),position = 'bottom',
                       breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
    geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
    geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
    #geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
    coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
             xlim = panel_extent$x,
             ylim = panel_extent$y,
             #lim = c(panel_extent$x[1]+200000,panel_extent$x[2]),
             #ylim = c(panel_extent$y[1],panel_extent$y[2]-200000),
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
          legend.key.width= unit(25, 'points'),axis.title = element_blank(),legend.position = c(0.12,0.47), #c(0.12,0.47)
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
          colour = '#800000', linetype='solid', fill = NA,linewidth=0.7) +
    annotate("text", x = -2300000, y = 2600000, label = "Russia",parse=TRUE,size=4)+
    annotate("text", x = 0, y = -1000000, label = "italic('Pacific Ocean')",parse=TRUE,size=4)+
    annotate("text", x = +2000000, y = 1400000, label = "Canada",parse=TRUE,size=4)+
    annotate("text", x = +3000000, y = 0, label = "USA",parse=TRUE,size=4)+
    annotation_scale()

#save plot
ragg::agg_png(paste0('./figures/map_bering.png'), width = 7, height = 7, units = "in", res = 300)
grid.newpage()
vp_b <- viewport(width = 1, height = 1, x = 0.5, y = 0.5)  # the larger map
vp_a <- viewport(width = 0.4, height = 0.3, x = 0.219, y = 0.846)  # the inset in upper left
print(zoomin , vp = vp_b)
print(zoomout , vp = vp_a)
dev.off()


