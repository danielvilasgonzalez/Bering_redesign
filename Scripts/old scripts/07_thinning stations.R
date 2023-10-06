####################################################################
####################################################################
##
##    Sampling scenarios maps
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu)
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c('rasterVis','rgeos','scales','rnaturalearth','grid','ggplot2',"SamplingStrata")

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
1 #for dvilasg@uw.edu

#####################################
# GET DEPTH RASTER
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
# BERING SHAPEFILES
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
bs_sh1<-terra::union(EBSshelf_sh,NBS_sh)
bs_sh<-terra::union(bs_sh1,EBSslope_sh)

#####################################
# GET EEZ
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
ak_bathy_3 <- crop(ak_bathy_2, extent(bs_sh))
ak_bathy_4 <- mask(ak_bathy_3, bs_sh)

#positive values equal to zero and convert negative to positive values
ak_bathy_4[ak_bathy_4>0]<-0 
ak_bathy_4<--ak_bathy_4 

#####################################
# GET Stations 
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
st_EBS1<-spTransform(st_EBS,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))

#convert to dataframe
st_EBS<-as.data.frame(st_EBS)
st_EBS$survey_definition_id<-as.factor(st_EBS$survey_definition_id)

#convert to dataframe
st_EBS1<-as.data.frame(st_EBS1)
st_EBS1$survey_definition_id<-as.factor(st_EBS1$survey_definition_id)
st_EBS1<-st_EBS1[which(st_EBS1$survey_definition_id=='98'),]

#remove some polygons of EEZ object
eez_sh11 <- eez_sh1[eez_sh1$AREA_KM2 == 5193061,] 
eez_sh22 <- eez_sh2[eez_sh2$AREA_KM2 == 5193061,]  #"5193061"  "24614858" "8521"    

#join both polygons without inner line
eez_sh3<-aggregate(rbind(eez_sh11,eez_sh22),dissolve=T)
eez_sh33<-rgeos::gUnaryUnion(eez_sh3)

#corner stations
st_corners1<-st_EBS1[which(nchar(st_EBS1$stationid)>=6 & st_EBS1$stationid!='AZ0504'),]
st_EBS2<-st_EBS1[which(nchar(st_EBS1$stationid)<=5 | st_EBS1$stationid=='AZ0504'),]

#########################
#current design
#########################

ggplot()+
  geom_point(data=st_EBS2,aes(x=longitude,y=latitude),shape=4,size=1)+
  geom_point(data=st_corners1,aes(x=longitude,y=latitude),color='red',shape=20,size=1)+
  geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
  geom_polygon(data=eez_sh33,aes(x=long,y=lat,group=group),fill=NA,color='grey40')+
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
  geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           xlim = panel_extent$x,
           ylim = panel_extent$y,
           label_axes = "-NE-")+
  
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

#############################################
# VARIANCE SBT STATIONS
############################################

#extract station EBS bottom trawl
st_EBS<-read.csv('./additional//ebs_nbs_temperature_full_area.csv')
st_EBS2019<-subset(st_EBS,year==2019 & survey_definition_id ==98 )

#unique station ID
st_all<-unique(st_EBS$stationid)

#st from NBS and EBS
st_EBS1<-st_EBS[which(st_EBS$survey_definition_id=='98'),]
st_NBS1<-st_EBS[which(st_EBS$survey_definition_id=='143'),]
st_NBS2<-st_EBS[which(st_EBS$survey_definition_id=='143'),'stationid']

#st crab corner
st_crab<-st_EBS1[which(nchar(st_EBS1$stationid)>=6 & st_EBS1$stationid!='AZ0504'),'stationid']

#EBS stations without crab corner
st_EBS2<-st_EBS1[which(st_EBS1$survey_definition_id=='98'),'stationid']
st_EBS2<-setdiff(st_EBS2,st_crab)
st_EBS11<-st_EBS1[which(st_EBS1$stationid %in% st_EBS2),]

#get mean lon and lat
Lon<-aggregate(st_EBS$longitude, by=list(st_EBS$stationid), FUN=mean)
#Lon<-Lon[which(Lon$Group.1 %in% st_EBS2),]
Lat<-aggregate(st_EBS$latitude, by=list(st_EBS$stationid), FUN=mean)
#Lat<-Lat[which(Lat$Group.1 %in% st_EBS2),]

#get variance
var<-aggregate(st_EBS$gear_temperature, by=list('stationid'=st_EBS$stationid), FUN=var,na.rm=TRUE)  
#var1<-var[which(var$stationid %in% st_EBS2),]
df_var<-data.frame('stationid'=Lon$Group.1,'Lon'=Lon$x,'Lat'=Lat$x)
df_var<-merge(df_var,var,by='stationid')
df_var<-df_var[complete.cases(df_var),]

#convert to spatial object to change lat lon
coordinates(df_var)<- ~Lon + Lat
proj4string(df_var) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
#change CRS
df_var<-spTransform(df_var,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))

#convert to df
df_var<-as.data.frame(df_var)

#plot
ggplot()+
  geom_point(data=df_var,aes(x=Lon,y=Lat,color=x),size=3)+
  #geom_point(data=st_corners1,aes(x=longitude,y=latitude),color='red',shape=20,size=1)+
  geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
  geom_polygon(data=eez_sh33,aes(x=long,y=lat,group=group),fill=NA,color='grey40')+
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
  geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           xlim = panel_extent$x,
           ylim = panel_extent$y,
           label_axes = "-NE-")+
  scale_color_viridis_c('SBT variance')+
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

#############################################
# SBT cluster STRATIFICATION
############################################

#get mean sbt
st_EBScurrent<-subset(st_EBS,stationid %in% st_EBS2019$stationid)

mean<-aggregate(st_EBS$gear_temperature, by=list('stationid'=st_EBS$stationid), FUN=mean,na.rm=TRUE)  
#var1<-var[which(var$stationid %in% st_EBS2),]
df_var<-merge(df_var,mean,by='stationid')
df_var<-df_var[complete.cases(df_var),]

#get mean depth
mean<-aggregate(st_EBS$bottom_depth, by=list('stationid'=st_EBS$stationid), FUN=mean,na.rm=TRUE)  
#var1<-var[which(var$stationid %in% st_EBS2),]
df_var<-merge(df_var,mean,by='stationid')
df_var<-df_var[complete.cases(df_var),]

#rename cols
colnames(df_var)<-c('stationid','sbt_var','Lon','Lat','sbt_mean','depth')

x<-st_EBScurrent[,c("year","stationid","gear_temperature")]
x<-x[!duplicated(x[,c('year','stationid')]),]

library(reshape2)
xx<-dcast(x, formula = stationid ~ year, value="gear_temperature", fill=0)
row.names(xx)<-xx[,1]
xx<-xx[,-1]
#colnames(xx)<-NULL

#x<-as.data.frame((df_var[,c("depth")]))
xxx<-as.matrix(dist(xx))
library(rioja)
clust <- rioja::chclust(dist(xx))
bstick(clust)
cut_avg <- as.data.frame(cutree(clust, k = 2))
y<-kmeans(x = as.matrix(dist(xx)),centers = 4) 
cut_avg1<-data.frame('stationid'=rownames(cut_avg),'cut1'=cut_avg$`cutree(clust, k = 2)`,'cut2'=y$cluster)
plot(clust)
rect.hclust(clust , k = 2, border = 2:6)

df_var1<-merge(df_var,cut_avg1,by='stationid')
df_var1$cut2<-as.factor(df_var1$cut2)

#plot
ggplot()+
  geom_point(data=df_var1,aes(x=Lon,y=Lat,color=cut2),size=3)+
  #geom_point(data=st_corners1,aes(x=longitude,y=latitude),color='red',shape=20,size=1)+
  geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
  geom_polygon(data=eez_sh33,aes(x=long,y=lat,group=group),fill=NA,color='grey40')+
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
  geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           xlim = panel_extent$x,
           ylim = panel_extent$y,
           label_axes = "-NE-")+
  scale_color_viridis_d('clusters')+
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



#############################################
# SBT cluster STRATIFICATION
############################################

df<-st_EBS[which(st_EBS$stationid %in% st_EBS2),c('stationid','gear_temperature','year','longitude','latitude')]
xx<-dcast(x, formula = stationid ~ year, value="gear_temperature", fill=0)
df_anom<-data.frame('stationid'=xx$stationid)

for (y in 1983:2019) {
  #y<-1983
  
  anom<-xx[,paste0(y-1)]-xx[,paste0(y)]
  anom<-data.frame('stationid'=xx$stationid,anom)
  df_anom1<-merge(df_anom,anom,by='stationid',all.x=TRUE)
  df_anom<-cbind(df_anom,df_anom1$anom)
  
  
  
  # #df_anom<-
  # df1<-subset(df,year==y)
  # df1<-df1[order(df1$stationid),]
  # df1_st<-unique(df1$stationid)
  # df2<-subset(df,year==y-1)
  # df2<-df2[order(df2$stationid),]
  # df2_st<-unique(df2$stationid)
  # 
  # st<-intersect(df1_st,df2_st)
  # 
  # df11<-df1[df1$stationid %in% st,]
  # df11 <- df11[!duplicated(df11$stationid), ]
  # df22<-df2[df2$stationid %in% st,]
  # df22 <- df22[!duplicated(df22$stationid), ]
  # #setdiff(df11$stationid,df22$stationid)
  # 
  # anom<-df11$gear_temperature-df22$gear_temperature
  # #length(anom)
  # #length(st)
  # anom<-data.frame('stationid'=st,anom)
  # df_anom1<-merge(df_anom,anom,by='stationid',all.x=TRUE)
  # 
  # df_anom<-cbind(df_anom,df_anom1$anom)
  # 
  
}


colnames(df_anom)<-c('stationid',1983:2019)
df_anom[,4:ncol(df_anom)]<-abs(df_anom[,4:ncol(df_anom)])
#View(df_anom)
#df_anom1<-df_anom[which(df_anom$stationid %in% st_EBScurrent$stationid),]



xx<-(df_anom1[,2:ncol(df_anom1)])
#xx<-xx[,-c(42,41)]

y<-kmeans(x = as.matrix(dist(xx)),centers = 5) 
cut_avg1<-data.frame('stationid'=rownames(cut_avg),'cut3'=cut_avg$`cutree(clust, k = 2)`,'cut4'=y$cluster)
plot(clust)
rect.hclust(clust , k = 2, border = 2:6)

df_var1<-merge(df_var,cut_avg1,by='stationid')
df_var1$cut4<-as.factor(df_var1$cut4)
#df_anom$total<-rowMeans(df_anom[,2:ncol(df_anom)],na.rm=TRUE)
#df_anom1<-df_anom[order(df_anom$total,decreasing = TRUE),]


#plot
ggplot()+
  geom_point(data=df_var1,aes(x=Lon,y=Lat,color=cut4),size=3)+
  #geom_point(data=st_corners1,aes(x=longitude,y=latitude),color='red',shape=20,size=1)+
  geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
  geom_polygon(data=eez_sh33,aes(x=long,y=lat,group=group),fill=NA,color='grey40')+
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
  geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           xlim = panel_extent$x,
           ylim = panel_extent$y,
           label_axes = "-NE-")+
  scale_color_viridis_d('clusters')+
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









#st_anom<-data.frame('stationid'=Lon$Group.1,'Lon'=Lon$x,'Lat'=Lat$x)
#st_anom$total<-df_anom$total

coordinates(st_anom)<- ~Lon + Lat
proj4string(st_anom) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 

#change CRS
st_anom<-spTransform(st_anom,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))

st_anom<-as.data.frame(st_anom)

#get mean sbt
st_EBScurrent<-subset(st_EBS,stationid %in% st_EBS2019$stationid)

mean<-aggregate(st_EBS$gear_temperature, by=list('stationid'=st_EBS$stationid), FUN=mean,na.rm=TRUE)  
#var1<-var[which(var$stationid %in% st_EBS2),]
df_var<-merge(df_var,mean,by='stationid')
df_var<-df_var[complete.cases(df_var),]

#get mean depth
mean<-aggregate(st_EBS$bottom_depth, by=list('stationid'=st_EBS$stationid), FUN=mean,na.rm=TRUE)  
#var1<-var[which(var$stationid %in% st_EBS2),]
df_var<-merge(df_var,mean,by='stationid')
df_var<-df_var[complete.cases(df_var),]

#rename cols
colnames(df_var)<-c('stationid','sbt_var','Lon','Lat','sbt_mean','depth')

x<-st_EBScurrent[,c("year","stationid","gear_temperature")]
x<-x[!duplicated(x[,c('year','stationid')]),]

library(reshape2)
xx<-dcast(x, formula = stationid ~ year, value="gear_temperature", fill=0)
row.names(xx)<-xx[,1]
xx<-xx[,-1]
#colnames(xx)<-NULL

#x<-as.data.frame((df_var[,c("depth")]))
xxx<-as.matrix(dist(xx))
library(rioja)
clust <- rioja::chclust(dist(xx))
bstick(clust)
cut_avg <- as.data.frame(cutree(clust, k = 2))
y<-kmeans(x = as.matrix(dist(xx)),centers = 4) 
cut_avg1<-data.frame('stationid'=rownames(cut_avg),'cut1'=cut_avg$`cutree(clust, k = 2)`,'cut2'=y$cluster)
plot(clust)
rect.hclust(clust , k = 2, border = 2:6)

df_var1<-merge(df_var,cut_avg1,by='stationid')
df_var1$cut2<-as.factor(df_var1$cut2)

#plot
ggplot()+
  geom_point(data=df_var1,aes(x=Lon,y=Lat,color=cut2),size=3)+
  #geom_point(data=st_corners1,aes(x=longitude,y=latitude),color='red',shape=20,size=1)+
  geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
  geom_polygon(data=eez_sh33,aes(x=long,y=lat,group=group),fill=NA,color='grey40')+
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
  geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           xlim = panel_extent$x,
           ylim = panel_extent$y,
           label_axes = "-NE-")+
  scale_color_viridis_d('clusters')+
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





#abline(h = 3, col = 'red')

clust$call
bstick.chclust(clust)

# Basic diagram
plot(clust, hang=-1)
# Rotated through 90 degrees
plot(clust, hang=-1, horiz=TRUE)
# Rotated and observations plotted according to sample depth.
plot(clust, xvar=x$gear_temperature, hang=-1, horiz=TRUE, x.rev=TRUE)

y<-kmeans(x = as.matrix(dist(xx)),centers = 5) 

plot(xx, col = y$cluster)




#build df for strata
frame1 <- buildFrameDF(df = df_var,
                       id = "stationid",
                       X = c("Lat","Lon"),
                       Y = c("depth"),
                       domainvalue = "REG")




















######################
#thinning 50km
library(BiodiversityR)
#######################

#using CRS for thinning
st_EBS12<-st_EBS[which(nchar(st_EBS$stationid)<=5 | st_EBS$stationid=='AZ0504'),]
st_EBS12<-st_EBS12[which(st_EBS12$survey_definition_id=='98'),]
lonlat.df<-data.frame('Lon'=st_EBS12$longitude,'Lat'=st_EBS12$latitude)
dim(lonlat.df)
#spatial thinning 50km
xy50<-ensemble.spatialThin(lonlat.df, thin.km = 50, 
                         runs = 100, silent = FALSE, verbose = FALSE, 
                         return.notRetained = FALSE)


dim(xy50)
xy50<-data.frame('longitude'=xy50[,1],'latitude'=xy50[,2])
#convert to spatial object
coordinates(xy50)<-c('longitude','latitude')
proj4string(xy50) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 

#change CRS
xy501<-spTransform(xy50,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))

#convert to dataframe
xy501<-as.data.frame(xy501)

#plot
ggplot()+
  geom_point(data=xy501,aes(x=longitude,y=latitude),shape=4,size=1)+
  geom_point(data=st_corners1,aes(x=longitude,y=latitude),color='red',shape=20,size=1)+
  geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
  geom_polygon(data=eez_sh33,aes(x=long,y=lat,group=group),fill=NA,color='grey40')+
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
  geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           xlim = panel_extent$x,
           ylim = panel_extent$y,
           label_axes = "-NE-")+
  
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

######################
#thinning 40km
library(BiodiversityR)
#######################

#spatial thinning 40km
xy40<-ensemble.spatialThin(lonlat.df, thin.km = 40, 
                           runs = 100, silent = FALSE, verbose = FALSE, 
                           return.notRetained = FALSE)

dim(xy40)
xy40<-data.frame('longitude'=xy40[,1],'latitude'=xy40[,2])
#convert to spatial object
coordinates(xy40)<-c('longitude','latitude')
proj4string(xy40) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 

#change CRS
xy401<-spTransform(xy40,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))

#convert to dataframe
xy401<-as.data.frame(xy401)

#plot
ggplot()+
  geom_point(data=xy401,aes(x=longitude,y=latitude),shape=4,size=1)+
  geom_point(data=st_corners1,aes(x=longitude,y=latitude),color='red',shape=20,size=1)+
  geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
  geom_polygon(data=eez_sh33,aes(x=long,y=lat,group=group),fill=NA,color='grey40')+
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
  geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           xlim = panel_extent$x,
           ylim = panel_extent$y,
           label_axes = "-NE-")+
  
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

######################
#thinning 35km
library(BiodiversityR)
#######################

#spatial thinning 35km
xy35<-ensemble.spatialThin(lonlat.df, thin.km = 35, 
                           runs = 100, silent = FALSE, verbose = FALSE, 
                           return.notRetained = FALSE)

dim(xy35)
xy35<-data.frame('longitude'=xy35[,1],'latitude'=xy35[,2])
#convert to spatial object
coordinates(xy35)<-c('longitude','latitude')
proj4string(xy35) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 

#change CRS
xy351<-spTransform(xy35,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))

#convert to dataframe
xy351<-as.data.frame(xy351)

#plot
ggplot()+
  geom_point(data=xy351,aes(x=longitude,y=latitude),shape=4,size=1)+
  geom_point(data=st_corners1,aes(x=longitude,y=latitude),color='red',shape=20,size=1)+
  geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
  geom_polygon(data=eez_sh33,aes(x=long,y=lat,group=group),fill=NA,color='grey40')+
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
  geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           xlim = panel_extent$x,
           ylim = panel_extent$y,
           label_axes = "-NE-")+
  
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


######################
#thinning random
library(BiodiversityR)
#######################

points<-sample(1:350,35,replace=FALSE)
points<-sample(1:350,18,replace=FALSE)

st_EBS21<-st_EBS2[-points,]

#plot
ggplot()+
  geom_point(data=st_EBS21,aes(x=longitude,y=latitude),shape=4,size=1)+
  geom_point(data=st_corners1,aes(x=longitude,y=latitude),color='red',shape=20,size=1)+
  geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
  geom_polygon(data=eez_sh33,aes(x=long,y=lat,group=group),fill=NA,color='grey40')+
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
  geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           xlim = panel_extent$x,
           ylim = panel_extent$y,
           label_axes = "-NE-")+
  
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


######################
#thinning by depth
#######################

st_EBS2nbs<-st_EBS2[which(st_EBS2$survey_definition_id==143),]
st_EBS2ebs<-st_EBS2[which(st_EBS2$survey_definition_id==98),]
st_EBS2ebs<-st_EBS2ebs[order(st_EBS2ebs$bottom_depth,decreasing = TRUE),]



st_EBS21ebs<-st_EBS2ebs[-c(1:50),]

st_EBS21<-rbind(st_EBS21ebs,st_EBS2nbs)
#plot
ggplot()+
  geom_point(data=st_EBS21,aes(x=longitude,y=latitude),shape=4,size=1)+
  geom_point(data=st_corners1,aes(x=longitude,y=latitude),color='red',shape=20,size=1)+
  geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
  geom_polygon(data=eez_sh33,aes(x=long,y=lat,group=group),fill=NA,color='grey40')+
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
  geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           xlim = panel_extent$x,
           ylim = panel_extent$y,
           label_axes = "-NE-")+
  
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


#read raster
r<-raster('./bathymetry/gebco_2022_n70.0_s50.0_w-180.0_e-155.0.asc')
bs_sh12<-spTransform(xy35,CRSobj = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") )
r1<-crop(r, extent(bs_sh12))
r2 <- mask(r1, bs_sh12)
r2[r2>0]<-NA

dim(lonlat)

xy50<-ensemble.environmentalThin(lonlat.df, predictors.stack = r2, 
                               extracted.data=NULL, thin.n = 400,
                               runs = 100, pca.var = 0.95, silent = FALSE, verbose = FALSE,
                               return.notRetained = FALSE)


xy50<-data.frame('longitude'=xy50[,1],'latitude'=xy50[,2])
#convert to spatial object
coordinates(xy50)<-c('longitude','latitude')
proj4string(xy50) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 

#change CRS
xy501<-spTransform(xy50,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))

#convert to dataframe
xy501<-as.data.frame(xy501)
#plot
ggplot()+
  geom_point(data=xy501,aes(x=longitude,y=latitude),shape=4,size=1)+
  geom_point(data=st_corners1,aes(x=longitude,y=latitude),color='red',shape=20,size=1)+
  geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
  geom_polygon(data=eez_sh33,aes(x=long,y=lat,group=group),fill=NA,color='grey40')+
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
  geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           xlim = panel_extent$x,
           ylim = panel_extent$y,
           label_axes = "-NE-")+
  
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

##############################
# environmental spatial 
###############################

#extract station EBS bottom trawl
st_EBS<-read.csv('./additional//ebs_nbs_temperature_full_area.csv')

#unique station ID
st_all<-unique(st_EBS$stationid)

#convert to dataframe
st_EBS1<-st_EBS1[which(st_EBS1$survey_definition_id=='98'),]
st_NBS1<-st_EBS[which(st_EBS$survey_definition_id=='143'),]

#corner stations
st_crab<-st_EBS1[which(nchar(st_EBS1$stationid)>=6 & st_EBS1$stationid!='AZ0504'),'stationid']

#EBS stations without crab
st_EBS2<-st_EBS1[which(st_EBS1$survey_definition_id=='98'),'stationid']
st_EBS2<-setdiff(st_EBS2,st_crab)
st_EBS11<-st_EBS1[which(st_EBS1$stationid %in% st_EBS2),]

#get mean lon and lat
Lon<-aggregate(st_EBS$longitude, by=list(st_EBS$stationid), FUN=mean)
Lon<-Lon[which(Lon$Group.1 %in% st_EBS2),]
Lat<-aggregate(st_EBS$latitude, by=list(st_EBS$stationid), FUN=mean)
Lat<-Lat[which(Lat$Group.1 %in% st_EBS2),]

#NBS stations
st_NBS2<-st_EBS[which(st_EBS$survey_definition_id=='143'),'stationid']

df<-st_EBS[which(st_EBS$stationid %in% st_EBS2),c('stationid','gear_temperature','year','longitude','latitude')]

var<-aggregate(st_EBS$gear_temperature, by=list('stationid'=st_EBS$stationid), FUN=var,na.rm=TRUE)  
var1<-var[which(var$stationid %in% st_EBS2),]
df_var<-data.frame('stationid'=Lon$Group.1,'Lon'=Lon$x,'Lat'=Lat$x)
df_var<-merge(df_var,var1,by='stationid')

df_anom<-data.frame('stationid'=Lon$Group.1,'Lon'=Lon$x,'Lat'=Lat$x)

for (y in 1983:2022) {
  #y<-1983
  
  #df_anom<-
  
    
    
    
  df1<-subset(df,year==y)
  df1<-df1[order(df1$stationid),]
  df1_st<-unique(df1$stationid)
  df2<-subset(df,year==y-1)
  df2<-df2[order(df2$stationid),]
  df2_st<-unique(df2$stationid)
  
  st<-intersect(df1_st,df2_st)
 
  df11<-df1[df1$stationid %in% st,]
  df11 <- df11[!duplicated(df11$stationid), ]
  df22<-df2[df2$stationid %in% st,]
  df22 <- df22[!duplicated(df22$stationid), ]
  #setdiff(df11$stationid,df22$stationid)
  
  anom<-df11$gear_temperature-df22$gear_temperature
  #length(anom)
  #length(st)
  anom<-data.frame('stationid'=st,anom)
  df_anom1<-merge(df_anom,anom,by='stationid',all.x=TRUE)
  
  df_anom<-cbind(df_anom,df_anom1$anom)
  
  
}


colnames(df_anom)<-c('stationid',1983:2022)
df_anom[,2:ncol(df_anom)]<-abs(df_anom[,2:ncol(df_anom)])
df_anom$total<-rowMeans(df_anom[,2:ncol(df_anom)],na.rm=TRUE)
df_anom1<-df_anom[order(df_anom$total,decreasing = TRUE),]
st_anom<-data.frame('stationid'=Lon$Group.1,'Lon'=Lon$x,'Lat'=Lat$x)
st_anom$total<-df_anom$total

coordinates(st_anom)<- ~Lon + Lat
proj4string(st_anom) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 

#change CRS
st_anom<-spTransform(st_anom,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))

st_anom<-as.data.frame(st_anom)


ggplot()+
  geom_point(data=st_anom,aes(x=Lon,y=Lat,color=total),size=3)+
  #geom_point(data=st_corners1,aes(x=longitude,y=latitude),color='red',shape=20,size=1)+
  geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
  geom_polygon(data=eez_sh33,aes(x=long,y=lat,group=group),fill=NA,color='grey40')+
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
  geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           xlim = panel_extent$x,
           ylim = panel_extent$y,
           label_axes = "-NE-")+
  scale_color_viridis_c('SBT anomaly')+
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

coordinates(df_var)<- ~Lon + Lat
proj4string(df_var) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 

#change CRS
df_var<-spTransform(df_var,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))

df_var<-as.data.frame(df_var)

ggplot()+
  geom_point(data=df_var,aes(x=Lon,y=Lat,color=x),size=3)+
  #geom_point(data=st_corners1,aes(x=longitude,y=latitude),color='red',shape=20,size=1)+
  geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
  geom_polygon(data=eez_sh33,aes(x=long,y=lat,group=group),fill=NA,color='grey40')+
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
  geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           xlim = panel_extent$x,
           ylim = panel_extent$y,
           label_axes = "-NE-")+
  scale_color_viridis_c('SBT variance')+
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


df_var<-df_var[order(df_var$x,decreasing=TRUE),]
df_var1<-df_var[-(1:18),]

ggplot()+
  geom_point(data=df_var1,aes(x=Lon,y=Lat,color=x),size=3)+
  geom_point(data=st_corners1,aes(x=longitude,y=latitude),color='red',shape=20,size=1)+
  geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
  geom_polygon(data=eez_sh33,aes(x=long,y=lat,group=group),fill=NA,color='grey40')+
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
  geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           xlim = panel_extent$x,
           ylim = panel_extent$y,
           label_axes = "-NE-")+
  scale_color_viridis_c('SBT variance')+
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

install.packages("spThin")

library(spThin)

#sbt ebs
sbt_brick<-coldpool::ebs_bottom_temperature

#names
years<-sort(as.numeric(gsub(".*?([0-9]+).*", "\\1", names(sbt_brick))))

#list to store values
plot_list<-list()

#loop over years, coldpool only has data between 1982 and 2022
for (y in min(years):max(years)) {
  
  y<-2021
  
  #print year to check progress
  cat(paste("    ----- year", y, "-----\n"))  
  
  #2020 empty, so all NAs
  # if (y == 2020) {
  #   r<-sbt_brick[[grep(y-1,names(sbt_brick))]]
  #   r[r]<-NA
  # } else {
    r<-sbt_brick[[grep(y,names(sbt_brick))]]
    r
  #}
  
    #spatial clustering of SBT
    r[r>=-2 & r<=-1]<-'-2'
     r[r>=-1 & r<=0]<-'-1'
     r[r>=0 & r<=+1]<-'0'
     r[r>=+1 & r<=+2]<-'+1'
     r[r>=+2]<-'+2'

coordinates(st_EBS2ebs)<- ~ longitude + latitude
proj4string(r) <- CRS('+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs ')

st_EBS2ebs$temp<-raster::extract(r,st_EBS2ebs)
st_EBS2ebs[which(st_EBS2ebs$temp)]

gplot(r)+
  geom_tile(aes(fill=value))+
  #geom_point(data=xy501,aes(x=longitude,y=latitude),shape=4,size=1)+
  #geom_point(data=st_corners1,aes(x=longitude,y=latitude),color='red',shape=20,size=1)+
  geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
  geom_polygon(data=eez_sh33,aes(x=long,y=lat,group=group),fill=NA,color='grey40')+
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
  geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           xlim = panel_extent$x,
           ylim = panel_extent$y,
           label_axes = "-NE-")+
  scale_fill_binned(limits=c(-2,2),na.value = 'white',low = "#007c9b", high = "#c1f2fe",name='SBT (°C)',guide = guide_colorbar(frame.colour = 'black',ticks.colour = NA))+
  
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



