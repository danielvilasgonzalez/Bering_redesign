####################################################################
####################################################################
##
##    Run sampling optimization based on predicted densities from VAST model
##    and get samples from each strata for each sampling design
##
##    by best stratification, we mean the stratification that ensures the minimum sample cost,
##    sufficient to satisfy precision constraints set on the accuracy of the estimates of the survey target variables Y’s
##    constraints expressed as maximum allowable coefficients of variation in the different domains of interest
##
##    https://cran.r-project.org/web/packages/SamplingStrata/vignettes/SamplingStrata.html
##    (using devtools::install_github(repo = "zoyafuso-NOAA/SamplingStrata"))
##    Daniel Vilas (daniel.vilas@noaa.gov/dvilasg@uw.edu)
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE))
#free up memrory and report the memory usage
gc()

#libraries from cran to call or install/load
pack_cran<-c("splines",'SamplingStrata','wesanderson','dplyr','sp',
             'sf','maptools','parallel','rasterVis','rgeos','scales',
             'rnaturalearth','grid','ggplot2','spatstat','raster')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#install akgfmaps to extract shapefile of Alaska
if (!('akgfmaps' %in% installed.packages())) {
  devtools::install_github('afsc-gap-products/akgfmaps')};library(akgfmaps)

#install coldpool to extract SBT for the EBS
if (!('VAST' %in% installed.packages())) {
  devtools::install_github("james-thorson/VAST@main", INSTALL_opts="--no-staged-install")};library(VAST)

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#setwd
out_dir<-'C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/'
out_dir<-'/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

#version VAST (cpp)
version<-'VAST_v14_0_1'

#selected species
spp<-c('Limanda aspera',
       'Gadus chalcogrammus',
       'Gadus macrocephalus',
       'Atheresthes stomias',
       'Reinhardtius hippoglossoides',
       'Lepidopsetta polyxystra',
       'Hippoglossoides elassodon',
       'Pleuronectes quadrituberculatus',
       'Hippoglossoides robustus',
       'Boreogadus saida',
       'Eleginus gracilis',
       'Anoplopoma fimbria',
       'Chionoecetes opilio',
       'Paralithodes platypus',
       'Paralithodes camtschaticus')

##############################
#OBJECTS FOR PLOTTING
#############################

# Define plot extent (through trial end error) units km
panel_extent <- data.frame(x = c(-1716559.21, -77636.05), #x = c(-1326559.21, -87636.05),
                           y = c(483099.5, 2194909.7)) #y = c(533099.5, 1894909.7))

#Alaska layer
ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")
ak_sppoly<-as(ebs_layers$akland, 'Spatial')

#create directory
dir.create('./shapefiles/',showWarnings = FALSE)

#name shapefiles
shfiles<-c('EBSshelfThorson','NBSThorson','EBSslopeThorson')

#get files from google drive and set up
files<-googledrive::drive_find()
3 #for dvilasg@uw.edu

#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='Shapefiles'),'id']

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)

#loop over shapefiles
for (i in shfiles) {
  
  #i=shfiles[1]
  
  #identify files
  id.data<-files.1[which(grepl(i,files.1$name)),]
  
  #loop over files
  for (j in 1:nrow(id.data)) {
    
    #if not file, download
    if (!(id.data$name[j] %in% list.files('./shapefiles/'))) {
      #download data
      googledrive::drive_download(file=id.data$id[j],
                                  path = paste0('./shapefiles/',id.data$name[j]),
                                  overwrite = TRUE)}
    
  }
  
  #shapefile EBS
  sh<-rgdal::readOGR(dsn='./shapefiles/',layer = i)
  
  #if slope reproject
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
nbs_sh<-NBS_sh
ebs_sh<-EBSshelf_sh

#color palette
pal<-wesanderson::wes_palette('Zissou1',21,type='continuous')

#####################################
# GET EEZ
#####################################

#create directory
dir.create('./shapefiles/',showWarnings = FALSE)

#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='EEZ'),'id']

#list of files and folder
id.data<-googledrive::drive_ls(id.bering.folder$id)

#loop over files
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

#####################################
# Get current ebs and NBS stations
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

#change projection of spatial object
coordinates(st_EBS)<- ~ longitude + latitude

#reproject shapefile
proj4string(st_EBS) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
st_EBS<-spTransform(st_EBS,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))

#to dataframe
st_EBS<-as.data.frame(st_EBS)

###################################
# LOAD GRID EBS (remember to keep the same order as in fit_model if multiple grids)
###################################

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)
#add col and row number
x1<-grid[,c('Lon','Lat','cell')]
names(x1)<-c('x','y','z')
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
colnames(x5)<-c('Lat','Lon','cell','optional','col','row')
grid<-x5[,c('Lat','Lon','cell','col','row')]

###################################
# BASELINE/CURRENT SAMPLING DESIGN
###################################

load('./output/baseline_strata.RData') #baseline_strata

###################################
# SCENARIOS
###################################

#sampling scenarios
samp_df<-expand.grid(strat_var=c('depth + SBT variance','SBT variance','depth'),
                     target_var=c('sumDensity'), #,'sqsumDensity'
                     n_samples=c(520), #c(300,500) 520 (EBS+NBS+CRAB);26 (CRAB); 350 (EBS-CRAB); 494 (NBS-CRAB)
                     n_strata=c(15)) #c(5,10,15)

#add scenario number
samp_df$samp_scn<-paste0(paste0('scn',1:nrow(samp_df)))

#########################
# RUN LOOP SPECIES
#########################



#loop over species
for (sp in spp) {
  
  sp<-'Gadus macrocephalus'
  
  #load optimization data
  load(paste0('./output/species/',sp,'/optimization data/optimization_static_data.RData')) #D6
  #load(paste0('./output/species/',sp,'/projection_data.RData')) #temp_dens_vals
  
  #load fit OM
  load(paste0('./shelf EBS NBS VAST/',sp,'/fit-001.RData'))
  
  #removed cells because of depth
  rem_cells<-D6[which(D6$include==FALSE),'cell']
  ok_cells<-D6[which(D6$include==1),'cell']
  
  #load data_geostat file
  data_geostat<-readRDS(paste0('./data processed/species/',sp,'/','data_geostat_temp.rds')) #fit
  
  ###################################
  # CREATE SPATIAL OBJECT BASED ON CELLS STRATA
  ###################################
  
  plot_list_nsamples<-list()
  plot_list_sampleprop<-list()
  
  for (s in 1:nrow(samp_df)) {
  
      
  #s<-3
  
  #load solutions
  #load(file=paste0('./output/species/',sp,'/optimization data/optimization_results_',samp_df[s,'samp_scn'],'.RData')) #result_list
  
  load(file = paste0("./output/ms_optim_allocations_",samp_df[s,'samp_scn'],".RData")) #all
  
  strata<-rbind(all$result_list$solution$indices,
                data.frame(ID=rem_cells,X1=NA))
  colnames(strata)<-c('cell','Strata')
  
  #n cells by strata
  strata_sum<-aggregate(cell ~ Strata, data = strata, FUN = length)
  names(strata_sum)[2]<-'total_cell'
  #merge
  strata<-merge(strata,all$samples_strata,by.x='Strata',by.y='strata')
  strata<-merge(strata,strata_sum,by='Strata')
  strata$prop<-strata$n_samples/strata$total_cell
  dim(strata)
  
  D8<-merge(D6,strata,by='cell')
  D8$Strata[is.na(D8$Strata)]<-99

  #df to spatialpoint df
  coordinates(D8) <- ~ Lon + Lat
  crs(D8)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  #reproject coordinates for plotting purposes
  D8_1<-spTransform(D8,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
  D8_2<-data.frame(D8_1)
  
  #x and y cells
  xycells<-as.integer(sqrt(dim(D8_1)[1]))
  
  # create a template raster
  r1 <- raster(ext=extent(D8_1),ncol=xycells, nrow=xycells) #c(15800,15800) 7000
  
  #create raster
  r2<-rasterize(D8_1, r1 ,field=c('Strata','n_samples','prop'))
  crs(r2) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
  #plot(r2)
  
  r2[r2==99]<-NA
  
  r3<-as.data.frame(r2,xy=TRUE)
  r4<-rasterToPolygons(r2$Strata,dissolve=TRUE,)
  
  pal <- wes_palette("Zissou1", 15, type = "continuous")
  
  r3<-r3[complete.cases(r3$Strata),] 
  
  r3$prop<-as.factor(r3$prop)
  r3$n_samples<-as.factor(r3$n_samples)
  # Create a named vector for scale_colour_manual
  color_scale <- setNames(as.character(pal), sort(unique(r3$n_samples)))
  
  p1<-
    ggplot()+
    geom_raster(data=r3,aes(x=x,y=y,fill=n_samples))+
    geom_polygon(data=r4,aes(x=long,y=lat,group=group), colour="black", fill=NA)+
    #geom_point(data=D8_2, aes(Lon, Lat, fill=Strata, group=NULL),size=2, stroke=0,shape=21)+
      #scale_fill_gradientn(colours=c("#ea5545", "#f46a9b", "#ef9b20", "#edbf33", "#ede15b", "#bdcf32", "#87bc45", "#27aeef", "#b33dc6"))+
      scale_fill_manual(values = color_scale)+ 
      #scale_fill_gradientn(colors=pal,values = sort(unique(r3$prop)))+ 
      scale_x_continuous(expand = c(0,0),position = 'bottom',
                         breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
      coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
               xlim = panel_extent$x-c(000000,100000),
               ylim = panel_extent$y-c(0,300000),
               label_axes = "-NE-")+
      theme(aspect.ratio = 1,panel.grid.major = element_blank(),
            panel.background = element_rect(fill = NA),panel.ontop = TRUE,
            legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(20, 'points'),
            legend.key.width= unit(20, 'points'),axis.title = element_blank(),legend.position = 'none',
            panel.border = element_rect(fill = NA, colour = NA),legend.key = element_rect(color="black"),
            legend.spacing.y = unit(8, 'points'),
            axis.text=element_blank(),axis.ticks = element_blank(),
            plot.margin = margin(0.01,0.01,0.01,0.01), 
            axis.ticks.length = unit(-5,"points"),plot.title = element_text(size=18,hjust = 0.5,vjust=-5, face="bold"))+
      scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())+
      guides(fill = guide_legend(order=2,override.aes=list(shape = 22,size=8)),
             color = guide_legend(order=1,override.aes=list(size=8)),
             shape = guide_legend(order=1),override.aes=list(size=8))+
    labs(title=paste0(gsub('_',' + ',samp_df[s,'strat_var'])),fill='')
    

    
    plot_list_nsamples[[s]]<-p1
    
    color_scale <- setNames(as.character(pal), sort(unique(r3$prop)))
    
    p2<-
      ggplot()+
      geom_raster(data=r3,aes(x=x,y=y,fill=prop))+
      geom_polygon(data=r4,aes(x=long,y=lat,group=group), colour="black", fill=NA)+
      #geom_point(data=D8_2, aes(Lon, Lat, fill=Strata, group=NULL),size=2, stroke=0,shape=21)+
      #scale_fill_gradientn(colours=c("#ea5545", "#f46a9b", "#ef9b20", "#edbf33", "#ede15b", "#bdcf32", "#87bc45", "#27aeef", "#b33dc6"))+
      scale_fill_manual(values = color_scale)+ 
      #scale_fill_gradientn(colors=pal,values = sort(unique(r3$prop)))+ 
      scale_x_continuous(expand = c(0,0),position = 'bottom',
                         breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
      coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
               xlim = panel_extent$x-c(000000,100000),
               ylim = panel_extent$y-c(0,300000),
               label_axes = "-NE-")+
      theme(aspect.ratio = 1,panel.grid.major = element_blank(),
            panel.background = element_rect(fill = NA),panel.ontop = TRUE,
            legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(20, 'points'),
            legend.key.width= unit(20, 'points'),axis.title = element_blank(),legend.position = 'none',
            panel.border = element_rect(fill = NA, colour = NA),legend.key = element_rect(color="black"),
            legend.spacing.y = unit(8, 'points'),
            axis.text=element_blank(),axis.ticks = element_blank(),
            plot.margin = margin(0.01,0.01,0.01,0.01), 
            axis.ticks.length = unit(-5,"points"),plot.title = element_text(size=18,hjust = 0.5,vjust=-5, face="bold"))+
      scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())+
      guides(fill = guide_legend(order=2,override.aes=list(shape = 22,size=8)),
             color = guide_legend(order=1,override.aes=list(size=8)),
             shape = guide_legend(order=1),override.aes=list(size=8))+
      labs(title=paste0(gsub('_',' + ',samp_df[s,'strat_var'])),fill='')  


    plot_list_sampleprop[[s]]<-p2

    
  
  #save plot
# ragg::agg_png(paste0('./figures/species/',sp,'/optimized_stratification_',samp_df[s,'samp_scn'],'.png'), width = 7, height = 7, units = "in", res = 300)
# print(p)
# dev.off()
  }
}
#}


##################3baseline

load('./output/baseline_strata.RData')

baseline_strata$cell_strata<-as.data.frame(baseline_strata$cell_strata)
baseline_strata$cell_strata<-merge(baseline_strata$cell_strata,baseline_strata$n_samples,by.x='Stratum',by.y='stratum')

#n cells by strata
strata_sum<-aggregate(cell ~ Stratum, data = baseline_strata$cell_strata, FUN = length)
names(strata_sum)[2]<-'total_cell'
#merge
baseline_strata$cell_strata<-merge(baseline_strata$cell_strata,strata_sum,by='Stratum')
baseline_strata$cell_strata$prop<-baseline_strata$cell_strata$scnbase/baseline_strata$cell_strata$total_cell
names(baseline_strata$cell_strata)[7]<-'n_samples'

#df to spatialpoint df
coordinates(baseline_strata$cell_strata) <- ~ Lon + Lat
crs(baseline_strata$cell_strata)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#reproject coordinates for plotting purposes
D8_1<-spTransform(baseline_strata$cell_strata,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
D8_2<-data.frame(D8_1)

#x and y cells
xycells<-as.integer(sqrt(dim(D8_1)[1]))

# create a template raster
r1 <- raster(ext=extent(D8_1),ncol=xycells, nrow=xycells) #c(15800,15800) 7000

#create raster
r2<-rasterize(D8_1, r1 ,field=c('Stratum','n_samples','prop'))
crs(r2) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
#plot(r2)

r2[r2==99]<-NA
r3<-as.data.frame(r2,xy=TRUE)
r4<-rasterToPolygons(r2,dissolve=TRUE)

#aggregate(r3$layer,by=list(r3$layer),FUN=length)
# 
# p<-
#   ggplot()+
#   geom_raster(data=r3,aes(x=x,y=y,fill=as.character(layer)))+
#   geom_polygon(data=r4,aes(x=long,y=lat,group=group), colour="black", fill=NA)+
#   scale_fill_gradientn(colours=c("#ea5545", "#f46a9b", "#ef9b20", "#edbf33", "#ede15b", "#bdcf32", "#87bc45", "#27aeef", "#b33dc6"),
#                            guide = guide_legend(frame.colour = "black", 
#                                                 ticks.colour = "black"),breaks=sort(unique(D8_2$Strata)),na.value = "transparent",labels=paste0(sort(unique(D8_2$Strata))," (n=",result_list$sample_allocations,')'))+
#       scale_x_continuous(expand = c(0,0),position = 'bottom',
#                          breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
#       coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
#                xlim = panel_extent$x-c(100000,100000),
#                ylim = panel_extent$y-c(0,300000),
#                label_axes = "-NE-")+
#       theme(aspect.ratio = 1,panel.grid.major = element_blank(),
#             panel.background = element_rect(fill = NA),panel.ontop = TRUE,
#             legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(20, 'points'),
#             legend.key.width= unit(20, 'points'),axis.title = element_blank(),legend.position = 'none',
#             panel.border = element_rect(fill = NA, colour = NA),legend.key = element_rect(color="black"),
#             legend.spacing.y = unit(8, 'points'),
#            axis.text=element_blank(),axis.ticks = element_blank(),
#             plot.margin = margin(0.01,0.01,0.01,0.01), 
#             axis.ticks.length = unit(-5,"points"),plot.title = element_text(size=14,hjust = 0.5,vjust=-10, face="bold"))+
#       scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())+
#       guides(fill = guide_legend(order=2,override.aes=list(shape = 22,size=8)),
#              color = guide_legend(order=1,override.aes=list(size=8)),
#              shape = guide_legend(order=1),override.aes=list(size=8))+
#       labs(title=paste0(gsub('_',' + ',samp_df[s,'strat_var'])),fill='')


r3<-r3[complete.cases(r3$Stratum),] 
r3$Stratum<-as.character(r3$Stratum)
#as.vector(unique(g$data[[1]]["fill"]))

r3$prop<-as.factor(r3$prop)
r3$n_samples<-as.factor(r3$n_samples)
# Create a named vector for scale_colour_manual
color_scale <- setNames(as.character(pal), sort(unique(r3$prop)))

p2<-
  ggplot()+
  geom_raster(data=r3,aes(x=x,y=y,fill=prop))+
  geom_polygon(data=r4,aes(x=long,y=lat,group=group), colour="black", fill=NA)+
  # scale_fill_manual(values=c('82'="#86BA59",'81'= "#4CA5EB",'90'= "#87BC45",'71'= "#629CE7",'70'= "#B33DC6",
  #                            '61'= "#B8CD34",'41'= "#EDC237",'62'= "#B3CB37",'43'= "#EDC940",'20'= "#F46A9B",
  #                            '10'= "#EA5545",'42'= "#EDC63C",'31'= "#EF9F22",'32'= "#EFA224",'50'= "#EDE15B"))+
  scale_fill_manual(values = color_scale)+ 
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           xlim = panel_extent$x-c(000000,100000),
           ylim = panel_extent$y-c(0,300000),
           label_axes = "-NE-")+
  theme(aspect.ratio = 1,panel.grid.major = element_blank(),
        panel.background = element_rect(fill = NA),panel.ontop = TRUE,
        legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(20, 'points'),
        legend.key.width= unit(20, 'points'),axis.title = element_blank(),
        panel.border = element_rect(fill = NA, colour = NA),legend.key = element_rect(color="black"),
        legend.spacing.y = unit(8, 'points'),
        axis.text=element_blank(),axis.ticks = element_blank(),
        plot.margin = margin(0.01,0.01,0.01,0.01), legend.position = 'none',
        axis.ticks.length = unit(-5,"points"),plot.title = element_text(size=18,hjust = 0.5,vjust=-5, face="bold"))+
  scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())+
  guides(fill = guide_legend(order=2,override.aes=list(shape = 22,size=8)),
         color = guide_legend(order=1,override.aes=list(size=8)),
         shape = guide_legend(order=1),override.aes=list(size=8))+
  labs(title='baseline',fill='')

pal1 <- wes_palette("Zissou1", length(unique(r3$n_samples)), type = "continuous")

color_scale <- setNames(as.character(pal1), sort(unique(r3$n_samples)))

p1<-
  ggplot()+
  geom_raster(data=r3,aes(x=x,y=y,fill=n_samples))+
  geom_polygon(data=r4,aes(x=long,y=lat,group=group), colour="black", fill=NA)+
  # scale_fill_manual(values=c('82'="#86BA59",'81'= "#4CA5EB",'90'= "#87BC45",'71'= "#629CE7",'70'= "#B33DC6",
  #                            '61'= "#B8CD34",'41'= "#EDC237",'62'= "#B3CB37",'43'= "#EDC940",'20'= "#F46A9B",
  #                            '10'= "#EA5545",'42'= "#EDC63C",'31'= "#EF9F22",'32'= "#EFA224",'50'= "#EDE15B"))+
  scale_fill_manual(values = color_scale)+ 
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           xlim = panel_extent$x-c(000000,100000),
           ylim = panel_extent$y-c(0,300000),
           label_axes = "-NE-")+
  theme(aspect.ratio = 1,panel.grid.major = element_blank(),
        panel.background = element_rect(fill = NA),panel.ontop = TRUE,
        legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(20, 'points'),
        legend.key.width= unit(20, 'points'),axis.title = element_blank(),
        panel.border = element_rect(fill = NA, colour = NA),legend.key = element_rect(color="black"),
        legend.spacing.y = unit(8, 'points'),
        axis.text=element_blank(),axis.ticks = element_blank(),
        plot.margin = margin(0.01,0.01,0.01,0.01), legend.position = 'none',
        axis.ticks.length = unit(-5,"points"),plot.title = element_text(size=18,hjust = 0.5,vjust=-5, face="bold"))+
  scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())+
  guides(fill = guide_legend(order=2,override.aes=list(shape = 22,size=8)),
         color = guide_legend(order=1,override.aes=list(size=8)),
         shape = guide_legend(order=1),override.aes=list(size=8))+
  labs(title='baseline',fill='')

  g <- ggplot_build(p1)
  unique(g$data[[1]]["fill"])
  
  
plot_list_nsamples[['baseline']]<-p1
plot_list_sampleprop[['baseline']]<-p2

cowplot::plot_grid(plot_list_nsamples[['baseline']],plot_list_nsamples[[3]],plot_list_nsamples[[2]],plot_list_nsamples[[1]])






  ragg::agg_png(paste0('./figures/sampling designs nsamples.png'), width = 10, height = 10, units = "in", res = 300)
  cowplot::plot_grid(plot_list_nsamples[['baseline']],plot_list_nsamples[[3]],plot_list_nsamples[[2]],plot_list_nsamples[[1]])
  dev.off()

  ragg::agg_png(paste0('./figures/sampling designs propsamples.png'), width = 10, height = 10, units = "in", res = 300)
  cowplot::plot_grid(plot_list_sampleprop[['baseline']],plot_list_sampleprop[[3]],plot_list_sampleprop[[2]],plot_list_sampleprop[[1]])
  dev.off()
  
    