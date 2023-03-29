####################################################################
####################################################################
##
##    Run sampling optimization based on predicted densities from VAST model
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
pack_cran<-c("splines",'SamplingStrata','wesanderson','dplyr','sp','rasterVis','rgeos','scales','rnaturalearth','grid','ggplot2')

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
setwd(out_dir)

#version VAST (cpp)
version<-'VAST_v13_1_0'

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
1 #for dvilasg@uw.edu

#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='Shapefiles'),'id']

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)

#loop over shapefiles
for (i in shfiles) {
  
  #i=shfiles[1]
  
  id.data<-files.1[which(grepl(i,files.1$name)),]
  
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
# ALL EBS+NBS points
#########################

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

#load grid data
#https://github.com/James-Thorson-NOAA/FishStatsUtils/tree/main/data
#https://github.com/danielvilasgonzalez/Bering_redesign/blob/main/Scripts/04_Bering10K_data.R
load('./data processed/lastversion_grid_EBS.RData') #grid.ebs_year

#remove slope grid
grid.ebs_year1<-grid.ebs_year[which(grid.ebs_year$region!='EBSslope'),]

###################################
# SCENARIOS
###################################

df_scn<-expand.grid(strat_var=c('Lat_LonE','Lat_Depth','Lat_varTemp','Lat_meanTempF','Depth_meanTempF','Depth_varTemp'),
                    target_var=c('sumDensity'), #,'sqsumDensity'
                    n_samples=c(300,500))

###################################
# LOAD FIT OBJECT (from VAST::fit_model()) or VAST::project_model()
###################################

#load fit file
load('./shelf EBS NBS VAST/Gadus macrocephalus/temp3d/b0_19822022fit.RData')

#project model example
#p1<-project_model(x = fit,n_proj = n_proj)

#################################################
# SIMULATE DATA
#################################################

#reload DLL file from model
#fit<-reload_model(fit)

#simulate data
#sim.data<-simulate_data(fit,
#                        type=1) #1 measurement error or conditional #2 unconditional #3 new fixed and random #4 random effects from MLE

#################################################
# ARRANGE SIMULATED or PREDICTED DATA DATA
#################################################

## Years to use
year_set <- as.numeric(fit$year_labels)
years_included <- 1:length(year_set)
n_years <- length(years_included)

## Scientific and common names used in optimization
spp<-'Gadus macrocephalus'
sp<-'Gadus macrocephalus'

#get predictions for category 1
temp_dens_vals <- array(NA,
                        dim = c(nrow(fit$extrapolation_list$Data_Extrap),length(spp),length(years_included)),
                        dimnames = list(1:nrow(fit$extrapolation_list$Data_Extrap),spp,years_included))
temp_dens_vals[, sp, ] <- fit$Report$D_gct[, 1, years_included]

density_input<-temp_dens_vals

D_gt<-fit$Report$D_gct[,1,]
#dim(D_gt)
#D_gt_proj<-D_gt[,paste0(project_yrs)]

#drop units
D_gt<-drop_units(D_gt)

#dataframe of cells with predictions
D_gt<-data.frame('cell'=c(1:fit$spatial_list$n_g),D_gt)

#remove duplicates just in case
D_gt<-D_gt[!duplicated(as.list(D_gt))]

#years simulated
yrs<-as.numeric(fit$year_labels)

#rename years predictions 
colnames(D_gt)<-c('cell',yrs) #,project_yrs

#reshape
D_gt1<-reshape2::melt(D_gt,id=c('cell'))

#get map info
mdl <- make_map_info(Region = fit$settings$Region,
                     spatial_list = fit$spatial_list,
                     Extrapolation_List = fit$extrapolation_list)

#merge cells and predictions
D <- merge(D_gt1, mdl$PlotDF, by.x='cell', by.y='x2i')

#rename columns
colnames(D)<-c('cell','Year','Density','Lat','Lon','Include')

#get grid data for yrs in the simulation
grid.ebs_year2<-grid.ebs_year1[which(grid.ebs_year1$Year %in% yrs),]

#merge grid data and predictions
D1<-merge(D,grid.ebs_year2,by=c('Lat','Lon','Year'))

#static sampling so, we want to aggregate annual predictions: mean density, mean temp, and temp var
D2<-aggregate(cbind(Temp,Density) ~ Lat+Lon+cell+Depth, data = D1, FUN = mean, na.rm = TRUE)
D3<-aggregate(cbind(Temp,Density) ~ Lat+Lon+cell+Depth, data = D1, FUN = var, na.rm = TRUE)
D4<-aggregate(cbind(Density) ~ Lat+Lon+cell+Depth, data = D1, FUN = sum, na.rm = TRUE)
colnames(D2)[5:6]<-paste0('mean',colnames(D2)[5:6])
colnames(D3)[5:6]<-paste0('var',colnames(D3)[5:6])
colnames(D4)[5]<-paste0('sum',colnames(D4)[5])

#merge aggregate values
D5<-merge(D2,D3,by=c('Lat','Lon','cell','Depth'))
D6<-merge(D5,D4,by=c('Lat','Lon','cell','Depth'))

#add squared density
D6$sqsumDensity<-(D6$sumDensity)^2

#duplicate object to remove 
D7<-D6[which(D6$Depth>0),]

#cells including positive depth
cells<-unique(D7$cell)

#convert SBT into F to get positive values only
D7$meanTempF<-(9/5)*D7$meanTemp + 32

#add longitude on eastings to get positive values
D7$LonE<-D7$Lon+180+180

#removed cells because of depth<0
rem_cells<-D6[which(D6$Depth<0),'cell']

#################################################
# STRATIFICATION OPTIMIZATION SETTINGS
#################################################

#n_years
n_years<-length(unique(D$Year))

#number of domains
n_dom<-1

#number of cells
n_cells<-length(cells)

#domain_input
domain_input<-rep(1, n_cells)

#create a list to store results
l<-list()
plot_l<-list()

#########################
# RUN LOOP SCENARIOS
#########################

#loop through scenarios
for (scn in 1:nrow(df_scn)) {
  
  #scn<-1
  
  #print species to check progress
  cat(paste(" #############  Scenario number", scn, " #############\n"))
  
  #stratification variables 
  stratum_var_input<-data.frame(X1 = D7[,paste0(sub("\\_.*", "", df_scn[scn,'strat_var']))],
                                X2 = D7[,paste0(sub(".*_", "", df_scn[scn,'strat_var']))]) #Xspp #set different scenarios and spp ############ TO CHECK
  
  #target variables
  target_var_input<-data.frame(Y1 = D7$sumDensity,
                               Y1_SQ_SUM = D7$sqsumDensity) #D7$sqsumDensity #Ynspp #set different scenarios and spp ############ TO CHECK
  
  #weights
  #in case add weights based on observed years
  
  #create df
  frame <- data.frame(domainvalue = domain_input,
                      id = cells,
                      stratum_var_input,
                      WEIGHT=n_years,
                      target_var_input) 
  
  ###################################
  # SIMPLE RANDOM SAMPLING CV CONSTRAINTS
  ###################################
  
  #Initiate CVs to be those calculated under simple random sampling (SRS)
  srs_stats <- SamplingStrata::buildStrataDF(dataset = cbind( frame[, -grep(x = names(frame), pattern = "X")],X1 = 1))
  
  #number of samples 
  #520 in 2022 (NBS and EBSshelf)
  #376 in 2018 (EBSshelf)
  n_samples <- df_scn[scn,'n_samples']
  
  #number of samples
  srs_n <- as.numeric(n_samples * table(frame$domainvalue) / n_cells)
  
  ## SRS statistics
  srs_var <- srs_stats$S1^2 * (1 - srs_n / n_cells) / srs_n
  srs_cv <- sqrt(srs_var) / srs_stats$M1
  
  #create cv object for constraints
  cv <- list()
  cv[["CV1"]] <- srs_cv
  cv[["DOM"]] <- 1:n_dom
  cv[["domainvalue"]] <- 1:n_dom
  cv <- as.data.frame(cv)
  
  ###################################
  # STRATAS
  ###################################
  
  #number of stratas 
  no_strata<-10
  
  #get n_strata from kmean suggestion
  kmean<-KmeansSolution2(frame=frame,
                         errors=cv,
                         maxclusters = 20)
  
  #number strata from kmean
  no_strata<-tapply(kmean$suggestions,
                    kmean$domainvalue,
                    FUN=function(x) length(unique(x)))
  
  ###################################
  # RUN OPTIMIZATION
  ###################################
  
  #run optimization
  solution <- optimStrata(method = "continuous",
                          errors = cv, 
                          framesamp = frame,
                          iter = 50, #300
                          pops = 10, #100
                          elitism_rate = 0.1,
                          mut_chance = 1 / (no_strata[1] + 1),
                          nStrata = no_strata,
                          showPlot = T,
                          writeFiles = T)
  
  ###################################
  # STORE SOLUTIONS
  ###################################
  
  #results
  framenew<-solution$framenew
  outstrata<-solution$aggr_strata
  ss<-summaryStrata(framenew,outstrata)
  #head(ss)
  
  #plot strata 2D
  #plotStrata2d(framenew,outstrata,domain=1,vars=c('X1','X2'),labels = c('VarTemp','Depth'))
  
  ## Organize result outputs
  solution$aggr_strata$STRATO <- as.integer(solution$aggr_strata$STRATO)
  solution$aggr_strata <- 
    solution$aggr_strata[order(solution$aggr_strata$DOM1,
                               solution$aggr_strata$STRATO), ]
  
  sum_stats <- summaryStrata(solution$framenew,
                             solution$aggr_strata,
                             progress=FALSE)
  sum_stats$stratum_id <- 1:nrow(sum_stats)
  sum_stats$Population <- sum_stats$Population / n_years
  sum_stats$wh <- sum_stats$Allocation / sum_stats$Population
  sum_stats$Wh <- sum_stats$Population / n_cells
  sum_stats <- cbind(sum_stats,
                     subset(x = solution$aggr_strata,
                            select = -c(STRATO, N, COST, CENS, DOM1, X1)))
  
  #store summary results
  l[[scn]]<-sum_stats
  
  ###################################
  # CREATE SPATIAL OBJECT BASED ON CELLS STRATA
  ###################################
  
  strata<-rbind(solution$indices,
                data.frame(ID=rem_cells,X1=NA))
  colnames(strata)<-c('cell','Strata')
  
  dim(strata)
  
  D8<-merge(D6,strata,by='cell')
  
  #df to spatialpoint df
  coordinates(D8) <- ~ Lon + Lat
  crs(D8)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  #reproject coordinates for plotting purposes
  D8_1<-spTransform(D8,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
  D8_2<-data.frame(D8_1)
  
  ###################################
  # MULTIVARIATE OPTIMAL ALLOCATION
  ###################################
  
  #create df to extract
  temp_stratif <- solution$aggr_strata
  temp_stratif$N <- temp_stratif$N / length(yrs)
  temp_stratif$DOM1 <- 1
  
  #run multivariate allocation
  temp_bethel <- SamplingStrata::bethel(
    errors = cv,
    stratif = temp_stratif, 
    realAllocation = T, 
    printa = T)
  temp_n <- sum(ceiling(temp_bethel))
  
  #number of samples per strata
  allocations<-as.integer(temp_bethel)
  
  #strata per cell
  temp_ids<-solution$indices
  
  sample_vec <- c()
  
  #random sample for each strata
  for(istrata in 1:length(allocations)) {
    sample_vec <- c(sample_vec,
                    sample(x = temp_ids[which(temp_ids$X1==istrata),'ID'], #which(temp_ids == istrata)
                           size = allocations[istrata]) )
  }
  
  #points dataframe
  points<-data.frame(cell=sample_vec,
                     strata=rep(1:length(temp_bethel),allocations))
  
  #merge
  points1<-merge(points,D7,by='cell',all.x=TRUE)
  
  #change projection of spatial object 
  coordinates(points1)<- ~ Lon + Lat
  
  #reproject shapefile
  proj4string(points1) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
  points1<-spTransform(points1,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))
  
  #to dataframe
  points1<-as.data.frame(points1)
  
  #########################
  # JOIN POINTS FOR LEGEND PURPOSES
  #########################
  
  df<-rbind(data.frame(Lat=points1$Lat,Lon=points1$Lon,Stations='optimization'),
            data.frame(Lat=st_EBS$latitude,Lon=st_EBS$longitude,Stations='current design'),
            data.frame(Lat=st_corners1$latitude,Lon=st_corners1$longitude,Stations='corner crab'))
  
  #########################
  # MAP POINTS
  #########################

  p<-
    ggplot()+
        geom_point(data=D8_2, aes(Lon, Lat, fill=Strata, group=NULL),size=2, stroke=0,shape=21)+
        scale_fill_gradientn(colours=c("#ea5545", "#f46a9b", "#ef9b20", "#edbf33", "#ede15b", "#bdcf32", "#87bc45", "#27aeef", "#b33dc6"),
                             guide = guide_legend(),breaks=sort(unique(D8_2$Strata)),labels=paste0(sort(unique(D8_2$Strata))," (n=",allocations,')'))+
        geom_point(data=df,aes(x=Lon,y=Lat,color=Stations,shape=Stations),fill='white',color='black',size=2)+
        scale_shape_manual(values = c('optimization'=21,
                                      'current design'=4,
                                      'corner crab'=8),
                           breaks=unique(df$Stations),
                           labels=paste0(unique(df$Stations)," (n=",c(nrow(points1),nrow(st_EBS),nrow(st_corners1)),')'))+
        #geom_point(data=st_EBS,aes(x=longitude,y=latitude),shape=4,size=1)+
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
        theme(aspect.ratio = 1,panel.grid.major = element_line(color = rgb(0, 0, 0,20, maxColorValue = 285), linetype = 'dashed', linewidth =  0.5),
              panel.background = element_rect(fill = NA),panel.ontop = TRUE,text = element_text(size=14),
              legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(20, 'points'),
              legend.key.width= unit(20, 'points'),axis.title = element_blank(),legend.position = c(0.12,0.47),
              panel.border = element_rect(fill = NA, colour = 'black'),legend.key = element_rect(color="black"),
              axis.text = element_text(color='black'),legend.spacing.y = unit(8, 'points'),
              axis.text.y.right = element_text(hjust= 0.1 ,margin = margin(-7,0,0,-30, unit = 'points'),color='black'),
              axis.text.x = element_text(vjust = 6, margin = margin(-7,0,0,-30, unit = 'points'),color='black'),
              axis.ticks.length = unit(-5,"points"),plot.title = element_text(size=12,vjust = -18, hjust=0.95,face="bold"),
              plot.margin=margin(c(-25,0,0,-10)))+
        annotate("text", x = -256559, y = 1354909, label = "Alaska",parse=TRUE,size=7)+
        annotate("text", x = -1296559, y = 2049090, label = "Russia",parse=TRUE,size=7)+
        scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())+
        annotate("text", x = -1296559, y = 744900, label = "italic('Bering Sea')",parse=TRUE,size=9)+
        guides(fill = guide_legend(order=2,override.aes=list(shape = 22,size=8)),
               color = guide_legend(order=1,override.aes=list(size=8)),
               shape = guide_legend(order=1),override.aes=list(size=8))+
        labs(title=paste0('Scenario\n',df_scn[scn,'strat_var'],' n=',df_scn[scn,'n_samples']))
  
  #save image
  print(p)
  ggsave(paste0('./figures/optimization_',sp,"_",scn,'.png'), width = 9, height = 9)
  
  #store into the list
  plot_l[[scn]]<-p
  
}

#save results list
save(l,file=paste0('./shelf EBS NBS VAST/Gadus macrocephalus/temp3d/optimization_summary_',sp,'.RData'))

#save plot list
save(plot_l,file=paste0('./shelf EBS NBS VAST/Gadus macrocephalus/temp3d/optimization_plot_',sp,'.RData'))
