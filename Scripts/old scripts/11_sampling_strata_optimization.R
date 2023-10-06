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
             'rnaturalearth','grid','ggplot2','spatstat')

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
out_dir<-'/Users/daniel/Work/Adapting to a Changing Seascape/'

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
2 #for dvilasg@uw.edu

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
samp_df<-expand.grid(strat_var=c('Depth_varTemp','varTemp','Depth'),
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
  load(paste0('./shelf EBS NBS VAST/',sp,'/fit.RData'))
  
  #removed cells because of depth
  rem_cells<-D6[which(D6$include==FALSE),'cell']
  ok_cells<-D6[which(D6$include==1),'cell']
    
  #load data_geostat file
  data_geostat<-readRDS(paste0('./data processed/species/',sp,'/','data_geostat_temp.rds')) #fit
  
  #n cells
  n_cells<-length(ok_cells)
  
  #how manyt yrs were encounter this species
  xall<-data_geostat %>% count(year)
  data_geostat[which(is.na(data_geostat$cpue_kgha)),'cpue_kgha']<-0
  x0<-subset(data_geostat,cpue_kgha==0) %>% count(year)
  xpct<-x0[,2]/xall[,2]*100
  n_years<-length(xpct)-length(0 %in% xpct)
  
  #number of domains
  n_dom<-1
  
  #domain_input
  domain_input<-rep(1, n_cells)

  #df summary
  sp_sum_stats<-data.frame(matrix(nrow = 0,ncol=18))
  names(sp_sum_stats)<- c("Domain","Stratum","Population","Allocation","SamplingRate","Lower_X1","Upper_X1","Lower_X2","Upper_X2",
                           "stratum_id","wh","Wh","M1","S1","SOLUZ","samp_scn","sp")

  #subset cells with appropiate depth
  static_df1<-subset(D6,cell %in% ok_cells)

  #########################
  # sampling designs
  #########################

  #loop through sampling designs
  for (s in 1:nrow(samp_df)) {
      
    #s<-1
      
    #print scenario to check progress
    cat(paste(" #############   Species", sp, match(sp,spp), 'out of',length(spp),  "  #############\n",
              #" #############  SBT Scenario", sbt_scn, " #############\n",
              " #############  Sampling Scenario", samp_df[s,"samp_scn"], " #############\n"))
    
    #if scenario includes two stratifying factors
    if (grepl('_',samp_df[s,'strat_var'])) {
      #stratification variables 
      stratum_var_input<-data.frame(X1 = static_df1[,paste0(sub("\\_.*", "", samp_df[s,'strat_var']))],
                                      X2 = static_df1[,paste0(sub(".*_", "", samp_df[s,'strat_var']))]) 
    } else {
      stratum_var_input<-data.frame(X1 = static_df1[,paste0(sub("\\_.*", "", samp_df[s,'strat_var']))])
    }

    #target variables
    target_var_input<-data.frame(Y1 = static_df1$sumDensity,
                                 Y1_SQ_SUM = static_df1$sqsumDensity) #D7$sqsumDensity #Ynspp #set different scenarios and spp ############ TO CHECK
      
    #create df
    frame <- data.frame(domainvalue = domain_input, #domain
                        id = static_df1$cell, #id as cells
                        stratum_var_input, #Stratification variables
                        WEIGHT=n_years, #weight for spp depending on the presence of years
                        target_var_input) #target variables 
      
    ###################################
    # Simple random sampling design CV constraints
    ###################################
      
    #Initiate CVs to be those calculated under simple random sampling (SRS)
    srs_stats <- SamplingStrata::buildStrataDF(dataset = cbind( frame[, -grep(x = names(frame), pattern = "X")],X1 = 1))
      
    #number of samples (maximum)
    #520 in 2022 (NBS and EBSshelf)
    #376 in 2018 (EBSshelf)
    n_samples <- samp_df[s,'n_samples']
    
    #number of samples
    srs_n <- as.numeric(n_samples * table(frame$domainvalue) / n_cells)
      
    #SRS statistics
    srs_var <- srs_stats$S1^2 * (1 - srs_n / n_cells) / srs_n #M1<- mean; S1<-SD
    srs_cv <- sqrt(srs_var) / srs_stats$M1
      
    #create cv object for constraints
    cv <- list()
    cv[["CV1"]] <- srs_cv
    cv[["DOM"]] <- 1:n_dom
    cv[["domainvalue"]] <- 1:n_dom
    cv <- as.data.frame(cv)
      
    ###################################
    # Strata
    ###################################
      
    #number of stratas 
    no_strata<-samp_df[s,'n_strata']
      
    ###################################
    # Run optimization
    ###################################
      
    #run optimization
    solution <- optimStrata(method = "continuous", #continous variables
                            errors = cv,  #precision level - maximum allowable coefficient of variation set by the simple random sampling 
                            framesamp = frame, #df of input variables
                            #iter = 300, #300 #aximum number of iterations
                            #pops = 100, #100  #dimension of each generations
                            elitism_rate = 0.1, #0.1
                            mut_chance = 1 / (no_strata[1] + 1), #mutation chance
                            nStrata = no_strata, #maximum strata
                            showPlot = TRUE, #FALSE
                            writeFiles = FALSE)
      
    ###################################
    # Store solutions from optimizations
    ###################################
      
    #results from optimization
    framenew<-solution$framenew
    outstrata<-solution$aggr_strata
    ss<-summaryStrata(framenew,outstrata)
      
    #organize result outputs
    solution$aggr_strata$STRATO <- as.integer(solution$aggr_strata$STRATO)
    solution$aggr_strata <- 
    solution$aggr_strata[order(solution$aggr_strata$DOM1,
                               solution$aggr_strata$STRATO), ]
    
    #save optimization stats  
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
      
    #add scn and sp
    sum_stats$samp_scn<-samp_df[s,'samp_scn']
    sum_stats$sp<-sp
    
    #if one stratifying factor add columns  
    if (!grepl('_',samp_df[s,'strat_var'])) {
        
      sum_stats<-data.frame(sum_stats[,c("Domain","Stratum","Population","Allocation","SamplingRate","Lower_X1","Upper_X1")],
                           "Lower_X2"=NA,"Upper_X2"=NA,
                           sum_stats[,c("stratum_id","wh","Wh","M1",'M2',"S1",'S2',"SOLUZ","samp_scn","sp")])
    }
    
    #append stat results  
    sp_sum_stats<-rbind(sp_sum_stats,sum_stats)
      
    ###################################
    # Create spatial objects
    ###################################

    #strata by cell  
    strata<-solution$indices
    colnames(strata)<-c('cell','Strata')
      
    #add a strata value to each cell
    D8<-merge(D6,strata,by='cell',all.x=TRUE)
    D8<-D8[,c("cell","Lat","Lon","Strata")]
    D8$Strata<-as.numeric(D8$Strata)
    D8$Strata<-ifelse(is.na(D8$Strata),999,D8$Strata)
      
    #df to spatialpoint df
    coordinates(D8) <- ~ Lon + Lat
    crs(D8)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
      
    #reproject coordinates for plotting purposes
    D8_1<-spTransform(D8,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
    D8_2<-merge(as.data.frame(D8_1),x5,by='cell')
    D8_2<-D8_2[,c('cell','Strata','Lon.x','Lat.x','col','row')]
    names(D8_2)[c(2,3,4)]<-c('strata','Lon','Lat')
    
    #x and y cells
    xycells<-as.integer(sqrt(dim(D8_1)[1]))
      
    # create a template raster
    r1 <- raster(ext=extent(D8_1),ncol=xycells, nrow=xycells) #c(15800,15800) 7000
      
    #create raster
    r2<-rasterize(D8_1, r1 ,field='Strata')
    crs(r2) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
    r2[r2==999] <- NA
    #plot(r2)
      
    ###################################
    # Multimariate optimal allocation
    ###################################
      
    #df for optimization
    temp_frame <- frame
    temp_frame$domainvalue <- n_dom
    srs_stats <- SamplingStrata::buildStrataDF(
    dataset = cbind( temp_frame[, -grep(x = names(temp_frame), 
                                          pattern = "X")],
                       X1 = 1))
    srs_n <- as.numeric(n_samples * table(temp_frame$domainvalue) / n_cells)
    
    #srs statistics
    srs_var <- srs_stats$S1^2 * (1 - srs_n / n_cells) / srs_n
    srs_cv <- sqrt(srs_var) / srs_stats$M1
    
    #error
    error_df <- data.frame("DOM" = "DOM1",
                           "CV1" = srs_cv,
                           "domainvalue"  = 1)
    
    #df to run bethel algorithm
    temp_stratif <- solution$aggr_strata
    temp_stratif$N <- temp_stratif$N / n_years
    temp_stratif$DOM1 <- 1
    
    #run multivariate allocation
    temp_bethel <- SamplingStrata::bethel(errors = error_df,
                                          stratif = temp_stratif,
                                          realAllocation = T,
                                          printa = T)
      
      
    temp_n <- sum(ceiling(temp_bethel))
    error_df$CV1 <- as.numeric(attributes(temp_bethel)$outcv[, "ACTUAL CV"])
      
    #run multivariate allocation
    temp_bethel <- SamplingStrata::bethel(stratif = temp_stratif,
                                          errors = error_df, 
                                          printa = TRUE)
      
    #number of samples per strata
    allocations<-as.integer(temp_bethel)
      
    #cv
    cv_temp <- data.frame(PLANNED_CV=as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "]),
                          ACTUAL_CV=as.numeric(attributes(temp_bethel)$outcv[, "ACTUAL CV"]))
      
    #strata per cell
    temp_ids<-solution$indices
      
    #iterations for samples
    n_iter<-100
      
    #######################
    # allocations of samples per strata, year, sampling design
    # this piece of code can be run through the #11A script
    #######################
    
    # #to store samples
    # all_points<-array(NA,
    #                   dim = list(sum(allocations$n_samples),4,length(1982:2027),3),
    #                   dimnames = list(c(1:sum(allocations$n_samples)),c('Lon','Lat','cell','strata'),1:length(1982:2027),c('current','buffer','random')))
    # 
    # for (y in 1:length(1982:2027)) {
    #   
    #   #y<-1
    #   
    #   #print scenario to check progress
    #   cat(paste(" #############   year", y, 'of',length(1982:2027),  "  #############\n"))
    #   
    #   #to store points
    #   dfcurrent<-data.frame(matrix(NA,nrow=0,ncol=4))
    #   colnames(dfcurrent)<-c('Lon','Lat','cell','strata')
    #   dfbuffer<-dfrandom<-dfcurrent
    #   
    #   #for while purposes
    #   flag<-TRUE
    #   
    #   #random sample for each strata wit distance constrains
    #   for(n_istrata in 1:nrow(allocations)) {
    #     
    #     #n_istrata<-10
    #     
    #     istrata<-allocations[n_istrata,'Strata']
    #     
    #     #subset cells for strata
    #     df<-subset(as.data.frame(D8_2),strata==istrata)
    #     df1<-subset(as.data.frame(D8_2),strata==istrata & cell %in% baseline_strata$locations$cell)
    #     n_i<-allocations[n_istrata,'n_samples']
    #     
    #     #while error, keep running
    #     while(flag){
    #       
    #       #keep running if error
    #       tryCatch({
    #         
    #         #print loop state
    #         cat(paste(" -- strata", n_istrata, 'of',nrow(allocations),  " -- CURRENT  --\n"))
    #         
    #         ##############################
    #         # CURRENT STRATIFIED APPROACH
    #         ##############################
    #         
    #         if(samp_df[s,'samp_scn']=='scnbase'){
    #           pointsc<-baseline_strata$locations[which(baseline_strata$locations$stratum==istrata),c('longitude','latitude','cell','stratum')]
    #           names(pointsc)<-c('Lon' ,'Lat', 'cell','strata')
    #         }else if (samp_df[s,'samp_scn']=='scnbase_bis'){
    #           pointsc<-baseline_strata$locations[which(baseline_strata$locations$corner==FALSE & baseline_strata$locations$stratum==istrata),c('longitude','latitude','cell','stratum')]
    #           names(pointsc)<-c('Lon' ,'Lat',  'cell','strata')
    #         } else {
    #           
    #           #if more required samples than available
    #           if (nrow(df1)<n_i) {
    #             
    #             #vectors to store results
    #             dropcell<-c()
    #             selcell<-c()
    #             
    #             #duplicate df strata
    #             dff<-df
    #             
    #             #df removing available samples from current design
    #             dff<-subset(dff, !(cell %in% df1$cell))
    #             
    #             #cells to complete the required cells
    #             ii<-n_i-nrow(df1)
    #             
    #             #loop over the required samples using buffer
    #             for (iii in rep(1,times=ii)) {
    #               
    #               #ii<-1
    #               
    #               #get random cell
    #               cell_i<-sample(dff$cell,iii)
    #               #get row
    #               #row<-dff[which(dff$cell==cell_i),c('row','col')][1,1] 
    #               row<-dff[which(dff$cell %in% c(cell_i,df1$cell)),c('row','col')][1,1]
    #               #get col
    #               #col<-dff[which(dff$cell==cell_i),c('row','col')][1,2] 
    #               col<-dff[which(dff$cell %in% c(cell_i,df1$cell)),c('row','col')][1,2]
    #               #get adjacent cells
    #               adj_cells<-expand.grid(row=c(row,row+(1:100),row-(1:100)),
    #                                      col=c(col,col+(1:100),col-(1:100)))
    #               adj_cells1<-merge(adj_cells,df,by=c('row', 'col'))[,'cell']
    #               #remove cells from available for sampling
    #               dropcell_i<-c(cell_i,adj_cells1)
    #               
    #               #if no more samples available to sample stop
    #               if (nrow(subset(dff, !(cell %in% dropcell_i)))==0) {flag<-TRUE; stop("-- no samples to select on current approach",call. = FALSE)}else {
    #                 dff<-subset(dff, !(cell %in% dropcell_i))
    #                 selcell<-c(selcell,cell_i)
    #                 dropcell<-c(dropcell,dropcell_i)}
    #             }
    #             
    #             #subset points if equal to required number of samples, if not ERROR and rerun iteration
    #             if (nrow(subset(df,cell %in% c(selcell,df1$cell))) != n_i) {flag<-TRUE; stop("-- different number of points on current approach",call. = FALSE)}else{
    #               pointsc<-subset(df,cell %in% c(selcell,df1$cell))[,c('Lon','Lat','cell','strata')]
    #               names(pointsc)<-c('Lon','Lat','cell','strata')}
    #             
    #             #else there are enough samples to get from the current sampling design    
    #           } else {
    #             ss<-sample(1:nrow(df1),size = n_i,replace = FALSE)
    #             pointsc<-df1[ss,c('Lon','Lat','cell','strata')]}}
    #         
    #         ##############################
    #         # STRATIFIED + BUFFER APPROACH
    #         ##############################
    #         
    #         cat(paste(" -- strata", n_istrata, 'of',nrow(allocations),  " -- BUFFER  --\n"))
    #         
    #         #replicate df of strata     
    #         dff<-df
    #         
    #         #create vectors to store cells
    #         dropcell<-c()
    #         selcell<-c()
    #         
    #         #loop over required samples
    #         for (iii in rep(1,times=n_i)) {
    #           
    #           #iii<-1
    #           
    #           #random sample
    #           cell_i<-sample(dff$cell,iii)
    #           #get row of selected sample
    #           row<-dff[which(dff$cell==cell_i),c('row','col')][1,1]
    #           #get col of selected sample
    #           col<-dff[which(dff$cell==cell_i),c('row','col')][1,2]
    #           #get adjacent cells of selected sample
    #           adj_cells<-expand.grid(row=c(row,row+(1:100),row-(1:100)),
    #                                  col=c(col,col+(1:100),col-(1:100)))
    #           adj_cells1<-merge(adj_cells,df,by=c('row', 'col'))[,'cell']
    #           #drop cells is equal to selected cell and adjacent cells
    #           dropcell_i<-c(cell_i,adj_cells1)
    #           
    #           #store cells if there are still samples on the filtered df
    #           if (nrow(subset(dff, !(cell %in% dropcell_i)))==0) {flag<-TRUE;stop("-- no samples to select on buffer approach",call. = FALSE)} else {
    #             dff<-subset(dff, !(cell %in% dropcell_i))
    #             selcell<-c(selcell,cell_i)
    #             dropcell<-c(dropcell,dropcell_i)}
    #           
    #         }
    #         
    #         #subset points if equal to required number of samples, if not ERROR and rerun iteration
    #         if (nrow(subset(df,cell %in% selcell)) != n_i ) {flag<-TRUE;stop("-- different number of points on buffer approach",call. = FALSE)}else{
    #           pointsb<-subset(df,cell %in% selcell)[,c('Lon','Lat','cell','strata')]
    #           names(pointsb)<-c('Lon','Lat','cell','strata')}
    #         
    #         ##############################
    #         # STRATIFIED RANDOM APPROACH
    #         ##############################
    #         
    #         cat(paste(" -- strata", n_istrata, 'of',nrow(allocations),  " -- RANDOM  --\n"))
    #         
    #         #duplicate df strata
    #         dff<-df
    #         
    #         #random selection of samples
    #         selcell<-sample(dff$cell,n_i)
    #         pointsr<-subset(df,cell %in% selcell)[,c('Lon','Lat','cell','strata')]
    #         names(pointsr)<-c('Lon','Lat','cell','strata')
    #         
    #         #append data if buffer and current df have equal to required number of samples
    #         if (nrow(pointsc) == n_i & nrow(pointsb) == n_i & nrow(pointsr) == n_i) {
    #           dfcurrent<-rbind(dfcurrent,pointsc)
    #           dfbuffer<-rbind(dfbuffer,pointsb)
    #           dfrandom<-rbind(dfrandom,pointsr)
    #           flag<-FALSE}
    #       },
    #       
    #       #error message
    #       error=function(e) {
    #         message(paste0("ERROR RETRYING",':\n'),e)})
    #       if (!flag) next
    #     }
    #     
    #     #to start while again
    #     flag<-TRUE
    #   }
    #   
    #   #append points
    #   all_points[,,y,'current']<-unlist(dfcurrent)
    #   all_points[,,y,'buffer']<-unlist(dfbuffer)
    #   all_points[,,y,'random']<-unlist(dfrandom)
    # }
    # 
    # #save sample selection
    # save(all_points,file=paste0('./output/species/',sp,'/optimization data/samples_optimization_',samp_df[s,'samp_scn'],'_dynamic.RData'))
  
  
      #save list results por SBTscn and Sampscn
      result_list <- list(solution = solution,
                          sum_stats = sum_stats,
                          cvs = cv_temp,
                          sample_allocations = allocations,
                          sol_by_cell = temp_ids)
                          #str_cell = points1)
      
      #save plot list
      save(result_list,file=paste0('./output/species/',sp,'/optimization data/optimization_results_',samp_df[s,'samp_scn'],'.RData'))
      

   }
    
  #save locations
  save(sp_sum_stats,file = paste0('./output/species/',sp,'/optimization data/optimization_summary_stats.RData'))

}  
