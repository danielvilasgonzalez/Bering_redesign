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

#call function distance points buffer
source('C:/Users/Daniel.Vilas/Work/GitHub/Bering_redesign/Scripts/genRandomPnts.R')

#setwd
out_dir<-'C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

#version VAST (cpp)
version<-'VAST_v13_1_0'

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

#sampling scenarios
samp_df<-expand.grid(strat_var=c('Lat_varTemp','Lat_meanTempF','Depth_meanTempF','Depth_varTemp','meanTempF_varTemp','meanTempF','varTemp','Depth'),
                    target_var=c('sumDensity'), #,'sqsumDensity'
                    n_samples=c(350), #c(300,500) 520 (EBS+NBS+CRAB);26 (CRAB); 350 (EBS-CRAB); 494 (NBS-CRAB)
                    n_strata=c(10)) #c(5,10,15)

#add scenario number
samp_df$samp_scn<-paste0(paste0('scn',1:nrow(samp_df)))

#########################
# RUN LOOP SPECIES
#########################

#loop over species
for (sp in spp) {
  
  sp<-'Gadus macrocephalus'
  
  #load optimization data
  load(paste0('./output/species/',sp,'/optimization_static_data.RData')) #D6
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
  
  #df locations
  sp_loc<-data.frame(matrix(nrow = 0,ncol=5))
  names(sp_loc)<- c("Lat","Lon","Stations","samp_scn","sp")
  
  #subset cells with appropiate depth
  static_df1<-subset(D6,cell %in% ok_cells)
    
  #create a list to store results
  plot_list<-list()
    
  #########################
  # RUN LOOP SAMPLING SCENARIOS
  #########################

  #loop through sampling scenarios
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
    # SIMPLE RANDOM SAMPLING CV CONSTRAINTS
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
    # STRATAS
    ###################################
      
    #get n_strata from kmean suggestion
    # kmean<-KmeansSolution2(frame=frame,
    #                        errors=cv,
    #                        maxclusters = 20,
    #                        showPlot = F)
    # 
    # #number strata from kmean
    # no_strata<-tapply(kmean$suggestions,
    #                   kmean$domainvalue,
    #                   FUN=function(x) length(unique(x)))
    # 
      
    #number of stratas 
    no_strata<-samp_df[s,'n_strata']
      
    ###################################
    # RUN OPTIMIZATION
    ###################################
      
    #run optimization
    solution <- optimStrata(method = "continuous", #continous variables
                            errors = cv,  #precision level - maximum allowable coefficient of variation set by the simple random sampling 
                            framesamp = frame, #df of input variables
                            iter = 50, #300 #aximum number of iterations
                            pops = 20, #100  #dimension of each generations
                            elitism_rate = 0.2, #0.1
                            mut_chance = 1 / (no_strata[1] + 1), #mutation chance
                            nStrata = no_strata,
                            showPlot = FALSE,
                            writeFiles = FALSE)
      
    ###################################
    # STORE SOLUTIONS
    ###################################
      
    #results from optimization
    framenew<-solution$framenew
    outstrata<-solution$aggr_strata
    ss<-summaryStrata(framenew,outstrata)
    #head(ss)
      
    #plot strata 2D
    #plotStrata2d(framenew,outstrata,domain=1,vars=c('X1','X2'),labels = c('VarTemp','Depth'))
      
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
                           sum_stats[,c("stratum_id","wh","Wh","M1","S1","SOLUZ","samp_scn","sp")])
    }
    
    #append stat results  
    sp_sum_stats<-rbind(sp_sum_stats,sum_stats)
      
    ###################################
    # CREATE SPATIAL OBJECT BASED ON CELLS STRATA
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
    #D8_2<-data.frame(D8_1)
      
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
    # MULTIVARIATE OPTIMAL ALLOCATION
    ###################################
      
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
    
    error_df <- data.frame("DOM" = "DOM1",
                           "CV1" = srs_cv,
                           "domainvalue"  = 1)
    
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
      
    #to store samples
    all_points<-array(NA,dim = list(sum(allocations),4,n_iter),
                        dimnames = list(c(1:sum(allocations)),c('Lon','Lat','cell','strata'),c(1:n_iter)))
      
    #loop over iterations
    for (iter in 1:n_iter) {
      
      #iter<-1
        
      #to store points
      dfpoints<-data.frame(matrix(NA,nrow=0,ncol=4))
      colnames(dfpoints)<-c('Lon','Lat','cell','strata')
        
      #random sample for each strata wit distance constrains
      for(istrata in 1:length(allocations)) {
          #istrata<-1
          
          #subset cells for strata
          df<-subset(as.data.frame(D8_1),Strata==istrata)
          
          #ratio of available cells and samples to take
          ratio<-dim(df)[1]/allocations[istrata]
          
          #select samples with buffer
          xy.buff<-buffer.f(data.frame(x=df$Lon,y=df$Lat,cell=df$cell),
                                       buffer = ratio*100, #30000
                                       reps = 1,
                                       n=allocations[istrata])
          colnames(xy.buff)[1:2]<-c("Lon",'Lat')
          xy.buff1<-data.frame(xy.buff,strata=istrata)
          
          dfpoints<-rbind(dfpoints,xy.buff1)
            
          #sample_vec<-c(sample_vec,xy.buff$cell)
          # sample_vec <- c(sample_vec,
          #                 sample(x = temp_ids[which(temp_ids$X1==istrata),'ID'], #which(temp_ids == istrata)
          #                        size = allocations[istrata]) )
        }
        
        #append points
        all_points[,,iter]<-unlist(dfpoints)
        
      }
      
      #save sample selection
      save(all_points,file=paste0('./output/species/',sp,'/samples_optimization_',samp_df[s,'samp_scn'],'.RData'))

      #change projection of spatial object 
      coordinates(dfpoints)<- ~ Lon + Lat
      
      #reproject shapefile
      proj4string(dfpoints) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
      #points1<-spTransform(dfpoints,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))
      
      #to dataframe
      points1<-as.data.frame(dfpoints)
      
      #add scn 
      points1$samp_scn<- s
      points1$sp<- sp

      #save list results por SBTscn and Sampscn
      result_list <- list(solution = solution,
                          sum_stats = sum_stats,
                          cvs = cv_temp,
                          sample_allocations = allocations,
                          sol_by_cell = temp_ids,
                          str_cell = points1)
      
      #save plot list
      save(result_list,file=paste0('./output/species/',sp,'/optimization_results_',samp_df[s,'samp_scn'],'.RData'))
      
      #########################
      # JOIN POINTS FOR LEGEND PURPOSES
      #########################
      
      loc_sur<-rbind(data.frame(Lat=points1$Lat,Lon=points1$Lon,Stations='optimization'),
                      data.frame(Lat=st_EBS$latitude,Lon=st_EBS$longitude,Stations='current design'),
                      data.frame(Lat=st_corners1$latitude,Lon=st_corners1$longitude,Stations='corner crab'))
      
      loc_sur$samp_scn<-samp_df[s,'samp_scn']
      loc_sur$sp<-sp
      
      loc_sur1<-subset(loc_sur,Stations=='optimization')
      
      sp_loc<-rbind(sp_loc,loc_sur1)
      
      #########################
      # MAP POINTS
      #########################
    
      p<-
        ggplot()+
            geom_raster(data=as.data.frame(r2, xy = TRUE),aes(x=x,y=y,fill=layer))+
            #geom_point(data=D8_2, aes(Lon, Lat, fill=Strata, group=NULL),size=2, stroke=0,shape=21)+
            scale_fill_gradientn(colours=c("#ea5545", "#f46a9b", "#ef9b20", "#edbf33", "#ede15b", "#bdcf32", "#87bc45", "#27aeef", "#b33dc6"),
                                 #guide = guide_legend(),na.value = 'white',breaks=sort(unique(D8_1$Strata)),
                                 #labels=paste0(sort(unique(D8_1$Strata))," (n=",allocations,')'))+ #,,
                                  guide = guide_legend(),na.value = 'white',breaks=sort(na.omit(unique(values(r2)))),
                                  labels=paste0(sort(na.omit(unique(values(r2))))," (n=",allocations,')'))+ #,,
            #geom_point(data=st_EBS,aes(x=longitude,y=latitude),shape=4,size=1)+
            #geom_point(data=st_corners1,aes(x=longitude,y=latitude),color='red',shape=20,size=1)+
            geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
            #geom_polygon(data=eez_sh33,aes(x=long,y=lat,group=group),fill=NA,color='grey40',linetype='dashed')+
            scale_x_continuous(expand = c(0,0),position = 'bottom',
                               breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
            geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
            geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
            geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
            geom_point(data=loc_sur,aes(x=Lon,y=Lat,color=Stations,shape=Stations),fill='white',color='black',size=1.2)+
            scale_shape_manual(values = c('optimization'=21,
                                          'current design'=4,
                                          'corner crab'=8),
                               breaks=unique(loc_sur$Stations),
                               labels=paste0(unique(loc_sur$Stations)," (n=",c(nrow(points1),nrow(st_EBS),nrow(st_corners1)),')'))+
            coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
                     xlim = panel_extent$x,
                     ylim = c(453099.5,2004909.7),
                         label_axes = "-NE-")+
            theme(aspect.ratio = 1,panel.grid.major = element_line(color = rgb(0, 0, 0,20, maxColorValue = 285), linetype = 'dashed', linewidth =  0.5),
                  panel.background = element_rect(fill = NA),panel.ontop = TRUE,text = element_text(size=12),
                  legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(14, 'points'),
                  legend.key.width= unit(12, 'points'),axis.title = element_blank(),
                  legend.position = c(0, 1), legend.justification = c(0, 1),
                  panel.border = element_rect(fill = NA, colour = 'black'),legend.key = element_rect(color="black"),
                  axis.text = element_text(color='black'),legend.spacing.y = unit(8, 'points'),
                  axis.text.y.right = element_text(hjust= 0.1 ,margin = margin(0,0,0,-28, unit = 'points'),color='black'),
                  axis.text.x.bottom = element_text(vjust = 6, margin = margin(-5,0,0,0, unit = 'points'),color='black'),
                  axis.ticks.length = unit(-5,"points"),plot.title = element_text(size=11,vjust = -10, hjust=0.95,face="bold"),
                  plot.margin=margin(c(-10,0,0,-10)),legend.title = element_text(size = 12))+
                #annotate("text", x = -256559, y = 1354909, label = "Alaska",parse=TRUE,size=7)+
                #annotate("text", x = -1296559, y = 2049090, label = "Russia",parse=TRUE,size=7)+
                scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())+
                #annotate("text", x = -1296559, y = 744900, label = "italic('Bering Sea')",parse=TRUE,size=9)+
                guides(fill = guide_legend(order=2,override.aes=list(size=4),title="Strata"),
                       color = guide_legend(order=1,override.aes=list(size=8)),
                       shape = guide_legend(order=1),override.aes=list(size=8))+
                labs(title=paste0(samp_df[s,'strat_var'],' n=',samp_df[s,'n_samples']))

        #store into the list
        plot_list[[s]]<-p  
   }
    
    #save multiplot
    mp<-cowplot::plot_grid(plotlist = plot_list,nrow = 3,ncol = 3)
    ragg::agg_png(paste0('./figures/species/',sp,'/optimization_sampling.png'), width = 14, height = 14, units = "in", res = 300)
    print(mp)
    dev.off()

    #save plot list
    save(plot_list,file=paste0('./output/species/',sp,'/optimization_plots.RData'))
  
  
  #save locations
  save(sp_loc,file = paste0('./output/species/',sp,'/optimization_locations.RData'))
  
  #save locations
  save(sp_sum_stats,file = paste0('./output/species/',sp,'/optimization_summary_stats.RData'))

}  
