library(VAST)

#check errors in models

#ssp with problems
spp<-c('Reinhardtius hippoglossoides', #missbehave model - year 2019 peack of biomass
       'Anoplopoma fimbria', #may be the reality for sablefish
       'Lepidopsetta polyxystra') #missbehave model - from 1982-1995 very low biomass (almost 0)

sp<-spp[1]
sp<-'Gadus macrocephalus'
load(paste0('./shelf EBS NBS VAST/',sp,'/fit-001.RData'))
load(paste0('./shelf EBS NBS VAST/',sp,'/fit.RData'))
plot(fit$data_frame$t_i,fit$data_frame$b_i)

y<-readRDS(paste0('./shelf EBS NBS VAST/',sp,'/data_geostat_temp.RDS'))
plot(y$Year,y$CPUE_kg)
subset(y,Year==2020)

#save data_geostat with SBT
yy<-readRDS(paste0('./data processed/species/',sp,'/data_geostat_temp.rds'))
subset(yy,year==2020)


z<-readRDS(paste0('./data processed/species/',sp,'/data_geostat.rds'))
subset(z,year==2020)
plot(z$year,z$cpue_kgha)
tapply(z$bottom_temp_c, z$year, summary)
zz<-subset(z,year==2019)


hist(zz$bottom_temp_c,breaks=20)

#plot(fit,working_dir ='./test/')
plot(x=1982:2022,y=drop_units(fit$Report$Index_ctl[,,1]))

plot(fit$data_frame$t_i,fit$data_frame$b_i)

xx<-subset(fit$data_frame,t_i>=2010)


x<-readRDS('./shelf EBS NBS VAST/',sp,'/data_geostat_temp.rds')
summary(x)
xx<-x[which(x$Year==2019),]
summary(xx$CPUE_kg)


######################################################
######################################################

#VARIANCE AMONG PROJECTIONS 

# Calculate Mean: For each time point, calculate the mean (average) value across all the projections. This will give you a baseline reference point.
# Calculate Variance: For each time point, calculate the squared difference between each projection's value and the mean value calculated in step 3. Sum up these squared differences for all projections and divide by the number of projections minus 1 (to get an unbiased estimate) to calculate the variance at that time point.
# Variance = Σ(projection_value - mean_value)^2 / (number_of_projections - 1)
# Repeat for All Time Points: Repeat step 4 for every time point in your time series to obtain a variance value for each time point.
# Interpret Results: The calculated variances at different time points represent the variability or spread of projections around the mean. Larger variance values indicate more disagreement or uncertainty among the projections, while smaller variance values suggest greater agreement.



####################################################################
####################################################################
##    
##    simulate data from OM for historical and projected years
##    danielvilasgonzalez@gmail.com/dvilasg@uw.edu
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 

#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c('ggplot2','units','splines','raster','sp')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#install VAST if it is not
if (!('VAST' %in% installed.packages())) {
  devtools::install_github("james-thorson/VAST@main", INSTALL_opts="--no-staged-install")};library(VAST)

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#setwd
out_dir<-'C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/'  
out_dir<-'/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

#list of sp
spp<-list.dirs('./data processed/species/',full.names = FALSE,recursive = FALSE)

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

#yrs
yrs<-1982:2022

#how manyt projected years we want
n_proj<-5

#project_yrs
project_yrs<-((yrs[length(yrs)])+1):(yrs[length(yrs)]+n_proj)

###################################
# GRID NBS AND EBS
###################################

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)

#load grid
load('./extrapolation grids/lastversion_grid_EBS.RData')
yrs<-1982:2022
grid_ebs<-grid.ebs_year[which(grid.ebs_year$region != 'EBSslope' & grid.ebs_year$Year %in% yrs),]
dim(grid_ebs)

#load baseline strata and specify corner stations
load('./output/baseline_strata.RData')
baseline_strata$locations[grep('GF|HG|JI|IH|ON|QP|PO',baseline_strata$locations$stationid),]
baseline_strata$locations$corner<-ifelse(grepl('GF|HG|JI|IH|ON|QP|PO',baseline_strata$locations$stationid),'TRUE','FALSE')

###################################
# Sampling designs (from script #11) 
###################################

#sampling scenarios
samp_df<-expand.grid(strat_var=c('Depth_varTemp','varTemp','Depth'),
                     target_var=c('sumDensity'), #,'sqsumDensity'
                     n_samples=c(520), #c(300,500) 520 (EBS+NBS+CRAB);26 (CRAB); 350 (EBS-CRAB); 494 (NBS-CRAB)
                     n_strata=c(15),
                     stringsAsFactors = FALSE) #c(5,10,15)

#add scenario number
samp_df$samp_scn<-paste0(paste0('scn',1:nrow(samp_df)))
samp_df<-rbind(samp_df,c('baseline','current',520,15,'scnbase'),
               c('baseline w/o corner','current',494,15,'scnbase_bis'))


###################################
# SBT projections
###################################

#save SBT table
load('./tables/SBT_projection.RData')#df_sbt

#name scenario
df_sbt$sbt<-paste0('SBT',df_sbt$sbt_n)
df_sbt$sbt2<-paste0(df_sbt$sbt,'_',df_sbt$Scenario)

df_sbt1<-subset(df_sbt,sbt %in% paste0('SBT',c(1,2,3,5,7))) #status quo, warm moderate variation, warm very high variation, gradually warm, severe gradually warm

#number of simulations
n_sim<- 30

#store index and dens
dens_index_proj_OM<-list()

#array to store simulated densities/CPUE
# sim_proj_dens_spp<-array(NA,
#                          dim=c(nrow(grid),length(project_yrs),1,length(1:8),length(spp)),
#                          dimnames=list(1:nrow(grid),project_yrs,1,1:8,spp))


sp<-"Gadus macrocephalus"
  
  #create folder simulation data
  dir.create(paste0('./output/species/',sp,'/'))
  
  #get list of fit data
  ff<-list.files(paste0('./shelf EBS NBS VAST/',sp),'fit',recursive = TRUE)
  
  #load fit file
  load(paste0('./shelf EBS NBS VAST/',sp,'/',ff)) #fit
  #getLoadedDLLs() #if check loaded DLLs
  
  #reload model
  #fit<-
  #  reload_model(x = fit)
  data_geostat<-readRDS(paste0('./shelf EBS NBS VAST/',sp,'/data_geostat_temp.rds'))
  
  #predictions
  #d_i<-fit$Report$D_i #nrow(fit$data_frame)
  #length(d_i)==nrow(data_geostat)
  
  #################
  # get predTF (required argument to get predictions on grid when simulating data)
  #################
  
  #read data_geostat_temp file
  df1<-readRDS(paste0('./data processed/species/',sp,'/data_geostat_temp.rds'))
  df2<-subset(df1,year %in% yrs)
  
  #select rows and rename
  df3<-df2[,c("lat_start","lon_start","year",'scientific_name','weight_kg','effort','depth_m','LogDepth',"ScaleLogDepth",'Scalebottom_temp_c','bottom_temp_c','survey_name')]
  colnames(df3)<-c('Lat','Lon','Year','Species','CPUE_kg','Effort','Depth','LogDepth','ScaleLogDepth','ScaleBotTemp','BotTemp','Region')
  
  #data geostat
  df4<-subset(df3,Region %in% c("Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey",
                                "Northern Bering Sea Crab/Groundfish Survey - Eastern Bering Sea Shelf Survey Extension"))
  
  data_geostat<-df4[complete.cases(df4[,c('CPUE_kg')]),]
  
  #covariate data - filter by year and complete cases for env variables
  #covariate_data<-subset(df2,Year>=yrs_region[1] & Year<=yrs_region[2])
  covariate_data<-df3[complete.cases(df3[,c('BotTemp')]),] #,'ScaleLogDepth'
  
  #add grid to get prediction for simulate data on each cell of the grid (sim$b_i)
  grid_df<-data.frame(Lat=grid_ebs$Lat,
                      Lon=grid_ebs$Lon,
                      Year=grid_ebs$Year,
                      Species=rep(sp,times=nrow(grid_ebs)),
                      CPUE_kg=mean(data_geostat$CPUE_kg),
                      Effort=grid_ebs$Area_in_survey_km2,
                      Depth=grid_ebs$Depth,
                      BotTemp=grid_ebs$Temp,
                      Region=grid_ebs$region,
                      stringsAsFactors = T)
  
  #ha to km2
  data_geostat$Effort<-data_geostat$Effort/100
  
  #rbind grid and data_geostat to get prediction into grid values when simulating data
  data_geostat1<-rbind(data_geostat[,c("Lat","Lon","Year","Species","CPUE_kg","Effort","Depth","BotTemp","Region")],
                       grid_df)
  
  #to get predictions in locations but not influencing fit
  pred_TF <- rep(1, nrow(data_geostat1))
  pred_TF[1:nrow(data_geostat)] <- 0
  
  ######################
  # PROJECTED DATA
  ######################
  
  #get raster stack
  stack_files<-list.files('./data processed/SBT projections/')
  
  #create folder simulation data
  dir.create(paste0('./output/species/',sp,'/simulated projected data/'))
  
  #loop over scenarios
  for (sbt in unique(df_sbt1$sbt_n)) {
    
    #sbt<-unique(df_sbt$sbt_n)[1]
    
    #print scenario to check progress
    cat(paste(" #############     PROJECTING    #############\n",
              " #############   Species", sp, match(sp,spp), 'out of',length(spp),  "  #############\n",
              " #############  SBT", sbt, " #############\n"))
    
    #open stack of rasters
    st<-stack_files[grepl(paste0('SBT_',sbt),stack_files)][1]
    st<-stack(paste0('./data processed/SBT projections/',st))
    #plot(st)
    #title(main='x')
    
    #reproject shapefile
    #proj4string(st) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs') 
    
    #raster to points
    points<-data.frame(rasterToPoints(st))
    
    #load fit file
    #load(paste0('./shelf EBS NBS VAST/',sp,'/',ff))
    
    #create a df to store
    points3<-data.frame(matrix(nrow = 0,ncol = ncol(fit$covariate_data)))
    names(points3)<-names(fit$covariate_data)
    
    for (y in project_yrs) {
      
      #y<-project_yrs[1]
      
      #get points for year
      points1<-points[,c('x','y',paste0('y',y))]
      names(points1)<-c('Lon',"Lat",'BotTemp')
      
      #reproject df
      coordinates(points1)<- ~ Lon + Lat
      proj4string(points1) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs') 
      points1<-data.frame(spTransform(points1,CRSobj = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")))
      names(points1)<-c('BotTemp','Lon',"Lat")
      
      #create a new df
      points2<-data.frame(cbind(Year=y,Lat=points1$Lat,Lon=points1$Lon,ScaleLogDepth=NA,LogDepth=NA,ScaleBotTemp=NA,BotTemp=points1$BotTemp,CPUE_kg=NA))
      
      #add year
      points3<-rbind(points3,points2)
      
    }
    
    #add year to covariate data from fit
    #cov_list[[sbt]]<-points3
    
    #add covariate data
    new_data<-rbind(fit$covariate_data,points3)
    
    #project model example
    pm<-VAST::project_model(x = fit,
                            working_dir = paste0('./shelf EBS NBS VAST/',sp,'/'),
                            n_proj = n_proj,
                            n_samples = n_sim, #n_sim?
                            new_covariate_data = new_data,
                            historical_uncertainty = 'none')
    
    #remove fit
    #rm(fit)
    #dyn.unload('C:/Program Files/R/R-4.2.2/library/VAST/executables/VAST_v13_1_0_TMBad.dll')
    save(pm, file = paste0("./output/species/",sp,'/simulated projected data/test_fit_projection_SBT',sbt,'.RData'))
    
    #store
    #sim_proj_dens_spp[,,1,sbt,sp]<-pm$D_gct[,1,42:46]
    
    #store index and dens
    #index<-pm$Index_ctl
    #dens<-pm$D_gct[,1,]
    #dens_index_proj_OM[[sbt]]<-pm

    rm(pm)
    gc()
  }


  #######################
  # VARIANCE PROJECTIONS - COMPARE/RANK AMONG SBT SCENARIOS
  #######################
  
  
  #VARIANCE AMONG PROJECTIONS 
  
  # Calculate Mean: For each time point, calculate the mean (average) value across all the projections. This will give you a baseline reference point.
  # Calculate Variance: For each time point, calculate the squared difference between each projection's value and the mean value calculated in step 3. Sum up these squared differences for all projections and divide by the number of projections minus 1 (to get an unbiased estimate) to calculate the variance at that time point.
  # Variance = Σ(projection_value - mean_value)^2 / (number_of_projections - 1)
  # Repeat for All Time Points: Repeat step 4 for every time point in your time series to obtain a variance value for each time point.
  # Interpret Results: The calculated variances at different time points represent the variability or spread of projections around the mean. Larger variance values indicate more disagreement or uncertainty among the projections, while smaller variance values suggest greater agreement.
  
  
  