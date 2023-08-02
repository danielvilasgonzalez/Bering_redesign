####################################################################
####################################################################
##    
##    simulate data and survey for historical and projected years
##    prepare estimates to compute design-based indices
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
project_yrs<-(last(yrs)+1):(last(yrs)+n_proj)

###################################
# # GRID NBS AND EBS
# ###################################
# 
# #load grid of NBS and EBS
# load('./extrapolation grids/northern_bering_sea_grid.rda')
# load('./extrapolation grids/eastern_bering_sea_grid.rda')
# grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
# grid$cell<-1:nrow(grid)
# 
# #load grid
# load('./extrapolation grids/lastversion_grid_EBS.RData')
# yrs<-1982:2022
# grid_ebs<-grid.ebs_year[which(grid.ebs_year$region != 'EBSslope' & grid.ebs_year$Year %in% yrs),]
# dim(grid_ebs)
# 
# #load baseline strata and specify corner stations
# load('./output/baseline_strata.RData')
# baseline_strata$locations[grep('GF|HG|JI|IH|ON|QP|PO',baseline_strata$locations$stationid),]
# baseline_strata$locations$corner<-ifelse(grepl('GF|HG|JI|IH|ON|QP|PO',baseline_strata$locations$stationid),'TRUE','FALSE')
# 
# ###################################
# # Sampling designs (from script #11) 
# ###################################
# 
# #sampling scenarios
# samp_df<-expand.grid(strat_var=c('Depth_varTemp','varTemp','Depth'),
#                      target_var=c('sumDensity'), #,'sqsumDensity'
#                      n_samples=c(520), #c(300,500) 520 (EBS+NBS+CRAB);26 (CRAB); 350 (EBS-CRAB); 494 (NBS-CRAB)
#                      n_strata=c(15),
#                      stringsAsFactors = FALSE) #c(5,10,15)
# 
# #add scenario number
# samp_df$samp_scn<-paste0(paste0('scn',1:nrow(samp_df)))
# samp_df<-rbind(samp_df,c('baseline','current',520,15,'scnbase'),
#                c('baseline w/o corner','current',494,15,'scnbase_bis'))
# #number of simulations

n_sim<- 100

#loop over spp
for (sp in spp) {
  
  #sp<-"Gadus macrocephalus"
  
  #get list of fit data
  ff<-list.files(paste0('./shelf EBS NBS VAST/',sp),'fit',recursive = TRUE)
  
  #load fit file
  load(paste0('./shelf EBS NBS VAST/',sp,'/',ff)) #fit
  #getLoadedDLLs() #if check loaded DLLs
  
  #reload model
  fit<-
    reload_model(x = fit)
  
  #check observation and predicted densities at each obs
  #observations
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
  # HISTORICAL DATA
  ######################
  
  #array to store simulated densities/CPUE
  sim_dens<-array(NA,
                  dim=c(nrow(grid),length(unique(yrs)),n_sim),
                  dimnames=c(1:nrow(grid),unique(yrs),1:n_sim))
  
  #create folder simulation data
  dir.create(paste0('./output/species/',sp,'/simulated historical data/'))
  
  #loop over simulations
  for (isim in 1:n_sim) { #simulations
    
    #isim<-1
    
    #print simulation to check progress
    cat(paste(" #############   Species", sp, match(sp,spp), 'out of',length(spp),  "  #############\n",
              " #############  historical simulation", isim, "of",n_sim, " #############\n"))

    #simulate data from OM
    Sim1 <- FishStatsUtils::simulate_data(fit = fit, #kg/km2
                                          type = 1, 
                                          random_seed = isim)
    
    #select simulated data that belong to grid points
    sim_bio <-matrix(data = Sim1$b_i[pred_TF == 1], #kg
                     nrow = nrow(grid), 
                     ncol = length(unique(yrs)))
    
    
    #biomass (kg) to CPUE (kg/km2)
    sim_dens[,,isim]<-sim_bio/grid$Area_in_survey_km2 
    
  }
  
  #save data
  save(sim_dens, file = paste0("./output/species/",sp,'/simulated historical data/sim_dens.RData')) 
  
  ######################
  # PROJECTED DATA
  ######################
  
  
}
