####################################################################
####################################################################
##    
##    simulate datafrom fit object, simulate data and get design sample index
##    get bias
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

#install coldpool to extract SBT for the EBS
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

###################################
# GRID NBS AND EBS
###################################

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)


###################################
# SAMPLING SCENARIOS from script #11 
###################################

#sampling scenarios
samp_df<-expand.grid(strat_var=c('Lat_varTemp','Lat_meanTempF','Depth_meanTempF','Depth_varTemp','meanTempF_varTemp','meanTempF','varTemp','Depth'),
                     target_var=c('sumDensity'), #,'sqsumDensity'
                     n_samples=c(350), #c(300,500) 520 (EBS+NBS+CRAB);26 (CRAB); 350 (EBS-CRAB); 494 (NBS-CRAB)
                     n_strata=c(10)) #c(5,10,15)

#add scenario number
samp_df$samp_scn<-paste0(paste0('scn',1:nrow(samp_df)))


#loop over spp
for (sp in spp) {
  
  sp<-"Gadus macrocephalus"
  
  #load fit file
  #getLoadedDLLs()
  load(paste0('./shelf EBS NBS VAST/',sp,'/fit.RData')) #fit

  Sys.time()
  #reload model
  fit<-
    reload_model(x = fit)
  Sys.time()
  
  #check observation and predicted densities at each obs
  data_geostat<-readRDS(paste0('./shelf EBS NBS VAST/',sp,'/data_geostat_temp.rds'))
  d_i<-fit$Report$D_i #nrow(fit$data_frame)
  length(d_i)==nrow(data_geostat)
  
    
  # fit_sim <-
  #  
  #     fit_model(settings=fit$settings,
  #               Lat_i=data_geostat$Lat, 
  #               Lon_i=data_geostat$Lon,
  #               t_i=data_geostat$Year,
  #               b_i=data_geostat$CPUE_kg,
  #               c_iz = as.numeric(factor(data_geostat$Species))-1,
  #               a_i=data_geostat$Effort/100,
  #               #input_grid=grid.ebs,
  #               getJointPrecision = TRUE,
  #               test_fit=FALSE,
  #               create_strata_per_region = TRUE,  
  #               covariate_data = fit$covariate_data, 
  #               X1_formula =  fit$X1_formula,
  #               X2_formula = fit$X2_formula, 
  #               newtonsteps = 1,
  #               Parameters=fit$ParHat,
  #               #X_gtp = X_gtp,
  #               working_dir = './shelf EBS NBS VAST/Gadus macrocephalus/test/')
  
           
  ##################################################
  ####   Simulate 1000 iterations of data
  ##################################################
  
  sim_data <- array(data = NA, dim = c(fit$spatial_list$n_g,
                                       length(fit$year_labels),
                                       1000))
  
  sim_data <- array(data = NA, dim = c(fit$spatial_list$n_g,
                                       length(fit$year_labels),
                                       1000))
  
  for (isim in 1:1000) {
    
    isim<-10
    
    Sim1 <- FishStatsUtils::simulate_data(fit = fit, 
                                          type = 1, 
                                          random_seed = isim)
    Sim1< - fit$D
    
    #PredTF_i = 1 just excludes a given datum from being used in the joint likelihoodPredTF_i	
    #OPTIONAL, whether each observation i is included in the likelihood, PredTF_i[i]=0, or in the predictive probability, PredTF_i[i]=1
    sim_data[, , isim] <- data.frame(Sim1$D_gct[,1,])
    sim_surv_data[,,isim,iter]
    
    
    for (samp in unique(samp_df$samp_scn)) {
      
      samp<-unique(samp_df$samp_scn)[1]
      
    
    
    #get station locations for each sampling design
    load(file=paste0('./output/species/',sp,'/samples_optimization_',samp,'.RData')) #all_points
    
    for (i in dimnames(all_points)[[3]]) {
      
      i<-dimnames(all_points)[[3]][1]
                                   
         points<-data.frame(all_points[,,i])                          
               
         coordinate                     
      dimnames(all_points)
      
    }
    
    #if(isim%%100 == 0) print(paste("Done with", ispp, "Iteration", isim))
    }
  }
  