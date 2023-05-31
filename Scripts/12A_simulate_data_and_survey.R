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

#load grid
load('./extrapolation grids/lastversion_grid_EBS.RData')
yrs<-1982:2022
grid_ebs<-grid.ebs_year[which(grid.ebs_year$region != 'EBSslope' & grid.ebs_year$Year %in% yrs),]
dim(grid_ebs)

#load baseline strata
load('./output/baseline_strata.RData')
baseline_strata$locations[grep('GF|HG|JI|IH|ON|QP|PO',baseline_strata$locations$stationid),]
baseline_strata$locations$corner<-ifelse(grepl('GF|HG|JI|IH|ON|QP|PO',baseline_strata$locations$stationid),'TRUE','FALSE')


###################################
# SAMPLING SCENARIOS from script #11 
###################################

#sampling scenarios
samp_df<-expand.grid(strat_var=c('Depth_varTemp','varTemp','Depth'),
                     target_var=c('sumDensity'), #,'sqsumDensity'
                     n_samples=c(520), #c(300,500) 520 (EBS+NBS+CRAB);26 (CRAB); 350 (EBS-CRAB); 494 (NBS-CRAB)
                     n_strata=c(15)) #c(5,10,15)

#add scenario number
samp_df$samp_scn<-paste0(paste0('scn',1:nrow(samp_df)))

#number of simulations
n_sim<- 100

#loop over spp
for (sp in spp) {
  
  sp<-"Gadus macrocephalus"
  
  #load fit file
  #getLoadedDLLs()
  load(paste0('./shelf EBS NBS VAST/',sp,'/fit.RData')) #fit

  # Sys.time()
  # #reload model
   fit<-
     reload_model(x = fit)
  # Sys.time()
  
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
  
  #################
  #get predTF
  #################
  
  #read data_geostat_temp file
  df1<-readRDS(paste0('./data processed/species/',sp,'/data_geostat_temp.rds'))
  #df1[which(df1$year==2020),'bottom_temp_c']<-NA
  df2<-subset(df1,year %in% 1982:2022)
  
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
           
  ##################################################
  ####   Simulate 100 iterations of historical data
  ##################################################
  
  #array to store results
  # sim_data <- array(data = NA, dim = c(fit$spatial_list$n_g,
  #                                      length(fit$year_labels),
  #                                      n_sim))
  # 
  
  #create folder simulation data
  dir.create(paste0('./output/species/',sp,'/simulated historical data/'))
  
  #create folder simulation number
  dir.create(paste0('./output/species/',sp,'/survey simulation historical data/'))
  
  # sim_data <- array(data = NA, dim = c(fit$spatial_list$n_g,
  #                                      length(fit$year_labels),
  #                                      1000))
  # 
  for (isim in 1:n_sim) { #simulations
    
    #isim<-1
    
    #print scenario to check progress
    cat(paste(" #############   Species", sp, match(sp,spp), 'out of',length(spp),  "  #############\n",
              " #############  historical simulation", isim, "of",n_sim, " #############\n"))
    
    #simulation folder
    sim_fold<-formatC(isim, width = 4, format = "d", flag = "0")
    
    #create folder simulation number
    dir.create(paste0('./output/species/',sp,'/survey simulation historical data/',sim_fold))
    
    #create folder simulation number
    dir.create(paste0('./output/species/',sp,'/simulated historical data/',sim_fold))
    
    #simulate data from OM
    Sim1 <- FishStatsUtils::simulate_data(fit = fit, 
                                          type = 1, 
                                          random_seed = isim)
    
    sim_dens <-matrix(data = Sim1$b_i[pred_TF == 1], 
                      nrow = nrow(grid), 
                      ncol = length(unique(yrs)))
    
    
    #PredTF_i = 1 just excludes a given datum from being used in the joint likelihoodPredTF_i	
    #OPTIONAL, whether each observation i is included in the likelihood, PredTF_i[i]=0, or in the predictive probability, PredTF_i[i]=1
    
    #append simulated densities
    #sim_dens <- unlist(Sim1$D_gct[,1,])
    #append simulated index   
    #sim_ind<- unlist(Sim1$Index_ctl)
    
    #list of densities and index
    #sim_data<-list(sim_dens=sim_dens,sim_ind=sim_ind)
    
    #save data
    save(sim_dens, file = paste0("./output/species/",sp,'/simulated historical data/',sim_fold,'/sim_dens.RData'))  
    
    for (samp in c('scnbase','scnbase_bis',unique(samp_df$samp_scn))) { #sampling designs
      
      #samp<-"scn1"
      
      if (!(samp %in% unique(samp_df$samp_scn))) {
        
        #conditions on baseline scenarios
        if (samp == 'scnbase') {
          points<-baseline_strata$locations
        } else if (samp == 'scnbase_bis') {
          points<-baseline_strata$locations[which(baseline_strata$locations$corner=='FALSE'),]
        } 
        
        #to store results
        sim_survey <- array(data = NA, dim = c(nrow(points),
                                               length(fit$year_labels)+1, #add strata number
                                               1),
                            dimnames = list(c(1:nrow(points)),
                                            c(fit$year_labels,'strata'),
                                            1))
        
        #store results
        sim_survey[,,1]<-cbind(sim_dens[points$cell,],points$stratum)
        
        #store results
        save(sim_survey, file = paste0("./output/species/",sp,'/survey simulation historical data/',sim_fold,'/sim_survey_',samp,'.RData'))  
        
      } else {
      
      #get station locations for each sampling design
      load(file=paste0('./output/species/',sp,'/optimization data/samples_optimization_',samp,'.RData')) #all_points
      
      #to store results
      sim_survey <- array(data = NA, dim = c(length(dimnames(all_points)[[1]]),
                                             length(fit$year_labels)+1, #add strata number
                                             length(dimnames(all_points)[[3]]),
                                             2),
                          dimnames = list(dimnames(all_points)[[1]],
                                          c(fit$year_labels,'strata'),
                                          dimnames(all_points)[[3]],
                                          c('current','buffer')))
      
      for (i in as.integer(dimnames(all_points)[[3]])) { #iters
        
        #i<-dimnames(all_points)[[3]][1]
        
        #get points iterations
        pointsc<-data.frame(unlist(all_points[,,i,'current']))
        pointsb<-data.frame(unlist(all_points[,,i,'buffer']))                          
        
        
        #append survey densities for each iteration and simulated data
        sim_survey[,,i,'current']<-pointsc
        sim_survey[,,i,'buffer']<-pointsb
        
      }
      
      #store results
      save(sim_survey, file = paste0("./output/species/",sp,'/survey simulation historical data/',sim_fold,'/sim_survey_',samp,'.RData'))  
   }
  }
 }
}

rm(fit)
rm(Sim1)

##################################################
####   Simulate 100 iterations of projected data
##################################################

#loop over spp
for (sp in spp) {
  
  sp<-"Gadus macrocephalus"
  
  #create folder simulation number
  dir.create(paste0('./output/species/',sp,'/survey simulation projected data/'))
 
  
  
  for (sbt in paste0('SBT',1:12)) {
      
    #sbt<-'SBT1'
    
    #print scenario to check progress
    cat(paste(" #############   Species", sp, match(sp,spp), 'out of',length(spp),  "  #############\n",
              " #############  projected simulation", sbt, "of", 'SBT12' , " #############\n"))
  
  
      #load sbt data
      load(paste0('./output/species/',sp,'/simulated projected data/fit_projection_',sbt,'.RData'))
      
      
      for (isim in 1:n_sim) { #simulations
        
        #isim<-1
        
        #simulation folder
        sim_fold<-formatC(isim, width = 4, format = "d", flag = "0")
        
        #create folder simulation number
        dir.create(paste0('./output/species/',sp,'/simulated projected data/',sim_fold),showWarnings = FALSE)
        
        #create folder simulation number
        dir.create(paste0('./output/species/',sp,'/survey simulation projected data/',sim_fold),showWarnings = FALSE)
        
        #simulated data
        sim_dens<-drop_units(unlist(pm[[isim]]$D_gct[,1,42:46]))
        sim_ind<-drop_units(unlist(pm[[isim]]$Index_ctl[,42:46,]))
        
        #list of results of sbt
        sim_data<-list(sim_dens=sim_dens,sim_ind=sim_ind)
        
        #save data
        save(sim_data, file = paste0("./output/species/",sp,'/simulated projected data/',sim_fold,'/sim_data_',sbt,'.RData'))  
        
        for (samp in c('scnbase','scnbase_bis',unique(samp_df$samp_scn))) { #sampling designs
          
          #samp<-"scn1"
          
          if (!(samp %in% unique(samp_df$samp_scn))) {
            
           
            #conditions on baseline scenarios
            if (samp == 'scnbase') {
              points<-baseline_strata$locations
            } else if (samp == 'scnbase_bis') {
              points<-baseline_strata$locations[which(baseline_strata$locations$corner=='FALSE'),]
            } 
            
            #to store results
            sim_survey <- array(data = NA, dim = c(nrow(points),
                                                   5+1, #add strata number
                                                   1),
                                dimnames = list(c(1:nrow(points)),
                                                c(2023:2027,'strata'),
                                                1))
            
            #store results
            sim_survey[,,1]<-cbind(sim_data$sim_dens[points$cell,],points$stratum)
            
            save(sim_survey, file = paste0("./output/species/",sp,'/survey simulation projected data/',sim_fold,'/sim_survey_',samp,'_',sbt,'.RData'))  
            
          } else {
            
            #get station locations for each sampling design
            load(file=paste0('./output/species/',sp,'/optimization data/samples_optimization_',samp,'.RData')) #all_points
            
            #to store results
            sim_survey <- array(data = NA, dim = c(length(dimnames(all_points)[[1]]),
                                                   5+1, #add strata number
                                                   length(dimnames(all_points)[[3]]),
                                                   2),
                                dimnames = list(dimnames(all_points)[[1]],
                                                c(2023:2027,'strata'),
                                                dimnames(all_points)[[3]],
                                                c('current','buffer')))
            
            for (i in as.integer(dimnames(all_points)[[3]])) { #iters
              
              #get points iterations
              pointsc<-data.frame(unlist(all_points[,,i,'current']))
              pointsb<-data.frame(unlist(all_points[,,i,'buffer']))                          
              
              #append survey densities for each iteration and simulated data
              sim_survey[,,i,'current']<-pointsc
              sim_survey[,,i,'buffer']<-pointsb
            }
            
            #store results
            save(sim_survey, file = paste0("./output/species/",sp,'/survey simulation projected data/',sim_fold,'/sim_survey_',samp,'_',sbt,'.RData'))  
          }
        }
      }
      rm(pm)
  }
}
    

  