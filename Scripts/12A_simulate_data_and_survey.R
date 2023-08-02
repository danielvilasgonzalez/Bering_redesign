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
#number of simulations
n_sim<- 100

#loop over spp
for (sp in spp) {
  
  sp<-"Gadus macrocephalus"
  
  #load fit file
  load(paste0('./shelf EBS NBS VAST/',sp,'/fit.RData')) #fit
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
           
  ######################
  # HISTORICAL DATA
  ######################
  
  #create folder simulation data
  dir.create(paste0('./output/species/',sp,'/simulated historical data/'))
  
  #create folder simulation number
  dir.create(paste0('./output/species/',sp,'/survey simulation historical data/'))
  
  #loop over simulations
  for (isim in 1:n_sim) { #simulations
    
    #isim<-1
    
    #print simulation to check progress
    cat(paste(" #############   Species", sp, match(sp,spp), 'out of',length(spp),  "  #############\n",
              " #############  historical simulation", isim, "of",n_sim, " #############\n"))
    
    #simulation folder
    sim_fold<-formatC(isim, width = 4, format = "d", flag = "0")
    
    #create folder simulation number
    dir.create(paste0('./output/species/',sp,'/survey simulation historical data/',sim_fold))
    
    #create folder simulation number
    dir.create(paste0('./output/species/',sp,'/simulated historical data/',sim_fold))
    
    #simulate data from OM
    Sim1 <- FishStatsUtils::simulate_data(fit = fit, #kg/km2
                                           type = 1, 
                                           random_seed = isim)
     
    #select simulated data that belong to grid points
    sim_bio <-matrix(data = Sim1$b_i[pred_TF == 1], #kg
                     nrow = nrow(grid), 
                     ncol = length(unique(yrs)))
     
     
    #biomass (kg) to CPUE (kg/km2)
    sim_dens<-sim_bio/grid$Area_in_survey_km2 
    
    #save data
    save(sim_dens, file = paste0("./output/species/",sp,'/simulated historical data/',sim_fold,'/sim_dens.RData')) 
    #load(file = paste0("./output/species/Gadus macrocephalus//simulated historical data/0001//sim_dens.RData"))  
    
    #checking units ---- kg/km2
    # x<-readRDS('./data processed/species/Gadus macrocephalus/data_geostat_temp.rds') #kg/km2
    # summary(x$cpue_kgkm2)
    # sim_dens1<-reshape2::melt(sim_dens/grid$Area_in_survey_km2)
    # summary(sim_dens/grid$Area_in_survey_km2)
    # summary(sim_dens1$value)
    
    #loop over sampling designs
    for (samp in unique(samp_df$samp_scn)) { #sampling designs
      
      #get station locations for each sampling design
      load(file=paste0('./output/species/',sp,'/optimization data/samples_optimization_',samp,'_dynamic.RData')) #all_points
      
      #to store results
      sim_survey <- array(data = NA, dim = c(length(dimnames(all_points)[[1]]),
                                             length(fit$year_labels)+1, #add strata number
                                             #length(dimnames(all_points)[[3]]),
                                             2),
                          dimnames = list(dimnames(all_points)[[1]],
                                          c(fit$year_labels,'strata'),
                                          #dimnames(all_points)[[3]],
                                          c('current','random')))
      
      for (y in 1:length(1982:2022)) { #years
        
        #y<-dimnames(all_points)[[3]][1]
        
        #year
        yy<-c(1982:2022)[y]
        
        #get points iterations
        pointsc<-data.frame(unlist(all_points[,,y,'current']))
        #pointsb<-data.frame(unlist(all_points[,,y,'buffer']))                          
        pointsr<-data.frame(unlist(all_points[,,y,'random']))   
        
        #append survey densities for each iteration and simulated data
        sim_survey[,,'current']<-cbind(sim_dens[pointsc$cell,],pointsc$strata)
        #sim_survey[,,'buffer']<-cbind(sim_dens[pointsb$cell,],pointsb$strata)
        sim_survey[,,'random']<-cbind(sim_dens[pointsr$cell,],pointsr$strata)
        
      }
      
      #store results
      save(sim_survey, file = paste0("./output/species/",sp,'/survey simulation historical data/',sim_fold,'/sim_survey_',samp,'_dynamic.RData'))  
    }
  }
}
#}

#remove objects
rm(fit);rm(Sim1)

######################
# PROJECTED DATA
######################

#loop over spp
for (sp in spp) {
  
  sp<-"Gadus macrocephalus"
  
  #create folder simulation number
  dir.create(paste0('./output/species/',sp,'/survey simulation projected data/'))
  
  #loop over SBT scenarios
  for (sbt in paste0('SBT',1:12)) {
      
    #sbt<-'SBT1'
    
    #print scenario to check progress
    cat(paste(" #############   Species", sp, match(sp,spp), 'out of',length(spp),  "  #############\n",
              " #############  projected simulation", sbt, "of", 'SBT12' , " #############\n"))
  
  
    #load sbt data
    load(paste0('./output/species/',sp,'/simulated projected data/fit_projection_',sbt,'.RData'))
    
    #loop over simulations  
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
        
        for (samp in unique(samp_df$samp_scn)) { #sampling designs
          
          #samp<-"scn1"
          
          #get station locations for each sampling design
          load(file=paste0('./output/species/',sp,'/optimization data/samples_optimization_',samp,'_dynamic.RData')) #all_points
            
          #to store results
          sim_survey <- array(data = NA, dim = c(length(dimnames(all_points)[[1]]),
                                                 5+1, #add strata number
                                                 2),
                              dimnames = list(dimnames(all_points)[[1]],
                                              c(2023:2027,'strata'),
                                              c('current','random')))
            
          for (y in 1:length(2022:2027)) { #years
            
            #y<-1  
              
            #years
            yy<-c(2022:2027)[y]
              
            #get points iterations
            pointsc<-data.frame(unlist(all_points[,,y,'current']))
            #pointsb<-data.frame(unlist(all_points[,,y,'buffer']))   
            pointsr<-data.frame(unlist(all_points[,,y,'random']))                          
            
            #append survey densities for each iteration and simulated data
            sim_survey[,,'current']<-cbind(sim_dens[pointsc$cell,],pointsc$strata)
            #sim_survey[,,'buffer']<-cbind(sim_dens[pointsb$cell,],pointsb$strata)
            sim_survey[,,'random']<-cbind(sim_dens[pointsr$cell,],pointsr$strata)

          }
            
          #store results
          save(sim_survey, file = paste0("./output/species/",sp,'/survey simulation projected data/',sim_fold,'/sim_survey_',samp,'_',sbt,'_dynamic.RData'))  
          gc();rm(sim_survey)
        }
      }
   gc();rm(pm)
   }
}
    
######################
# PROJECTED DATA TEST
######################

#project just one 


  