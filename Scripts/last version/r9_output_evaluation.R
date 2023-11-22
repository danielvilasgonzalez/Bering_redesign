####################################################################
####################################################################
##    
##    plot comparison evaluation for spp
##    
##    danielvilasgonzalez@gmail.com/dvilasg@uw.edu
##
####################################################################
####################################################################

######################
# SETTINGS
######################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c('raster','ggplot2')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

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
       'Paralithodes camtschaticus',
       'Chionoecetes bairdi')

#remove Anoploma and Reinhardtius because habitat preference reasons
spp<-setdiff(spp, c('Anoplopoma fimbria','Reinhardtius hippoglossoides'))

#years
yrs<-c(1982:2022)
n_yrs<-length(yrs)


#how manyt projected years we want
n_proj<-5

#project_yrs
project_yrs<-((yrs[length(yrs)])+1):(yrs[length(yrs)]+n_proj)

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)

#load baseline strata data
load('./output/baseline_strata.RData') #baseline_strata

#add percent of total area per strata
baseline_strata$strata_areas$pct<-baseline_strata$strata_areas$Area_in_survey_km2/sum(baseline_strata$strata_areas$Area_in_survey_km2)

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
samp_df<-rbind(samp_df,c('baseline','systematic',520,15,'scnbase'),
               c('baseline w/o corner','systematic',494,15,'scnbase_bis'))

###################################
# SBT scenarios
###################################

#load SBT scenarios table
load('./tables/SBT_projection.RData')#df_sbt

#name SBT scenarios
df_sbt$sbt<-paste0('SBT',df_sbt$sbt_n)
df_sbt$sbt2<-paste0(df_sbt$sbt,'_',df_sbt$Scenario)

#number of historical simulations and projected simulations
n_sim<- 100

######################
# HISTORICAL 
######################

#INDEX TRUE (MODEL$BASED)
#load(file = paste0("./output/species/dens_index_hist_OM.RData"))  #dens_index_hist_OM, 

#array to store
index_hist<-array(NA,
                  dim=list(length(spp),3,n_yrs,2,n_sim,nrow(samp_df)),
                  dimnames=list(spp,c('STRS_mean','STRS_var','CV_sim'),paste0('y',yrs),c('systematic','random'),1:n_sim,samp_df$samp_scn))

#loop over sampling designs
for (samp in samp_df$samp_scn) { #sampling designs
  
  #samp<-unique(samp_df$samp_scn)[1]
  
  s<-match(samp,samp_df$samp_scn)
  
  #store results
  load(file = paste0('./output/species/sim_hist_survey_',samp,'.RData'))  #sim_survey

  #when base sampling other files
  if (grepl('base',samp)) {
    
    #conditions on baseline scenarios
    if (samp == 'scnbase') {
      baseline_strata$locations2<-baseline_strata$locations
    } else if (samp == 'scnbase_bis') {
      baseline_strata$locations2<-baseline_strata$locations[which(baseline_strata$locations$corner=='FALSE'),]
    } 
    
    #sort by strata
    baseline_strata$strata_areas<-baseline_strata$strata_areas[order(baseline_strata$strata_areas$X1),]
    baseline_strata$locations2<-baseline_strata$locations2[order(baseline_strata$locations2$stratum),]
    
    #area by strata
    strata_areas <- baseline_strata$strata_areas
    
    #strata data
    survey_detail <- data.frame("Stratum" = baseline_strata$strata_areas$X1, #strata
                                'Nh' = baseline_strata$strata_areas$pct*53464, #number of cells
                                "nh" = data.frame(table(baseline_strata$locations2$stratum))[,c('Freq')]) #number of sample allocations
    
    #weight of strata for each
    survey_detail$Wh <- survey_detail$Nh / sum(survey_detail$Nh)
    survey_detail$wh <- with(survey_detail, nh/Nh)
    
  } else {
    
    #load optimization results
    load(paste0("./output/ms_optim_allocations_",samp_df[s,'samp_scn'],".RData")) #all
    
    #area
    area_cell<-merge(all$result_list$solution$indices, grid, by.x='ID',by.y='cell')
    
    #area by strata
    strata_areas <- aggregate(Area_in_survey_km2 ~ X1, 
                              FUN = sum,
                              data = area_cell)
    
    #strata data
    survey_detail <- data.frame("Stratum" = all$samples_strata$strata, #strata
                                'Nh' = as.integer(table(all$result_list$solution$indices$X1)), #number of cells
                                "nh" = all$samples_strata$n_samples) #number of sample allocations
    
    #weight of strata for each
    survey_detail$Wh <- survey_detail$Nh / sum(survey_detail$Nh)
    survey_detail$wh <- with(survey_detail, nh/Nh)
    
  }   
  
  #loop over simulation (n_sim=100)
  for (isim in 1:n_sim) {
    
    #isim<-1
    
    #print scenario to check progress
    cat(paste(" #############  sampling design -", samp, " #############\n",
              " #############  simulations -", isim, " #############\n"))

    #loop over approaches  
    for (apr in dimnames(sim_survey)[[4]]) {
      
      #apr<-dimnames(sim_survey)[[4]][1]
      
      #loop over approaches  
      for (year in dimnames(sim_survey)[[3]]) {
          
        #get strata mean density over year  
        dens<-unlist(sim_survey[,,year,apr,isim])
        
        y<-data.frame(dens[,c(spp,'strata')],check.names = FALSE)
        yy<-reshape2::melt(y,id.vars=c('strata'))
        
        #mean, sum and var by strata and year (variable)
        yyy<-aggregate(x=yy$value,
                       by=list(strata=yy$strata,sp=yy$variable),
                       FUN = function(x) c('mean' = mean(x,na.rm=T), 'sum' = sum(x),'var' = var(x,na.rm=T) ))
        
        #create df
        zzz<-data.frame('strata'=yyy$strata,'sp'=yyy$sp,'mean'=yyy$x[,c('mean')],'var'=yyy$x[,c('var')]) #/length(yy$value)
        zzzz<-merge(zzz,strata_areas,by.x='strata',by.y='X1',all.x=TRUE)
        zzzz<-merge(zzzz,survey_detail,by.x='strata',by.y='Stratum',all.x=TRUE)
        
        #add index strata for sum to compute index (mean strata density * area of strata) kg!
        zzzz$index_strata<-zzzz$mean*zzzz$Area_in_survey_km2
        
        #add strata var 
        zzzz$strs_var<-zzzz$var*(zzzz$Area_in_survey_km2^2)/zzzz$nh #sum(survey_detail$Nh) 
        
        #sum of strata var and mean density across years (kg/km2)
        zzzz1 <- aggregate(zzzz[,c('strs_var','index_strata')], by= list(zzzz$sp),FUN = sum)
        
        #get CV across years
        zzzz1$cv<- sqrt(zzzz1$strs_var) / zzzz1$index_strata
        
        #mean CV 
        mean(zzzz1$cv,na.rm=TRUE)
        
        #get outputs
        STRS_mean <- zzzz1$index_strata
        STRS_var <- zzzz1$strs_var
        CV <- sqrt(STRS_var) / STRS_mean
        
        #store outputs
        index_hist[,'STRS_mean',paste0('y',year),apr,isim,samp]<-STRS_mean
        index_hist[,'STRS_var',paste0('y',year),apr,isim,samp]<-STRS_var
        index_hist[,'CV_sim',paste0('y',year),apr,isim,samp]<-CV
      }
    }
  }
}
      
save(index_hist, file = paste0("./output/species/historical design estimates.RData" ))
  
######################
# PROJECTED
######################

#INDEX TRUE (MODEL$BASED)
#load(file = paste0("./output/species/dens_index_project_OM.RData"))  #dens_index_hist_OM, 

n_sim<-1
n_sur<-100

#array to store
index_proj<-array(NA,
                  dim=list(length(spp),3,n_proj,2,length(1:n_sim),length(1:n_sur),nrow(df_sbt),nrow(samp_df)),
                  dimnames=list(spp,c('STRS_mean','STRS_var','CV_sim'),paste0('y',project_yrs),c('systematic','random'),1:n_sim,1:n_sur,df_sbt$sbt,samp_df$samp_scn))

 
#loop over sampling designs
for (samp in samp_df$samp_scn) {
  
  #samp<-samp_df$samp_scn[3]
  
  s<-match(samp,samp_df$samp_scn)
  
  #print scenario to check progress
  cat(paste(" #############  sampling design -", samp, " #############\n"))
  
  #store results
  load(file = paste0('./output/species/sim_proj_survey_',samp,'.RData'))  #sim_survey
  
  #when base sampling other files
  if (grepl('base',samp)) {
    
    #conditions on baseline scenarios
    if (samp == 'scnbase') {
      baseline_strata$locations2<-baseline_strata$locations
    } else if (samp == 'scnbase_bis') {
      baseline_strata$locations2<-baseline_strata$locations[which(baseline_strata$locations$corner=='FALSE'),]
    } 
  
    #sort by strata
    baseline_strata$strata_areas<-baseline_strata$strata_areas[order(baseline_strata$strata_areas$X1),]
    baseline_strata$locations2<-baseline_strata$locations2[order(baseline_strata$locations2$stratum),]
    
    #area by strata
    strata_areas <- baseline_strata$strata_areas
    
    #strata data
    survey_detail <- data.frame("Stratum" = baseline_strata$strata_areas$X1, #strata
                                'Nh' = baseline_strata$strata_areas$pct*53464, #number of cells
                                "nh" = data.frame(table(baseline_strata$locations2$stratum))[,c('Freq')]) #number of sample allocations
    
    #weight of strata for each
    survey_detail$Wh <- survey_detail$Nh / sum(survey_detail$Nh)
    survey_detail$wh <- with(survey_detail, nh/Nh)
    
  } else {
  
    #load optimization results
    load(paste0("./output/ms_optim_allocations_",samp_df[s,'samp_scn'],".RData")) #all
    
    #area
    area_cell<-merge(all$result_list$solution$indices, grid, by.x='ID',by.y='cell')
    
    #area by strata
    strata_areas <- aggregate(Area_in_survey_km2 ~ X1, 
                              FUN = sum,
                              data = area_cell)
    
    #strata data
    survey_detail <- data.frame("Stratum" = all$samples_strata$strata, #strata
                                'Nh' = as.integer(table(all$result_list$solution$indices$X1)), #number of cells
                                "nh" = all$samples_strata$n_samples) #number of sample allocations
    
    #weight of strata for each
    survey_detail$Wh <- survey_detail$Nh / sum(survey_detail$Nh)
    survey_detail$wh <- with(survey_detail, nh/Nh)
  }  
  
  #loop over sbt scenarios
  for (sbt in df_sbt$sbt) {
    
    #print scenario to check progress
    cat(paste(" #############  SBT projection -", sbt, " #############\n"))
    
    #sbt<-'SBT1'
    #loop over simulation
    for (isim in 1:n_sim) {
      
    #loop over simulation
    for (isur in 1:n_sur) {
      
      #isim<-1
    
      #loop over approaches
      for (apr in dimnames(sim_survey)[[4]]) {
        
        #apr<-dimnames(sim_survey)[[4]][1]
        
        #loop over sbt scenarios
        for (year in project_yrs) {
          
          #year<-'2023'
          
          #get strata mean over year  
          dens<-unlist(sim_survey[,,as.character(year),apr,sbt,isim,isur])
          y<-data.frame(dens[,c(spp,'strata')],check.names = FALSE)
          yy<-reshape2::melt(y,id.vars=c('strata'))
          
          #mean, sum and var by strata and year (variable)
          yyy<-aggregate(x=yy$value,
                         by=list(strata=yy$strata,sp=yy$variable),
                         FUN = function(x) c('mean' = mean(x,na.rm=T), 'sum' = sum(x),'var' = var(x,na.rm=T) ))
          
          
          
          
          #create df
          zzz<-data.frame('strata'=yyy$strata,'sp'=yyy$sp,'mean'=yyy$x[,c('mean')],'var'=yyy$x[,c('var')]) #/length(yy$value)
          zzzz<-merge(zzz,strata_areas,by.x='strata',by.y='X1',all.x=TRUE)
          zzzz<-merge(zzzz,survey_detail,by.x='strata',by.y='Stratum',all.x=TRUE)
          
          #add index strata for sum to compute index (mean strata density * area of strata) kg!
          zzzz$index_strata<-zzzz$mean*zzzz$Area_in_survey_km2
          
          #add strata var 
          zzzz$strs_var<-zzzz$var*(zzzz$Area_in_survey_km2^2)/zzzz$nh #sum(survey_detail$Nh) 
          
          #sum of strata var and mean density across years (kg/km2)
          zzzz1 <- aggregate(zzzz[,c('strs_var','index_strata')], by= list(zzzz$sp),FUN = sum)
          
          #get CV across years
          zzzz1$cv<- sqrt(zzzz1$strs_var) / zzzz1$index_strata
          
          #mean CV - to check 
          mean(zzzz1$cv,na.rm=TRUE) 
          
          #get outputs
          STRS_mean <- zzzz1$index_strata #absolute strata year biomass
          STRS_var <- zzzz1$strs_var #strata variance 
          CV <- sqrt(STRS_var) / STRS_mean #simulated CV
          
          #append results
          index_proj[,'STRS_mean',paste0('y',year),apr,isim,isur,sbt,samp]<-STRS_mean
          index_proj[,'STRS_var',paste0('y',year),apr,isim,isur,sbt,samp]<-STRS_var
          index_proj[,'CV_sim',paste0('y',year),apr,isim,isur,sbt,samp]<-CV
        }
        }
      }  
    }
  }
}
  
#save object
save(index_proj, file = paste0("./output/species/projected design estimates.RData" ))

#comparison estimates proj and hist
#CV (STRS_var/index)
mean(index_hist[,'CV_sim',,,,],na.rm=TRUE)
mean(index_proj[,'CV_sim',,,,,,],na.rm=TRUE)
#index
mean(index_hist[,'STRS_mean',,,,])
mean(index_proj[,'STRS_mean',,,,,,])
