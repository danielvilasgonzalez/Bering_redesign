####################################################################
####################################################################
##    
##    simulate survey by extracting predicted densities for
##    each location under each sampling design for each SBT projection
##    get design-based estimates for each year and iteration
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
pack_cran<-c('raster')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

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

#years
yrs<-c(1982:2022)
n_yrs<-length(yrs)

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)

#load baseline strata data
load('./output/baseline_strata.RData') #baseline_strata

#add percent of total area per strata
baseline_strata$strata_areas$pct<-baseline_strata$strata_areas$Area_in_survey_km2/sum(baseline_strata$strata_areas$Area_in_survey_km2)

######################
# DESIGN-BASED INDEX HISTORICAL DATA
######################

#loop over species
for (sp in spp) {
  
  sp<-"Gadus macrocephalus"
  
  #create folder simulation number
  dir.create(paste0('./output/species/',sp,'/'))
  
  #get samples scenarios
  lf<-list.files(path = paste0('./output/species/',sp,'/survey simulation historical data/0001/'),pattern = 'sim_survey')
  lf<-lf[grep('_dynamic',lf)]
  lf1<-gsub('sim_survey_','',lf)
  lf2<-gsub('_dynamic.RData','',lf1)
  
  #folder simulations
  ld<-list.dirs(path = paste0('./output/species/',sp,'/survey simulation historical data/'),full.names = FALSE,recursive = FALSE)
  
  #load stats optimization summary for weighting and area
  load(paste0('./output/species/',sp,'/optimization data/optimization_summary_stats.RData')) #sp_sum_stats,
  
  #array to store
  index_hist<-array(NA,
                    dim=list(4,41,2,length(ld),length(lf2)),
                    dimnames=list(c('STRS_mean','STRS_var','cv','index'),paste0('y',1982:2022),c('current','random'),ld,lf2))
  
  
  #loop over simulation (n_sim=100)
  for (sim in ld) {
    
    #sim<-ld[1]
    
    #loop over sampling designs
    for (samp in lf2) {
      
      #samp<-lf2[4]
      
      #print scenario to check progress
      cat(paste(" #############  ", sim, " #############\n",
                " #############  ", samp, " #############\n"))
      
      #load survey samples from historical
      load( file = paste0("./output/species/",sp,'/survey simulation historical data/',sim,'/sim_survey_',samp,'_dynamic.RData')) #sim_survey
      
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
        load(paste0('./output/species/',sp,'/optimization data/optimization_results_',samp,'.RData')) #result_list
        
        #area
        area_cell<-merge(result_list$solution$indices, grid, by.x='ID',by.y='cell')
        
        #area by strata
        strata_areas <- aggregate(Area_in_survey_km2 ~ X1, 
                                  FUN = sum,
                                  data = area_cell)
        
        #strata data
        survey_detail <- data.frame("Stratum" = result_list$solution$aggr_strata$STRATO, #strata
                                    'Nh' = as.integer(table(result_list$solution$indices$X1)), #number of cells
                                    "nh" = result_list$sample_allocations) #number of sample allocations
        
        #weight of strata for each
        survey_detail$Wh <- survey_detail$Nh / sum(survey_detail$Nh)
        survey_detail$wh <- with(survey_detail, nh/Nh)
        
      }   
      
      #loop over approaches  
      for (apr in dimnames(sim_survey)[[3]]) {
        
        #apr<-dimnames(sim_survey)[[3]][1]
        
        #get strata mean density over year  
        dens<-unlist(sim_survey[,,apr])
        y<-data.frame(dens)
        rowMeans(y[,-ncol(y)])
        yy<-reshape2::melt(y,id.vars=c('strata'))
        #yy$value<-yy$value/1000 #t/km²
        
        #mean, sum and var by strata and year (variable)
        yyy<-aggregate(x=yy$value,
                       by=list(strata=yy$strata,year=yy$variable),
                       FUN = function(x) c('mean' = mean(x,na.rm=T), 'sum' = sum(x),'var' = var(x,na.rm=T) ))
        
        #yyy$x[,c('mean')] #strata_mean
        #yyy$x[,c('var')]/length(yy$value) #strata var
        
        #create df
        zzz<-data.frame('strata'=yyy$strata,'year'=yyy$year,'mean'=yyy$x[,c('mean')],'var'=yyy$x[,c('var')]) #/length(yy$value)
        zzzz<-merge(zzz,strata_areas,by.x='strata',by.y='X1',all.x=TRUE)
        zzzz<-merge(zzzz,survey_detail,by.x='strata',by.y='Stratum',all.x=TRUE)
        #add index strata for sum to compute index (mean strata density * area of strata) t!
        zzzz$index_strata<-zzzz$mean*zzzz$Area_in_survey_km2
        #add strata var 
        zzzz$strs_var<-zzzz$var*(zzzz$Area_in_survey_km2^2)/zzzz$nh #sum(survey_detail$Nh) 
        #sum of strata var and mean density across years (kg/km2)
        zzzz1 <- aggregate(zzzz[,c('strs_var','index_strata')], by= list(zzzz$year) , 
                           FUN = sum)
        #get CV across years
        zzzz1$cv<- sqrt(zzzz1$strs_var) / zzzz1$index_strata
        #mean CV 
        mean(zzzz1$cv)
        
        #get outputs
        STRS_mean <- zzzz1$index_strata
        STRS_var <- zzzz1$strs_var
        CV <- sqrt(STRS_var) / STRS_mean
        index<- zzzz1$index_strata
        
        index_hist['STRS_mean',,apr,sim,samp]<-STRS_mean
        index_hist['STRS_var',,apr,sim,samp]<-STRS_var
        index_hist['cv',,apr,sim,samp]<-CV
        index_hist['index',,apr,sim,samp]<-index
      } 
    }
  }
  #save object
  save(index_hist, file = paste0("./output/species/",sp,'/historical design-based indices/indices.RData'))
}

######################
# DESIGN-BASED INDEX PROJECTED DATA
######################

#loop over species
for (sp in spp) {
  
  sp<-"Gadus macrocephalus"
  
  #create folder simulation number
  dir.create(paste0('./output/species/',sp,'/'))
  
  #get samples scenarios
  lf<-list.files(path = paste0('./output/species/',sp,'/survey simulation projected data/0001/'),pattern = 'sim_survey')
  lf<-lf[grepl('_dynamic',lf)]
  lf0<-lf[grepl('SBT',lf)]
  lf1<-gsub('sim_survey_','',lf0)
  lf2<-unique(sub("\\_.*", "", lf1))
  lf2<-c('scn1','scn2','scn3','scnbase','scnbase_bis')
  
  #folder simulations
  ld<-list.dirs(path = paste0('./output/species/',sp,'/survey simulation projected data/'),full.names = FALSE,recursive = FALSE)
  
  #load stats optimization summary for weighting and area
  load(paste0('./output/species/',sp,'/optimization data/optimization_summary_stats.RData')) #sp_sum_stats,
  
  #array to store
  index_proj<-array(NA,
                    dim=list(4,5,2,length(ld),length(lf2),length(paste0('SBT',1:12))),
                    dimnames=list(c('STRS_mean','STRS_var','cv','index'),paste0('y',2023:2027),c('current','random'),ld,lf2,paste0('SBT',1:12)))
  
  #loop over simulation
  for (sim in ld) {
    
    #sim<-ld[1]
    
    #loop over sbt scenarios
    for (sbt in paste0('SBT',1:12)) {
      
      #sbt<-'SBT1'
      
      #loop over sampling designs
      for (samp in lf2) {
        
        #samp<-lf2[1]
        #samp<-lf2[5]
        
        #print scenario to check progress
        cat(paste(" #############  ", sim, " #############\n",
                  " #############  ", sbt, " #############\n",
                  " #############  ", samp, " #############\n"))
        
        #load survey samples from historical
        load( file = paste0("./output/species/",sp,'/survey simulation projected data/',sim,'/sim_survey_',samp,'_',sbt,'_dynamic.RData')) #sim_survey
        
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
          load(paste0('./output/species/',sp,'/optimization data/optimization_results_',samp,'.RData')) #result_list
          
          #area
          area_cell<-merge(result_list$solution$indices, grid, by.x='ID',by.y='cell')
          
          #area by strata
          strata_areas <- aggregate(Area_in_survey_km2 ~ X1, 
                                    FUN = sum,
                                    data = area_cell)
          
          #strata data
          survey_detail <- data.frame("Stratum" = result_list$solution$aggr_strata$STRATO, #strata
                                      'Nh' = as.integer(table(result_list$solution$indices$X1)), #number of cells
                                      "nh" = result_list$sample_allocations) #number of sample allocations
          
          #weight of strata for each
          survey_detail$Wh <- survey_detail$Nh / sum(survey_detail$Nh)
          survey_detail$wh <- with(survey_detail, nh/Nh)
        }  
        
        #loop over approaches
        for (apr in dimnames(sim_survey)[[3]]) {
          
          #apr<-dimnames(sim_survey)[[3]][1]
          
          #get strata mean over year  
          dens<-unlist(sim_survey[,,apr])
          y<-data.frame(dens)
          yy<-reshape2::melt(y,id.vars=c('strata'))
          #yy$value<-yy$value/1000 #t/km²
          
          #mean, sum and var by strata and year (variable)
          yyy<-aggregate(x=yy$value,
                         by=list(strata=yy$strata,year=yy$variable),
                         FUN = function(x) c(mean = mean(x,na.rm=T), sum = sum(x),var = var(x,na.rm=T) ))
          
          
          #yyy$x[,c('mean')] #strata_mean
          #yyy$x[,c('var')]/length(yy$value) #strata var
          
          #create df
          zzz<-data.frame('strata'=yyy$strata,'year'=yyy$year,'mean'=yyy$x[,c('mean')],'var'=yyy$x[,c('var')]) #/length(yy$value)
          zzzz<-merge(zzz,strata_areas,by.x='strata',by.y='X1',all.x=TRUE)
          zzzz<-merge(zzzz,survey_detail,by.x='strata',by.y='Stratum',all.x=TRUE)
          #add index strata for sum to compute index (mean strata density * area of strata) kg!
          zzzz$index_strata<-zzzz$mean*zzzz$Area_in_survey_km2
          #add strata var 
          zzzz$strs_var<-zzzz$var*(zzzz$Area_in_survey_km2^2)/zzzz$nh #sum(survey_detail$Nh) 
          #sum of strata var and mean density across years (kg/km2)
          zzzz1 <- aggregate(zzzz[,c('strs_var','index_strata')], by= list(zzzz$year) , 
                             FUN = sum)
          #get CV across years
          zzzz1$cv<- sqrt(zzzz1$strs_var) / zzzz1$index_strata
          #mean CV 
          mean(zzzz1$cv)
          
          #get outputs
          STRS_mean <- zzzz1$index_strata
          STRS_var <- zzzz1$strs_var
          CV <- sqrt(STRS_var) / STRS_mean
          index<- zzzz1$index_strata
          
          #append results
          index_proj['STRS_mean',,apr,sim,samp,sbt]<-STRS_mean
          index_proj['STRS_var',,apr,sim,samp,sbt]<-STRS_var
          index_proj['cv',,apr,sim,samp,sbt]<-CV
          index_proj['index',,apr,sim,samp,sbt]<-index
          
        }  
      }
    }
  }
  
  #save object
  save(index_proj, file = paste0("./output/species/",sp,'/projected design-based indices/indices.RData'))
  
}

#comparison estimates proj and hist
#CV (STRS_var/index)
mean(index_hist['cv',,,,])
mean(index_proj['cv',,,,,])
#STRS_var
mean(index_hist['STRS_var',,,,])
mean(index_proj['STRS_var',,,,,])
#index
mean(index_hist['index',,,,])
mean(index_proj['index',,,,,])
