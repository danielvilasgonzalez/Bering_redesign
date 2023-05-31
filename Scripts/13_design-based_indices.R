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
  lf1<-gsub('sim_survey','',lf)
  lf2<-gsub('.RData','',lf1)
  
  #folder simulations
  ld<-list.dirs(path = paste0('./output/species/',sp,'/survey simulation historical data/'),full.names = FALSE,recursive = FALSE)
  
  #load stats optimization summary for weighting and area
  load(paste0('./output/species/',sp,'/optimization data/optimization_summary_stats.RData')) #sp_sum_stats,

  #array to store
  index_hist<-array(NA,
                     dim=list(100,4,41,length(ld),length(lf2)),
                     dimnames=list(1:100,c('STRS_mean','STRS_var','cv','index'),paste0('y',1982:2022),ld,lf2))
  
  
  #loop over simulation
  for (sim in ld) {
    
    #sim<-ld[1]
   
    #print scenario to check progress
    #cat(paste(" #############  ", sim, " #############\n"))
    
    #loop over sampling designs
    for (samp in lf2) {
      
      #samp<-lf2[49]
      
      #print scenario to check progress
      cat(paste(" #############  ", sim, " #############\n",
                " #############  ", samp, " #############\n"))
      
      #load survey samples from historical
      load( file = paste0("./output/species/",sp,'/survey simulation historical data/',sim,'/sim_survey',samp,'.RData')) #sim_survey
      
      #when base scenario there is no replicates
      if (grepl('base',samp)) {
        
        #replicates
        repls<-1
        
        #in case scn base need to remove stations
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
        
        #replicates
        repls<-1:100
        
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
        
        #subset by sampling design
        sum_stats<-subset(sp_sum_stats,samp_scn==samp)
        
      }  
      
      #loop over replicates
      for (repl in repls) { #:100
        
        #repl<-1
        
        #get strata mean over year  
        dens<-unlist(sim_survey[,,repl])
        y<-data.frame(dens)
        yy<-reshape2::melt(y,id.vars=c('strata'))
        
        #mean
        yyy<-aggregate(x=yy$value,
                       by=list(strata=yy$strata,year=yy$variable),
                       FUN = function(x) c(mean = mean(x,na.rm=T), var = var(x,na.rm=T) ))
        
        #rename and dataframe
        yyy$year<-gsub('X','',yyy$year) 
        yyyy<-do.call(data.frame, yyy)
        names(yyyy)<-c('strata','year','mean','var')
        
        #merge weight and strata
        yyyyy<-merge(yyyy,survey_detail[,c("Stratum","Wh")],by.x='strata',by.y='Stratum',all.x=TRUE)
        yyyyy<-merge(yyyyy,strata_areas,by.x='strata',by.y='X1',all.x=TRUE)
        yyyyyy<-merge(yyyyy,data.frame('strata'=unique(survey_detail$Stratum),'v'=with(survey_detail, Wh^2 * (1 - wh) / nh)),by='strata',all.x=TRUE)
        
        #get average mean
        yyyyyy$mean_w<-yyyyyy$mean*yyyyyy$Wh
        yyyyyy$index<-yyyyyy$mean*yyyyyy$Area_in_survey_km2
        yyyyyy$var_w<-yyyyyy$var*yyyyyy$v
        
        #area by strata
        yyyyyyy <- aggregate(yyyyyy[,c('mean_w','var_w','index')], by= list(yyyyyy$year) , 
                             FUN = sum)
        
        #get outputs
        STRS_mean <- yyyyyyy$mean_w
        STRS_var <- yyyyyyy$var_w
        CV <- sqrt(STRS_var) / STRS_mean
        index<- yyyyyyy$index/1000
        
        index_hist[repl,'STRS_mean',,sim,samp]<-STRS_mean
        index_hist[repl,'STRS_var',,sim,samp]<-STRS_var
        index_hist[repl,'cv',,sim,samp]<-CV
        index_hist[repl,'index',,sim,samp]<-index
        
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
  lf0<-lf[grepl('SBT',lf)]
  lf1<-gsub('sim_survey_','',lf0)
  lf2<-unique(sub("\\_.*", "", lf1))
  
  #folder simulations
  ld<-list.dirs(path = paste0('./output/species/',sp,'/survey simulation projected data/'),full.names = FALSE,recursive = FALSE)
  
  #load stats optimization summary for weighting and area
  load(paste0('./output/species/',sp,'/optimization data/optimization_summary_stats.RData')) #sp_sum_stats,

  #array to store
  index_proj<-array(NA,
                    dim=list(100,4,5,length(ld),length(lf2),length(paste0('SBT',1:12))),
                    dimnames=list(1:100,c('STRS_mean','STRS_var','cv','index'),paste0('y',2023:2027),ld,lf2,paste0('SBT',1:12)))
  
  #loop over simulation
  for (sim in ld) {
    
    #sim<-ld[1]
    
    #print scenario to check progress
    #cat(paste(" #############  ", sim, " #############\n"))
    
    #loop over sbt scenarios
    for (sbt in paste0('SBT',1:12)) {
    
      #sbt<-'SBT1'
      
      #loop over sampling designs
      for (samp in lf2) {
        
        #samp<-lf2[48]
        
        #print scenario to check progress
        cat(paste(" #############  ", sim, " #############\n",
                  " #############  ", sbt, " #############\n",
                  " #############  ", samp, " #############\n"))
        
        #load survey samples from historical
        load( file = paste0("./output/species/",sp,'/survey simulation projected data/',sim,'/sim_survey_',samp,'_',sbt,'.RData')) #sim_survey

        #when base scenario there is no replicates
        if (grepl('base',samp)) {
          
          #replicates
          repls<-1
          
          if (samp == "scnbase") {
            baseline_strata$locations2<-baseline_strata$locations
          } else if (samp == "scnbase25") { #remove 25 stations
            baseline_strata$locations2<-baseline_strata$locations[which(baseline_strata$locations$rm25==1),]
          } else if (samp == "scnbase9") { #remove 9 stations
            baseline_strata$locations2<-baseline_strata$locations[which(baseline_strata$locations$rm9==1),] }
          
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
          
          #replicates
          repls<-1:100
          
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
          
          #subset by sampling design
          sum_stats<-subset(sp_sum_stats,samp_scn==samp)
          
        }  
        
        #loop over replicates
        for (repl in repls) { #:100
          
          #repl<-1
          
          #get strata mean over year  
          dens<-unlist(sim_survey[,,repl])
          y<-data.frame(dens)
          yy<-reshape2::melt(y,id.vars=c('strata'))
          
          #mean
          yyy<-aggregate(x=yy$value,
                         by=list(strata=yy$strata,year=yy$variable),
                         FUN = function(x) c(mean = mean(x,na.rm=T), var = var(x,na.rm=T) ))
          
          #rename and dataframe
          yyy$year<-gsub('X','',yyy$year) 
          yyyy<-do.call(data.frame, yyy)
          names(yyyy)<-c('strata','year','mean','var')
          
          #merge weight and strata
          yyyyy<-merge(yyyy,survey_detail[,c("Stratum","Wh")],by.x='strata',by.y='Stratum',all.x=TRUE)
          yyyyy<-merge(yyyyy,strata_areas,by.x='strata',by.y='X1',all.x=TRUE)
          yyyyyy<-merge(yyyyy,data.frame('strata'=unique(survey_detail$Stratum),'v'=with(survey_detail, Wh^2 * (1 - wh) / nh)),by='strata',all.x=TRUE)
          
          #get average mean
          yyyyyy$mean_w<-yyyyyy$mean*yyyyyy$Wh
          yyyyyy$index<-yyyyyy$mean*yyyyyy$Area_in_survey_km2
          yyyyyy$var_w<-yyyyyy$var*yyyyyy$v
          
          #area by strata
          yyyyyyy <- aggregate(yyyyyy[,c('mean_w','var_w','index')], by= list(yyyyyy$year) , 
                              FUN = sum)

          #get outputs
          STRS_mean <- yyyyyyy$mean_w
          STRS_var <- yyyyyyy$var_w
          CV <- sqrt(STRS_var) / STRS_mean
          index<- yyyyyyy$index/1000
          
          #append results
          index_proj[repl,'STRS_mean',,sim,samp,sbt]<-STRS_mean
          index_proj[repl,'STRS_var',,sim,samp,sbt]<-STRS_var
          index_proj[repl,'cv',,sim,samp,sbt]<-CV
          index_proj[repl,'index',,sim,samp,sbt]<-index
      
            
          }
        }
      }
    }
  
  #save object
  save(index_proj, file = paste0("./output/species/",sp,'/projected design-based indices/indices.RData'))
  
}
  