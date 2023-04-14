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
  
######################
# GET ANNUAL PREDICTIONS (2023-2027) SBT-
######################

#loop over species
for (sp in spp) {
  
  sp<-"Gadus macrocephalus"
  
  #get sampling designs
  lf<-list.files(path = paste0('./output/species/',sp,'/'),pattern = 'samples_optimization_')
  lf1<-gsub('samples_optimization_','',lf)
  lf2<-gsub('.RData','',lf1)
  
  #load projection fit
  load( file = paste0("./output/species/",sp,'/fit_projection.RData')) #pr_list
  
  #loop over sampling designs
  for (samp in lf2) {
    
    #samp<-lf2[1]
    
    #to store results
    samp_proj<-list()
    
    #get station locations for each sampling design
    load(file=paste0('./output/species/',sp,'/samples_optimization_',samp,'.RData')) #all_points
    
    #loop over SBT scenarios
    for (sbt in names(pr_list)) {
      
      #sbt<-names(pr_list)[1]
      
      #print scenario to check progress
      cat(paste(" #############     PROJECTING    #############\n",
                " #############   Species", sp, match(sp,spp), 'out of',length(spp),  "  #############\n",
                " #############  Sampling ", samp, " #############\n",
                " #############  SBT ", sbt, " #############\n"))
   
      
      #get raster stack for each SBT scenario
      bio_proj<-stack(paste0('./output/species/',sp,'/biomass_projection_',sbt,'.grd'))
      names(bio_proj)<-paste0('y',c(2023:2027))
      
      #to store points
      points2<-data.frame(matrix(nrow = 0,ncol = 10))
      colnames(points2)<-c(dimnames(all_points)[[2]],names(bio_proj),'iter')
      
      #loop over iterations
      for (iter in dimnames(all_points)[[3]]) {
        
        #iter<-'50'
        
        #subset by iteration
        points<-data.frame(all_points[,,iter])
  
        #convert to spatial object
        coordinates(points)<- ~ Lon + Lat
        
        #extract sp biomass projections for each sampling design and iteration
        points1<-data.frame(data.frame(all_points[,,iter]),raster::extract(bio_proj,points))
        points1$iter<-iter
        points2<-rbind(points2,points1)
      }
    
      #append biomass samples
      samp_proj[[paste0('SBT',sbt)]]<-points2
    }
    #save biomass samples list
    save(samp_proj, file = paste0("./output/species/",sp,'/projections_survey_',samp,'.RData'))
  }
}

#array to store
index_array<-array(NA,
                  dim=list(100,4,5,length(names(samp_proj)),length(lf2),length(spp)),
                  dimnames=list(1:100,c('STRS_mean','STRS_var','cv','index'),paste0('y',2023:2027),names(samp_proj),lf2,spp))

#loop over species
for (sp in spp) {
  
  sp<-"Gadus macrocephalus"
  
  #get samples scenarios
  lf<-list.files(path = paste0('./output/species/',sp,'/'),pattern = 'samples_optimization_')
  lf1<-gsub('samples_optimization_','',lf)
  lf2<-gsub('.RData','',lf1)

  #load stats optimization summary for weighting and area
  load(paste0('./output/species/',sp,'/optimization_summary_stats.RData')) #sp_sum_stats,
  
  #loop over sampling designs
  for (samp in lf2) {
    
    #samp<-lf2[1]
    
    #load survey samples
    load( file = paste0("./output/species/",sp,'/projections_survey_',samp,'.RData')) #samp_proj
    
    #load optimization results
    load(paste0('./output/species/',sp,'/optimization_results_',samp,'.RData')) #result_list

    #subset by sampling design
    sum_stats<-subset(sp_sum_stats,samp_scn==samp)
    
    #loop over sbt scenarios  
    for (sbt in names(samp_proj)) {
      
      #sbt<-names(samp_proj)[1]
      
      #print scenario to check progress
      cat(paste(" #############   Species", sp, match(sp,spp), 'out of',length(spp),  "  #############\n",
                " #############  Sampling Scenario", samp, " #############\n",
                " #############  SBT Scenario", sbt, " #############\n"))
      
      #get densities from sbt scenario
      x<-samp_proj[[sbt]]
      
      #merge biomass samples of survey with grid
      y<-merge(x,grid,by='cell')
      
      # yy<-subset(y,iter==1)
      # 
      # ggplot()+
      #   geom_point(data=yy,aes(x=Lon.x,y=Lat.x,size=y2023,color=y2023))
      # 
      
      #loop over years
      for (year in paste0('y',2023:2027)) {
      
        #year<-'y2023'
        
        #catch to CPUE
        #cpue<-y[,year]/y[,"Area_in_survey_km2"] 
        y$cpue<-y[,year]*y[,"Area_in_survey_km2"] 
        cpue<-y[,year]*y[,"Area_in_survey_km2"] 
        
         # ggplot()+
        #   geom_point(data=yy,aes(x=Lon.x,y=Lat.x,size=cpue,color=cpue))
        
        
        #strata data
        survey_detail <- data.frame("Stratum" = result_list$solution$aggr_strata$STRATO, #strata
                                    'Nh' = as.integer(table(result_list$solution$indices$X1)), #number of cells
                                    "nh" = result_list$sample_allocations) #number of sample allocations
        
        #weight of strata for each
        survey_detail$Wh <- survey_detail$Nh / sum(survey_detail$Nh)
        survey_detail$wh <- with(survey_detail, nh/Nh)
  
        #mean by strata
        strata_means<-aggregate(x=cpue,
                                by=list(iter=y$iter,strata=y$strata),
                                FUN=mean, na.rm=T)
        #var by strata
        strata_vars<-aggregate(x=cpue,
                                by=list(iter=y$iter,strata=y$strata),
                                FUN=var, na.rm=T)
        
        #area by strata
        strata_areas <- aggregate(Area_in_survey_km2 ~ strata + iter, 
                                  FUN = sum,
                                  data = y)
      
        #loop over iterations
        for (iters in sort(as.integer(unique(strata_areas$iter)))) {
          
          #iters<-1
          
          #subset by iter
          strata_means1<-subset(strata_means,iter==iters)
          strata_areas1<-subset(strata_areas,iter==iters)
          strata_vars1<-subset(strata_vars,iter==iters)
          
          #calculate STRS mean and variance of density and CV
          index_array[iters,'STRS_mean',year,sbt,samp,sp] <- STRS_mean <- sum(strata_means1$x * survey_detail$Wh)
          index_array[iters,'STRS_var',year,sbt,samp,sp] <- STRS_var <- sum(strata_vars1$x * with(survey_detail, Wh^2 * (1 - wh) / nh) )
          index_array[iters,'cv',year,sbt,samp,sp] <- sqrt(STRS_var) / STRS_mean
          
          # Calculate total index
          index_array[iters,'index',year,sbt,samp,sp] <- sum(strata_areas1$Area_in_survey_km2 * strata_means1$x) * 0.001
        }
      }
    }
  }
  #save object
  sp_index<-unlist(index_array[,,,,,sp])
  save(sp_index, file = paste0("./output/species/",sp,'/indices.RData'))
}
    
#save object
save(index_array, file = paste0("./output/indices.RData"))
