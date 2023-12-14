####################################################################
####################################################################
##    
##    simulate data and survey for historical and projected years
##    prepare estimates to compute design-based indices
##    danielvilasgonzalez@gmail.com/dvilasg@uw.edu
##    
##    systematic / spatially balanced / random
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 

#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c('ggplot2','units','splines','raster','sp','Spbsampling','sf','doParallel','foreach')

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
       'Paralithodes camtschaticus',
       'Chionoecetes bairdi')


#remove Anoploma and Reinhardtius because habitat preference reasons
spp<-setdiff(spp, c('Anoplopoma fimbria','Reinhardtius hippoglossoides'))

#yrs
yrs<-setdiff(1982:2022,2020)

#how manyt projected years we want
n_proj<-5

#project_yrs
project_yrs<-((yrs[length(yrs)])+1):(yrs[length(yrs)]+n_proj)

###################################
# LOAD GRID EBS (remember to keep the same order as in fit_model if multiple grids)
###################################

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)
grid2<-grid
#add col and row number
x1<-grid[,c('Lon','Lat','cell')]
names(x1)<-c('x','y','z')
coordinates(x1)=~x + y
crs(x1)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
x2<-spTransform(x1,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
x3<-data.frame(x2)
x3$x<-as.integer(x3$coords.x1)
x3$y<-as.integer(x3$coords.x2)
lon<-sort(unique(x3$x),decreasing = FALSE) #1556
lat<-sort(unique(x3$y),decreasing = TRUE) #1507
lons<-data.frame(x=lon,col=1:length(lon))
lats<-data.frame(y=lat,row=1:length(lat))
x4<-merge(x3,lons,by='x',all.x=TRUE)
x5<-merge(x4,lats,by='y',all.x=TRUE)
x5<-x5[,c('y','x','z','col','row')]
colnames(x5)<-c('Lat','Lon','cell','col','row')
grid<-x5[,c('Lat','Lon','cell','col','row')]

###################################
# YEARS AND BASELINE
###################################

#load grid
load('./data processed/grid_EBS_NBS.RData')
#yrs<-setdiff(1982:2022,2020)
grid_ebs<-grid.ebs_year[which(grid.ebs_year$region != 'EBSslope' & grid.ebs_year$Year %in% yrs),]
dim(grid_ebs)

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
# BASELINE STRATA
###################################

load('./output/baseline_strata.RData')
#add percent of total area per strata
baseline_strata$strata_areas$pct<-baseline_strata$strata_areas$Area_in_survey_km2/sum(baseline_strata$strata_areas$Area_in_survey_km2)

###################################
# SBT scenarios
###################################

#load SBT scenarios table
load('./tables/SBT_projection.RData')#df_sbt

#name SBT scenarios
df_sbt$sbt<-paste0('SBT',df_sbt$sbt_n)
df_sbt$sbt2<-paste0(df_sbt$sbt,'_',df_sbt$Scenario)

#number of historical simulations and projected simulations
n_sim_hist<- 100
n_sim_proj<- 100
 
#number of surveys
n_sur<-100

################
# HISTORICAL
################
 
#create a df of 40years x 100sur/sim (it changes on the projected since it simulations have different indeces)
n_sur<-100
sur_df<-data.frame(num=1:(length(yrs)*n_sur),
           year=rep(yrs,times=n_sur),
           sur=rep(1:n_sur,each=length(yrs)))

#simulated densities
load(file = paste0('./output/species/ms_sim_dens.RData'))  #sim_dens1

#ms_sim_survey folder
dir.create('./output/ms_sim_survey/')

#create folder
dir.create(paste0('./output/ms_sim_survey_hist/'))

#loop over n combinations of simulated
for (sim in 1:n_sim_hist) {
  
  #sim<-1
  
  # Convert 0 to '001'
  sim_fol <- sprintf("%03d", sim)
  
  #create folder
  dir.create(paste0('./output/ms_sim_survey_hist/sim',sim_fol))
  
  #array to store
  index_hist<-array(NA,
                    dim=list(length(spp),3,length(yrs),3,n_sur,nrow(samp_df)),
                    dimnames=list(spp,c('STRS_mean','STRS_var','CV_sim'),paste0('y',yrs),c('sys','rand','sb'),1:n_sur,samp_df$samp_scn))
  
  #loop over sampling design
    for (samp in samp_df$samp_scn)  {
      
      #samp<-'scnbase_bis'
      #start_time_parallel <- Sys.time()
      
      #number of sampling design
      s<-match(samp,samp_df$samp_scn)
      
      #when base sampling other files
      if (grepl('base',samp)) {
        
        #conditions on baseline scenarios
        if (samp == 'scnbase') {
          baseline_strata$locations2<-baseline_strata$locations
        } else if (samp == 'scnbase_bis') {
          baseline_strata$locations2<-baseline_strata$locations[which(baseline_strata$locations$corner=='FALSE'),]
        } 
        
        #sort by strata
        #baseline_strata$strata_areas<-baseline_strata$strata_areas[order(baseline_strata$strata_areas$X1),]
        #baseline_strata$locations2<-baseline_strata$locations2[order(baseline_strata$locations2$stratum),]
        names(baseline_strata$strata_areas)[1]<-'X1'
        names(baseline_strata$locations2)[ncol(baseline_strata$locations2)]<-'stratum'
        
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
        area_cell<-merge(all$result_list$solution$indices, grid2, by.x='ID',by.y='cell')
        
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
      
      #array to store simulated densities/CPUE
      alloc<-ifelse(samp=='scnbase_bis',494,520)
      
      #load survey allocations by sampling design
      load(file = paste0('./output/survey_allocations_',samp,'.RData')) #scn_allocations
      dimnames(scn_allocations)[[3]]<-c('sys','rand','sb')
      
        #array to store results  
        sim_survey <- array(NA,
                            dim = c(alloc, length(spp)+2, n_sur, length(unique(yrs)),length(c('sys','rand','sb'))),
                            dimnames = list(1:alloc, c('cell','strata',spp), 1:n_sur,unique(yrs), c('sys','rand','sb')))
        
        
        
        #loop over n combinations of simulated
        for (n in sur_df$num) {
    
          #n<-sur_df$num[1]
          
          #year of simulation
          y<-as.character(sur_df[which(sur_df$num==n),'year'])
          #isurvey of simulation
          sur<-sur_df[which(sur_df$num==n),'sur']
          
          #simulated densities of survey and year
          sim_dens2<-sim_dens1[,,y,sim]
          
          #loop over station allocation approac
          for (apr in c('sys','rand','sb')) {
            
            #apr<-'sys'
    
            #print process        
            cat(paste(" #############  ",samp,'- simdata',sim,'- survey',sur, '- year',y ,'- allocation',apr," #############\n"))
            
            #get densities based on station allocations
            sim_survey<-data.frame(cbind(strata=scn_allocations[scn_allocations[,'sur',apr]==n,c('strata'),apr],
                                 dens=sim_dens2[scn_allocations[scn_allocations[,'sur',apr]==n,c('cell'),apr],]),check.names = FALSE)
            
            sim_survey1<-reshape2::melt(sim_survey,id.vars=c('strata'))
            
            #mean, sum and var by strata and year (variable)
            sim_survey2<-aggregate(x=sim_survey1$value,
                           by=list(strata=sim_survey1$strata,sp=sim_survey1$variable),
                           FUN = function(x) c('mean' = mean(x,na.rm=T), 'sum' = sum(x),'var' = var(x,na.rm=T) ))
            
            #create df
            zzz<-data.frame('strata'=sim_survey2$strata,'sp'=sim_survey2$sp,'mean'=sim_survey2$x[,c('mean')],'var'=sim_survey2$x[,c('var')]) #/length(yy$value)
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
            index_hist[,'STRS_mean',paste0('y',y),apr,sur,samp]<-STRS_mean
            index_hist[,'STRS_var',paste0('y',y),apr,sur,samp]<-STRS_var
            index_hist[,'CV_sim',paste0('y',y),apr,sur,samp]<-CV
        
          }
        }
  }
  
  save(index_hist, file = paste0('./output/ms_sim_survey_hist/sim',sim_fol,'/index_hist.RData'))  

}

################
# PROJECTED
################

#dir create sim survey projected
dir.create('./output/ms_sim_survey_proj')

#create a df of 5years x 100sur/sim x x 8sbt scenarios (it changes on the projected since it simulations have different indeces)
n_sur<-100

sur_df<-cbind(expand.grid('year'=project_yrs,'sur'=1:n_sim_proj,'sbt'=df_sbt$sbt_n),'num'=1:4000)


#simulated densities
#load(file = paste0('./output/species/ms_sim_dens.RData'))  #sim_dens1

#ms_sim_survey folder
#dir.create('./output/ms_sim_survey/')

#loop over n combinations of simulated
for (sim in 1:n_sim_proj) {
  
  #sim<-1
  
  # Convert 0 to '001'
  sim_fol <- sprintf("%03d", sim)
  
  #create folder
  dir.create(paste0('./output/ms_sim_survey_proj/sim',sim_fol))

  
#loop over projections
for (sbt in df_sbt$sbt_n) {
  
  #sbt<-df_sbt$sbt_n[1]
    
  #num surveys
  nums<-sur_df[which(sur_df$sbt==sbt),'num']
  
  #sbt<-df_sbt$sbt_n[1]
  #load densities projections
  load(paste0('./output/species/SBT',sbt,' ms_sim_proj_dens.RData'))
  
  #array to store
  index_proj<-array(NA,
                    dim=list(length(spp),3,length(project_yrs),3,n_sur,nrow(samp_df)),
                    dimnames=list(spp,c('STRS_mean','STRS_var','CV_sim'),paste0('y',project_yrs),c('sys','rand','sb'),1:n_sur,samp_df$samp_scn))
  
  
  #array to store results  
  # sim_survey <- array(NA,
  #                     dim = c(alloc, length(spp)+2, n_sur, length(unique(project_yrs)),length(samp_df$samp_scn),length(c('sys','rand','sb'))),
  #                     dimnames = list(1:alloc, c('cell','strata',spp), 1:n_sur,unique(project_yrs), samp_df$samp_scn, c('sys','rand','sb')))
  
  
  #loop over sampling design
  for (samp in samp_df$samp_scn)  {
    
    #samp<-'scnbase'
    #start_time_parallel <- Sys.time()
    
    #number of sampling design
    s<-match(samp,samp_df$samp_scn)
    
    #when base sampling other files
    if (grepl('base',samp)) {
      
      #conditions on baseline scenarios
      if (samp == 'scnbase') {
        baseline_strata$locations2<-baseline_strata$locations
      } else if (samp == 'scnbase_bis') {
        baseline_strata$locations2<-baseline_strata$locations[which(baseline_strata$locations$corner=='FALSE'),]
      } 
      
      #sort by strata
      #baseline_strata$strata_areas<-baseline_strata$strata_areas[order(baseline_strata$strata_areas$X1),]
      #baseline_strata$locations2<-baseline_strata$locations2[order(baseline_strata$locations2$stratum),]
      names(baseline_strata$strata_areas)[1]<-'X1'
      names(baseline_strata$locations2)[ncol(baseline_strata$locations2)]<-'stratum'
      
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
      area_cell<-merge(all$result_list$solution$indices, grid2, by.x='ID',by.y='cell')
      
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
    
    #array to store simulated densities/CPUE
    alloc<-ifelse(samp=='scnbase_bis',494,520)
    
    #load survey allocations by sampling design
    load(file = paste0('./output/survey_allocations_',samp,'.RData')) #scn_allocations
    dimnames(scn_allocations)[[3]]<-c('sys','rand','sb')
      
    #loop over n combinations of simulated
    for (n in nums) {
      
      #n<-nums[1]
      
      #year of simulation
      y<-as.character(sur_df[which(sur_df$num==n),'year'])
      #isurvey of simulation
      sur<-sur_df[which(sur_df$num==n),'sur']
      
      #simulated densities of survey and year
      sim_dens2<-simdata[,,y,sim]
      
      #loop over station allocation approac
      for (apr in c('sys','rand','sb')) {
        
        #apr<-'sys'
        
        #print process        
        cat(paste(" #############  ",samp,'- simdata',sim,'- sbt',sbt,'- survey',sur, '- year',y ,'- allocation',apr," #############\n"))
        
        #get densities based on station allocations
        sim_survey<-data.frame(cbind(strata=scn_allocations[scn_allocations[,'sur',apr]==n,c('strata'),apr],
                                     dens=sim_dens2[scn_allocations[scn_allocations[,'sur',apr]==n,c('cell'),apr],]),check.names = FALSE)
        
        sim_survey1<-reshape2::melt(sim_survey,id.vars=c('strata'))
        
        #mean, sum and var by strata and year (variable)
        sim_survey2<-aggregate(x=sim_survey1$value,
                               by=list(strata=sim_survey1$strata,sp=sim_survey1$variable),
                               FUN = function(x) c('mean' = mean(x,na.rm=T), 'sum' = sum(x),'var' = var(x,na.rm=T) ))
        
        #create df
        zzz<-data.frame('strata'=sim_survey2$strata,'sp'=sim_survey2$sp,'mean'=sim_survey2$x[,c('mean')],'var'=sim_survey2$x[,c('var')]) #/length(yy$value)
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
        index_proj[,'STRS_mean',paste0('y',y),apr,sur,samp]<-STRS_mean
        index_proj[,'STRS_var',paste0('y',y),apr,sur,samp]<-STRS_var
        index_proj[,'CV_sim',paste0('y',y),apr,sur,samp]<-CV
        
      }
     }
    }
  save(index_proj, file = paste0('./output/ms_sim_survey_proj/sim',sim_fol,'/SBT',sbt,' index_proj.RData'))  
  
  } 
}



# ################
# # get simulated densities from simulated survey PROJECTED
# ################
# 
# #all 4000 surveys (year5, sur100, sbt8)
# sur_df<-cbind(num=1:nrow(expand.grid(project_yrs,df_sbt$sbt,1:n_sur)),expand.grid(year=project_yrs,sbt=df_sbt$sbt_n,sur=1:n_sur))
# 
# #load sbt scenarios
# load('./tables/SBT_projection.RData') #df_sbt
# 
# #loop over species
# for (sp in spp) {
#   
#   #sp<-spp[1]
#   
#   #loop over 8 temperature scenarios
#   for (sbt in unique(df_sbt$sbt_n)) {
#     
#     #sbt<-df_sbt$sbt_n[1]
#     
#     cat(paste(" #############  spp",sp,'- sbt',sbt ," #############\n"))
#     
#     
#     load(file = paste0("./output/species/",sp,"/simulated projected data/SBT",sbt," dens_index_proj_OM_50.RData")) #dens_index_proj_OM
#     dens_index_proj_OM_50<-dens_index_proj_OM
#     load(file = paste0("./output/species/",sp,"/simulated projected data/SBT",sbt," dens_index_proj_OM_100.RData")) #dens_index_proj_OM
#     dens_index_proj_OM_100<-dens_index_proj_OM
#     rm(dens_index_proj_OM)
#     
#       #for each sampling design
#       for (samp in samp_df$samp_scn) {
#         
#         #samp<-samp_df$samp_scn[1]
#         
#         #load survey allocations by sampling design
#         load(file = paste0('./output/species/ms_sim_survey/survey_allocations_',samp,'.RData')) #scn_allocations
#         
#         #to store survey samples
#         sim_survey <- array(data = NA, dim = c(dim(scn_allocations)[1]/n_sur/length(yrs),
#                                                length(c("Lon","Lat","cell","strata","sur"))+1,
#                                                n_sim_proj,
#                                                length(c('sys','rand','sb')),
#                                                n_sur,
#                                                length(project_yrs)),
#                             dimnames = list(c(1:(dim(scn_allocations)[1]/n_sur/length(yrs))),
#                                             c("Lon","Lat","cell","strata","sur",'CPUE_kgkm'),
#                                             1:n_sim_proj,
#                                             c('sys','rand','sb'),
#                                             1:n_sur,
#                                             as.character(project_yrs)))
#       
#       #loop over each of 100 simulated projections
#       for (isim in 1:n_sim_proj) {
#       
#         cat(paste(" #############  simulated data", isim," #############\n"))
#         
#         #isim<-1
#         
#         if (isim %in% 1:50) {
#           dens_index_proj_OM1<-data.frame(drop_units(dens_index_proj_OM_50[[paste0('sim',isim)]]$dens),check.names = FALSE)
#         } else if (isim %in% 51:100) {
#           dens_index_proj_OM1<-data.frame(drop_units(dens_index_proj_OM_100[[paste0('sim',isim-50)]]$dens),check.names = FALSE)
#         }
#         
#         # #for 100 simulated surveys in each year
#         # for (y in project_yrs) {
#         #   
#         #   y<-2023
#         #   dd<-dens_index_proj_OM1[,as.character(y)]
#           
#         dim(dens_index_proj_OM1)
#         dd<-dens_index_proj_OM1[,as.character(project_yrs)]
#     
#           #for each allocation sampling approach
#           for (ap in c('sys','rand','sb')) {
#             
#             #ap<-'systematic'
#               
#             #get allocations by sampling allocation approach
#             scn_allocations1<-data.frame(scn_allocations[,,ap])
#             
#                 #loop over 100 surveys
#                 for (isur in 1:n_sur) {
#                   
#                   #get survey number 1:4000 (5y * 8 sbt * 100sur)
#                   sur_num<-sur_df[which(sur_df$sbt==sbt & sur_df$sur==isur),] #'num'
#                   
#                   for (y in as.character(project_yrs)) {
#                     
#                     #y<-project_yrs[1]
#                     
#                     scn_allocations2<-scn_allocations1[which(scn_allocations1$sur == sur_num[which(sur_num$year==y),'num']),]
#                     dd1<-dd[scn_allocations2$cell,as.character(y)]
# 
#                     sim_survey[,,isim,ap,isur,y]<-cbind(scn_allocations2,CPUE_kgkm=dd1)
#                     
#                   }
#                 }
#               }
#             }
#         #store results
#         save(sim_survey, file = paste0('./output/species/',sp,'/simulated survey project/sim_proj_survey_sbt',sbt,'_',samp,'.RData'))  
#         #save(all_points, file = paste0('./output/species/allocations_hist_survey_',samp,'.RData'))     
#         
#         rm(sim_survey)
#       }
#     }
#   }              
   
# ################
# # PROJECTED
# ################
# 
# #for each sampling design
# for (samp in samp_df$samp_scn) {
#     
#     #samp<-samp_df$samp_scn[1]
#     
#     #load survey allocations by sampling design
#     load(file = paste0('./output/species/ms_sim_survey/survey_allocations_',samp,'.RData')) #scn_allocations
#     
#     #to store survey samples
#     sim_survey <- array(data = NA, dim = c(dim(scn_allocations)[1]/n_sur/length(yrs),
#                                            length(c("Lon","Lat","cell","strata","sur"))+1,
#                                            n_sim_proj,
#                                            length(c('sys','rand','sb')),
#                                            n_sur,
#                                            length(project_yrs)),
#                         dimnames = list(c(1:(dim(scn_allocations)[1]/n_sur/length(yrs))),
#                                         c("Lon","Lat","cell","strata","sur",'CPUE_kgkm'),
#                                         1:n_sim_proj,
#                                         c('sys','rand','sb'),
#                                         1:n_sur,
#                                         as.character(project_yrs)))
#     
#     #loop over each of 100 simulated projections
#     for (isim in 1:n_sim_proj) {
#       
#       cat(paste(" #############  simulated data", isim," #############\n"))
#       
#       #isim<-1
#       
#       if (isim %in% 1:50) {
#         dens_index_proj_OM1<-data.frame(drop_units(dens_index_proj_OM_50[[paste0('sim',isim)]]$dens),check.names = FALSE)
#       } else if (isim %in% 51:100) {
#         dens_index_proj_OM1<-data.frame(drop_units(dens_index_proj_OM_100[[paste0('sim',isim-50)]]$dens),check.names = FALSE)
#       }
#       
#       # #for 100 simulated surveys in each year
#       # for (y in project_yrs) {
#       #  
#       #   y<-2023
#       #   dd<-dens_index_proj_OM1[,as.character(y)]
#       
#       dim(dens_index_proj_OM1)
#       dd<-dens_index_proj_OM1[,as.character(project_yrs)]
#       
#       #for each allocation sampling approach
#       for (ap in c('sys','rand','sb')) {
#         
#         #ap<-'systematic'
#         
#         #get allocations by sampling allocation approach
#         scn_allocations1<-data.frame(scn_allocations[,,ap])
#         
#         #loop over 100 surveys
#         for (isur in 1:n_sur) {
#           
#           #get survey number 1:4000 (5y * 8 sbt * 100sur)
#           sur_num<-sur_df[which(sur_df$sbt==sbt & sur_df$sur==isur),] #'num'
#           
#           for (y in as.character(project_yrs)) {
#             
#             #y<-project_yrs[1]
#             
#             scn_allocations2<-scn_allocations1[which(scn_allocations1$sur == sur_num[which(sur_num$year==y),'num']),]
#             dd1<-dd[scn_allocations2$cell,as.character(y)]
#             
#             sim_survey[,,isim,ap,isur,y]<-cbind(scn_allocations2,CPUE_kgkm=dd1)
#             
#           }
#         }
#       }
#     }
#     #store results
#     save(sim_survey, file = paste0('./output/species/',sp,'/simulated survey project/sim_proj_survey_sbt',sbt,'_',samp,'.RData'))  
#     #save(all_points, file = paste0('./output/species/allocations_hist_survey_',samp,'.RData'))    
#     
#     rm(sim_survey)
#   }
                  