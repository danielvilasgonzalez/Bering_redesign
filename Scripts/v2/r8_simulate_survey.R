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
# LOAD GRID EBS (remember to keep the same order as in fit_model if multiple grids)
###################################

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)
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
# SIMULATED DENSITIES
###################################

#load data
load(file = paste0("./output/species/sim_hist_dens_spp.RData"))    #sim_hist_dens_spp
sim_hist_dens<-sim_hist_dens_spp
rm(sim_hist_dens_spp)

#load data
load(file = paste0("./output/species/sim_proj_dens_spp.RData"))  #sim_proj_dens_spp
sim_proj_dens<-sim_proj_dens_spp
rm(sim_proj_dens_spp)

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

    #loop over sampling designs
    for (samp in unique(samp_df$samp_scn)[4:5]) { #sampling designs
      
      #samp<-unique(samp_df$samp_scn)[1]
      
      s<-match(samp,samp_df$samp_scn)
      
      if (grepl('base',samp)) { #if it contains base is baseline scenario
        
        strata<-as.data.frame(baseline_strata$cell_strata)[,c('cell','Lat','Lon','Stratum')]
        names(strata)[4]<-c('strata')
        D8<-strata
        
        #allocations<-
        if(samp_df[s,'samp_scn']=='scnbase'){
          allocations<-data.frame('strata'=baseline_strata$n_samples$stratum,
                                  'n_samples'=baseline_strata$n_samples$scnbase)
        }else{
          allocations<-data.frame('strata'=baseline_strata$n_samples$stratum,
                                  'n_samples'=baseline_strata$n_samples$scnbase_bis)
        }
        
      } else {
        
        #load results_optimization
        load(file=paste0('./output/ms_optim_allocations_',samp,'.RData')) #list = c('result_list','ss_sample_allocations','ms_sample_allocations','samples_strata','cv_temp')
        load(file=paste0('./output/species/multisp_optimization_static_data.RData')) #df
        df<-df[,c("Lat",'Lon','cell')]
        
        #cell - strata
        strata<-all$result_list$solution$indices #all$result_list$sol_by_cell
        colnames(strata)<-c('cell','strata')
        strata<-strata[order(strata$cell),]

        #add a strata value to each cell
        D8<-merge(df,strata,by='cell',all.x=TRUE)
        D8<-D8[,c("cell","Lat","Lon","strata")]
        names(D8)[4]<-'strata'
        D8$strata<-as.numeric(D8$strata)
        D8$strata<-ifelse(is.na(D8$strata),999,D8$strata)
        
        #strata - n_samples
        allocations<-data.frame('strata'=all$samples_strata$strata,
                                'n_samples'=all$samples_strata$n_samples)
        
      }
      
      #to store all points
      all_points<-array(NA,
                        dim = list(sum(allocations$n_samples),4,length(yrs),2,n_sim),
                        dimnames = list(c(1:sum(allocations$n_samples)),c('Lon','Lat','cell','strata'),1:length(yrs),c('systematic','random'),1:n_sim))
      
      #to store survey samples
      sim_survey <- array(data = NA, dim = c(length(dimnames(all_points)[[1]]),
                                             length(spp)+2,
                                             length(yrs), #add strata number
                                             #length(dimnames(all_points)[[3]]),
                                             2,
                                             n_sim),
                          dimnames = list(dimnames(all_points)[[1]],
                                          c(spp,'strata','cell'),
                                          c(yrs),
                                          #dimnames(all_points)[[3]],
                                          c('systematic','random'),
                                          1:n_sim))
      
      # #loop over simulations
      for (isim in 1:n_sim) { #simulations
        #   
        #isim<-1
          
        #print simulation to check progress
        cat(paste(" #############  historical - sampling design",samp ,'- simulation ', isim, "of",n_sim, " #############\n"))
           
        #simulation folder
        #sim_fold<-formatC(isim, width = 4, format = "d", flag = "0")

      for (y in 1:length(1982:2022)) { #years
        
        #y<-1
        
        #year
        yy<-c(1982:2022)[y]
        
        #to store points
        dfcurrent<-data.frame(matrix(NA,nrow=0,ncol=4))
        colnames(dfcurrent)<-c('Lon','Lat','cell','strata')
        dfrandom<-dfcurrent
        
        #for while purposes
        flag<-TRUE
        
        #random sample for each strata wit distance constrains
        for(n_istrata in 1:nrow(allocations)) {
          
          #n_istrata<-1
          
          istrata<-allocations[n_istrata,'strata']
          
          #subset cells for strata
          df<-subset(as.data.frame(D8),strata==istrata)
          df1<-subset(as.data.frame(D8),strata==istrata & cell %in% baseline_strata$locations$cell)
          n_i<-allocations[n_istrata,'n_samples']
          
          #while error, keep running
          while(flag){
            
            #keep running if error
            tryCatch({
              
              #print loop state
              #cat(paste(" -- strata", n_istrata, 'of',nrow(allocations),  " -- CURRENT  --\n"))
              
              ##############################
              # CURRENT STRATIFIED APPROACH - SYSTEMATIC
              ##############################
              
              if(samp_df[s,'samp_scn']=='scnbase'){
                pointsc<-baseline_strata$locations[which(baseline_strata$locations$stratum==istrata),c('longitude','latitude','cell','stratum')]
                names(pointsc)<-c('Lon' ,'Lat', 'cell','strata')
              }else if (samp_df[s,'samp_scn']=='scnbase_bis'){
                pointsc<-baseline_strata$locations[which(baseline_strata$locations$corner==FALSE & baseline_strata$locations$stratum==istrata),c('longitude','latitude','cell','stratum')]
                names(pointsc)<-c('Lon' ,'Lat',  'cell','strata')
              } else {
                
                #if more required samples than available
                if (nrow(df1)<n_i) {
                  
                  #vectors to store results
                  #dropcell<-c()
                  #selcell<-c()
                  
                  #duplicate df strata
                  dff<-df
                  
                  #df removing available samples from current design
                  dff<-subset(dff, !(cell %in% df1$cell))
                  
                  #cells to complete the required cells
                  ii<-n_i-nrow(df1)
                  
                  #loop over the required samples using random
                  selcell<-sample(dff$cell,ii)
                  
                  
                  #subset points if equal to required number of samples, if not ERROR and rerun iteration
                  if (nrow(subset(df,cell %in% c(selcell,df1$cell))) != n_i) {flag<-TRUE; stop("-- different number of points on current approach",call. = FALSE)}else{
                    pointsc<-subset(df,cell %in% c(selcell,df1$cell))[,c('Lon','Lat','cell','strata')]
                    names(pointsc)<-c('Lon','Lat','cell','strata')}
                  
                  #else there are enough samples to get from the current sampling design    
                } else {
                  ss<-sample(1:nrow(df1),size = n_i,replace = FALSE)
                  pointsc<-df1[ss,c('Lon','Lat','cell','strata')]}}

              ##############################
              # STRATIFIED RANDOM APPROACH
              ##############################
              
              #cat(paste(" -- strata", n_istrata, 'of',nrow(allocations),  " -- RANDOM  --\n"))
              
              #duplicate df strata
              dff<-df
              
              #random selection of samples
              selcell<-sample(dff$cell,n_i)
              pointsr<-subset(df,cell %in% selcell)[,c('Lon','Lat','cell','strata')]
              names(pointsr)<-c('Lon','Lat','cell','strata')
              
              #append data if buffer and current df have equal to required number of samples
              if (nrow(pointsc) == n_i & nrow(pointsr) == n_i) {
                dfcurrent<-rbind(dfcurrent,pointsc)
                #dfbuffer<-rbind(dfbuffer,pointsb)
                dfrandom<-rbind(dfrandom,pointsr)
                flag<-FALSE}
            },
            
            #error message
            error=function(e) {
              message(paste0("ERROR RETRYING",':\n'),e)})
            if (!flag) next
          }
          
          #to start while again
          flag<-TRUE
        }
        
        #append points
        all_points[,,y,'systematic',isim]<-unlist(dfcurrent)
        #all_points[,,y,'buffer']<-unlist(dfbuffer)
        all_points[,,y,'random',isim]<-unlist(dfrandom)
        
        #get points iterations
        pointsc<-data.frame(unlist(all_points[,,y,'systematic',isim]))
        #pointsb<-data.frame(unlist(all_points[,,y,'buffer']))                          
        pointsr<-data.frame(unlist(all_points[,,y,'random',isim]))   

        #append survey densities for each iteration and simulated data
        sim_survey[,,as.character(yy),'systematic',isim]<-cbind(sim_hist_dens[pointsc$cell,y,isim,],'strata'=pointsc$strata,'cell'=pointsc$cell)
        #sim_survey[,,'buffer']<-cbind(sim_dens[pointsb$cell,],pointsb$strata)
        sim_survey[,,as.character(yy),'random',isim]<-cbind(sim_hist_dens[pointsr$cell,y,isim,],'strata'=pointsr$strata,'cell'=pointsr$cell)
        
      }
    }
    #store results
    save(sim_survey, file = paste0('./output/species/sim_hist_survey_',samp,'.RData'))  
    save(all_points, file = paste0('./output/species/allocations_hist_survey_',samp,'.RData'))     
      
  }

######################
# PROJECTED DATA
######################

    #loop over sampling designs
    for (samp in unique(samp_df$samp_scn)) { #sampling designs
      
      #samp<-unique(samp_df$samp_scn)[1]
      
      s<-match(samp,samp_df$samp_scn)
      
      if (grepl('base',samp)) { #if it contains base is baseline scenario
        
        strata<-as.data.frame(baseline_strata$cell_strata)[,c('cell','Lat','Lon','Stratum')]
        names(strata)[4]<-c('strata')
        D8<-strata
        
        #allocations<-
        if(samp_df[s,'samp_scn']=='scnbase'){
          allocations<-data.frame('strata'=baseline_strata$n_samples$stratum,
                                  'n_samples'=baseline_strata$n_samples$scnbase)
        }else{
          allocations<-data.frame('strata'=baseline_strata$n_samples$stratum,
                                  'n_samples'=baseline_strata$n_samples$scnbase_bis)
        }
        
      } else {
        
        #load results_optimization
        load(file=paste0('./output/ms_optim_allocations_',samp,'.RData')) #list = c('result_list','ss_sample_allocations','ms_sample_allocations','samples_strata','cv_temp')
        load(file=paste0('./output/species/multisp_optimization_static_data.RData')) #df
        df<-df[,c("Lat",'Lon','cell')]
        
        #cell - strata
        strata<-all$result_list$solution$indices #all$result_list$sol_by_cell
        colnames(strata)<-c('cell','strata')
        strata<-strata[order(strata$cell),]
        
        #add a strata value to each cell
        D8<-merge(df,strata,by='cell',all.x=TRUE)
        D8<-D8[,c("cell","Lat","Lon","strata")]
        names(D8)[4]<-'strata'
        D8$strata<-as.numeric(D8$strata)
        D8$strata<-ifelse(is.na(D8$strata),999,D8$strata)
        
        #strata - n_samples
        allocations<-data.frame('strata'=all$samples_strata$strata,
                                'n_samples'=all$samples_strata$n_samples)
        
      }
      
      #to store all points
      all_points<-array(NA,
                        dim = list(sum(allocations$n_samples),4,length(project_yrs),2,n_sim),
                        dimnames = list(c(1:sum(allocations$n_samples)),c('Lon','Lat','cell','strata'),1:length(project_yrs),c('systematic','random'),1:n_sim))
      
      #to store survey samples
      sim_survey <- array(data = NA, dim = c(length(dimnames(all_points)[[1]]),
                                             length(spp)+2,
                                             length(project_yrs), #add strata number
                                             #length(dimnames(all_points)[[3]]),
                                             2,
                                             nrow(df_sbt),
                                             n_sim),
                          dimnames = list(dimnames(all_points)[[1]],
                                          c(spp,'strata','cell'),
                                          c(project_yrs),
                                          #dimnames(all_points)[[3]],
                                          c('systematic','random'),
                                          df_sbt$sbt,
                                          1:n_sim))

      
      # #loop over simulations
      for (isim in 1:n_sim) { #simulations
        #   
        #isim<-1
        
        #print simulation to check progress
        cat(paste(" #############  projected - sampling design",samp ,'- simulation ', isim, "of",n_sim, " #############\n"))
        
        #simulation folder
        #sim_fold<-formatC(isim, width = 4, format = "d", flag = "0")
        
        
        for (y in 1:length(project_yrs)) { #years
          
          #y=1
          #year
          yy<-c(project_yrs)[y]
          
          #to store points
          dfcurrent<-data.frame(matrix(NA,nrow=0,ncol=4))
          colnames(dfcurrent)<-c('Lon','Lat','cell','strata')
          dfrandom<-dfcurrent
          
          #for while purposes
          flag<-TRUE
          
          #random sample for each strata wit distance constrains
          for(n_istrata in 1:nrow(allocations)) {
            
            #n_istrata<-1
            
            istrata<-allocations[n_istrata,'strata']
            
            #subset cells for strata
            df<-subset(as.data.frame(D8),strata==istrata)
            df1<-subset(as.data.frame(D8),strata==istrata & cell %in% baseline_strata$locations$cell)
            n_i<-allocations[n_istrata,'n_samples']
            
            #while error, keep running
            while(flag){
              
              #keep running if error
              tryCatch({
                
                #print loop state
                #cat(paste(" -- strata", n_istrata, 'of',nrow(allocations),  " -- CURRENT  --\n"))
                
                ##############################
                # CURRENT STRATIFIED APPROACH
                ##############################
                
                if(samp_df[s,'samp_scn']=='scnbase'){
                  pointsc<-baseline_strata$locations[which(baseline_strata$locations$stratum==istrata),c('longitude','latitude','cell','stratum')]
                  names(pointsc)<-c('Lon' ,'Lat', 'cell','strata')
                }else if (samp_df[s,'samp_scn']=='scnbase_bis'){
                  pointsc<-baseline_strata$locations[which(baseline_strata$locations$corner==FALSE & baseline_strata$locations$stratum==istrata),c('longitude','latitude','cell','stratum')]
                  names(pointsc)<-c('Lon' ,'Lat',  'cell','strata')
                } else {
                  
                  #if more required samples than available
                  if (nrow(df1)<n_i) {
                    
                    #vectors to store results
                    #dropcell<-c()
                    #selcell<-c()
                    
                    #duplicate df strata
                    dff<-df
                    
                    #df removing available samples from current design
                    dff<-subset(dff, !(cell %in% df1$cell))
                    
                    #cells to complete the required cells
                    ii<-n_i-nrow(df1)
                    
                    #loop over the required samples using random
                    selcell<-sample(dff$cell,ii)
                    
                    
                    #subset points if equal to required number of samples, if not ERROR and rerun iteration
                    if (nrow(subset(df,cell %in% c(selcell,df1$cell))) != n_i) {flag<-TRUE; stop("-- different number of points on current approach",call. = FALSE)}else{
                      pointsc<-subset(df,cell %in% c(selcell,df1$cell))[,c('Lon','Lat','cell','strata')]
                      names(pointsc)<-c('Lon','Lat','cell','strata')}
                    
                    #else there are enough samples to get from the current sampling design    
                  } else {
                    ss<-sample(1:nrow(df1),size = n_i,replace = FALSE)
                    pointsc<-df1[ss,c('Lon','Lat','cell','strata')]}}
                
                
                ##############################
                # STRATIFIED RANDOM APPROACH
                ##############################
                
                #cat(paste(" -- strata", n_istrata, 'of',nrow(allocations),  " -- RANDOM  --\n"))
                
                #duplicate df strata
                dff<-df
                
                #random selection of samples
                selcell<-sample(dff$cell,n_i)
                pointsr<-subset(df,cell %in% selcell)[,c('Lon','Lat','cell','strata')]
                names(pointsr)<-c('Lon','Lat','cell','strata')
                
                #append data if buffer and current df have equal to required number of samples
                if (nrow(pointsc) == n_i & nrow(pointsr) == n_i) {
                  dfcurrent<-rbind(dfcurrent,pointsc)
                  #dfbuffer<-rbind(dfbuffer,pointsb)
                  dfrandom<-rbind(dfrandom,pointsr)
                  flag<-FALSE}
              },
              
              #error message
              error=function(e) {
                message(paste0("ERROR RETRYING",':\n'),e)})
              if (!flag) next
            }
            
            #to start while again
            flag<-TRUE
          }
          
          #append points
          all_points[,,y,'systematic',isim]<-unlist(dfcurrent)
          #all_points[,,y,'buffer']<-unlist(dfbuffer)
          all_points[,,y,'random',isim]<-unlist(dfrandom)
          
          #get points iterations
          pointsc<-data.frame(unlist(all_points[,,y,'systematic',isim]))
          #pointsb<-data.frame(unlist(all_points[,,y,'buffer']))                          
          pointsr<-data.frame(unlist(all_points[,,y,'random',isim]))   
          
          
          #loop over SBT scenarios
          for (sbt in paste0(df_sbt$sbt_n)) {
            
            #dimnames(sim_proj_dens)
            sim_proj_dens1<-sim_proj_dens[,,1,sbt,]
            
          #append survey densities for each iteration and simulated data
          sim_survey[,,as.character(yy),'systematic',paste0('SBT',sbt),isim]<-cbind(sim_proj_dens1[pointsc$cell,y,],'strata'=pointsc$strata,'cell'=pointsc$cell)
          #sim_survey[,,'buffer']<-cbind(sim_dens[pointsb$cell,],pointsb$strata)
          sim_survey[,,as.character(yy),'random',paste0('SBT',sbt),isim]<-cbind(sim_proj_dens1[pointsr$cell,y,],'strata'=pointsr$strata,'cell'=pointsr$cell)
          }
          
        }
      }
      #store results
      save(sim_survey, file = paste0('./output/species/sim_proj_survey_',samp,'.RData'))  
      save(all_points, file = paste0('./output/species/allocations_proj_survey_',samp,'.RData'))  
      
    }


