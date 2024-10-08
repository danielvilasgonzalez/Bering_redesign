
###################################
# Grid
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

    #parameters
    #grid <- 53464 #row grid
    spp <- spp
    yrs<-c(1982:2019,2021:2022)
    n_sim<-100
    
  
    #########################
    # reshape in parallel array of historical simulations
    #########################
  
  
    library(foreach)
    library(doParallel)
    
    # Initializing parallel backend
    cl <- makeCluster(detectCores()-1)  # Using all available cores
    registerDoParallel(cl)
    
    # Parallelizing the loop
    #foreach(samp = samp_df$samp_scn) %dopar% {
      
      #samp<-'scn1'
      #start_time_parallel <- Sys.time()
      
      #array to store simulated densities/CPUE
      sim_dens1 <- array(NA,
                         dim = c(nrow(grid), length(spp), length(unique(yrs)), n_sim),
                         dimnames = list(1:nrow(grid), spp, unique(yrs), 1:n_sim))
      
      
      foreach(sp = spp) %do% {
        
        #sp<-spp[1]
        
        load(paste0('./output/species/', sp, '/simulated historical data/sim_dens.RData'))
        
        #sim_dens[,as.character(yrs),]
        
        foreach(y = yrs) %:%
          foreach(sim = 1:n_sim) %do% {
            #y<-'1982';sim<-'1'
            
            sim_dens1[, sp, as.character(y), as.character(sim)] <- sim_dens[, as.character(y), as.character(sim)]
          }
      }
      #end_time_parallel <- Sys.time()
      #store results
    #}
      # Stopping the parallel backend
      stopCluster(cl)
      
      
    save(sim_dens1, file = paste0('./output/species/ms_sim_dens.RData'))  
    
    #########################
    # select historical dens based on allocation and simulations
    #########################
    
    
    library(foreach)
    library(doParallel)
    
    # Initializing parallel backend
    cl <- makeCluster(detectCores()-1)  # Using all available cores
    registerDoParallel(cl)
    
    load(file = paste0('./output/species/ms_sim_dens.RData'))  #sim_dens1, 
    
    # Parallelizing the loop
    foreach(samp = samp_df$samp_scn) %dopar% {
      
      #samp<-'scn1'
      #start_time_parallel <- Sys.time()
      
      #array to store simulated densities/CPUE
      alloc<-ifelse(samp=='scnbase_bis',494,520)
      sim_survey <- array(NA,
                         dim = c(alloc, length(spp)+2, length(unique(yrs)), n_sim,length(c('sys','rand','sb'))),
                         dimnames = list(1:alloc, c('cell','strata',spp), unique(yrs), 1:n_sim,c('sys','rand','sb')))
      
      #load survey allocations by sampling design
      load(file = paste0('./output/species/ms_sim_survey/survey_allocations_',samp,'.RData')) #scn_allocations
      dimnames(scn_allocations)[[3]]<-c('sys','rand','sb')
        
      foreach(sur = as.character(unique(sur_df$num)), .combine='c') %dopar% {
        
        #sur<-as.character(unique(sur_df$num)[1])
          
        #systematic
        sim_survey[,,as.character(sur_df[which(sur_df$num==sur),'year']),as.character(sur_df[which(sur_df$num==sur),'sur']),'sys'] <- 
            cbind(scn_allocations[scn_allocations[,'sur','sys']==sur,c('cell','strata'),'sys'],
                  dens=sim_dens1[scn_allocations[scn_allocations[,'sur','sys']==sur,c('strata'),'sys'],,
                                as.character(sur_df[which(sur_df$num==sur),'year']),
                                as.character(sur_df[which(sur_df$num==sur),'sur'])])
       
        #random
        sim_survey[,,as.character(sur_df[which(sur_df$num==sur),'year']),as.character(sur_df[which(sur_df$num==sur),'sur']),'rand'] <- 
          cbind(scn_allocations[scn_allocations[,'sur','rand']==sur,c('cell','strata'),'rand'],
                dens=sim_dens1[scn_allocations[scn_allocations[,'sur','rand']==sur,c('strata'),'rand'],,
                               as.character(sur_df[which(sur_df$num==sur),'year']),
                               as.character(sur_df[which(sur_df$num==sur),'sur'])])
        
        #sb
        sim_survey[,,as.character(sur_df[which(sur_df$num==sur),'year']),as.character(sur_df[which(sur_df$num==sur),'sur']),'sb'] <- 
          cbind(scn_allocations[scn_allocations[,'sur','sb']==sur,c('cell','strata'),'sb'],
                dens=sim_dens1[scn_allocations[scn_allocations[,'sur','sb']==sur,c('strata'),'sb'],,
                               as.character(sur_df[which(sur_df$num==sur),'year']),
                               as.character(sur_df[which(sur_df$num==sur),'sur'])])
            
        }
      save(sim_survey, file = paste0('./output/species/ms_sim_survey/hist_survey_',samp,'.RData'))  
      
      }

    # Stopping the parallel backend
    stopCluster(cl)
  
  
  
  
  
  
  
  
  
  
  
  
  





