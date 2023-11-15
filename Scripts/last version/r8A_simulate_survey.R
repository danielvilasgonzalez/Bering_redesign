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

# ###################################
# # SIMULATED DENSITIES
# ###################################
# 
# #load data
# load(file = paste0("./output/species/sim_hist_dens_spp.RData"))    #sim_hist_dens_spp
# sim_hist_dens<-sim_hist_dens_spp
# rm(sim_hist_dens_spp)
# 
# #load data
# load(file = paste0("./output/species/sim_proj_dens_spp.RData"))  #sim_proj_dens_spp
# sim_proj_dens<-sim_proj_dens_spp
# rm(sim_proj_dens_spp)
# 
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

###############
# get simulated allocations - survey simulation
################
 
#loop over sampling designs
for (samp in unique(samp_df$samp_scn)) { #sampling designs
      
      #samp<-unique(samp_df$samp_scn)[1]
      
      s<-match(samp,samp_df$samp_scn)
      
      n_strata<-samp_df$n_strata[s]
      n_samples<-samp_df$n_samples[s]
      
      if (grepl('base',samp)) { #if it contains base is baseline scenario

        strata<-as.data.frame(baseline_strata$cell_strata)[,c('cell','Lat','Lon','Stratum')]
        names(strata)[4]<-c('strata')
        aggregate(strata$cell,by=list(strata$strata),FUN=length)
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
        load(file=paste0('./output/multisp_optimization_static_data.RData')) #df
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
      
        #load baseline strata and specify corner stations
        load('./output/baseline_strata.RData')
        baseline_strata$locations[grep('GF|HG|JI|IH|ON|QP|PO',baseline_strata$locations$stationid),]
        baseline_strata$locations$corner<-ifelse(grepl('GF|HG|JI|IH|ON|QP|PO',baseline_strata$locations$stationid),'TRUE','FALSE')
        
        #get strata and cells from the current allocation sampling design (stratum current design - strata optimized stratification)
        current<-merge(baseline_strata$locations,strata,by='cell',all.x=TRUE)

        #create raster polygon of stratified design
        #df to spatialpoint df
        D11<-D8
        coordinates(D11) <- ~ Lon + Lat
        crs(D11)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
        
        #reproject coordinates for plotting purposes
        df_1<-spTransform(D11,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
        df_2<-data.frame(df_1)
        
        #x and y cells
        xycells<-as.integer(sqrt(dim(df_1)[1]))
        
        # create a template raster
        r1 <- raster(ext=extent(df_1),ncol=xycells, nrow=xycells) #c(15800,15800) 7000
        
        #create raster
        r2<-rasterize(df_1, r1 ,field=c('strata'))
        crs(r2) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
        r2[r2==999]<-NA
        
        #create polygon to get boundaries of each strata
        r3<-as.data.frame(r2,xy=TRUE)
        r4<-rasterToPolygons(r2$layer,dissolve=TRUE)
        r5<-spTransform(r4,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
        # plot(r4)
        # plot(r5)
        # r6<-spTransform(r5,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
        # 
        # plot(r2)
        # plot(r6,add=TRUE)
        # 
        #replicates n_surveys x n_years
        rep_sur<-n_sur*length(yrs)
        
        #create array to store results
        #to store survey samples
        scn_allocations <- array(data = NA, dim = c(as.numeric(n_samples)*rep_sur,
                                                    length(c('Lat','Lon','cell','strata','sur')),
                                                    #length(yrs), 
                                                    length(c('systematic','random','sb'))),
                                 dimnames = list(1:(as.numeric(n_samples)*rep_sur),
                                                 c('Lat','Lon','cell','strata','sur'),
                                                 #c(yrs),
                                                 c('sys','rand','sb')))
      
        

        
        #for (isur in 1:n_sur) {
        #  
        #  isur<-1
        #  
        #  for (y in 1:length(yrs)) { #years
        #    
        #    y<-1
        #    
        #    #year
        #    yy<-yrs[y]
        #    
        #    #to store points
            dfcurrent<-data.frame(matrix(NA,nrow=0,ncol=5))
            colnames(dfcurrent)<-c('Lon','Lat','cell','strata','sur')
            dfspb<-dfrandom<-dfcurrent
        #    
        #    #for while purposes
        #    flag<-TRUE
            
            str_alloc<-list()
            
          
            for (str in allocations$strata) {
              
              #str<-allocations$strata[1]
              
              #print simulation to check progress
              cat(paste(" #############  sampling design",samp ,'- strata',str ," #############\n"))
              
              #number of samples required by str
              n_i<-allocations[which(allocations$strata==str),'n_samples']
              
              #subset cells grid by strata
              D9<-subset(D8,strata==str)
              
              ########################
              ### RANDOM SAMPLING
              ########################
              
              # Initializing parallel backend
              cl <- makeCluster(detectCores()-1)
              registerDoParallel(cl)
              

              # Parallelizing the loop
              dfrandom <- foreach(sur = 1:rep_sur, .combine=rbind) %dopar% {
                #cat(paste(" #############  random sampling", '- survey', sur, " #############\n"))
                
                D10 <- D9[sample(nrow(D9), size = n_i, replace = FALSE), ]
                D10 <- D10[, c('Lon', 'Lat', 'cell', 'strata')]
                D10$sur <- sur
                
                return(D10)
              }
              
              # Stopping the parallel backend
              stopCluster(cl)
              
              ########################
              ### SPATIALLY BALANCE SAMPLING
              ########################
              
              #get distance among points
              dis_la <- as.matrix(dist(cbind(D9$Lon, D9$Lat)))
              D9_st_sf <- st_as_sf(D9, coords = c("Lon", "Lat"))
              dis_la_sf <- st_distance(D9_st_sf)
              
              #PWD
              con <- rep(0, nrow(dis_la))
              stand_dist_la_pwd <- Spbsampling::stprod(mat = dis_la, con = con)$mat

              # Initializing parallel backend
              cl <- makeCluster(detectCores()-1)
              registerDoParallel(cl)
              
              # Create an empty list to store the results
              dfspb <- foreach(sur = 1:rep_sur, .combine = rbind) %dopar% {
                
                #cat(paste(" #############  spatially-balanced sampling", '- survey', sur, " #############\n"))
                
                s_pwd_la <- Spbsampling::pwd(dis = stand_dist_la_pwd, n = n_i)$s
                D10 <- D9[s_pwd_la[1, ], ]
                
                D10$sur <- sur
                # spatially-balancing strata samples
                spbsamp <- D10[, c('Lon', 'Lat', 'cell', 'strata', 'sur')]
                
                return(spbsamp)
              }

              # Stopping the parallel backend
              stopCluster(cl)
              
              #######################
              ### SYSTEMATIC SAMPLING
              #######################
              
              # Initializing parallel backend
              cl <- makeCluster(detectCores()-1)
              registerDoParallel(cl)
              
              # Create an empty list to store the results
              dfcurrent <- foreach(sur = 1:rep_sur, .combine = rbind) %dopar% {
                
                #cat(paste(" #############  systematic sampling", '- survey', sur, " #############\n"))
                
                #if current sampling design
                if (samp_df[s, 'samp_scn'] == 'scnbase') {
                  pointsc <- current[which(current$stratum == str), c('longitude', 'latitude', 'cell', 'stratum')]
                  names(pointsc) <- c('Lon', 'Lat', 'cell', 'strata')
                  
                  #if current sampling design w/o corner
                } else if (samp_df[s, 'samp_scn'] == 'scnbase_bis') {
                  pointsc <- current[which(current$corner == FALSE & current$stratum == str), c('longitude', 'latitude', 'cell', 'stratum')]
                  names(pointsc) <- c('Lon', 'Lat', 'cell', 'strata')
                  
                  #if optimized sampling design
                } else {
                  
                  sys <- current[which(current$strata == str), c('longitude', 'latitude', 'cell', 'strata')]
                  names(sys) <- c('Lon', 'Lat', 'cell', 'strata')
                  
                  #if more required samples than available
                  if (nrow(sys) < n_i) {
                    
                    #duplicate df strata
                    dff <- D9
                    
                    #df removing available samples from current design
                    dff <- subset(dff, !(cell %in% sys$cell))
                    
                    #cells to complete the required cells
                    ii <- n_i - nrow(sys)
                    
                    #loop over the required samples using random
                    selcell <- sample(dff$cell, ii)
                    
                    if (nrow(subset(D9, cell %in% c(selcell, sys$cell))) == 0) {
                      
                      ss <- sample(1:nrow(D9), size = n_i, replace = FALSE)
                      pointsc <- D9[ss, c('Lon', 'Lat', 'cell', 'strata')]
                      
                      #subset points if equal to required number of samples, if not ERROR and rerun iteration
                    } else if (nrow(subset(D9, cell %in% c(selcell, sys$cell))) != n_i) {
                      flag <- TRUE
                      stop("-- different number of points on the current approach", call. = FALSE)
                    } else {
                      pointsc <- subset(D9, cell %in% c(selcell, sys$cell))[, c('Lon', 'Lat', 'cell', 'strata')]
                      names(pointsc) <- c('Lon', 'Lat', 'cell', 'strata')
                    }
                    
                    #else there are enough samples to get from the current sampling design    
                  } else {
                    ss <- sample(1:nrow(sys), size = n_i, replace = FALSE)
                    pointsc <- sys[ss, c('Lon', 'Lat', 'cell', 'strata')]
                  }
                }
                
                #systematic strata samples
                syssamp <- pointsc[, c('Lon', 'Lat', 'cell', 'strata')]
                syssamp$sur <- sur
                
                return(syssamp)
              }
              
              stopCluster(cl)
              
              #######################
              ### PLOT
              #######################
              
              #coordinates(rsamp)=~Lon + Lat
              #crs(rsamp)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
              #rsamp<-spTransform(rsamp,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
              #
              #coordinates(spbsamp)=~Lon + Lat
              #crs(spbsamp)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
              #spbsamp<-spTransform(spbsamp,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
              #
              #coordinates(syssamp)=~Lon + Lat
              #crs(syssamp)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
              #syssamp<-spTransform(syssamp,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
              #
              #pol<-r4[r4$layer==str,]
              #plot(pol)
              #points(rsamp,col='green',pch=19)
              #points(spbsamp,col='blue',pch=19)
              #points(syssamp,col='red',pch=19)
              
              #store into a list
              str_alloc[[str]]<-list(dfrandom = dfrandom, dfspb = dfspb, dfcurrent = dfcurrent)
              
            }
            
            #reshape
            cur<-dplyr::bind_rows(lapply(str_alloc, function(x) x[['dfcurrent']]))
            rand<-dplyr::bind_rows(lapply(str_alloc, function(x) x[['dfrandom']]))
            spb<-dplyr::bind_rows(lapply(str_alloc, function(x) x[['dfspb']]))

          #append results 
          scn_allocations[,,'sys'] <- unlist(cur)
          scn_allocations[,,'rand'] <- unlist(rand)
          scn_allocations[,,'sb'] <- unlist(spb)
        
          
          #remove
          rm(dfcurrent,dfrandom,dfspb,str_alloc,cur,rand,spb)
          
          #data.frame(table(dfcurrent[,c("strata",'sur')]))
          #data.frame(table(dfcurrent[,c('sur')]))    
              
          #dfcurrent<-dfcurrent[order(dfcurrent$sur),]  #520 100        
          #nrow(dfcurrent)==as.numeric(n_samples)*n_sur*length(yrs)    
          #dfcurrent$year<-rep(yrs,each=as.numeric(n_samples)*n_sur)

          save(scn_allocations, file = paste0('./output/survey_allocations_',samp,'.RData')) 
          
    }
      
      
################
# get simulated densities from simulated survey HISTORICAL
################
 
#samp<-samp_df$samp_scn[1];ap<-'sb';isim<-1;isur<-1;y<-'1986'
n_sur=100
sur_df<-
data.frame(num=1:(length(yrs)*n_sur),
           year=rep(yrs,times=n_sur),
           sur=rep(1:n_sur,each=length(yrs)))

#load simulated densities by spp
#load(paste0('./output/species/sim_hist_dens_spp.RData')) #sim_hist_dens_spp

# library(foreach)
# library(doParallel)

# # Initializing parallel backend
# cl <- makeCluster(detectCores()-1)  # Using all available cores
# registerDoParallel(cl)

#simulated densities
load(file = paste0('./output/species/ms_sim_dens.RData'))  #sim_dens1

#ms_sim_survey folder
dir.create('./output/ms_sim_survey/')

# Create folders with numbers from 1 to 100 and leading zeros from simulated data 
for (i in 1:100) {
  folder_name <- paste0("./output/ms_sim_survey/", formatC(i, width = 3, flag = "0"))
  dir.create(folder_name)
}

# Check the created folders
list_of_folders <- list.dirs("./output/ms_sim_survey/", recursive = FALSE)
print(list_of_folders)

# Parallelizing the loop
for (samp in samp_df$samp_scn)  {
  
  samp<-'scn1'
  #start_time_parallel <- Sys.time()
  
  #array to store simulated densities/CPUE
  alloc<-ifelse(samp=='scnbase_bis',494,520)
  
  #load survey allocations by sampling design
  load(file = paste0('./output/survey_allocations_',samp,'.RData')) #scn_allocations
  dimnames(scn_allocations)[[3]]<-c('sys','rand','sb')
  
  for (isim in 2:n_sim_hist) {

    fol<- list_of_folders[isim]
    
    sim_survey <- array(NA,
                        dim = c(alloc, length(spp)+2, length(unique(yrs)), n_sur,length(c('sys','rand','sb'))),
                        dimnames = list(1:alloc, c('cell','strata',spp), unique(yrs), 1:n_sur,c('sys','rand','sb')))
    
    
    for (y in as.character(yrs)) {
      
      #isim<-1;y<-as.character(yrs[1])
      
      cat(paste(" #############  ",samp,'- simdata',isim, '- year',y," #############\n"))
      
      
      sim_dens2<-sim_dens1[,,y,isim]
      
      surs<-sample(1:max(sur_df$num),size = n_sur)  
      
      for (sur in as.character(surs)) {
        
        #sur<-as.character(surs)[1]

        cat(paste(" #############  ",samp,'- simdata',isim, '- year',y,'- simsur',sur ," #############\n"))
      
        #systematic
        sim_survey[,,y, match(sur,surs),'sys'] <- 
          cbind(scn_allocations[scn_allocations[,'sur','sys']==sur,c('cell','strata'),'sys'],
                dens=sim_dens2[scn_allocations[scn_allocations[,'sur','sys']==sur,c('cell'),'sys'],])
        
        #random
        sim_survey[,,y, match(sur,surs),'rand'] <- 
          cbind(scn_allocations[scn_allocations[,'sur','rand']==sur,c('cell','strata'),'rand'],
                dens=sim_dens2[scn_allocations[scn_allocations[,'sur','rand']==sur,c('cell'),'rand'],])
        
        #sb
        sim_survey[,,y, match(sur,surs),'sb']  <- 
          cbind(scn_allocations[scn_allocations[,'sur','sb']==sur,c('cell','strata'),'sb'],
                dens=sim_dens2[scn_allocations[scn_allocations[,'sur','sb']==sur,c('cell'),'sb'],])
        
      }
      
    }
    
    save(sim_survey, file = paste0(fol,'/sim_survey_',samp,'.RData'))  
    
  }

  
}

# Stopping the parallel backend
stopCluster(cl)




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
   


################
# get simulated densities from simulated survey PROJECTED
################

#all 4000 surveys (year5, sur100, sbt8)
sur_df<-cbind(num=1:nrow(expand.grid(project_yrs,df_sbt$sbt,1:n_sur)),expand.grid(year=project_yrs,sbt=df_sbt$sbt_n,sur=1:n_sur))

#load sbt scenarios
load('./tables/SBT_projection.RData') #df_sbt

simdata<-array(NA,
               dim = c(nrow(grid), length(spp), length(unique(2023:2027)), n_sim_proj,length(df_sbt$sbt)),
               dimnames = list(1:nrow(grid),spp, as.character(2023:2027), 1:n_sim_proj,df_sbt$sbt))

#loop over species
for (sp in spp) {
  
  sp<-spp[1]
  
  #loop over 8 temperature scenarios
  for (sbt in unique(df_sbt$sbt_n)) {
    
    #sbt<-df_sbt$sbt_n[1]
    
    cat(paste(" #############  spp",sp,'- sbt',sbt ," #############\n"))
    
    
    load(file = paste0("./output/species/",sp,"/simulated projected data/SBT",sbt," dens_index_proj_OM_50.RData")) #dens_index_proj_OM
    dens_index_proj_OM_50<-dens_index_proj_OM
    load(file = paste0("./output/species/",sp,"/simulated projected data/SBT",sbt," dens_index_proj_OM_100.RData")) #dens_index_proj_OM
    dens_index_proj_OM_100<-dens_index_proj_OM
    rm(dens_index_proj_OM)
    
    for (i in 1:50) {
      
      #i<-1
      d1<-dens_index_proj_OM_50[[i]]$dens[,as.character(2023:2027)]
      d2<-dens_index_proj_OM_100[[i]]$dens[,as.character(2023:2027)]
      
      for (y in as.character(2023:2027)) {
        
        #y<-as.character(2023:2027)[1]
        
        simdata[,sp,y,i,sbt]<-d1[,y]
        simdata[,sp,y,i+50,sbt]<-d2[,y]
        
      }
    }
  }
  
  
  
  
  
  #for each sampling design
  for (samp in samp_df$samp_scn) {
    
    #samp<-samp_df$samp_scn[1]
    
    #load survey allocations by sampling design
    load(file = paste0('./output/species/ms_sim_survey/survey_allocations_',samp,'.RData')) #scn_allocations
    
    #to store survey samples
    sim_survey <- array(data = NA, dim = c(dim(scn_allocations)[1]/n_sur/length(yrs),
                                           length(c("Lon","Lat","cell","strata","sur"))+1,
                                           n_sim_proj,
                                           length(c('sys','rand','sb')),
                                           n_sur,
                                           length(project_yrs)),
                        dimnames = list(c(1:(dim(scn_allocations)[1]/n_sur/length(yrs))),
                                        c("Lon","Lat","cell","strata","sur",'CPUE_kgkm'),
                                        1:n_sim_proj,
                                        c('sys','rand','sb'),
                                        1:n_sur,
                                        as.character(project_yrs)))
    
    #loop over each of 100 simulated projections
    for (isim in 1:n_sim_proj) {
      
      cat(paste(" #############  simulated data", isim," #############\n"))
      
      #isim<-1
      
      if (isim %in% 1:50) {
        dens_index_proj_OM1<-data.frame(drop_units(dens_index_proj_OM_50[[paste0('sim',isim)]]$dens),check.names = FALSE)
      } else if (isim %in% 51:100) {
        dens_index_proj_OM1<-data.frame(drop_units(dens_index_proj_OM_100[[paste0('sim',isim-50)]]$dens),check.names = FALSE)
      }
      
      # #for 100 simulated surveys in each year
      # for (y in project_yrs) {
      #  
      #   y<-2023
      #   dd<-dens_index_proj_OM1[,as.character(y)]
      
      dim(dens_index_proj_OM1)
      dd<-dens_index_proj_OM1[,as.character(project_yrs)]
      
      #for each allocation sampling approach
      for (ap in c('sys','rand','sb')) {
        
        #ap<-'systematic'
        
        #get allocations by sampling allocation approach
        scn_allocations1<-data.frame(scn_allocations[,,ap])
        
        #loop over 100 surveys
        for (isur in 1:n_sur) {
          
          #get survey number 1:4000 (5y * 8 sbt * 100sur)
          sur_num<-sur_df[which(sur_df$sbt==sbt & sur_df$sur==isur),] #'num'
          
          for (y in as.character(project_yrs)) {
            
            #y<-project_yrs[1]
            
            scn_allocations2<-scn_allocations1[which(scn_allocations1$sur == sur_num[which(sur_num$year==y),'num']),]
            dd1<-dd[scn_allocations2$cell,as.character(y)]
            
            sim_survey[,,isim,ap,isur,y]<-cbind(scn_allocations2,CPUE_kgkm=dd1)
            
          }
        }
      }
    }
    #store results
    save(sim_survey, file = paste0('./output/species/',sp,'/simulated survey project/sim_proj_survey_sbt',sbt,'_',samp,'.RData'))  
    #save(all_points, file = paste0('./output/species/allocations_hist_survey_',samp,'.RData'))    
    
    rm(sim_survey)
  }
}
}                 
                  