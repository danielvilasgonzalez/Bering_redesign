####################################################################
####################################################################
##
##    Run sampling optimization based on predicted densities from VAST OM EBS+NBS
##    and calculate stratification boundaries and sample allocations for each sampling design
##
##    by best stratification, we mean the stratification that ensures the minimum sample cost, 
##    sufficient to satisfy precision constraints set on the accuracy of the estimates of the survey target variables Y’s
##    constraints expressed as maximum allowable coefficients of variation in the different domains of interest
##
##    https://cran.r-project.org/web/packages/SamplingStrata/vignettes/SamplingStrata.html
##    (using devtools::install_github(repo = "zoyafuso-NOAA/SamplingStrata"))
##    Daniel Vilas (daniel.vilas@noaa.gov/dvilasg@uw.edu)
##    Lewis Barnett
##    Zack Oyafuso (some code borrowed from his previous analysis in the GOA)
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c("splines",'SamplingStrata','wesanderson','dplyr','sp',
             'sf','maptools','rgeos','scales','raster',
             'rnaturalearth','grid','ggplot2','spatstat','ragg')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#install akgfmaps to extract shapefile of Alaska
if (!('akgfmaps' %in% installed.packages())) {
  devtools::install_github('afsc-gap-products/akgfmaps')};library(akgfmaps)

#install coldpool to extract SBT for the EBS
if (!('VAST' %in% installed.packages())) {
  devtools::install_github("james-thorson/VAST@main", INSTALL_opts="--no-staged-install")};library(VAST)

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#parallel
#cl <- parallel::makeCluster(parallel::detectCores() - 2)
#doParallel::registerDoParallel(cl)

#setwd - depends on computer using
#out_dir<-'C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/' #NOAA laptop  
out_dir<-'/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/' #mac
#out_dir<-'/Users/daniel/Work/VM' #VM
setwd(out_dir)

#version VAST (cpp)
version<-'VAST_v14_0_1'

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
       'Chionoecetes bairdi',
       'Sebastes alutus',
       #'Sebastes melanostictus',
       'Atheresthes evermanni')

#remove Anoploma and Reinhardtius because habitat preference reasons
#spp<-setdiff(spp, c('Anoplopoma fimbria','Reinhardtius hippoglossoides'))

spp1<-c('Yellowfin sole',
        'Alaska pollock',
        'Pacific cod',
        'Arrowtooth flounder',
        'Greenland turbot',
        'Northern rock sole',
        'Flathead sole',
        'Alaska plaice',
        'Bering flounder',
        'Arctic cod',
        'Saffon cod',
        'Sablefish',
        'Snow crab',
        'Blue king crab',
        'Red king crab',
        'Tanner crab',
        'Pacific ocean perch',
        'Kamchatka flounder')

spp_name<-data.frame('spp'=spp,
                     'common'=spp1) 


#df spp, number and target variables
df_spp<-data.frame('spp'=spp,
                   'n'=c(1:length(spp)),
                   'Y'=paste0('Y',c(1:length(spp))))

#number sp
n_spp<-length(spp)

###################################
# LOAD GRID EBS (remember to keep the same order as in fit_model if multiple grids)
###################################

#load slope grid
load('./extrapolation grids/bering_sea_slope_grid.rda')
dim(bering_sea_slope_grid)
names(bering_sea_slope_grid)[4]<-'Stratum'
bering_sea_slope_grid$Stratum<-999
#gridslope<-data.frame(bering_sea_slope_grid,region='SLP')

#load EBS+NBS grid
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),
                          data.frame(eastern_bering_sea_grid,region='EBS'),
                          data.frame(bering_sea_slope_grid,region='SLP')))
ncell_ebsnbs<-nrow(rbind(data.frame(northern_bering_sea_grid,region='NBS'),
                         data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)
grid$cell<-as.numeric(grid$cell)

#add col and row number
x1<-grid[,c('Lon','Lat','cell')]
names(x1)<-c('x','y','z')
coordinates(x1)=~x + y
#crs(x1)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
proj4string(x1)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
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
x5<-x5[,c('y','x','z','optional','col','row'),]
colnames(x5)<-c('Lat','Lon','cell','optional','col','row')
grid<-x5[,c('Lat','Lon','cell','col','row')]

# Define plot extent (through trial end error) units km (for plotting purposes)
panel_extent <- data.frame(x = c(-1716559.21, -77636.05), #x = c(-1326559.21, -87636.05),
                           y = c(483099.5, 2194909.7)) #y = c(533099.5, 1894909.7))

#####################################
# Polygon regions shapefiles (EBS, NBS and slope)
#####################################

#name shapefiles 
shfiles<-c('EBSshelfThorson','NBSThorson','EBSslopeThorson')

#loop over shapefiles
for (i in shfiles) {
  
  #shapefile EBS
  sh<-rgdal::readOGR(dsn='./shapefiles/',layer = i)
  
  if (i=='EBSslopeThorson') {
    
    #reproject shapefile
    proj4string(sh) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
    sh<-spTransform(sh,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))
    
  }
  
  #shapefile name
  shname<-paste0(gsub('Thorson','',i),'_sh')
  
  #assign shapefiles
  assign(shname,sh)
  
}

#merge shapefiles
bs_sh1<-raster::union(EBSshelf_sh,NBS_sh)
bs_sh<-raster::union(bs_sh1,EBSslope_sh)


###################################
# FXNs - extracted from ZO code
###################################

calc_expected_CV <- function (strata) {
  
  n_spp <- length(grep(x = names(strata), pattern = "M"))
  n_h <- strata$Allocation
  N_h <- strata$Population
  
  cv <- vector(length = n_spp)
  names(cv) <- paste0("Y", 1:n_spp)
  
  for (ispp in 1:n_spp) {
    S_h <- strata[, paste0("S", ispp)] #SD from 
    M_h <- strata[, paste0("M", ispp)] #mean from 
    
    Y_h <- N_h * M_h
    Var_h <- (N_h^2) * (1 - n_h/N_h) * ((S_h^2)/n_h)
    CV <- sqrt(sum(Var_h))/sum(Y_h)
    
    cv[ispp] <- CV
  }
  
  cv <- round(cv, 3)
  return(cv)
}

###################################
# Sampling designs
###################################

#sampling scenarios
samp_df<-expand.grid(strat_var=c('varTemp','Depth','varTemp_forced','Depth_forced'), #LonE and combinations
                     target_var=c('sumDensity'), #,'sqsumDensity'
                     n_samples=c(520), #c(300,500) 520 (EBS+NBS+CRAB);26 (CRAB); 350 (EBS-CRAB); 494 (NBS-CRAB)
                     n_strata=c(10),
                     domain=1) #c(5,10,15)

#add other columns
samp_df$idomain<-NA
samp_df1<-samp_df
samp_df1$n_strata<-5
samp_df1$domain<-2
samp_df1$idomain<-'region'

#scenario slope forced
# slope_for<-data.frame(strat_var=c('varTemp_forced','Depth_forced'),
#                        target_var='sumDensity',
#                        n_samples=520,
#                        n_strata=5,
#                        domain=2,
#                        idomain='region',stringsAsFactors = FALSE)


# slope_samp<-data.frame(strat_var=c('Lat','Depth'),
#                        target_var='sumDensity',
#                        n_samples=520,
#                        n_strata=5,
#                        domain=1,
#                        idomain='slope',stringsAsFactors = FALSE)

#rbind scenarios
samp_df<-rbind(samp_df,samp_df1)

#add scenario number
samp_df$samp_scn<-paste0(paste0('scn',1:nrow(samp_df)))


###############
# load ms data and settings
###############

  #load multispecies data
  load(paste0('./output/multisp_optimization_static_data_ebsnbs_slope.RData')) #df
  
  tar_var<-paste0(rep(df_spp$Y,each=2),c('','_SQ_SUM'))
  names(df)[((ncol(df)-length(tar_var))+1):ncol(df)]<-tar_var
  ispp<-n_spp

  #removed cells because of depth
  rem_cells<-df[which(df$include==FALSE),'cell']
  ok_cells<-df[which(df$include==1),'cell']
  
  #load data_geostat file
  #data_geostat<-readRDS(paste0('./data processed/species/',sp,'/','data_geostat_temp.rds')) 
  
  #years
  n_years<-length(1982:2022)
  
  #subset cells with appropiate depth
  static_df1<-subset(df,cell %in% ok_cells & Depth <= 400)
  
  #n cells
  n_cells<-nrow(static_df1)
  
  #domain_input
  domain_input<-rep(1, nrow(static_df1))
  domain_grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),
                                   data.frame(eastern_bering_sea_grid,region='EBS'),
                                   data.frame(bering_sea_slope_grid,region='SLP')))
  #cell
  domain_grid$cell<-1:nrow(domain_grid)
  domain_grid$cell<-as.numeric(domain_grid$cell)
  #domain based on the NBS-EBS
  domain_grid$domain_region<-ifelse(domain_grid$region=='NBS',1,2)
  domain_grid<-domain_grid[which(domain_grid$cell %in% static_df1$cell),]
  dim(domain_grid)
  
  #add forced attribute to force optimize separetely the slope
  slope_cells<-domain_grid[which(domain_grid$region=='SLP'),'cell']
  #nbs_cells<-domain_grid[which(domain_grid$region=='NBS'),'cell']
  static_df1$forced<-1
  #static_df1[which(static_df1$cell %in% slope_cells),'forced']<-99999
  static_df1[which(static_df1$cell %in% slope_cells),'forced']<-2
  #static_df1[which(static_df1$cell %in% nbs_cells),'forced']<-9
  summary(static_df1)
  static_all<-static_df1
  
  aggregate(cell ~ forced,static_all,FUN=length)
  
  #########################
  # loop over optimized sampling designs
  #########################
  
  #s<-5 cannot find a solution
  
  #loop through sampling designs
  for (s in c(1:nrow(samp_df))) { #nrow(samp_df)

    #s<- 6
      
    #domain
    dom<-samp_df[s,'domain']
    idom<-samp_df[s,'idomain']
    
    #static_all[which(static_all$cell>=ncell_ebsnbs+1),'Depth']<-1000
    
    if (is.na(idom)) {
      static_df1<-static_all
      tar_var<-paste0(rep(df_spp$Y,each=2),c('','_SQ_SUM'))
      #names(df)[((ncol(df)-length(tar_var))+1):ncol(df)]<-tar_var
      ispp<-n_spp<-18
    
    } else if(idom=='slope') {
      static_df1<-static_all[which(static_all$cell>=ncell_ebsnbs+1),]
      static_df1<-static_df1[,colSums(static_df1[,1:ncol(static_df1)]) > 0]
      colnames(static_df1)[9:(ncol(static_df1)-1)]<-
        sort(c(paste0('Y',1:(length(colnames(static_df1)[9:(ncol(static_df1)-1)])/2)),
          paste0('Y',1:(length(colnames(static_df1)[9:(ncol(static_df1)-1)])/2),c('_SQ_SUM'))))
      tar_var<-colnames(static_df1)[9:(ncol(static_df1)-1)]
      #names(df)[((ncol(df)-length(tar_var))+1):ncol(df)]<-tar_var
      ispp<-n_spp<-length(tar_var)/2
          
    } else {
      static_df1<-static_all
      tar_var<-paste0(rep(df_spp$Y,each=2),c('','_SQ_SUM'))
      #names(df)[((ncol(df)-length(tar_var))+1):ncol(df)]<-tar_var
      ispp<-n_spp<-18
    }
    
    
    #domain input
    if (dom==1) {
      domain_input<-rep(1, nrow(static_df1))
    } else if (idom=='region') {
      domain_input<-domain_grid$domain_region
    } else if (idom=='lat'){
      domain_input<-domain_grid$domain_lat
    }
    
    
    
    
    
    #print scenario to check progress
    cat(paste("\n #############  Sampling Scenario", samp_df[s,"samp_scn"], " #############\n"))
    
    #if scenario includes two stratifying factors
    if (grepl('_',samp_df[s,'strat_var'])) {
      #stratification variables 
      stratum_var_input<-data.frame(X1 = static_df1[,paste0(sub("\\_.*", "", samp_df[s,'strat_var']))],
                                    X2 = static_df1[,paste0(sub(".*_", "", samp_df[s,'strat_var']))]) 
    } else {
      stratum_var_input<-data.frame(X1 = static_df1[,paste0(sub("\\_.*", "", samp_df[s,'strat_var']))])
    }

    #test
    #stratum_var_input$X2<-rlnorm(length(stratum_var_input$X2),meanlog = 2,sdlog = 1)
    
    #target variables
    target_var_input<-static_df1[,tar_var]
    n_samples<-srs_n<-samp_df[s,'n_samples']
    
    #create df
    frame <- data.frame(domainvalue = domain_input, #domain
                        id = static_df1$cell, #id as cells
                        stratum_var_input, #Stratification variables
                        WEIGHT=n_years, #weight for spp depending on the presence of years
                        target_var_input) #target variables 
    
    #correct years
    # frame$WEIGHT<-ifelse(frame$id %in% slope_cells,6,
    #                      ifelse(frame$id %in% nbs_cells,5, 40))
    
    ###################################
    # Simple random sampling design CV constraints
    ###################################
    
    #Initiate CVs to be those calculated under simple random sampling (SRS)    
    #build strata frame (if 2 factors stratification or 1 factors) 
    if (grepl('_',samp_df[s,'strat_var'])) {
      #stratification variables 
      srs_stats <- buildStrataDF(dataset = cbind( frame[, -grep(x = names(frame), pattern = "X")],X1 = 1,X2=1))
    } else {
      srs_stats <- buildStrataDF(dataset = cbind( frame[, -grep(x = names(frame), pattern = "X")],X1 = 1))
    }
    
    ## SRS statistics
    srs_var <- srs_stats[, paste0("S", 1:ispp)]^2 * (1 -  srs_n/ n_cells) / srs_n
    srs_cv <- sqrt(srs_var) / srs_stats[, paste0("M", 1:ispp)]
    
    #CV
    cv_df <- list()
    cv_df[["DOM"]] <-c(1:dom)
    for (iispp in 1:ispp) cv_df[[paste0("CV", iispp)]] <- as.numeric(srs_cv[,iispp])
    cv_df[["domainvalue"]] <-c(1:dom)
    cv_df <- as.data.frame(cv_df)
    cv_df[is.na(cv_df)]<-1
    
    
    #store Cv
    cv_initial<-cv_df
    
    #number of samples ###520 in 2022 (NBS and EBSshelf) #376 in 2018 (EBSshelf)
    n_samples <- samp_df[s,'n_samples']
    
    #number of stratas 
    no_strata<-samp_df[s,'n_strata']
    
    #number of samples
    srs_n <- as.numeric(n_samples * table(frame$domainvalue) / n_cells)
    
    ###################################
    # Run optimization
    ###################################
    
    #for plot 
    load('./data processed/grid_EBS_NBS.RData')
    gridi<-
      unique(grid.ebs_year[,c('Lat','Lon','DepthGEBCO')])
    gridi1<-subset(gridi,DepthGEBCO <= 400)
    
    #count
    count<-0
    
    #while loop until the desired number of strata  
    flag<-TRUE
    while (flag) {
      
      count<-count+1
      
      #print state of spp[1] to get an approximation of the state
      cat(paste("\n #############   OPTIMIZING STRATA - CV of ",spp[1],' ', cv_df[,2],"  #############\n"))
      cat(paste("\n #############   COUNT", count ,"  #############\n"))
      
      #run optimization
      solution <- optimStrata(method = "continuous", #continous variables
                              errors = cv_df,  #precision level - maximum allowable coefficient of variation set by the simple random sampling 
                              framesamp = frame, #df of input variables
                              iter = 150, #30, #300 #maximum number of iterations
                              pops = 50,#10, #100  #dimension of each generations
                              elitism_rate = 0.1, #0.1
                              mut_chance = 1 / (no_strata[1] + 1), #mutation chance
                              nStrata = c(no_strata,no_strata), #maximum strata
                              showPlot = TRUE, #FALSE
                              writeFiles = FALSE)
      
      #plot strata
      solution$framenew$strata<-paste0(solution$framenew$DOMAINVALUE,'_',solution$framenew$STRATO)
      print(
      ggplot()+
        geom_point(data=cbind(static_df1[,c('Lat','Lon')],solution$framenew),aes(x=Lon,y=Lat,color=strata))
      )
      
      #solution$aggr_strata
      #flag to keep the loop
      flag<-ifelse(nrow(solution$aggr_strata)!=samp_df$n_strata[s]*dom,TRUE,FALSE)
      
       #if condition reduce or increase CV to achieve the objective
       if (nrow(solution$aggr_strata)<samp_df$n_strata[s]*dom) {
         cv_df[,c(2:(ispp+1))]<-cv_df[,c(2:(ispp+1))]-(cv_df[,c(2:(ispp+1))]*0.001) #reduce 0.1% CV ###0.001
       } else if (nrow(solution$aggr_strata)>samp_df$n_strata[s]*dom){
        cv_df[,c(2:(ispp+1))]<-cv_df[,c(2:(ispp+1))]+cv_df[,c(2:(ispp+1))]*0.10} #increase 1%
    }




    
    ###################################
    # Store solutions from optimizations
    ###################################
    
    #store strata CV
    cv_strata_final<-cv_df
    
    solution$aggr_strata$STRATO <- as.integer(solution$aggr_strata$STRATO)
    solution$aggr_strata <- 
      solution$aggr_strata[order(solution$aggr_strata$DOM1,
                                 solution$aggr_strata$STRATO), ]
    
    
    # plotStrata2d(solution$framenew,
    #               solution$aggr_strata,
    #              domain = 1,
    #              vars=c('X1','X2'),
    #              labels = c('Depth','Lat'))
    
    
    sum_stats <- SamplingStrata::summaryStrata(solution$framenew,
                                               solution$aggr_strata,
                                               progress=FALSE)
    sum_stats$stratum_id <- 1:nrow(sum_stats)
    sum_stats$wh <- sum_stats$Allocation / sum_stats$Population
    sum_stats$Wh <- sum_stats$Population / nrow(frame)
    sum_stats <- cbind(sum_stats,
                       subset(x = solution$aggr_strata,
                              select = -c(STRATO, N, COST, CENS, DOM1, X1)))
    
    plot_solution <- solution$indices
    
    #save optimization stats  
    sum_stats <- summaryStrata(solution$framenew,
                               solution$aggr_strata,
                               progress=FALSE)
    sum_stats$stratum_id <- 1:nrow(sum_stats)
    sum_stats$Population <- sum_stats$Population / n_years
    sum_stats$wh <- sum_stats$Allocation / sum_stats$Population
    sum_stats$Wh <- sum_stats$Population / n_cells
    sum_stats <- cbind(sum_stats,
                       subset(x = solution$aggr_strata,
                              select = -c(STRATO, N, COST, CENS, DOM1, X1)))
    
    #add scn and sp
    sum_stats$samp_scn<-samp_df[s,'samp_scn']
    
    #if one stratifying factor add columns  
    if (!grepl('_',samp_df[s,'strat_var'])) {
      
      sum_stats<-data.frame(sum_stats[,c("Domain","Stratum","Population","Allocation","SamplingRate","Lower_X1","Upper_X1")],
                            "Lower_X2"=NA,"Upper_X2"=NA,
                            sum_stats[,c("stratum_id","wh","Wh","M1",'S1',"SOLUZ","samp_scn")])
    }

    #store results
    result_list <- list(solution = solution,
                        sum_stats = sum_stats,
                        cvs = as.numeric(calc_expected_CV(sum_stats)),
                        n = sum(sum_stats$Allocation),
                        sol_by_cell = plot_solution)
    
    #######################
    ##   7) Single-Species Optimization ----
    ##   Calculate single-species CV subject to the initial stratification
    #######################

      ss_sample_allocations <- expand.grid(n = n_samples, species = spp)
    
    for (iispp in 1:ispp) {
      #iispp<-17
      
      temp_n <- result_list$n
      
      ## Subset density data for species ispp
      ss_df <- subset(x = frame, 
                      select = c("domainvalue", "id", "WEIGHT", "X1", #"X2", 
                                 paste0("Y", iispp), paste0("Y", iispp, "_SQ_SUM")))
      names(ss_df)[grep(x = names(ss_df), pattern = "Y")] <- c("Y1", "Y1_SQ_SUM")
      
      ## Create CV inputs to the Bethel algorithm; initialize at SRS CV
        error_df <- data.frame("DOM" = c(1:dom),
                              'CV1'=as.numeric(srs_cv[,iispp]),
                               "domainvalue"  = c(1:dom))
        #names(error_df)[2] <- "CV1"
        
        ## subset stratum stats for the species of interest as inputs to the 
        ## Bethel algorithm
        temp_stratif <- 
            solution$aggr_strata[, c("STRATO", "N", 
                                     paste0("M", iispp), paste0("S", iispp), 
                                     "COST", "CENS", "DOM1", "X1" ,"SOLUZ")]
        
        if (dom==2) {
          temp_stratif$DOM2<-temp_stratif$DOM1
        }
        
        
        temp_stratif$N <- temp_stratif$N / n_years
        #temp_stratif$DOM1 <- 1
        names(temp_stratif)[3:4] <- paste0(c("M", "S"), 1)
        
        # checkInput(errors = checkInput(errors = error_df, 
        #                                strata = temp_stratif, 
        #                                sampframe = ss_df))
        #> 
        
        error_df[is.na(error_df)]<-1
        
        ## run bethel at the SRS CV
        temp_bethel <- SamplingStrata::bethel(
          errors = error_df,
          stratif = temp_stratif, 
          realAllocation = T, 
          printa = T)
        
        ## Save the current n and cv constraint
        temp_n <- sum(ceiling(temp_bethel))
        updated_cv_constraint <-ifelse(dom==2,
                                       as.numeric(attributes(temp_bethel)$outcv[c(1,4), "PLANNED CV "]),
                                       as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "]))
        
          #as.numeric(attributes(temp_bethel)$outcv[c(1,4), "PLANNED CV "])
        
        ## modify CV, rerun bethel, save temp_n, run until temp_n == n_samples
        while (temp_n != n_samples){
          over_under <- temp_n > n_samples
          CV_adj <- ifelse(over_under == TRUE, 
                           yes = 1.01,
                           no = 0.999)
          
          updated_cv_constraint <- updated_cv_constraint * CV_adj
          
          error_df[, "CV1"] <- updated_cv_constraint
          
          temp_bethel <- SamplingStrata::bethel(stratif = temp_stratif,
                                                errors = error_df, 
                                                printa = TRUE)
          
          temp_n <- sum(as.numeric(temp_bethel))
          updated_cv_constraint <- ifelse(dom==2,
                                          as.numeric(attributes(temp_bethel)$outcv[c(1,4), "PLANNED CV "]),
                                          as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "]))
          print(paste0("n = ", temp_n, ", ", updated_cv_constraint) )
        }
        
        ## Save the CV and station allocations that corresponds to n_samples

          temp_idx <- ss_sample_allocations$n == n_samples & 
            ss_sample_allocations$species == spp[iispp]
        
        ss_sample_allocations[temp_idx, paste0('CV',1:dom)] <- ifelse(dom==2,
                                                                      as.numeric(attributes(temp_bethel)$outcv[c(1,4), "ACTUAL CV"]),
                                                                      as.numeric(attributes(temp_bethel)$outcv[, "ACTUAL CV"]))
        
        ss_sample_allocations[temp_idx, paste0("Str_", 1:length(temp_bethel))] <- 
          as.integer(temp_bethel)
        
      
    }
    
    ######################
    ##   8) Adjust MS solution ----
    ##   Optimize allocation across a range of sample sizes, given the original
    ##   stratification.
    #####################
    
    ms_sample_allocations <- expand.grid(n = rep(n_samples,length(unique(domain_input))))
    
    ## Subset lower limits of CVs from the ss cvs
    ss_cvs <- subset(ss_sample_allocations, n == n_samples)[,paste0('CV',1:dom)]
    ss_cvs[is.na(ss_cvs)]<-1  
    
    ## CV dataframe input to Bethel algorithm
    error_df <-  data.frame("DOM" = c(1:dom),
                            srs_cv,
                            "domainvalue"  = c(1:dom))
    
      
      names(error_df)[2:(1 + n_spp)] <- paste0("CV", 1:n_spp)
      
    
      ## Stratum statistics input to the Bethel algorithm
      temp_stratif <- solution$aggr_strata
      temp_stratif$N <- temp_stratif$N / n_years
      
      if (dom==1) {
        
        temp_stratif$DOM1 <- c(rep(1,nrow(temp_stratif)))
        
      } else if (dom==2) {
        
        temp_stratif$DOM1 <- c(rep(1,nrow(temp_stratif)/2),rep(2,nrow(temp_stratif)/2))
        temp_stratif$DOM2<-temp_stratif$DOM1
        
      }
      
      
      error_df[is.na(error_df)]<-1
      
      ## Run Bethel algorithm and save current n and cv constraints
      temp_bethel <- SamplingStrata::bethel(
        errors = error_df,
        stratif = temp_stratif, 
        realAllocation = T, 
        printa = T)
      
      temp_n <- sum(ceiling(temp_bethel))
      
      if (dom==2) {
        
        updated_cv_constraint <-  as.numeric(attributes(temp_bethel)$outcv[c(1:n_spp,((n_spp*3)+1):(n_spp*4)), "PLANNED CV "])
                                        
      } else if (dom==1) {
        
        updated_cv_constraint <-  as.numeric(attributes(temp_bethel)$outcv[c(1:n_spp), "PLANNED CV "])
        
      }
      

      ## Rerun Bethel algorithm, modifying the CVs relative to the distances
      ## between the SRS and SS CVs given n_samples stations. First we calculate 
      ## CVs calculated under SRS for each species with n_samples stations

        temp_srs_var <- srs_stats[, paste0("S", 1:n_spp)]^2 * (1 - n_samples / n_cells) / n_samples
        temp_srs_cv <- sqrt(temp_srs_var) / srs_stats[, paste0("M", 1:n_spp)]
      

      
      
      
      while (temp_n != n_samples) {
        over_under <- temp_n > n_samples
        
        ## If the current n is < n_samples, decrease the CV by a small amount
        ## relative to the distance between the current CV and the SS CV
        if (over_under == FALSE) {
          updated_cv_constraint <- updated_cv_constraint  - updated_cv_constraint * 0.001
        }
        
        ## If the current n is > n_samples, increase the CV by a small amount
        ## relative to the distance between the current CV and the SRS CV
        if(over_under == TRUE) {
          updated_cv_constraint <- updated_cv_constraint  + updated_cv_constraint * 0.01
        }
        
        ## Update the CV dataframe input with the updated_cv_constraint

          error_df[, paste0("CV", 1:n_spp)] <- updated_cv_constraint
        
        
      
        ## Rerun Bethel algorithm
        temp_bethel <- SamplingStrata::bethel(stratif = temp_stratif,
                                              errors = error_df, 
                                              printa = TRUE)
        
        temp_n <- sum(as.numeric(temp_bethel))
        
        ## Save sample size and CV constraint
        if (dom==2) {
          
          updated_cv_constraint <-  as.numeric(attributes(temp_bethel)$outcv[c(1:n_spp,((n_spp*3)+1):(n_spp*4)), "PLANNED CV "])
          
        } else if (dom==1) {
          
          updated_cv_constraint <-  as.numeric(attributes(temp_bethel)$outcv[c(1:n_spp), "PLANNED CV "])
          
        }
        
        
        ## Print out result to console
        print(paste0("n = ", temp_n) )
      }
      
      ## Save optimized CV 
      #temp_idx <- ms_sample_allocations$n == n_samples

        ms_sample_allocations[1, paste0("CV", 1:n_spp)] <- 
          as.numeric(attributes(temp_bethel)$outcv[1:n_spp, "ACTUAL CV"])
      
        if (dom==2) {
          ms_sample_allocations[2, paste0("CV", 1:n_spp)] <- 
            as.numeric(attributes(temp_bethel)$outcv[c(((n_spp*3)+1):(n_spp*4)), "ACTUAL CV"])

        }

        if (dom==1) {
          
          ms_sample_allocations[1, paste0("Str_", 1:(length(temp_bethel)))] <- 
            as.integer(temp_bethel)[1:(length(temp_bethel))]
          
        } else if (dom==2) {
          
          ms_sample_allocations[1, paste0("Str_", 1:(length(temp_bethel)/2))] <- 
            as.integer(temp_bethel)[1:(length(temp_bethel)/2)]
          ms_sample_allocations[2, paste0("Str_", 1:(length(temp_bethel)/2))] <- 
            as.integer(temp_bethel)[((length(temp_bethel)/2)+1):length(temp_bethel)]
          
        }
  
        
        

      
      #store number samples per strata
    samples_strata<-as.integer(temp_bethel)

    #store bethel CV
    cv_bethel_final<-error_df
    
    #cv
    cv_temp <- rbind(CV_random=cv_initial, #max
                          CV_strata=cv_strata_final, #min to achieve strata
                          CV_bethel=cv_bethel_final) #min to achieve naximum sample
    
    #results list to save
    all<-list(result_list=result_list,
              ss_sample_allocations=ss_sample_allocations,
              ms_sample_allocations=ms_sample_allocations,
              samples_strata=data.frame(strata=1:no_strata,n_samples=samples_strata),
              cv=cv_temp)
    
    #strata to plot
    dd<-all$result_list$solution$framenew
    #all$ms_sample_allocations
    effort<-all$ms_sample_allocations[,colnames(all$ms_sample_allocations)[grepl(pattern = 'Str',colnames(all$ms_sample_allocations))]]
    effort1<-reshape2::melt(effort)
    names(effort1)<-c('STRATO','effort')
    effort1$STRATO<-gsub('Str_','',effort1$STRATO)
    dd1<-merge(dd[,c('DOMAINVALUE','STRATO','ID')],grid,by.x='ID',by.y='cell')
    dd2<-merge(dd1,effort1,by='STRATO')
    dd2$strata<-paste0(dd2$DOMAINVALUE,'_',dd2$STRATO)
    
    #save list
      save(all,
           file = paste0("./output/ms_optim_allocations_ebsnbs_slope_",samp_df[s,'samp_scn'],'_','DOM',length(unique(domain_input)),".RData"))
      
      
      #save plot
      ragg::agg_png(paste0('./figures slope/',"ms_optim_allocations_ebsnbs_slope_",samp_df[s,'samp_scn'],'.png'),  width = 7, height = 7, units = "in", res = 300)
      print(
        ggplot()+
          geom_tile(data=dd2,aes(x=Lon,y=Lat,color=strata),size=1.5)+
          geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
          geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
          geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
          ggtitle(samp_df[s,'strat_var'])#+
          #scale_color_continuous(type='viridis')
      )
      dev.off()
      
      #save plot
      ragg::agg_png(paste0('./figures slope/',"ms_optim_allocations_ebsnbs_slope_",samp_df[s,'samp_scn'],'_effort.png'),  width = 7, height = 7, units = "in", res = 300)
      print(
        ggplot()+
          geom_tile(data=dd2,aes(x=Lon,y=Lat,color=effort),size=1.5)+
          geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
          geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
          geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
          ggtitle(samp_df[s,'strat_var'])+
          scale_color_continuous(type='viridis')
      )
      dev.off()
      
      
   }

  for (s in c(2,4,6)) { #nrow(samp_df)
  
    #2,4,6
    
    #s<-6
    
    dom<-samp_df[s,'domain']
    
    
    #load list
    load(
      file = paste0("./output/ms_optim_allocations_ebsnbs_slope_",samp_df[s,'samp_scn'],'_','DOM',dom,".RData")) #all
    
    #strata to plot
    dd<-all$result_list$solution$framenew
    #all$ms_sample_allocations
    effort<-all$ms_sample_allocations[,colnames(all$ms_sample_allocations)[grepl(pattern = 'Str',colnames(all$ms_sample_allocations))]]
    effort1<-reshape2::melt(effort)
    names(effort1)<-c('STRATO','effort')
    effort1$STRATO<-gsub('Str_','',effort1$STRATO)
 
    
    if (dom==2) {
      
      effort1$STRATO<-paste0(c(1,2,1,2,1,2,1,2,1,2),'_',effort1$STRATO)
      names(effort1)[1]<-'strata'
      dd1<-merge(dd[,c('DOMAINVALUE','STRATO','ID')],grid,by.x='ID',by.y='cell')
      dd1$strata<-paste0(dd1$DOMAINVALUE,'_',dd1$STRATO)
      dd2<-merge(dd1,effort1,by='strata')
      #dd2$strata<-paste0(dd2$DOMAINVALUE,'_',dd2$STRATO)
      
    } else {
      
    dd1<-merge(dd[,c('DOMAINVALUE','STRATO','ID')],grid,by.x='ID',by.y='cell')
    dd2<-merge(dd1,effort1,by='STRATO')
    dd2$strata<-paste0(dd2$DOMAINVALUE,'_',dd2$STRATO)
    }
    
    if (dom==2) {
      
    nbs_effort<-aggregate(ID ~ X1,all$result_list$sol_by_cell[1:15180,],FUN=length)
    nbs_effort$X1<-paste0('1_',nbs_effort$X1)
    
    ebs_effort<-aggregate(ID ~ X1,all$result_list$sol_by_cell[15181:nrow(all$result_list$sol_by_cell),],FUN=length)
    ebs_effort$X1<-paste0('2_',ebs_effort$X1)
  
    all_effort<-rbind(nbs_effort,ebs_effort)
    names(all_effort)<-c('strata','effort_area')
    dd2<-merge(dd2,all_effort,by=('strata'))
    dd2$effort_area<-dd2$effort_area/dd2$effort
    
    } else {
      
      all_effort<- aggregate(ID ~ X1,all$result_list$sol_by_cell,FUN=length)
      names(all_effort)<-c('strata','effort_area')
      all_effort$strata<-paste0('1_',all_effort$strata)
      dd2<-merge(dd2,all_effort,by=('strata'))
      dd2$effort_area<-dd2$effort_area/dd2$effort
    }
  #save plot
  #ragg::agg_png(paste0('./figures slope/',"ms_optim_allocations_ebsnbs_slope_",samp_df[s,'samp_scn'],'.png'),  width = 7, height = 7, units = "in", res = 300)
  print(
    ggplot()+
      geom_point(data=dd2,aes(x=Lon,y=Lat,color=strata),size=0.3)+
      geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
      geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
      geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
      ggtitle(samp_df[s,'strat_var'])#+
    #scale_color_discrete(type='viridis')
  )
  #dev.off()
  
  #save plot
  #ragg::agg_png(paste0('./figures slope/',"ms_optim_allocations_ebsnbs_slope_",samp_df[s,'samp_scn'],'_effort.png'),  width = 7, height = 7, units = "in", res = 300)
  print(
    ggplot()+
      geom_point(data=dd2,aes(x=Lon,y=Lat,color=effort_area),size=0.3)+
      geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
      geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
      geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
      ggtitle(samp_df[s,'strat_var'])+
      theme_minimal()+
      theme(text = element_blank())+
      scale_color_continuous(type='viridis',name='stations/area')
  )
  #dev.off()
  }
  
  
  
  ###################
  # Plot comparison sampling effort x strata for each species under singlesp or multisp allocation of samples
  ###################
    
  cvs<-data.frame(matrix(NA,nrow = 0,ncol = 4))
  names(cvs)<-c('sp','samp','ss','ms')
  
  for (s in 1:nrow(samp_df)) {
    
    #s<-1
    
    samp<-samp_df[s,'samp_scn']
    
    #load multispecies data
    load(paste0('./output/multisp_optimization_static_data_ebsnbs_slope.RData')) #df
    
    #load optimized stratification
    #load(file = paste0("./output/ms_optim_allocations_ebsnbs_slope_",samp_df[s,'samp_scn'],".RData")) #all
    load(paste0("./output/ms_optim_allocations_ebsnbs_slope_",samp_df[s,'samp_scn'],'_','DOM',samp_df[s,'domain'],".RData"))
    
    strata<-rbind(all$result_list$solution$indices,
                  data.frame(ID=rem_cells,X1=NA))
    colnames(strata)<-c('cell','Strata')
    
    #n cells by strata to calculate proportion of sampling effort as n samples/total cells
    strata_sum<-aggregate(cell ~ Strata, data = strata, FUN = length)
    names(strata_sum)[2]<-'total_cell'
    
    #merge with ss data
    ss<-all$ss_sample_allocations
    ss1<-data.frame(ss[,c(paste0('Str_',1:samp_df$n_strata[s]))],row.names = ss$species)
    colnames(ss1)<-1:samp_df$n_strata[s]
    ss2<-reshape2::melt(as.matrix(ss1))
    names(ss2)<-c('species','Strata','ss_samples')
    strata1<-merge(strata,ss2,by='Strata')
    
    #merge with ms data
    ms<-all$ms_sample_allocations
    ms_strata<-data.frame('Strata'=1:samp_df$n_strata[s],
                          'ms_samples'=as.numeric(ms[1,20:(20+samp_df$n_strata[s]-1)]))
    
    
    strata2<-merge(strata1,ms_strata,by='Strata')
    strata2<-merge(strata2,strata_sum,by='Strata')
    #strata2$prop<-strata2$n_samples/strata2$total_cell
    strata2$ratio<-log(strata2$ss_samples/strata2$ms_samples)
    dim(strata2)
    
    df1<-df[,c('Lat','Lon','cell')]
    df1<-merge(df1,strata2,by='cell')
    df1$Strata[is.na(df1$Strata)]<-999
    
    #to store plots
    plot_list_n<-list()
    plot_list_d<-list()
    
    for (isp in spp) {
      
    #isp<-spp[1]
    
    df2<-subset(df1,species==isp)
    
    #ms CV
    ms_cv<-
      ms[,paste0('CV',match(isp,spp))]
    
    #ms SCV
    ms_scv<-(ms_cv*100)^2
    
    #ss CV
    ss_cv<-
      ss$CV[match(isp,spp)]
    
    #ss SCV
    ss_scv<-(ss_cv*100)^2
    
    icvs<-data.frame('sp'=isp,
                     'samp'=samp,
                     'ss'=ss_cv,
                     'ms'=ms_cv)
    
    #append
    cvs<-rbind(cvs,icvs)
    
    #df to spatialpoint df
    coordinates(df2) <- ~ Lon + Lat
    proj4string(df2)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    
    #reproject coordinates for plotting purposes
    df_1<-spTransform(df2,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
    df_2<-data.frame(df_1)
    
    #x and y cells
    xycells<-as.integer(sqrt(dim(df_1)[1]))
    
    # create a template raster
    r1 <- raster(ext=extent(df_1),ncol=xycells, nrow=xycells) #c(15800,15800) 7000
    
    #create raster
    r2<-rasterize(df_1, r1 ,field=c('Strata','ms_samples','ss_samples','ratio'))
    crs(r2) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
    r2[r2==999]<-NA
    
    #create polygon to get boundaries of each strata
    r3<-as.data.frame(r2,xy=TRUE)
    r4<-rasterToPolygons(r2$Strata,dissolve=TRUE,digits = 1)
    
    #color palette
    pal <- wes_palette("Zissou1", length(sort(unique(r3$ss_samples))), type = "continuous")
    color_scale_ss <- setNames(as.character(pal), sort(unique(r3$ss_samples)))
    r3<-r3[complete.cases(r3$Strata),] 
    
    #as factors
    r3$ss_samples<-as.factor(r3$ss_samples)
    #r3$ratio<-as.factor(r3$ratio)
  
    #common name
    com<-spp_name[which(spp_name$spp==isp),'common']
    
    #plot by number of samples
    pn<-
    ggplot()+
      geom_raster(data=r3,aes(x=x,y=y,fill=ss_samples))+
      geom_polygon(data=r4,aes(x=long,y=lat,group=group), colour="black", fill=NA)+
      scale_fill_manual(values = color_scale_ss,name='sampling effort')+
      guides(fill=guide_colorbar(title.position = 'top', title.hjust = 0.5,ticks.colour = NA,frame.colour = 'black'))+
      scale_x_continuous(expand = c(0,0),position = 'bottom',
                         breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
      coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
               xlim = panel_extent$x-c(000000,100000),
               ylim = panel_extent$y-c(0,300000),
               label_axes = "-NE-")+
      theme(aspect.ratio = 1,panel.grid.major = element_blank(),
            panel.background = element_rect(fill = NA),panel.ontop = TRUE,
            legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(20, 'points'),
            legend.key.width= unit(20, 'points'),axis.title = element_blank(),legend.position = 'none',
            panel.border = element_rect(fill = NA, colour = NA),legend.key = element_rect(color="black"),
            legend.spacing.y = unit(8, 'points'),
            axis.text=element_blank(),axis.ticks = element_blank(),
            plot.margin = margin(0.01,0.01,0.01,0.01), 
            axis.ticks.length = unit(-5,"points"),plot.title = element_text(size=12,hjust = 0.5,vjust=-5, face="bold"))+
      scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())+
      #labs(title=paste0(com,'\n(msSCV=',round(ms_scv,digits = 3),'; ssSCV=',round(ss_scv,digits = 3),')'),fill='')
      labs(title=paste0(com,'\n(msCV=',round(ms_cv,digits = 3),'; ssCV=',round(ss_cv,digits = 3),')'),fill='')
    
    #plot by number of samples
    pd<-
    ggplot()+
      geom_raster(data=r3,aes(x=x,y=y,fill=ratio))+
      geom_polygon(data=r4,aes(x=long,y=lat,group=group), colour="black", fill=NA)+
      scale_fill_gradient2(midpoint = 0, low = "red", mid = "white",
                             high = "blue",name='ss samples - ms samples')+
      guides(fill=guide_colorbar(title.position = 'top', title.hjust = 0.5,ticks.colour = NA,frame.colour = 'black'))+
      scale_x_continuous(expand = c(0,0),position = 'bottom',
                         breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
      coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
               xlim = panel_extent$x-c(000000,100000),
               ylim = panel_extent$y-c(0,300000),
               label_axes = "-NE-")+
      theme(aspect.ratio = 1,panel.grid.major = element_blank(),
            panel.background = element_rect(fill = NA),panel.ontop = TRUE,
            legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(20, 'points'),
            legend.key.width= unit(20, 'points'),axis.title = element_blank(),legend.position = 'none',
            panel.border = element_rect(fill = NA, colour = NA),legend.key = element_rect(color="black"),
            legend.spacing.y = unit(8, 'points'),
            axis.text=element_blank(),axis.ticks = element_blank(),
            plot.margin = margin(0.01,0.01,0.01,0.01), 
            axis.ticks.length = unit(-5,"points"),plot.title = element_text(size=12,hjust = 0.5,vjust=-5, face="bold"))+
      scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())+
      #labs(title=paste0(com,'\n(msSCV=',round(ms_scv,digits = 3),'; ssSCV=',round(ss_scv,digits = 3),')'),fill='')
      labs(title=paste0(com,'\n(msCV=',round(ms_cv,digits = 3),'; ssCV=',round(ss_cv,digits = 3),')'),fill='')
    
    plot_list_n[[isp]]<-pn
    plot_list_d[[isp]]<-pd
    }
  
    
    #color palette
    pal <- wes_palette("Zissou1", length(sort(unique(r3$ms_samples))), type = "continuous")
    color_scale_ms <- setNames(as.character(pal), sort(unique(r3$ms_samples)))
    
    #as factors
    r3$ms_samples<-as.factor(r3$ms_samples)
    #r3$ratio<-as.factor(r3$ratio)
    
    #plot by number of samples
    pm<-
      ggplot()+
      geom_raster(data=r3,aes(x=x,y=y,fill=ms_samples))+
      geom_polygon(data=r4,aes(x=long,y=lat,group=group), colour="black", fill=NA)+
      scale_fill_manual(values = color_scale_ms,name='sampling effort')+
      guides(fill=guide_colorbar(title.position = 'top', title.hjust = 0.5,ticks.colour = NA,frame.colour = 'black'))+
      scale_x_continuous(expand = c(0,0),position = 'bottom',
                         breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
      coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
               xlim = panel_extent$x-c(000000,100000),
               ylim = panel_extent$y-c(0,300000),
               label_axes = "-NE-")+
      theme(aspect.ratio = 1,panel.grid.major = element_blank(),
            panel.background = element_rect(fill = NA),panel.ontop = TRUE,
            legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(20, 'points'),
            legend.key.width= unit(20, 'points'),axis.title = element_blank(),legend.position = 'none',
            panel.border = element_rect(fill = NA, colour = NA),legend.key = element_rect(color="black"),
            legend.spacing.y = unit(8, 'points'),
            axis.text=element_blank(),axis.ticks = element_blank(),
            plot.margin = margin(0.01,0.01,0.01,0.01), 
            axis.ticks.length = unit(-5,"points"),plot.title = element_text(size=14,hjust = 0.5,vjust=-5, face="bold"))+
      scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())+
      labs(title=paste0('multispecies\n',gsub('_',' + ',samp_df[s,'strat_var'])),fill='')  
      
  #create legends
  legend_d<-
    ggplot()+
    geom_raster(data=r3,aes(x=x,y=y,fill=as.numeric(ratio)))+
    scale_fill_gradient2(midpoint = mean(range(r3$ratio)), low = "red", mid = "white",
                         high = "blue",breaks=range(as.numeric(r3$ratio)),labels=c("Lower","Higher"),name='log(ss/ms samples)')+
    guides(fill=guide_colorbar(title.position = 'top', title.hjust = 0.5,ticks.colour = NA,frame.colour = 'black'))+
    theme(legend.position = 'bottom')+
    labs(fill='')
  
  legend_n<-
    ggplot()+
    geom_raster(data=r3,aes(x=x,y=y,fill=as.numeric(ss_samples)))+
    scale_fill_gradientn(colours = pal,breaks=range(as.numeric(r3$ss_samples)),labels=c("Low","High"),name='sampling effort (n samples)')+
    guides(fill=guide_colorbar(title.position = 'top', title.hjust = 0.5,ticks.colour = NA,frame.colour = 'black'))+
    theme(legend.position = 'bottom')+
    labs(fill='')
  
  legend1 <- cowplot::get_legend( 
    legend_d + 
      theme(legend.position = "bottom") 
  ) 
  
  legend2 <- cowplot::get_legend( 
    legend_n + 
      theme(legend.position = "bottom") 
  ) 
  
  pgrid1<-cowplot::plot_grid(plotlist = plot_list_d, nrow = 2)
  pgrid2<-cowplot::plot_grid(plotlist = plot_list_n, nrow = 2)
  
  #save plots
  ragg::agg_png(paste0('./figures/sampling designs ss ratio_',samp_df[s,'strat_var'],'.png'), width = 20, height = 7, units = "in", res = 300)
  print(cowplot::plot_grid(pgrid1, legend1, nrow = 2, rel_heights = c(1, .1)))
  dev.off()
  
  ragg::agg_png(paste0('./figures/sampling designs ss n_',samp_df[s,'strat_var'],'.png'), width = 20, height = 7, units = "in", res = 300)
  print(cowplot::plot_grid(pgrid2, legend2, nrow = 2, rel_heights = c(1, .1)))
  dev.off()
  
  ragg::agg_png(paste0('./figures/sampling designs ms_',samp_df[s,'strat_var'],'.png'), width = 5, height = 5, units = "in", res = 300)
  print(pm)
  dev.off()
  
}  
  
  save(cvs,file = './output/ss_ms_slope.RData')
  load(file = './output/ss_ms_slope.RData')
  spp_name$common<-gsub('_EBSNBS','',spp_name$common)
  cvs1<-merge(cvs,spp_name,by.x='sp',by.y='spp')
  cvs1<-cvs1[order(cvs1$common,decreasing = FALSE),]
  cvs1$common<-factor(cvs1$common,levels = rev(sort(unique(cvs1$common))))
  #cvs1$samp<-factor(cvs1$samp,
  #               levels = c('scn3','scn2','scn1'))
  
  levels(cvs1$samp)
  
  p1<-
  ggplot()+
    geom_linerange(data=cvs1,aes(xmin=ss,xmax=ms,x=ms,y=common,color=samp),linewidth=1,stat = "identity", position = position_dodge(width = 0.7))+ 
    geom_point(data=cvs1,aes(x=ss,y=common,group=samp),shape=4,stat = "identity", position = position_dodge(width = 0.7),size=2)+
    geom_point(data=cvs1,aes(x=ms,y=common,group=samp),shape=16,stat = "identity", position = position_dodge(width = 0.7),size=2)+
    scale_color_manual(values=c('scn1'='#4b7a99','scn2'='#679bc3','scn3'='#8db6c3'),
                      labels = c('opt depth','opt varSBT','opt depth + varSBT'),
                      limits=c('scn3','scn2','scn1'),
                      name='stratification')+
    # scale_color_manual(values=c('scn1'='#4b7a99','scn2'='#679bc3','scn3'='#8db6c3'),
    #                   labels = c('opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
    theme_bw()+
    theme(panel.grid.minor = element_line(linetype=2,color='grey90',),#strip.background = element_rect(fill='white'),
          legend.key.size = unit(12, 'points'),legend.direction = 'vertical',legend.text = element_text(size=9),legend.position=c(.78,.925),
          legend.title = element_text(size=10),legend.spacing.x = unit(0.05, "cm"),legend.box.spacing =  unit(0.01, "cm"), #,strip.text = element_text(size=12)
          strip.background = element_blank(), legend.box.background = element_rect(fill = "white", color = "black"),#legend.background = element_blank(),legend.box = 'horizontal',#legend.justification = 'right',legend.position='bottom',#
          strip.text = element_blank())+ #axis.text.x = element_text(angle=90,vjust=0.5),
    expand_limits(x = 0)+
    labs(y='',x='CV')+
    scale_x_continuous(limits = c(0,max(cvs1$ms)+max(cvs1$ms)*0.1),expand = c(NA,0)) #expand = c(NA,0.1),limits = c(0,NA)
  
  cvs1$samp<-factor(cvs1$samp,
                    levels = c('scn3','scn2','scn1'))
  
  p<- 
  ggplot()+
    #geom_linerange(data=cvs1,aes(xmin=ss,xmax=ms,x=ms,y=common,color=samp),linewidth=1,stat = "identity", position = position_dodge(width = 0.7))+ 
    geom_point(data=cvs1,aes(x=log(ms/ss),y=rev(common),group=samp,fill=samp),stat = "identity", position = position_dodge(width = 0.5),size=3,shape=21)+
    #geom_point(data=cvs1,aes(x=ms,y=common,group=samp),shape=16,stat = "identity", position = position_dodge(width = 0.7),size=2)+
      scale_fill_manual(values=c('scn1'='#4b7a99','scn2'='#679bc3','scn3'='#8db6c3'),
                        labels = c('opt depth','opt varSBT','opt depth + varSBT'),
                        limits=c('scn3','scn2','scn1'),
                        name='stratification')+
    # scale_color_manual(values=c('scn1'='#4b7a99','scn2'='#679bc3','scn3'='#8db6c3'),
    #                   labels = c('opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
    theme_bw()+
    theme(panel.grid.minor = element_line(linetype=2,color='grey90',),#strip.background = element_rect(fill='white'),
          legend.key.size = unit(12, 'points'),legend.direction = 'vertical',legend.text = element_text(size=9),legend.position=c(.78,.925),
          legend.title = element_text(size=10),legend.spacing.x = unit(0.05, "cm"),legend.box.spacing =  unit(0.01, "cm"), #,strip.text = element_text(size=12)
          strip.background = element_blank(), legend.box.background = element_rect(fill = "white", color = "black"),#legend.background = element_blank(),legend.box = 'horizontal',#legend.justification = 'right',legend.position='bottom',#
          strip.text = element_blank())+ #axis.text.x = element_text(angle=90,vjust=0.5),
    expand_limits(x = 0)+
    labs(y='',x='log(CVms/CVss)')+
    scale_x_continuous(limits = c(0,max(log(cvs1$ms/cvs1$ss))+max(log(cvs1$ms/cvs1$ss))*0.1),expand = c(NA,0)) #expand = c(NA,0.1),limits = c(0,NA)
  
  #save plot
  ragg::agg_png(paste0('./figures/CVss_CVms_ratio.png'),  width = 5, height = 7, units = "in", res = 300)
  p
  dev.off()
  
  #save plot
  ragg::agg_png(paste0('./figures/CVss_CVms.png'),  width = 5, height = 7, units = "in", res = 300)
  p1
  dev.off()
  
  
  ###################
  # Plot spatial random fields
  ###################
  plot_list_rf<-list()
  
  for (isp in spp) {
    
    #isp<-'Gadus macrocephalus'
    cat(paste(" #############  ",isp ," #############\n"))
    
    #common name
    com<-spp_name[which(spp_name$spp==isp),'common']
    
    #fit file
    ff<-list.files(paste0('./shelf EBS NBS VAST/',isp,'/'),'fit',recursive=TRUE)
    
    #load fit file
    load(paste0('./shelf EBS NBS VAST/',isp,'/',ff)) #fit
    
    #get spatial random fields
    spf<-fit$Report$Omega1_gc
    
    #get dataframe
    D_gt <- spf
    D_gt<-data.frame('cell'=c(1:fit$spatial_list$n_g),D_gt,check.names = FALSE)
    
    #merge with grid
    D_gt1<-merge(D_gt,grid,by='cell')
    
    #df to spatialpoint df
    coordinates(D_gt1) <- ~ Lon + Lat
    crs(D_gt1)<-c(crs='+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
    
    #reproject coordinates for plotting purposes
    #df_1<-spTransform(D_gt1,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
    df_2<-data.frame(D_gt1)
    
    #x and y cells
    xycells<-as.integer(sqrt(dim(D_gt1)[1]))
    
    # create a template raster
    r1 <- raster(ext=extent(df_1),ncol=xycells, nrow=xycells) #c(15800,15800) 7000
    
    #create raster
    r2<-rasterize(D_gt1, r1 ,field=c('1'))
    crs(r2) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
    r2[r2==999]<-NA
    
    #create polygon to get boundaries of each strata
    r3<-as.data.frame(r2,xy=TRUE)
    
    #Alaska layer
    ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")
    ak_sppoly<-as(ebs_layers$akland, 'Spatial')
    
    #plot
    prf<-
    ggplot() +
      geom_raster(data=r3,aes(x=x,y=y,fill=layer))+
      #geom_polygon(data=r4,aes(x=long,y=lat,group=group), colour="black", fill=NA)+
      guides(fill=guide_colorbar(title.position = 'top', title.hjust = 0.5,ticks.colour = NA,frame.colour = 'black'))+
      scale_x_continuous(expand = c(0,0),position = 'bottom',
                         breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
      coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
               xlim = panel_extent$x-c(000000,100000),
               ylim = panel_extent$y-c(0,300000),
               label_axes = "-NE-")+
      scale_fill_viridis_c(option = 'A',name=('Spatial Random\nField deviations'),
                           guide = guide_colorbar(  frame.colour = "black",ticks.colour = 'black'),na.value=rgb(1, 0, 0, 0))+
      #scale_fill_gradientn(colours = wesanderson::wes_palette("Zissou1", 21, type = "continuous"),name=('Spatial Random\nField deviations'),
      #                      guide = guide_colorbar(  frame.colour = "black",ticks.colour = 'black'),na.value=rgb(1, 0, 0, 0))+
      theme(aspect.ratio = 1,panel.grid.major = element_blank(),
            panel.background = element_rect(fill = NA),panel.ontop = TRUE,
            legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(20, 'points'),
            legend.key.width= unit(20, 'points'),axis.title = element_blank(),legend.position = 'none',
            panel.border = element_rect(fill = NA, colour = NA),legend.key = element_rect(color="black"),
            legend.spacing.y = unit(8, 'points'),
            axis.text=element_blank(),axis.ticks = element_blank(),
            plot.margin = margin(0.01,0.01,0.01,0.01), 
            axis.ticks.length = unit(-5,"points"),plot.title = element_text(size=12,hjust = 0.5,vjust=-5, face="bold"))+
      scale_y_continuous(expand = c(0,0),position = 'right',sec.axis = dup_axis())+
      labs(title=paste0(com),fill='')
  
    #store plot
    plot_list_rf[[isp]]<-prf
  } 
  
  #grid of plots
  pgrid1<-cowplot::plot_grid(plotlist = plot_list_rf, nrow = 2)
  
  #plot for common legend
  legend_rf<-
    ggplot()+
    geom_raster(data=r3,aes(x=x,y=y,fill=as.numeric(layer)))+
    scale_fill_viridis_c(breaks=range(as.numeric(r3$layer),na.rm = TRUE),labels=c("Low","High"),option = 'A',name=('Spatial Random\nField deviations'),
                         guide = guide_colorbar(  frame.colour = "black",ticks.colour = 'black'),na.value=rgb(1, 0, 0, 0))+
    guides(fill=guide_colorbar(title.position = 'top', title.hjust = 0.5,ticks.colour = NA,frame.colour = 'black'))+
    theme(legend.position = 'bottom')+
    labs(fill='')
  
  #legend
  legend1 <- cowplot::get_legend( 
    legend_rf + 
      theme(legend.position = "bottom") 
  ) 
  
  #save plot
  ragg::agg_png(paste0('./figures/SpatialRandomFields.png'),  width = 20, height = 7, units = "in", res = 300)
  print(cowplot::plot_grid(pgrid1, legend1, nrow = 2, rel_heights = c(1, .15)))
  dev.off()
  