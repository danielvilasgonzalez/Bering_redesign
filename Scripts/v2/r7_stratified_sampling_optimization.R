####################################################################
####################################################################
##
##    Run sampling optimization based on predicted densities from VAST model
##    and get samples from each strata for each sampling design
##
##    by best stratification, we mean the stratification that ensures the minimum sample cost, 
##    sufficient to satisfy precision constraints set on the accuracy of the estimates of the survey target variables Y’s
##    constraints expressed as maximum allowable coefficients of variation in the different domains of interest
##
##    https://cran.r-project.org/web/packages/SamplingStrata/vignettes/SamplingStrata.html
##    (using devtools::install_github(repo = "zoyafuso-NOAA/SamplingStrata"))
##    Daniel Vilas (daniel.vilas@noaa.gov/dvilasg@uw.edu)
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c("splines",'SamplingStrata','wesanderson','dplyr','sp',
             'sf','maptools','parallel','rasterVis','rgeos','scales',
             'rnaturalearth','grid','ggplot2','spatstat','parallel','doParallel')

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
cl <- parallel::makeCluster(parallel::detectCores() - 2)
doParallel::registerDoParallel(cl)

#setwd
out_dir<-'C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/'  #out_dir<-'/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/'
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
       'Paralithodes camtschaticus')

df_spp<-data.frame('spp'=spp,
                   'n'=c(1:length(spp)),
                   'Y'=paste0('Y',c(1:length(spp))))

tar_var<-paste0(rep(df_spp$Y,each=2),c('','_SQ_SUM'))

#number sp
n_spp<-length(spp)

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
x3$x<-as.integer(x3$x)
x3$y<-as.integer(x3$y)
lon<-sort(unique(x3$x),decreasing = FALSE) #1556
lat<-sort(unique(x3$y),decreasing = TRUE) #1507
lons<-data.frame(x=lon,col=1:length(lon))
lats<-data.frame(y=lat,row=1:length(lat))
x4<-merge(x3,lons,by='x',all.x=TRUE)
x5<-merge(x4,lats,by='y',all.x=TRUE)
colnames(x5)<-c('Lat','Lon','cell','optional','col','row')
grid<-x5[,c('Lat','Lon','cell','col','row')]

###################################
# FXNs
###################################

calc_expected_CV <- function (strata) {
  
  n_spp <- length(grep(x = names(strata), pattern = "M"))
  n_h <- strata$Allocation
  N_h <- strata$Population
  
  cv <- vector(length = n_spp)
  names(cv) <- paste0("Y", 1:n_spp)
  
  for (ispp in 1:n_spp) {
    S_h <- strata[, paste0("S", ispp)]
    M_h <- strata[, paste0("M", ispp)]
    
    Y_h <- N_h * M_h
    Var_h <- (N_h^2) * (1 - n_h/N_h) * ((S_h^2)/n_h)
    CV <- sqrt(sum(Var_h))/sum(Y_h)
    
    cv[ispp] <- CV
  }
  
  cv <- round(cv, 3)
  return(cv)
}

###################################
# BASELINE/CURRENT SAMPLING DESIGN
###################################

load('./output/baseline_strata.RData') #baseline_strata

###################################
# SCENARIOS
###################################

#sampling scenarios
samp_df<-expand.grid(strat_var=c('Depth_varTemp','varTemp','Depth'),
                     target_var=c('sumDensity'), #,'sqsumDensity'
                     n_samples=c(520), #c(300,500) 520 (EBS+NBS+CRAB);26 (CRAB); 350 (EBS-CRAB); 494 (NBS-CRAB)
                     n_strata=c(15)) #c(5,10,15)

#add scenario number
samp_df$samp_scn<-paste0(paste0('scn',1:nrow(samp_df)))

#########################
# RUN LOOP SPECIES
#########################

#load multispecies data
load(paste0('./output/species/multisp_optimization_static_data.RData')) #df
names(df)[((ncol(df)-length(tar_var))+1):ncol(df)]<-tar_var

#loop over species
#for (sp in spp) {
  
  #sp<-'Gadus macrocephalus'
  
  #load optimization data
  #load(paste0('./output/species/',sp,'/optimization data/optimization_static_data.RData')) #D6
  #load(paste0('./output/species/',sp,'/projection_data.RData')) #temp_dens_vals
  
  #load fit OM
  #load(paste0('./shelf EBS NBS VAST/',sp,'/fit.RData'))
  
  #removed cells because of depth
  rem_cells<-df[which(df$include==FALSE),'cell']
  ok_cells<-df[which(df$include==1),'cell']
  
  #load data_geostat file
  #data_geostat<-readRDS(paste0('./data processed/species/',sp,'/','data_geostat_temp.rds')) 
  
  #n cells
  n_cells<-length(ok_cells)
  
  #years
  n_years<-length(1982:2022)
  
  #number of domains
  n_dom<-1
  
  #domain_input
  domain_input<-rep(1, n_cells)
  
  # #df summary
  # sp_sum_stats<-data.frame(matrix(nrow = 0,ncol=18))
  # names(sp_sum_stats)<- c("Domain","Stratum","Population","Allocation","SamplingRate","Lower_X1","Upper_X1","Lower_X2","Upper_X2",
  #                         "stratum_id","wh","Wh","M1","S1","SOLUZ","samp_scn","sp")
  # 
  #subset cells with appropiate depth
  static_df1<-subset(df,cell %in% ok_cells)
  
  # #frame_df
  # frame_df <- data.frame(domainvalue = domain_input, #domain
  #                     id = static_df1$cell, #id as cells
  #                     stratum_var_input, #Stratification variables
  #                     WEIGHT=n_years, #weight for spp depending on the presence of years
  #                     target_var_input) #target variables 
  # 
  # srs_stats <- SamplingStrata::buildStrataDF(
  #   dataset = cbind( frame_df[, -grep(x = names(frame_df), pattern = "X")],
  #                    X1 = 1))
  # 
  # 
  # srs_stats <- SamplingStrata::buildStrataDF(
  #   dataset = cbind( frame_df[, -grep(x = names(frame_df), pattern = "X")],
  #                    X1 = 1))
  # ## SRS statistics
  # srs_var <- srs_stats[, paste0("S", 1:n_spp)]^2 * (1 - srs_n / n_cells) / srs_n
  # srs_cv <- sqrt(srs_var) / srs_stats[, paste0("M", 1:n_spp)]
  
  
  #########################
  # loop over optimized sampling designs
  #########################
  
  #loop through sampling designs
  for (s in 1:nrow(samp_df)) {
    
    #s<-2
    
    #print scenario to check progress
    cat(paste(" #############  Sampling Scenario", samp_df[s,"samp_scn"], " #############\n"))
    
    #if scenario includes two stratifying factors
    if (grepl('_',samp_df[s,'strat_var'])) {
      #stratification variables 
      stratum_var_input<-data.frame(X1 = static_df1[,paste0(sub("\\_.*", "", samp_df[s,'strat_var']))],
                                    X2 = static_df1[,paste0(sub(".*_", "", samp_df[s,'strat_var']))]) 
    } else {
      stratum_var_input<-data.frame(X1 = static_df1[,paste0(sub("\\_.*", "", samp_df[s,'strat_var']))])
    }
    
    # srs_stats <- SamplingStrata::buildStrataDF(
    #   dataset = cbind( static_df1[, -grep(x = names(static_df1), pattern = "X")],
    #                    X1 = 1))
    
    ## SRS statistics
    #srs_var <- srs_stats[, paste0("S", 1:n_spp)]^2 * (1 - srs_n / n_cells) / srs_n
    #srs_cv <- sqrt(srs_var) / srs_stats[, paste0("M", 1:n_spp)]
    
    
    #target variables
    target_var_input<-data.frame(Y1 = static_df1$sumDensity,
                                 Y1_SQ_SUM = static_df1$sumDensity_sq) #D7$sqsumDensity #Ynspp #set different scenarios and spp ############ TO CHECK
    
    target_var_input<-static_df1[,tar_var]
    
    #create df
    frame <- data.frame(domainvalue = domain_input, #domain
                        id = static_df1$cell, #id as cells
                        stratum_var_input, #Stratification variables
                        WEIGHT=n_years, #weight for spp depending on the presence of years
                        target_var_input) #target variables 
    
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
    srs_var <- srs_stats[, paste0("S", 1:n_spp)]^2 * (1 - srs_n / n_cells) / srs_n
    srs_cv <- sqrt(srs_var) / srs_stats[, paste0("M", 1:n_spp)]
    
    #CV
    cv_df <- list()
    cv_df[["DOM"]] <- 1
    for (ispp in 1:n_spp) cv_df[[paste0("CV", ispp)]] <- as.numeric(srs_cv[ispp])
    cv_df[["domainvalue"]] <- 1
    cv_df <- as.data.frame(cv_df)
    
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
    
    #while loop until the desired number of strata  
    flag<-TRUE
    while (flag) {
      
      #print state
      cat(paste(" #############   OPTIMIZING STRATA - CV of ",spp[1],' ', cv_df[,2],"  #############\n"))
      
      #run optimization
      solution <- optimStrata(method = "continuous", #continous variables
                              errors = cv_df,  #precision level - maximum allowable coefficient of variation set by the simple random sampling 
                              framesamp = frame, #df of input variables
                              iter = 20, #300 #aximum number of iterations
                              pops = 10, #100  #dimension of each generations
                              elitism_rate = 0.1, #0.1
                              mut_chance = 1 / (no_strata[1] + 1), #mutation chance
                              nStrata = no_strata, #maximum strata
                              showPlot = TRUE, #FALSE
                              writeFiles = FALSE)
      
      #flag to keep the loop
      flag<-ifelse(nrow(solution$aggr_strata)!=unique(samp_df$n_strata),TRUE,FALSE)
      
      #if condition reduce or increase CV to achieve the objective
      if (nrow(solution$aggr_strata)<unique(samp_df$n_strata)) {
        cv_df[,c(2:(n_spp+1))]<-cv_df[,c(2:(n_spp+1))]-(cv_df[,c(2:(n_spp+1))]*0.001) #reduce 0.1% CV
      } else if (nrow(solution$aggr_strata)>unique(samp_df$n_strata)){
        cv_df[,c(2:(n_spp+1))]<-cv_df[,c(2:(n_spp+1))]+cv_df[,c(2:(n_spp+1))]*0.01} #increase 1%
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
    
    sum_stats <- SamplingStrata::summaryStrata(solution$framenew,
                                               solution$aggr_strata,
                                               progress=FALSE)
    sum_stats$stratum_id <- 1:nrow(sum_stats)
    sum_stats$wh <- sum_stats$Allocation / sum_stats$Population
    sum_stats$Wh <- sum_stats$Population / nrow(frame)
    sum_stats <- cbind(sum_stats,
                       subset(x = solution$aggr_strata,
                              select = -c(STRATO, N, COST, CENS, DOM1, X1)))
    
    plot_solution <- solution$indices$X1
    
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
    
    # #append stat results  
    # sp_sum_stats<-rbind(sp_sum_stats,sum_stats)
    
    #store results
    result_list <- list(solution = solution,
                        sum_stats = sum_stats,
                        cvs = as.numeric(calc_expected_CV(sum_stats)),
                        n = sum(sum_stats$Allocation),
                        sol_by_cell = plot_solution)
    save(list = "result_list", file = "./output/ms_optim_strata_result_list_',",samp_df[s,'samp_scn'],".RData")
    
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##   7) Single-Species Optimization ----
    ##   Calculate single-species CV subject to the initial stratification
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

    ss_sample_allocations <- expand.grid(n = n_samples, species = spp)
    for (ispp in 1:n_spp) {
      #ispp<-5
      
      temp_n <- result_list$n
      
      #if one stratifying factor add columns  
      if (!grepl('_',samp_df[s,'strat_var'])) {
      ## Subset density data for species ispp
      ss_df <- subset(x = frame, 
                      select = c("domainvalue", "id", "WEIGHT", "X1", #"X2", 
                                 paste0("Y", ispp), paste0("Y", ispp, "_SQ_SUM")))
      } else {
        ss_df <- subset(x = frame, 
                        select = c("domainvalue", "id", "WEIGHT", "X1", "X2", 
                                   paste0("Y", ispp), paste0("Y", ispp, "_SQ_SUM"))) 
      }
      names(ss_df)[grep(x = names(ss_df), pattern = "Y")] <- c("Y1", "Y1_SQ_SUM")
      
      ## Create CV inputs to the Bethel algorithm; initialize at SRS CV
        error_df <- data.frame("DOM" = "DOM1",
                               as.numeric(srs_cv)[ispp],
                               "domainvalue"  = 1)
        names(error_df)[2] <- "CV1"
        
        ## subset stratum stats for the species of interest as inputs to the 
        ## Bethel algorithm
        #if one stratifying factor add columns  
        if (!grepl('_',samp_df[s,'strat_var'])) {
        temp_stratif <- 
          solution$aggr_strata[, c("STRATO", "N", 
                                   paste0("M", ispp), paste0("S", ispp), 
                                   "COST", "CENS", "DOM1", "X1" , "SOLUZ")]
        } else {
          temp_stratif <- 
            solution$aggr_strata[, c("STRATO", "N", 
                                     paste0("M", ispp), paste0("S", ispp), 
                                     "COST", "CENS", "DOM1", "X1" , "X2","SOLUZ")]
        } 
        temp_stratif$N <- temp_stratif$N / n_years
        temp_stratif$DOM1 <- 1
        names(temp_stratif)[3:4] <- paste0(c("M", "S"), 1)
        
        ## run bethel at the SRS CV
        temp_bethel <- SamplingStrata::bethel(
          errors = error_df,
          stratif = temp_stratif, 
          realAllocation = T, 
          printa = T)
        
        ## Save the current n and cv constraint
        temp_n <- sum(ceiling(temp_bethel))
        updated_cv_constraint <- 
          as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "])
        
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
          updated_cv_constraint <- 
            as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "])
          print(paste0("n = ", temp_n, ", ", updated_cv_constraint) )
        }
        
        ## Save the CV and station allocations that corresponds to n_samples
        temp_idx <- ss_sample_allocations$n == n_samples & 
          ss_sample_allocations$species == spp[ispp]
        
        ss_sample_allocations[temp_idx, "CV"] <- 
          as.numeric(attributes(temp_bethel)$outcv[, "ACTUAL CV"])
        
        ss_sample_allocations[temp_idx, paste0("Str_", 1:length(temp_bethel))] <- 
          as.integer(temp_bethel)
        
      
    }
    
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##   8) Adjust MS solution ----
    ##   Optimize allocation across a range of sample sizes, given the original
    ##   stratification.
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
    ms_sample_allocations <- expand.grid(n = n_samples)
    #temp_n <- result_list$n
    
    ## Subset lower limits of CVs from the ss cvs
    ss_cvs <- subset(ss_sample_allocations, n == n_samples)$CV
      
    ## CV dataframe input to Bethel algorithm
    error_df <-  data.frame("DOM" = "DOM1",
                            srs_cv,
                            "domainvalue"  = 1)
    names(error_df)[2:(1 + n_spp)] <- paste0("CV", 1:n_spp)
      
      ## Stratum statistics input to the Bethel algorithm
      temp_stratif <- solution$aggr_strata
      temp_stratif$N <- temp_stratif$N / n_years
      temp_stratif$DOM1 <- 1
      
      ## Run Bethel algorithm and save current n and cv constraints
      temp_bethel <- SamplingStrata::bethel(
        errors = error_df,
        stratif = temp_stratif, 
        realAllocation = T, 
        printa = T)
      
      temp_n <- sum(ceiling(temp_bethel))
      updated_cv_constraint <- 
        as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "])
      
      ## Rerun Bethel algorithm, modifying the CVs relative to the distances
      ## between the SRS and SS CVs given n_samples stations. First we calculate 
      ## CVs calculated under SRS for each species with n_samples stations
      temp_srs_var <- 
        srs_stats[, paste0("S", 1:n_spp)]^2 * (1 - n_samples / n_cells) / n_samples
      temp_srs_cv <- sqrt(temp_srs_var) / srs_stats[, paste0("M", 1:n_spp)]
      
      while (temp_n != n_samples) {
        over_under <- temp_n > n_samples
        
        ## If the current n is < n_samples, decrease the CV by a small amount
        ## relative to the distance between the current CV and the SS CV
        if (over_under == FALSE) {
          # CV_adj <- 0.95
          # updated_cv_constraint <- 
          #   updated_cv_constraint * (CV_adj) + ss_cvs * (1  - CV_adj)
          
          updated_cv_constraint <- updated_cv_constraint  - updated_cv_constraint * 0.001
        }
        
        ## If the current n is > n_samples, increase the CV by a small amount
        ## relative to the distance between the current CV and the SRS CV
        if(over_under == TRUE) {
          # CV_adj = .05
          # updated_cv_constraint <- 
          #   temp_srs_cv * (CV_adj) + updated_cv_constraint * (1  - CV_adj)
          
          updated_cv_constraint <- updated_cv_constraint  + updated_cv_constraint * 0.01
        }
        
        ## Update the CV dataframe input with the updated_cv_constraint
        error_df[, paste0("CV", 1:n_spp)] <- updated_cv_constraint
        
        ## Rerun Bethel algorithm
        temp_bethel <- SamplingStrata::bethel(stratif = temp_stratif,
                                              errors = error_df, 
                                              printa = TRUE)
        
        ## Save sample size and CV constraint
        temp_n <- sum(as.numeric(temp_bethel))
        updated_cv_constraint <- 
          as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "])
        
        ## Print out result to console
        print(paste0("n = ", temp_n) )
      }
      
      ## Save optimized CV 
      temp_idx <- ms_sample_allocations$n == n_samples
      
      ms_sample_allocations[temp_idx, paste0("CV", 1:n_spp)] <- 
        as.numeric(attributes(temp_bethel)$outcv[, "ACTUAL CV"])
      
      ms_sample_allocations[temp_idx, paste0("Str_", 1:length(temp_bethel))] <- 
        as.integer(temp_bethel)
      
    

    #store bethel CV
    cv_bethel_final<-error_df
    
    #cv
    cv_temp <- data.frame(CV_random=cv_initial, #max
                          CV_strata=cv_strata_final, #min to achieve strata
                          CV_bethel=cv_bethel_final) #min to achieve naximum sample
    
    save(list = c('ss_sample_allocations','ms_sample_allocations','cv_temp'),
         file = paste0("./output/ms_optim_allocations",samp_df[s,'samp_scn'],".RData"))
    
}  
