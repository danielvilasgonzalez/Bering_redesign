####################################################################
####################################################################
##
##    get raw slope EBS data
##    fit single VAST models
##    4 models for each sp (null, depth, sbt and depth+sbt)
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu)
##
####################################################################
####################################################################

#load libraries
library(googledrive)
library(dplyr)
library(ggplot2)
library(sp) 
library(sf)
library(raster)
library(TMB)
library(VAST)
library(effects)

#setwd
setwd('E:/UW/Adapting Monitoring to a Changing Seascape/')

#get files from google drive and set up
files<-googledrive::drive_find()
1 #depending on your email

#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='Bering redesign RWP project'),'id']

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)
id.data<-files.1[which(files.1$name=='data raw'),'id']
files.2<-googledrive::drive_ls(id.data$id)

#####################################
# SP CODE DATA
#####################################
#get species file
files.bering<-files.2$name #"metadata"           "ebs_slope_cpue.csv" "ebs_shelf_cpue.csv" "nbs_cpue.csv"       "species.csv" 
file<-files.bering[grep('species',files.bering)]
file.id<-files.2[which(files.2$name==file),'id']

#download file into temp folder
googledrive::drive_download(file=file.id$id,
                            path = paste0('./Data/',file),
                            overwrite = TRUE)

#read csv file
df.sp<-read.csv(paste0('./Data/',file))
df.sp1<-df.sp[,c('SPECIES_CODE','SPECIES_NAME')]

#sp list
sp.list<-c('Limanda aspera','Gadus chalcogrammus','Gadus macrocephalus','Atheresthes stomias','Reinhardtius hippoglossoides',
           'Lepidopsetta polyxystra','Hippoglossoides elassodon','Pleuronectes quadrituberculatus','Hippoglossoides robustus')

#####################################
# EBS SLOPE CPUE DATA
#####################################
#get cpue file
file<-files.bering[grep('slope',files.bering)]
file.id<-files.2[which(files.2$name==file),'id']

#create temp file to store file
#temp <- tempfile(fileext = ".zip")
temp<-tempdir()

#download file into temp folder
googledrive::drive_download(file=file.id$id,
                            path = paste0('./Data/',file),
                            overwrite = TRUE)

#read csv file
df.slope<-read.csv(paste0('./Data/',file))

#get species name by merging both df
df.slope1<-merge(x = df.slope,
                 y = df.sp1,
                 by = 'SPECIES_CODE',
                 all.x = TRUE)

#rename cols df slope
names(df.slope1)<-c('Sp_code','Year','Hauljoin','Strata','Lat','Lon','SST','SBT','Depth','CPUE_kg','CPUE_n','Sp_name')

#correct slope data (slope is in ha, while others and in km2)
df.slope1$CPUE_kg<-df.slope1$CPUE_kg/100

#load grids
grid.slope<-read.csv('./Resources/Extrapolation Grids/SlopeThorsonGrid.csv')
colnames(grid.slope)[10]<-'Area_km2'

##############################################
# fit VAST singlesp models
##############################################

#VAST version to use (code from ZO)
vast_cpp_version <- "VAST_v12_0_0"

#code from ZO
# n_g <- nrow(grid.slope) #number of grid cells
# n_t <- diff(range(df.slope1$YEAR)) + 1 #Number of total years
# n_p <- 2 #two density covariates
# 
# grid.slope$LOG_DEPTH<-log(grid.slope$ARDEMdepth)
# grid.slope$LOG_Q_DEPTH<-log((grid.slope$ARDEMdepth)^2)
# 
# X_gtp <- array(dim = c(n_g, n_t, n_p) )
# for (i in 1:n_t) {
#   X_gtp[, i, ] <- as.matrix(grid.slope[, c("LOG_DEPTH","LOG_Q_DEPTH")])
# }

#models to run: NULL, DEPTH, SBT, and DEPTH+SBT
models<-c('null','depth','sbt','depth_sbt')

#create folder
dir.create('./slope EBS VAST',showWarnings = FALSE)

#species list
sp.list1<-sp.list[sp.list %in% df.slope1$Sp_name]

#loop over species
for (sp in sp.list1) {
  #sp<-sp.list1[5]
  
  #windows progress bar
  py <- winProgressBar(title = paste0(sp, ' (',which(sp.list1 == sp),' out of ',length(sp.list1),')'), # Window title
                       label = "Percentage completed", # Window label
                       min = 0,      # Minimum value of the bar
                       max = 4, #(models) #length(sp.list1), # Maximum value of the bar
                       initial = 0,  # Initial value of the bar
                       width = 300L) # Width of the window   
  
  #create folder to store results
  dir.create(paste0(getwd(),'/slope EBS VAST/',sp),
             showWarnings = FALSE)
  
  #create df to store diagnostics
  df.diagnostics<-data.frame(row.names = models)
  
  #filter by sp
  df<-subset(df.slope1, Sp_name == sp)
  
  #remove rows with NAs in env data
  df<-df[complete.cases(df[c('SBT','Depth')]),]
  
  #covariate data
  covariate_data<-df[,c("Lat","Lon","Year",'CPUE_kg',"Depth",'SBT')]
  
  #VAST model settings
  settings <- make_settings(n_x=100, 
                            Region="bering_sea_slope", #'User'
                            purpose="index2", 
                            bias.correct=FALSE,
                            knot_method='grid',
                            use_anisotropy=TRUE,
                            Version = vast_cpp_version,
                            fine_scale=TRUE,
                            ObsModel = c(2,1), #c(2,1) #biomass
                            max_cells = Inf,
                            Options = c("Calculate_Range" = F, 
                                        "Calculate_effective_area" = F))
  
  #each model
  for (m in models) {
    #m<-models[3]
    
    #create model folder to store results
    dir.create(paste(getwd(),'slope EBS VAST',sp,m,sep='/'),
               showWarnings = FALSE)
    
    if (m=='null') {
     covariate_data1<- cbind(covariate_data[,c("Lat","Lon")],Year=NA)
    } else if (m=='depth'){
      covariate_data1<- cbind(covariate_data[,c("Lat","Lon","Depth")],Year=NA)
    } else if (m=='sbt'){
      covariate_data1<- cbind(covariate_data[,c("Lat","Lon","SBT")],Year=NA) #covariate_data[,c("Lat","Lon","SBT","Year")],
    } else if (m=='depth_sbt'){
      covariate_data1<- cbind(covariate_data[,c("Lat","Lon","Depth",'SBT')],Year=NA) #covariate_data[,c("Lat","Lon","SBT","Year")],
    }

    #fit
    tryCatch({
    fit <- fit_model(settings=settings,
                     Lat_i=df$Lat, 
                     Lon_i=df$Lon,
                     t_i=df$Year,
                     b_i=df$CPUE_kg,
                     c_iz = as.numeric(factor(df$Sp_name))-1,
                     a_i=rep(1,times=nrow(df)),
                     X1_formula = ifelse(m=='null','~0',
                                         ifelse(m=='depth','~log(Depth)+(log(Depth))^2',
                                                ifelse(m=='sbt','~log(SBT)+(log(SBT))^2',
                                                       '~ log(Depth)+(log(Depth))^2+log(SBT)+(log(SBT))^2'))),
                     X2_formula = ifelse(m=='null','~0',
                                         ifelse(m=='depth','~log(Depth)+(log(Depth))^2',
                                                ifelse(m=='sbt','~log(SBT)+(log(SBT))^2',
                                                       '~ log(Depth)+(log(Depth))^2+log(SBT)+(log(SBT))^2'))),
                     covariate_data = covariate_data1, #covariate_data[,c("Lat","Lon","SBT","Year")],
                     #input_grid=grid.slope,
                     #strata.limits=data.frame(STRATA = c("All_areas","slope","shelf"))
                     #getJointPrecision = TRUE,
                     #newtonsteps = 1,
                     category_name=sp,
                     test_fit=FALSE,
                     working_dir = paste(getwd(),'slope EBS VAST',sp,m,'/',sep='/'))#;beepr::beep(sound = 8)
    }, error=function(e){})
    
    try({
      if (class(fit$Report)!='list') {
       
        settings2<-settings
        settings2$ObsModel<-c(1,1)
        fit <- fit_model(settings=settings2,
                         Lat_i=df$Lat, 
                         Lon_i=df$Lon,
                         t_i=df$Year,
                         b_i=df$CPUE_kg,
                         c_iz = as.numeric(factor(df$Sp_name))-1,
                         a_i=rep(1,times=nrow(df)),
                         X1_formula = ifelse(m=='null','~0',
                                             ifelse(m=='depth','~log(Depth)+(log(Depth))^2',
                                                    ifelse(m=='sbt','~log(SBT)+(log(SBT))^2',
                                                           '~ log(Depth)+(log(Depth))^2+log(SBT)+(log(SBT))^2'))),
                         X2_formula = ifelse(m=='null','~0',
                                             ifelse(m=='depth','~log(Depth)+(log(Depth))^2',
                                                    ifelse(m=='sbt','~log(SBT)+(log(SBT))^2',
                                                           '~ log(Depth)+(log(Depth))^2+log(SBT)+(log(SBT))^2'))),
                         covariate_data = covariate_data1, #covariate_data[,c("Lat","Lon","SBT","Year")],
                         #input_grid=grid.slope,
                         #strata.limits=data.frame(STRATA = c("All_areas","slope","shelf"))
                         #getJointPrecision = TRUE,
                         #newtonsteps = 1,
                         category_name=sp,
                         test_fit=FALSE,
                         working_dir = paste(getwd(),'slope EBS VAST',sp,m,'/',sep='/'))#;beepr::beep(sound = 8)
      }
    })
  
    #save fit
    save(list = "fit", file = paste(getwd(),'slope EBS VAST',sp,m,'fit.RData',sep='/'))
    load(paste(getwd(),'slope EBS VAST',sp,m,'fit.RData',sep='/'))
    #plot
    if (class(fit$Report) =='list') {
    plot(fit,
         #plot_set = c(1:21),
         working_dir =  paste(getwd(),'slope EBS VAST',sp,m,'/',sep='/'))
    }
    
    #diagnostics
    if (class(fit$Report)!='list') {
      df.diagnostics[m,'ObsModel1']<-fit$data_list$ObsModel_ez[1,1]
      df.diagnostics[m,'ObsModel2']<-fit$data_list$ObsModel_ez[1,2]
      df.diagnostics[m,'Covergence']<-fit$Report
      df.diagnostics[m,c('AIC','max_gradient','deltaAIC','jnll','rmse','mae','depth_effect1','depth_effect2','sbt_effect1','sbt_effect2')] <-NA
      
    } else{
      df.diagnostics[m,'ObsModel1']<-fit$data_list$ObsModel_ez[1,1]
      df.diagnostics[m,'ObsModel2']<-fit$data_list$ObsModel_ez[1,2]
      df.diagnostics[m,'Covergence']<-fit$parameter_estimates$Convergence_check
      df.diagnostics[m,'AIC']<-fit$parameter_estimates$AIC
      df.diagnostics[m,'max_gradient']<-fit$parameter_estimates$max_gradient
      df.diagnostics[m,'jnll']<-fit$parameter_estimates$objective
      df.diagnostics[m,'rmse']<-sqrt(mean((fit$data_frame$b_i - fit$Report$D_i)^2))
      
      sum = 0      
      for (i in 1:length(fit$data_frame$b_i)){
        sum <- abs(fit$data_frame$b_i[i] - fit$Report$D_i[i]) + sum}
      
      df.diagnostics[m,'mae'] <- sum/length(fit$data_frame$b_i)
      
        if (m=='null') {
          df.diagnostics[m,c('depth_effect1','depth_effect2','sbt_effect1','sbt_effect2')]<-NA
        } else if (m %in% c('depth','sbt')) {
          df.diagnostics[m,paste0(m,'_effect1')]<-fit$ParHat$gamma1_cp
          df.diagnostics[m,paste0(m,'_effect2')]<-fit$ParHat$gamma2_cp
        } else if (m == 'depth_sbt'){
          df.diagnostics[m,'depth_effect1']<-fit$ParHat$gamma1_cp[1]
          df.diagnostics[m,'depth_effect2']<-fit$ParHat$gamma2_cp[1]
          df.diagnostics[m,'sbt_effect1']<-fit$ParHat$gamma1_cp[2]
          df.diagnostics[m,'sbt_effect2']<-fit$ParHat$gamma2_cp[2]
        }
      }
      
    #remove fit object
    rm(fit)
    
    #progress bar
    pctgy <- paste0(round(which(models == m)/length(models) *100, 0), "% completed")
    setWinProgressBar(py, which(models == m), label = pctgy) # The label will override the label set on the
    
    }    

  #add sp column
  df.diagnostics$sp<-sp
  
  #calculate deltaAIC
  for (m in models) {
    df.diagnostics[m,'deltaAIC']<-df.diagnostics['null','AIC']-df.diagnostics[m,'AIC']
  }

  #write df diagnostics
  write.csv(df.diagnostics, 
            file = paste(getwd(),'slope EBS VAST',sp,"table_diagnostics.csv",sep='/'),
            row.names = TRUE)
  
  #close window
  close(py)
  }
    
    
##############################
