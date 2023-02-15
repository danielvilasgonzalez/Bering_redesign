####################################################################
####################################################################
##    
##    fit single sp VAST models using depth and temp (SBT)
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu)
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c("splines",'effects')

#install pacman to use p_load function - call library and if not installed, then install
 if (!('pacman' %in% installed.packages())) {
   install.packages("pacman")}

#install coldpool to extract SBT for the EBS
if (!('VAST' %in% installed.packages())) {
  devtools::install_github("james-thorson/VAST@main", INSTALL_opts="--no-staged-install")};library(VAST)

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#setwd
out_dir<-'E:/UW/Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

#number of knots
knots<-'200'

#list of sp
splist<-list.dirs('./slope shelf EBS NBS VAST',full.names = FALSE,recursive = FALSE)
splist<-sort(splist[-1])

#list of models
models<-c('null',
          as.vector(outer(c('depth','temp'), c('2d','3d'), paste, sep="")),
          as.vector(outer(as.vector(outer(c('depth','temp'), c('2d','3d'), paste, sep=""))[c(1,3)],
                          as.vector(outer(c('depth','temp'), c('2d','3d'), paste, sep=""))[c(2,4)],
                          paste, sep="_")))
#models<-c('null','depth','temp','full')

#diagnostics df
diagnostics<-array(dim = c(length(models),9,length(splist)),
                   dimnames = list(models,c('status','maxgradient','aic','jnll','rmse','ScaleLogDepth1','ScaleLogDepth2','ScaleTemp1','ScaleTemp2'),splist))

#loop over species to fit models
#for (sp in sp.list) {

#Pcod example
sp<-'Gadus macrocephalus'

#check % of process  
#windows progress bar
py <- winProgressBar(title = paste0(sp, ' (',which(splist == sp),' out of ',length(splist),')'), # Window title
                     label = "Percentage completed", # Window label
                     min = 0,      # Minimum value of the bar
                     max = length(models), #(models) # Maximum value of the bar
                     initial = 0,  # Initial value of the bar
                     width = 300L) # Width of the window   

#read data_geostat_temp file
df1<-readRDS(paste0('./slope shelf EBS NBS VAST/',sp,'/data_geostat_temp.rds'))

#remove rows with NAs
df2<-df1[complete.cases(df1),]

#covariate data
covariate_data<-df1[,c("Lat","Lon","Year",'CPUE_kg',"ScaleLogDepth",'ScaleTemp')]

#regions (predefined in VAST)
region<-c("bering_sea_slope","eastern_bering_sea",'northern_bering_sea')

  #loop over models
  for (m in models) {
  
  #m<-models[2]
  
  #print year to check progress
  cat(paste("\n","    ----- ", sp, " -----\n","       - ", m, " model\n"))  
  
  #create folder to store results
  dir.create(paste0(getwd(),'/slope shelf EBS NBS VAST/',sp,'/',m,'/'),
             showWarnings = FALSE)
    
  #VAST model settings
  settings <- make_settings(n_x=knots, 
                            Region=region, #c("bering_sea_slope","eastern_bering_sea",'northern_bering_sea'
                            purpose="index2", 
                            bias.correct=FALSE,
                            knot_method='grid',
                            use_anisotropy=TRUE,
                            #FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID"),
                            #RhoConfig=c("Beta1"=3,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0), #RhoConfig=c("Beta1"=1,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0)
                            #Version = "VAST_v12_0_0",
                            #fine_scale=TRUE,
                            ObsModel = c(2,1), #c(1,1) #biomass
                            max_cells = Inf,
                            Options = c("Calculate_Range" =  T, 
                                        "Calculate_effective_area" = T)) 
  
  #Kmeans_knots-200
  if (!file.exists(paste0('./slope shelf EBS NBS VAST/',sp,'/',m,'/','Kmeans_knots-',knots,'.RData')) & m!=models[1]) {
    file.copy(paste0('./slope shelf EBS NBS VAST/',sp,'/',models[1],'/','Kmeans_knots-',knots,'.RData'),
              paste0('./slope shelf EBS NBS VAST/',sp,'/',m,'/','Kmeans_knots-',knots,'.RData'))}
  
  #formula for each model
  if(grepl('depth',m) & grepl('temp',m)){
    
    f1<-ifelse(grepl('depth2d',m),
               ' ~ bs(ScaleLogDepth, degree=2)',
               ' ~ bs(ScaleLogDepth, degree=3)')
    f2<-ifelse(grepl('temp2d',m),
               ' + bs(ScaleTemp, degree=2)',
               ' + bs(ScaleTemp, degree=3)')
    
    X1_formula<-paste(f1,f2)
  
  } else if (grepl('depth',m)){
    
    X1_formula<-ifelse(grepl('depth2d',m),
                       ' ~ bs(ScaleLogDepth, degree=2)',
                       ' ~ bs(ScaleLogDepth, degree=3)')
    
  } else if (grepl('temp',m)) {
    
    X1_formula<-ifelse(grepl('temp2d',m),
                       ' + bs(ScaleTemp, degree=2)',
                       ' + bs(ScaleTemp, degree=3)')
    
  } else {
    
    X1_formula<-'~0'
    
  }

  #formula for positive catch rates equal to presence/absence
  X2_formula<-X1_formula
  
  #modify settings to use 2 or 3 factors on the spatial and spatiotemporal variation
  # if (grepl('IID',m)) {
  #   settings$FieldConfig[c('Epsilon','Omega'),]<-'IID'
  # } else if (grepl('f2',m)) {
  #   settings$FieldConfig[c('Epsilon','Omega'),]<-2
  # } else if (grepl('f3',m)) {
  #   settings$FieldConfig[c('Epsilon','Omega'),]<-3
  # }
  
  #fit model #### ADD TryCatch{(),}
  fit <- fit_model(settings=settings,
                   Lat_i=df2$Lat, 
                   Lon_i=df2$Lon,
                   t_i=df2$Year,
                   b_i=df2$CPUE_kg,
                   c_iz = as.numeric(factor(df2$Species))-1,
                   a_i=rep(0.1,times=nrow(df2)),
                   #input_grid=grid.ebs,
                   getJointPrecision = TRUE,
                   test_fit=FALSE,
                   create_strata_per_region = TRUE,  
                   covariate_data = cbind(covariate_data[,c("Lat","Lon","ScaleLogDepth",'ScaleTemp','Year')]), 
                   X1_formula =  X1_formula,
                   X2_formula = X2_formula, 
                   newtonsteps = 0,
                   #X_gtp = X_gtp,
                   working_dir = paste0('./slope shelf EBS NBS VAST/',sp,'/',m,'/'))
  
  #save fit
  save(list = "fit", file = paste0('./slope shelf EBS NBS VAST/',sp,'/',m,'/fit.RData'))
  
  #convergence
  diagnostics[m,'status',sp]<-ifelse(test = is.null(fit) == T | is.null(fit$parameter_estimates$max_gradient),"no_convergence","check_gradient")
  #max gradient
  diagnostics[m,'maxgradient',sp]<-fit$parameter_estimates$max_gradient
  #AIC
  diagnostics[m,'aic',sp]<-round(fit$parameter_estimates$AIC[1],3)
  #JNLL
  diagnostics[m,'jnll',sp]<-round(fit$parameter_estimates$objective[1],3)
  #RMSE
  diagnostics[m,'rmse',sp]<-round(sqrt(mean((df2$CPUE_kg - fit$Report$D_i)^2)) / mean(df2$CPUE_kg),3)

  #depth effects
  if (grepl('depth|full',m)) {
    diagnostics[m,'ScaleLogDepth1',sp]<-round(fit$ParHat$gamma1_cp[,'ScaleLogDepth'],3)
    diagnostics[m,'ScaleLogDepth2',sp]<-round(fit$ParHat$gamma2_cp[,'ScaleLogDepth'],3)
  }
  #temp effects
  if (grepl('temp|full',m)) {
    diagnostics[m,'ScaleTemp1',sp]<-round(fit$ParHat$gamma1_cp[,'ScaleTemp'],3)
    diagnostics[m,'ScaleTemp2',sp]<-round(fit$ParHat$gamma2_cp[,'ScaleTemp'],3)
  }
  
  #progress bar
  pctgy <- paste0(round(which(models == m)/length(models) *100, 0), "% completed")
  setWinProgressBar(py, which(models == m), label = pctgy) # The label will override the label set on the

  }

close(py)

}