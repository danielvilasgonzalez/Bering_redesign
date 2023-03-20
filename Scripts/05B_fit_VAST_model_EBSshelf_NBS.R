####################################################################
####################################################################
##    
##    fit single sp VAST model for the EBS shelf and NBS using depth and temp (SBT)
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu)
##    evaluate parameters : https://github.com/James-Thorson-NOAA/VAST/blob/main/R/make_parameters.R
##    https://rdrr.io/github/James-Thorson/VAST/man/check_fit.html
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c("splines")

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#install coldpool to extract SBT for the EBS
if (!('VAST' %in% installed.packages())) {
  devtools::install_github("james-thorson/VAST@main", INSTALL_opts="--no-staged-install")};library(VAST)

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#setwd
out_dir<-'C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

#version VAST (cpp)
version<-'VAST_v13_1_0'

#number of knots
knots<-'500' #100

#years
yrs<-1982:2015

#list of sp
splist<-list.dirs('./data processed/',full.names = FALSE,recursive = FALSE)

#folder region
fol_region<-'shelf EBS NBS VAST'

#dir create for slope region results
dir.create(paste(out_dir,fol_region,sep='/'))

#dir create for splist
for (sp in splist) {
  dir.create(paste(out_dir,fol_region,sp,sep='/'))
}

#list of models
# models<-c('null',
#           as.vector(outer(c('depth','temp'), c('2d','3d'), paste, sep="")),
#           as.vector(outer(as.vector(outer(c('depth','temp'), c('2d','3d'), paste, sep=""))[c(1,3)],
#                           as.vector(outer(c('depth','temp'), c('2d','3d'), paste, sep=""))[c(2,4)],
#                           paste, sep="_")))

models<-c('null','temp2d','temp3d')

#diagnostics df
diagnostics<-array(dim = c(length(models),7,length(splist)),
                   dimnames = list(models,c('status','maxgradient','aic','jnll','rmse','dev','pct.dev'),splist))

# effects_df<-array(dim = c(length(models),12,length(splist)),
#                   dimnames = list(models,c('pres_depth1','pres_depth2','pres_depth3',
#                                            'pos_depth1','pos_depth2','pos_depth3',
#                                            'pres_temp1','pres_temp2','pres_temp3',
#                                            'pos_temp1','pos_temp2','pos_temp3'),splist))

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
df1<-readRDS(paste0('./data processed/',sp,'/data_geostat_temp.rds'))
#df1[which(df1$year==2020),'bottom_temp_c']<-NA
df2<-subset(df1,year %in% yrs)

#select rows and rename
df3<-df2[,c("lat_start","lon_start","year",'scientific_name','weight_kg','effort','depth_m','LogDepth',"ScaleLogDepth",'Scalebottom_temp_c','bottom_temp_c','survey_name')]
colnames(df3)<-c('Lat','Lon','Year','Species','CPUE_kg','Effort','Depth','LogDepth','ScaleLogDepth','ScaleBotTemp','BotTemp','Region')

#data geostat
df4<-subset(df3,Region %in% c("Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey",
                              "Northern Bering Sea Crab/Groundfish Survey - Eastern Bering Sea Shelf Survey Extension"))

data_geostat<-df4[complete.cases(df4[,c('CPUE_kg')]),]

#covariate data - filter by year and complete cases for env variables
#covariate_data<-subset(df2,Year>=yrs_region[1] & Year<=yrs_region[2])
covariate_data<-df3[complete.cases(df3[,c('BotTemp')]),] #,'ScaleLogDepth'

#save data
saveRDS(data_geostat,paste(out_dir,fol_region,sp,'data_geostat_temp.rds',sep='/'))

#regions (predefined in VAST)
region<-c("northern_bering_sea","eastern_bering_sea")

#loop over models
#for (m in models) {
  
  m<-models[3]
  
  
  #print year to check progress
  cat(paste("\n","    ----- ", sp, " -----\n","       - ", m, " model\n"))  
  
  #create folder to store results
  dir.create(paste(out_dir,fol_region,sp,m,sep='/'),
             showWarnings = FALSE)
  
  #VAST model settings
  settings <- make_settings(n_x=knots, 
                            Region=region, #c("bering_sea_slope","eastern_bering_sea",'northern_bering_sea'
                            purpose="index2", 
                            bias.correct=FALSE,
                            knot_method='grid',
                            use_anisotropy=TRUE,
                            #FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID"),
                            RhoConfig=c("Beta1"=2,"Beta2"=2,"Epsilon1"=4,"Epsilon2"=4), 
                            Version = version,
                            #fine_scale=TRUE,
                            ObsModel = c(2,1), #c(2,3) if encounter 100% year
                            max_cells = Inf,
                            Options = c("Calculate_Range" =  F, 
                                        "Calculate_effective_area" = F)) 
  
  
  # if (m=='null') {
  #   
  #   # Make settings (turning off bias.correct to save time for example)
  #   settings$FieldConfig = matrix( c(0,0,0,0,"IID","IID"), byrow=TRUE, ncol=2 )
  #   settings$RhoConfig[c("Beta1","Beta2")] = 3
  #   
  # }
  
  #Kmeans_knots-200
  if (!file.exists(paste(out_dir,fol_region,sp,m,'Kmeans_knots-',knots,'.RData',sep='/')) & m!=models[1]) {
    file.copy(paste(out_dir,fol_region,sp,models[1],'Kmeans_knots-',knots,'.RData',sep='/'),
              paste(out_dir,fol_region,sp,m,'Kmeans_knots-',knots,'.RData',sep='/'))}
  
  #formula and predictors settings for each model
  if(grepl('depth',m) & grepl('temp',m)){
    
    f1<-ifelse(grepl('depth2d',m),
               ' ~ bs(ScaleLogDepth, degree=2, intercept=FALSE)',
               ' ~ bs(ScaleLogDepth, degree=3, intercept=FALSE)')
    f2<-ifelse(grepl('temp2d',m),
               ' + bs(ScaleBotTemp, degree=2, intercept=FALSE)',
               ' + bs(ScaleBotTemp, degree=3, intercept=FALSE)')
    
    X1_formula<-paste(f1,f2)
    
    X1config_cp = array( c(1,1), dim=c(1,2))
    
  } else if (grepl('depth',m)){
    
    X1_formula<-ifelse(grepl('depth2d',m),
                       ' ~ bs(ScaleLogDepth, degree=2, intercept=FALSE)',
                       ' ~ bs(ScaleLogDepth, degree=3, intercept=FALSE)')
    
    X1config_cp = array( c(1,1), dim=c(1,1))
    
  } else if (grepl('temp',m)) {
    
    X1_formula<-ifelse(grepl('temp2d',m),
                       ' ~ bs(BotTemp, degree=2, intercept=FALSE)',
                       ' ~ bs(BotTemp, degree=3, intercept=FALSE)')
    
    X1config_cp = array( c(1,1), dim=c(1,1))
    
  } else {
    
    X1_formula<-' ~ 0'
    
    X1config_cp = NULL
    
  }
  
  #formula for positive catch rates equal to presence/absence
  X2_formula<-X1_formula
  
  #predictor settings
  X2config_cp = X1config_cp
  
  #modify settings to use 2 or 3 factors on the spatial and spatiotemporal variation
  # if (grepl('IID',m)) {
  #   settings$FieldConfig[c('Epsilon','Omega'),]<-'IID'
  # } else if (grepl('f2',m)) {
  #   settings$FieldConfig[c('Epsilon','Omega'),]<-2
  # } else if (grepl('f3',m)) {
  #   settings$FieldConfig[c('Epsilon','Omega'),]<-3
  # }
  
  #fit model #### ADD TryCatch{(),}
  fit <- tryCatch( {fit_model(settings=settings,
                              Lat_i=data_geostat$Lat, 
                              Lon_i=data_geostat$Lon,
                              t_i=data_geostat$Year,
                              b_i=data_geostat$CPUE_kg,
                              c_iz = as.numeric(factor(data_geostat$Species))-1,
                              a_i=data_geostat$Effort/100,
                              #input_grid=grid.ebs,
                              getJointPrecision = TRUE,
                              test_fit=FALSE,
                              create_strata_per_region = TRUE,  
                              covariate_data = covariate_data[,c('Year',"Lat","Lon","ScaleLogDepth","LogDepth",'ScaleBotTemp','BotTemp',"CPUE_kg")], 
                              X1_formula =  X1_formula,
                              X2_formula = X2_formula, 
                              newtonsteps = 1,
                              #X_gtp = X_gtp,
                              working_dir = paste(out_dir,fol_region,sp,m,'/',sep='/'))},
                   error = function(cond) {
                     message("Did not converge. Here's the original error message:")
                     message(cond)
                     cat(paste(" ERROR MODEL",m,"IS NOT CONVERGING "))  
                     # Choose a return value in case of error
                     return(NULL)
                   })
  
  #plot(fit,
  #     working_dir=paste(out_dir,fol_region,sp,m,'/',sep='/'))
  
  #save fit
  save(list = "fit", file = paste(out_dir,fol_region,sp,m,'b2_19822015fit.RData',sep='/')) #paste(yrs_region,collapse = "")
  
  #convergence
  diagnostics[m,'status',sp]<-ifelse(test = is.null(fit) == T | is.null(fit$parameter_estimates$max_gradient),"no_convergence","check_gradient")
  #max gradient
  
  if (diagnostics[m,'status',sp]=="check_gradient") {
    
    diagnostics[m,'maxgradient',sp]<-fit$parameter_estimates$max_gradient
    #AIC
    diagnostics[m,'aic',sp]<-round(fit$parameter_estimates$AIC[1],3)
    #JNLL
    diagnostics[m,'jnll',sp]<-round(fit$parameter_estimates$objective[1],3)
    #RMSE
    diagnostics[m,'rmse',sp]<-round(sqrt(mean((data_geostat$CPUE_kg - fit$Report$D_i)^2,na.rm = TRUE)) / mean(data_geostat$CPUE_kg,na.rm = TRUE),3)
    #DEVIANCE EXPLAINED
    diagnostics[m,'dev',sp]<-round(fit$Report$deviance,3)
    
    #depth effects
    if (grepl('depth',m)) {
      effects_df[m,c('pres_depth1','pres_depth2'),sp]<-round(fit$ParHat$gamma1_cp[ ,grepl( "ScaleLogDepth",names(data.frame(fit$ParHat$gamma1_cp)))][c(1:2)],3)
      effects_df[m,c('pos_depth1','pos_depth2'),sp]<-round(fit$ParHat$gamma2_cp[ ,grepl( "ScaleLogDepth",names(data.frame(fit$ParHat$gamma2_cp)))][c(1:2)],3)
      #if 3 degrees
      if (grepl(3,m)) {
        effects_df[m,'pres_depth3',sp]<-round(fit$ParHat$gamma1_cp[ ,grepl( "ScaleLogDepth",names(data.frame(fit$ParHat$gamma1_cp)))][3],3)
        effects_df[m,'pos_depth3',sp]<-round(fit$ParHat$gamma2_cp[ ,grepl( "ScaleLogDepth",names(data.frame(fit$ParHat$gamma2_cp)))][3],3)
      }
    }
    
    #temp effects
    if (grepl('temp',m)) {
      effects_df[m,c('pres_temp1','pres_temp2'),sp]<-round(fit$ParHat$gamma1_cp[ ,grepl( "ScaleBotTemp",names(data.frame(fit$ParHat$gamma1_cp)))][c(1:2)],3)
      effects_df[m,c('pos_temp1','pos_temp2'),sp]<-round(fit$ParHat$gamma2_cp[ ,grepl( "ScaleBotTemp",names(data.frame(fit$ParHat$gamma2_cp)))][c(1:2)],3)
      #if 3 degrees
      if (grepl(3,m)) {
        effects_df[m,'pres_temp3',sp]<-round(fit$ParHat$gamma1_cp[ ,grepl( "ScaleBotTemp",names(data.frame(fit$ParHat$gamma1_cp)))][3],3)
        effects_df[m,'pos_temp3',sp]<-round(fit$ParHat$gamma2_cp[ ,grepl( "ScaleBotTemp",names(data.frame(fit$ParHat$gamma2_cp)))][3],3)
      }
    }
  }
  
  #progress bar
  pctgy <- paste0(round(which(models == m)/length(models) *100, 0), "% completed")
  setWinProgressBar(py, which(models == m), label = pctgy) # The label will override the label set on the
  
  
#}


#percent of dev explained
diagnostics[,'pct.dev',sp]<-(1 - as.numeric(diagnostics[,'dev',sp])/as.numeric(diagnostics['null','dev',sp]))*100

#df diagnostics
df.diagnostics<-data.frame(diagnostics[,,sp])
df.diagnostics[,'pct.dev']<-(1 - as.numeric(df.diagnostics[,'dev'])/as.numeric(df.diagnostics['null','dev']))*100
#df.diagnostics$pct.dev<-(1 - as.numeric(diagnostics[,'dev',sp])/as.numeric(diagnostics['null','dev',sp]))*100

#df effects
df.effects<-data.frame(effects_df[,,sp])

#save RDS effects and diagnostics - move out of the loop when end testing pcod
saveRDS(data.frame(effects_df[,,sp]),paste(out_dir,fol_region,sp,'effects.RData',sep='/'))
saveRDS(data.frame(df.diagnostics[,,sp]),paste(out_dir,fol_region,sp,'diagnostics.RData',sep='/'))

#save csv file
write.csv(df.diagnostics,paste(out_dir,fol_region,sp,'diagnostics.csv',sep='/'))

#close process window
close(py)

#}

##################

cf<-check_fit(fit$parameter_estimates)
data.frame(fit$parameter_estimates$par)


