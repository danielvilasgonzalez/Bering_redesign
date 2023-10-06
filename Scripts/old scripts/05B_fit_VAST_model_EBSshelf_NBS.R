####################################################################
####################################################################
##    
##    fit single sp VAST model for the EBS shelf and NBS using temp (SBT, BotTemp) with cubic effect 
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
pack_cran<-c("splines",'dplyr')

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
out_dir<-'/Users/daniel/Work/Adapting to a Changing Seascape/'
setwd(out_dir)

#version VAST (cpp)
version<-"VAST_v14_0_1" #if using "VAST_v13_1_0" follow covariate values

#number of knots
knots<-'500' #500, 1000 

#years
yrs<-1982:2022

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

#folder region
fol_region<-c('shelf EBS NBS VAST')

#load grid
load('./extrapolation grids/lastversion_grid_EBS.RData')
grid_ebs<-grid.ebs_year[which(grid.ebs_year$region != 'EBSslope' & grid.ebs_year$Year %in% yrs),]
dim(grid_ebs)

#dir create for slope region results
dir.create(paste(out_dir,fol_region,sep='/'))

#dir create for splist
for (sp in spp) {
  dir.create(paste(out_dir,fol_region,sp,sep='/'))
}

#diagnostics df
diagnostics<-array(dim = c(1,6,length(spp)),
                   dimnames = list('BotTemp3d',c('status','maxgradient','aic','jnll','rmse','dev'),spp))

#loop over species to fit models
#for (sp in spp) {

#Pcod example
sp<-'Paralithodes platypus'

#'Limanda aspera',
#'Gadus chalcogrammus',
#'Gadus macrocephalus',
#'Hippoglossoides elassodon'
#''Atheresthes stomias'
#''Reinhardtius hippoglossoides'
#''Lepidopsetta polyxystra' #not convergenced
#'Pleuronectes quadrituberculatus'
#''Boreogadus saida'
#''Eleginus gracilis'
#''Anoplopoma fimbria'
#''Chionoecetes opilio',
#''Paralithodes platypus' #not converged
#''Paralithodes camtschaticus'
#'
#print year to check progress
cat(paste("\n","    ----- ", sp, " -----\n"))  

#check % of process  
#windows progress bar
# py <- winProgressBar(title = paste0(sp, ' (',which(spp == sp),' out of ',length(spp),')'), # Window title
#                      label = "Percentage completed", # Window label
#                      min = 0,      # Minimum value of the bar
#                      max = length(spp), #(species) # Maximum value of the bar
#                      initial = 0,  # Initial value of the bar
#                      width = 300L) # Width of the window   

#read data_geostat_temp file
df1<-readRDS(paste0('./data processed/species/',sp,'/data_geostat_temp.rds'))
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

#add grid to get prediction for simulate data on each cell of the grid (sim$b_i)
grid_df<-data.frame(Lat=grid_ebs$Lat,
                    Lon=grid_ebs$Lon,
                    Year=grid_ebs$Year,
                    Species=rep(sp,times=nrow(grid_ebs)),
                    CPUE_kg=mean(data_geostat$CPUE_kg),
                    Effort=grid_ebs$Area_in_survey_km2,
                    Depth=grid_ebs$Depth,
                    BotTemp=grid_ebs$Temp,
                    Region=grid_ebs$region,
                    stringsAsFactors = T)

#ha to km2
data_geostat$Effort<-data_geostat$Effort/100

#rbind grid and data_geostat to get prediction into grid values when simulating data
data_geostat1<-rbind(data_geostat[,c("Lat","Lon","Year","Species","CPUE_kg","Effort","Depth","BotTemp","Region")],
                     grid_df)

#to get predictions in locations but not influencing fit
pred_TF <- rep(1, nrow(data_geostat1))
pred_TF[1:nrow(data_geostat)] <- 0

#save data
saveRDS(data_geostat1,paste(out_dir,fol_region,sp,'data_geostat_temp.rds',sep='/'))

#regions (predefined in VAST)
region<-c("northern_bering_sea","eastern_bering_sea")

#get percentage of encounters
xall<-data_geostat %>% count(Year)
x0<-subset(data_geostat,CPUE_kg!=0) %>% count(Year)
xpct<-x0[,2]/xall[,2]*100

#any year with 100%encounters or 0%encounters
enc100<-ifelse(100 %in% xpct,TRUE,FALSE)
enc0<-ifelse(0 %in% xpct,TRUE,FALSE)

#set settings based on enc100
# if(enc100==TRUE){
#   obs <- c(2,3)
# }
# if(enc100==FALSE){
#   obs <- c(2,1)
# }

#VAST model settings
settings <- make_settings(n_x=knots, 
                          Region=region, #c("bering_sea_slope","eastern_bering_sea",'northern_bering_sea'
                          purpose="index2", 
                          bias.correct=FALSE,
                          knot_method='grid',
                          use_anisotropy=TRUE,
                          # FieldConfig = matrix( c("IID","IID",0,"Identity", "IID","IID",0,"Identity"), 
                          #                       ncol=2, 
                          #                       nrow=4, 
                          #                       dimnames=list(c("Omega","Epsilon","Beta","Epsilon_year"),c("Component_1","Component_2"))),
                           #FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID",'Beta1'=0,'Beta2'=0),
                          RhoConfig=c("Beta1"=2,"Beta2"=2,"Epsilon1"=2,"Epsilon2"=4), 
                          Version = version,
                          #fine_scale=TRUE,
                          ObsModel = c(2,1), #obs
                          max_cells = Inf,
                          Options = c("Calculate_Range" =  F, 
                                      "Calculate_effective_area" = F)) 

#Kmeans_knots
# if (!file.exists(paste(out_dir,fol_region,sp,'Kmeans_knots-',knots,'.RData',sep='/'))) {
#   file.copy(paste(out_dir,fol_region,sp,models[1],'Kmeans_knots-',knots,'.RData',sep='/'),
#             paste(out_dir,fol_region,sp,'Kmeans_knots-',knots,'.RData',sep='/'))}

#formula and predictors settings for each model
X1_formula<-' ~ bs(BotTemp, degree=3, intercept=FALSE)'
X1config_cp = array( c(1,1), dim=c(1,1))

#formula for positive catch rates equal to presence/absence
X2_formula<-X1_formula

#predictor settings
X2config_cp = X1config_cp

#fit model #### ADD TryCatch{(),}
fit <- tryCatch( {fit_model(settings=settings,
                            Lat_i=data_geostat1$Lat, 
                            Lon_i=data_geostat1$Lon,
                            t_i=data_geostat1$Year,
                            b_i=data_geostat1$CPUE_kg,
                            c_iz = as.numeric(factor(data_geostat1$Species))-1,
                            a_i=data_geostat1$Effort,
                            #input_grid=grid.ebs,
                            getJointPrecision = TRUE,
                            test_fit=FALSE,
                            create_strata_per_region = TRUE,  
                            covariate_data = covariate_data[,c('Year',"Lat","Lon","ScaleLogDepth","LogDepth",'ScaleBotTemp','BotTemp',"CPUE_kg")], 
                            X1_formula =  X1_formula,
                            X2_formula = X2_formula, 
                            newtonsteps = 1,
                            PredTF_i = pred_TF,
                            #X_gtp = X_gtp,
                            working_dir = paste(out_dir,fol_region,sp,'/',sep='/'))},
                 error = function(cond) {
                   message("Did not converge. Here's the original error message:")
                   message(cond)
                   # Choose a return value in case of error
                   return(NULL)
                 })

#check fit
check_fit(fit$parameter_estimates)

#plot(fit,
#     working_dir=paste(out_dir,fol_region,sp,m,'/',sep='/'))

#save fit
save(list = "fit", file = paste(out_dir,fol_region,sp,'fit.RData',sep='/')) #paste(yrs_region,collapse = "")

#convergence
diagnostics['BotTemp3d','status',sp]<-ifelse(test = is.null(fit) == T | is.null(fit$parameter_estimates$max_gradient),"no_convergence","check_gradient")
#max gradient

if (diagnostics['BotTemp3d','status',sp]=="check_gradient") {
  
  diagnostics['BotTemp3d','maxgradient',sp]<-fit$parameter_estimates$max_gradient
  #AIC
  diagnostics['BotTemp3d','aic',sp]<-round(fit$parameter_estimates$AIC[1],3)
  #JNLL
  diagnostics['BotTemp3d','jnll',sp]<-round(fit$parameter_estimates$objective[1],3)
  #RMSE
  diagnostics['BotTemp3d','rmse',sp]<-round(sqrt(mean((data_geostat$CPUE_kg - fit$Report$D_i)^2,na.rm = TRUE)) / mean(data_geostat$CPUE_kg,na.rm = TRUE),3)
  #DEVIANCE EXPLAINED
  diagnostics['BotTemp3d','dev',sp]<-round(fit$Report$deviance,3)
  
}

#progress bar
pctgy <- paste0(round(which(spp == sp)/length(spp) *100, 0), "% completed")
setWinProgressBar(py, which(spp == sp), label = pctgy) # The label will override the label set on the

#}


#percent of dev explained (only we can calculated if NULL model fitted)
# diagnostics[,'pct.dev',sp]<-(1 - as.numeric(diagnostics[,'dev',sp])/as.numeric(diagnostics['null','dev',sp]))*100
# 
# #df diagnostics
# df.diagnostics<-data.frame(diagnostics[,,sp])
# df.diagnostics[,'pct.dev']<-(1 - as.numeric(df.diagnostics[,'dev'])/as.numeric(df.diagnostics['null','dev']))*100
#df.diagnostics$pct.dev<-(1 - as.numeric(diagnostics[,'dev',sp])/as.numeric(diagnostics['null','dev',sp]))*100


#save RDS effects and diagnostics - move out of the loop when end testing pcod
saveRDS(data.frame(df.diagnostics[,,sp]),paste(out_dir,fol_region,sp,'diagnostics.RData',sep='/'))

#close process window
gc()
#}

close(py)

