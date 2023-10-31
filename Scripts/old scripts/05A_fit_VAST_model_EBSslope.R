####################################################################
####################################################################
##    
##    fit single sp VAST model for the slope BS using depth
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu)
##    evaluate parameters : https://github.com/James-Thorson-NOAA/VAST/blob/main/R/make_parameters.R
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c("splines",'ggplot2')

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
out_dir<- '/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

#version VAST (cpp)
version<-'VAST_v14_0_1'

#number of knots
knots<-'300' #200

#list of sp
splist<-list.dirs('./data processed/',full.names = FALSE,recursive = FALSE)

#folder region
fol_region<-'slope EBS VAST'

#load grid
load('./extrapolation grids/lastversion_grid_EBS.RData')
load('./data processed/grid_EBS_NBS.RData')

#add grid to get prediction for simulate data on each cell of the grid (sim$b_i)
load('./extrapolation grids/eastern_bering_sea_grid.rda')
head(eastern_bering_sea_grid)
dim(eastern_bering_sea_grid)
eastern_bering_sea_grid<-subset(as.data.frame(eastern_bering_sea_grid),Stratum %in% c(50,61))
load('./extrapolation grids/bering_sea_slope_grid.rda')
names(bering_sea_slope_grid)[4]<-'Stratum'
bering_sea_slope_grid$Stratum<-99
grids<-rbind(eastern_bering_sea_grid,bering_sea_slope_grid)  
names(grids)[3]<-'Area_km2'
grids

#dir create for slope region results
dir.create(paste(out_dir,fol_region,sep='/'))

#dir create for splist
for (sp in splist) {
  dir.create(paste(out_dir,fol_region,sp,sep='/'))
}

#loop over species to fit models
#for (sp in sp.list) {

#Pcod example
sp<-'Atheresthes stomias'#'Reinhardtius hippoglossoides'#'Gadus macrocephalus'
# spp<-c('Limanda aspera',
#        'Gadus chalcogrammus',
#        'Gadus macrocephalus',
#        'Atheresthes stomias',
#        'Reinhardtius hippoglossoides',
#        'Lepidopsetta polyxystra',
#        'Hippoglossoides elassodon',
#        'Pleuronectes quadrituberculatus',
#        'Hippoglossoides robustus',
#        'Boreogadus saida',
#        'Eleginus gracilis',
#        'Anoplopoma fimbria',
#        'Chionoecetes opilio',
#        'Paralithodes platypus',
#        'Paralithodes camtschaticus')

#read data_geostat_temp file
df1<-readRDS(paste0('./data processed/species/',sp,'/data_geostat_temp.rds'))

#for slope data
df11<-subset(df1,survey_name== "Eastern Bering Sea Slope Bottom Trawl Survey" | survey_name== "Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey" & depth_m>=100)

#yrs
yrs<-c(2002,2004,2008,2010,2012,2016)

#df1<-readRDS(paste0('./data processed/',sp,'/data_geostat_temp.rds'))
#select rows and rename
df2<-df1[,c("lat_start","lon_start","year",'scientific_name','weight_kg','effort','depth_m','LogDepth',"ScaleLogDepth",'Scalebottom_temp_c','bottom_temp_c','survey_name')]
colnames(df2)<-c('Lat','Lon','Year','Species','CPUE_kg','Effort','Depth','LogDepth','ScaleLogDepth','ScaleBotTemp','BotTemp','Region')

#data geostat
df3<-subset(df2,Region %in% c("Eastern Bering Sea Slope Bottom Trawl Survey",'Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey') & Year %in% yrs & Depth>=100)
yrs_region<-range(df3$Year)
data_geostat<-df3[complete.cases(df3$CPUE_kg),]
summary(data_geostat)

#covariate data - filter by year and complete cases for env variables
#covariate_data<-subset(df2,Year>=yrs_region[1] & Year<=yrs_region[2])
covariate_data<-df2[complete.cases(df2[,c('ScaleLogDepth')]),]
covariate_data$Year<-NA

#add grid to get prediction for simulate data on each cell of the grid (sim$b_i)
grid_ebs<-grid.ebs_year[which(grid.ebs_year$region == 'EBSslope' & grid.ebs_year$Year %in% yrs),]
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

#ha to km2 ------ so kg/km2
data_geostat$Effort<-data_geostat$Effort/100

#rbind grid and data_geostat to get prediction into grid values when simulating data
data_geostat1<-rbind(data_geostat[,c("Lat","Lon","Year","Species","CPUE_kg","Effort","Depth","BotTemp","Region")],
                     grid_df)

data_geostat1<-data_geostat1[,c("Lat","Lon","Year","Species","CPUE_kg","Effort","Depth","BotTemp","Region")]

#to get predictions in locations but not influencing fit
pred_TF <- rep(1, nrow(data_geostat1))
pred_TF[1:nrow(data_geostat)] <- 0

#create folder
dir.create(paste0('./slope EBS VAST/',sp))

#save data
saveRDS(data_geostat1,paste(out_dir,fol_region,sp,'data_geostat_temp.rds',sep='/'))

##################
#explore data

ggplot()+
  geom_point(data=data_geostat1,aes(x=Lon,y=Lat))+
  facet_wrap(~Year)


###################

#regions (predefined in VAST)
#region<-c("bering_sea_slope")
region<-'user'
  #loop over models
  #for (m in models) {

  #m<-models[1]

  
  #print year to check progress
  #cat(paste("\n","    ----- ", sp, " -----\n","       - ", m, " model\n"))  
  
  #create folder to store results
  dir.create(paste(out_dir,fol_region,sp,sep='/'),
             showWarnings = FALSE)
    
  #get percentage of encounters
  xall<-aggregate(data_geostat1$Year,by=list(data_geostat1$Year),FUN=length)
  x0<-subset(data_geostat1,CPUE_kg!=0)
  x0<-aggregate(x0$Year,by=list(x0$Year),FUN=length)
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
  settings <- make_settings(n_x=knots,#knots, 
                            Region=region, #c("bering_sea_slope","eastern_bering_sea",'northern_bering_sea'
                            purpose="index2", 
                            bias.correct=FALSE,
                            knot_method='grid',
                            use_anisotropy=TRUE, #TRUE
                            RhoConfig=c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0), #RhoConfig=c("Beta1"=2,"Beta2"=2,"Epsilon1"=4,"Epsilon2"=4), 
                              #FieldConfig = matrix( c("IID","IID",'IID',"Identity","IID","IID",'IID',"Identity"), #c("IID","IID",0,"Identity", "IID","IID",0,"Identity"), 
                              #                       ncol=2, 
                              #                       nrow=4, 
                              #                       dimnames=list(c("Omega","Epsilon","Beta","Epsilon_year"),c("Component_1","Component_2"))),
                            #FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID"),
                            #RhoConfig=c("Beta1"=4,"Beta2"=4,"Epsilon1"=0,"Epsilon2"=0), # Change Beta1 to AR1, to allow linear covariate effect
                            Version = version,
                            #fine_scale=TRUE,
                            ObsModel = c(2,1),#c(2,1), #c(1,1) #biomass
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
  if (!file.exists(paste(out_dir,fol_region,sp,'Kmeans_knots-',knots,'.RData',sep='/')) ) {
    file.copy(paste(out_dir,fol_region,sp,'Kmeans_knots-',knots,'.RData',sep='/'),
              paste(out_dir,fol_region,sp,'Kmeans_knots-',knots,'.RData',sep='/'))}
  
  #formula and predictors settings for each model
  X1_formula<-' ~ bs(ScaleLogDepth, degree=3, intercept=FALSE)'
  X1config_cp = array( c(1,1), dim=c(1,1))
  
  #X1_formula<-' ~ 0'
  #X1config_cp = NULL
  
  #formula for positive catch rates equal to presence/absence
  X2_formula<-X1_formula
  
  #predictor settings
  X2config_cp = X1config_cp
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
  st<-Sys.time()
  #fit model #### ADD TryCatch{(),}
  fit <- tryCatch( {fit_model(settings=settings,
                               Lat_i=data_geostat1$Lat, 
                               Lon_i=data_geostat1$Lon,
                               t_i=data_geostat1$Year,
                               b_i=data_geostat1$CPUE_kg,
                               c_iz = as.numeric(factor(data_geostat1$Species))-1,
                               a_i=data_geostat1$Effort,
                               input_grid=grids,
                               getJointPrecision = TRUE,
                               test_fit=FALSE,
                               create_strata_per_region = TRUE,  
                               covariate_data = covariate_data[,c('Year',"Lat","Lon","ScaleLogDepth","LogDepth",'ScaleBotTemp','BotTemp',"CPUE_kg")], 
                               X1_formula =  X1_formula,
                               X2_formula = X2_formula, 
                               #newtonsteps = 1,
                               PredTF_i = pred_TF,
                               #X_gtp = X_gtp,
                               working_dir = paste(out_dir,fol_region,sp,'/',sep='/'))},
      error = function(cond) {
      message("Did not converge. Here's the original error message:")
      message(cond)
      cat(paste(" ERROR MODEL",m,"IS NOT CONVERGING "))  
      # Choose a return value in case of error
      return(NULL)
    })
  
  check_fit(fit$parameter_estimates)
  plot(fit)
  fit$Report$D_gct[1,,]
  
  end<-Sys.time()
  
#}

  load('./extrapolation grids/bering_sea_slope_grid.rda')
  grid<-as.data.frame(rbind(data.frame(bering_sea_slope_grid,region='SLOPE')))
  grid$cell<-1:nrow(grid)
  
  
x<-check_fit(parameter_estimates = fit$parameter_estimates)
fit$parameter_estimates$par
fit$parameter_estimates$Convergence_check
Sim1<-simulate_data(fit)
sim_bio <-matrix(data = Sim1$b_i[pred_TF == 1], #kg
                                    nrow = nrow(grid), 
                                    ncol = length(unique(yrs)))
#project_model example
project_model(x = fit,n_proj = 5)



