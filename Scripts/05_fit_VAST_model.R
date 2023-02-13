####################################################################
####################################################################
##    
##    fit single sp VAST models using depth and SBT
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu)
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
#pack_cran<-c("TMB")

#install pacman to use p_load function - call library and if not installed, then install
# if (!('pacman' %in% installed.packages())) {
#   install.packages("pacman")}

#install coldpool to extract SBT for the EBS
if (!('VAST' %in% installed.packages())) {
  devtools::install_github("james-thorson/VAST@main", INSTALL_opts="--no-staged-install")};library(VAST)

#load/install packages
#pacman::p_load(pack_cran,character.only = TRUE)

#setwd
out_dir<-'E:/UW/Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

#range years of data
sta_y<-2002
end_y<-2022

#number of knots
knots<-'200'

#list of sp
splist<-list.dirs('./slope shelf EBS NBS VAST',full.names = FALSE)
splist<-splist[-1]

#list of models
models<-c('null',
          'depth','sbt','sp','st',
          'depth_sbt','sp_st','depth_sbt_sp','depth_sbt_st',
          'full')

#loop over species to fit models
#for (sp in sp.list) {

#Pcod example
sp<-splist[3]

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
colnames(df1)[8]<-'Temp'

#remove rows with NAs
df2<-df1[complete.cases(df1),]

#covariate data
covariate_data<-df1[,c("Lat","Lon","Year",'CPUE_kg',"Depth",'Temp')]

#regions (predefined in VAST)
region<-c("bering_sea_slope","eastern_bering_sea",'northern_bering_sea')

  #loop over models
  #for (m in models) {
  
  m<-models[2]
  
  #print year to check progress
  cat(paste("    ----- ", sp, " -----\n","       - ", m, " model\n"))  
  
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
                            #strata.limits=data.frame(STRATA = c('All_areas',"shelf","slope")),
                            #fine_scale=TRUE,
                            ObsModel = c(2,1), #c(1,1) #biomass
                            max_cells = Inf,
                            Options = c("Calculate_Range" =  T, 
                                        "Calculate_effective_area" = T)) 
  
  #Kmeans_knots-200
  if (!file.exists(paste0('./slope shelf EBS NBS VAST/',sp,'/',m,'/','Kmeans_knots-',knots,'.RData'))) {
    file.copy(paste0('./slope shelf EBS NBS VAST/',splist[1],'/',models[1],'/','Kmeans_knots-',knots,'.RData'),
              paste0('./slope shelf EBS NBS VAST/',sp,'/',m,'/','Kmeans_knots-',knots,'.RData'))
  }
  
  #formula for each model
  X1_formula<-ifelse(m %in% c('full','depth_sbt','depth_sbt_sp','depth_sbt_st'), '~log(Depth)+(log(Depth))^2+Temp',
                     ifelse(m=='depth','~log(Depth)+(log(Depth))^2',
                            ifelse(m=='sbt','~Temp',
                                   ifelse(m %in% c('null','sp','st','sp_st'),'~0'))))

  #formula for positive catch rates equal to presence/absence
  X2_formula<-X1_formula
  
  #modify settings
  #FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID"),
  #FieldConfig = c("Omega1"=2, "Epsilon1"=2, "Omega2"=2, "Epsilon2"=2), #2 factors
  #FieldConfig = c("Omega1"=3, "Epsilon1"=3, "Omega2"=3, "Epsilon2"=3), #3 factors
  settings$
  
  
  #check factors
  
  
  
  
  #fit model
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
                   covariate_data = cbind(covariate_data[,c("Lat","Lon","Depth",'Temp','Year')]), 
                   X1_formula =  X1_formula,
                   X2_formula = X2_formula, 
                   #X_gtp = X_gtp,
                   working_dir = paste0('./slope shelf EBS NBS VAST/',sp,'/',m,'/'))
  
  
  #progress bar
  pctgy <- paste0(round(which(models == m)/length(models) *100, 0), "% completed")
  setWinProgressBar(py, which(models == m), label = pctgy) # The label will override the label set on the

  }

close(py)

}