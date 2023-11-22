########################
##
## slope - shelf combination and exploration
##
########################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c("splines",'ggplot2','dplyr')

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
knots<-'200' #200

#list of sp
splist<-list.dirs('./data processed/',full.names = FALSE,recursive = FALSE)

#folder region
fol_region<-'slope EBS VAST'

#load grid
load('./extrapolation grids/lastversion_grid_EBS.RData')
load('./data processed/grid_EBS_NBS.RData')

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
       'Sebastes melanostictus')
 
splist<-list() 
 
for (sp in spp) {

 #example
 #sp<-'Reinhardtius hippoglossoides'
 #sp<-'Atheresthes stomias'
 
#read data_geostat_temp file
df1<-readRDS(paste0('./data processed/species/',sp,'/data_geostat.rds'))

#for slope data
df11<-subset(df1,survey_name== "Eastern Bering Sea Slope Bottom Trawl Survey")
df11$survey_name[df11$survey_name == 'Eastern Bering Sea Slope Bottom Trawl Survey'] <- 'slope'
#yrs only for slope
yrs<-unique(df11$year)

#for shelf data
df12<-subset(df1,survey_name== "Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey")
df12$survey_name[df12$survey_name == 'Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey'] <- 'EBS shelf'

#for shelf data deeper than 106 (3rd quartile)
df13<-subset(df1,survey_name== "Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey" & depth_m>=100)
df13$survey_name[df13$survey_name == 'Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey'] <- '>100 EBS shelf'

#rbind region specific df
df1<-rbind(df11,df12,df13)
df1<-subset(df1, year %in% yrs)

#store df
splist[[sp]]<-df1

} 
 
#rbind list dfs
df2<-dplyr::bind_rows(splist, .id = "column_label")

ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")
ebs_layers$survey.strata <- sf::st_transform(ebs_layers$survey.strata, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")#'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')

#plot colored deeper than 100 meters
ggplot()+
  geom_point(data=subset(df12,year==2016),aes(x=lon_start,y=lat_start,group=year))+
  geom_point(data=subset(df11,year==2016),aes(x=lon_start,y=lat_start,group=year),col='blue',alpha=0.8)+
  geom_sf(data=ebs_layers$survey.strata,fill = NA)+
  geom_point(data=subset(df12,year==2016 & depth_m >100),aes(x=lon_start,y=lat_start,group=year),col='green',alpha=0.8)+
  facet_wrap(~year)

#plot boxplot by region and species
df2$survey_name <- factor(df2$survey_name, levels = c('EBS shelf','>100 EBS shelf','slope'))
ggplot()+
  geom_boxplot(data=df2,aes(x=as.factor(year),y=log(cpue_kgha+1),group=interaction(year,survey_name,scientific_name),color=survey_name))+
  facet_wrap(~scientific_name,scales='free',nrow=5)

#df1<-readRDS(paste0('./data processed/',sp,'/data_geostat_temp.rds'))
#select rows and rename
df2<-df2[,c("lat_start","lon_start","year",'scientific_name','weight_kg','effort','depth_m','survey_name')]
colnames(df2)<-c('Lat','Lon','Year','Species','CPUE_kg','Effort','Depth','Region')

#data geostat
df3<-subset(df2,Region %in%  c(">100 EBS shelf",'slope'))
yrs_region<-unique(df3$Year)
df3<-df3[complete.cases(df3$CPUE_kg),]

for (sp in spp) {

#example
sp<-spp[5]  
  
#filter by sp
data_geostat<-subset(df3,Species==sp)

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

# grid_ebs<-grid.ebs_year[which(grid.ebs_year$region == 'EBSslope' & grid.ebs_year$Year %in% yrs),]
# grid_df<-data.frame(Lat=grid_ebs$Lat,
#                     Lon=grid_ebs$Lon,
#                     Year=grid_ebs$Year,
#                     Species=rep(sp,times=nrow(grid_ebs)),
#                     CPUE_kg=mean(data_geostat$CPUE_kg),
#                     Effort=grid_ebs$Area_in_survey_km2,
#                     Depth=grid_ebs$Depth,
#                     #BotTemp=grid_ebs$Temp,
#                     Region=grid_ebs$region,
#                     stringsAsFactors = T)
# 

#ha to km2 ------ so kg/km2
data_geostat$Effort<-data_geostat$Effort/100

# 
# #rbind grid and data_geostat to get prediction into grid values when simulating data

# data_geostat1<-rbind(data_geostat[,c("Lat","Lon","Year","Species","CPUE_kg","Effort","Depth","Region")],
#                      grids)

data_geostat1<-rbind(data_geostat[,c("Lat","Lon","Year","Species","CPUE_kg","Effort","Depth","Region")])
data_geostat1$ScaleLogDepth<-scale(log(data_geostat1$Depth))

#covariate data - filter by year and complete cases for env variables
#covariate_data<-subset(df2,Year>=yrs_region[1] & Year<=yrs_region[2])
covariate_data<-data_geostat1[complete.cases(data_geostat1[,c('Depth')]),]
covariate_data$Year<-NA

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
  geom_point(data=subset(data_geostat1,CPUE_kg!=0),aes(x=Lon,y=Lat,size=CPUE_kg,fill=CPUE_kg),color='transparent',shape=21)+
  facet_wrap(~Year)+
  theme_bw()


# check percent of zeros
ggplot(data = data_geostat, aes(CPUE_kg)) + 
  geom_histogram(bins =20,
                 aes(y = after_stat(density))) +
  facet_wrap(~Year) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw()



# Calculate the percentage of zeros for each group
percent_zeros <- data_geostat %>%
  group_by(Year) %>%
  summarize(percentage_zeros = mean(CPUE_kg == 0) * 100)

# Print the results
print(percent_zeros)


###################

#any year with 100%encounters or 0%encounters
enc100<-ifelse(0 %in% percent_zeros$percentage_zeros,TRUE,FALSE)
enc0<-ifelse(100 %in% percent_zeros$percentage_zeros,TRUE,FALSE)

#set settings based on enc100
 if(enc100==TRUE){
   obs <- c(2,3)
 }
 if(enc100==FALSE){
   obs <- c(2,1)
 }
#Specify observation model and config settings based on species encounter probability (i.e. presence of 100% or 0% encounter probability)
if(enc0==TRUE){
  rho_c <- c("Beta1"=1,"Beta2"=1,"Epsilon1"=1,"Epsilon2"=1)
}
if(enc0==FALSE){
  rho_c <- c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0)
}

#regions (predefined in VAST)
region<-'User'#c("bering_sea_slope")
#region<-'bering_sea_slope'#c("bering_sea_slope")

#create folder to store results
dir.create(paste(out_dir,fol_region,sp,sep='/'),
           showWarnings = FALSE)

#VAST model settings
settings <- make_settings(n_x=knots,#knots, 
                          Region=region, #c("bering_sea_slope","eastern_bering_sea",'northern_bering_sea'
                          purpose="index2", 
                          bias.correct=FALSE,
                          knot_method='grid',
                          use_anisotropy=TRUE, #TRUE
                          #RhoConfig=rho_c, 
                          RhoConfig=c("Beta1"=2,"Beta2"=2,"Epsilon1"=0,"Epsilon2"=0), 
                          #FieldConfig = matrix( c("IID","IID",'IID',"Identity","IID","IID",'IID',"Identity"), #c("IID","IID",0,"Identity", "IID","IID",0,"Identity"), 
                          #                       ncol=2, 
                          #                       nrow=4, 
                          #                       dimnames=list(c("Omega","Epsilon","Beta","Epsilon_year"),c("Component_1","Component_2"))),
                          #FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID"),
                          #RhoConfig=c("Beta1"=4,"Beta2"=4,"Epsilon1"=0,"Epsilon2"=0), # Change Beta1 to AR1, to allow linear covariate effect
                          Version = version,
                          #fine_scale=TRUE,
                          ObsModel = obs,#c(2,1), #c(1,1) #biomass
                          max_cells = Inf,
                          Options = c("Calculate_Range" =  F, 
                                      "Calculate_effective_area" = F)) 

# ####################
# ####create region grid
# ####################
# 
# library(sp) # 1.4.4
# library(sf) # 0.9.6
# 
# ### Method 1: use a set of lat/lon coordinates which define the
# ### outer edge of the region. For instance you might want to plot
# ### your data and simply create a region that captures it. The
# ### locator() function can be useful for this as shown
# ### below. Here we use a subset of the Eastern Bering Sea.
# 
# ### Use this to draw points around your data
# data_geostat_grid<-subset(data_geostat,Year==2016)
# plot(data_geostat_grid$Lon, data_geostat_grid$Lat)
# LL <- locator()
# saveRDS(LL, 'extent_LL.rds')
# 
# ## Take a data.frame of coordinates in longitude/latitude that
# ## define the outer limits of the region (the extent).
# LL <- readRDS('extent_LL.rds')
# region_extent <- data.frame(long=LL$x, lat=LL$y)
# str(region_extent)
# ## > 'data.frame':	42 obs. of  2 variables:
# ## $ long: num  -166 -166 -165 -165 -164 ...
# ## $ lat : num  53.9 54.1 54.2 54.6 55 ...
# 
# #### Turn it into a spatial polygon object
# ## Need to duplicate a point so that it is connected
# region_extent <- rbind(region_extent, region_extent[1,])
# ## https://www.maths.lancs.ac.uk/~rowlings/Teaching/Sheffield2013/cheatsheet.html
# poly <- Polygon(region_extent)
# polys <- Polygons(list(poly), ID='all')
# sps <- SpatialPolygons(list(polys))
# ## I think the F_AREA could be dropped here
# sps <- SpatialPolygonsDataFrame(sps, data.frame(Id=factor('all'), F_AREA=1, row.names='all'))
# proj4string(sps)<- CRS("+proj=longlat +datum=WGS84")
# sps <- spTransform(sps, CRS("+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "))
# ### Get UTM zone for conversion to UTM projection
# ## retrieves spatial bounding box from spatial data [,1] is
# ## longitude
# lon <- sum(bbox(sps)[1,])/2
# ## convert decimal degrees to utm zone for average longitude, use
# ## for new CRS
# utmzone <- floor((lon + 180)/6)+1
# crs_LL <- CRS('+proj=longlat +ellps=WGS84 +no_defs')
# sps@proj4string <- crs_LL


#Kmeans_knots-200
if (!file.exists(paste(out_dir,fol_region,sp,'Kmeans_knots-',knots,'.RData',sep='/')) ) {
  file.copy(paste(out_dir,fol_region,sp,'Kmeans_knots-',knots,'.RData',sep='/'),
            paste(out_dir,fol_region,sp,'Kmeans_knots-',knots,'.RData',sep='/'))}

#formula and predictors settings for each model
X1_formula<-' ~ bs(ScaleLogDepth, degree=2, intercept=FALSE)'
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
                            covariate_data = covariate_data[,c('Year',"Lat","Lon","ScaleLogDepth","Depth","CPUE_kg")], 
                            X1_formula =  X1_formula,
                            X2_formula = X2_formula, 
                            #newtonsteps = 3,
                            #PredTF_i = pred_TF,
                            #X_gtp = X_gtp,
                            working_dir = paste(out_dir,fol_region,sp,'/',sep='/'))},
                 error = function(cond) {
                   message("Did not converge. Here's the original error message:")
                   message(cond)
                   cat(paste(" ERROR MODEL IS NOT CONVERGING "))  
                   # Choose a return value in case of error
                   return(NULL)
                 })

# #change lower
# fit$tmb_list$Lower
# 
# Lower = fit$tmb_list$Lower
# Lower["logkappa1"]<--10
# 
# fit <- tryCatch( {fit_model(settings=settings,
#                             Lat_i=data_geostat1$Lat, 
#                             Lon_i=data_geostat1$Lon,
#                             t_i=data_geostat1$Year,
#                             b_i=data_geostat1$CPUE_kg,
#                             c_iz = as.numeric(factor(data_geostat1$Species))-1,
#                             a_i=data_geostat1$Effort,
#                             input_grid=grids,
#                             getJointPrecision = TRUE,
#                             test_fit=FALSE,
#                             lower=Lower,
#                             create_strata_per_region = TRUE,  
#                             covariate_data = covariate_data[,c('Year',"Lat","Lon","ScaleLogDepth","Depth","CPUE_kg")], 
#                             X1_formula =  X1_formula,
#                             X2_formula = X2_formula, 
#                             #newtonsteps = 3,
#                             #PredTF_i = pred_TF,
#                             #X_gtp = X_gtp,
#                             working_dir = paste(out_dir,fol_region,sp,'/',sep='/'))},
#                  error = function(cond) {
#                    message("Did not converge. Here's the original error message:")
#                    message(cond)
#                    cat(paste(" ERROR MODEL IS NOT CONVERGING "))  
#                    # Choose a return value in case of error
#                    return(NULL)
#                  })



#check_fit(fit$parameter_estimates)
plot(fit,working_dir = paste(out_dir,fol_region,sp,'/',sep='/'))

#save fit
save(list = "fit", file = paste(out_dir,fol_region,sp,'fit.RData',sep='/')) #paste(yrs_region,collapse = "")

tryCatch({
cat(paste('#######################',fit$parameter_estimates$Convergence_check,'for',sp,'#######################'))},
error = function(cond) {
  message("Did not converge. Here's the original error message:")
  message(cond)
  #cat(paste(" ERROR MODEL IS NOT CONVERGING "))  
  # Choose a return value in case of error
  return(NULL)
})


}




