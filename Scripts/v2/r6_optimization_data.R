####################################################################
####################################################################
##
##    Get true index from the OM, prepare data for optimization
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
pack_cran<-c("splines",'SamplingStrata','wesanderson','dplyr','sp','rgeos','scales','rnaturalearth','grid','ggplot2')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#install coldpool to extract SBT for the EBS
if (!('VAST' %in% installed.packages())) {
  devtools::install_github("james-thorson/VAST@main", INSTALL_opts="--no-staged-install")};library(VAST)

#install akgfmaps to extract shapefile of Alaska
if (!('akgfmaps' %in% installed.packages())) {
  devtools::install_github('afsc-gap-products/akgfmaps')};library(akgfmaps)

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

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

###################################
# LOAD GRID EBS (remember to keep the same order as in fit_model if multiple grids)
###################################

#load grid data
#https://github.com/James-Thorson-NOAA/FishStatsUtils/tree/main/data
#https://github.com/danielvilasgonzalez/Bering_redesign/blob/main/Scripts/04_Bering10K_data.R
load('./data processed/lastversion_grid_EBS.RData') #grid.ebs_year

#remove slope grid
#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)

#annually grid data
grid.ebs_year1<-grid.ebs_year[which(grid.ebs_year$region!='EBSslope'),]
ncells<-nrow(grid.ebs_year1[which(grid.ebs_year1$Year==1982),])
yrs<-1982:2022

#build array for temporal array to store results
temp_dens_vals <- array(NA,
                        dim = c(ncells,
                                length(yrs),
                                length(spp)),
                        dimnames = list(1:ncells,yrs,spp))

#build array for static array
static_dens_vals <- array(NA,
                          dim = c(ncells,
                                  length(c('Lat','Lon','cell','Depth','meanTemp','meanDensity','varTemp','varDensity','sumDensity','sqsumDensity',"include","meanTempF","LonE")),
                                  length(spp)),
                          dimnames = list(1:ncells,
                                          c('Lat','Lon','cell','Depth','meanTemp','meanDensity','varTemp','varDensity','sumDensity','sqsumDensity',"include","meanTempF","LonE"),
                                          spp))
#array to store indices
true_index<-array(NA,
                  dim = list(length(yrs),3,length(spp)),
                  dimnames = list(yrs,c('EBS_NBS','NBS','EBS'),spp))

# ###################################
# # LOOP OVER SPECIES
# ###################################

#loop over species
for (sp in spp) {
  
  #sp<-'Gadus macrocephalus'
  
  #print scenario to check progress
  cat(paste(" #############   Species", sp, match(sp,spp), 'out of',length(spp),  "  #############\n"))
  
  #create folder optimization data
  dir.create(paste0('./output/species/',sp))
  
  #create folder optimization data
  dir.create(paste0('./output/species/',sp,'/optimization data/'))
  
  ###################################
  # LOAD FIT OBJECTS (fit<-from VAST::fit_model()) and pr_list<-VAST::project_model()
  ###################################
  
  #fit file
  ff<-list.files(paste0('./shelf EBS NBS VAST/',sp,'/'),'fit',recursive=TRUE)
  
  #load fit file
  load(paste0('./shelf EBS NBS VAST/',sp,'/',ff)) #fit
  
  #################################################
  # ARRANGE PREDICTED DENSITIES FROM OM
  #################################################
  
  #get predicted densities for sp
  temp_dens_vals[,,sp] <- unlist(fit$Report$D_gct[, 1, as.character(yrs)]) #[kg/km2]
  
  #density_input<-temp_dens_vals
  D_gt<-unlist(fit$Report$D_gct[, 1, as.character(yrs)])
  
  #get true index for NBS_EBS, NBS and EBS
  true_index[,,sp]<-fit$Report$Index_ctl[1, , ]
  
  #dataframe of cells with predictions
  D_gt<-data.frame('cell'=c(1:fit$spatial_list$n_g),D_gt)
  
  #years simulated
  yrs<-as.numeric(yrs)
  
  #rename years predictions 
  colnames(D_gt)<-c('cell',yrs) #,project_yrs
  
  #reshape
  D_gt1<-reshape2::melt(D_gt,id=c('cell'))
  
  #merge cells and predictions
  D_gt2 <- merge(D_gt1, grid, by=c('cell'))
  
  #get map info
  mdl <- make_map_info(Region = fit$settings$Region, 
                       spatial_list = fit$spatial_list,
                       Extrapolation_List = fit$extrapolation_list)
  
  #merge cells and predictions
  D <- merge(D_gt2, mdl$PlotDF, by.x=c('cell',"Lat",'Lon'), by.y=c('x2i',"Lat",'Lon'))
  
  #rename columns
  colnames(D)<-c('cell','Lat','Lon','Year','Density','Area','Stratum','region','Include')
  
  #create biomass
  D$Biomass<-D$Density*D$Area
  D$Biomass<-drop_units(D$Biomass)
  
  #get grid data for yrs in the simulation
  grid.ebs_year2<-grid.ebs_year1[which(grid.ebs_year1$Year %in% yrs),]
  
  #merge grid data and predictions
  D1<-merge(D,grid.ebs_year2,by=c('Lat','Lon','Year'))
  D1$Biomass_sq<-(D1$Biomass)^2
  D1$Density_sq<-(D1$Density)^2
  #subset by year (maybe to change to get for the forecasted ones)
  
  #static sampling so, we want to aggregate annual predictions: mean density, mean temp, and temp var
  D2<-aggregate(cbind(Temp,Density) ~ Lat+Lon+cell+Depth, data = D1, FUN = mean, na.rm = TRUE)
  D3<-aggregate(cbind(Temp,Density) ~ Lat+Lon+cell+Depth, data = D1, FUN = var, na.rm = TRUE)
  D4<-aggregate(cbind(Density) ~ Lat+Lon+cell+Depth, data = D1, FUN = sum, na.rm = TRUE)
  D41<-aggregate(cbind(Density_sq) ~ Lat+Lon+cell+Depth, data = D1, FUN = sum, na.rm = TRUE)
  colnames(D2)[5:6]<-paste0('mean',colnames(D2)[5:6])
  colnames(D3)[5:6]<-paste0('var',colnames(D3)[5:6])
  colnames(D4)[5]<-paste0('sum',colnames(D4)[5])
  colnames(D41)[5]<-paste0('sum',colnames(D41)[5])
  
  #merge aggregate values
  D5<-merge(D2,D3,by=c('Lat','Lon','cell','Depth'))
  D51<-merge(D5,D41,by=c('Lat','Lon','cell','Depth'))
  D6<-merge(D51,D4,by=c('Lat','Lon','cell','Depth'))
  
  #keep only cells with positive cells
  D6$include<-ifelse(D6$Depth>0,TRUE,FALSE)
  
  #convert SBT into F to get positive values only
  D6$meanTempF<-(9/5)*D6$meanTemp + 32
  
  #add longitude on eastings to get positive values
  D6$LonE<-D6$Lon+180+180
  
  #get predictions for sp
  static_dens_vals[,,sp] <- unlist(D6)
  
  #only one species
  tdf<-temp_dens_vals[,,sp]
  
  #list OM true CPUE and true index
  CPUE_index<-list('CPUE'=D1,
                   'true_index'=true_index[,,sp])
  
  #save results list
  save(D6,file=paste0('./output/species/',sp,'/optimization data/optimization_static_data.RData'))
  
  #save results list
  save(CPUE_index,file=paste0('./output/species/',sp,'/optimization data/OM_CPUE_index.RData'))
  
  #save results list
  save(tdf,file=paste0('./output/species/',sp,'/optimization data/fit_temporal_data.RData'))
}

#join optimmization data into a one single 
for (sp in spp) {
  
  #sp<-spp[2]
  
  if (sp==spp[1]) {
    load(paste0('./output/species/',sp,'/optimization data/optimization_static_data.RData'))
    df<-D6[,c("Lat","Lon","cell","Depth","meanTemp","varTemp","sumDensity_sq","sumDensity","include","meanTempF","LonE")]  
  }
  
  load(paste0('./output/species/',sp,'/optimization data/optimization_static_data.RData'))
  dens<-data.frame(D6$sumDensity,D6$sumDensity_sq)
  names(dens)<-c(paste0(sp,'_sumDensity'),paste0(sp,'_sumDensity_sq'))
  df<-cbind(df,dens)
}

save(df,file=paste0('./output/species/multisp_optimization_static_data.RData'))

#################################################
# CREATE DATA SAMPLING SCENARIO BASELINE
#################################################

#EBS and NBS layer
ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")

x<-ebs_layers$survey.strata
y<-x['OBJECTID']
z<-y[which(y$OBJECTID==14),]
zz<-x[which(x$OBJECTID==14),]

#sf to sp object
#z1<-as(z, 'Spatial')
z1<-as(y, 'Spatial')

#get stations
#extract station EBS bottom trawl
st_EBS<-read.csv('./additional//ebs_nbs_temperature_full_area.csv')
st_EBS<-subset(st_EBS,year %in% c(2000:2022))
st2022<-subset(st_EBS,year==2022)
st2022<-unique(st2022$stationid)

#get mean latitude and longitude
st<-aggregate(cbind(latitude,longitude) ~ stationid + stratum, 
              FUN = mean,
              data = st_EBS)

#subset stations only appearing in 2022
st<-subset(st, stationid %in% st2022)

#select n stations on the corner
# st1<-st[which(st$stratum =='70' & st$latitude<=61.5 & st$longitude>=-168.5),] #if 9 3x3
# st$rm9<-ifelse(st$stationid %in% st1$stationid,0,1)
# st1<-st[which(st$stratum =='70' & st$latitude<=61.7 & st$longitude>=-169),] #if 16 4x4
# st$rm16<-ifelse(st$stationid %in% st1$stationid,0,1)
# st1<-st[which(st$stratum =='70' & st$latitude<=62.1 & st$longitude>=-169.6),] #if 25 5x5
# st$rm25<-ifelse(st$stationid %in% st1$stationid,0,1)
# nrow(st1) #25 stations to remove
corner<- c('GF','HG','IH','QP','JI','ON','PO')
st.corner<-paste(corner,collapse = '|')

#for plot purposes reproject
coordinates(st)<- ~ longitude + latitude
proj4string(st) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
st1<-spTransform(st,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))
st2<-as.data.frame(st1)
st2$corner<-ifelse(grepl(st.corner,st2$stationid),TRUE,FALSE)

#check location station removed
ggplot()+
  geom_sf(data=ebs_layers$survey.strata,fill = NA)+
  geom_point(data=st2,aes(x=longitude,y=latitude,color=corner))

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)

#check grid area
sum(grid$Area_in_survey_km2)
#check shapefile area
sum(x$Precise_Ar)/1000000
#sum(x$Shape_Area)

#df to spatialpoint df
coordinates(grid) <- ~ Lon + Lat
crs(grid)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#reproject coordinates for plotting purposes
D2_1<-grid
D2_2<-data.frame(D2_1)

#x and y cells
xycells<-as.integer(sqrt(dim(D2_1)[1]))

# create a template raster
r1 <- raster(ext=extent(D2_1),ncol=xycells, nrow=xycells) #c(15800,15800) 7000

#create raster
r2<-rasterize(D2_1, r1 ,field='cell')
plot(r2)

#get streatum area base on current sampling design
area_baseline<-as.data.frame(x[,c("Stratum","Precise_Ar","SURVEY")])
area_baseline$Precise_Ar<-area_baseline$Precise_Ar/1000000

#dataframe stratum and area
strata_areas <- data.frame('X1'=x$Stratum,'Area_in_survey_km2'=x$Precise_Ar/1000000)

#locations of stations
locations <- as.data.frame(st2)
head(locations)

#cell
locations$cell<-extract(r2,st)

#number of samples per strata for random sampling
y<-aggregate(locations$cell,by=list(locations$stratum),length)
yc<-aggregate(subset(locations,corner!=TRUE)[,'cell'],by=list(subset(locations,corner!=TRUE)[,'stratum']),length)
n_samples<-data.frame('stratum'=yc$Group.1,'scnbase'=y$x,'scnbase_bis'=yc$x)

# grid1<-grid
# coordinates(grid1)<- ~ Lon + Lat
# crs(grid1)<-'+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs'
# 
# xx<- as(x, 'Spatial')
# #grid over polygon to get samples grid which current baseline strata
# cell_strata<-data.frame(as.data.frame(grid1,'stratum'=over(grid1,xx)[,'Stratum']))

#list baseline strata
baseline_strata<-list(strata_areas=strata_areas,locations=locations,n_samples=n_samples,cell_strata=grid)

#save data
save(baseline_strata,file='./output/baseline_strata.RData')
