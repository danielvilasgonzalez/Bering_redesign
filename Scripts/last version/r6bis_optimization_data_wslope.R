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



#####################
# check the slope model that converged
#####################

#fitfiles<-list.files('./slope EBS VAST/',recursive = TRUE,pattern = 'fit.RData')
#spp<-gsub('/fit.RData','',fitfiles)
df_conv<-data.frame(spp=c(spp),conv=NA)

#prepare dataframe optimization
df_conv$EBS_NBS<-NA

for (sp in spp) {
  
  #sp<-spp[2]
  
  cat(paste0('#####  ',sp,'  #######\n'))
  
  #f<-fitfiles[1]
  if (length(list.files(paste0('./slope EBS VAST/',sp,'/'),pattern = 'fit.RData'))!=0) {
    load(paste0('./slope EBS VAST/',sp,'/fit.RData'))
  }
  
  if (length(list.files(paste0('./slope EBS VAST/',sp,'/'),pattern = 'fit.RData'))==0) {
    df_conv[which(df_conv$spp==sp),'conv']<-'non convergence'
  } else if (is.null(fit)) {
    df_conv[which(df_conv$spp==sp),'conv']<-'non convergence'
  } else if (is.null(fit$parameter_estimates$Convergence_check)) {
    df_conv[which(df_conv$spp==sp),'conv']<-fit$Report
  }else{
    df_conv[which(df_conv$spp==sp),'conv']<-fit$parameter_estimates$Convergence_check
  }
  
  #EBS+NBS fit
  if (file.exists(paste0('./shelf EBS NBS VAST//',sp,'/fit.RData'))) {
    
    #load fit file
    load(paste0('./shelf EBS NBS VAST//',sp,'/fit.RData'))
    
    #dimensions and check fit
    #dim(fit$Report$D_gct) #53464
    #check_fit(fit$parameter_estimates)
    
    if (is.null(fit)) {
      df_conv[which(df_conv$spp==sp),'EBS_NBS']<-'non convergence'
    } else if (is.null(fit$parameter_estimates$Convergence_check)) {
      df_conv[which(df_conv$spp==sp),'EBS_NBS']<-fit$Report
    }else{
      df_conv[which(df_conv$spp==sp),'EBS_NBS']<-fit$parameter_estimates$Convergence_check
    }
    
  } else {
    
    df_conv[which(df_conv$spp==sp),'EBS_NBS']<-'no model'
    
  }
  
}

#write table
write.csv(df_conv,'./tables/slope_ebsnbs_convspp.csv')

conv_spp<-
  df_conv[which(df_conv$conv=='There is no evidence that the model is not converged'),'spp']


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
grid$cell<-1:nrow(grid)
grid$cell<-as.numeric(grid$cell)

gridslope<-grid[which(grid$region=='SLP'),]

ggplot()+
  geom_point(data=grid,aes(x=Lon,y=Lat,col=cell))

#load grid all years
load(file = './data processed/grid_EBS_NBS.RData')

#yrs slope
yrs_slope<-c(2002,2004,2008,2010,2012,2016)

#remove grids deeper than XXX because they cannot be sampled
grid.ebs_year1<-grid.ebs_year[which(grid.ebs_year$Year %in% yrs_slope & grid.ebs_year$region =='EBSslope'),]
hist(unique(grid.ebs_year1$depth_m))

gridslope ##########
#400m

grid.ebs_year2<-grid.ebs_year1[which(grid.ebs_year1$Year==yrs_slope[1]),]
nrow(grid.ebs_year2)
nrow(grid.ebs_year2[which(grid.ebs_year2$depth_m > 1000 ),])
nrow(grid.ebs_year2[which(grid.ebs_year2$depth_m <= 400 ),])
grid.ebs_year2[which(grid.ebs_year2$depth_m <= 400  ),]


ggplot()+
  geom_point(data=grid.ebs_year2[which(grid.ebs_year2$depth_m <= 400  ),],aes(x=Lon,y=Lat,col=depth_m))+
  scale_fill_continuous(limits=c(0,2600),high = "#132B43", low = "#56B1F7")

ggplot()+
  geom_point(data=grid.ebs_year2,aes(x=Lon,y=Lat,col=depth_m))+
  scale_color_continuous(limits=c(0,2600),high = "#132B43", low = "#56B1F7")




#converged species by region
spp_conv_ebsnbs<-df_conv[which(df_conv$EBS_NBS=='There is no evidence that the model is not converged'),'spp']
spp_conv_slope<-df_conv[which(df_conv$conv=='There is no evidence that the model is not converged'),'spp']

#years by region
yrs_slope<-c(2002,2004,2008,2010,2012,2016)
yrs<-c(1982:2019,2021:2022)

#build array for temporal array to store results
temp_dens_vals <- array(NA,
                        dim = c(nrow(gridslope),
                                length(yrs_slope),
                                length(spp)),
                        dimnames = list(1:nrow(gridslope),yrs_slope,spp))

#array to store indices
true_index<-array(NA,
                  dim = list(length(yrs_slope),1,length(spp)),
                  dimnames = list(yrs_slope,'slope',spp))


# ###################################
# # LOOP OVER SPECIES
# ###################################

#loop over species
for (sp in spp) {
  
  cat(paste('#######################',sp,'#######################\n'))
  
  dir.create(paste0('./output/species/',sp),showWarnings = FALSE)
  dir.create(paste0('./output/species/',sp,'/optimization data/'),showWarnings = FALSE)

  
  
  if (sp %in% spp_conv_slope) {
    
    
    #load slope fit
    load(paste0('./slope EBS VAST/',sp,'/fit.RData'))
    
    #get predicted densities for sp
    temp_dens_vals[,,sp] <- unlist(fit$Report$D_gct[, 1, as.character(yrs_slope)]) #[kg/km2]
    
    #density_input<-temp_dens_vals
    D_gt<-unlist(fit$Report$D_gct[, 1, as.character(yrs_slope)])
    
    
    #get true index for NBS_EBS, NBS and EBS
    true_index[,,sp]<-fit$Report$Index_ctl[1,paste0(yrs_slope),1 ]
    
    #get map info
    mdl <- make_map_info(Region = fit$settings$Region, 
                         spatial_list = fit$spatial_list,
                         Extrapolation_List = fit$extrapolation_list)
    mdl$PlotDF$x2i<-mdl$PlotDF$x2i+53464
    
  } else {
    
    D_gt<-data.frame(matrix(0,nrow = nrow(gridslope),ncol = length(yrs_slope)))
    names(D_gt)<-as.character(yrs_slope)
    
    #get true index for NBS_EBS, NBS and EBS
    true_index[,,sp]<-0
    
    #load slope fit
    load(paste0('./slope EBS VAST/',spp_conv_slope[1],'/fit.RData'))
    
    #get map info
    mdl <- make_map_info(Region = fit$settings$Region, 
                         spatial_list = fit$spatial_list,
                         Extrapolation_List = fit$extrapolation_list)
    mdl$PlotDF$x2i<-mdl$PlotDF$x2i+53464

  }  
  
  #nrow(gridslope)
  
  #dataframe of cells with predictions
  D_gt<-data.frame('cell'=c(min(gridslope$cell):max(gridslope$cell)),D_gt)
  
  #years simulated
  #yrs<-as.numeric(yrs)
  
  #rename years predictions 
  colnames(D_gt)<-c('cell',yrs_slope) #,project_yrs
  
  #reshape
  D_gt1<-reshape2::melt(D_gt,id=c('cell'))
  
  #merge cells and predictions
  D_gt2 <- merge(D_gt1, gridslope, by=c('cell'))
  
  #merge cells and predictions
  D <- merge(D_gt2, mdl$PlotDF, by.x=c('cell','Lat','Lon'), by.y=c('x2i','Lat','Lon'))
  
  #rename columns
  colnames(D)<-c('cell','Lat','Lon','Year','Density','Area','Stratum','region','Include')
  
  #create biomass
  D$Biomass<-D$Density*D$Area
  
  if (sp %in% spp_conv_slope) {
    D$Biomass<-drop_units(D$Biomass)
  }
  
  #get grid data for yrs in the simulation
  grid.ebs_year1<-grid.ebs_year[which(grid.ebs_year$Year %in% yrs_slope & grid.ebs_year$region =='EBSslope'),]
  grid.ebs_year1<-grid.ebs_year1[,c("Lat","Lon","Area_in_survey_km2","DepthGEBCO","Temp","Year")]
  names(grid.ebs_year1)[4]<-'Depth'
  
  #merge grid data and predictions
  D1<-merge(D,grid.ebs_year1,by=c('Lat','Lon','Year'))
  D1$Biomass_sq<-(D1$Biomass)^2
  D1$Density_sq<-(D1$Density)^2
  #subset by year (maybe to change to get for the forecasted ones)
  
  if (sp %in% c('Atheresthes stomias','Atheresthes evermanni')){
    D1<-subset(D1,Year %in% c(1991:2022))
  }
  
  
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
  #static_dens_vals[,,sp] <- unlist(D6)
  #D6[which(D6$Depth<=400),]
  
  #only one species
  tdf<-temp_dens_vals[,,sp]
  
  #list OM true CPUE and true index
  CPUE_index<-list('CPUE'=D1,
                   'true_index'=true_index[,,sp])
  
  #save results list
  save(D6,file=paste0('./output/species/',sp,'/optimization data/optimization_static_data_slope.RData'))
  
  #save results list
  save(CPUE_index,file=paste0('./output/species/',sp,'/optimization data/OM_CPUE_index_slope.RData'))
  
  #save results list
  save(tdf,file=paste0('./output/species/',sp,'/optimization data/fit_temporal_data_slope.RData'))
  
}

#join optimmization data into a one single 
for (sp in spp) {
  
  #sp<-spp_conv[2]
  
  if (sp==spp[1]) {
    load(paste0('./output/species/',sp,'/optimization data/optimization_static_data_slope.RData'))
    df<-D6[,c("Lat","Lon","cell","Depth","meanTemp","varTemp","sumDensity_sq","sumDensity","include","meanTempF","LonE")]  
  }
  
  load(paste0('./output/species/',sp,'/optimization data/optimization_static_data_slope.RData'))
  dens<-data.frame(D6$sumDensity,D6$sumDensity_sq)
  names(dens)<-c(paste0(sp,'_sumDensity'),paste0(sp,'_sumDensity_sq'))
  df<-cbind(df,dens)
}

save(df,file=paste0('./output/multisp_optimization_static_data_slope.RData'))


######################
# COMBINE DATA FROM SLOPE AND EBS_NBS FOR OPTIMIZATION
#####################

load(paste0('./output/multisp_optimization_static_data_ebsnbs.RData'))
head(df)
dim(df)
df_ebs_nbs<-df
df_ebs_nbs<-df_ebs_nbs[order(df_ebs_nbs$cell),]
load(paste0('./output/multisp_optimization_static_data_slope.RData'))
head(df)
dim(df)
df_slope<-df
df_slope<-df_slope[order(df_slope$cell),]
df<-rbind(df_ebs_nbs,df_slope)

save(df,file=paste0('./output/multisp_optimization_static_data_ebsnbs_slope.RData'))
load(paste0('./output/multisp_optimization_static_data_ebsnbs_slope.RData'))

