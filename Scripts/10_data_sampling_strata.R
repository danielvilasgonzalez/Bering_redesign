####################################################################
####################################################################
##
##    Run sampling optimization based on predicted densities from VAST model
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

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#setwd
out_dir<-'C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

#version VAST (cpp)
version<-'VAST_v13_1_0'

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

#fit file
ff<-'fit.RData'

###################################
# LOAD GRID EBS (remember to keep the same order as in fit_model if multiple grids)
###################################

#load grid data
#https://github.com/James-Thorson-NOAA/FishStatsUtils/tree/main/data
#https://github.com/danielvilasgonzalez/Bering_redesign/blob/main/Scripts/04_Bering10K_data.R
load('./data processed/lastversion_grid_EBS.RData') #grid.ebs_year

#remove slope grid
grid.ebs_year1<-grid.ebs_year[which(grid.ebs_year$region!='EBSslope'),]

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

# ###################################
# # LOOP OVER SPECIES
# ###################################

#loop over species
for (sp in spp) {
  
  sp<-'Gadus macrocephalus'
  
  #print scenario to check progress
  cat(paste(" #############   Species", sp, match(sp,spp), 'out of',length(spp),  "  #############\n"))
  
  ###################################
  # LOAD FIT OBJECTS (fit<-from VAST::fit_model()) and pr_list<-VAST::project_model()
  ###################################
  #load fit file
  load(paste0('./shelf EBS NBS VAST/',sp,'/',ff)) #fit

  #load list of projections
  #load(paste0('./output/species/',sp,'/fit_projection.RData')) #pr_list
  
  #n cells
  ncells<-fit$spatial_list$n_g
  
  #how manyt projected years we want
  yrs<-as.integer(dimnames(fit$Report$D_gct)[[3]])
  
      #scn<-names(pr_list)[1]
      
      #years fit+projected
      #yrs<-as.integer(dimnames(fit$Report$D_gct)[[3]])

      #################################################
      # SIMULATE DATA
      #################################################
      
      #reload DLL file from model
      #fit<-reload_model(fit)
      
      #simulate data
      #sim.data<-simulate_data(fit,
      #                        type=1) #1 measurement error or conditional #2 unconditional #3 new fixed and random #4 random effects from MLE
      
      #################################################
      # ARRANGE SIMULATED or PREDICTED DATA DATA
      #################################################
      
      #get predictions for sp
      temp_dens_vals[,,sp] <- unlist(fit$Report$D_gct[, 1, as.character(yrs)])
      
      #density_input<-temp_dens_vals
      D_gt<-unlist(fit$Report$D_gct[, 1, as.character(yrs)])
      #dim(D_gt)
      #D_gt_proj<-D_gt[,paste0(project_yrs)]
      
      #drop units
      #D_gt<-drop_units(D_gt)
      
      #dataframe of cells with predictions
      D_gt<-data.frame('cell'=c(1:fit$spatial_list$n_g),D_gt)
      
      #years simulated
      yrs<-as.numeric(yrs)
      
      #rename years predictions 
      colnames(D_gt)<-c('cell',yrs) #,project_yrs
      
      #reshape
      D_gt1<-reshape2::melt(D_gt,id=c('cell'))
      
      #get map info
      mdl <- make_map_info(Region = fit$settings$Region, 
                           spatial_list = fit$spatial_list,
                           Extrapolation_List = fit$extrapolation_list)
      
      #merge cells and predictions
      D <- merge(D_gt1, mdl$PlotDF, by.x='cell', by.y='x2i')
      
      #rename columns
      colnames(D)<-c('cell','Year','Density','Lat','Lon','Include')
      
      #get grid data for yrs in the simulation
      grid.ebs_year2<-grid.ebs_year1[which(grid.ebs_year1$Year %in% yrs),]
      
      #merge grid data and predictions
      D1<-merge(D,grid.ebs_year2,by=c('Lat','Lon','Year'))
      
      #subset by year (maybe to change to get for the forecasted ones)
      
      #static sampling so, we want to aggregate annual predictions: mean density, mean temp, and temp var
      D2<-aggregate(cbind(Temp,Density) ~ Lat+Lon+cell+Depth, data = D1, FUN = mean, na.rm = TRUE)
      D3<-aggregate(cbind(Temp,Density) ~ Lat+Lon+cell+Depth, data = D1, FUN = var, na.rm = TRUE)
      D4<-aggregate(cbind(Density) ~ Lat+Lon+cell+Depth, data = D1, FUN = sum, na.rm = TRUE)
      colnames(D2)[5:6]<-paste0('mean',colnames(D2)[5:6])
      colnames(D3)[5:6]<-paste0('var',colnames(D3)[5:6])
      colnames(D4)[5]<-paste0('sum',colnames(D4)[5])
      
      #merge aggregate values
      D5<-merge(D2,D3,by=c('Lat','Lon','cell','Depth'))
      D6<-merge(D5,D4,by=c('Lat','Lon','cell','Depth'))
      
      #add squared density
      D6$sqsumDensity<-(D6$sumDensity)^2
      
      #keep only cells with positive cells
      D6$include<-ifelse(D6$Depth>0,TRUE,FALSE)
  
      #convert SBT into F to get positive values only
      D6$meanTempF<-(9/5)*D6$meanTemp + 32
      
      #add longitude on eastings to get positive values
      D6$LonE<-D6$Lon+180+180
      
      #get predictions for sp
      static_dens_vals[,,sp] <- unlist(D6)
    
      tdf<-temp_dens_vals[,,sp]
  
    #save results list
    save(D6,file=paste0('./output/species/',sp,'/optimization_static_data.RData'))
    
    #save results list
    save(tdf,file=paste0('./output/species/',sp,'/fit_temporal_data.RData'))
}
