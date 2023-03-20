####################################################################
####################################################################
##
##    Sampling strata
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu)
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c("splines",'SamplingStrata')

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

##############################
#OBJECTS TO PLOT
#############################
# Define plot extent (through trial end error) units km
panel_extent <- data.frame(x = c(-1716559.21, -77636.05), #x = c(-1326559.21, -87636.05),
                           y = c(483099.5, 2194909.7)) #y = c(533099.5, 1894909.7))

#Alaska layer
ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")
ak_sppoly<-as(ebs_layers$akland, 'Spatial')

#create directory
dir.create('./shapefiles/',showWarnings = FALSE)

#name shapefiles 
shfiles<-c('EBSshelfThorson','NBSThorson','EBSslopeThorson')

#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='Shapefiles'),'id']

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)

#loop over shapefiles
for (i in shfiles) {
  
  #i=shfiles[1]
  
  id.data<-files.1[which(grepl(i,files.1$name)),]
  
  for (j in 1:nrow(id.data)) {
    
    #download data
    googledrive::drive_download(file=id.data$id[j],
                                path = paste0('./shapefiles/',id.data$name[j]),
                                overwrite = TRUE)
    
  }
  
  #shapefile EBS
  sh<-rgdal::readOGR(dsn='./shapefiles/',layer = i)
  
  if (i=='EBSslopeThorson') {
    
    #reproject shapefile
    proj4string(sh) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
    sh<-spTransform(sh,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))
    
  }
  
  #shapefile name
  shname<-paste0(gsub('Thorson','',i),'_sh')
  
  #assign shapefiles
  assign(shname,sh)
  
}

#merge shapefiles
bs_sh1<-terra::union(EBSshelf_sh,NBS_sh)
bs_sh<-terra::union(bs_sh1,EBSslope_sh)
nbs_sh<-NBS_sh
ebs_sh<-EBSshelf_sh

#color palette
pal<-wesanderson::wes_palette('Zissou1',21,type='continuous')

##################################################
####   Our df will have fields for:
####   domain: only one domain so the value is just 1
####   id: unique ID for each sampling cell
####   X1: strata variable 2: depth of cell (m) 
####   X2: strata variable 1: longitude in eastings (km). Because the 
####       optimization does not read in negative values, I shift the 
####       values so that the lowest value is 0
####
####   Variables used to more efficiently calcualte stratum variance 
####
####   WEIGHT: number of observed years 
####   Y1, Y2, ... : density for a given cell summed across observed years
####   Y1_SQ_SUM, Y2_SQ_SUM, ... : density-squared for a given cell, 
####           summed across observed years
##################################################

#project model example
#p1<-project_model(x = fit,n_proj = n_proj)


#get predictions
D_gt<-p1$D_gct[,1,]
#D_gt_proj<-D_gt[,paste0(project_yrs)]
D_gt<-drop_units(D_gt)
D_gt<-data.frame('cell'=c(1:fit$spatial_list$n_g),D_gt)
D_gt<-D_gt[!duplicated(as.list(D_gt))]
colnames(D_gt)<-c('cell',c(fit$year_labels,project_yrs))
D_gt1<-reshape2::melt(D_gt,id=c('cell'))

mdl <- make_map_info(Region = fit$settings$Region,
                     spatial_list = fit$spatial_list,
                     Extrapolation_List = fit$extrapolation_list)

D <- merge(D_gt1, mdl$PlotDF, by.x='cell', by.y='x2i')
#x<-unique(D[c("Lat","Lon","cell")])

D1<-aggregate(D$value, by=list(cell=D$cell,Lat=D$Lat,Lon=D$Lon), FUN=mean)
colnames(D1)[ncol(D1)]<-'density'

#read raster
r<-raster('./bathymetry/gebco_2022_n70.0_s50.0_w-180.0_e-155.0.asc')

coordinates(D1) <- ~ Lon + Lat
depth<-raster::extract(r,D1)
D1$depth<-depth
D1<-as.data.frame(D1)
D1[which(D1$depth>0),'depth']<--1
D1$depth<--D1$depth


# buildFrameDF(df = D5,
#              id = 1:500,
#              X='depth',Y='value')

n_cells<-max(D1$cell)
cell_idx <- rep(x = TRUE, times = max(D1$cell))
domain_input<-rep(1, sum(cell_idx))
stratum_var_input<-data.frame(X1 = D1$depth)

#create df
frame <- cbind(
  data.frame(domainvalue = domain_input,
             id = (1:n_cells)[cell_idx],
             stratum_var_input),
  
  matrix(data = D1$density, #D5$value
         ncol = 1,
         dimnames = list(NULL, paste0("Y1")))
  )



srs_stats <- SamplingStrata::buildStrataDF(frame)

srs_n <- as.numeric(1 * table(frame$domainvalue) / nrow(frame))
srs_var <- as.matrix(srs_stats[, paste0("S", 1:1)])^2

srs_var <- sweep(x = srs_var, 
                 MARGIN = 1, 
                 STATS = (1 - srs_n / n_cells) / srs_n, 
                 FUN = "*")

srs_cv <- sqrt(srs_var) / srs_stats[, paste0("M", 1:1)]

cv_constraints <- srs_cv

## Create CV constraint df
cv <- list()
for (spp in 1:ns_opt) 
  cv[[paste0("CV", 1)]] <- 
  as.numeric(cv_constraints[, 1])
cv[["DOM"]] <- 1:1
cv[["domainvalue"]] <- 1:1
cv <- as.data.frame(cv)



srs_n <- as.numeric(5 * table(frame$domainvalue) / n_cells)

## SRS statistics
srs_var <- srs_stats$S1^2 * (1 - 5 / 500) / 5 #srs_var <- srs_stats$S1^2 * (1 - srs_n / n_cells) / srs_n

srs_cv <- sqrt(srs_var) / srs_stats$M1

## cv is a data input to the SamplingStrata package, assign the initial 
## cv constraints 
cv2<-as.data.frame(list(DOM=rep('DOM1',1),
                        CV2=rep(0.05,1),
                        domainvalue=c(1:1)))


cv <- list()
cv[["CV1"]] <- srs_cv
cv[["DOM"]] <- 1:1
cv[["domainvalue"]] <- 1:1
cv <- as.data.frame(cv)

solution <- optimStrata(method = "continuous",
                        errors = cv2, 
                        framesamp = frame,
                        iter = 300,
                        pops = 100,
                        elitism_rate = 0.1,
                        #mut_chance = 1 / (10 + 1),
                        mut_chance = 1 / (rep(10, 1) + 1),
                        nStrata = rep(10, 1),
                        showPlot = T,
                        writeFiles = T)

## Set wd for output files, create a directory if it doesn"t exist yet
setwd(result_dir)


#https://github.com/James-Thorson-NOAA/FishStatsUtils/tree/main/data
 load('./extrapolation grids/eastern_bering_sea_grid.rda')
 dim(eastern_bering_sea_grid)
 load('./extrapolation grids/northern_bering_sea_grid.rda')
 dim(northern_bering_sea_grid)
 
 rbind(eastern_bering_sea_grid,northern_bering_sea_grid)
 
 
 summary(D)
 
 fit$spatial_list$PolygonList$NN_Extrap$nn.idx #knots for each spatial grid
 fit$spatial_list$PolygonList$NN_Extrap$nn.idx

grid<-readRDS('./data processed//grid_slope_shelf_EBS_NBS_covariate_data.rds')
grid1<-grid[,c("Year","Lon","Lat","Temp","STRATA")]
grid2<-subset(grid1,STRATA!='EBSslope' & Year %in% project_yrs)
