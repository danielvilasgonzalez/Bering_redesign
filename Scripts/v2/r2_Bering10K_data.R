####################################################################
####################################################################
##
##    extract sea bottom temperature (SBT) from netcdf of Bering 10K ROMS
##    netcdf downloaded from https://data.pmel.noaa.gov/aclim/thredds/
##    incorporate Temp (SBT) to Bering Sea grid
##    incorporate Temp (SBT) to data_geostat file from Bering Sea 
##    (shelf Eastern Bering Sea, slope Eastern Bering Sea, northern Bering Sea)
##    save data_geostat_temp file to fit VAST
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu)
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c('ncdf4','raster','FNN','lubridate')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#setwd
out_dir<-'C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

#range years of data
sta_y<-1982
end_y<-2022

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

#get files from google drive and set up
files<-googledrive::drive_find()
2 #for dvilasg@uw.edu

#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='Bering redesign RWP project'),'id']

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)
id.data<-files.1[which(files.1$name=='data raw'),'id']
files.2<-googledrive::drive_ls(id.data$id)

#####################################
# GET STATIONS (HAUL)
#####################################

#create directory
dir.create('./data raw/',showWarnings = FALSE)

#get haul (stations) data
file<-files.2[grep('haul',files.2$name),]
#file.id<-files.2[which(files.2$name %in% file),]

#download file
googledrive::drive_download(file=file$id,
                            path = paste0('./data raw/',file$name),
                            overwrite = TRUE)

#read csv file
haul<-readRDS(paste0('./data raw/',file$name))
dim(haul);length(unique(haul$hauljoin))

#get year and month from haul
haul$month<-month(as.POSIXlt(haul$date, format="%d/%m/%Y"))
haul$year<-year(as.POSIXlt(haul$date, format="%d/%m/%Y"))

#####################################
# BERING SEA GRIDS 
#####################################
#https://github.com/James-Thorson-NOAA/FishStatsUtils/tree/main/data
load('./extrapolation grids/eastern_bering_sea_grid.rda')
dim(eastern_bering_sea_grid)
load('./extrapolation grids/northern_bering_sea_grid.rda')
dim(northern_bering_sea_grid)
load('./extrapolation grids/bering_sea_slope_grid.rda')
dim(bering_sea_slope_grid)

eastern_bering_sea_grid<-as.data.frame(eastern_bering_sea_grid)
eastern_bering_sea_grid$region<-'EBSshelf'
northern_bering_sea_grid<-as.data.frame(northern_bering_sea_grid)
northern_bering_sea_grid$region<-'NBS'
bering_sea_slope_grid<-as.data.frame(bering_sea_slope_grid)
bering_sea_slope_grid$Stratum<-'NA'
bering_sea_slope_grid$region<-'EBSslope'

EBSgrid<-rbind(northern_bering_sea_grid,
               eastern_bering_sea_grid,
               bering_sea_slope_grid[,c("Lat","Lon","Area_in_survey_km2",'Stratum',"region")])

#####################################
# ADD DEPTH from GEBCO (just in case)
#####################################

#create directory
dir.create('./bathymetry/',showWarnings = FALSE)

#get raster depth from gebco (https://download.gebco.net/)
#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='Bathymetry'),'id']

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)
id.data<-files.1[which(files.1$name=='gebco_2022_n70.0_s50.0_w-180.0_e-155.0.asc'),]

googledrive::drive_download(file=id.data$id,
                            path = paste0('./bathymetry/',id.data$name),
                            overwrite = TRUE)

#read raster
r<-raster('./bathymetry/gebco_2022_n70.0_s50.0_w-180.0_e-155.0.asc')

#extract depth values for each station of grid
rr<-extract(r, SpatialPoints(cbind(grid.ebs$Lon,grid.ebs$Lat)))
grid.ebs$DepthGEBCO<--rr
grid.ebs$depth_m <- ifelse(is.na(grid.ebs$Depth), grid.ebs$DepthGEBCO, grid.ebs$Depth)

#same but last version of the grid
coordinates(EBSgrid)<- ~Lon+Lat
proj4string(EBSgrid)<-crs('+proj=longlat +datum=WGS84 +no_defs ')
EBSgrid$Depth<--raster::extract(r,EBSgrid)

# #####################################
# # LOOP OVER YEARS for the old grid
# #####################################
# 
# #create directory
# dir.create('./bering 10k roms/',showWarnings = FALSE)
# 
# #get id shared folder from google drive
# id.roms.folder<-files[which(files$name=='Bering 10K ROMS'),'id']
# 
# #list of files and folder of netcd files
# files.1<-googledrive::drive_ls(id.roms.folder$id)
# id.data<-files.1[which(files.1$name=='netcdf_historical'),'id']
# nc_histfiles<-googledrive::drive_ls(id.data$id)
# id.data<-files.1[which(files.1$name=='netcdf_forecast'),'id']
# nc_forfiles<-googledrive::drive_ls(id.data$id)
# files.hist<-sort(nc_histfiles$name)
# files.for<-sort(nc_forfiles$name)
# 
# #create df to store results
# grid.ebs_year<-data.frame(matrix(nrow = 0,
#                                  ncol = ncol(grid.ebs)+2))
# colnames(grid.ebs_year)<-c(colnames(grid.ebs),
#                            'Temp',
#                            'Year')
# 
# #create df to store results
# st_year<-data.frame(matrix(nrow = 0,
#                            ncol = ncol(haul)+1))
# colnames(st_year)<-c(colnames(haul),
#                      'Temp')
# 
# #loop over years to incorporate values into the Bering Sea grid
# for (y in sta_y:end_y) {
#   
#   #y<-2019
#   
#   #print year to check progress
#   cat(paste("    ----- year", y, "-----\n"))  
#   
#   #open netcdf file for each year
#   if (y %in% c(1980:1984)) {
#     f<-files.hist[1]
#     file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
#   } else  if (y %in% c(1985:1989)) {
#     f<-files.hist[2]
#     file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
#   } else  if (y %in% c(1990:1994)) {
#     f<-files.hist[3]
#     file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
#   } else  if (y %in% c(1995:1999)) {
#     f<-files.hist[4]
#     file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
#   } else if (y %in% c(2000:2004)) {
#     f<-files.hist[5]
#     file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
#   } else if (y %in% c(2005:2009)) {
#     f<-files.hist[6]
#     file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
#   } else if (y %in% c(2010:2014)) {
#     f<-files.hist[7]
#     file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
#   } else if (y %in% c(2015:2019)) {
#     f<-files.hist[8]
#     file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
#   } else if (y %in% c(2020)) {
#     f<-files.hist[9]
#     file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
#   } else if (y %in% c(2021:2022)) {
#     f<-files.for[1]
#     file.id<-nc_forfiles[which(nc_forfiles$name==f),'id']} #for year>2020 have to select projection
#   
#   #if not file, download
#   if (!(f %in% list.files('./bering 10k roms/'))) {
#     
#     #download file into temp folder
#     googledrive::drive_download(file=file.id$id,
#                                 path = paste0('./bering 10k roms/',f),
#                                 overwrite = TRUE)
#   }
#   
#   #open netcdf
#   nc<-nc_open(paste0('./bering 10k roms/',f))
#   
#   #dimensions netcdf files
#   #258 rows
#   #182 cols
#   #46956 cells
#   #259 time steps
#   
#   #get variables
#   #names(nc$var)
#   
#   #get latitude
#   lats <- ncvar_get(nc,"lat_rho")
#   
#   #get longitude
#   lons <- ncvar_get(nc,"lon_rho")
#   
#   #get SBT
#   temp<-ncvar_get(nc,'temp')
#   
#   #get time
#   t_axis<-ncvar_get(nc,"ocean_time")
#   
#   #convert time
#   time_axis <- as.POSIXct(t_axis, origin = "1900-01-01", tz = "GMT") 
#   
#   #get weekly temp slices from specific year y
#   nc_y<-ncvar_get(nc, "temp")[,,which(grepl(paste0(y,'-'),time_axis))]
#   
#   #get mean matrix for this year
#   mean_nc<-apply(nc_y,c(1,2),mean,na.rm=TRUE)
#   
#   #create dataframe with lats, lons and mean year SBT
#   df_nc<-data.frame(Lat=as.vector(lats),
#                     Lon=as.vector(lons),
#                     temp=as.vector(mean_nc))
#   
#   #longitude are in eastern. get SBT for the western hemisphere (Bering Sea). longitude greater than 180 degrees
#   df_nc1<-df_nc[which(df_nc$Lon>=180),]
#   
#   #convert eastern longitude to western values (higher). longitude should be negative
#   df_nc1$Lon<-df_nc1$Lon-180-180
#   
#   #filter values from the grid 
#   df_nc2<-subset(df_nc1,Lat >= min(grid.ebs$Lat) & Lat <= max(grid.ebs$Lat) & Lon >= min(grid.ebs$Lon) & Lon <= max(grid.ebs$Lon))
#   
#   #remove rows with NAs
#   df_nc3<-df_nc2[complete.cases(df_nc2),]
#   
#   #create spatial object from df
#   coordinates(df_nc3) <- ~ Lon + Lat
#   
#   ########################################
#   # for GRID duplicate (EBSgrid is the last version and grid.ebs is the old)
#   ########################################
#   
#   #create spatial object from grid
#   spg <- grid.ebs
#   coordinates(spg) <- ~ Lon + Lat
#   
#   #get the nearests points from one df to other df
#   nn<-get.knnx(coordinates(df_nc3),coordinates(spg),1)
#   nc_index<-nn$nn.index[,1]
#   
#   #get SBT
#   temps<-as.data.frame(df_nc3)$temp[nc_index]
#   grid.ebs$Temp<-temps
#   grid.ebs$Year<-y
#   
#   #incorporate SBT
#   grid.ebs_year<-rbind(grid.ebs_year,grid.ebs)
#   
#   ########################################
#   # for STATIONS
#   ########################################
#   
#   #create spatial object from grid
#   haul1<-subset(haul,year==y)
#   
#   if (nrow(haul1)==0) {
#     
#     haul2<-data.frame(matrix(nrow = nrow(grid.ebs),ncol = ncol(haul1)))
#     colnames(haul2)<-colnames(haul1)
#     haul2$survey_name<-grid.ebs$STRATA
#     haul2$depth_m<-grid.ebs$depth_m
#     haul2$lat_end<-grid.ebs$Lat
#     haul2$lon_end<-grid.ebs$Lon
#     haul2$lat_start<-grid.ebs$Lat
#     haul2$lon_start<-grid.ebs$Lon
#     haul2$year<-grid.ebs$Year
#     haul2$Temp<-grid.ebs$Temp
#     
#     #incorporate df with SBT
#     st_year<-rbind(st_year,haul2)
#     
#   } else {
#     
#     spg <- haul1
#     coordinates(spg) <- ~ lon_start + lat_start
#     
#     #get the nearests points from one df to other df
#     nn<-get.knnx(coordinates(df_nc3),coordinates(spg),1)
#     nc_index<-nn$nn.index[,1]
#     
#     #get SBT
#     temps<-as.data.frame(df_nc3)$temp[nc_index]
#     haul1$Temp<-temps
#     
#     #incorporate df with SBT
#     st_year<-rbind(st_year,haul1)}
#   
#   #close netcdf file
#   nc_close(nc)
#   
# }
# 
# #save grid Bering Sea with SBT and depth as dataframe
# saveRDS(grid.ebs_year,'./data processed/grid_slope_shelf_EBS_NBS_covariate_data.rds')
# 
# #modify survey labels
# st_year$survey_name[st_year$survey_name == "EBSshelf"] <- "Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey"
# st_year$survey_name[st_year$survey_name == "EBSslope"] <- "Eastern Bering Sea Slope Bottom Trawl Survey"
# st_year$survey_name[st_year$survey_name == "NBS"] <- "Northern Bering Sea Crab/Groundfish Survey - Eastern Bering Sea Shelf Survey Extension"
# 
# #save grid Bering Sea with SBT and depth as dataframe
# saveRDS(st_year,'./data processed/hauls_slope_shelf_EBS_NBS_covariate_data.rds')

#####################################
# LOOP OVER YEARS for the last version grid
#####################################

#create directory
dir.create('./bering 10k roms/',showWarnings = FALSE)

#get id shared folder from google drive
id.roms.folder<-files[which(files$name=='Bering 10K ROMS'),'id']

#list of files and folder of netcd files
files.1<-googledrive::drive_ls(id.roms.folder$id)
id.data<-files.1[which(files.1$name=='netcdf_historical'),'id']
nc_histfiles<-googledrive::drive_ls(id.data$id)
id.data<-files.1[which(files.1$name=='netcdf_forecast'),'id']
nc_forfiles<-googledrive::drive_ls(id.data$id)
files.hist<-sort(nc_histfiles$name)
files.for<-sort(nc_forfiles$name)

#create df to store results
grid.ebs_year<-data.frame(matrix(nrow = 0,
                                 ncol = ncol(EBSgrid)+2))
colnames(grid.ebs_year)<-c(colnames(EBSgrid),
                           'Temp',
                           'Year')

# #create df to store results
# st_year<-data.frame(matrix(nrow = 0,
#                            ncol = ncol(haul)+1))
# colnames(st_year)<-c(colnames(haul),
#                      'Temp')

#loop over years to incorporate values into the Bering Sea grid
for (y in 1982:2024) {
  
  #y<-2019
  
  #print year to check progress
  cat(paste("    ----- year", y, "-----\n"))  
  
  #open netcdf file for each year
  if (y %in% c(1980:1984)) {
    f<-files.hist[1]
    file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
  } else  if (y %in% c(1985:1989)) {
    f<-files.hist[2]
    file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
  } else  if (y %in% c(1990:1994)) {
    f<-files.hist[3]
    file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
  } else  if (y %in% c(1995:1999)) {
    f<-files.hist[4]
    file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
  } else if (y %in% c(2000:2004)) {
    f<-files.hist[5]
    file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
  } else if (y %in% c(2005:2009)) {
    f<-files.hist[6]
    file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
  } else if (y %in% c(2010:2014)) {
    f<-files.hist[7]
    file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
  } else if (y %in% c(2015:2019)) {
    f<-files.hist[8]
    file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
  } else if (y %in% c(2020)) {
    f<-files.hist[9]
    file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
  } else if (y %in% c(2021:2024)) {
    f<-files.for[1]
    file.id<-nc_forfiles[which(nc_forfiles$name==f),'id']} #for year>2020 have to select projection
  
  #if not file, download
  if (!(f %in% list.files('./bering 10k roms/'))) {
    
    #download file into temp folder
    googledrive::drive_download(file=file.id$id,
                                path = paste0('./bering 10k roms/',f),
                                overwrite = TRUE)
  }
  
  #open netcdf
  nc<-nc_open(paste0('./bering 10k roms/',f))
  
  #dimensions netcdf files
  #258 rows
  #182 cols
  #46956 cells
  #259 time steps
  
  #get variables
  #names(nc$var)
  
  #get latitude
  lats <- ncvar_get(nc,"lat_rho")
  
  #get longitude
  lons <- ncvar_get(nc,"lon_rho")
  
  #get SBT
  temp<-ncvar_get(nc,'temp')
  
  #get time
  t_axis<-ncvar_get(nc,"ocean_time")
  
  #convert time
  time_axis <- as.POSIXct(t_axis, origin = "1900-01-01", tz = "GMT") 
  
  #get weekly temp slices from specific year y
  nc_y<-ncvar_get(nc, "temp")[,,which(grepl(paste0(y,'-'),time_axis))]
  
  #get mean matrix for this year
  mean_nc<-apply(nc_y,c(1,2),mean,na.rm=TRUE)
  
  #create dataframe with lats, lons and mean year SBT
  df_nc<-data.frame(Lat=as.vector(lats),
                    Lon=as.vector(lons),
                    temp=as.vector(mean_nc))
  
  #longitude are in eastern. get SBT for the western hemisphere (Bering Sea). longitude greater than 180 degrees
  df_nc1<-df_nc[which(df_nc$Lon>=180),]
  
  #convert eastern longitude to western values (higher). longitude should be negative
  df_nc1$Lon<-df_nc1$Lon-180-180
  
  #filter values from the grid 
  df_nc2<-subset(df_nc1,Lat >= min(EBSgrid$Lat) & Lat <= max(EBSgrid$Lat) & Lon >= min(EBSgrid$Lon) & Lon <= max(EBSgrid$Lon))
  
  #remove rows with NAs
  df_nc3<-df_nc2[complete.cases(df_nc2),]
  
  #create spatial object from df
  coordinates(df_nc3) <- ~ Lon + Lat
  
  ########################################
  # for GRID duplicate (EBSgrid is the last version and grid.ebs is the old)
  ########################################
  
  #create spatial object from grid
  spg <- EBSgrid
  #coordinates(spg) <- ~ Lon + Lat
  
  #get the nearests points from one df to other df
  nn<-get.knnx(coordinates(df_nc3),coordinates(spg),1)
  nc_index<-nn$nn.index[,1]
  
  #get SBT
  temps<-as.data.frame(df_nc3)$temp[nc_index]
  EBSgrid$Temp<-temps
  EBSgrid$Year<-y
  
  #incorporate SBT
  grid.ebs_year<-rbind(grid.ebs_year,as.data.frame(EBSgrid))
  
  #   ########################################
  #   # for STATIONS
  #   ########################################
  #   
  #   #create spatial object from grid
  #   haul1<-subset(haul,year==y)
  #   
  #   if (nrow(haul1)==0) {
  #     
  #     haul2<-data.frame(matrix(nrow = nrow(grid.ebs),ncol = ncol(haul1)))
  #     colnames(haul2)<-colnames(haul1)
  #     haul2$survey_name<-grid.ebs$STRATA
  #     haul2$depth_m<-grid.ebs$depth_m
  #     haul2$lat_end<-grid.ebs$Lat
  #     haul2$lon_end<-grid.ebs$Lon
  #     haul2$lat_start<-grid.ebs$Lat
  #     haul2$lon_start<-grid.ebs$Lon
  #     haul2$year<-grid.ebs$Year
  #     haul2$Temp<-grid.ebs$Temp
  #     
  #     #incorporate df with SBT
  #     st_year<-rbind(st_year,haul2)
  #     
  #   } else {
  #     
  #     spg <- haul1
  #     coordinates(spg) <- ~ lon_start + lat_start
  #     
  #     #get the nearests points from one df to other df
  #     nn<-get.knnx(coordinates(df_nc3),coordinates(spg),1)
  #     nc_index<-nn$nn.index[,1]
  #     
  #     #get SBT
  #     temps<-as.data.frame(df_nc3)$temp[nc_index]
  #     haul1$Temp<-temps
  #     
  #     #incorporate df with SBT
  #     st_year<-rbind(st_year,haul1)}
  #   
  #   #close netcdf file
  #   nc_close(nc)
  #   
}

#save grid Bering Sea with SBT and depth as dataframe
save(grid.ebs_year,file = './data processed/grid_EBS_NBS.RData')

#####################################
# LOOP OVER SPP
#####################################

#get species name
spp<-list.dirs('./data processed/species/',full.names = FALSE,recursive = FALSE)

#loop over species to add SBT to data_geostat
for (sp in spp) {
  
  #sp<-spp[3]
  
  #print species to check progress
  cat(paste(" ############# ", sp, " #############\n"))
  
  #open data_geostat
  df1<-readRDS(paste0('./data processed/species/',sp,'/data_geostat.rds'))
  
  #create df to store results
  df1_temp<-data.frame(matrix(nrow=0,
                              ncol=ncol(df1)+4))
  colnames(df1_temp)<-c(colnames(df1),"Temp")
  
  for (y in sta_y:end_y) {
    
    #y<-2019
    
    #print year to check progress
    cat(paste("    ---- year", y, "----\n"))
    
    #subset df by year
    df2<-subset(df1,year==y)
    
    #if no data for that year, use stations file
    if (nrow(df2)==0) {
      
      #filter year and remove negative depth values
      st_year1<-subset(grid.ebs_year,Year==y & Depth > 0)
      df3<-data.frame(matrix(nrow = nrow(st_year1),ncol = ncol(df2)))
      colnames(df3)<-colnames(df2)
      df3$scientific_name<-sp
      df3$year<-y
      df3$lat_start<-st_year1$Lat
      df3$lat_end<-st_year1$Lat
      df3$lon_start<-st_year1$Lon
      df3$lon_end<-st_year1$Lon
      df3$depth_m<-st_year1$Depth
      df3$Temp<-st_year1$Temp
      df3$bottom_temp_c<-st_year1$Temp
      df2<-df3
      
    }
    
    #open netcdf file for each year
    if (y %in% c(1980:1984)) {
      f<-files.hist[1]
      file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
    } else  if (y %in% c(1985:1989)) {
      f<-files.hist[2]
      file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
    } else  if (y %in% c(1990:1994)) {
      f<-files.hist[3]
      file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
    } else  if (y %in% c(1995:1999)) {
      f<-files.hist[4]
      file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
    } else if (y %in% c(2000:2004)) {
      f<-files.hist[5]
      file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
    } else if (y %in% c(2005:2009)) {
      f<-files.hist[6]
      file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
    } else if (y %in% c(2010:2014)) {
      f<-files.hist[7]
      file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
    } else if (y %in% c(2015:2019)) {
      f<-files.hist[8]
      file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
    } else if (y %in% c(2020)) {
      f<-files.hist[9]
      file.id<-nc_histfiles[which(nc_histfiles$name==f),'id']
    } else if (y %in% c(2021:2022)) {
      f<-files.for[1]
      file.id<-nc_forfiles[which(nc_forfiles$name==f),'id']} #for year>2020 have to select projection
    
    #if not file, download
    # if (!(f %in% list.files('./Data/Bering 10K ROMS/'))) {
    #   
    #   #download file into temp folder
    #   googledrive::drive_download(file=file.id$id,
    #                               path = paste0('./Data/Bering 10K ROMS/',f),
    #                               overwrite = TRUE)
    # }
    
    #open netcdf
    nc<-nc_open(paste0('./bering 10k roms/',f))
    
    #dimensions netcdf file
    #258 rows
    #182 cols
    #46956 cells
    #259 time steps
    
    #get variables
    #names(nc$var)
    
    #get latitude
    lats <- ncvar_get(nc,"lat_rho")
    
    #get longitude
    lons <- ncvar_get(nc,"lon_rho")
    
    #get SBT
    temp<-ncvar_get(nc,'temp')
    
    #get time
    t_axis<-ncvar_get(nc,"ocean_time")
    
    #convert time
    time_axis <- as.POSIXct(t_axis, origin = "1900-01-01", tz = "GMT") 
    
    #get weekly temp slices from specific year y
    nc_y<-ncvar_get(nc, "temp")[,,which(grepl(paste0(y,'-'),time_axis))]
    
    #get mean matrix for this year
    mean_nc<-apply(nc_y,c(1,2),mean,na.rm=TRUE)
    
    #create dataframe with lats, lons and mean year SBT
    df_nc<-data.frame(Lat=as.vector(lats),
                      Lon=as.vector(lons),
                      temp=as.vector(mean_nc))
    
    #longitude are in eastern. get SBT for the western hemisphere (Bering Sea). longitude greater than 180 degrees
    df_nc1<-df_nc[which(df_nc$Lon>=180),]
    
    #convert eastern longitude to western values (higher). longitude should be negative
    df_nc1$Lon<-df_nc1$Lon-180-180
    
    #filter values from the grid 
    df_nc2<-subset(df_nc1,Lat >= min(df1$lat_start) & Lat <= max(df1$lat_start) & Lon >= min(df1$lon_start) & Lon <= max(df1$lon_start))
    
    #remove NA rows
    df_nc3<-df_nc2[complete.cases(df_nc2),]
    
    #create a spatial object
    coordinates(df_nc3) <- ~ Lon + Lat
    #plot(df_nc3)
    #as.data.frame(df_nc)$temp[nc_index]
    spg <- df2
    coordinates(spg) <- ~ lon_start + lat_start
    
    #get the nearests points from one df to other df
    nn<-get.knnx(coordinates(df_nc3),coordinates(spg),1)
    nc_index<-nn$nn.index[,1]
    
    #get SBT
    temps<-as.data.frame(df_nc3)$temp[nc_index]
    df2$Temp<-temps
    
    #add results
    df1_temp<-rbind(df1_temp,df2)
    
  }
  
  #scale covariates
  df1_temp$Scalebottom_temp_c<-scale(df1_temp$bottom_temp_c)
  df1_temp$ScaleTemp<-scale(df1_temp$Temp)
  df1_temp$LogDepth<-log(df1_temp$depth_m)
  df1_temp$ScaleLogDepth<-scale(df1_temp$LogDepth)
  
  #xx<-df1_temp[which(is.na(df1_temp$cpue_noha)),]
  #summary(xx$year)
  
  #save data_geostat with SBT
  saveRDS(df1_temp,
          paste0('./data processed/species/',sp,'/data_geostat_temp.rds'))
  
}


