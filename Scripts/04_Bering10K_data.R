####################################################################
####################################################################
##
##    extract sea bottom temperature (SBT) from netcdf of Bering 10K ROMS
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
pack_cran<-c('ncdf4','raster','FNN')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#setwd
out_dir<-'E:/UW/Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

#range years of data
sta_y<-2002
end_y<-2022

#####################################
# BERING SEA GRIDS
#####################################
 
#create directory
dir.create('./Data/Extrapolation Grids/',showWarnings = FALSE)
 
#get files from google drive and set up
files<-googledrive::drive_find()
1 #for dvilasg@uw.edu

#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='Extrapolation Grids'),'id']

#grid names
grfiles<-c('SlopeThorsonGrid.csv','EBSThorsonGrid.csv','NBSThorsonGrid.csv')

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)
id.data<-files.1[which(files.1$name %in% grfiles),]

  #loop over grid files
  for (i in 1:nrow(id.data)) {
    
    #i=2
   googledrive::drive_download(file=id.data$id[i],
                               path = paste0('./Data/Extrapolation Grids/',id.data$name[i]),
                               overwrite = TRUE)
    
   gr<-read.csv(paste0('./Data/Extrapolation Grids/',id.data$name[i]))
   colnames(gr)[10]<-'Area_km2'
   
   #shapefile name
   if (id.data$name[i]=='SlopeThorsonGrid.csv') {
     grname<-'EBSslope_gr'
   } else if (id.data$name[i]=='EBSThorsonGrid.csv') {
     grname<-'EBSshelf_gr'     
   } else if (id.data$name[i]=='NBSThorsonGrid.csv') {
     grname<-'NBS_gr' 
   }
  
   #name strata
   gr$STRATA<-gsub('_gr','',grname)
   
   #assign shapefiles
   assign(grname,gr)
   
  }
 
#join grids
grid.ebs<-rbind(EBSshelf_gr,EBSslope_gr,NBS_gr)

#####################################
# ADD DEPTH
#####################################

#create directory
dir.create('./Data/Bathymetry/',showWarnings = FALSE)

#get raster depth from gebco (https://download.gebco.net/)
#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='Bathymetry'),'id']

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)
id.data<-files.1[which(files.1$name=='gebco_2022_n70.0_s50.0_w-180.0_e-155.0.asc'),]

googledrive::drive_download(file=id.data$id,
                            path = paste0('./Data/Bathymetry/',id.data$name),
                            overwrite = TRUE)

#read raster
r<-raster('./Data/Bathymetry/gebco_2022_n70.0_s50.0_w-180.0_e-155.0.asc')

#extract depth values for each station
rr<-extract(r, SpatialPoints(cbind(grid.ebs$Lon,grid.ebs$Lat)))
grid.ebs$DepthGEBCO<-rr

#####################################
# LOOP OVER YEARS
#####################################

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
                                 ncol = ncol(grid.ebs)+2))
colnames(grid.ebs_year)<-c(colnames(grid.ebs),
                           'Temp',
                           'Year')

#loop over years to incorporate values into the Bering Sea grid
for (y in sta_y:end_y) {
  
  #y<-2020
  
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
  } else if (y %in% c(2021:2022)) {
    f<-files.for[1]
    file.id<-nc_forfiles[which(nc_forfiles$name==f),'id']} #for year>2020 have to select projection
  
  #if not file, download
  if (!(f %in% list.files('./Data/Bering 10K ROMS/'))) {
    
    #download file into temp folder
    googledrive::drive_download(file=file.id$id,
                                path = paste0('./Data/Bering 10K ROMS/',f),
                                overwrite = TRUE)
  }
  
  #open netcdf
  nc<-nc_open(paste0('./Data/Bering 10K ROMS/',f))
  
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
  df_nc2<-subset(df_nc1,Lat >= min(grid.ebs$Lat) & Lat <= max(grid.ebs$Lat) & Lon >= min(grid.ebs$Lon) & Lon <= max(grid.ebs$Lon))
  
  #remove rows with NAs
  df_nc3<-df_nc2[complete.cases(df_nc2),]
  
  #create spatial object from df
  coordinates(df_nc3) <- ~ Lon + Lat
  
  #create spatial object from grid
  spg <- grid.ebs
  coordinates(spg) <- ~ Lon + Lat
  
  #get the nearests points from one df to other df
  nn<-get.knnx(coordinates(df_nc3),coordinates(spg),1)
  nc_index<-nn$nn.index[,1]
  
  #get SBT
  temps<-as.data.frame(df_nc3)$temp[nc_index]
  grid.ebs$Temp<-temps
  grid.ebs$Year<-y
  
  #incorporate SBT
  grid.ebs_year<-rbind(grid.ebs_year,grid.ebs)
  
  #close netcdf file
  nc_close(nc)
  
}

#save grid Bering Sea with SBT and depth as dataframe
saveRDS(grid.ebs_year,'./slope shelf EBS NBS VAST/grid_ebs_covariate_data.rds')

#####################################
# LOOP OVER SPP
#####################################

#get species name
splist<-list.dirs('./slope shelf EBS NBS VAST/',full.names = FALSE,recursive = FALSE)
splist<-splist[-1]

#loop over species to add SBT to data_geostat
for (sp in splist) {
  
  #sp<-splist[3]
  
  #print species to check progress
  cat(paste(" ############# ", sp, " #############\n"))
  
  #open data_geostat
  df1<-readRDS(paste0('./slope shelf EBS NBS VAST/',sp,'/data_geostat.rds'))
  
  #create df to store results
  df1_temp<-data.frame(matrix(nrow=0,
                              ncol=ncol(df1)+2))
  colnames(df1_temp)<-c("Species","Year","Lat","Lon","CPUE_kg","Survey",'Depth','LogDepth',"ScaleLogDepth","Temp","ScaleTemp")

  for (y in sta_y:end_y) {
    
    #y<-2020
    
    #print year to check progress
    cat(paste("    ---- year", y, "----\n"))
    
    #subset df by year
    df2<-subset(df1,Year==y)
    
    #if no data for that year, add NA in biomass estimates because covariate data needs SBT for all years
    if (nrow(df2)==0) {
      df2<-subset(df1,Year==y-1)
      df2$Year<-y
      df2$CPUE_kg<-NA
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
    nc<-nc_open(paste0('./Data/Bering 10K ROMS/',f))
    
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
    df_nc2<-subset(df_nc1,Lat >= min(df1$Lat) & Lat <= max(df1$Lat) & Lon >= min(df1$Lon) & Lon <= max(df1$Lon))
    
    #remove NA rows
    df_nc3<-df_nc2[complete.cases(df_nc2),]
    
    #create a spatial object
    coordinates(df_nc3) <- ~ Lon + Lat
    #plot(df_nc3)
    #as.data.frame(df_nc)$temp[nc_index]
    spg <- df2
    coordinates(spg) <- ~ Lon + Lat
    
    #get the nearests points from one df to other df
    nn<-get.knnx(coordinates(df_nc3),coordinates(spg),1)
    nc_index<-nn$nn.index[,1]
    
    #get SBT
    temps<-as.data.frame(df_nc3)$temp[nc_index]
    df2$Temp<-temps
    df2$ScaleTemp<-scale(df2$Temp)
    
    #add results
    df1_temp<-rbind(df1_temp,df2)
    
  }
  
  #save data_geostat with SBT
  saveRDS(df1_temp,
          paste0('./slope shelf EBS NBS VAST/',sp,'/data_geostat_temp.rds'))
  
}

