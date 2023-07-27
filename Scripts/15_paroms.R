#open netcdf file for specific location
#extract values for specific locations
#getting the value X, drop stations in specific strata (middle shelf)
#convert values sea ice cover on cold-pool ice extension on the following year (month)?
#PAROMS

####################################################################
####################################################################
##    
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 

#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c('sp','FNN','ggplot2','ncdf4')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#setwd
out_dir<-'C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

#####################################
# LAST VERSION OF GRIDS 
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
# GET STATIONS (HAUL)
#####################################

library(lubridate)

#read csv file
haul<-readRDS('C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/data raw/afsc_haul_raw_2023_2_21.rds')
dim(haul);length(unique(haul$hauljoin))

#get year and month from haul
haul$month<-month(as.POSIXlt(haul$date, format="%d/%m/%Y"))
haul$year<-year(as.POSIXlt(haul$date, format="%d/%m/%Y"))

#only for one year
summary(haul)

#create df to store results
st_year_hist<-data.frame(matrix(nrow = 0,
                           ncol = ncol(haul)+1))
colnames(st_year_hist)<-c(colnames(haul),
                     'aice')

st_year_proj<-st_year_hist

#####################################
# GET ENVIR VAR
#####################################

#setwd
out_dir<-'C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/paroms/'
setwd(out_dir)


#list files
lf<-list.files('./',pattern = '.nc')

#load package
library()

i<-1

#get grid name
grid<-lf[grepl('grid',lf)]

#open netcdf of grid
grid_paroms<-nc_open(grid)

#get latitude
lats <- ncvar_get(grid_paroms,"lat_rho")

#get longitude
lons <- ncvar_get(grid_paroms,"lon_rho")

#data nc
data_lf<-setdiff(lf,grid)

#open netcdf of seaice
nc<-nc_open(paste0('./',lf[i]))

nc_t<-gsub('arctic4_','',data_lf)
nc_t<-gsub('.nc','',nc_t)
nc_y<-sub("\\_.*", "", nc_t)
nc_m<-sub('.*_','',nc_t)
  
#vasr with ice
names(nc$var)
ices<-names(nc$var)[grep('ice',names(nc$var))]
nc$var$aice

#get sea ice coverage (% ice cover)
var<-ncvar_get(nc,'aice')

#create dataframe with lats, lons and mean year SBT
df_nc<-data.frame(Lat=as.vector(lats),
                  Lon=as.vector(lons),
                  var=as.vector(var))

dim(df_nc)

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

plot(df_nc3)

########################################
# get var at haul locations
#########################################


#for historic time to fit the OM

  #loop over years
  for (y in 1982:2022) {
    
    y<-2022
    
    #subset by year
    haul1<-subset(haul,year==y) 
    
    #get month
    mm<-unique(haul1$month)
    
    #loop over months
    for (m in mm) {
      
      m<-mm[1]
      
      #subset by month
      haul2<-subset(haul1,month==m)
    
      #create spatial object from grid
      spg <- haul2
      coordinates(spg) <- ~ lon_start + lat_start
      
      #get the nearests points from one df to other df
      nn<-get.knnx(coordinates(df_nc3),coordinates(spg),1)
      nc_index<-nn$nn.index[,1]
      
      #get SBT
      var<-as.data.frame(df_nc3)$var[nc_index]
      haul1$aice<-var
      
      #plot
      ggplot()+
        geom_point(data=haul1,aes(x=lon_start,y=lat_start,color=aice))
      
      #incorporate df with SBT
      st_year_hist<-rbind(st_year_hist,haul1)
      
      #close netcdf file
      nc_close(nc)
    }
  }


  #for projection data,
  #get sea-ice prior March, so if want to predict for 2023, we will pull out sea-ice proportion for 2022-03 

  #get lat_lon from the 
  haul2022<-subset(haul,year==y) 

  #loop over projected years
  for (y in 2023:2027) {
    
    y<-2023
    
    #get previous year march netcdf file
    #nc_file<-data_lf[grepl(paste0(y-1,'_03.'),data_lf)]
    nc_file<-data_lf[grepl(paste0(y-1,'_03.'),data_lf)]
    
    #get sea ice coverage (% ice cover)
    aice<-ncvar_get(nc,'aice')
    
    #create dataframe with lats, lons and mean year SBT
    df_nc<-data.frame(Lat=as.vector(lats),
                      Lon=as.vector(lons),
                      aice=as.vector(aice))
    
    #longitude are in eastern. get SBT for the western hemisphere (Bering Sea). longitude greater than 180 degrees
    df_nc1<-df_nc[which(df_nc$Lon>=180),]
    
    #convert eastern longitude to western values (higher). longitude should be negative
    df_nc1$Lon<-df_nc1$Lon-180-180
    
    #filter values from the grid 
    df_nc2<-subset(df_nc1,Lat >= min(EBSgrid$Lat) & Lat <= max(EBSgrid$Lat) & Lon >= min(EBSgrid$Lon) & Lon <= max(EBSgrid$Lon))
    
    #remove rows with NAs
    df_nc3<-df_nc2[complete.cases(df_nc2),]
    
    #create spatial object from df
    #coordinates(df_nc3) <- ~ Lon + Lat

    #plot
    #ggplot()+
    #  geom_point(data=df_nc3,aes(x=Lon,y=Lat,color=aice))
    
    
    
    
    
  }



# if (nrow(haul1)==0) {
#   
#   haul2<-haul1[,c('lon_start','lat_start')]
#   
#   
#   haul2<-data.frame(matrix(nrow = nrow(grid.ebs),ncol = ncol(haul1)))
#   colnames(haul2)<-colnames(haul1)
#   haul2$survey_name<-grid.ebs$STRATA
#   haul2$depth_m<-grid.ebs$depth_m
#   haul2$lat_end<-grid.ebs$Lat
#   haul2$lon_end<-grid.ebs$Lon
#   haul2$lat_start<-grid.ebs$Lat
#   haul2$lon_start<-grid.ebs$Lon
#   haul2$year<-grid.ebs$Year
#   haul2$Temp<-grid.ebs$Temp
#   
#   #incorporate df with SBT
#   st_year<-rbind(st_year,haul2)
#   
# } else {
  

  #create spatial object from grid
  haul1<-subset(haul,year==2022) #remove if we want past sea ice data (these data are not in PAROMS)

  spg <- haul1
  coordinates(spg) <- ~ lon_start + lat_start
  
  #get the nearests points from one df to other df
  nn<-get.knnx(coordinates(df_nc3),coordinates(spg),1)
  nc_index<-nn$nn.index[,1]
  
  #get SBT
  var<-as.data.frame(df_nc3)$var[nc_index]
  haul1$aice<-var
  
  #plot
  ggplot()+
    geom_point(data=haul1,aes(x=lon_start,y=lat_start,color=aice))
  
  #incorporate df with SBT
  st_year<-rbind(st_year,haul1)

  #close netcdf file
  nc_close(nc)




########################################
# get var at haul locations
#########################################

