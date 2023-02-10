####################################################################
####################################################################
##
##    get raw data from the EBS and NBS 
##    (shelf Eastern Bering Sea, slope Eastern Bering Sea, northern Bering Sea)
##    add sp names
##    include depth from GEBCO
##    create data_geostat file to fit VAST
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu)
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c('googledrive','raster')

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

#get files from google drive and set up
files<-googledrive::drive_find()
1 #for dvilasg@uw.edu

#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='Bering redesign RWP project'),'id']

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)
id.data<-files.1[which(files.1$name=='data raw'),'id']
files.2<-googledrive::drive_ls(id.data$id)

#####################################
# SP CODE
#####################################

#create directory
dir.create('./Data/Additional/',showWarnings = FALSE)

#get species file
files.bering<-files.2$name
file<-files.bering[grep('species',files.bering)]
file.id<-files.2[which(files.2$name==file),'id']

#download file into temp folder
googledrive::drive_download(file=file.id$id,
                            path = paste0('./Data/Additional/',file),
                            overwrite = TRUE)

#read csv file
df.sp<-read.csv(paste0('./Data/Additional/',file))
df.sp1<-df.sp[,c('SPECIES_CODE','SPECIES_NAME')]

#sp list
sp.list<-c('Limanda aspera','Gadus chalcogrammus','Gadus macrocephalus','Atheresthes stomias','Reinhardtius hippoglossoides',
           'Lepidopsetta polyxystra','Hippoglossoides elassodon','Pleuronectes quadrituberculatus','Hippoglossoides robustus')

#####################################
# EBS SLOPE CPUE
#####################################

#create directory
dir.create('./Data/Surveys/',showWarnings = FALSE)

#get cpue file
file<-files.bering[grep('slope',files.bering)]
file.id<-files.2[which(files.2$name==file),'id']

#download file into temp folder
googledrive::drive_download(file=file.id$id,
                            path = paste0('./Data/Surveys/',file),
                            overwrite = TRUE)

#read csv file
df.slope<-read.csv(paste0('./Data/Surveys/',file))

df.slope1<-merge(x = df.slope,
                 y = df.sp1,
                 by = 'SPECIES_CODE',
                 all.x = TRUE)

#####################################
# EBS SHELF
#####################################

#get cpue file
file<-files.bering[grep('shelf',files.bering)]
file.id<-files.2[which(files.2$name==file),'id']

#download file into temp folder
googledrive::drive_download(file=file.id$id,
                            path = paste0('./Data/Surveys/',file),
                            overwrite = TRUE)

#read csv file
df.shelf<-read.csv(paste0('./Data/Surveys/',file))

#####################################
# NBS SHELF
#####################################

#get cpue file
file<-files.bering[grep('nbs',files.bering)]
file.id<-files.2[which(files.2$name==file),'id']

#download file into temp folder
googledrive::drive_download(file=file.id$id,
                            path = paste0('./Data/Surveys/',file),
                            overwrite = TRUE)

#read csv file
df.nbs<-read.csv(paste0('./Data/Surveys/',file))

#####################################
# MERGE DATA
#####################################

#by species
df.shelf1<-subset(df.shelf, SPECIES_NAME %in% sp.list) 
df.nbs1<-subset(df.nbs, SPECIES_NAME %in% sp.list)
df.slope1<-subset(df.slope1,SPECIES_NAME %in% sp.list)

#select columns
df.shelf2<-df.shelf1[,c("SPECIES_NAME","YEAR","LATITUDE","LONGITUDE","CPUE_KGHA")]
df.slope2<-df.slope1[,c("SPECIES_NAME","YEAR","LATITUDE","LONGITUDE","WGTCPUE")]
colnames(df.slope2)[5]<-'CPUE_KGHA'
df.nbs2<-df.nbs1[,c("SPECIES_NAME","YEAR","LATITUDE","LONGITUDE","CPUE_KGHA")]

#correct slope data (slope is in ha, while others and in km2)
df.slope2$CPUE_KGHA<-df.slope2$CPUE_KGHA/100

#add survey
df.slope2$SURVEY<-'slope'
unique(df.slope2[,c("SPECIES_NAME",'YEAR')])
df.shelf2$SURVEY<-'shelf'
unique(df.shelf2[,c("SPECIES_NAME",'YEAR')])
df.nbs2$SURVEY<-'nbs'
unique(df.nbs2[,c("SPECIES_NAME",'YEAR')])

#bind dfs
df<-rbind(df.slope2,df.shelf2,df.nbs2)

#####################################
# ADD DEPTH
#####################################

#get raster depth from gebco (https://download.gebco.net/)
#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='Bathymetry'),'id']

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)
id.data<-files.1[which(files.1$name=='gebco_2022_n70.0_s50.0_w-180.0_e-155.0.asc'),]

googledrive::drive_download(file=id.data$id,
                            path = paste0('./Data/Additional/',id.data$name),
                            overwrite = TRUE)

r<-raster('./Data/Additional/gebco_2022_n70.0_s50.0_w-180.0_e-155.0.asc')

#extract depth values for each station
rr<-extract(r, SpatialPoints(cbind(df$LONGITUDE,df$LATITUDE)))
df$DEPTH<-rr

#rename cols
colnames(df)<-c('Species','Year',"Lat","Lon",'CPUE_kg','Survey','Depth')

#####################################
# BERING SEA GRIDS
#####################################

#get raster depth from gebco (https://download.gebco.net/)
#get id shared folder from google drive
id.bering.folder<-files[which(files$name=='Extrapolation Grids'),'id']

#list of files and folder
files.1<-googledrive::drive_ls(id.bering.folder$id)
id.data<-files.1[which(files.1$name=='gebco_2022_n70.0_s50.0_w-180.0_e-155.0.asc'),]

googledrive::drive_download(file=id.data$id,
                            path = paste0('./Data/Additional/',id.data$name),
                            overwrite = TRUE)


grid.slope<-read.csv('./Resources/Extrapolation Grids/SlopeThorsonGrid.csv')
colnames(grid.slope)[10]<-'Area_km2'
grid.slope$STRATA<-'slope'
grid.shelf<-read.csv('./Resources/Extrapolation Grids/EBSThorsonGrid.csv')
colnames(grid.shelf)[10]<-'Area_km2'
grid.shelf$STRATA<-'shelf'
grid.nbs<-read.csv('./Resources/Extrapolation Grids/NBSThorsonGrid.csv')
colnames(grid.nbs)[10]<-'Area_km2'
grid.nbs$STRATA<-'nbs'

#join grids
grid.ebs<-rbind(grid.slope,grid.shelf,grid.nbs)
head(grid.ebs)


